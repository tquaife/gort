#!/usr/bin/env python3
"""
Run gortt and FLIGHT over the BRDF angle grid and write tidy CSVs.

FLIGHT does one view direction per process -- MULT_VIEWS is compiled out and
marked "Not yet properly implemented" -- so the angular sweep is one FLIGHT
process per (scene, realisation, solar zenith, view direction).  That is the
entire cost of this comparison: about 6,900 processes at roughly 4.5 s each,
which is a couple of hours across four cores.

gortt, by contrast, takes the whole angle list on stdin in one run, so it is
evaluated twice: once on exactly FLIGHT's directions, for the difference panels
and the agreement statistics, and once on a much finer grid for smooth contours.

The work is checkpointed per (scene, realisation) into work/brdf/, so an
interrupted sweep resumes rather than restarting.  Pass --force to redo.

    ./run_brdf.py --dry-run     # how many runs, and how long
    ./run_brdf.py --quick       # one scene, one realisation, coarse -- smoke test
    ./run_brdf.py               # the full sweep
"""

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc

SECONDS_PER_FLIGHT_RUN = 4.5        # measured: im_size 150, 2 bands, 500 crowns


def task_rows(scn, realisation, directions, szas, im_size, bands, workroot):
    """
    One (scene, realisation): build the scene once, then walk the angles.

    The crowns, spectra and reverse.data are written once and left alone; only
    in_flight.data is rewritten per direction, which is all that changes.
    """
    wd = os.path.join(workroot, "%s-r%d" % (scn.name, realisation))
    sc.make_workdir(wd)
    crowns = sc.tile_ghosts(scn, sc.crown_positions(scn, realisation))
    sc.write_crowns(os.path.join(wd, "crowns.data"), scn, crowns)
    sc.write_spec(wd, bands)
    sc.write_reverse(os.path.join(wd, "reverse.data"),
                     config.BRDF_BOUNCE, config.BRDF_BRANCH, im_size)

    inflight = os.path.join(wd, "in_flight.data")
    spectrum = os.path.join(wd, "RESULTS", "out-spectrum.csv")

    rows, diag = [], None
    for sza in szas:
        for vza, raa in directions:
            sc.write_in_flight(inflight, scn, "r", sza, vza, raa, len(bands))
            out = sc.run_flight(wd)
            if diag is None:
                diag = sc.parse_flight_stdout(out)
                diag["ncrowns_incl_ghosts"] = len(crowns)
            brf = sc.read_spectrum(spectrum)
            for band, value in zip(bands, brf):
                rows.append((scn.name, realisation, sza, vza, raa, band.name, value))
    return rows, diag


def _worker(payload):
    name, realisation, directions, szas, im_size, workroot, checkpoint = payload
    scn = config.scene_by_name(name)
    t0 = time.time()
    rows, diag = task_rows(scn, realisation, directions, szas, im_size,
                           config.BRDF_BANDS, workroot)
    with open(checkpoint, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "realisation", "sza", "vza", "raa", "band", "brf"])
        for r in rows:
            w.writerow(["%s" % r[0], r[1], "%g" % r[2], "%g" % r[3], "%g" % r[4],
                        r[5], "%.8g" % r[6]])
    return name, realisation, len(rows), time.time() - t0, diag


def run_gortt(scenes, szas, flight_dirs, out_path):
    """
    gortt on both grids.  One process per (scene, band, solar zenith), because
    the leaf and soil overrides are spectrally flat so a band is a separate run,
    and because gortt takes a whole angle list at once.
    """
    dense = [(v, r) for v in config.GORTT_DENSE_VZA for r in config.GORTT_DENSE_RAA]
    with open(out_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "sza", "vza", "raa", "band", "grid", "brf"])
        for scn in scenes:
            for band in config.BRDF_BANDS:
                for sza in szas:
                    for grid_name, grid in (("flight", flight_dirs), ("dense", dense)):
                        dirs = [(v, sza, r) for v, r in grid]
                        brf = sc.gortt_brdf(scn, band, dirs)
                        for (v, r), value in zip(grid, brf):
                            w.writerow([scn.name, "%g" % sza, "%g" % v, "%g" % r,
                                        band.name, grid_name, "%.8g" % value])
            sys.stderr.write("run_brdf: gortt %s done\n" % scn.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--realisations", type=int, default=config.BRDF_REALISATIONS,
                    help="independent scene draws to average (default: %(default)s)")
    ap.add_argument("--im-size", type=int, default=config.BRDF_IM_SIZE,
                    help="FLIGHT image side; it traces im_size^2 paths "
                         "(default: %(default)s)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count(),
                    help="parallel FLIGHT workers (default: %(default)s)")
    ap.add_argument("--scenes", nargs="*", help="limit to these scenes")
    ap.add_argument("--quick", action="store_true",
                    help="one scene, one realisation, principal plane only")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the run count and time estimate, then stop")
    ap.add_argument("--force", action="store_true", help="ignore checkpoints")
    ap.add_argument("--gortt-only", action="store_true",
                    help="redo just the gortt side, which takes seconds")
    args = ap.parse_args()

    scenes = config.scenes()
    if args.scenes:
        scenes = [config.scene_by_name(n) for n in args.scenes]
    szas = config.SZAS
    directions = config.flight_directions()
    realisations = args.realisations
    im_size = args.im_size

    if args.quick:
        scenes = scenes[:1]
        realisations = 1
        im_size = 60
        szas = [30.0]
        directions = [(v, r) for v, r in directions if r in (0.0, 180.0)][::2]

    os.makedirs(config.RESULTS_DIR, exist_ok=True)
    workroot = os.path.join(config.WORK_DIR, "brdf")
    os.makedirs(workroot, exist_ok=True)
    ckpt_dir = os.path.join(workroot, "checkpoints")
    os.makedirs(ckpt_dir, exist_ok=True)

    nruns = len(scenes) * realisations * len(szas) * len(directions)
    est = nruns * SECONDS_PER_FLIGHT_RUN / max(1, args.jobs)
    sys.stderr.write("run_brdf: %d scenes x %d realisations x %d sza x %d directions "
                     "= %d FLIGHT runs\n"
                     % (len(scenes), realisations, len(szas), len(directions), nruns))
    sys.stderr.write("run_brdf: about %.0f min on %d workers\n" % (est / 60.0, args.jobs))
    if args.dry_run:
        return

    sys.stderr.write("run_brdf: gortt ...\n")
    run_gortt(scenes, szas, directions, os.path.join(config.RESULTS_DIR, "brdf_gortt.csv"))
    if args.gortt_only:
        return

    if not os.path.exists(config.flight_binary()):
        sys.exit("run_brdf: no FLIGHT binary -- run build_flight.py first")

    payloads, done = [], []
    for scn in scenes:
        for rz in range(realisations):
            ck = os.path.join(ckpt_dir, "%s-r%d.csv" % (scn.name, rz))
            if os.path.exists(ck) and not args.force:
                done.append(ck)
                continue
            payloads.append((scn.name, rz, directions, szas, im_size, workroot, ck))

    if done:
        sys.stderr.write("run_brdf: %d task(s) already checkpointed, skipping\n" % len(done))

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(_worker, p): p for p in payloads}
        for i, fut in enumerate(as_completed(futures), 1):
            name, rz, n, dt, diag = fut.result()
            done.append(futures[fut][6])
            sys.stderr.write("run_brdf: [%d/%d] %s r%d  %d rows  %.0f s  "
                             "(cover %.3f, LAI %.2f)\n"
                             % (i, len(payloads), name, rz, n, dt,
                                diag.get("frac_cover", float("nan")),
                                diag.get("scene_lai", float("nan"))))

    # merge the checkpoints
    out = os.path.join(config.RESULTS_DIR, "brdf_flight.csv")
    with open(out, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "realisation", "sza", "vza", "raa", "band", "brf"])
        for ck in sorted(done):
            with open(ck) as cf:
                rd = csv.reader(cf)
                next(rd)
                for row in rd:
                    w.writerow(row)
    sys.stderr.write("run_brdf: wrote %s in %.0f min\n" % (out, (time.time() - t0) / 60.0))


if __name__ == "__main__":
    main()
