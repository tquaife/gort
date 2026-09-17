#!/usr/bin/env python3
"""
A single-scattering FLIGHT BRDF, as a diagnostic for the -rt_* options.

run_brdf.py runs FLIGHT at BOUNCE_THRESHOLD 8, so results/brdf_flight.csv holds
the sum over scattering orders 1 to 8 and there is no way to take it apart.  That
is the right reference for "is gortt's BRDF correct", but it cannot answer the
question the -rt_* options are really about, which is how much of the near
infrared is multiple scattering and whether a recollision probability predicts it.

This runs the same scenes at BOUNCE_THRESHOLD 1 on the principal plane only.
Then, per (scene, band):

    BRF(all orders) / BRF(single) = 1 / ( 1 - omega * p )

is the spectral-invariant prediction, which inverts to an observed recollision
probability

    p = ( 1 - BRF(single)/BRF(all) ) / omega

that can be held against the one gortt computes from the crown's own geometry
(gortt_recollision_scale in gortt_brdf.c).  That is the check the -rt_recollision
option needs and cannot get from the data already on disk.

It writes results/brdf_flight_single.csv and touches nothing else.  In
particular it does not go near brdf_flight.csv or lidar_flight.csv -- the latter
has no checkpoints and is the one irreplaceable artefact here.

    ./run_flight_single.py --dry-run
    ./run_flight_single.py
    ./run_flight_single.py --report          # just re-read the CSV and tabulate
"""

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc

OUT_NAME = "brdf_flight_single.csv"

# The principal plane only.  The point is the single/multiple ratio, which is a
# smooth function of direction, so the polar grid would buy little for four
# times the cost.
SZAS = [30.0]


def directions():
    seen, out = set(), []
    for vza in config.PP_VZA:
        for raa in config.PP_RAA:
            key = (round(vza, 4), 0.0 if vza == 0.0 else round(raa, 4))
            if key not in seen:
                seen.add(key)
                out.append(key)
    return sorted(out)


def task_rows(scn, dirs, im_size, workroot):
    """One scene: build it once, then walk the angles.  See run_brdf.task_rows."""
    wd = os.path.join(workroot, scn.name)
    sc.make_workdir(wd)
    crowns = sc.tile_ghosts(scn, sc.crown_positions(scn, 0))
    sc.write_crowns(os.path.join(wd, "crowns.data"), scn, crowns)
    sc.write_spec(wd, config.BRDF_BANDS)
    # PRI_BAND_NO is -1 in write_reverse, which matters here and nowhere else:
    # with it set, FLIGHT silently raises a bounce threshold of 1 to 2
    # (flight.c:7381) and this whole run would be measuring the wrong thing.
    sc.write_reverse(os.path.join(wd, "reverse.data"), 1, config.BRDF_BRANCH, im_size)

    inflight = os.path.join(wd, "in_flight.data")
    spectrum = os.path.join(wd, "RESULTS", "out-spectrum.csv")

    rows = []
    for sza in SZAS:
        for vza, raa in dirs:
            sc.write_in_flight(inflight, scn, "r", sza, vza, raa,
                               len(config.BRDF_BANDS))
            sc.run_flight(wd)
            for band, value in zip(config.BRDF_BANDS, sc.read_spectrum(spectrum)):
                rows.append((scn.name, sza, vza, raa, band.name, value))
    return rows


def _worker(payload):
    name, dirs, im_size, workroot = payload
    t0 = time.time()
    rows = task_rows(config.scene_by_name(name), dirs, im_size, workroot)
    return name, rows, time.time() - t0


def _black_worker(payload):
    """One black-soil nadir run.  Module level so it can be pickled to a pool."""
    name, bounce, im_size, workroot = payload
    scn = config.scene_by_name(name)
    wd = os.path.join(workroot, "%s-b%d" % (name, bounce))
    sc.make_workdir(wd)
    sc.write_crowns(os.path.join(wd, "crowns.data"), scn,
                    sc.tile_ghosts(scn, sc.crown_positions(scn, 0)))
    sc.write_spec(wd, config.BRDF_BANDS, rsoil=0.0)
    sc.write_reverse(os.path.join(wd, "reverse.data"), bounce,
                     config.BRDF_BRANCH, im_size)
    sc.write_in_flight(os.path.join(wd, "in_flight.data"), scn, "r",
                       SZAS[0], 0.0, 0.0, len(config.BRDF_BANDS))
    sc.run_flight(wd)
    return name, bounce, sc.read_spectrum(os.path.join(wd, "RESULTS",
                                                       "out-spectrum.csv"))


def black_soil_p(im_size, jobs):
    """
    A clean measurement of the recollision probability.

    The report below takes p from the ratio of the all-orders BRF to the
    single-scattering one, and that ratio is contaminated: a photon scattered
    down onto the soil and back up counts as "multiple scattering" but is not a
    leaf recollision at all.  It shows: the p it infers differs between red and
    near infrared for the same canopy, and p is supposed to be spectrally
    invariant, so the number cannot be what it claims to be.

    With a black soil there is no such path, and

        BRF(all orders) / BRF(single) = 1 / ( 1 - omega p )

    holds for the canopy alone.  Nadir only, one realisation, orders 1 and 8:
    24 FLIGHT runs, about a minute.  Both bands are run and must give the same
    p -- that agreement is the test that this estimator is the honest one.
    """
    workroot = os.path.join(config.WORK_DIR, "blacksoil")
    os.makedirs(workroot, exist_ok=True)

    jobs_l = [(s.name, b, im_size, workroot) for s in config.scenes() for b in (1, 8)]
    got = {}
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for fut in as_completed([pool.submit(_black_worker, j) for j in jobs_l]):
            name, bounce, brf = fut.result()
            got[(name, bounce)] = brf

    print("Recollision probability with a BLACK SOIL, so the ratio is the canopy")
    print("alone.  p should not depend on the waveband: that is the check.\n")
    print("%-22s %9s %8s %8s %8s %8s"
          % ("scene", "crown_lai", "p gortt", "p red", "p nir", "p mean"))
    for scn in config.scenes():
        if (scn.name, 1) not in got or (scn.name, 8) not in got:
            continue
        ps = []
        for i, band in enumerate(config.BRDF_BANDS):
            s1, s8 = got[(scn.name, 1)][i], got[(scn.name, 8)][i]
            w = band.rleaf + band.tleaf
            ps.append((1.0 - s1 / s8) / w if s8 > 0 else float("nan"))
        cl = scn.favd * 4.0 * scn.b / 3.0
        print("%-22s %9.3f %8.3f %8.3f %8.3f %8.3f"
              % (scn.name, cl, 0.5 * cl / (1.0 + 0.5 * cl), ps[0], ps[1],
                 float(np.mean(ps))))


def report(path):
    """
    The observed recollision probability, per scene and band.

    Needs brdf_flight.csv for the all-orders side.  Both are averaged over the
    directions the two have in common at the shared solar zenith, which keeps
    the ratio a scene property rather than a directional one -- p is defined as
    an average over the canopy and has no direction of its own.
    """
    single, allo = {}, {}
    for r in csv.DictReader(open(path)):
        single.setdefault((r["scene"], r["band"], float(r["sza"])), {})[
            (float(r["vza"]), float(r["raa"]))] = float(r["brf"])
    fp = os.path.join(config.RESULTS_DIR, "brdf_flight.csv")
    if not os.path.exists(fp):
        sys.exit("run_flight_single: %s not found" % fp)
    for r in csv.DictReader(open(fp)):
        allo.setdefault((r["scene"], r["band"], float(r["sza"])), {}).setdefault(
            (float(r["vza"]), float(r["raa"])), []).append(float(r["brf"]))

    print("Observed recollision probability from FLIGHT, sza %g, principal plane."
          % SZAS[0])
    print("p = (1 - BRF_single/BRF_all) / omega, the spectral-invariant relation.\n")
    print("%-22s %6s %8s %9s %9s %7s %7s"
          % ("scene", "band", "omega", "BRF 1", "BRF 1-8", "ratio", "p"))
    for scn in config.scenes():
        for band in config.BRDF_BANDS:
            key = (scn.name, band.name, SZAS[0])
            if key not in single or key not in allo:
                continue
            keys = sorted(set(single[key]) & set(allo[key]))
            if not keys:
                continue
            s = float(np.mean([single[key][k] for k in keys]))
            a = float(np.mean([np.mean(allo[key][k]) for k in keys]))
            w = band.rleaf + band.tleaf
            print("%-22s %6s %8.2f %9.5f %9.5f %7.3f %7.3f"
                  % (scn.name, band.name, w, s, a, a / s,
                     (1.0 - s / a) / w if a > 0 else float("nan")))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--im-size", type=int, default=config.BRDF_IM_SIZE)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--scenes", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report", action="store_true",
                    help="skip the run, just tabulate an existing CSV")
    ap.add_argument("--black-soil", action="store_true",
                    help="instead, measure the recollision probability with a "
                         "black soil at nadir (24 FLIGHT runs, ~1 min)")
    args = ap.parse_args()

    out = os.path.join(config.RESULTS_DIR, OUT_NAME)
    if args.black_soil:
        return black_soil_p(args.im_size, args.jobs)
    if args.report:
        return report(out)

    scenes = config.scenes()
    if args.scenes:
        scenes = [config.scene_by_name(n) for n in args.scenes]
    dirs = directions()

    n = len(scenes) * len(SZAS) * len(dirs)
    sys.stderr.write("run_flight_single: %d FLIGHT processes, ~%.0f min at %d jobs\n"
                     % (n, n * 4.5 / max(1, args.jobs) / 60.0, args.jobs))
    if args.dry_run:
        return

    if not os.path.exists(config.flight_binary()):
        sys.exit("run_flight_single: no FLIGHT binary -- run build_flight.py first")

    os.makedirs(config.RESULTS_DIR, exist_ok=True)
    workroot = os.path.join(config.WORK_DIR, "single")
    os.makedirs(workroot, exist_ok=True)

    t0 = time.time()
    with open(out, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "sza", "vza", "raa", "band", "brf"])
        payloads = [(s.name, dirs, args.im_size, workroot) for s in scenes]
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_worker, p) for p in payloads]
            for i, fut in enumerate(as_completed(futures), 1):
                name, rows, dt = fut.result()
                for r in rows:
                    w.writerow([r[0], "%g" % r[1], "%g" % r[2], "%g" % r[3],
                                r[4], "%.8g" % r[5]])
                sys.stderr.write("run_flight_single: %2d/%d %-22s %5.0f s\n"
                                 % (i, len(payloads), name, dt))
    sys.stderr.write("run_flight_single: wrote %s in %.1f min\n"
                     % (out, (time.time() - t0) / 60.0))
    report(out)


if __name__ == "__main__":
    main()
