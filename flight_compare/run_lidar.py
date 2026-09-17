#!/usr/bin/env python3
"""
Run the nadir lidar waveform in both models and write tidy CSVs.

gortt's waveform (Ni-Meister et al. 2001) is first order by construction: the
canopy return is J0 rho_v dP/dz, one scattering event per photon.  FLIGHT is
told to match that with BOUNCE_THRESHOLD 1 in reverse.data -- the field in
flightline.data that looks like it selects single scattering ('t'/'s'/'m') is
read into g_lidar_mode at flight.c:6465 and then never referenced again.  A
multiple-scattering pass at BOUNCE_THRESHOLD 8 is run as well, at fewer photons,
to show what gortt is leaving out; FLIGHT retraces the whole path for every
scattering order separately (flight.c:6878), so order 8 costs about 36 times
order 1 per photon.

Both models are nadir only.  gortt's -lidar does not read stdin at all and
exits before the per-angle loop, and FLIGHT's lidar code carries the comment
"Currently 'hardwired' to cope only with nadir and distant imagery".

Heights: FLIGHT records one-way range from path_zero, which sits REC_ALT_ST
above the ground target, so height = REC_ALT_ST - range.  gortt reports height
above ground directly.  The two are put on the same grid by construction --
scene.lidar_height_binsize() makes gortt's bins equal to FLIGHT's, and
scene.lidar_pulse_height_sigma() converts FLIGHT's pulse width in ns into the
metres of height that gortt's -lidar_pulse expects.

    ./run_lidar.py --dry-run
    ./run_lidar.py --quick
    ./run_lidar.py
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


def flight_waveform(scn, realisation, bounce, nphotons, workroot):
    """One FLIGHT pulse at nadir.  Returns (z, waveform) for the lidar band."""
    wd = os.path.join(workroot, "%s-r%d-b%d" % (scn.name, realisation, bounce))
    sc.make_workdir(wd)
    bands = [config.LIDAR_BAND]
    crowns = sc.tile_ghosts(scn, sc.crown_positions(scn, realisation))
    sc.write_crowns(os.path.join(wd, "crowns.data"), scn, crowns)
    sc.write_spec(wd, bands)
    sc.write_lidar(wd, scn, bands)
    sc.write_reverse(os.path.join(wd, "reverse.data"), bounce,
                     config.BRDF_BRANCH, config.BRDF_IM_SIZE)
    sc.write_in_flight(os.path.join(wd, "in_flight.data"), scn, "l",
                       0.0, 0.0, 0.0, len(bands), nphotons=nphotons)
    out = sc.run_flight(wd)
    z, w = sc.read_waveform(wd)
    return z, w[:, 0], sc.parse_flight_stdout(out)


def _worker(payload):
    name, realisation, bounce, nphotons, workroot = payload
    scn = config.scene_by_name(name)
    t0 = time.time()
    z, w, diag = flight_waveform(scn, realisation, bounce, nphotons, workroot)
    return name, realisation, bounce, z, w, diag, time.time() - t0


def gortt_waveforms(scenes, out_path):
    """
    gortt -lidar -lidar_waveform on the same height grid and pulse width as
    FLIGHT, plus the per-layer profile, whose P(z) column is the pure geometry.
    """
    with open(out_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "z", "return_density", "rho_v", "rho_g", "Rg"])
        for scn in scenes:
            head, rows = sc.gortt_lidar(scn, config.LIDAR_BAND, waveform=True)
            _, rho_v, rho_g = head[0], head[2], head[3]
            # the profile run carries the ground return on its header line
            phead, _ = sc.gortt_lidar(scn, config.LIDAR_BAND, waveform=False)
            Rg = phead[5]
            for z, d in rows:
                w.writerow([scn.name, "%.6g" % z, "%.8g" % d,
                            "%.8g" % rho_v, "%.8g" % rho_g, "%.8g" % Rg])
            sys.stderr.write("run_lidar: gortt %s done\n" % scn.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--realisations", type=int, default=config.LIDAR_REALISATIONS,
                    help="scene draws for the single-scattering pass (default: %(default)s)")
    ap.add_argument("--multi-realisations", type=int, default=3,
                    help="scene draws for the multiple-scattering pass, which is "
                         "~36x dearer per photon (default: %(default)s)")
    ap.add_argument("--photons", type=int, default=config.LIDAR_PHOTONS,
                    help="photons per pulse, single scattering (default: %(default)s)")
    ap.add_argument("--multi-photons", type=int, default=40000,
                    help="photons per pulse, multiple scattering (default: %(default)s)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--scenes", nargs="*")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--gortt-only", action="store_true")
    args = ap.parse_args()

    scenes = config.scenes()
    if args.scenes:
        scenes = [config.scene_by_name(n) for n in args.scenes]
    nreal, nmulti = args.realisations, args.multi_realisations
    photons, mphotons = args.photons, args.multi_photons
    if args.quick:
        scenes, nreal, nmulti = scenes[:1], 1, 1
        photons, mphotons = 20000, 10000

    os.makedirs(config.RESULTS_DIR, exist_ok=True)
    workroot = os.path.join(config.WORK_DIR, "lidar")
    os.makedirs(workroot, exist_ok=True)

    jobs = []
    for scn in scenes:
        for rz in range(nreal):
            jobs.append((scn.name, rz, config.LIDAR_BOUNCE_SINGLE, photons, workroot))
        for rz in range(nmulti):
            jobs.append((scn.name, rz, config.LIDAR_BOUNCE_MULTI, mphotons, workroot))

    sys.stderr.write("run_lidar: %d FLIGHT pulses (%d single, %d multiple)\n"
                     % (len(jobs), len(scenes) * nreal, len(scenes) * nmulti))
    if args.dry_run:
        return

    sys.stderr.write("run_lidar: gortt ...\n")
    gortt_waveforms(scenes, os.path.join(config.RESULTS_DIR, "lidar_gortt.csv"))
    if args.gortt_only:
        return

    if not os.path.exists(config.flight_binary()):
        sys.exit("run_lidar: no FLIGHT binary -- run build_flight.py first")

    out = os.path.join(config.RESULTS_DIR, "lidar_flight.csv")
    t0 = time.time()
    with open(out, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["scene", "realisation", "bounce", "z", "energy"])
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_worker, p) for p in jobs]
            for i, fut in enumerate(as_completed(futures), 1):
                name, rz, bounce, z, wave, diag, dt = fut.result()
                keep = wave > 0.0           # most of the range window is empty
                for zz, ee in zip(z[keep], wave[keep]):
                    w.writerow([name, rz, bounce, "%.6g" % zz, "%.8g" % ee])
                sys.stderr.write("run_lidar: [%d/%d] %s r%d bounce %d  %d bins  %.0f s\n"
                                 % (i, len(jobs), name, rz, bounce, int(keep.sum()), dt))
    sys.stderr.write("run_lidar: wrote %s in %.0f min\n" % (out, (time.time() - t0) / 60.0))


if __name__ == "__main__":
    main()
