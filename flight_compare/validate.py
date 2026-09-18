#!/usr/bin/env python3
"""
Check that gortt and FLIGHT have been told about the same canopy, before an
hour is spent sweeping angles with them.

Four checks, each of which has caught something real:

1. Fractional cover.  An independent Monte Carlo of our own crown list, against
   FLIGHT's find_frac_cover(), against GORT's 1 - exp(-lambda pi r^2).  This is
   what showed that FLIGHT does not replicate crowns across its periodic walls
   (see scene.tile_ghosts).

2. Nadir gap probability.  With black leaves (reflectance = transmittance = 0)
   and a white Lambertian soil, viewed and illuminated at nadir, FLIGHT's
   scene-mean BRF is exactly the probability that a vertical ray reaches the
   ground: rays that touch foliage return nothing, rays that reach the soil
   return 1, and at sza = vza = 0 the downward and upward paths are the same
   line so the ground point is always sunlit.  That is directly comparable with
   p_n0[0][0] + epgap[0][0] from `gortt -W`, and it is the decisive test that
   favd, lambda, r, b, h1 and h2 mean the same thing in both models.

   It also settles a question the LAI diagnostics cannot: FLIGHT reports scene
   LAI both with and without a crown-overlap correction, and GORT's elai
   assumes overlapping crowns attenuate twice.  The gap probability is the
   quantity that actually depends on which is true.

3. Hotspot azimuth.  FLIGHT at vza = sza must be brighter at raa = 0 than at
   raa = 180.  This pins the azimuth convention against the model rather than
   against a reading of the source.

4. The same, for gortt, so that a convention error in either is attributed to
   the right one.

Nothing here writes to /home/tquaife/src/flight.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc


def mc_cover(scn, crowns, n=400000, seed=1):
    """Our own estimate of the projected crown cover inside the scene box."""
    rng = np.random.default_rng(seed)
    hw = scn.half_width
    px = rng.uniform(-hw, hw, n)
    py = rng.uniform(-hw, hw, n)
    hit = np.zeros(n, dtype=bool)
    r2 = scn.r ** 2
    for cx, cy, _ in crowns:
        hit |= ((px - cx) ** 2 + (py - cy) ** 2) <= r2
    return hit.mean()


def gortt_nadir_gap(scn):
    """
    P(0) at nadir from `gortt -W`: the "p 0 0" row of the probability file,
    which is the ground layer at theta = 0, summed over its p_n0 and epgap
    columns.  See the format comment on gortt_write_prob_file.
    """
    out = sc.gortt_probabilities(scn)
    for line in out.splitlines():
        f = line.split()
        if len(f) == 5 and f[0] == "p" and f[1] == "0" and f[2] == "0":
            return float(f[3]) + float(f[4])
    sys.exit("validate: no \"p 0 0\" row in -W output")


def flight_run(scn, workdir, mode, sza, vza, raa, realisation=0,
               rleaf=None, tleaf=None, rsoil=None, im_size=None, bounce=None):
    sc.make_workdir(workdir)
    crowns = sc.tile_ghosts(scn, sc.crown_positions(scn, realisation))
    sc.write_crowns(os.path.join(workdir, "crowns.data"), scn, crowns)
    sc.write_spec(workdir, config.BANDS, rleaf=rleaf, tleaf=tleaf, rsoil=rsoil)
    sc.write_reverse(os.path.join(workdir, "reverse.data"),
                     bounce or config.BRDF_BOUNCE, config.BRDF_BRANCH,
                     im_size or config.BRDF_IM_SIZE)
    sc.write_in_flight(os.path.join(workdir, "in_flight.data"), scn, mode,
                       sza, vza, raa, len(config.BANDS))
    out = sc.run_flight(workdir)
    brf = sc.read_spectrum(os.path.join(workdir, "RESULTS", "out-spectrum.csv"))
    return crowns, sc.parse_flight_stdout(out), brf


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenes", nargs="*", default=["sparse-low-sphere",
                                                    "dense-high-sphere",
                                                    "dense-high-prolate"],
                    help="scenes to check (default: %(default)s)")
    ap.add_argument("--im-size", type=int, default=200,
                    help="FLIGHT image side for the gap test; noise goes as "
                         "1/im_size (default: %(default)s)")
    ap.add_argument("--workdir", default=os.path.join(config.WORK_DIR, "validate"))
    args = ap.parse_args()

    if not os.path.exists(config.GORTT_BIN):
        sys.exit("validate: %s not found -- run make first" % config.GORTT_BIN)
    if not os.path.exists(config.flight_binary()):
        sys.exit("validate: no FLIGHT binary -- run build_flight.py first")

    failures = []

    print("%-22s %10s %10s %10s   %10s %10s %8s"
          % ("scene", "cover_MC", "cover_FLT", "cover_GORT",
             "gap_FLIGHT", "gap_gortt", "rel.diff"))
    print("-" * 92)

    for name in args.scenes:
        scn = config.scene_by_name(name)
        wd = os.path.join(args.workdir, name)

        # --- 1 & 2: black leaves, white soil, nadir view and nadir sun -------
        crowns, diag, brf = flight_run(scn, wd, "r", 0.0, 0.0, 0.0,
                                       rleaf=0.0, tleaf=0.0, rsoil=1.0,
                                       im_size=args.im_size)
        cover_mc = mc_cover(scn, crowns)
        gap_flight = brf[0]                 # identical in every band by construction
        gap_gortt = gortt_nadir_gap(scn)
        rel = (gap_flight - gap_gortt) / gap_gortt

        print("%-22s %10.4f %10.4f %10.4f   %10.4f %10.4f %7.1f%%"
              % (name, cover_mc, diag["frac_cover"], scn.nadir_cover,
                 gap_flight, gap_gortt, 100.0 * rel))

        if abs(cover_mc - scn.nadir_cover) > 0.01:
            failures.append("%s: cover %.4f vs analytic %.4f"
                            % (name, cover_mc, scn.nadir_cover))
        if abs(rel) > 0.05:
            failures.append("%s: nadir gap probability differs by %.1f%% "
                            "(FLIGHT %.4f, gortt %.4f)"
                            % (name, 100.0 * rel, gap_flight, gap_gortt))
        if max(brf) - min(brf) > 1e-6:
            failures.append("%s: black-leaf BRF varies between bands (%s)" % (name, brf))

    # --- 3 & 4: the hotspot is at relative azimuth zero ----------------------
    print()
    scn = config.scene_by_name(args.scenes[0])
    sza = 30.0
    wd = os.path.join(args.workdir, "hotspot")
    _, _, back = flight_run(scn, wd, "r", sza, sza, 0.0, im_size=args.im_size)
    _, _, fwd = flight_run(scn, wd, "r", sza, sza, 180.0, im_size=args.im_size)
    g_back = sc.gortt_brdf(scn, config.BRDF_BANDS[1], [(sza, sza, 0.0)])[0]
    g_fwd = sc.gortt_brdf(scn, config.BRDF_BANDS[1], [(sza, sza, 180.0)])[0]

    print("hotspot check at sza = vza = %g, %s band" % (sza, config.BRDF_BANDS[1].name))
    print("  FLIGHT  raa=0 %.4f   raa=180 %.4f   %s"
          % (back[1], fwd[1], "OK" if back[1] > fwd[1] else "WRONG WAY ROUND"))
    print("  gortt   raa=0 %.4f   raa=180 %.4f   %s"
          % (g_back, g_fwd, "OK" if g_back > g_fwd else "WRONG WAY ROUND"))
    if back[1] <= fwd[1]:
        failures.append("FLIGHT hotspot is not at raa = 0")
    if g_back <= g_fwd:
        failures.append("gortt hotspot is not at raa = 0")

    print()
    if failures:
        for f in failures:
            print("FAIL: %s" % f)
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
