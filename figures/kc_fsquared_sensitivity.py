#!/usr/bin/env python3
"""
How much did dropping the F-squared from gortt_kc change the model?

gortt_kc interpolates the Li & Strahler (1992) f-ratio between the two halves
of the principal plane.  Up to commit e109ff0 (tag pre-kc-fsquared) it read

    f = (1.-frac)*f0*F0 + frac*f180*F180 ;

but gortt_kc_fFbeta already returns f with its factor of F in it -- eq. (17)
of the paper -- so that carried F squared.  Eq. (25) interpolates f itself:

    f = (1.-frac)*f0 + frac*f180 ;

This script runs an OLD and a NEW gortt binary on identical inputs and reports
the difference in Kc, in the BRF, and in the hemispherically integrated
albedo, over a set of canopies and sun angles.

  usage:  ./kc_fsquared_sensitivity.py --old /path/to/gortt_old [--new ../gortt]
                                       [-o out.png] [--csv data.csv]

To get an OLD binary:

    git worktree add /tmp/gort-old pre-kc-fsquared && make -C /tmp/gort-old

It also checks three things that must hold if the change is what it claims to
be, and exits non-zero if any of them fails:

  1. -beta 0 is byte-identical between the binaries.  f is weighted by beta,
     so at beta = 0 ("random" crown heights, f = F) the change cannot show.
  2. Kc is unchanged at the hotspot (vza = sza, raa = 0), where F0 = 1.

It then reports, for OLD and NEW, how many geometries hit the Kt clamp in
gortt_rsurf (Kt = MAX(0, 1 - Kc - Kg - Kz)), i.e. have Kc > 1 - Kg - Kz, and
the smallest forward-scatter principal-plane vza at which that starts.  This
is reported rather than failed on: Kc can only grow under the change (F <= 1),
and it does hit the clamp at large forward-scatter view zeniths.  That comes
from f in eq. (17) growing like sec(vza') with the P terms as gortt evaluates
them, which the extra F used to hide -- see OPEN-QUESTIONS.md.
"""

import argparse
import csv
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gortt_run as gr

HERE = os.path.dirname(os.path.abspath(__file__))

SZAS = (15.0, 30.0, 45.0, 60.0)

# Principal plane, signed as in gortt_run.principal_plane: + is backscatter.
PP_VZA = np.arange(-80.0, 80.0 + 1e-9, 2.0)

# The principal cone: fixed view zenith, relative azimuth swept 0..180.
CONE_VZA = (20.0, 45.0, 70.0)
CONE_RAA = np.arange(0.0, 180.0 + 1e-9, 15.0)

# Canopies: each is appended to gortt_run.SOBS, so the later option wins.
# name, extra args, line style.  Colours are assigned in this fixed order.
CANOPIES = [
    ("SOBS",            [],                                   "-"),
    ("sparse λ=0.1",    ["-lambda", "0.1"],                   "--"),
    ("dense λ=1.0",     ["-lambda", "1.0"],                   "-."),
    ("SOBS, β=1",       ["-beta", "1"],                       ":"),
    ("HB2.5 BR1.5 PCC.6", ["-HB", "2.5", "-BR", "1.5", "-PCC", "0.6"], (0, (5, 1, 1, 1, 1, 1))),
]
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

# Table II of Ni et al. (1999), with the red background at 0.117 -- see the
# fig67_ni1999.py notes in README.md.
BANDS = [
    ("red", 670.0, ["-refl_leaf", "0.09", "-trans_leaf", "0.002", "-alb_soil", "0.117"]),
    ("nir", 865.0, ["-refl_leaf", "0.515", "-trans_leaf", "0.305", "-alb_soil", "0.29"]),
]


def run_raw(binary, rows, wavelength, args):
    """Run gortt at one wavelength with -prnprop; return (stdout, parsed rows)."""
    header = "%d 1 %g\n" % (len(rows), wavelength)
    body = "".join("%.17g %.17g %.17g %.17g\n" % r for r in rows)
    env = dict(os.environ, GORTT_PREC="17")
    argv = [binary, "-prnprop"] + gr.SOBS + list(args)
    proc = subprocess.run(argv, input=header + body, capture_output=True,
                          text=True, env=env)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("sensitivity: %s exited %d" % (binary, proc.returncode))

    out = []
    for line in proc.stdout.splitlines()[1:]:
        pre, rest = line.split("[", 1)
        props, post = rest.split("]", 1)
        pre = pre.split()
        out.append(dict(brf=float(pre[4]),
                        K=[float(x) for x in props.split()],
                        energy=[float(x) for x in post.split()]))
    if len(out) != len(rows):
        sys.exit("sensitivity: asked for %d angles, parsed %d" % (len(rows), len(out)))
    return proc.stdout, out


def cone_rows(sza):
    return [(v, float(r), sza, 0.0) for v in CONE_VZA for r in CONE_RAA]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="gortt built before the change")
    ap.add_argument("--new", default=gr.default_binary(), help="gortt built after it")
    ap.add_argument("-o", "--output", default=os.path.join(HERE, "kc_fsquared_sensitivity.png"))
    ap.add_argument("--csv", help="also write every old/new pair to this file")
    a = ap.parse_args()
    for b in (a.old, a.new):
        gr.check_binary(b)

    failures = []
    records = []      # for --csv
    summary = []      # one row per canopy x sza
    curves = {}       # (canopy, sza) -> (vza, dKc, dBRF_nir)

    # ---- check 1: beta = 0 is byte-identical --------------------------------
    for bname, wl, bargs in BANDS:
        for canopy, extra, _ in CANOPIES:
            if "-beta" in extra:
                continue
            rows = []
            for s in SZAS:
                rows += gr.principal_plane(PP_VZA, s) + cone_rows(s)
            args = extra + bargs + ["-beta", "0"]
            if run_raw(a.old, rows, wl, args)[0] != run_raw(a.new, rows, wl, args)[0]:
                failures.append("beta=0 not byte-identical: %s, %s" % (canopy, bname))

    # ---- the sweep ----------------------------------------------------------
    for canopy, extra, _ in CANOPIES:
        for s in SZAS:
            pp = gr.principal_plane(PP_VZA, s)
            cone = cone_rows(s)
            hot = [(s, 0.0, s, 0.0)]
            rows = pp + cone + hot
            npp, ncone = len(pp), len(cone)

            row = dict(canopy=canopy, sza=s)
            for bname, wl, bargs in BANDS:
                old = run_raw(a.old, rows, wl, extra + bargs)[1]
                new = run_raw(a.new, rows, wl, extra + bargs)[1]

                kc_o = np.array([r["K"][0] for r in old])
                kc_n = np.array([r["K"][0] for r in new])
                brf_o = np.array([r["brf"] for r in old])
                brf_n = np.array([r["brf"] for r in new])

                # check 2: hotspot
                if kc_o[-1] != kc_n[-1]:
                    failures.append("hotspot Kc moved: %s sza=%g %r -> %r"
                                    % (canopy, s, kc_o[-1], kc_n[-1]))
                # the Kt clamp: Kc + Kg + Kz > 1
                for tag, res in (("old", old), ("new", new)):
                    K = np.array([r["K"] for r in res])
                    over = K[:, 0] + K[:, 1] + K[:, 3] > 1.0 + 1e-12
                    row["clamp_" + tag] = int(over.sum())
                    fwd = [abs(v) for v, o in zip(PP_VZA, over[:npp]) if o and v < 0]
                    row["clamp_vza_" + tag] = min(fwd) if fwd else None

                dkc = kc_n - kc_o
                dbrf = brf_n - brf_o
                row["dKc"] = np.max(np.abs(dkc))
                row["Kc_old_at_max"] = kc_o[np.argmax(np.abs(dkc))]
                row["dBRF_" + bname] = np.max(np.abs(dbrf))
                row["rBRF_" + bname] = np.max(np.abs(dbrf) / np.abs(brf_o))
                if bname == "nir":
                    curves[(canopy, s)] = (PP_VZA, dkc[:npp], dbrf[:npp])

                for (v, vaa, ss, saa), o, n, ko, kn in zip(rows, brf_o, brf_n, kc_o, kc_n):
                    records.append((canopy, bname, ss, v, (vaa - saa) % 360.0, ko, kn, o, n))

                # hemispherical terms: independent of view, so one row at nadir
                eo = run_raw(a.old, [(0.0, 0.0, s, 0.0)], wl, extra + bargs + ["-energy"])[1][0]["energy"]
                en = run_raw(a.new, [(0.0, 0.0, s, 0.0)], wl, extra + bargs + ["-energy"])[1][0]["energy"]
                row["alb_old_" + bname] = eo[0]
                row["dalb_" + bname] = en[0] - eo[0]
                row["dfaveg_" + bname] = en[1] - eo[1]
            summary.append(row)

    # ---- report -------------------------------------------------------------
    print("Max over the principal plane (vza -80..80) and principal cone "
          "(vza 20/45/70, raa 0..180).")
    print("dX = new - old.  Albedo is the white-sky/black-sky mix gortt_energy "
          "returns at that sza.\n")
    hdr = ("%-18s %4s | %7s | %8s %6s | %8s %6s | %6s %8s %8s | %6s %8s %8s"
           % ("canopy", "sza", "max|dKc|", "|dBRF|r", "rel", "|dBRF|n", "rel",
              "alb_r", "dalb_r", "dfveg_r", "alb_n", "dalb_n", "dfveg_n"))
    print(hdr)
    print("-" * len(hdr))
    for r in summary:
        print("%-18s %4.0f | %7.4f | %8.5f %5.1f%% | %8.5f %5.1f%% | %6.4f %+8.5f %+8.5f | %6.4f %+8.5f %+8.5f"
              % (r["canopy"], r["sza"], r["dKc"],
                 r["dBRF_red"], 100 * r["rBRF_red"], r["dBRF_nir"], 100 * r["rBRF_nir"],
                 r["alb_old_red"], r["dalb_red"], r["dfaveg_red"],
                 r["alb_old_nir"], r["dalb_nir"], r["dfaveg_nir"]))

    print("\nGeometries hitting the Kt clamp (Kc > 1 - Kg - Kz), of %d per canopy x sza;"
          % (len(PP_VZA) + len(CONE_VZA) * len(CONE_RAA) + 1))
    print("'from' is the smallest forward-scatter principal-plane vza where it starts.\n")
    print("%-18s %4s | %5s %5s | %5s %5s" % ("canopy", "sza", "old", "from", "new", "from"))
    fmt_v = lambda v: "-" if v is None else "%.0f" % v
    for r in summary:
        print("%-18s %4.0f | %5d %5s | %5d %5s"
              % (r["canopy"], r["sza"], r["clamp_old"], fmt_v(r["clamp_vza_old"]),
                 r["clamp_new"], fmt_v(r["clamp_vza_new"])))

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["canopy", "band", "sza", "vza", "raa",
                        "Kc_old", "Kc_new", "BRF_old", "BRF_new"])
            w.writerows(records)
        print("\nwrote", a.csv)

    plot(curves, a.output)
    print("wrote", a.output)

    print()
    if failures:
        for f in failures:
            print("CHECK FAILED:", f)
        sys.exit(1)
    print("checks passed: beta=0 byte-identical, hotspot Kc unchanged")


def plot(curves, path):
    ink, muted, grid = "#1f1f1e", "#6b6a64", "#e4e3dc"
    fig, axes = plt.subplots(2, len(SZAS), figsize=(13, 6.4), sharex=True, sharey="row")
    for j, s in enumerate(SZAS):
        for i, (label, k) in enumerate((("ΔKc", 1), ("ΔBRF, 865 nm", 2))):
            ax = axes[i, j]
            ax.axhline(0.0, color=muted, lw=0.8)
            ax.axvline(s, color=grid, lw=6, zorder=0)   # hotspot band
            for c, (canopy, _, ls) in enumerate(CANOPIES):
                v = curves[(canopy, s)]
                ax.plot(v[0], v[k], color=COLOURS[c], ls=ls, lw=2, label=canopy)
            ax.grid(True, color=grid, lw=0.6)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(muted)
            ax.tick_params(colors=muted, labelsize=8)
            if i == 0:
                ax.set_title("sza = %g°" % s, color=ink, fontsize=10)
            if j == 0:
                ax.set_ylabel(label + "  (new − old)", color=ink, fontsize=9)
            if i == 1:
                ax.set_xlabel("signed vza (°), + = backscatter", color=muted, fontsize=8)
    axes[0, 0].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.suptitle("gortt_kc: effect of dropping the extra F (principal plane; grey band = hotspot)",
                 color=ink, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


if __name__ == "__main__":
    main()
