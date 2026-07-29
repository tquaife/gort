#!/usr/bin/env python3
"""
Reproduce Figures 6 and 7 of Ni et al. (1999) using gortt.

Figures 6 and 7 are the same pair of panels for two wavebands -- Fig. 6 in the
visible ("RED"), Fig. 7 in the near infrared -- for the SOBS site:

  (a) the four component signatures C, G, T and Z as a function of solar
      zenith angle, with vertical bars showing the range of C over viewing
      directions at that solar zenith;
  (b) C at a solar zenith of 45 degrees as a function of viewing zenith angle
      in the principal plane.

These are the "{ C G T Z }" values that "gortt -prnspec" prints.

Bypassing PROSPECT and the Price soil
-------------------------------------
The figures are drawn for the fixed optical properties of the paper's Table II,
not for a PROSPECT leaf or a Price soil, so both are switched off from the
command line:

    -refl_leaf R   leaf reflectance   } together these turn PROSPECT off
    -trans_leaf T  leaf transmittance }
    -alb_soil A    background albedo, which turns the Price EOFs off

Table II:

    band     reflectance  transmittance   background albedo
    540 nm      0.2          0.06              0.069
    670 nm      0.09         0.002             0.29   (see below)
    865 nm      0.515        0.305             0.29

Reflectance and transmittance are given separately because the radiative
transfer needs both: the single-scattering albedo is their sum, but the
scattering asymmetry is built from their difference,

    g_asym = -(4/9)*( rleaf - tleaf )/w

so Table II's asymmetric leaves (g_asym = -0.43 at 670 nm, -0.11 at 865 nm)
can now be given exactly.  The older -alb_leaf could not: it set r = t = W/2
and so always gave g_asym = 0.

The red background albedo used here is 0.117 rather than Table II's 0.29: G is
identically the background albedo, and the published Fig. 6 draws it at 0.117.
See BANDS below and README.md.

  usage:  ./fig67_ni1999.py                 # both figures
          ./fig67_ni1999.py --band red      # just Fig. 6
"""

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gortt_run as gr

# Bands, as the paper's Table II: leaf reflectance and transmittance separately.
#
# n.b. the red background albedo is 0.117, not Table II's 0.29.  G is
# identically the background albedo (gortt.c: Gd = Gf = rsoil, so the fd
# mixture collapses), and the dashed G line in the published Fig. 6 measures
# 0.117 off the page, not 0.29.  Table II gives 0.29 for both 670 and 865 nm,
# which looks like the NIR row duplicated: 0.117 also brings panel (b) closer
# to the published curve at both the nadir and the hotspot.  --soil 0.29
# restores the tabulated value.
BANDS = {
    "red": dict(figure=6, wavelength=670.0, refl=0.09, trans=0.002, soil=0.117,
                label="RED (670 nm)", ymax_a=0.4, ymax_b=0.2),
    "nir": dict(figure=7, wavelength=865.0, refl=0.515, trans=0.305, soil=0.29,
                label="NIR (865 nm)", ymax_a=1.0, ymax_b=0.6),
}

# Component signature order, as gortt_rsurf composes the scene.  Line styles
# are the ones the published legend uses.
COMPONENTS = [
    ("C", "sunlit crown",      "-"),
    ("G", "sunlit background", "--"),
    ("T", "shaded crown",      "-."),
    ("Z", "shaded background", ":"),
]

COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

SZA_PANEL_B = 45.0   # panel (b) is drawn at a solar zenith of 45 degrees


def band_args(band, extra):
    return ["-refl_leaf", str(band["refl"]),
            "-trans_leaf", str(band["trans"]),
            "-alb_soil", str(band["soil"])] + list(extra)


def panel_a_curves(binary, band, szas, extra):
    """
    The four component signature curves against solar zenith angle.

    The published solid C curve rises with solar zenith and peaks near 65-70
    degrees, which is C *at the hotspot* -- the paper says C "has a maximum
    value at the hotspot viewing direction".  C at a fixed nadir view does the
    opposite and falls slightly, so the hotspot is what is plotted here.

    G, T and Z are taken at nadir; the paper notes they are almost constant
    with viewing direction, and G is exactly the background albedo.
    """
    args = band_args(band, extra)
    hotspot = gr.run(binary, [(s, 0.0, s, 0.0) for s in szas],
                     band["wavelength"], args, want="spec")
    nadir = gr.run(binary, [(0.0, 0.0, s, 0.0) for s in szas],
                   band["wavelength"], args, want="spec")
    return hotspot[:, 0], nadir


def panel_a_bars(binary, band, bar_szas, vza_grid, extra):
    """
    Range of C over viewing direction at each of a few solar zeniths, which is
    what the published caption says the vertical bars are.  Drawn on a coarser
    grid than the curves, as in the paper.
    """
    args = band_args(band, extra)
    lo, hi = [], []
    for s in bar_szas:
        c = gr.run(binary, gr.principal_plane(vza_grid, s),
                   band["wavelength"], args, want="spec")[:, 0]
        lo.append(c.min())
        hi.append(c.max())
    return np.array(lo), np.array(hi)


def panel_b_data(binary, band, vza, extra):
    """C against viewing zenith in the principal plane, at 45 degrees solar."""
    sig = gr.run(binary, gr.principal_plane(vza, SZA_PANEL_B),
                 band["wavelength"], band_args(band, extra), want="spec")
    return sig[:, 0]


def plot(band, szas, c_hot, nadir, bar_szas, lo, hi, vza, c_b, style, title):
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.0, 7.5))

    # ---- panel (a) --------------------------------------------------------
    for i, (name, longname, dash) in enumerate(COMPONENTS):
        colour = COLOURS[i] if style == "colour" else "black"
        width = 1.8 if style == "colour" else 1.1
        y = c_hot if i == 0 else nadir[:, i]
        ax.plot(szas, y, dash, color=colour, linewidth=width,
                label="%s  %s" % (name, longname))

    # Range of C over viewing direction, as in the published figure.
    for s, a, b in zip(bar_szas, lo, hi):
        ax.plot([s, s], [a, b],
                "-", color=COLOURS[0] if style == "colour" else "black",
                linewidth=0.6, alpha=0.7)

    ax.set_xlim(0, 80)
    ax.set_ylim(0, band["ymax_a"])
    ax.set_xlabel("solar zenith angle")
    ax.set_ylabel("component signature")
    ax.text(0.03, 0.90, "(a)  %s" % band["label"], transform=ax.transAxes, fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2,
              frameon=False, fontsize=9)

    # ---- panel (b) --------------------------------------------------------
    bx.plot(vza, c_b, "-",
            color=COLOURS[0] if style == "colour" else "black",
            linewidth=1.8 if style == "colour" else 1.1)
    bx.set_xlim(-90, 90)
    bx.set_ylim(0, band["ymax_b"])
    bx.set_xticks(np.arange(-80, 81, 20))
    bx.set_xlabel("viewing zenith angle")
    bx.set_ylabel("component signature")
    bx.text(0.03, 0.88, "(b)  C, szn=%g degree" % SZA_PANEL_B,
            transform=bx.transAxes, fontsize=10)

    for a in (ax, bx):
        for spine in a.spines.values():
            spine.set_linewidth(0.8)
        a.tick_params(width=0.8)

    if title:
        fig.suptitle(title, y=0.99, fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def main():
    here = os.path.dirname(os.path.abspath(__file__))

    ap = argparse.ArgumentParser(
        description="Reproduce Figures 6 and 7 of Ni et al. (1999) using gortt."
    )
    ap.add_argument("--band", choices=("red", "nir", "both"), default="both",
                    help="red is Fig. 6, nir is Fig. 7 (default: both)")
    ap.add_argument("--outdir", default=here,
                    help="where to write the images (default: %(default)s)")
    ap.add_argument("--gortt", default=gr.default_binary(),
                    help="path to the gortt binary (default: ../gortt)")
    ap.add_argument("--style", choices=("paper", "colour"), default="paper",
                    help="paper: black lines as published; colour: same dashes")
    ap.add_argument("--refl", type=float,
                    help="override the leaf reflectance")
    ap.add_argument("--trans", type=float,
                    help="override the leaf transmittance")
    ap.add_argument("--soil", type=float,
                    help="override the background albedo")
    ap.add_argument("--step", type=float, default=1.0,
                    help="viewing zenith step for panel (b) (default: %(default)s)")
    ap.add_argument("--limit", type=float, default=89.0,
                    help="largest |viewing zenith| (default: %(default)s); both "
                         "the direct and -P paths return -nan at exactly 90")
    ap.add_argument("--csv", metavar="PREFIX",
                    help="also write PREFIX-<band>-a.csv and -b.csv")
    ap.add_argument("--q08", action="store_true",
                    help="use -q08_pn_kopen (faster, only strictly valid for h=0)")
    ap.add_argument("--title", default=None,
                    help="figure title ('' for none)")
    args = ap.parse_args()

    gr.check_binary(args.gortt)

    extra = ["-q08_pn_kopen"] if args.q08 else []
    szas = np.arange(0.0, 80.1, 2.0)
    bar_szas = np.arange(0.0, 80.1, 5.0)
    vza = np.arange(-args.limit, args.limit + args.step / 2.0, args.step)

    names = ("red", "nir") if args.band == "both" else (args.band,)

    for name in names:
        band = dict(BANDS[name])
        if args.refl is not None:
            band["refl"] = args.refl
        if args.trans is not None:
            band["trans"] = args.trans
        if args.soil is not None:
            band["soil"] = args.soil

        sys.stderr.write("fig%d: running gortt for %s...\n"
                         % (band["figure"], band["label"]))

        c_hot, nadir = panel_a_curves(args.gortt, band, szas, extra)
        lo, hi = panel_a_bars(args.gortt, band, bar_szas, vza, extra)
        c_b = panel_b_data(args.gortt, band, vza, extra)

        title = args.title
        if title is None:
            title = ("Ni et al. (1999) Fig. %d, reproduced with gortt "
                     "-- SOBS site, %s" % (band["figure"], band["label"]))

        fig = plot(band, szas, c_hot, nadir, bar_szas, lo, hi,
                   vza, c_b, args.style, title)

        out = os.path.join(args.outdir, "fig%d_ni1999.png" % band["figure"])
        fig.savefig(out, dpi=150)
        plt.close(fig)
        sys.stderr.write("fig%d: wrote %s\n" % (band["figure"], out))

        if args.csv:
            with open("%s-%s-a.csv" % (args.csv, name), "w", newline="") as fp:
                w = csv.writer(fp)
                w.writerow(["sza", "C_hotspot", "G", "T", "Z"])
                for s, c, row in zip(szas, c_hot, nadir):
                    w.writerow(["%g" % s, "%.17g" % c] +
                               ["%.17g" % x for x in row[1:]])
            with open("%s-%s-b.csv" % (args.csv, name), "w", newline="") as fp:
                w = csv.writer(fp)
                w.writerow(["vza", "C"])
                for v, c in zip(vza, c_b):
                    w.writerow(["%g" % v, "%.17g" % c])
            sys.stderr.write("fig%d: wrote %s-%s-{a,b}.csv\n"
                             % (band["figure"], args.csv, name))


if __name__ == "__main__":
    main()
