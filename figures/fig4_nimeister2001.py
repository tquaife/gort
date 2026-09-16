#!/usr/bin/env python3
"""
Reproduce Figure 4, case B, of Ni-Meister et al. (2001) using gortt.

The figure is three stacked panels against height, for one canopy:

    P(z)        the total nadir gap probability
    dP/dz       the lidar waveform, up to a constant
    dlogP/dz    the apparent foliage profile

which are exactly the first four columns that "gortt -lidar" prints.  The
published curves are overlaid from fig4_nimeister2001_digitised.csv, traced
off a 400 dpi render of the scanned page; see that file's header and the
"Digitising" section below.

Table I of the paper gives case B as

    R = 1.2 m,  b = 3.6 m,  lambda = 0.24 /m2,  h1 = 6 m,  h2 = 12 m,  Fa = 0.44 /m

-- but the figure was not drawn with h2 = 12.
----------------------------------------------------------------------------
All three panels carry dotted reference lines labelled z1, h1, h2 and z2.
Traced off the page they sit at

    z1 = 2.40    h1 = 6.00    h2 = 10.00    z2 = 13.60   (m)

consistently across the three panels to better than 0.03 m.  GORT defines
z1 = h1 - b and z2 = h2 + b, so z1 = 2.40 with h1 = 6 forces b = 3.60, which is
Table I's value exactly; and z2 = 13.60 with that b then forces h2 = 10.00, not
12.  There is no (h1, b) pair that gives both z1 = 2.4 and z2 = 13.6 with
h2 = 12, and h2 = 12 would put z2 at 15.6, off the top of the figure's own
0-15 m axis.

So the figure and Table I disagree, and this script follows the figure.
--h2 12 runs the tabulated value instead, and --h2 both overlays the two, which
is the quickest way to see that it is the figure's geometry gortt reproduces.

Digitising
----------
The curves were traced by locating the axis spines in the render, then walking
one image row at a time.  For dP/dz and dlogP/dz there is a single curve per
panel, so the longest run of dark pixels in the row is it.  The P panel carries
three curves -- P solid, P(n=0) dash-dot and P(n>0) dashed -- and the solid one
is the rightmost, except above z = 11.5 where it merges with the dash-dot; the
csv holds NaN there and a second, independent P column obtained by integrating
the dP/dz panel down from P(z2) = 1.  The two agree to 0.01-0.03 wherever both
exist, which is the honest digitisation error, and dlogP/dz reproduces
dPdz/P to 0.003 over the whole profile, which says the three panels were traced
consistently.

The wobble in the derivative panels
-----------------------------------
gortt's dP/dz and dlogP/dz are visibly ragged where the published curves are
smooth.  That is not the differencing: it is already in the gap probability
tables, and differencing only makes it visible.  Measured as the rms second
difference at 80 layers, the between-crown term p_n0 is 0.05% of its peak and
the within-crown term epgap 0.21% of its.  Both come from quadratures that step
with the layer thickness -- gortt_crown_proj_volume integrates over crown-centre
heights in steps of dz_p, and gortt_s_to_index bins path lengths onto ds = dz --
so both shrink with --nlayers, but only about as fast as sqrt(dz): the rms
wobble in dP/dz is 16% of its peak at 20 layers, 9% at 40, 6% at 80 and 4% at
120.  Raising --nlayers is expensive, and the physical answer for a waveform is
-lidar_pulse, which convolves it away.
"""

import argparse
import csv
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gortt_run as gr

# Table I, case B.  h2 is the one the figure and the table disagree about, so
# it is not in here -- it comes from the command line.
CASE_B = ["-r", "1.2", "-b", "3.6", "-lambda", "0.24", "-h1", "6", "-favd", "0.44"]

H2_FIGURE = 10.0    # measured off the dotted reference lines
H2_TABLE = 12.0     # what Table I lists

# (column in gortt's -lidar output, axis label, published csv column, x limit)
PANELS = [
    (1, "Gap probability (P)", "P_solid", 1.0),
    (2, "dP/dz", "dPdz", 0.3),
    (3, "dlogP/dz", "dlogPdz", 0.3),
]


def run_lidar(binary, h2, nlayers, extra=()):
    """gortt -lidar for case B; returns (nlayers, 5) of z, P, dP/dz, Fapp, Rv."""
    argv = ([binary, "-lidar", "-h2", "%g" % h2, "-nlayers", "%d" % nlayers]
            + CASE_B + list(extra))
    env = dict(os.environ, GORTT_PREC="17")

    try:
        proc = subprocess.run(argv, input="", capture_output=True, text=True, env=env)
    except OSError as err:
        sys.exit("fig4: cannot run %s: %s" % (binary, err))

    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("fig4: %s exited %d" % (binary, proc.returncode))

    # first line is the header: nlayers, wavelength, rho_v, rho_g, J0, Rg
    lines = proc.stdout.splitlines()
    rows = [[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()]
    if len(rows) != nlayers:
        sys.exit("fig4: asked for %d layers, parsed %d" % (nlayers, len(rows)))
    return np.array(rows)


def read_published(path):
    """The digitised csv; '#' comment lines and empty fields (NaN) allowed."""
    if not os.path.exists(path):
        sys.exit("fig4: %s not found" % path)

    cols, data = None, []
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if cols is None:
                cols = line.split(",")
                continue
            data.append([float(x) if x else np.nan for x in line.split(",")])

    arr = np.array(data)
    return {name: arr[:, i] for i, name in enumerate(cols)}


def panel(ax, pub, models, key, xlabel, xmax, h2_lines, style):
    published = "0.55" if style == "paper" else "tab:red"

    for z_ref, label in h2_lines:
        ax.axhline(z_ref, color="0.75", lw=0.6, ls=":")
        ax.text(xmax * 0.985, z_ref, label, fontsize=7, color="0.45",
                ha="right", va="bottom")

    ax.plot(pub[key], pub["z"], color=published, lw=3.0, solid_capstyle="round",
            label="Ni-Meister et al. (2001), fig. 4B", zorder=1)

    for (h2, arr), (colour, dash) in zip(models, CURVE_STYLES):
        ax.plot(arr[:, PANELS_IDX[key]], arr[:, 0], color=colour, lw=1.4, ls=dash,
                label="gortt, h2 = %g m" % h2, zorder=2)

    ax.set_xlim(0.0, xmax)
    ax.set_ylim(0.0, 16.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("height (m)")
    ax.tick_params(direction="in", top=True, right=True)


PANELS_IDX = {"P_solid": 1, "dPdz": 2, "dlogPdz": 3}
CURVE_STYLES = [("k", "-"), ("tab:blue", "--")]


def main():
    here = os.path.dirname(os.path.abspath(__file__))

    ap = argparse.ArgumentParser(
        description="Reproduce Fig. 4 (case B) of Ni-Meister et al. (2001) using gortt."
    )
    ap.add_argument("--h2", choices=("10", "12", "both"), default="10",
                    help="upper bound of crown centres (m): 10 is what the figure "
                         "was drawn with, 12 is what Table I lists (default: %(default)s)")
    ap.add_argument("--nlayers", type=int, default=80,
                    help="canopy height layers (default: %(default)s). The gap "
                         "probabilities are the expensive part; 80 takes a few seconds.")
    ap.add_argument("--outdir", default=here,
                    help="where to write the image (default: %(default)s)")
    ap.add_argument("--gortt", default=gr.default_binary(),
                    help="path to the gortt binary (default: ../gortt)")
    ap.add_argument("--style", choices=("paper", "colour"), default="paper",
                    help="paper: greys as published; colour: coloured model curves")
    ap.add_argument("--csv", help="also write the model numbers to this file")
    args = ap.parse_args()

    gr.check_binary(args.gortt)

    pub = read_published(os.path.join(here, "fig4_nimeister2001_digitised.csv"))

    h2s = [H2_FIGURE] if args.h2 == "10" else \
          [H2_TABLE] if args.h2 == "12" else [H2_FIGURE, H2_TABLE]

    models = []
    for h2 in h2s:
        sys.stderr.write("fig4: running gortt -lidar, h2 = %g m, %d layers ...\n"
                         % (h2, args.nlayers))
        models.append((h2, run_lidar(args.gortt, h2, args.nlayers)))

    # the paper's own reference lines, for the geometry the figure used
    b, h1 = 3.6, 6.0
    h2_lines = [(h1 - b, "z1"), (h1, "h1"), (H2_FIGURE, "h2"), (H2_FIGURE + b, "z2")]

    if args.style == "colour":
        styles = CURVE_STYLES
    else:
        styles = [("k", "-"), ("k", "--")]
    CURVE_STYLES[:] = styles

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.2))
    for ax, (_, xlabel, key, xmax) in zip(axes, PANELS):
        panel(ax, pub, models, key, xlabel, xmax, h2_lines, args.style)

    axes[0].legend(loc="upper left", fontsize=7, frameon=False)
    fig.suptitle("Ni-Meister et al. (2001) fig. 4, case B: "
                 "R=1.2  b=3.6  $\\lambda$=0.24  h1=6  Fa=0.44", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    out = os.path.join(args.outdir, "fig4_nimeister2001.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    sys.stderr.write("fig4: wrote %s\n" % out)

    # how well does it do?  Compare on the published grid, over the canopy.
    for h2, arr in models:
        m = np.isfinite(pub["P_solid"])
        gp = np.interp(pub["z"][m], arr[::-1, 0], arr[::-1, 1])
        gd = np.interp(pub["z"], arr[::-1, 0], arr[::-1, 2])
        gl = np.interp(pub["z"], arr[::-1, 0], arr[::-1, 3])
        def stat(a, b):
            d = np.abs(a - b)
            return np.sqrt(np.nanmean(d ** 2)), np.nanmax(d)

        sys.stderr.write(
            "fig4: h2 = %5.1f m vs published:  P rms %.3f max %.3f | "
            "dP/dz rms %.4f max %.4f | dlogP/dz rms %.4f max %.4f\n"
            % ((h2,) + stat(gp, pub["P_solid"][m])
                     + stat(gd, pub["dPdz"])
                     + stat(gl, pub["dlogPdz"])))

    if args.csv:
        with open(args.csv, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["h2", "z", "P", "dPdz", "dlogPdz", "Rv"])
            for h2, arr in models:
                for row in arr:
                    w.writerow(["%g" % h2] + ["%.17g" % x for x in row])
        sys.stderr.write("fig4: wrote %s\n" % args.csv)


if __name__ == "__main__":
    main()
