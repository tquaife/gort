#!/usr/bin/env python3
"""
Reproduce Figure 5 of Ni et al. (1999) using gortt.

Figure 5 shows the modelled areal proportions of the four scene components --
sunlit crown (Kc), sunlit background (Kg), shaded crown (Kt) and shaded
background (Kz) -- in the principal plane, for the BOREAS Old Black Spruce
(SOBS) site, at solar zenith angles of 30 and 60 degrees.

Those four proportions are exactly what "gortt -prnprop" prints, so the whole
figure is one gortt run plus a plot.  The tree parameters are Table I of the
paper:

    R = 0.76 m,  b = 2.7 m,  lambda = 0.41 /m2,
    h1 = 3.0 m,  h2 = 8.5 m,  FAVD = 0.858 /m

which are (near enough) gortt's own built-in defaults -- the defaults in
gortt_set_defaults are the SOBS site.  They are passed explicitly here anyway,
so that the figure does not silently change if the defaults ever do.

Sign convention on the x axis: the paper puts the hotspot at *positive* viewing
zenith, so positive vza is the backscatter direction (relative azimuth 0) and
negative vza is forward scatter (relative azimuth 180).

  usage:  ./fig5_ni1999.py [-o out.png] [--style paper|colour] [--csv data.csv]

Run it from anywhere; it finds ../gortt relative to this file unless --gortt
says otherwise.

The reproduction is close but not exact: Kt, Kz and Kg match the published
curves well, and Kc matches at the hotspot and at nadir, but Kc is too low away
from the hotspot in both directions.  See README.md in this directory -- that
discrepancy looks like the F-versus-F-squared question recorded in
OPEN-QUESTIONS.md, and this figure is a fairly direct test of it.
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

# ---------------------------------------------------------------------------
# Table I of Ni et al. (1999): input tree parameters for the SOBS site.
# ---------------------------------------------------------------------------
SOBS = [
    "-r", "0.76",       # horizontal crown radius (m)
    "-b", "2.7",        # vertical crown radius (m)
    "-lambda", "0.41",  # crown count density (1/m2)
    "-h1", "3.0",       # lower boundary of crown centres (m)
    "-h2", "8.5",       # upper boundary of crown centres (m)
    "-favd", "0.858",   # foliage area volume density (1/m)
]

# Solar zenith angles of the two panels, top to bottom.
SZA_PANELS = (30.0, 60.0)

# The wavelength is irrelevant -- only the viewed proportions are plotted --
# but gortt requires at least one, and it must be within 400-2500 nm.
WAVELENGTH = 800.0

# Component order as gortt prints it inside [ ].  Line styles are the paper's.
COMPONENTS = [
    ("Kc", "sunlit crown",      "-"),
    ("Kg", "sunlit background", "--"),
    ("Kt", "shaded crown",      "-."),
    ("Kz", "shaded background", ":"),
]

# Categorical slots 1-4 of the reference palette, in fixed order.  Only used by
# --style colour; the default reproduction is the paper's black and white.  The
# dash patterns above are carried in both styles, so the four curves are never
# distinguished by colour alone.
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def build_angles(vza, sza):
    """
    Principal-plane angles for one solar zenith, as (vza, vaa, sza, saa) rows.

    gortt takes a non-negative viewing zenith plus azimuths, so a negative
    plotting angle becomes |vza| with the sun moved to the opposite azimuth.
    """
    rows = []
    for v in vza:
        if v < 0.0:
            rows.append((abs(v), 0.0, sza, 180.0))  # forward scatter
        else:
            rows.append((v, 0.0, sza, 0.0))         # backscatter, hotspot side
    return rows


def run_gortt(binary, rows, extra_args=()):
    """
    Run gortt once over every angle and return the four proportions per row.

    Returns an (nrows, 4) array of Kc, Kg, Kt, Kz.
    """
    header = "%d 1 %g\n" % (len(rows), WAVELENGTH)
    body = "".join("%g %g %g %g\n" % r for r in rows)

    argv = [binary, "-prnprop"] + SOBS + list(extra_args)

    # gortt prints %f -- six decimals -- unless GORTT_PREC asks for more.  At
    # the default the four proportions only sum to one to about 1e-6, which is
    # visible if you difference them.
    env = dict(os.environ, GORTT_PREC="17")

    try:
        proc = subprocess.run(
            argv, input=header + body, capture_output=True, text=True, env=env
        )
    except OSError as err:
        sys.exit("fig5: cannot run %s: %s" % (binary, err))

    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("fig5: %s exited %d" % (binary, proc.returncode))

    out = []
    for line in proc.stdout.splitlines()[1:]:  # first line echoes the header
        if "[" not in line:
            continue
        inside = line.split("[", 1)[1].split("]", 1)[0]
        vals = [float(x) for x in inside.split()]
        if len(vals) != 4:
            sys.exit("fig5: expected 4 proportions, got %d: %s" % (len(vals), line))
        out.append(vals)

    if len(out) != len(rows):
        sys.exit("fig5: asked for %d angles, parsed %d" % (len(rows), len(out)))

    return np.array(out)


def plot(vza, panels, style, title):
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.0))

    for ax, (sza, k) in zip(axes, panels):
        for i, (name, longname, dash) in enumerate(COMPONENTS):
            if style == "colour":
                colour, width = COLOURS[i], 1.8
            else:
                colour, width = "black", 1.1
            ax.plot(
                vza, k[:, i],
                dash, color=colour, linewidth=width,
                label="$K_%s$  %s" % (name[1], longname),
            )

        ax.set_xlim(-90, 90)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(-80, 81, 20))
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.set_xlabel("viewing zenith angle")
        ax.set_ylabel("K")

        # Panel label, as in the paper.  szn=60 sits lower to stay clear of the
        # shaded-crown curve, which runs close to 1 on the left of that panel.
        ax.text(
            0.04, 0.88 if sza == 30 else 0.72, "szn=%g" % sza,
            transform=ax.transAxes, fontsize=11,
        )

        # Recessive axes: keep the frame but drop the emphasis.
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.tick_params(width=0.8)

    # One legend for both panels; the paper carries this in its caption.
    axes[0].legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.20),
        ncol=2, frameon=False, fontsize=9,
    )

    if title:
        fig.suptitle(title, y=0.99, fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def main():
    here = os.path.dirname(os.path.abspath(__file__))

    ap = argparse.ArgumentParser(
        description="Reproduce Figure 5 of Ni et al. (1999) using gortt."
    )
    ap.add_argument(
        "-o", "--output", default=os.path.join(here, "fig5_ni1999.png"),
        help="output image (default: %(default)s); extension picks the format",
    )
    ap.add_argument(
        "--gortt", default=os.path.join(here, os.pardir, "gortt"),
        help="path to the gortt binary (default: ../gortt)",
    )
    ap.add_argument(
        "--style", choices=("paper", "colour"), default="paper",
        help="paper: black lines as published; colour: same dashes, coloured",
    )
    ap.add_argument(
        "--step", type=float, default=0.5,
        help="viewing zenith step in degrees (default: %(default)s)",
    )
    ap.add_argument(
        "--limit", type=float, default=89.0,
        help="largest |viewing zenith| in degrees (default: %(default)s). "
             "Both the direct and -P paths return -nan at exactly 90, so the "
             "endpoint is excluded by default",
    )
    ap.add_argument(
        "--csv", metavar="FILE",
        help="also write the plotted numbers to FILE",
    )
    ap.add_argument(
        "--title", default="Ni et al. (1999) Fig. 5, reproduced with gortt "
                           "-- SOBS site, principal plane",
        help="figure title ('' for none)",
    )
    ap.add_argument(
        "--q08", action="store_true",
        help="use the Quaife et al. (2008) gap probability approximation "
             "(-q08_pn_kopen); much faster, only strictly valid for h=0",
    )
    args = ap.parse_args()

    if not os.path.exists(args.gortt):
        sys.exit("fig5: %s not found -- run make first" % args.gortt)

    vza = np.arange(-args.limit, args.limit + args.step / 2.0, args.step)
    extra = ["-q08_pn_kopen"] if args.q08 else []

    panels = []
    for sza in SZA_PANELS:
        sys.stderr.write("fig5: running gortt for szn=%g (%d angles)...\n"
                         % (sza, len(vza)))
        panels.append((sza, run_gortt(args.gortt, build_angles(vza, sza), extra)))

    fig = plot(vza, panels, args.style, args.title)
    fig.savefig(args.output, dpi=150)
    sys.stderr.write("fig5: wrote %s\n" % args.output)

    if args.csv:
        with open(args.csv, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["sza", "vza", "Kc", "Kg", "Kt", "Kz"])
            for sza, k in panels:
                for v, row in zip(vza, k):
                    w.writerow(["%g" % sza, "%g" % v] + ["%.17g" % x for x in row])
        sys.stderr.write("fig5: wrote %s\n" % args.csv)


if __name__ == "__main__":
    main()
