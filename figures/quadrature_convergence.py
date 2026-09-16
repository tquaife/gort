#!/usr/bin/env python3
"""
Show that the gap probabilities no longer depend on -nlayers, and that the
quadratures behind them are converged at the default -nquad.

Four quadratures used to take their step size from nlayers, which is the
model's *output* grid: how finely the lidar profile is reported, and the h
index of the probability tables.  Accuracy therefore moved when you asked for
a different resolution, and the BRF -- which reads only the ground layer --
was silently inheriting whatever error the default 15 layers happened to give.
They now each have their own resolution, scaled by -nquad.

This script draws the evidence:

  (a) k_openep against -nlayers.  Flat is the point.  With --ref it also
      plots an older binary, whose curve drifts as O(1/n) towards the same
      value the current code returns at every resolution.
  (b) k_openep against -nquad, i.e. is the default actually converged.
  (c) p_n0[0][t] against its closed form.  At the ground layer Vgamma has an
      exact solution at every zenith angle -- see below -- so this panel is a
      check against arithmetic rather than against a finer grid.
  (d) how ragged the lidar dP/dz profile is, against -nlayers.

The closed form in (c).  GORT puts the lowest foliage at z1 = h1 - b, and in
the sphere-transformed space b becomes r, so z1_p = h1_p - r: the ground layer
sits exactly one crown radius below the lowest crown centre.  Every crown is
then entirely above it, the crown projection cross section is the constant
pi r^2 sec(theta'), and

    Vgamma(theta', z1_p) = pi r^2 sec(theta') ( h2_p - h1_p )
    p_n0[0][t]           = exp( -lv' Vgamma )

with no quadrature at all.  gortt matches this to about 1 ulp.

Usage
-----
    ./quadrature_convergence.py
    ./quadrature_convergence.py --ref /path/to/older/gortt
    ./quadrature_convergence.py --csv out.csv --outdir .
"""

import argparse
import math
import os
import subprocess
import sys

import numpy as np

import gortt_run


# Canopies to sweep.  The default is gortt's own (the SOBS stand); the other
# two bracket it, because a quadrature that is converged for one crown density
# is not automatically converged for another.
CANOPIES = {
    "default": [],
    "dense": ["-LAI", "6.0", "-lambda", "0.6"],
    "sparse": ["-LAI", "1.5", "-lambda", "0.05"],
}

# The geometry the closed form in panel (c) is written for.  These are
# gortt_set_defaults' values; passing them explicitly keeps the analytic
# comparison honest if a default ever moves.
#
# n.b. b must be passed to full precision.  gortt_set_defaults writes
# b = 3.55263 * r, which is 2.699998..., and rounding that to 2.700 shifts
# theta' by about a part in 10^6 -- invisible at nadir, but p_n0 is 1e-65 at
# 89 degrees and an exponent is not a forgiving place to lose six digits.
GEOM_R, GEOM_H1, GEOM_H2, GEOM_LAM = 0.76, 3.0, 8.5, 0.405
GEOM_B = 3.55263 * GEOM_R
GEOM = ["-r", "%.17g" % GEOM_R, "-b", "%.17g" % GEOM_B,
        "-h1", "%.17g" % GEOM_H1, "-h2", "%.17g" % GEOM_H2,
        "-lambda", "%.17g" % GEOM_LAM, "-favd", "0.858"]


def _run(binary, args, stdin=""):
    env = dict(os.environ, GORTT_PREC="17")
    try:
        proc = subprocess.run([binary] + list(args), input=stdin,
                              capture_output=True, text=True, env=env)
    except OSError as err:
        sys.exit("quadrature: cannot run %s: %s" % (binary, err))
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("quadrature: %s exited %d" % (binary, proc.returncode))
    return proc.stdout


def probabilities(binary, args):
    """
    -W output: p_n0[0][t] and epgap[0][t] for t = 0..90, plus k_open[0] and
    k_openep[0] off the final "-1" row.
    """
    pn0, epg, kopen, kopenep = [], [], None, None
    for line in _run(binary, ["-W"] + list(args)).splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "-1":
            kopen, kopenep = float(f[1]), float(f[2])
        else:
            pn0.append(float(f[1]))
            epg.append(float(f[2]))
    if kopen is None:
        sys.exit("quadrature: no k_open row in -W output")
    return np.array(pn0), np.array(epg), kopen, kopenep


def profile(binary, args):
    """-lidar profile: the z, P(z) and dP/dz columns, top-down."""
    rows = [l.split() for l in _run(binary, ["-lidar"] + list(args)).splitlines()]
    data = np.array([[float(x) for x in r[:3]] for r in rows[1:]])
    return data[:, 0], data[:, 1], data[:, 2]


def raggedness(dpdz):
    """
    rms second difference of dP/dz as a percentage of its peak.

    A second difference has zero gain at zero frequency and maximum gain at
    the Nyquist of the layer grid, so it measures exactly the layer-to-layer
    wobble and ignores the shape of the profile underneath it.
    """
    d2 = dpdz[2:] - 2.0 * dpdz[1:-1] + dpdz[:-2]
    return 100.0 * math.sqrt(float(np.mean(d2 ** 2))) / float(np.max(dpdz))


def pn0_ground_exact(theta_deg):
    """p_n0[0][t] in closed form -- see the module docstring."""
    ell = GEOM_B / GEOM_R
    h1p, h2p = GEOM_H1 / ell, GEOM_H2 / ell
    lvp = (GEOM_LAM / (GEOM_H2 - GEOM_H1)) * ell

    out = []
    for td in theta_deg:
        # gortt clamps 90 degrees to 89; sec(90) is infinite and the
        # probability tables have to stay finite.
        th = math.radians(min(td, 89.0))
        # gortt clamps theta' the same way, but atan() never reaches 90
        # degrees so that clamp is dead code and is not reproduced here.
        thp = math.atan(math.tan(th) * ell)
        vg = math.pi * GEOM_R ** 2 / math.cos(thp) * (h2p - h1p)
        out.append(math.exp(-lvp * vg))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gortt", default=gortt_run.default_binary())
    ap.add_argument("--ref", default=None,
                    help="second binary to overlay, e.g. a build of an older commit")
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--csv", default=None)
    ap.add_argument("--nlayers", type=int, nargs="+",
                    default=[15, 20, 30, 45, 60, 90, 120, 180, 240])
    ap.add_argument("--nquad", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    args = ap.parse_args()

    gortt_run.check_binary(args.gortt)
    if args.ref:
        gortt_run.check_binary(args.ref)

    rows = []

    # ---- (a) against -nlayers -------------------------------------------
    print("k_openep against -nlayers (flat is the point):")
    sweep_n = {}
    for name, opts in CANOPIES.items():
        vals = []
        for n in args.nlayers:
            _, _, _, ke = probabilities(args.gortt, opts + ["-nlayers", str(n)])
            vals.append(ke)
            rows.append(("nlayers", name, n, ke))
        sweep_n[name] = np.array(vals)
        spread = (vals[-1] - vals[0]) / vals[0] if vals[0] else 0.0
        print("  %-8s %.12g  spread over nlayers %d..%d: %.2e relative"
              % (name, vals[0], args.nlayers[0], args.nlayers[-1], abs(spread)))

    ref_n = {}
    if args.ref:
        print("reference binary, same sweep:")
        for name, opts in CANOPIES.items():
            vals = [probabilities(args.ref, opts + ["-nlayers", str(n)])[3]
                    for n in args.nlayers]
            ref_n[name] = np.array(vals)
            print("  %-8s %.12g -> %.12g" % (name, vals[0], vals[-1]))

    # ---- (b) against -nquad ---------------------------------------------
    print("\nk_openep against -nquad (is the default converged):")
    sweep_q = {}
    for name, opts in CANOPIES.items():
        vals = []
        for q in args.nquad:
            _, _, _, ke = probabilities(args.gortt, opts + ["-nquad", str(q)])
            vals.append(ke)
            rows.append(("nquad", name, q, ke))
        sweep_q[name] = np.array(vals)
        err = abs(vals[0] - vals[-1]) / abs(vals[-1])
        print("  %-8s default is %.2e relative from -nquad %d"
              % (name, err, args.nquad[-1]))

    # ---- (c) the ground layer against its closed form --------------------
    theta = np.arange(91.0)
    pn0, _, _, _ = probabilities(args.gortt, GEOM)
    exact = pn0_ground_exact(theta)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(pn0 - exact) / np.abs(exact)
    print("\np_n0[0][t] against its closed form:")
    # Past about 80 degrees p_n0 is below 1e-20 and a relative comparison is
    # really a comparison of exponents, so quote the bulk separately.
    print("  t <=  80: max relative difference %.2e  (at t = %d)"
          % (np.max(rel[:81]), int(np.argmax(rel[:81]))))
    print("  t <=  90: max relative difference %.2e  (at t = %d), where p_n0 is %.1e"
          % (np.max(rel), int(np.argmax(rel)), exact[int(np.argmax(rel))]))

    # ---- (d) raggedness --------------------------------------------------
    print("\nlidar dP/dz raggedness, rms 2nd difference as % of peak:")
    rag_n = [20, 40, 80, 160]
    rag = [raggedness(profile(args.gortt, ["-nlayers", str(n)])[2]) for n in rag_n]
    rag_ref = None
    if args.ref:
        rag_ref = [raggedness(profile(args.ref, ["-nlayers", str(n)])[2]) for n in rag_n]
    for i, n in enumerate(rag_n):
        line = "  nlayers=%-4d  %6.2f%%" % (n, rag[i])
        if rag_ref:
            line += "   (reference %6.2f%%)" % rag_ref[i]
        print(line)
        rows.append(("raggedness", "default", n, rag[i]))

    if args.csv:
        with open(args.csv, "w") as fp:
            fp.write("sweep,canopy,x,value\n")
            for r in rows:
                fp.write("%s,%s,%g,%.17g\n" % r)
        print("\nwrote %s" % args.csv)

    # ---- plot ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not available; numbers above are the result)")
        return

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    for name in CANOPIES:
        ax[0][0].plot(args.nlayers, sweep_n[name] / sweep_n[name][0],
                      "o-", label=name)
        if name in ref_n:
            ax[0][0].plot(args.nlayers, ref_n[name] / ref_n[name][0],
                          "x--", alpha=0.5, label=name + " (ref)")
    ax[0][0].set_xlabel("-nlayers")
    ax[0][0].set_ylabel("k_openep, normalised to its first value")
    ax[0][0].set_title("(a) accuracy must not depend on the output grid")
    ax[0][0].set_xscale("log")
    ax[0][0].legend(fontsize=7)
    ax[0][0].grid(alpha=0.3)

    # The finest run is the yardstick, so it is not plotted against itself --
    # that point is zero by construction and would read as machine precision.
    for name in CANOPIES:
        final = sweep_q[name][-1]
        err = np.abs(sweep_q[name][:-1] - final) / abs(final)
        ax[0][1].plot(args.nquad[:-1], err, "o-", label=name)
    ax[0][1].set_xlabel("-nquad")
    ax[0][1].set_ylabel("relative difference from -nquad %d" % args.nquad[-1])
    ax[0][1].set_title("(b) the default orders are converged")
    ax[0][1].set_xscale("log")
    ax[0][1].set_yscale("log")
    ax[0][1].legend(fontsize=7)
    ax[0][1].grid(alpha=0.3)

    ax[1][0].semilogy(theta, np.maximum(rel, 1e-18), "-")
    ax[1][0].axhline(2.3e-16, ls=":", c="k", lw=0.8)
    ax[1][0].text(2, 3e-16, "1 ulp", fontsize=7)
    ax[1][0].set_xlabel("zenith angle (degrees)")
    ax[1][0].set_ylabel("|gortt - closed form| / closed form\n(floored at 1e-18; the floor is an exact match)")
    ax[1][0].set_title("(c) p_n0 at the ground against its exact solution")
    ax[1][0].grid(alpha=0.3)

    ax[1][1].plot(rag_n, rag, "o-", label="now")
    if rag_ref:
        ax[1][1].plot(rag_n, rag_ref, "x--", alpha=0.6, label="reference")
    ax[1][1].set_xlabel("-nlayers")
    ax[1][1].set_ylabel("rms 2nd difference of dP/dz, % of peak")
    ax[1][1].set_title("(d) raggedness of the lidar profile")
    ax[1][1].set_xscale("log")
    ax[1][1].set_yscale("log")
    ax[1][1].legend(fontsize=7)
    ax[1][1].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(args.outdir, "quadrature_convergence.png")
    fig.savefig(out, dpi=140)
    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
