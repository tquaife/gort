#!/usr/bin/env python3
"""
Lidar waveform plots, gortt beside FLIGHT.

One figure per scene, three panels:

  (a) the waveforms normalised to unit total energy, against height.  This is
      the comparison that means something.  gortt's canopy return is
      J0 rho_v dP/dz, so its *shape* is dP/dz alone -- the vertical profile of
      the nadir gap probability -- and is completely independent of rho_v.
      Normalising puts FLIGHT's photon-traced waveform on the same footing and
      tests exactly that profile.

  (b) the same over the canopy only, with the ground return excluded and each
      curve renormalised, so that the canopy shape is not squashed by a ground
      return that is several times larger.

  (c) where the energy goes: the fraction returned from the ground against the
      fraction returned from the canopy.  This is the panel the backscatter
      coefficients control, and the one that caught the error in the published
      model -- see the note below.

Why the waveforms are normalised anyway
---------------------------------------
Panels (a) and (b) are normalised because the shape is dP/dz alone and carries
no information about the backscatter coefficients either way.  Panel (c) is
where those coefficients show up, and it is now a real comparison.

It did not start that way.  gortt originally used

    rho_v = 2 * t_leaf / 3          rho_g = 0.5 * r_soil

following Ni-Meister et al. (2001).  rho_v built from the leaf *transmittance*
looked wrong for a backscattered return, and panel (c) showed it: gortt's
ground fraction came out at 0.68-0.87 of FLIGHT's in every scene.  That was
recorded as a finding rather than worked around, and it was right.  Ni-Meister
et al. (2018), RSE 204, 229-243 corrects both coefficients -- the hotspot is at
scattering angle xi = pi, not 0, so rho_v takes the reflection lobe of the
phase function and not the transmission one, and rho_g is the soil albedo
outright:

    rho_v = 2 * r_leaf / 3          rho_g = r_soil

gortt now uses those (gortt_lidar.c), and panel (c) agrees with FLIGHT to
within 7% in every scene.  The change was made in the model on the paper's
authority; nothing here is tuned against FLIGHT.

What "orders 1-8" is
--------------------
FLIGHT's lidar loop runs the scattering order from 1 to BOUNCE_THRESHOLD and
accumulates a contribution for each (flight.c:6879), so the curve is the
*cumulative* waveform over orders 1 to 8, single scattering included -- not
order 8 alone.  It is the same quantity as the single-scattering curve with
more physics in it, which is why the two are plotted on the same axes.

Its visible effect is on the ground return.  Multiply scattered photons take a
longer route for the same reflecting height, and a longer route reads as a
lower height, so the ground peak grows a tail *downwards*: the standard
deviation of the sub-1 m return goes from 0.18 m at single scattering to
0.8-1.8 m over orders 1-8, and 1-9% of the total return arrives below z = -1 m,
running all the way to the floor of the recording window.  gortt has no
counterpart -- Ni-Meister et al. (2001) is a single-scattering waveform -- which
is the point of showing it.

Why the ground return is not centred on zero
--------------------------------------------
It sits at -0.198 m in every scene, and that is two separate things:

* -0.15 m of range walk, which is real.  FLIGHT's source is a point 1000 m up,
  not a collimated beam, so a photon entering the 25 m footprint off-axis
  travels sqrt(range^2 + rho^2) rather than range.  Averaged over a uniform
  disc that is E[rho^2]/(2*range) = 0.156 m of extra one-way path, and a longer
  path reads as a lower height.  Shrinking the footprint to 1 m removes it
  (measured: -0.198 -> -0.052 m), and so does moving the source to 10 km.  It
  is skewed, not symmetric, which is why the peak has a longer tail below it
  than above.

* -0.05 m of fixed offset in FLIGHT's range origin, which is not.  With a
  pinhole footprint, a near-zero pulse and no crowns the whole return lands in
  one bin whose true height is -0.125 +/- 0.025 m, i.e. about 0.25 m of two-way
  path.  It is invariant to canopy top (tested at 5, 20 and 50 m), footprint,
  sensor range, pulse width and bin count, so it is a bookkeeping constant, not
  physics; it is the same size as the 0.205 m gap FLIGHT leaves between the
  photon launch plane (canopy_top + 0.01) and its sky plane (canopy_top +
  0.215).  Part of it is cancelled by the bin index being truncated rather than
  rounded, which biases the reported height up by half a bin.

Both are well inside one 0.1 m bin of each other and neither affects a shape
comparison, but it is why the two models' ground spikes do not coincide to the
bin.  gortt's own ground return is centred on zero exactly: its grid is
anchored on z2 + 3 sigma so the delta rarely falls on a bin centre, but the
linear deposition added in gort commit 0ce47c1 splits it between the two bins
either side in proportion, which preserves the centroid.
"""

import argparse
import collections
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc

GROUND_TOP = 1.0        # m: everything below this is the ground return.  Safe
                        # for every scene here -- foliage starts at z1 = 2.0 m.


def read_gortt(path):
    """(variant, scene) -> (z, return density), plus the per-scene rho_v/rho_g/Rg."""
    out = collections.defaultdict(list)
    meta = {}
    with open(path) as fp:
        for r in csv.DictReader(fp):
            key = (r["variant"], r["scene"])
            out[key].append((float(r["z"]), float(r["return_density"])))
            meta[key] = (float(r["rho_v"]), float(r["rho_g"]), float(r["Rg"]))
    return {k: np.array(sorted(v)).T for k, v in out.items()}, meta


VARIANT_COLOURS = ["#2a78d6", "#1baf7a", "#8e44ad", "#c0392b", "#16a085"]


def read_flight(path):
    """(scene, bounce) -> {realisation: {z: energy}}"""
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    with open(path) as fp:
        for r in csv.DictReader(fp):
            out[(r["scene"], int(r["bounce"]))][int(r["realisation"])][
                round(float(r["z"]), 4)] = float(r["energy"])
    return out


def stack(realisations):
    """
    Put every realisation on one height axis.

    Only non-empty bins are stored, and which bins are empty differs between
    realisations, so the union of the heights is taken and the gaps filled with
    zeros -- an absent bin really did receive no energy.
    """
    zs = sorted(set().union(*[set(d) for d in realisations.values()]))
    arr = np.zeros((len(realisations), len(zs)))
    for i, k in enumerate(sorted(realisations)):
        d = realisations[k]
        for j, z in enumerate(zs):
            arr[i, j] = d.get(z, 0.0)
    return np.array(zs), arr


def normalise(z, w):
    """Unit total energy.  Bins are uniform, so the sum is the integral / db."""
    tot = w.sum()
    return w / tot if tot > 0 else w


def partition(z, w):
    """(ground fraction, canopy fraction) of the total returned energy."""
    tot = w.sum()
    if tot <= 0:
        return float("nan"), float("nan")
    g = w[z < GROUND_TOP].sum()
    return g / tot, 1.0 - g / tot


def plot_scene(scn, gortt, gmeta, flight, outdir, style, variants=None, suffix=""):
    """
    variants is the list of gortt formulations drawn on the same axes, baseline
    first.  In practice it is the baseline alone: no -rt_* option reaches the
    lidar path, which run_lidar.gortt_waveforms explains and checks.  The
    machinery is here so that one which did would need no changes.
    """
    variants = variants or [config.BASELINE]
    cg, cf, cm = ("#2a78d6", "#eb6834", "#1baf7a") if style == "colour" \
        else ("k", "0.35", "0.6")

    gz, gw = gortt[(variants[0], scn.name)]
    gn = normalise(gz, gw)

    single = flight.get((scn.name, config.LIDAR_BOUNCE_SINGLE))
    multi = flight.get((scn.name, config.LIDAR_BOUNCE_MULTI))

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.6))

    for panel, canopy_only in ((axes[0], False), (axes[1], True)):
        def prep(z, w):
            m = z >= GROUND_TOP if canopy_only else np.ones_like(z, dtype=bool)
            return z[m], normalise(z[m], w[m])

        for n, variant in enumerate(variants):
            if (variant, scn.name) not in gortt:
                continue
            vz, vw = gortt[(variant, scn.name)]
            zz, ww = prep(vz, vw)
            colour = (VARIANT_COLOURS[n % len(VARIANT_COLOURS)]
                      if style == "colour" else ("k" if n == 0 else "0.25"))
            panel.plot(ww, zz, color=colour, lw=1.5,
                       ls="-" if n == 0 else "--",
                       label="gortt, %s" % variant)

        if single:
            fz, fa = stack(single)
            fn = np.array([normalise(fz, row) for row in fa])
            zz, mu = prep(fz, fn.mean(axis=0))
            if len(fn) > 1:
                _, sd = prep(fz, fn.std(axis=0, ddof=1))
                # clipped at zero: a return cannot be negative, and the spread is
                # wide because each realisation is a single 25 m footprint
                panel.fill_betweenx(zz, np.maximum(mu - sd, 0.0), mu + sd,
                                    color=cf, alpha=0.25, lw=0)
            panel.plot(mu, zz, color=cf, lw=1.3,
                       label="FLIGHT, single scattering")
        if multi:
            mz, ma = stack(multi)
            mn = np.array([normalise(mz, row) for row in ma])
            zz, mu = prep(mz, mn.mean(axis=0))
            panel.plot(mu, zz, color=cm, lw=1.1, ls="--",
                       label="FLIGHT, orders 1-8")
            if not canopy_only:
                deep = mn.mean(axis=0)[mz < -1.0].sum()

        for zr, lab in ((scn.z1, "z1"), (scn.h1, "h1"), (scn.h2, "h2"), (scn.z2, "z2")):
            panel.axhline(zr, color="0.75", lw=0.6, ls=":")
            panel.text(panel.get_xlim()[1], zr, lab, fontsize=6, color="0.45",
                       ha="right", va="bottom")
        # The lower limit used to be -1 m, which cropped the whole multiple
        # scattering tail: photons that bounce inside the canopy before or
        # after reaching the ground arrive late, and a late arrival reads as a
        # height below the ground.  It runs to the floor of the recording
        # window, so the full waveform panel shows most of it.
        if canopy_only:
            panel.set_ylim(GROUND_TOP - 0.5, scn.z2 + 2.0)
        else:
            panel.set_ylim(-6.0, scn.z2 + 2.0)
            if multi:
                panel.text(0.97, 0.02,
                           "orders 1-8: %.1f%% of the return below $-$1 m,\n"
                           "tail continues to the window floor at %g m"
                           % (100.0 * deep, config.LIDAR_REC_LOW),
                           transform=panel.transAxes, fontsize=6, color=cm,
                           ha="right", va="bottom")
        panel.set_xlabel("normalised return per bin")
        panel.set_ylabel("height above ground (m)")
        panel.tick_params(direction="in", top=True, right=True)
        panel.set_title("canopy only, renormalised" if canopy_only
                        else "full waveform, unit total energy", fontsize=9)
        panel.axhline(0.0, color="0.55", lw=0.7)

    axes[0].legend(fontsize=7, frameon=False, loc="upper right")

    # (c) where the energy goes
    ax = axes[2]
    labels, ground, canopy = [], [], []
    for variant in variants:
        if (variant, scn.name) not in gortt:
            continue
        vz, vw = gortt[(variant, scn.name)]
        g, c = partition(vz, vw)
        labels.append("gortt\n%s" % variant)
        ground.append(g)
        canopy.append(c)
    if single:
        fz, fa = stack(single)
        parts = np.array([partition(fz, row) for row in fa])
        labels.append("FLIGHT\nsingle")
        ground.append(parts[:, 0].mean())
        canopy.append(parts[:, 1].mean())
    if multi:
        mz, ma = stack(multi)
        parts = np.array([partition(mz, row) for row in ma])
        labels.append("FLIGHT\norders 1-8")
        ground.append(parts[:, 0].mean())
        canopy.append(parts[:, 1].mean())

    x = np.arange(len(labels))
    ax.bar(x, ground, 0.55, label="ground (z < %g m)" % GROUND_TOP, color=cg)
    ax.bar(x, canopy, 0.55, bottom=ground, label="canopy", color=cf)
    for xi, g in zip(x, ground):
        ax.text(xi, g / 2.0, "%.2f" % g, ha="center", va="center",
                fontsize=8, color="w")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of returned energy")
    ax.set_title("energy partition\n(set by rho_g/rho_v -- see docstring)",
                 fontsize=8)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.tick_params(direction="in", right=True)

    rho_v, rho_g, _ = gmeta[(variants[0], scn.name)]
    fig.suptitle("%s   nadir lidar waveform at %d nm\n"
                 "r=%.1f b=%.1f m, $\\lambda$=%.4f m$^{-2}$, favd=%.2f m$^{-1}$, "
                 "PCC=%.2f, target LAI=%.1f   |   gortt $\\rho_v$=%.3f $\\rho_g$=%.3f"
                 % (scn.name, config.LIDAR_BAND.wl, scn.r, scn.b, scn.lam,
                    scn.favd, scn.pcc, scn.target_lai, rho_v, rho_g),
                 fontsize=9)
    fig.tight_layout()
    out = os.path.join(outdir, "%s-lidar%s.png" % (scn.name, suffix))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def statistics(scenes, gortt, flight, path, variants=None):
    """
    How close are the normalised shapes?  Compared on the canopy only, where
    the two models are claiming the same thing, and summarised by the mean
    return height as well as by an rms difference.

    One row per (variant, scene).  FLIGHT is the same single-scattering Monte
    Carlo in every row: the -rt_* options are gortt-side formulations and the
    reference does not move.
    """
    variants = variants or [config.BASELINE]
    with open(path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["variant", "scene", "ground_frac_gortt", "ground_frac_flight",
                    "mean_height_gortt", "mean_height_flight",
                    "canopy_shape_rmsd", "peak_z_gortt", "peak_z_flight"])
        for variant in variants:
            for scn in scenes:
                if (variant, scn.name) not in gortt:
                    continue
                gz, gw = gortt[(variant, scn.name)]
                single = flight.get((scn.name, config.LIDAR_BOUNCE_SINGLE))
                if not single:
                    continue
                fz, fa = stack(single)
                fw = fa.mean(axis=0)

                gg, _ = partition(gz, gw)
                fg, _ = partition(fz, fw)

                gm = gz >= GROUND_TOP
                fm = fz >= GROUND_TOP
                gc = normalise(gz[gm], gw[gm])
                fc = normalise(fz[fm], fw[fm])
                # FLIGHT's grid is a subset of the same 0.1 m bins, so compare
                # by interpolating gortt's canopy shape onto FLIGHT's heights.
                gi = np.interp(fz[fm], gz[gm][::-1], gc[::-1], left=0.0, right=0.0)
                rmsd = float(np.sqrt(np.mean((gi - fc) ** 2)))

                w.writerow([variant, scn.name, "%.4f" % gg, "%.4f" % fg,
                            "%.3f" % float((gz[gm] * gc).sum()),
                            "%.3f" % float((fz[fm] * fc).sum()),
                            "%.5f" % rmsd,
                            "%.2f" % gz[gm][gc.argmax()],
                            "%.2f" % fz[fm][fc.argmax()]])
    sys.stderr.write("plot_lidar: wrote %s\n" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=config.RESULTS_DIR,
                    help="directory holding the CSVs (default: %(default)s)")
    ap.add_argument("--outdir", default=config.FIGS_DIR)
    ap.add_argument("--style", choices=("colour", "paper"), default="colour")
    ap.add_argument("--scenes", nargs="*")
    ap.add_argument("--variants", nargs="+", metavar="NAME",
                    default=[config.BASELINE],
                    help="gortt variants to draw (default: just the baseline, "
                         "since no -rt_* option reaches the lidar path). Any "
                         "name not present in the CSV is skipped.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gortt, gmeta = read_gortt(os.path.join(args.results, "lidar_gortt.csv"))
    flight = read_flight(os.path.join(args.results, "lidar_flight.csv"))

    scenes = config.scenes()
    if args.scenes:
        scenes = [config.scene_by_name(n) for n in args.scenes]
    have = {k[1] for k in gortt}
    scenes = [s for s in scenes if s.name in have]
    if not scenes:
        sys.exit("plot_lidar: no results found -- run run_lidar.py first")

    args.variants = [v for v in args.variants if any(v == k[0] for k in gortt)]
    for variant in args.variants:
        suffix = "" if variant == config.BASELINE else "-" + variant
        pair = ([config.BASELINE] if variant == config.BASELINE
                else [config.BASELINE, variant])
        for scn in scenes:
            sys.stderr.write("plot_lidar: %s\n"
                             % plot_scene(scn, gortt, gmeta, flight, args.outdir,
                                          args.style, variants=pair, suffix=suffix))

    statistics(scenes, gortt, flight,
               os.path.join(args.results, "lidar_stats.csv"),
               variants=args.variants)


if __name__ == "__main__":
    main()
