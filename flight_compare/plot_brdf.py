#!/usr/bin/env python3
"""
Polar BRDF plots and principal-plane cross sections, gortt beside FLIGHT.

Two figures per scene:

  <scene>-polar-<band>.png   3 rows (solar zenith) x 3 columns (gortt, FLIGHT,
                             gortt - FLIGHT).  Radius is view zenith, angle is
                             relative azimuth with the hotspot at the top, so
                             the backscatter half is the upper semicircle and
                             the forward-scatter half the lower one.  The two
                             model columns share a colour scale within a row so
                             that they can be read against each other.

  <scene>-pp.png             the principal plane, 2 bands x 3 solar zeniths.
                             Signed view zenith: positive is backscatter
                             (relative azimuth 0), following the convention of
                             figures/gortt_run.py and of the GORT papers, which
                             put the hotspot at positive view zenith.  gortt is
                             a line on its own fine grid; FLIGHT is points with
                             the standard error across scene realisations.

The polar panels contour FLIGHT only on the regular (view zenith x relative
azimuth) grid it was run on.  The extra principal-plane directions are left out
of the contouring -- they are there to resolve the hotspot in the line plots,
and feeding a partly irregular grid to a contour routine would invent structure
that was never sampled.
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


def read_gortt(path):
    """(variant, scene, band, sza, grid) -> {(vza, raa): brf}"""
    out = collections.defaultdict(dict)
    with open(path) as fp:
        for r in csv.DictReader(fp):
            key = (r["variant"], r["scene"], r["band"], float(r["sza"]), r["grid"])
            out[key][(float(r["vza"]), float(r["raa"]))] = float(r["brf"])
    return out


# Colours for the gortt curves when more than one variant is drawn. The
# baseline keeps the blue it has always had and FLIGHT keeps its orange, so a
# variant page reads the same way as the main one with an extra line on it.
VARIANT_COLOURS = ["#2a78d6", "#1baf7a", "#8e44ad", "#c0392b", "#16a085"]
FLIGHT_COLOUR = "#eb6834"


def read_flight(path):
    """(scene, band, sza) -> {(vza, raa): [brf per realisation]}"""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(path) as fp:
        for r in csv.DictReader(fp):
            key = (r["scene"], r["band"], float(r["sza"]))
            out[key][(float(r["vza"]), float(r["raa"]))].append(float(r["brf"]))
    return out


def polar_field(lookup, vzas, raas, reduce=np.mean):
    """
    Values on the (vza x raa) grid, mirrored into a full circle.

    raa 0 is backscatter and sits at the top; the canopy is statistically
    azimuthally symmetric about the principal plane, so raa and -raa are the
    same direction and the mirrored half is not new information -- it is drawn
    because a half-disc is much harder to read.
    """
    z = np.full((len(vzas), len(raas)), np.nan)
    for i, v in enumerate(vzas):
        for j, a in enumerate(raas):
            key = (v, 0.0 if v == 0.0 else a)
            vals = lookup.get(key)
            if vals is not None:
                z[i, j] = reduce(vals) if isinstance(vals, list) else vals

    raas = np.asarray(raas, dtype=float)
    full_raa = np.concatenate([-raas[::-1][:-1], raas])
    full_z = np.concatenate([z[:, ::-1][:, :-1], z], axis=1)
    return np.asarray(vzas, dtype=float), full_raa, full_z


def polar_panel(ax, vzas, raas, z, levels, cmap, title):
    """
    One polar panel.  Relative azimuth 0 -- the hotspot -- is at the top and
    azimuth increases clockwise, so the upper half of the disc is backscatter
    and the lower half forward scatter.  The radius is view zenith.
    """
    th, rr = np.meshgrid(np.radians(raas), vzas)
    # levels already span the full range of both panels, so no extend arrows --
    # they would collide with the end tick labels on the colorbar.
    cs = ax.contourf(th, rr, z, levels=levels, cmap=cmap, extend="neither")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, max(vzas))
    ax.set_rticks([20, 40, 60])
    ax.set_yticklabels([])           # the radial scale is in the panel caption
    ax.set_xticks(np.radians([0, 90, 180, 270]))
    # 90 and 270 are the same plane by symmetry; labelling only the principal
    # plane keeps neighbouring panels from running their tick labels together.
    ax.set_xticklabels(["hotspot", "", "forward", ""], fontsize=6, color="0.35")
    ax.grid(color="0.6", lw=0.4, alpha=0.6)
    ax.set_title(title, fontsize=8, pad=6)
    return cs


def add_colorbar(fig, cs, ax, label):
    cb = fig.colorbar(cs, ax=ax, fraction=0.05, pad=0.06)
    cb.ax.tick_params(labelsize=6)
    cb.locator = matplotlib.ticker.MaxNLocator(5)
    cb.formatter = matplotlib.ticker.FormatStrFormatter("%.3f")
    cb.update_ticks()
    cb.set_label(label, fontsize=6)
    return cb


def plot_polar(scn, band, gortt, flight, outdir, style, variant=None, suffix=""):
    """
    Three polar panels per solar zenith: gortt, FLIGHT, and the difference.

    Only one gortt variant fits here -- the layout is a fixed three columns on a
    shared colour scale, and a fourth field would need its own scale to be
    readable. A variant page gets its own file via `suffix` instead, and the
    difference panel is the one to read: it is the same FLIGHT either way.
    """
    variant = variant or config.BASELINE
    cmap = "viridis" if style == "colour" else "Greys"
    vzas, raas = config.POLAR_VZA, config.POLAR_RAA

    fig, axes = plt.subplots(len(config.SZAS), 3,
                             figsize=(10.5, 3.2 * len(config.SZAS)),
                             subplot_kw={"projection": "polar"},
                             layout="constrained")
    if len(config.SZAS) == 1:
        axes = np.array([axes])

    for row, sza in enumerate(config.SZAS):
        gk = (variant, scn.name, band.name, sza, "dense")
        fk = (scn.name, band.name, sza)
        if gk not in gortt or fk not in flight:
            continue

        gv, gr, gz = polar_field(gortt[gk], config.GORTT_DENSE_VZA, config.GORTT_DENSE_RAA)
        fv, fr, fz = polar_field(flight[fk], vzas, raas)
        # gortt sampled on FLIGHT's own directions, so the difference is like for like
        dv, dr, dz = polar_field(gortt[(variant, scn.name, band.name, sza, "flight")],
                                 vzas, raas)
        diff = dz - fz

        lo = np.nanmin([np.nanmin(gz), np.nanmin(fz)])
        hi = np.nanmax([np.nanmax(gz), np.nanmax(fz)])
        levels = np.linspace(lo, hi, 21)
        dmax = np.nanmax(np.abs(diff)) or 1e-6
        dlevels = np.linspace(-dmax, dmax, 21)

        polar_panel(axes[row][0], gv, gr, gz, levels, cmap,
                    "gortt   sza %g$^\\circ$" % sza)
        c2 = polar_panel(axes[row][1], fv, fr, fz, levels, cmap,
                         "FLIGHT   sza %g$^\\circ$" % sza)
        c3 = polar_panel(axes[row][2], dv, dr, diff, dlevels, "RdBu_r",
                         "gortt $-$ FLIGHT")

        # Columns 1 and 2 share a scale so they can be read against each other;
        # the difference gets its own, centred on zero.
        add_colorbar(fig, c2, axes[row][1], "BRF")
        add_colorbar(fig, c3, axes[row][2], "$\\Delta$BRF")

    fig.suptitle("%s   %s (%d nm)   BRF\n"
                 "r=%.1f b=%.1f m, $\\lambda$=%.4f m$^{-2}$, favd=%.2f m$^{-1}$, "
                 "PCC=%.2f, target LAI=%.1f"
                 % (scn.name, band.name.upper(), band.wl, scn.r, scn.b, scn.lam,
                    scn.favd, scn.pcc, scn.target_lai)
                 + "\nrings at view zenith 20, 40, 60$^\\circ$",
                 fontsize=9)
    out = os.path.join(outdir, "%s-polar-%s.png" % (scn.name, band.name))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def signed_pp(lookup, reduce=None):
    """
    (signed vza, value) along the principal plane.  Positive is backscatter
    (raa 0), negative is forward scatter (raa 180).
    """
    pts = []
    for (vza, raa), val in lookup.items():
        if raa not in (0.0, 180.0):
            continue
        v = reduce(val) if reduce else val
        pts.append((vza if raa == 0.0 else -vza, v))
    pts.sort()
    x, y = zip(*pts) if pts else ([], [])
    return np.array(x), np.array(y)


def plot_pp(scn, gortt, flight, outdir, style, variants=None, suffix=""):
    """
    The principal plane, one panel per (band, solar zenith).

    variants is a list of gortt variant names drawn on the same axes; the
    default is the baseline alone, which is what the main page uses. A variant
    page passes [baseline, <one option>] so that the change has something to be
    a change from -- FLIGHT is the same Monte Carlo in both cases and is never
    re-run.
    """
    variants = variants or [config.BASELINE]
    fcolour = FLIGHT_COLOUR if style == "colour" else "0.45"
    fig, axes = plt.subplots(len(config.BRDF_BANDS), len(config.SZAS),
                             figsize=(3.6 * len(config.SZAS), 3.1 * len(config.BRDF_BANDS)),
                             squeeze=False)

    for i, band in enumerate(config.BRDF_BANDS):
        for j, sza in enumerate(config.SZAS):
            ax = axes[i][j]
            fk = (scn.name, band.name, sza)
            for n, variant in enumerate(variants):
                gk = (variant, scn.name, band.name, sza, "dense")
                if gk not in gortt:
                    continue
                x, y = signed_pp(gortt[gk])
                colour = (VARIANT_COLOURS[n % len(VARIANT_COLOURS)]
                          if style == "colour" else ("k" if n == 0 else "0.25"))
                ax.plot(x, y, color=colour, lw=1.4,
                        ls="-" if n == 0 else "--",
                        label="gortt, %s" % variant)
            if fk in flight:
                x, y = signed_pp(flight[fk], reduce=np.mean)
                _, e = signed_pp(flight[fk],
                                 reduce=lambda v: np.std(v, ddof=1) / np.sqrt(len(v))
                                 if len(v) > 1 else 0.0)
                ax.errorbar(x, y, yerr=e, fmt="o", ms=3.0, lw=0.8, capsize=2,
                            color=fcolour, label="FLIGHT")
            ax.axvline(sza, color="0.7", lw=0.6, ls=":")
            ax.annotate("hotspot", xy=(sza, 0.97), xycoords=("data", "axes fraction"),
                        fontsize=6, color="0.45", ha="left", va="top", rotation=90)
            ax.set_title("%s, sza %g$^\\circ$" % (band.name.upper(), sza), fontsize=8)
            ax.tick_params(direction="in", top=True, right=True, labelsize=7)
            if i == len(config.BRDF_BANDS) - 1:
                ax.set_xlabel("view zenith (deg)\n$+$ backscatter    $-$ forward", fontsize=7)
            if j == 0:
                ax.set_ylabel("BRF", fontsize=8)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, frameon=False)

    fig.suptitle("%s   principal plane\n"
                 "r=%.1f b=%.1f m, $\\lambda$=%.4f m$^{-2}$, favd=%.2f m$^{-1}$, "
                 "PCC=%.2f, target LAI=%.1f"
                 % (scn.name, scn.r, scn.b, scn.lam, scn.favd, scn.pcc, scn.target_lai),
                 fontsize=9)
    fig.tight_layout()
    out = os.path.join(outdir, "%s-pp%s.png" % (scn.name, suffix))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_overview(scenes, gortt, flight, outdir, style, variant=None, suffix=""):
    """
    One picture of the whole comparison: every direction in every scene, gortt
    against FLIGHT, and the relative RMSD broken down by scene.

    The scatter is the honest summary -- it shows the spread as well as the
    slope, and it separates the two bands, which behave quite differently.
    """
    cols = {"red": "#c1443c", "nir": "#1baf7a"} if style == "colour" \
        else {"red": "0.2", "nir": "0.6"}

    variant = variant or config.BASELINE

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8),
                             gridspec_kw={"width_ratios": [1.0, 1.45]})

    ax = axes[0]
    lim = 0.0
    for band in config.BRDF_BANDS:
        gx, fy = [], []
        for scn in scenes:
            for sza in config.SZAS:
                gk = (variant, scn.name, band.name, sza, "flight")
                fk = (scn.name, band.name, sza)
                if gk not in gortt or fk not in flight:
                    continue
                for k in sorted(set(gortt[gk]) & set(flight[fk])):
                    gx.append(gortt[gk][k])
                    fy.append(np.mean(flight[fk][k]))
        if not gx:
            continue
        lim = max(lim, max(gx), max(fy))
        ax.plot(fy, gx, ".", ms=2.2, alpha=0.35, color=cols[band.name],
                label="%s (%d nm)" % (band.name.upper(), band.wl))
    lim *= 1.06
    ax.plot([0, lim], [0, lim], "-", lw=0.8, color="0.4")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("FLIGHT BRF")
    ax.set_ylabel("gortt BRF")
    ax.set_title("every direction, every scene", fontsize=9)
    ax.tick_params(direction="in", top=True, right=True)
    leg = ax.legend(fontsize=8, frameon=False, loc="upper left")
    for h in leg.legend_handles:
        h.set_markersize(7)
        h.set_alpha(1.0)

    ax = axes[1]
    names = [s.name for s in scenes]
    width = 0.38
    for i, band in enumerate(config.BRDF_BANDS):
        vals = []
        for scn in scenes:
            d, f = [], []
            for sza in config.SZAS:
                gk = (variant, scn.name, band.name, sza, "flight")
                fk = (scn.name, band.name, sza)
                if gk not in gortt or fk not in flight:
                    continue
                for k in sorted(set(gortt[gk]) & set(flight[fk])):
                    d.append(gortt[gk][k] - np.mean(flight[fk][k]))
                    f.append(np.mean(flight[fk][k]))
            vals.append(np.sqrt(np.mean(np.square(d))) / np.mean(f) if d else np.nan)
        ax.barh(np.arange(len(names)) + (i - 0.5) * width, vals, width,
                color=cols[band.name], label=band.name.upper())
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels(names, fontsize=8, family="monospace")
    ax.invert_yaxis()
    ax.set_xlabel("RMSD / mean FLIGHT BRF")
    ax.set_title("disagreement by scene", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax.tick_params(direction="in", top=True)
    ax.grid(axis="x", color="0.85", lw=0.5)
    ax.set_axisbelow(True)

    fig.suptitle("gortt against FLIGHT: all scenes, all view directions, "
                 "solar zeniths %s$^\\circ$"
                 % ", ".join("%g" % s for s in config.SZAS), fontsize=10)
    fig.tight_layout()
    out = os.path.join(outdir, "overview%s.png" % suffix)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def statistics(scenes, gortt, flight, path, variants=None):
    """
    Agreement between the two models on FLIGHT's own directions, per variant.

    FLIGHT is averaged over its realisations and is identical for every variant,
    so a row is a fair comparison of one gortt formulation against one fixed
    reference. The variant column is what the comparison pages rank on.
    """
    variants = variants or [v[0] for v in config.GORTT_VARIANTS]
    with open(path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["variant", "scene", "band", "sza", "n", "mean_flight",
                    "mean_gortt", "bias", "rmsd", "rel_rmsd", "r"])
        for variant in variants:
            for scn in scenes:
                for band in config.BRDF_BANDS:
                    for sza in config.SZAS:
                        gk = (variant, scn.name, band.name, sza, "flight")
                        fk = (scn.name, band.name, sza)
                        if gk not in gortt or fk not in flight:
                            continue
                        keys = sorted(set(gortt[gk]) & set(flight[fk]))
                        g = np.array([gortt[gk][k] for k in keys])
                        f = np.array([np.mean(flight[fk][k]) for k in keys])
                        d = g - f
                        rms = float(np.sqrt(np.mean(d ** 2)))
                        r = (float(np.corrcoef(g, f)[0, 1]) if len(g) > 2
                             else float("nan"))
                        w.writerow([variant, scn.name, band.name, "%g" % sza,
                                    len(keys), "%.5f" % f.mean(), "%.5f" % g.mean(),
                                    "%.5f" % d.mean(), "%.5f" % rms,
                                    "%.4f" % (rms / f.mean()), "%.4f" % r])
    sys.stderr.write("plot_brdf: wrote %s\n" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=config.RESULTS_DIR,
                    help="directory holding the CSVs (default: %(default)s)")
    ap.add_argument("--outdir", default=config.FIGS_DIR,
                    help="where to write the images (default: %(default)s)")
    ap.add_argument("--style", choices=("colour", "paper"), default="colour",
                    help="colour: viridis polar panels (default: %(default)s)")
    ap.add_argument("--scenes", nargs="*", help="limit to these scenes")
    ap.add_argument("--variants", nargs="+", metavar="NAME",
                    default=[v[0] for v in config.GORTT_VARIANTS],
                    help="gortt variants to draw (default: all of them). The "
                         "baseline keeps the unsuffixed filenames; every other "
                         "variant gets its own -<name> set, drawn against the "
                         "baseline.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gortt = read_gortt(os.path.join(args.results, "brdf_gortt.csv"))
    flight = read_flight(os.path.join(args.results, "brdf_flight.csv"))

    scenes = config.scenes()
    if args.scenes:
        scenes = [config.scene_by_name(n) for n in args.scenes]
    have = {k[0] for k in flight}
    scenes = [s for s in scenes if s.name in have]
    if not scenes:
        sys.exit("plot_brdf: no FLIGHT results found -- run run_brdf.py first")

    # The baseline figures keep their unsuffixed names, so the main page does
    # not move. Each variant adds its own -<name> set beside them and is drawn
    # against the baseline, since the interesting thing is the change.
    for variant in args.variants:
        suffix = "" if variant == config.BASELINE else "-" + variant
        pair = ([config.BASELINE] if variant == config.BASELINE
                else [config.BASELINE, variant])
        for scn in scenes:
            # Polar contours are baseline only. The layout is a fixed three
            # columns on one colour scale, so it can show one formulation at a
            # time anyway, and a full set per variant would be 144 more PNGs
            # and about 36 MB in a repository that already carries 13. The
            # principal-plane panels put both curves on the same axes, which is
            # what a variant page actually wants.
            if variant == config.BASELINE:
                for band in config.BRDF_BANDS:
                    sys.stderr.write("plot_brdf: %s\n"
                                     % plot_polar(scn, band, gortt, flight,
                                                  args.outdir, args.style,
                                                  variant=variant, suffix=suffix))
            sys.stderr.write("plot_brdf: %s\n"
                             % plot_pp(scn, gortt, flight, args.outdir, args.style,
                                       variants=pair, suffix=suffix))
        sys.stderr.write("plot_brdf: %s\n"
                         % plot_overview(scenes, gortt, flight, args.outdir,
                                         args.style, variant=variant, suffix=suffix))

    statistics(scenes, gortt, flight,
               os.path.join(args.results, "brdf_stats.csv"),
               variants=args.variants)


if __name__ == "__main__":
    main()
