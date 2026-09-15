#!/usr/bin/env python3
"""
The experiments behind KC-INVESTIGATION.md: gortt's geometric-optical kernels,
and the alternatives to them, scored against a Monte Carlo reference.

  usage:  ./kc_investigation.py [--outdir DIR] [--jobs 4] [--refresh]

Monte Carlo results are cached in DIR/mc_cache.csv (keyed by canopy and
geometry, fixed seeds), so only the first run is slow -- about five minutes
on four cores.  --refresh throws the cache away.  Writes, into DIR:

  tables.md          every table quoted in KC-INVESTIGATION.md
  fig3_P_terms.png   Pv and Pi at the parameters of L&S Fig. 3/4
  nesting.png        MC PvMv, PiMi, Po against min(PvMv, PiMi)
  kg_overlap.png     Kg against MC for the overlap height choices
  kc_variants.png    Kc against MC for the P-term / beta / F choices
  fig5.png           Ni et al. (1999) Fig. 5, with MC

See goms_kernels.py for the variant switches and goms_mc.py for the MC.
"""

import argparse
import csv
import dataclasses
import os
import sys
import zlib
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import goms_kernels as gk
import goms_mc as mc

R = np.radians
MC_RAYS, MC_ELEM, MC_SCENES = 40000, 16000, 8

CANOPIES = {
    "SOBS":        gk.SOBS,
    "same-height": dataclasses.replace(gk.SOBS, h1=5.75, h2=5.75),
    "sparse":      dataclasses.replace(gk.SOBS, lam=0.1),
    "dense":       dataclasses.replace(gk.SOBS, lam=1.0),
    "wide-spread": dataclasses.replace(gk.SOBS, h1=3.0, h2=14.0),
}
# L&S (1992) Figs 3-4: 30 crowns of 3 m radius per 900 m2, b/R = 1.5, b/h = 0.18,
# h = 25 m, and (h2-h1)/h = 0, 0.2, 0.4
LS92 = {
    "L&S (a) h2-h1=0": gk.Canopy(r=3.0, b=4.5, lam=30.0 / 900.0, h1=25.0, h2=25.0),
    "L&S (b) h2-h1=5": gk.Canopy(r=3.0, b=4.5, lam=30.0 / 900.0, h1=22.5, h2=27.5),
    "L&S (c) h2-h1=10": gk.Canopy(r=3.0, b=4.5, lam=30.0 / 900.0, h1=20.0, h2=30.0),
}
ALL_CANOPIES = {**CANOPIES, **LS92}

SZAS = (15.0, 30.0, 45.0, 60.0)
PP = np.arange(-85.0, 85.1, 5.0)

BASE = gk.GORTT.replace(overlap_h="avg")
KC_VARIANTS = [
    ("gortt, F squared (master)",      gk.GORTT.replace(fsq="old")),
    ("gortt (this branch)",            gk.GORTT),
    ("+ overlap h averaged",           BASE),
    ("+ Po = min(PvMv, PiMi)",         BASE.replace(po="nest")),
    ("  ... beta = eq. 23",            BASE.replace(po="nest", beta="eq23")),
    ("  ... beta = 1",                 BASE.replace(po="nest", beta=1.0)),
    ("lune model",                     BASE.replace(lune=True)),
    ("eqs 21-22 as printed",           BASE.replace(pimi="eq21", pvmv="eq22")),
    ("gortt_kc_PP_only as coded",      BASE.replace(pimi="eq21raw", pvmv="ppcoded")),
    ("F squared + h avg + Po nest",    BASE.replace(po="nest", fsq="old")),
]
KG_VARIANTS = [(h, az) for h in ("depth", "centre", "avg") for az in ("ambrals", "pp", "eq7", "wenge")]

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]


# ---------------------------------------------------------------------------
# geometry grid and the Monte Carlo cache
# ---------------------------------------------------------------------------

def geometries():
    """(canopy, vza, sza, raa) in degrees, gortt conventions, for every MC run."""
    out = []
    for name in CANOPIES:
        for s in SZAS:
            for v in PP:
                out.append((name, abs(v), s, 180.0 if v < 0 else 0.0))
            for ra in np.arange(10.0, 171.0, 10.0):
                out.append((name, s, s, ra))
            for ra in np.arange(15.0, 166.0, 15.0):
                out.append((name, 45.0, s, ra))
    for name, c in LS92.items():
        for vp in np.arange(-50.0, 60.1, 5.0):
            v = np.degrees(np.arctan(np.tan(R(abs(vp))) * c.r / c.b))
            out.append((name, round(v, 10), 30.0, 180.0 if vp < 0 else 0.0))
    return sorted(set(out))


def _mc_job(g):
    name, v, s, ra = g
    seed = zlib.crc32(("%s|%r|%r|%r" % (name, v, s, ra)).encode()) % 100000
    o = mc.simulate(ALL_CANOPIES[name], R(v), R(s), R(ra), MC_RAYS, MC_ELEM, seed=seed,
                    n_scenes=MC_SCENES)
    return g, o


def load_mc(path, jobs, refresh):
    cache = {}
    if os.path.exists(path) and not refresh:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                key = (row["canopy"], float(row["vza"]), float(row["sza"]), float(row["raa"]))
                cache[key] = {k: float(x) for k, x in row.items()
                              if k not in ("canopy", "vza", "sza", "raa")}
    todo = [g for g in geometries() if g not in cache]
    if todo:
        print("Monte Carlo: %d geometries to run (%d cached)" % (len(todo), len(cache)))
        with Pool(jobs) as pool:
            for k, (g, o) in enumerate(pool.imap_unordered(_mc_job, todo, chunksize=4)):
                cache[g] = o
                if (k + 1) % 100 == 0:
                    print("  %d/%d" % (k + 1, len(todo)))
        keys = sorted(next(iter(cache.values())).keys())
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["canopy", "vza", "sza", "raa"] + keys)
            for g in sorted(cache):
                w.writerow(list(g) + [repr(float(cache[g][k])) for k in keys])
    return cache


def select(cache, canopy, pred=lambda v, s, ra: True):
    keys = [g for g in sorted(cache) if g[0] == canopy and pred(*g[1:])]
    v = np.array([g[1] for g in keys])
    s = np.array([g[2] for g in keys])
    ra = np.array([g[3] for g in keys])
    return keys, v, s, ra


def col(cache, keys, name):
    return np.array([cache[k][name] for k in keys])


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def md_table(header, rows, fmt):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" if i == 0 else "---:" for i in range(len(header))) + "|"]
    for r in rows:
        out.append("| " + " | ".join([str(r[0])] + [f % x if isinstance(x, (float, np.floating))
                                                      else str(x) for f, x in zip(fmt, r[1:])]) + " |")
    return "\n".join(out)


def rmse(e):
    return float(np.sqrt(np.mean(np.square(e))))


def paper_P(c, vza, sza, raa):
    """L&S (1992) Pv, Pi by case: eqs 19/21/22 with Pv = 1 in (c) and Pi = 1 in (b)."""
    ip, vp = gk.prime(c, sza), gk.prime(c, vza)
    Mi, Mv = gk.mutual(c, ip), gk.mutual(c, vp)
    tMi, tMv = np.arccos(1 - 2 * Mi), np.arccos(1 - 2 * Mv)
    u = np.where(raa > np.pi / 2, -vp, vp)
    phase = np.cos(vp) * np.cos(ip) + np.sin(vp) * np.sin(ip) * np.cos(raa)
    Gv = np.pi * c.rr / np.cos(vp)
    Gc = Gv * 0.5 * (1 + phase)
    pv19 = (Mv * Gv - (Gv - Gc)) / (Mv * Gv)
    arg21 = tMi - ip + u
    pi21 = np.where(arg21 > 0, (1 - np.cos(arg21)) / (1 - np.cos(tMi)), 0.0)
    arg22 = tMv - ip + u
    pv22 = np.where(arg22 >= 0, (1 - np.cos(arg22)) / (1 - np.cos(tMv)), 0.0)
    Pv = np.where(u >= ip, pv19, np.where(u >= 0, 1.0, pv22))
    Pi = np.where(u >= ip, 1.0, pi21)
    return Pv, Pi


def experiment_kg(cache):
    lines = ["### Kg: overlap height and azimuthal form against Monte Carlo",
             "",
             "RMSE of Kg over every geometry for the canopy (principal plane, the principal cone "
             "vza = sza, and vza = 45 across azimuth; sza 15-60). In brackets, the principal cone "
             "alone, which is where the azimuthal form of the hotspot matters.", ""]
    rows = []
    for h, az in KG_VARIANTS:
        row = ["h=%s, %s" % (h, az)]
        for name, c in CANOPIES.items():
            keys, v, s, ra = select(cache, name)
            cone = (v == s) & (ra > 0) & (ra < 180)
            kg = gk.kg(c, gk.GORTT.replace(overlap_h=h, overlap_az=az), R(v), R(s), R(ra))
            e = kg - col(cache, keys, "Kg")
            bad = ~np.isfinite(e)
            e = np.where(bad, 0.0, e)
            ok = ~bad
            row.append("%.4f (%.4f)%s" % (rmse(e[ok]), rmse(e[cone & ok]),
                                          " [%d NaN]" % bad.sum() if bad.any() else ""))
        rows.append(row)
    lines.append(md_table(["variant"] + list(CANOPIES), rows, ["%s"] * len(CANOPIES)))
    se = max(np.max(col(cache, select(cache, n)[0], "Kg_se")) for n in CANOPIES)
    lines += ["", "Largest MC standard error on Kg in these runs: %.4f. `wenge` is NaN at the "
              "hotspot itself, where it divides by D = 0; its RMSE is over the other geometries."
              % se, ""]
    return lines


def experiment_kc(cache):
    lines = ["### Kc: P terms, Po, beta and F against Monte Carlo", "",
             "RMSE of Kc per canopy over the same geometries (sza 15-60 pooled). "
             "`clamp` counts geometries where Kc > 1 - Kg - Kz, i.e. where gortt_rsurf "
             "would clamp Kt to zero; the MC never does.", ""]
    rows = []
    for label, var in KC_VARIANTS:
        row = [label]
        allerr, nclamp = [], 0
        for name, c in CANOPIES.items():
            keys, v, s, ra = select(cache, name)
            K = gk.proportions(c, var, R(v), R(s), R(ra))
            e = K["Kc"] - col(cache, keys, "Kc")
            allerr.append(e)
            nclamp += int(np.sum(K["Kt_raw"] < -1e-12))
            row.append(rmse(e))
        e = np.concatenate(allerr)
        row += [rmse(e), float(np.max(np.abs(e))), nclamp]
        rows.append(row)
    lines.append(md_table(["variant"] + list(CANOPIES) + ["all", "max abs", "clamp"], rows,
                          ["%.3f"] * (len(CANOPIES) + 2) + ["%d"]))

    lines += ["", "Split by direction for SOBS, variant `+ Po = min(PvMv, PiMi)` against "
              "gortt on this branch: RMSE in backscatter (principal plane, vza > sza), near "
              "nadir (|vza| <= 20) and forward scatter (principal plane, vza >= 25).", ""]
    rows = []
    c = CANOPIES["SOBS"]
    for s0 in SZAS:
        keys, v, s, ra = select(cache, "SOBS", lambda v, s, ra, s0=s0: s == s0 and ra in (0.0, 180.0))
        mcK = col(cache, keys, "Kc")
        signed = np.where(ra > 90, -v, v)
        groups = [("back", signed > s0), ("nadir", np.abs(signed) <= 20), ("fwd", signed <= -25)]
        row = ["sza %g" % s0]
        for var in (gk.GORTT, BASE.replace(po="nest")):
            K = gk.proportions(c, var, R(v), R(s), R(ra))["Kc"]
            row += [rmse((K - mcK)[m]) for _, m in groups]
        rows.append(row)
    lines.append(md_table(["", "gortt back", "gortt nadir", "gortt fwd",
                           "nest back", "nest nadir", "nest fwd"], rows, ["%.3f"] * 6))
    return lines


def experiment_azimuth(cache):
    lines = ["### Azimuthal interpolation", "",
             "gortt_kc (and L&S eq. 25) interpolate f between raa = 0 and pi. The alternatives "
             "interpolate PvMv + PiMi - Po, or evaluate the P terms at the actual azimuth. All "
             "with `+ Po = min(PvMv, PiMi)`; RMSE of Kc over off-principal-plane geometries.", ""]
    rows = []
    N = BASE.replace(po="nest")
    for label, var in (("interpolate f (gortt)", N), ("interpolate P", N.replace(az_interp="P")),
                       ("P at actual azimuth", N.replace(az_interp="direct"))):
        row, allerr, ncl = [label], [], 0
        for name, c in CANOPIES.items():
            keys, v, s, ra = select(cache, name, lambda v, s, ra: ra not in (0.0, 180.0))
            K = gk.proportions(c, var, R(v), R(s), R(ra))
            e = K["Kc"] - col(cache, keys, "Kc")
            allerr.append(e)
            ncl += int(np.sum(K["Kt_raw"] < -1e-12))
            row.append(rmse(e))
        rows.append(row + [rmse(np.concatenate(allerr)), ncl])
    lines.append(md_table(["", *CANOPIES, "all", "clamp"], rows, ["%.3f"] * (len(CANOPIES) + 1) + ["%d"]))
    return lines


def experiment_beta(cache):
    lines = ["### Effective beta", "",
             "Kc = beta*Kc(uniform height) + (1-beta)*Kc(random), so the MC implies "
             "beta_eff = (Kc_MC - Kc_random)/(Kc_uniform - Kc_random), with both extremes from "
             "`+ Po = min(PvMv, PiMi)`. Medians over the principal plane by direction; only "
             "geometries where the two extremes differ by more than 0.05 are used.", ""]
    rows = []
    for name in ("SOBS", "dense", "sparse", "wide-spread", "same-height"):
        c = CANOPIES[name]
        for s0 in (15.0, 30.0, 60.0):
            keys, v, s, ra = select(cache, name, lambda v, s, ra, s0=s0: s == s0 and ra in (0.0, 180.0))
            signed = np.where(ra > 90, -v, v)
            mcK = col(cache, keys, "Kc")
            uni = gk.proportions(c, BASE.replace(po="nest", beta=1.0), R(v), R(s), R(ra))["Kc"]
            rnd = gk.proportions(c, BASE.replace(po="nest", beta=0.0), R(v), R(s), R(ra))["Kc"]
            ok = np.abs(uni - rnd) > 0.05
            be = np.where(ok, (mcK - rnd) / np.where(ok, uni - rnd, 1.0), np.nan)
            cb = gk.fFbeta(c, BASE, R(v), R(s), R(ra), gk.kg(c, BASE, R(v), R(s), R(ra)))["beta"]
            code_b = "%.2f-%.2f" % (np.min(cb), np.max(cb))
            dh = c.h2 - c.h1
            eq23 = (1 - dh / (4 * c.b)) ** 2 if dh <= 4 * c.b else 0.0

            def med(m):
                x = be[m & ok]
                return float(np.median(x)) if len(x) else float("nan")
            rows.append(["%s, sza %g" % (name, s0), med(signed > s0), med(np.abs(signed) <= 20),
                         med((signed <= -25) & (signed >= -60)), med(signed <= -70), code_b, eq23])
    lines.append(md_table(["", "back", "nadir", "fwd 25-60", "fwd 70-85", "gortt beta (range)",
                           "eq. 23"], rows, ["%.2f"] * 4 + ["%s", "%.2f"]))
    lines += ["", "gortt's beta depends on sza and, through Gamma_i = Gamma_v, on vza -- but "
              "symmetrically in vza, so it cannot be larger in backscatter than in forward "
              "scatter, which is what the Monte Carlo asks for."]
    return lines


def experiment_P(cache):
    lines = ["### P terms at the parameters of L&S Figs 3-4 (sza 30)", "",
             "View-projected fractions of a crown's projected area. `nest` is min(PvMv, PiMi). "
             "Signed vza' (+ = backscatter).", ""]
    rows = []
    for name, c in LS92.items():
        keys, v, s, ra = select(cache, name)
        vp = np.degrees(gk.prime(c, R(v))) * np.where(ra > 90, -1, 1)
        pick = np.isin(np.round(vp), [-50, -30, -10, 0, 20, 40, 60])
        d = gk.fFbeta(c, gk.GORTT, R(v), R(s), R(ra), gk.kg(c, gk.GORTT, R(v), R(s), R(ra)))
        for k in np.flatnonzero(pick)[np.argsort(vp[pick])]:
            mv = [cache[keys[k]][n] for n in ("PvMv", "PiMi", "Po")]
            rows.append(["%s, vza' %+.0f" % (name[:7], vp[k]), mv[0], d["PvMv"][k], mv[1],
                         d["PiMi"][k], mv[2], min(mv[0], mv[1]), d["Po"][k]])
    lines.append(md_table(["", "PvMv MC", "PvMv gortt", "PiMi MC", "PiMi gortt", "Po MC",
                           "nest(MC)", "Po gortt"], rows, ["%.3f"] * 7))

    # nesting across the whole grid
    po, nest, n = [], [], 0
    for name in CANOPIES:
        keys = select(cache, name)[0]
        po.append(col(cache, keys, "Po"))
        nest.append(np.minimum(col(cache, keys, "PvMv"), col(cache, keys, "PiMi")))
    po, nest = np.concatenate(po), np.concatenate(nest)
    lines += ["", "Over all %d MC geometries of the five canopies: Po / min(PvMv, PiMi) has "
              "median %.3f, 5th percentile %.3f; |Po - min| has RMSE %.4f."
              % (len(po), np.median(po / np.maximum(nest, 1e-9)),
                 np.percentile(po / np.maximum(nest, 1e-9), 5), rmse(po - nest)), ""]
    return lines


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def mc_points(ax, x, y, se, label="Monte Carlo"):
    ax.errorbar(x, y, yerr=2 * se, fmt="o", ms=3.5, color=INK, ecolor=MUTED, elinewidth=0.8,
                zorder=5, label=label)


def plot_fig3(cache, path):
    fig, axes = plt.subplots(3, 3, figsize=(14, 9.5), sharex=True)
    for i, (name, c) in enumerate(LS92.items()):
        keys, v, s, ra = select(cache, name)
        sgn = np.where(ra > 90, -1, 1)
        vp = np.degrees(gk.prime(c, R(v))) * sgn
        o = np.argsort(vp)
        keys, v, s, ra, vp = [keys[k] for k in o], v[o], s[o], ra[o], vp[o]
        fine = np.linspace(-50, 60, 441)
        fv = np.degrees(np.arctan(np.tan(R(np.abs(fine))) * c.r / c.b))
        fra = np.where(fine < 0, np.pi, 0.0)
        Kg = gk.kg(c, gk.GORTT, R(fv), R(30.0) + 0 * fv, fra)
        code = gk.fFbeta(c, gk.GORTT, R(fv), R(30.0) + 0 * fv, fra, Kg)
        lune = gk.fFbeta(c, gk.GORTT.replace(lune=True), R(fv), R(30.0) + 0 * fv, fra, Kg)
        pv_p, pi_p = paper_P(c, R(fv), R(30.0) + 0 * fv, fra)
        for j, (P, PM, M) in enumerate((("Pv", "PvMv", "Mv"), ("Pi", "PiMi", "Mi"))):
            ax = axes[i, j]
            style(ax)
            m = col(cache, keys, M)
            mc_points(ax, vp, col(cache, keys, PM) / m, col(cache, keys, PM + "_se") / m)
            ax.plot(fine, code[PM] / code[M], color=SERIES[0], lw=2, label="gortt (IGARSS forms)")
            ax.plot(fine, (pv_p, pi_p)[j], color=SERIES[1], lw=2, ls="--",
                    label="L&S eqs 19-22 by case")
            ax.plot(fine, lune[PM] / lune[M], color=SERIES[2], lw=2, ls=":", label="lune model")
            ax.set_ylim(-0.25, 1.25)
            ax.set_title("%s   %s" % (P, name), color=INK, fontsize=9)
        # the only combination f depends on
        ax = axes[i, 2]
        style(ax)
        dP = lambda d: d["PvMv"] + d["PiMi"] - d["Po"]
        mdp = col(cache, keys, "PvMv") + col(cache, keys, "PiMi") - col(cache, keys, "Po")
        mse = np.sqrt(col(cache, keys, "PvMv_se") ** 2 + col(cache, keys, "PiMi_se") ** 2
                      + col(cache, keys, "Po_se") ** 2)
        mc_points(ax, vp, mdp, mse)
        nest = gk.fFbeta(c, gk.GORTT.replace(po="nest"), R(fv), R(30.0) + 0 * fv, fra, Kg)
        ax.plot(fine, dP(code), color=SERIES[0], lw=2, label="gortt")
        ax.plot(fine, np.where(fine >= np.degrees(gk.prime(c, R(30.0))),
                               code["PvMv"], pi_p * code["Mi"]),
                color=SERIES[1], lw=2, ls="--", label="L&S by case")
        ax.plot(fine, dP(lune), color=SERIES[2], lw=2, ls=":", label="lune model")
        ax.plot(fine, dP(nest), color=SERIES[3], lw=2, ls="-.", label="Po = min(PvMv, PiMi)")
        ax.set_title("ΔP = PvMv + PiMi − Po   %s" % name, color=INK, fontsize=9)
        ax.legend(fontsize=7, frameon=False, loc="upper left") if i == 0 else None
    for j in range(3):
        axes[2, j].set_xlabel("signed view zenith θv′ (°), + = backscatter", color=MUTED, fontsize=8)
    axes[0, 0].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.suptitle("Mutual-shadowing conditional probabilities, L&S (1992) Fig. 3/4 parameters, "
                 "θi = 30°", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_nesting(cache, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, name in zip(axes, ("SOBS", "dense", "wide-spread")):
        keys, v, s, ra = select(cache, name, lambda v, s, ra: s == 30.0 and ra in (0.0, 180.0))
        x = np.where(ra > 90, -v, v)
        o = np.argsort(x)
        keys, x = [keys[k] for k in o], x[o]
        style(ax)
        for k, (n, lab) in enumerate((("PvMv", "PvMv"), ("PiMi", "PiMi"))):
            ax.plot(x, col(cache, keys, n), color=SERIES[k], lw=2, label="MC " + lab)
        ax.plot(x, np.minimum(col(cache, keys, "PvMv"), col(cache, keys, "PiMi")), color=MUTED,
                lw=5, alpha=0.35, label="min(PvMv, PiMi)")
        mc_points(ax, x, col(cache, keys, "Po"), col(cache, keys, "Po_se"), label="MC Po")
        ax.set_title("%s, sza 30°" % name, color=INK, fontsize=10)
        ax.set_xlabel("signed vza (°), + = backscatter", color=MUTED, fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.suptitle("Po is the intersection of the two mutually-shadowed sets: it follows the smaller one",
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_kg(cache, path):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    c = CANOPIES["SOBS"]
    for j, s0 in enumerate((15.0, 30.0, 60.0)):
        keys, v, s, ra = select(cache, "SOBS", lambda v, s, ra, s0=s0: s == s0 and ra in (0.0, 180.0))
        x = np.where(ra > 90, -v, v)
        o = np.argsort(x)
        keys, v, s, ra, x = [keys[k] for k in o], v[o], s[o], ra[o], x[o]
        fine = np.linspace(-85, 85, 681)
        fv, fs, fra = gk.signed_pp(fine, s0)
        ax = axes[0, j]
        style(ax)
        mc_points(ax, x, col(cache, keys, "Kg"), col(cache, keys, "Kg_se"))
        for k, h in enumerate(("depth", "centre", "avg")):
            ax.plot(fine, gk.kg(c, gk.GORTT.replace(overlap_h=h), fv, fs, fra), color=SERIES[k],
                    lw=2, ls=("-", "--", ":")[k], label="h = " + h)
        ax.set_title("SOBS principal plane, sza %g°" % s0, color=INK, fontsize=9)
        ax.set_xlabel("signed vza (°)", color=MUTED, fontsize=8)

        keys, v, s, ra = select(cache, "SOBS", lambda v, s, ra, s0=s0: s == s0 and v == s0)
        o = np.argsort(ra)
        keys, ra = [keys[k] for k in o], ra[o]
        fra = np.linspace(0, 180, 361)
        ax = axes[1, j]
        style(ax)
        mc_points(ax, ra, col(cache, keys, "Kg"), col(cache, keys, "Kg_se"))
        for k, (h, az) in enumerate((("depth", "ambrals"), ("avg", "ambrals"), ("avg", "eq7"), ("avg", "wenge"))):
            ax.plot(fra, gk.kg(c, gk.GORTT.replace(overlap_h=h, overlap_az=az), R(s0), R(s0), R(fra)),
                    color=SERIES[k], lw=2, ls=("-", "--", ":", "-.")[k], label="h=%s, %s" % (h, az))
        ax.set_title("SOBS principal cone vza = sza = %g°" % s0, color=INK, fontsize=9)
        ax.set_xlabel("relative azimuth (°)", color=MUTED, fontsize=8)
    axes[0, 0].legend(fontsize=7.5, frameon=False)
    axes[1, 0].legend(fontsize=7.5, frameon=False)
    axes[0, 0].set_ylabel("Kg", color=INK)
    axes[1, 0].set_ylabel("Kg", color=INK)
    fig.suptitle("Sunlit viewed background: overlap-function choices against Monte Carlo",
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_kc(cache, path):
    show = [KC_VARIANTS[0], KC_VARIANTS[1], KC_VARIANTS[3], KC_VARIANTS[5]]
    names = ("SOBS", "same-height", "wide-spread")
    fig, axes = plt.subplots(len(names), 4, figsize=(15, 10), sharex=True)
    for i, name in enumerate(names):
        c = CANOPIES[name]
        for j, s0 in enumerate(SZAS):
            keys, v, s, ra = select(cache, name, lambda v, s, ra, s0=s0: s == s0 and ra in (0.0, 180.0))
            x = np.where(ra > 90, -v, v)
            o = np.argsort(x)
            keys, x = [keys[k] for k in o], x[o]
            fine = np.linspace(-85, 85, 681)
            fv, fs, fra = gk.signed_pp(fine, s0)
            ax = axes[i, j]
            style(ax)
            mc_points(ax, x, col(cache, keys, "Kc"), col(cache, keys, "Kc_se"))
            for k, (label, var) in enumerate(show):
                ax.plot(fine, gk.proportions(c, var, fv, fs, fra)["Kc"], color=SERIES[k], lw=2,
                        ls=("-", "--", "-", ":")[k], label=label)
            ax.set_ylim(-0.02, 1.25)
            ax.set_title("%s, sza %g°" % (name, s0), color=INK, fontsize=9)
            if i == len(names) - 1:
                ax.set_xlabel("signed vza (°), + = backscatter", color=MUTED, fontsize=8)
    axes[0, 0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("Sunlit viewed crown Kc in the principal plane: variants against Monte Carlo",
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# Ni et al. (1999) Fig. 5, read off the page (sza 30): see figures/README.md
FIG5_PUBLISHED = [(30.0, 0.81), (44.0, 0.79), (60.0, 0.80), (0.0, 0.15), (-40.0, 0.13), (-60.0, 0.15)]


def plot_fig5(cache, path):
    c = CANOPIES["SOBS"]
    fig, axes = plt.subplots(2, 1, figsize=(9, 8.5), sharex=True)
    for ax, s0 in zip(axes, (30.0, 60.0)):
        keys, v, s, ra = select(cache, "SOBS", lambda v, s, ra, s0=s0: s == s0 and ra in (0.0, 180.0))
        x = np.where(ra > 90, -v, v)
        o = np.argsort(x)
        keys, x = [keys[k] for k in o], x[o]
        fine = np.linspace(-89, 89, 713)
        fv, fs, fra = gk.signed_pp(fine, s0)
        style(ax)
        mc_points(ax, x, col(cache, keys, "Kc"), col(cache, keys, "Kc_se"), label="Kc Monte Carlo")
        for k, (label, var) in enumerate((KC_VARIANTS[0], KC_VARIANTS[1], KC_VARIANTS[3])):
            ax.plot(fine, gk.proportions(c, var, fv, fs, fra)["Kc"], color=SERIES[k], lw=2,
                    ls=("-", "--", "-")[k], label="Kc " + label)
        if s0 == 30.0:
            px, py = zip(*FIG5_PUBLISHED)
            ax.plot(px, py, "D", ms=7, mfc="none", mec=SERIES[4], mew=2,
                    label="Ni et al. (1999) Fig. 5, read off")
        ax.set_ylim(-0.02, 1.25)
        ax.set_ylabel("Kc", color=INK)
        ax.set_title("SOBS, sza %g°" % s0, color=INK, fontsize=10)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    axes[1].set_xlabel("signed vza (°), + = backscatter", color=MUTED, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(HERE, "kc_investigation"))
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    cache = load_mc(os.path.join(a.outdir, "mc_cache.csv"), a.jobs, a.refresh)

    lines = ["<!-- generated by figures/kc_investigation.py; MC: %d rays, %d elements, %d tiles "
             "per geometry -->" % (MC_RAYS, MC_ELEM, MC_SCENES), ""]
    for exp in (experiment_P, experiment_kg, experiment_kc, experiment_azimuth, experiment_beta):
        lines += exp(cache) + [""]
    with open(os.path.join(a.outdir, "tables.md"), "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))

    plot_fig3(cache, os.path.join(a.outdir, "fig3_P_terms.png"))
    plot_nesting(cache, os.path.join(a.outdir, "nesting.png"))
    plot_kg(cache, os.path.join(a.outdir, "kg_overlap.png"))
    plot_kc(cache, os.path.join(a.outdir, "kc_variants.png"))
    plot_fig5(cache, os.path.join(a.outdir, "fig5.png"))
    print("wrote", a.outdir)


if __name__ == "__main__":
    main()
