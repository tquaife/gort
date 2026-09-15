#!/usr/bin/env python3
"""
A numpy port of gortt's geometric-optical kernels, with every contested choice
behind a switch.

This is an investigation tool, not a second implementation of the model.  With
the switches at their defaults it reproduces what gortt prints for -prnprop
(Kc, Kg, Kt, Kz) -- "--selftest" checks that against a built binary.
GORTT_129ED80 and GORTT_MASTER reproduce the two earlier versions the
investigation started from, and the other settings are the alternatives found
while checking gortt_brdf.c against
Li & Strahler (1992), IEEE TGRS 30(2) 276-292 ("L&S" below).  See
KC-INVESTIGATION.md for what each one is and why it is here.

Everything is in radians and works elementwise on arrays.  Angles are gortt's
own convention: vza and sza non-negative, raa folded into [0, pi], and raa = 0
is backscatter (the viewer on the sun's side).

The switches (a Variant):

  fsq         "new"   f = (1-frac)*f0 + frac*f180            (L&S eq. 25), gortt
              "old"   f = (1-frac)*f0*F0 + frac*f180*F180    (F squared, master)

  overlap_h   "avg"    the overlap averaged over centre heights in [h1, h2],
                       which is what the O-bar in L&S eq. 3 asks for; gortt
              "centre" h = (h1 + h2)/2, the crown-centre height L&S eqs 4-6 use
              "depth"  h = 2b + h2 - h1, the canopy depth (gortt up to 59121b9)

  overlap_az  "ambrals" t2 = sqrt(D^2 + (tan i' tan v' sin raa)^2), gortt
              "pp"      t2 = |tan i' - tan v' cos raa|, L&S eq. 6 (removed from gortt)
              "eq7"     L&S eq. 7: linear in raa from O(0) to O(pi) at the cutoff
                        Phi = 4R/(h(tan v + tan i))
              "wenge"   a t1 "from Wenge's code" (removed from gortt)

  pimi        "igarss"   (1 - cos(tMi*(1 - (i' - v' cos raa)/pi)))/2, gortt
              "eq21"     (1 - cos(tMi - i' + v' cos raa))/2, clamped to 0 when the
                         argument is negative, as the L&S text says
              "eq21raw"  the same without the clamp, as the removed
                         gortt_kc_PP_only coded it

  pvmv        "igarss"   Mv - (1 - cos(v' cos raa - i'))/2, gortt
              "clamped"  the same, floored at 0
              "eq22"     (1 - cos(tMv - i' + v' cos raa))/2 or 0, L&S eq. 22 as printed
              "ppcoded"  the removed gortt_kc_PP_only: condition on
                         tMv - i' + v' cos raa, value from tMv - v' + v' cos raa

  po          "nest"     min(PvMv, PiMi): the smaller mutually-shadowed set lies
                         inside the larger (what the Monte Carlo shows); gortt
              "code"     PvMv in forward scatter; PiMi if vza > sza; else PvMv
                         (gortt up to 0204222)
              "zerofwd"  as "code" but 0 in forward scatter

  lune        False, or True: ignore pimi/pvmv/po and compute all three from the
              exact interval geometry on the primed sphere in the principal
              plane (KC-INVESTIGATION.md, "the lune model")

  az_interp   "f"        interpolate f between raa = 0 and pi, as gortt_kc does
              "P"        interpolate PvMv + PiMi - Po instead, and build f from
                         F, M and Gamma at the actual azimuth (L&S sec. H
                         suggests interpolating the P's)
              "direct"   no interpolation: the P terms at the actual azimuth

  beta        "code"     gortt's (IGARSS '92) formula with D = r cot(i'/2)
              "eq23"     (1 - (h2-h1)/4b)^2, 0 if h2-h1 > 4b   (L&S eq. 23)
              a float    fixed, like gortt -beta
"""

import argparse
import dataclasses
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gortt_run as gr

HALF_PI = 0.5 * np.pi


@dataclasses.dataclass(frozen=True)
class Canopy:
    r: float = 0.76
    b: float = 2.7
    lam: float = 0.41
    h1: float = 3.0
    h2: float = 8.5

    @property
    def rr(self):
        return self.r * self.r

    def args(self):
        """The gortt options that describe this canopy."""
        return ["-r", repr(self.r), "-b", repr(self.b), "-lambda", repr(self.lam),
                "-h1", repr(self.h1), "-h2", repr(self.h2)]


SOBS = Canopy()


@dataclasses.dataclass(frozen=True)
class Variant:
    fsq: str = "new"
    overlap_h: str = "avg"
    overlap_az: str = "ambrals"
    pimi: str = "igarss"
    pvmv: str = "igarss"
    po: str = "nest"
    lune: bool = False
    beta: object = "code"
    az_interp: str = "f"

    def replace(self, **kw):
        return dataclasses.replace(self, **kw)

    def label(self):
        base = Variant()
        diff = ["%s=%s" % (f.name, getattr(self, f.name))
                for f in dataclasses.fields(self)
                if getattr(self, f.name) != getattr(base, f.name)]
        return ", ".join(diff) if diff else "gortt"


GORTT = Variant()
# the two versions the investigation started from
GORTT_129ED80 = Variant(po="code", overlap_h="depth")
GORTT_MASTER = GORTT_129ED80.replace(fsq="old")


def sec(x):
    return 1.0 / np.cos(x)


def prime(c, za):
    return np.arctan((c.b / c.r) * np.tan(za))


# ---------------------------------------------------------------------------
# overlap and Kg
# ---------------------------------------------------------------------------

def _overlap_at_h(c, v, ip, vp, raa, hb):
    """gortt_overlap at a given h/b, for the ambrals/pp/wenge azimuth forms."""
    ti, tv = np.tan(ip), np.tan(vp)
    d = ti * ti + tv * tv - 2.0 * ti * tv * np.cos(raa)
    D = np.sqrt(np.maximum(0.0, d))

    if v.overlap_az == "pp":
        t2 = np.abs(ti - tv * np.cos(raa))
    else:
        t2 = np.sqrt(D * D + (ti * tv * np.sin(raa)) ** 2)

    if v.overlap_az == "wenge":
        with np.errstate(divide="ignore", invalid="ignore"):
            a = (ti * np.sin(raa) / D) ** 2 + ((ti * np.cos(raa) - tv) / D) ** 2 * sec(ip) ** 2
            bb = (tv * np.sin(raa) / D) ** 2 + ((tv * np.cos(raa) - ti) / D) ** 2 * sec(vp) ** 2
        t1 = np.sqrt(np.maximum(0.0, a)) + np.sqrt(np.maximum(0.0, bb))
    else:
        t1 = sec(ip) + sec(vp)

    with np.errstate(divide="ignore", invalid="ignore"):
        cos_t = hb * t2 / t1
    cos_t = np.clip(cos_t, -1.0, 1.0)
    t = np.arccos(cos_t)
    return np.maximum(0.0, (t - np.sin(t) * cos_t) * (sec(ip) + sec(vp)) / np.pi)


def _overlap_h(c, v, ip, vp, raa, h, vza, sza):
    if v.overlap_az == "eq7":
        pp = v.replace(overlap_az="pp")
        o0 = _overlap_at_h(c, pp, ip, vp, 0.0 * raa, h / c.b)
        o180 = _overlap_at_h(c, pp, ip, vp, np.pi + 0.0 * raa, h / c.b)
        with np.errstate(divide="ignore"):
            phi = 4.0 * c.r / (h * (np.tan(vza) + np.tan(sza)))
        phi = np.minimum(phi, np.pi)
        w = np.clip(raa / phi, 0.0, 1.0)
        return o0 + w * (o180 - o0)
    return _overlap_at_h(c, v, ip, vp, raa, h / c.b)


NH_OVERLAP = 24


def gauleg(n, eps=3.0e-11):
    """gauleg() from gortt_albedo.c on [-1, 1], iteration for iteration, so the
    overlap average uses the same nodes to the last bit or two."""
    x, w = np.zeros(n), np.zeros(n)
    for i in range((n + 1) // 2):
        z = np.cos(3.141592654 * (i + 0.75) / (n + 0.5))
        while True:
            p1, p2 = 1.0, 0.0
            for j in range(1, n + 1):
                p3, p2 = p2, p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) <= eps:
                break
        x[i], x[n - 1 - i] = -z, z
        w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp)
    return x, w


def overlap(c, v, vza, sza, raa):
    """gortt_overlap, with the h and azimuth switches."""
    ip, vp = prime(c, sza), prime(c, vza)
    if v.overlap_h == "depth":
        return _overlap_h(c, v, ip, vp, raa, 2.0 * c.r * (c.b / c.r) + c.h2 - c.h1, vza, sza)
    if v.overlap_h == "centre":
        return _overlap_h(c, v, ip, vp, raa, 0.5 * (c.h1 + c.h2), vza, sza)
    if v.overlap_h == "avg":
        x, w = gauleg(NH_OVERLAP)
        hs = 0.5 * (c.h2 - c.h1) * x + 0.5 * (c.h1 + c.h2)
        return sum(0.5 * wk * _overlap_h(c, v, ip, vp, raa, hk, vza, sza)
                   for hk, wk in zip(hs, w))
    raise ValueError(v.overlap_h)


def kg(c, v, vza, sza, raa):
    ip, vp = prime(c, sza), prime(c, vza)
    o = overlap(c, v, vza, sza, raa)
    return np.exp(-(c.lam * c.r ** 2 * np.pi * (sec(ip) + sec(vp) - o)))


# ---------------------------------------------------------------------------
# the f-ratio
# ---------------------------------------------------------------------------

def mutual(c, za_p):
    """Mi or Mv, L&S eqs 10-11."""
    x = c.lam * np.pi * c.rr * sec(za_p)
    return 1.0 - (1.0 - np.exp(-x)) / x


def _strip(a, b, u):
    """View-projected fraction of the lune psi in (a, b), clipped to the visible half."""
    xa = np.clip(a - u, -HALF_PI, HALF_PI)
    xb = np.clip(b - u, -HALF_PI, HALF_PI)
    return np.where(xb > xa, 0.5 * (np.sin(xb) - np.sin(xa)), 0.0)


def lune_terms(ip, vp, raa, tMi, tMv):
    """
    PvMv, PiMi, Po from the interval geometry of the principal plane.

    psi is the angle of a surface normal from the zenith, in the principal
    plane, positive towards the sun.  u is the signed view direction.  Each
    set below is a lune -- all normals whose psi lies in an interval -- and
    the view-projected area of a lune is (sin(x_b) - sin(x_a))/2 as a fraction
    of the crown's projected disc, with x = psi - u clipped to the visible
    half.  Off the principal plane this uses cos(raa) to place the viewer,
    which is only exact at 0 and pi -- the only places gortt_kc evaluates f.
    """
    u = vp * np.sign(np.cos(raa) + 1e-300)
    sun_lo, sun_hi = ip - HALF_PI, ip + HALF_PI
    mi_lo, mi_hi = ip + HALF_PI - tMi, ip + HALF_PI
    # lower limb on the viewer's side; nadir goes with the sun side
    mv_lo = np.where(u >= 0.0, u + HALF_PI - tMv, u - HALF_PI)
    mv_hi = np.where(u >= 0.0, u + HALF_PI, u - HALF_PI + tMv)

    pvmv = _strip(np.maximum(mv_lo, sun_lo), np.minimum(mv_hi, sun_hi), u)
    pimi = _strip(mi_lo, mi_hi, u)
    po = _strip(np.maximum(mi_lo, mv_lo), np.minimum(mi_hi, mv_hi), u)
    return pvmv, pimi, po


def fFbeta(c, v, vza, sza, raa, Kg):
    """gortt_kc_fFbeta: f, F and beta at one relative azimuth."""
    ip, vp = prime(c, sza), prime(c, vza)
    o = overlap(c, v, vza, sza, raa)
    phase = np.cos(vp) * np.cos(ip) + np.sin(vp) * np.sin(ip) * np.cos(raa)

    Mi, Mv = mutual(c, ip), mutual(c, vp)
    Gamma = np.pi * c.rr * (sec(ip) + sec(vp) - o)
    Gamma_c = np.pi * c.rr * sec(vp) * 0.5 * (1.0 + phase)
    Gamma_v = np.pi * c.rr * sec(vp)
    F = Gamma_c / Gamma
    M = 1.0 - (1.0 - Kg) / (c.lam * Gamma)

    tMi = np.arccos(1.0 - 2.0 * Mi)
    tMv = np.arccos(1.0 - 2.0 * Mv)

    if v.lune:
        PvMv, PiMi, Po = lune_terms(ip, vp, raa, tMi, tMv)
    else:
        if v.pimi == "igarss":
            PiMi = (1 - np.cos(tMi * (1 - (ip - vp * np.cos(raa)) / np.pi))) / 2.0
        elif v.pimi in ("eq21", "eq21raw"):
            arg = tMi - ip + vp * np.cos(raa)
            PiMi = (1.0 - np.cos(arg)) / 2.0
            if v.pimi == "eq21":
                PiMi = np.where(arg > 0.0, PiMi, 0.0)
        else:
            raise ValueError(v.pimi)

        if v.pvmv in ("igarss", "clamped"):
            PvMv = Mv - (1.0 - np.cos(vp * np.cos(raa) - ip)) / 2.0
            if v.pvmv == "clamped":
                PvMv = np.maximum(0.0, PvMv)
        elif v.pvmv in ("eq22", "ppcoded"):
            cond = tMv - ip + vp * np.cos(raa)
            val = cond if v.pvmv == "eq22" else tMv - vp + vp * np.cos(raa)
            PvMv = np.where(cond >= 0.0, (1.0 - np.cos(val)) / 2.0, 0.0)
        else:
            raise ValueError(v.pvmv)

        forward = (raa < np.radians(270.0)) & (raa > np.radians(90.0))
        Po = np.where(forward, PvMv, np.where(np.abs(vza) > np.abs(sza), PiMi, PvMv))
        if v.po == "zerofwd":
            Po = np.where(forward, 0.0, Po)
        elif v.po == "nest":
            Po = np.minimum(PvMv, PiMi)
        elif v.po != "code":
            raise ValueError(v.po)

    if isinstance(v.beta, str) and v.beta == "code":
        with np.errstate(divide="ignore", invalid="ignore"):
            D = c.r * (1.0 / np.tan(ip / 2.0))
            Gi = Gamma_v
            dh = c.h2 - c.h1
            beta = (c.lam * Gi) / (c.lam * Gi + dh / D) \
                * (1.0 - np.exp(-c.lam * Gi - dh / D)) / (1.0 - np.exp(-c.lam * Gi))
        beta = np.where(ip < 0.000000001, 0.0, beta)
    elif isinstance(v.beta, str) and v.beta == "eq23":
        dh = c.h2 - c.h1
        beta = np.full_like(ip, (1.0 - dh / (4.0 * c.b)) ** 2 if dh <= 4.0 * c.b else 0.0)
    else:
        beta = np.full_like(ip, float(v.beta))

    f = F * (1.0 - Gamma_v * (PvMv + PiMi - Po) / Gamma_c) / (1.0 - M)
    return dict(f=f, F=F, beta=beta, PvMv=PvMv, PiMi=PiMi, Po=Po, Mi=Mi, Mv=Mv,
                Gamma=Gamma, Gamma_c=Gamma_c, Gamma_v=Gamma_v, M=M)


def kc(c, v, vza, sza, raa, Kg):
    """gortt_kc."""
    at = fFbeta(c, v, vza, sza, raa, Kg)
    z = 0.0 * raa
    Kg0 = kg(c, v, vza, sza, z)
    p0 = fFbeta(c, v, vza, sza, z, Kg0)
    Kg180 = kg(c, v, vza, sza, z + np.pi)
    p180 = fFbeta(c, v, vza, sza, z + np.pi, Kg180)

    frac = raa / np.pi
    frac = np.where(frac > 1.0, 2.0 - frac, frac)

    if v.az_interp in ("P", "direct"):
        if v.fsq != "new":
            raise ValueError("az_interp=%s needs fsq=new" % v.az_interp)
        dP = lambda d: d["PvMv"] + d["PiMi"] - d["Po"]
        dPx = dP(at) if v.az_interp == "direct" else (1. - frac) * dP(p0) + frac * dP(p180)
        f = at["F"] * (1.0 - at["Gamma_v"] * dPx / at["Gamma_c"]) / (1.0 - at["M"])
    elif v.fsq == "new":
        f = (1. - frac) * p0["f"] + frac * p180["f"]
    elif v.fsq == "old":
        f = (1. - frac) * p0["f"] * p0["F"] + frac * p180["f"] * p180["F"]
    else:
        raise ValueError(v.fsq)

    beta = at["beta"]
    f = beta * f + (1.0 - beta) * at["F"]
    return f * (1.0 - Kg)


def proportions(c, v, vza, sza, raa):
    """(Kc, Kg, Kt, Kz) as gortt_rsurf computes them, plus the unclamped Kt."""
    vza, sza, raa = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in (vza, sza, raa)))
    Kg = kg(c, v, vza, sza, raa)
    Kc = kc(c, v, vza, sza, raa, Kg)
    vp = prime(c, vza)
    Kz = np.exp(-(c.lam * np.pi * c.r ** 2) / np.cos(vp)) - Kg
    Kt_raw = 1.0 - Kc - Kz - Kg
    return dict(Kc=Kc, Kg=Kg, Kt=np.maximum(0.0, Kt_raw), Kz=Kz, Kt_raw=Kt_raw)


def signed_pp(vza_signed, sza):
    """Signed principal-plane vza (+ = backscatter) -> (vza, sza, raa) in radians."""
    v = np.radians(np.asarray(vza_signed, dtype=float))
    return np.abs(v), np.radians(sza) + 0.0 * v, np.where(v < 0.0, np.pi, 0.0)


# ---------------------------------------------------------------------------
# self-test against gortt
# ---------------------------------------------------------------------------

def selftest(binary, variant, tol=1e-12):
    canopies = [SOBS, dataclasses.replace(SOBS, lam=0.1), dataclasses.replace(SOBS, lam=1.0),
                Canopy(r=1.0, b=3.0, h1=2.0, h2=9.0, lam=0.3)]
    betas = ["code", 1.0, 0.0]
    worst = 0.0
    n = 0
    for c in canopies:
        for beta in betas:
            rows = []
            for s in (0.0, 15.0, 30.0, 45.0, 60.0):
                rows += gr.principal_plane(np.arange(-88.0, 88.1, 2.0), s)
                rows += [(vz, float(ra), s, 0.0) for vz in (5.0, 20.0, 45.0, 70.0)
                         for ra in np.arange(0.0, 360.0, 15.0)]
            extra = c.args() + (["-beta", repr(beta)] if beta != "code" else [])
            # gortt_run.run prepends SOBS; ours come after, so they win
            got = gr.run(binary, rows, 800.0, extra, want="prop")

            vza = np.radians([r[0] for r in rows])
            sza = np.radians([r[2] for r in rows])
            raa = np.radians([(r[1] - r[3]) % 360.0 for r in rows])
            raa = np.where(raa > np.pi, 2.0 * np.pi - raa, raa)
            K = proportions(c, variant.replace(beta=beta), vza, sza, raa)
            mine = np.stack([K["Kc"], K["Kg"], K["Kt"], K["Kz"]], axis=1)
            err = np.max(np.abs(mine - got))
            worst = max(worst, err)
            n += len(rows)
            if not err <= tol:
                i = np.unravel_index(np.argmax(np.abs(mine - got)), mine.shape)
                print("MISMATCH %s beta=%s row %r K%s: gortt %r port %r"
                      % (c, beta, rows[i[0]], "cgtz"[i[1]], got[i], mine[i]))
    return worst, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true",
                    help="compare the default Variant with gortt -prnprop")
    ap.add_argument("--new", default=gr.default_binary(), help="current gortt")
    ap.add_argument("--at129ed80", help="gortt built at 129ed80 (optional)")
    ap.add_argument("--master", help="gortt built at e109ff0 / pre-kc-fsquared (optional)")
    a = ap.parse_args()
    if not a.selftest:
        ap.print_help()
        return
    ok = True
    for name, var, b in (("current", GORTT, a.new), ("129ed80", GORTT_129ED80, a.at129ed80),
                         ("master", GORTT_MASTER, a.master)):
        if b is None:
            continue
        gr.check_binary(b)
        worst, n = selftest(b, var)
        print("selftest %s: %d geometries x 4 proportions, max |port - gortt| = %.3g"
              % (name, n, worst))
        ok = ok and worst <= 1e-12
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
