#!/usr/bin/env python3
"""
Monte Carlo reference for the geometric-optical proportions Kc, Kg, Kt, Kz and
for the mutual-shadowing terms of Li & Strahler (1992) eqs 10-17.

The scene is the one GOMS assumes: solid spheroidal crowns of horizontal
radius r and vertical half-axis b, centres uniformly random on the ground at
density lambda (a fixed count per periodic tile, i.e. a binomial field -- the
Poisson field of the model in the limit of a large tile), centre heights
uniform in [h1, h2].  (Within-crown gaps are a
matter for the signatures in gortt, not for these proportions.)

Everything is done in "sphere space": scaling z by r/b turns every crown into
a sphere of radius r and every zenith angle theta into theta' = atan((b/r)
tan theta).  Ground-projected areas are unchanged by a vertical scaling,
straight lines stay straight, and the sign of n.d for a surface normal n and
direction d is preserved, so being sunlit, being viewed, and every projected
area fraction are the same in both spaces.

Two estimators per geometry:

  rays       Parallel view rays through uniformly sampled ground points.  The
             first crown hit (or the ground) is classified sunlit if it faces
             the sun (crowns only) and the sun ray from it clears every other
             crown.  Gives Kc, Kt, Kg, Kz directly.

  elements   Uniform samples over the view-projected disc of randomly chosen
             crowns, mapped onto the crown's viewer-facing surface -- i.e.
             crown surface elements weighted by their area projected in the
             view direction, as L&S eq. 16 measures them.  Tallies
               Gc/Gv  facing the sun
               Mv     hidden from the viewer by other crowns
               PvMv   facing the sun and hidden from the viewer
               PiMi   facing the sun and shadowed from the sun by other crowns
               Po     all three
             and, from a separate sample over the sun-projected disc, Mi.  Then
             Kc = lambda*pi*r^2*sec(v') * (Gc/Gv - PvMv - PiMi + Po) exactly,
             which ties the P terms to the rays estimate.

  usage:  ./goms_mc.py [--vza 30 --sza 30 --raa 0] [--rays N] [--elements N]
                       [-r .. -b .. --lambda .. --h1 .. --h2 ..] [--seed S]
"""

import argparse
import dataclasses
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import goms_kernels as gk

EPS = 1e-9


def unit(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


class Scene:
    """A periodic tile of spheres, with enough image copies to trace from it."""

    def __init__(self, c, seed=0, n_target=1500, extent=0.0):
        rng = np.random.default_rng(seed)
        self.c = c
        self.k = c.r / c.b
        self.L = np.sqrt(n_target / c.lam)
        # a fixed count (a binomial rather than Poisson field): a Poisson count
        # makes each tile's effective lambda wander by 1/sqrt(N), which moved Kg
        # by far more than its sampling error
        N = int(round(c.lam * self.L ** 2))
        self.L = np.sqrt(N / c.lam)
        xy = rng.uniform(0.0, self.L, (N, 2))
        z = rng.uniform(c.h1 * self.k, c.h2 * self.k, N)
        self.base = np.column_stack([xy, z])
        self.zmax = c.h2 * self.k + c.r
        self.rng = rng
        self._images(extent)

    def _images(self, extent):
        n = int(np.ceil((extent + 2.0 * self.c.r) / self.L))
        shifts = np.array([(i, j) for i in range(-n, n + 1) for j in range(-n, n + 1)])
        N = len(self.base)
        self.centres = np.concatenate([self.base + [sx * self.L, sy * self.L, 0.0]
                                       for sx, sy in shifts])
        # the same crown in every copy shares an identity, for self-exclusion
        self.ident = np.tile(np.arange(N), len(shifts))
        self.central = np.flatnonzero((shifts[:, 0] == 0) & (shifts[:, 1] == 0))[0] * N \
            + np.arange(N)
        self._trees = {}

    def _frame(self, d):
        """Map 3-D points to 2-D so that a line along d hits a sphere iff within r."""
        dxy = np.hypot(d[0], d[1])
        e1 = np.array([d[0], d[1]]) / dxy if dxy > 0 else np.array([1.0, 0.0])
        e2 = np.array([-e1[1], e1[0]])

        def to2d(P):
            proj = P[:, :2] - (P[:, 2:3] / d[2]) * d[:2]
            return np.column_stack([(proj @ e1) * d[2], proj @ e2])
        return to2d

    def tree(self, d):
        key = tuple(np.round(d, 15))
        if key not in self._trees:
            to2d = self._frame(d)
            self._trees[key] = (cKDTree(to2d(self.centres)), to2d)
        return self._trees[key]

    def first_hits(self, P, d, exclude=None):
        """
        For points P and direction d: the nearest-to-source crown each line
        P + t d (t > 0) enters, as (hit mask, sphere index, t).  With exclude,
        ignores the crown whose identity matches (for rays leaving a surface).
        Returns the largest t, i.e. the first surface met coming in from +d.
        """
        tree, to2d = self.tree(d)
        pts = cKDTree(to2d(P))
        pairs = pts.sparse_distance_matrix(tree, self.c.r, output_type="ndarray")
        i = pairs["i"].astype(np.int64)
        j = pairs["j"].astype(np.int64)
        w = P[i] - self.centres[j]
        b = w @ d
        cc = np.einsum("ij,ij->i", w, w) - self.c.r ** 2
        disc = b * b - cc
        tp = -b + np.sqrt(np.maximum(disc, 0.0))
        ok = (disc > 0.0) & (tp > EPS * self.c.r)
        if exclude is not None:
            ok &= self.ident[j] != self.ident[exclude[i]]
        i, j, tp = i[ok], j[ok], tp[ok]

        hit = np.zeros(len(P), bool)
        best = np.full(len(P), -1, np.int64)
        tbest = np.full(len(P), -np.inf)
        if len(i):
            order = np.lexsort((tp, i))
            i, j, tp = i[order], j[order], tp[order]
            last = np.r_[i[1:] != i[:-1], True]
            hit[i[last]] = True
            best[i[last]] = j[last]
            tbest[i[last]] = tp[last]
        return hit, best, tbest

    def blocked(self, P, d, exclude=None):
        return self.first_hits(P, d, exclude)[0]


def extent_for(c, vza, sza):
    k = c.r / c.b
    zmax = c.h2 * k + c.r
    return zmax * (np.tan(gk.prime(c, vza)) + np.tan(gk.prime(c, sza)))


def _front_samples(scene, d, n, rng):
    """Uniform over the d-projected discs of central crowns, lifted to the d-facing surface."""
    r = scene.c.r
    j = scene.central[rng.integers(0, len(scene.central), n)]
    a = np.cross(d, [0.0, 0.0, 1.0]) if abs(d[2]) < 1.0 - 1e-12 else np.array([1.0, 0.0, 0.0])
    e1 = a / np.linalg.norm(a)
    e2 = np.cross(d, e1)
    rho = r * np.sqrt(rng.uniform(0.0, 1.0, n))
    ang = rng.uniform(0.0, 2.0 * np.pi, n)
    off = rho[:, None] * (np.cos(ang)[:, None] * e1 + np.sin(ang)[:, None] * e2)
    q = scene.centres[j] + off + np.sqrt(np.maximum(r * r - rho * rho, 0.0))[:, None] * d
    return q, j


def simulate(c, vza, sza, raa, n_rays=40000, n_elem=40000, seed=0, n_target=1500, n_scenes=8):
    """
    One geometry (radians, gortt conventions), averaged over n_scenes independent
    tiles.  Standard errors come from the spread between tiles, which includes
    the scene-to-scene variance as well as the sampling error.
    """
    runs = [_simulate_one(c, vza, sza, raa, n_rays // n_scenes, n_elem // n_scenes,
                          seed * 1000 + k, n_target) for k in range(n_scenes)]
    out = {}
    for key in runs[0]:
        if key.endswith("_se"):
            continue
        vals = np.array([r[key] for r in runs], dtype=float)
        out[key] = vals.mean()
        out[key + "_se"] = vals.std(ddof=1) / np.sqrt(n_scenes) if n_scenes > 1 else 0.0
    return out


def _simulate_one(c, vza, sza, raa, n_rays, n_elem, seed, n_target):
    vp, ip = gk.prime(c, vza), gk.prime(c, sza)
    v = unit(vp, raa)       # towards the viewer; the sun is at azimuth 0
    s = unit(ip, 0.0)
    scene = Scene(c, seed=seed, n_target=n_target, extent=extent_for(c, vza, sza))
    rng = np.random.default_rng(seed + 1)
    out = {}

    # ---- rays ------------------------------------------------------------
    g = np.column_stack([rng.uniform(0.0, scene.L, (n_rays, 2)), np.zeros(n_rays)])
    hit, j, t = scene.first_hits(g, v)
    q = g + np.where(hit, t, 0.0)[:, None] * v
    n = np.zeros_like(q)
    n[hit] = (q[hit] - scene.centres[j[hit]]) / c.r
    facing = np.where(hit, n @ s > 0.0, True)
    excl = np.where(hit, j, 0)
    lit = facing.copy()
    idx_c = np.flatnonzero(hit & facing)
    if len(idx_c):
        lit[idx_c] = ~scene.blocked(q[idx_c], s, exclude=excl[idx_c])
    idx_g = np.flatnonzero(~hit)
    if len(idx_g):
        lit[idx_g] = ~scene.blocked(q[idx_g], s)

    for name, m in (("Kc", hit & lit), ("Kt", hit & ~lit), ("Kg", ~hit & lit), ("Kz", ~hit & ~lit)):
        p = m.mean()
        out[name] = p
        out[name + "_se"] = np.sqrt(p * (1.0 - p) / n_rays)

    # ---- elements --------------------------------------------------------
    qv, jv = _front_samples(scene, v, n_elem, rng)
    nv = (qv - scene.centres[jv]) / c.r
    fac = nv @ s > 0.0
    hid = scene.blocked(qv, v, exclude=jv)
    sh = np.zeros(n_elem, bool)
    fi = np.flatnonzero(fac)
    sh[fi] = scene.blocked(qv[fi], s, exclude=jv[fi])
    qs, js = _front_samples(scene, s, n_elem, rng)
    shs = scene.blocked(qs, s, exclude=js)

    for name, m in (("GcGv", fac), ("Mv", hid), ("PvMv", fac & hid),
                    ("PiMi", fac & sh), ("Po", fac & hid & sh), ("Mi", shs)):
        p = m.mean()
        out[name] = p
        out[name + "_se"] = np.sqrt(p * (1.0 - p) / n_elem)

    out["Kc_from_P"] = c.lam * np.pi * c.rr / np.cos(vp) * \
        (out["GcGv"] - out["PvMv"] - out["PiMi"] + out["Po"])
    out["n_crowns"] = len(scene.base)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vza", type=float, default=30.0)
    ap.add_argument("--sza", type=float, default=30.0)
    ap.add_argument("--raa", type=float, default=0.0)
    ap.add_argument("-r", type=float, default=gk.SOBS.r)
    ap.add_argument("-b", type=float, default=gk.SOBS.b)
    ap.add_argument("--lambda", dest="lam", type=float, default=gk.SOBS.lam)
    ap.add_argument("--h1", type=float, default=gk.SOBS.h1)
    ap.add_argument("--h2", type=float, default=gk.SOBS.h2)
    ap.add_argument("--rays", type=int, default=40000)
    ap.add_argument("--elements", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tile", type=int, default=1500, help="crowns per tile")
    ap.add_argument("--scenes", type=int, default=8, help="independent tiles to average")
    a = ap.parse_args()
    c = gk.Canopy(r=a.r, b=a.b, lam=a.lam, h1=a.h1, h2=a.h2)
    o = simulate(c, np.radians(a.vza), np.radians(a.sza), np.radians(a.raa),
                 a.rays, a.elements, a.seed, a.tile, a.scenes)
    K = gk.proportions(c, gk.GORTT, np.radians(a.vza), np.radians(a.sza), np.radians(a.raa))
    for name in ("Kc", "Kt", "Kg", "Kz"):
        print("%-6s MC %.4f +- %.4f   gortt %.4f" % (name, o[name], o[name + "_se"], K[name]))
    print("sum    %.6f" % sum(o[k] for k in ("Kc", "Kt", "Kg", "Kz")))
    for name in ("GcGv", "Mv", "Mi", "PvMv", "PiMi", "Po"):
        print("%-6s MC %.4f +- %.4f" % (name, o[name], o[name + "_se"]))
    print("Kc from P terms %.4f   (%d crowns per tile)" % (o["Kc_from_P"], o["n_crowns"]))


if __name__ == "__main__":
    main()
