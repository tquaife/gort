#!/usr/bin/env python3
"""
Image a single crown with FLIGHT and measure how its brightness varies across
the projected disc.

This is section 3 of CROWN-BRDF-NEXT.md, and it is meant to be done before any
model code is written.  GORT treats the sunlit crown as a disc of *uniform*
radiance: `Kc` is the area of the sunlit-and-viewed part of the envelope and
`C` is one number for the whole of it.  The candidate replacements disagree
about exactly that -- a Lambert sphere is limb-darkened as sqrt(1 - u^2) at
zero phase, Lommel-Seeliger is perfectly flat -- so the profile across the disc
discriminates between them immediately, and far more sharply than the
whole-scene BRDF, which convolves the crown law with the Li-Strahler geometry
and the background.

The scene is one crown at the origin over a black soil, so the only thing in
the image is the crown and there are no neighbours to scatter into it.

    ./run_crown_image.py --dry-run     # how many runs, and how long
    ./run_crown_image.py --quick       # one scene, coarse -- smoke test
    ./run_crown_image.py               # the full sweep

Two CSVs come out:

    results/crown_profile.csv   ring-averaged radiance factor against u
    results/crown_phase.csv     the disc mean against view angle, i.e. the
                                crown's own phase function

Geometry.  FLIGHT's reverse mode launches one ray per pixel from the plane
z = canopy_top + 0.01 at (x, y) on a regular grid over the whole box, in the
view direction (flight.c:7649-7690); off_x and off_y are zeroed two lines
above, so the grid is not shifted with view zenith.  A pixel is therefore a
known straight line through the scene, and the crown is a known spheroid, so
the impact parameter is analytic -- nothing is fitted or registered.  Scaling
z by r/b turns the spheroid into a sphere of radius r, which is GORT's own
prime transform, so the u that comes out is the u of the chord average in
gortt_crown_reflectance.

That mapping is checked rather than assumed: every run reports how many lit
pixels fall outside the predicted silhouette, and a single one is fatal.  It is
the test that would catch a wrong canopy top, a transposed image, or a photon
wrapped round the scene wall and back onto the crown -- and it did catch the
last of those, which is why the box is sized per view zenith rather than fixed
(see half_width_for).
"""

import argparse
import csv
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc

# How far outside u = 1 a lit pixel may sit before it counts as a fault rather
# than a grazing ray; see the check in image_one.
LIMB_TOL = 0.01


def crown_centre_height(scn):
    """
    The one crown sits at the middle of the height range the full scenes draw
    from.  Nothing here depends on it -- the soil is black and the crown is
    alone -- but it keeps the canopy the same shape as its namesake.
    """
    return 0.5 * (scn.h1 + scn.h2)


def canopy_top_in(scn):
    """
    Where FLIGHT launches the reverse rays from.

    flight.c:3710-3713: g_canopy_top is the highest crown top, floored at
    understory_ht + 0.2 (understory_ht is 0 here), and g_canopy_top_in is that
    plus 0.01.  FLIGHT never prints it, hence the reimplementation -- and hence
    the silhouette check, which is what would notice if this drifted.
    """
    return max(crown_centre_height(scn) + scn.b, 0.2) + 0.01


def view_vector(vza, raa):
    """
    The direction a reverse ray travels, pointing down into the scene.

    flight.c:7677-7686 sets theta = pi - vza and phi = raa - pi, then
    (x, y, z) = (sin(theta) sin(phi), sin(theta) cos(phi), cos(theta)) -- note
    x takes the sine and y the cosine, which is FLIGHT's convention, not the
    usual one.
    """
    v, a = math.radians(vza), math.radians(raa)
    return np.array([-math.sin(v) * math.sin(a),
                     -math.sin(v) * math.cos(a),
                     -math.cos(v)])


def phase_angle(sza, vza, raa):
    """
    Angle between the view direction and the sun direction, both pointing up
    and away from the scene.  raa = 0 is backscatter, so this is |vza - sza|
    in the principal plane on the hotspot side.
    """
    s, v, a = math.radians(sza), math.radians(vza), math.radians(raa)
    c = (math.cos(v) * math.cos(s) + math.sin(v) * math.sin(s) * math.cos(a))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def impact_parameter(scn, im_size, half_width, vza, raa):
    """
    u = (perpendicular distance of the pixel's ray from the crown centre) / r,
    evaluated after the ellipsoid -> sphere transform, for every pixel.

    u < 1 is exactly the projected silhouette of the crown, whatever the view
    zenith, which is what makes it the right coordinate: it is the radius
    across the projected disc, normalised so that the limb is always at 1.  The
    projected-area weight is 2u du.
    """
    step = 2.0 * half_width / im_size
    axis = -half_width + 0.5 * step + step * np.arange(im_size)
    # flight.c:7649 loops yst outside xst, so rows are y and columns are x.
    x, y = np.meshgrid(axis, axis)

    zc, squeeze = crown_centre_height(scn), scn.r / scn.b

    origin = np.stack([x, y, np.full_like(x, (canopy_top_in(scn) - zc) * squeeze)], axis=-1)
    d = view_vector(vza, raa)
    d = np.array([d[0], d[1], d[2] * squeeze])
    d = d / np.linalg.norm(d)

    return np.linalg.norm(np.cross(origin, d), axis=-1) / scn.r


def half_width_for(scn, vza, floor, margin=1.0):
    """
    How wide the box has to be for the crown to be genuinely isolated at this
    view zenith.

    FLIGHT is periodic: replace_if_out (flight.c:943) puts a photon that leaves
    a wall back on the opposite wall with the same direction and the same
    height.  A wrapped ray re-enters at least (hw - r) from the crown and
    descends (hw - r)/tan(vza) getting back to it, while the most height it has
    to spend is z_top - (zc - b) = 2b.  So it can reach the crown again only if

        hw < r + 2b tan(vza)

    and the box has to clear that.  At hw = 12 m this bites for the prolate
    crowns beyond about 45 degrees: b is 5 m, and the measured failure was 299
    lit pixels out to u = 2.5 at vza 60, a clean second image of the crown.

    Sizing it per direction rather than at the worst case matters for the cost,
    which goes as hw^2: a fixed 30 m box would put every run at the price of
    the most expensive one.
    """
    reach = scn.r + (2.0 * scn.b + 0.01) * math.tan(math.radians(vza))
    return max(floor, reach + margin)


def read_float_image(path, nbands):
    """
    FLIGHT's floating-point image: one text line per image row, and within a
    row "%7.5f " per waveband per pixel (flight.c:7766, 7787).

    In reverse mode every pixel is written -- the guard at flight.c:7743 is
    `(g_mode == 'r') || ...`, so the "dud photon" branch that would make the
    raster ragged cannot be taken -- which is why this can be reshaped rather
    than having to be walked with the mask file.
    """
    rows = [line.split() for line in open(path) if line.strip()]
    if not rows:
        sys.exit("run_crown_image: %s is empty" % path)
    widths = set(len(r) for r in rows)
    if len(widths) != 1:
        sys.exit("run_crown_image: ragged image in %s (row lengths %s)"
                 % (path, sorted(widths)))
    a = np.array([[float(v) for v in r] for r in rows], dtype=float)
    if a.shape[1] != a.shape[0] * nbands:
        sys.exit("run_crown_image: %s is %d x %d, not square in %d bands"
                 % (path, a.shape[0], a.shape[1], nbands))
    return a.reshape(a.shape[0], a.shape[0], nbands)


def image_one(scn, wd, sza, vza, raa, im_size, half_width, bands, check):
    """One FLIGHT run: returns the image and the per-pixel impact parameter."""
    # The box, and so the image, is resized per direction, which means
    # crowns.data and reverse.data are rewritten here rather than once per
    # scene.  One crown at the origin over a black soil: the image holds
    # nothing else.
    sc.write_crowns(os.path.join(wd, "crowns.data"), scn,
                    np.array([[0.0, 0.0, crown_centre_height(scn)]]),
                    half_width=half_width)
    sc.write_reverse(os.path.join(wd, "reverse.data"),
                     config.BRDF_BOUNCE, config.BRDF_BRANCH, im_size)
    sc.write_in_flight(os.path.join(wd, "in_flight.data"), scn, "r",
                       sza, vza, raa, len(bands), output_mode=2)
    sc.run_flight(wd)

    img = read_float_image(os.path.join(wd, "RESULTS", "out.txt"), len(bands))
    if img.shape[0] != im_size:
        sys.exit("run_crown_image: asked for %d pixels, got %d" % (im_size, img.shape[0]))
    u = impact_parameter(scn, im_size, half_width, vza, raa)

    if check:
        # The tolerance is for rays that graze the limb, where whether the
        # ellipsoid is hit at all comes down to the last few digits: FLIGHT
        # intersects in float and carries 0.01 m offsets on the canopy top, so
        # a pixel whose centre this computes at u = 1.0003 can still be lit.
        # One such pixel in 218,058 was what prompted the tolerance.  It is far
        # below anything a real fault produces -- the wrapped second image of
        # the crown that this check was written for put 299 pixels out at
        # u = 2.5 -- so the ghost detection survives it intact.
        lit = img.max(axis=-1) > 0.0
        graze = int((lit & (u > 1.0) & (u <= 1.0 + LIMB_TOL)).sum())
        outside = int((lit & (u > 1.0 + LIMB_TOL)).sum())
        print("    check: %7d lit, %7d in silhouette, %d grazing, %d outside, "
              "max u lit %.4f" % (lit.sum(), (u <= 1.0).sum(), graze, outside,
                                  u[lit].max() if lit.any() else float("nan")))
        if outside:
            sys.exit("run_crown_image: %d lit pixels lie more than %g outside the "
                     "predicted silhouette (max u %.3f) -- the geometry mapping is "
                     "wrong" % (outside, LIMB_TOL, u[lit].max()))
    return img, u


def profile_rows(scn, band_names, img, u, sza, vza, raa, nu):
    """
    Ring-average the radiance factor in nu equal-width bins of u.

    Equal width in u, not in u^2: the interesting behaviour is at the limb and
    the outer rings hold the most pixels anyway, so equal width already puts
    the sampling where it is needed.
    """
    inside = u <= 1.0
    edges = np.linspace(0.0, 1.0, nu + 1)
    idx = np.clip(np.digitize(u[inside], edges) - 1, 0, nu - 1)

    rows = []
    for k, name in enumerate(band_names):
        v = img[:, :, k][inside]
        for j in range(nu):
            sel = idx == j
            n = int(sel.sum())
            if not n:
                continue
            vals = v[sel]
            rows.append((scn.name, name, sza, vza, raa, phase_angle(sza, vza, raa),
                         0.5 * (edges[j] + edges[j + 1]), n,
                         vals.mean(), vals.std(ddof=1) / math.sqrt(n) if n > 1 else 0.0,
                         float((vals > 0.0).mean())))
    return rows


def phase_rows(scn, band_names, img, u, sza, vza, raa):
    """
    The disc mean: the crown's own radiance factor, with no soil and no
    neighbours.  This is what GORT's Kc*C + Kt*T has to add up to over the
    projected crown area, so it is directly comparable rather than being a
    shape only.
    """
    inside = u <= 1.0
    rows = []
    for k, name in enumerate(band_names):
        v = img[:, :, k][inside]
        rows.append((scn.name, name, sza, vza, raa, phase_angle(sza, vza, raa),
                     int(inside.sum()), v.mean(),
                     v.std(ddof=1) / math.sqrt(inside.sum()),
                     float((v > 0.0).mean())))
    return rows


def directions(vzas, raas):
    """The (vza, raa) pairs; nadir is one direction whatever the azimuth."""
    out = []
    for raa in raas:
        for vza in vzas:
            if vza == 0.0 and raa != raas[0]:
                continue
            out.append((vza, raa))
    return out


def sweep(scenes, szas, dirs, pixel, floor, nu, workroot, check):
    bands = config.BRDF_BANDS
    band_names = [b.name for b in bands]
    prof, phase = [], []

    for scn in scenes:
        wd = os.path.join(workroot, scn.name)
        sc.make_workdir(wd)
        sc.write_spec(wd, bands, rsoil=0.0)

        for sza in szas:
            for vza, raa in dirs:
                hw = half_width_for(scn, vza, floor)
                im_size = int(round(2.0 * hw / pixel))
                t0 = time.time()
                img, u = image_one(scn, wd, sza, vza, raa, im_size, hw, bands, check)
                prof += profile_rows(scn, band_names, img, u, sza, vza, raa, nu)
                phase += phase_rows(scn, band_names, img, u, sza, vza, raa)
                sys.stderr.write("run_crown_image: %-20s sza %2g vza %2g raa %3g  "
                                 "hw %5.1f m  %4d px  (%.0f s)\n"
                                 % (scn.name, sza, vza, raa, hw, im_size,
                                    time.time() - t0))
    return prof, phase


def write_csv(path, header, rows, fmt):
    with open(path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(header)
        for r in rows:
            w.writerow([f % v if f else v for f, v in zip(fmt, r)])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenes", nargs="+", default=config.CROWN_IMAGE_SCENES,
                    help="scene names (default: %(default)s)")
    ap.add_argument("--szas", nargs="+", type=float, default=config.CROWN_IMAGE_SZAS)
    ap.add_argument("--pixel", type=float, default=config.CROWN_IMAGE_PIXEL,
                    help="pixel size, m (default: %(default)s)")
    ap.add_argument("--half-width", type=float, default=config.CROWN_IMAGE_HALFWIDTH,
                    help="minimum scene half width, m; raised per direction "
                         "by half_width_for (default: %(default)s)")
    ap.add_argument("--nu", type=int, default=config.CROWN_IMAGE_NU)
    ap.add_argument("--workdir", default=os.path.join(config.WORK_DIR, "crown_image"))
    ap.add_argument("--outdir", default=config.RESULTS_DIR)
    ap.add_argument("--quick", action="store_true",
                    help="one scene, 0.1 m pixels, three view zeniths")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-check", action="store_true",
                    help="skip the silhouette check (it is cheap; do not)")
    args = ap.parse_args()

    vzas, raas = config.CROWN_IMAGE_VZAS, config.CROWN_IMAGE_RAAS
    if args.quick:
        args.scenes = args.scenes[:1]
        args.pixel = 0.1
        vzas, raas = [0.0, 30.0, 60.0], [0.0, 180.0]

    scenes = [config.scene_by_name(n) for n in args.scenes]
    dirs = directions(vzas, raas)
    nruns = len(scenes) * len(args.szas) * len(dirs)

    # The box, and so the image, is sized per direction, so the cost has to be
    # summed rather than multiplied out.  16 s per 1200^2 image was measured on
    # a black single-crown scene; every ray that misses the crown ends at the
    # soil on its first intersection, so the rate is flat across the scenes.
    px = [int(round(2.0 * half_width_for(s, v, args.half_width) / args.pixel))
          for s in scenes for _ in args.szas for v, _ in dirs]
    secs = sum(16.0 * (n / 1200.0) ** 2 for n in px)

    print("run_crown_image: %d scenes x %d sza x %d directions = %d FLIGHT runs"
          % (len(scenes), len(args.szas), len(dirs), nruns))
    print("                 %.3f m pixels, %.0f across a crown radius; images "
          "%d to %d px" % (args.pixel, config.CROWN_RADIUS / args.pixel,
                           min(px), max(px)))
    print("                 about %.0f min" % (secs / 60.0))
    if args.dry_run:
        return

    if min(px) < 40:
        sys.exit("run_crown_image: %d pixels across the whole box is too coarse "
                 "to resolve a crown" % min(px))

    os.makedirs(args.outdir, exist_ok=True)
    prof, phase = sweep(scenes, args.szas, dirs, args.pixel, args.half_width,
                        args.nu, args.workdir, not args.no_check)

    write_csv(os.path.join(args.outdir, "crown_profile.csv"),
              ["scene", "band", "sza", "vza", "raa", "phase", "u", "n",
               "brf", "sem", "frac_lit"], prof,
              [None, None, "%g", "%g", "%g", "%.4f", "%.4f", None,
               "%.8g", "%.4g", "%.5f"])
    write_csv(os.path.join(args.outdir, "crown_phase.csv"),
              ["scene", "band", "sza", "vza", "raa", "phase", "n",
               "brf", "sem", "frac_lit"], phase,
              [None, None, "%g", "%g", "%g", "%.4f", None,
               "%.8g", "%.4g", "%.5f"])
    print("run_crown_image: wrote crown_profile.csv (%d rows) and "
          "crown_phase.csv (%d rows)" % (len(prof), len(phase)))


if __name__ == "__main__":
    main()
