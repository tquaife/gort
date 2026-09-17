"""
Single source of truth for the GORT/FLIGHT intercomparison: the scene matrix,
the optical properties, the angle grids and the paths.

Everything the two models are told about a scene is derived here, so the one
thing the comparison must not get wrong -- that both models are being shown the
same canopy -- is decided in one place and checked by validate.py.

Conventions
-----------
Crowns are spheroids of horizontal radius r and vertical radius b, centres
uniformly distributed in [h1, h2] and Poisson distributed in plan at lambda
crowns per square metre.  Foliage fills the crown at a constant area volume
density favd (1/m).  That is GORT's model exactly (Ni et al. 1999); FLIGHT is
given the same crowns as an explicit list.

Relative azimuth raa = 0 is the backscatter direction (the hotspot) in both
models: gortt computes raa = saa - vaa, and FLIGHT's reverse-mode view vector
negates to (VIEW_ZENITH, VIEW_AZIMUTH) in the same frame as its solar vector, so
the retro-reflection direction is VIEW_AZIMUTH == SOLAR_AZIMUTH.  We therefore
run FLIGHT with SOLAR_AZIMUTH 0 and VIEW_AZIMUTH raa, and gortt with vaa 0 and
saa raa.  validate.py checks this against the models rather than trusting it.

Why favd and not LAI
--------------------
gortt's -LAI converts with favd = LAI / (lambda (4/3) pi r^2 b) (gortt.c:1435),
which counts crown volume without correcting for overlap.  Specifying favd
directly on both sides makes the foliage density identical by construction and
leaves the overlap treatment as something this comparison measures instead of
something it assumes.  Scene LAI is reported as a diagnostic from both models.
"""

import math
import os

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
GORT_ROOT = os.path.dirname(HERE)
GORTT_BIN = os.path.join(GORT_ROOT, "gortt")
FLIGHT_SRC = "/home/tquaife/src/flight"      # read only, never written to
FLIGHT_DATA = os.path.join(FLIGHT_SRC, "DATA")

RESULTS_DIR = os.path.join(HERE, "results")
FIGS_DIR = os.path.join(HERE, "figs")
WORK_DIR = os.path.join(HERE, "work")        # gitignored, per-run scratch


def flight_binary():
    """
    The patched, lidar-enabled build; see build_flight.py.

    build_flight.py writes the path it used into .flight_binary, because the
    build lands in a scratch directory that cannot be recomputed from here.
    """
    pointer = os.path.join(HERE, ".flight_binary")
    if os.path.exists(pointer):
        path = open(pointer).read().strip()
        if path and os.path.exists(path):
            return path
    import build_flight
    return os.path.join(build_flight.default_builddir(), "flight")


# ---------------------------------------------------------------------------
# spectral properties
#
# Three bands, one grid shared by every file: FLIGHT requires SPEC/*.spec and
# LIDAR/*.spec to agree wavelength for wavelength, and truncates the wavelength
# to an integer, so these must be whole numbers.  Values are plausible broadleaf
# and mid-brightness soil, chosen rather than modelled -- PROSPECT and the Price
# EOFs are both switched off so the two models get identical numbers.
# ---------------------------------------------------------------------------

class Band:
    def __init__(self, wl, name, rleaf, tleaf, rsoil, role):
        self.wl = wl            # nm, integer
        self.name = name
        self.rleaf = rleaf
        self.tleaf = tleaf
        self.rsoil = rsoil
        self.role = role        # "brdf" or "lidar"

    def __repr__(self):
        return "Band(%d, %r)" % (self.wl, self.name)


BANDS = [
    Band(660, "red", 0.08, 0.04, 0.12, "brdf"),
    Band(860, "nir", 0.45, 0.45, 0.22, "brdf"),
    Band(1064, "lidar", 0.50, 0.44, 0.25, "lidar"),
]

BRDF_BANDS = [b for b in BANDS if b.role == "brdf"]
LIDAR_BAND = [b for b in BANDS if b.role == "lidar"][0]

# Spherical leaf angle distribution: the fraction of leaf normals in each 10 deg
# zenith bin is cos(lower) - cos(upper).  gortt is hardwired to G = 0.5 (LAD_05,
# gortt.c:870) and this is the distribution that gives it.  Same numbers as the
# shipped in_flight.data.
SPHERICAL_LAD = [math.cos(math.radians(10 * i)) - math.cos(math.radians(10 * (i + 1)))
                 for i in range(9)]

LEAF_SIZE = 0.05        # m, FLIGHT's within-crown hotspot scale.  gortt has no
                        # equivalent -- its hotspot is purely geometric-optical.

# ---------------------------------------------------------------------------
# scene matrix: 2 covers x 2 leaf areas x 3 crown shapes
# ---------------------------------------------------------------------------

CROWN_RADIUS = 2.0          # r, m.  Fixed, so cover and shape vary independently.
CENTRE_RANGE = 4.0          # h2 - h1, m
BASE_CLEARANCE = 2.0        # z1 = h1 - b, m: foliage always starts 2 m up
MAX_CROWNS = 500            # FLIGHT's MAX_NOSTANDS (flight.c); hard limit
REAL_CROWNS = 430           # crowns actually drawn.  The rest of the budget goes
                            # on the periodic images of the crowns that straddle a
                            # wall -- see scene.tile_ghosts() for why they matter.

SHAPES = [                  # (name, b/r)
    ("oblate", 0.5),
    ("sphere", 1.0),
    ("prolate", 2.5),
]
COVERS = [                  # (name, PCC = lambda pi r^2, the nadir cover before overlap)
    ("sparse", 0.20),
    ("dense", 0.65),
]
LAIS = [                    # (name, target scene LAI, m2/m2, ignoring crown overlap)
    ("low", 1.0),
    ("high", 4.0),
]


class Scene:
    """One canopy, in the terms both models need."""

    def __init__(self, cover_name, pcc, lai_name, lai, shape_name, br):
        self.cover_name = cover_name
        self.lai_name = lai_name
        self.shape_name = shape_name
        self.name = "%s-%s-%s" % (cover_name, lai_name, shape_name)

        self.r = CROWN_RADIUS
        self.b = br * self.r
        self.br = br
        self.pcc = pcc
        self.lam = pcc / (math.pi * self.r ** 2)        # crowns per m2

        self.h1 = self.b + BASE_CLEARANCE
        self.h2 = self.h1 + CENTRE_RANGE
        self.z1 = self.h1 - self.b                      # = BASE_CLEARANCE
        self.z2 = self.h2 + self.b

        self.crown_volume = (4.0 / 3.0) * math.pi * self.r ** 2 * self.b
        self.target_lai = lai
        # favd = LAI / (lambda * crown volume); identical on both sides.
        self.favd = lai / (self.lam * self.crown_volume)

        # The crown count fixes the scene size at the given density.
        self.ncrowns = REAL_CROWNS
        area = self.ncrowns / self.lam
        self.half_width = 0.5 * math.sqrt(area)         # X_DIM = Y_DIM

    @property
    def nadir_cover(self):
        """Projected crown cover at nadir, corrected for overlap."""
        return 1.0 - math.exp(-self.pcc)

    def gortt_geometry(self):
        """The crown geometry arguments, shared by every gortt call."""
        return ["-r", "%.10g" % self.r,
                "-b", "%.10g" % self.b,
                "-h1", "%.10g" % self.h1,
                "-h2", "%.10g" % self.h2,
                "-lambda", "%.10g" % self.lam,
                "-favd", "%.10g" % self.favd,
                "-diffuse", "0"]        # direct beam only; matches AER_OPT < 0

    def __repr__(self):
        return "Scene(%s)" % self.name


def scenes():
    out = []
    for cover_name, pcc in COVERS:
        for lai_name, lai in LAIS:
            for shape_name, br in SHAPES:
                out.append(Scene(cover_name, pcc, lai_name, lai, shape_name, br))
    return out


def scene_by_name(name):
    for s in scenes():
        if s.name == name:
            return s
    raise KeyError(name)


# ---------------------------------------------------------------------------
# angle grids
#
# FLIGHT runs one view direction per process (MULT_VIEWS is compiled out), so
# the grid is the cost.  The polar grid is what the contour plots are drawn
# from; the principal plane is sampled more finely on top of it because the
# hotspot is narrow and a 10 degree grid smears it.
# ---------------------------------------------------------------------------

SZAS = [15.0, 30.0, 45.0]

POLAR_VZA = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
POLAR_RAA = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]
PP_VZA = [5.0 * i for i in range(15)]           # 0..70 in 5 deg steps
PP_RAA = [0.0, 180.0]

# gortt is cheap, so the contour panels use a much finer grid of its own.
GORTT_DENSE_VZA = [2.5 * i for i in range(29)]  # 0..70
GORTT_DENSE_RAA = [5.0 * i for i in range(37)]  # 0..180


def flight_directions():
    """
    The (vza, raa) pairs FLIGHT is run at: the polar grid plus a finer principal
    plane.  vza = 0 is a single direction whatever the azimuth, so it appears
    once.
    """
    seen, out = set(), []
    for vza in POLAR_VZA:
        for raa in POLAR_RAA:
            key = (round(vza, 4), 0.0 if vza == 0.0 else round(raa, 4))
            if key not in seen:
                seen.add(key)
                out.append(key)
    for vza in PP_VZA:
        for raa in PP_RAA:
            key = (round(vza, 4), 0.0 if vza == 0.0 else round(raa, 4))
            if key not in seen:
                seen.add(key)
                out.append(key)
    return sorted(out)


# ---------------------------------------------------------------------------
# Monte Carlo settings
# ---------------------------------------------------------------------------

BRDF_IM_SIZE = 150          # reverse mode traces im_size^2 paths
BRDF_BRANCH = 2             # diffuse sub-samples per pixel at the first bounce
BRDF_BOUNCE = 8             # scattering orders
BRDF_REALISATIONS = 3       # independent scene draws, averaged

LIDAR_PHOTONS = 200000
LIDAR_REALISATIONS = 10     # a 25 m footprint holds ~30 crowns when sparse
LIDAR_BOUNCE_SINGLE = 1     # the comparison gortt's first-order waveform matches
LIDAR_BOUNCE_MULTI = 8      # overlay showing what multiple scattering adds

# lidar instrument.  Airborne-ish: close enough that the footprint is a sensible
# fraction of the scene, far enough that the beam is effectively collimated.
LIDAR_RANGE = 1000.0        # m, sensor to scene centre
LIDAR_FOV_RADIUS = 25.0     # m on the ground; must be < scene half width
LIDAR_PULSE_NS = 1.0        # ns, Gaussian sigma (= 0.3 m, or 0.15 m in height)
LIDAR_REC_HIGH = 40.0       # m above the ground target; recording window top
LIDAR_REC_LOW = -10.0       # m; window bottom, below ground to catch the tail
LIDAR_NBINS = 500           # bins across the window (0.1 m of height each)

# gortt's reporting grid for -lidar.  Set from a target layer thickness rather
# than as a fixed layer count, because gortt spaces its nlayers layers over
# z1..z2 (gortt.c:987) and that span is 6 m for the oblate crowns but 14 m for
# the prolate ones -- a fixed count gives them 2.3x coarser layers.
#
# This comparison found that gortt's waveform binning aliased: it dropped the
# whole of each layer-to-layer increment of P into the nearest range bin, so
# once dz exceeded the bin size the bins alternated between receiving one layer
# and receiving none and the waveform came out as a staircase.  Because gortt
# spreads its layers over z1..z2, that showed up in the 14 m prolate canopies
# and not the 6 m oblate ones, and looked like a property of prolate crowns.
# Fixed upstream in gort commit 0ce47c1, which splits each contribution between
# the two bins it falls between; the prolate raggedness at 120 layers went from
# 0.037 to 0.007 (rms second difference of the normalised canopy return over its
# peak).
#
# LIDAR_DZ remains, for a different reason: it gives every scene the same layer
# thickness, so the residual numerical roughness is the same in all twelve and
# the per-scene RMSDs against FLIGHT can be compared with each other.  A fixed
# layer count cannot, since z2 - z1 is 6 m for the oblate crowns and 14 m for
# the prolate ones.  0.02 m is five layers per 0.1 m bin; the cost is linear in
# the layer count and the whole gortt lidar pass is about three minutes.
LIDAR_DZ = 0.02             # m, target gortt layer thickness

SEED = 20260916             # our own RNG; FLIGHT's is fixed and cannot be set
