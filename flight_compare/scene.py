"""
Build matched scenes for gortt and FLIGHT, run both, and parse what they print.

FLIGHT's FPATH is "./" (flight.c:209), so every input and output path is
relative to the *current directory* rather than to the binary.  A run therefore
needs its own working directory holding in_flight.data, reverse.data,
crowns.data, SPEC/, LIDAR/, DATA/ and RESULTS/.  make_workdir() builds one;
DATA/ is symlinked to FLIGHT's own, which is the only thing this module reads
out of the FLIGHT tree, and it is never written to.

Crowns are written explicitly rather than left to FLIGHT's forest_gen3(), which
places a hardcoded 499 crowns and derives the scene size from the requested
cover.  An explicit crowns.data gets us GORT's canopy exactly: spheroids of
semi-axes (r, r, b), centres uniform in [h1, h2] and uniform in plan, with the
in-crown foliage area density in column 7.  That column is only honoured when
TOTAL_LAI in in_flight.data is negative (flight.c:7499), which is what we do,
so both models are given favd rather than a scene LAI.

FLIGHT's own RNG is seeded fixed unless argv[1] is "random", and "random" is
srand(time(NULL)), which collides between processes started in the same second.
Ensemble members are therefore different *scenes* drawn from our own seeded
numpy Generator, with FLIGHT left on its fixed seed.
"""

import csv
import glob
import hashlib
import math
import os
import shutil
import subprocess
import sys

import numpy as np

import config


# ---------------------------------------------------------------------------
# scene realisation
# ---------------------------------------------------------------------------

def crown_positions(scene, realisation):
    """
    (x, y, z) for every crown: uniform in plan over the scene, and with the
    centre height uniform in [h1, h2].  Seeded on (scene, realisation) so a
    realisation is reproducible and independent of the order things are run in.
    """
    # hashlib, not hash(): Python randomises str hashing per process, so hash()
    # would give a different scene in every worker.
    key = "%d|%s|%d" % (config.SEED, scene.name, realisation)
    seed = int(hashlib.sha256(key.encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    hw = scene.half_width
    x = rng.uniform(-hw, hw, scene.ncrowns)
    y = rng.uniform(-hw, hw, scene.ncrowns)
    z = rng.uniform(scene.h1, scene.h2, scene.ncrowns)
    return np.column_stack([x, y, z])


def tile_ghosts(scene, crowns):
    """
    Add the periodic images of every crown that straddles a wall.

    FLIGHT wraps *photons* at the scene walls (replace_if_out) but does not
    replicate *crowns*: the part of a crown lying outside the box is simply not
    there.  With 2 m crowns in a 46 m half-scene that loses a couple of percent
    of the foliage, which reads as a systematic deficit against GORT's infinite
    canopy rather than as noise.  Measured on dense-low-sphere before this was
    added, FLIGHT's fractional cover came out at 0.465 against the 0.478 that
    1 - exp(-lambda pi r^2) predicts; with the ghosts it lands on it.

    A ghost is a full crown at x +/- 2*half_width (and the same in y), kept
    only when its envelope actually reaches into the box.  Corners need both
    offsets, so the nine-cell neighbourhood is walked.
    """
    hw, r = scene.half_width, scene.r
    out = [crowns]
    for dx in (-2.0 * hw, 0.0, 2.0 * hw):
        for dy in (-2.0 * hw, 0.0, 2.0 * hw):
            if dx == 0.0 and dy == 0.0:
                continue
            g = crowns.copy()
            g[:, 0] += dx
            g[:, 1] += dy
            keep = (np.abs(g[:, 0]) - r < hw) & (np.abs(g[:, 1]) - r < hw)
            if keep.any():
                out.append(g[keep])
    all_crowns = np.vstack(out)
    if len(all_crowns) > config.MAX_CROWNS:
        sys.exit("scene: %s realisation needs %d crowns including ghosts, but "
                 "FLIGHT's MAX_NOSTANDS is %d -- lower config.REAL_CROWNS"
                 % (scene.name, len(all_crowns), config.MAX_CROWNS))
    return all_crowns


def write_crowns(path, scene, crowns, half_width=None):
    """
    crowns.data: "X_DIM Y_DIM NOSTANDS" then one row per crown,
    "type x y z Exy Ez density LAD leaf_size DBH" (read_canopy_data,
    flight.c:3562).  type 1 is an ellipsoid with semi-axes (Exy, Exy, Ez).
    The LAD and leaf_size columns are read and then discarded by FLIGHT --
    the global values win -- so they are written as zeros.

    half_width overrides the scene's own box, which run_crown_image.py needs:
    a single crown in the full 46 m box would be a handful of pixels across.
    """
    hw = scene.half_width if half_width is None else half_width
    with open(path, "w") as fp:
        fp.write("%.6f %.6f %d\n" % (hw, hw, len(crowns)))
        for x, y, z in crowns:
            fp.write("1 %.6f %.6f %.6f %.6f %.6f %.6f 0 %.4f 0.0\n"
                     % (x, y, z, scene.r, scene.b, scene.favd, config.LEAF_SIZE))


# ---------------------------------------------------------------------------
# FLIGHT input files
# ---------------------------------------------------------------------------

def write_in_flight(path, scene, mode, sza, vza, raa, nbands,
                    nphotons=10000, outname="out", leaf_size=None,
                    output_mode=0):
    """
    The 23 positional lines of in_flight.data (read_struct_data, flight.c:3777).
    Line order *is* the schema, so this is written out in full rather than
    edited in place.

    mode: "r" for the reverse-mode BRDF, "l" for the waveform lidar.
    TOTAL_LAI is negative so the per-crown density column of crowns.data is
    used.  AER_OPT is negative, which makes diffuse_est() return exactly zero
    (flight.c:2180) and the illumination purely direct -- gortt's -diffuse 0.

    output_mode is FLIGHT's OUTPUT_MODE: 0 writes no image, 1 an ENVI BIP
    binary and 2 that plus a floating-point text image (flight.c:7262).  It is
    0 for every run that only wants the scene-average spectrum; 2 is what
    run_crown_image.py uses to get the per-pixel crown brightness.
    """
    lad = " ".join("%.6f" % v for v in config.SPHERICAL_LAD)
    lines = [
        mode,                                           # 1  MODE
        "0",                                            # 2  ONED_FLAG (3D)
        "%.6f %.6f" % (sza, vza),                       # 3  SOLAR_ZENITH VIEW_ZENITH
        "%.6f %.6f" % (0.0, raa),                       # 4  SOLAR_AZIMUTH VIEW_AZIMUTH
        "%d" % nbands,                                  # 5  NO_WVBANDS
        "%d" % nphotons,                                # 6  NO_PHOTONS
        "-1 0",                                         # 7  TOTAL_LAI UNDERSTORY_LAI
        "1 0 0",                                        # 8  FRAC_GRN FRAC_SEN FRAC_BARK
        lad,                                            # 9  LAD
        "0",                                            # 10 SOILROUGH (Lambertian)
        "-1",                                           # 11 AER_OPT (< 0: direct only)
        "%.4f" % (config.LEAF_SIZE if leaf_size is None else leaf_size),   # 12 LF_SIZE
        "%.6f" % scene.nadir_cover,                     # 13 FRAC_COV (overwritten by FLIGHT)
        "0",                                            # 14 GROUND_SLOPE
        "f",                                            # 15 CROWN_SHAPE: read crowns.data
        "%.6f %.6f" % (scene.r, scene.b),               # 16 EXY EZ (unused; kept consistent)
        "%.6f %.6f" % (scene.z1, scene.z1 + config.CENTRE_RANGE),   # 17 MIN_HT MAX_HT (unused)
        "0",                                            # 18 DBH (no trunks: bark.spec absent)
        lad,                                            # 19 UNDERSTORY_LAD (unused)
        "0",                                            # 20 UNDERSTORY_HT
        "0.05",                                         # 21 UNDERSTORY_LF_SIZE
        outname,                                        # 22 OUTPUT_FILENAME
        "%d" % output_mode,                             # 23 OUTPUT_MODE
    ]
    with open(path, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def write_reverse(path, bounce, branch, im_size):
    """
    reverse.data, 10 positional lines (read_reverse_data, flight.c:4634).
    PRI_BAND_NO is -1: with it set, FLIGHT silently raises a bounce threshold
    of 1 to 2 (flight.c:7381), which would defeat the single-scattering runs.
    """
    lines = ["%d" % bounce, "%d" % branch, "%d" % im_size,
             "0", "-1", "2", "2000.0", "2000.0", "200.0", "10.0"]
    with open(path, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def write_spec(workdir, bands, rleaf=None, tleaf=None, rsoil=None):
    """
    SPEC/leaf.spec and SPEC/soil.spec.  One line per band, wavelength first and
    truncated to an integer by FLIGHT, and every spectral file must agree
    wavelength for wavelength or the run aborts.

    rleaf/tleaf/rsoil override the band values, which validate.py uses to build
    a black-leaf, white-soil scene.
    """
    spec = os.path.join(workdir, "SPEC")
    os.makedirs(spec, exist_ok=True)
    with open(os.path.join(spec, "leaf.spec"), "w") as fp:
        for b in bands:
            fp.write("%d\t%.6f\t%.6f\n" % (b.wl,
                                           b.rleaf if rleaf is None else rleaf,
                                           b.tleaf if tleaf is None else tleaf))
    with open(os.path.join(spec, "soil.spec"), "w") as fp:
        for b in bands:
            fp.write("%d\t%.6f\n" % (b.wl, b.rsoil if rsoil is None else rsoil))


def write_lidar(workdir, scene, bands):
    """
    LIDAR/source.data, flightline.data, source.spec and solar.spec.

    The beam entry point is sampled uniformly over a disc of radius
    range*tan(ifov/2) (flight.c:6890), so the footprint is the sensor IFOV and
    is uniformly illuminated -- the closest FLIGHT gets to GORT's statistically
    uniform canopy.  read_active_source() exits if that radius exceeds the
    scene half width, so it is checked here where the message can be useful.

    flightline.data's scattering-mode field is read into g_lidar_mode
    (flight.c:6465) and then never used anywhere, so single scattering is
    selected with BOUNCE_THRESHOLD in reverse.data instead.
    """
    lid = os.path.join(workdir, "LIDAR")
    os.makedirs(lid, exist_ok=True)

    fov_rad = config.LIDAR_FOV_RADIUS
    if fov_rad >= scene.half_width:
        sys.exit("scene: lidar footprint %.1f m does not fit in a %.1f m half-scene (%s)"
                 % (fov_rad, scene.half_width, scene.name))
    half_angle = math.atan(fov_rad / config.LIDAR_RANGE)

    with open(os.path.join(lid, "source.data"), "w") as fp:
        fp.write("0.0 0.0 %.4f\n" % config.LIDAR_RANGE)     # nadir, range to centre
        fp.write("0.0 0.0\n")                               # ground target
        fp.write("%.8f\n" % half_angle)                     # beam angular sd (bounds check only)
        fp.write("%.4f\n" % config.LIDAR_PULSE_NS)          # pulse sigma, ns
        fp.write("%.8f\n" % (2.0 * half_angle))             # FULL cone IFOV; code halves it
        fp.write("0.5\n")                                   # aperture radius, m
        fp.write("%.4f %.4f\n" % (config.LIDAR_REC_HIGH, config.LIDAR_REC_LOW))
        fp.write("%d\n" % config.LIDAR_NBINS)
        fp.write("0.0\n")                                   # dead time
        fp.write("y\n")                                     # afterpulse ON -- see below

    # One pulse at scene centre.  The waveform file is reopened "w" per
    # flightline sample (flight.c:6991), so more than one would overwrite.
    with open(os.path.join(lid, "flightline.data"), "w") as fp:
        fp.write("0 0\n0 0\n0.01 0\n15\nt\nn\n")

    with open(os.path.join(lid, "source.spec"), "w") as fp:
        for b in bands:
            # wavelength, energy (mJ), QE, bandwidth, atmos transmittance, optics
            fp.write("%d 1.0 1.0 1.0 1.0 1.0\n" % b.wl)
    with open(os.path.join(lid, "solar.spec"), "w") as fp:
        for b in bands:
            fp.write("%d 0.0 0.0\n" % b.wl)     # no solar background

    # Afterpulsing is switched ON, with an afterpulse of exactly zero energy.
    #
    # This is a workaround for a double free in FLIGHT, not a modelling choice.
    # read_active_source() closes the solar.spec handle at flight.c:6422 and
    # then closes the same FILE* again at flight.c:6444; with afterpulsing off
    # nothing has reopened it in between, so every lidar run aborts with
    # "free(): double free detected in tcache 2" before tracing a single
    # photon.  Turning afterpulsing on makes the second fclose() land on
    # afterpulse.spec instead, which is a file we control.
    #
    # A zero energy fraction does NOT switch the afterpulse off.  flight.c:6954
    # draws the selector as `r = rand() * RAND_MAX`, an int*int product that
    # overflows, so r is negative about half the time and `r < 0.0` at
    # flight.c:6955 is true about half the time whatever the fraction says.
    # Measured: with the pulse width set to 1e-4 ns and the afterpulse sd left
    # at 1 ns, 43% of the return picks up the afterpulse spread; with the
    # afterpulse sd also 1e-4 ns the whole return collapses into a single bin.
    #
    # So the afterpulse parameters must be written to be indistinguishable from
    # the real pulse rather than to be inert: sd equal to LIDAR_PULSE_NS and
    # delay zero.  Whichever branch fires, the photon then gets
    # gaussrand() * 0.3 * LIDAR_PULSE_NS metres of path and the waveform is the
    # one asked for -- checked against a run with the afterpulse sd set tiny,
    # which gives the same 0.15 m ground-return sd.  The file must exist:
    # with the flag on and the file missing, fpspec is NULL and the fclose()
    # above crashes again.
    with open(os.path.join(lid, "afterpulse.spec"), "w") as fp:
        for b in bands:
            # wavelength, energy frac, sd (ns), delay (ns).  sd tracks the pulse
            # width deliberately -- see above.
            fp.write("%d 0.0 %.4f 0.0\n" % (b.wl, config.LIDAR_PULSE_NS))


def make_workdir(path):
    """A clean FLIGHT working directory.  DATA/ is symlinked, never copied."""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(os.path.join(path, "RESULTS"))
    data = os.path.join(path, "DATA")
    if not os.path.exists(data):
        os.symlink(config.FLIGHT_DATA, data)
    return path


# ---------------------------------------------------------------------------
# running FLIGHT
# ---------------------------------------------------------------------------

class FlightError(RuntimeError):
    pass


def run_flight(workdir, binary=None, timeout=1800):
    """
    Run FLIGHT in workdir and return its stdout.

    FLIGHT exits 0 on its own validation errors (every check_data() failure is
    printf + exit(0)), so the return code is not a reliable signal.  The stdout
    is scanned for "ERROR" instead.
    """
    binary = binary or config.flight_binary()
    proc = subprocess.run([binary], cwd=workdir, capture_output=True,
                          text=True, timeout=timeout)
    out = proc.stdout
    if proc.returncode != 0 or "ERROR" in out or "not supported" in out:
        raise FlightError("FLIGHT failed in %s\n%s\n%s"
                          % (workdir, out[-3000:], proc.stderr[-2000:]))
    return out


def parse_flight_stdout(out):
    """Pull the scene diagnostics FLIGHT prints, for the record and the checks."""
    d = {}
    for line in out.splitlines():
        if line.startswith("Frac. cover (accounting for overlap):"):
            d["frac_cover"] = float(line.split(":")[1])
        elif line.startswith("Scene LAI (m2/m2) (accounting for overlap):"):
            d["scene_lai"] = float(line.split(":")[1])
        elif line.startswith("Scene LAI (m2/m2) (not accounting for overlap):"):
            d["scene_lai_nooverlap"] = float(line.split(":")[1])
        elif line.startswith("No. trees:"):
            d["ntrees"] = int(line.split(":")[1])
        elif line.startswith("Total LAI (m2/m2):"):
            parts = line.replace("Total LAI (m2/m2):", "").split("density (m2/m3):")
            d["total_lai"] = float(parts[0])
            if len(parts) > 1:
                d["density"] = float(parts[1])
    return d


def read_spectrum(path):
    """
    RESULTS/<name>-spectrum.csv.  Returns the "-1, mean" row: the scene-average
    BRF per band, which is the same number as the Rf column of the log file.
    The per-component rows above it are the sunlit/shaded soil and leaf
    breakdown, deliberately unused here.
    """
    with open(path) as fp:
        rows = list(csv.reader(fp))
    for row in rows:
        if row and row[0].strip() == "-1":
            return [float(v) for v in row[3:] if v.strip()]
    raise FlightError("no mean row in %s" % path)


def read_waveform(workdir):
    """
    RESULTS/flt-lidar_waveform_sens-*.out.  Columns are
    bin_index, one_way_range_m, then one value per band.

    The range is measured from path_zero, which sits at REC_ALT_ST above the
    ground target (flight.c:6845), so height above ground is
    REC_ALT_ST - range.  The filename embeds the derived LAI and cover, so it
    is globbed rather than predicted.
    """
    hits = glob.glob(os.path.join(workdir, "RESULTS", "flt-lidar_waveform_sens-*.out"))
    if len(hits) != 1:
        raise FlightError("expected one waveform file in %s, found %d"
                          % (workdir, len(hits)))
    arr = np.loadtxt(hits[0])
    z = config.LIDAR_REC_HIGH - arr[:, 1]
    return z, arr[:, 2:]


# ---------------------------------------------------------------------------
# running gortt
# ---------------------------------------------------------------------------

def gortt_env():
    """gortt prints %f -- four significant figures -- unless asked for more."""
    return dict(os.environ, GORTT_PREC="17")


def _gortt(argv, stdin=""):
    try:
        proc = subprocess.run(argv, input=stdin, capture_output=True,
                              text=True, env=gortt_env())
    except OSError as err:
        sys.exit("scene: cannot run %s: %s" % (argv[0], err))
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("scene: gortt exited %d\n  %s" % (proc.returncode, " ".join(argv)))
    return proc.stdout


def gortt_spectral_args(band, rleaf=None, tleaf=None, rsoil=None):
    """
    Switch PROSPECT and the Price soil off and hand gortt the same constants
    FLIGHT is given.  -refl_leaf and -trans_leaf must be supplied together
    (gortt.c:1418); both are spectrally flat, so one gortt run per band.
    """
    return ["-refl_leaf", "%.10g" % (band.rleaf if rleaf is None else rleaf),
            "-trans_leaf", "%.10g" % (band.tleaf if tleaf is None else tleaf),
            "-alb_soil", "%.10g" % (band.rsoil if rsoil is None else rsoil)]


def gortt_brdf(scene, band, directions, extra=(), binary=None):
    """
    One gortt process for a whole list of (vza, raa) directions at one solar
    zenith.  Returns a list of BRFs in the order given.

    gortt takes a non-negative view zenith and a pair of azimuths, and forms
    raa = saa - vaa (gortt.c:545), so the relative azimuth goes in as the solar
    azimuth with the view azimuth at zero.  raa = 0 is backscatter.
    """
    binary = binary or config.GORTT_BIN
    rows = [(vza, 0.0, sza, raa) for vza, sza, raa in directions]
    header = "%d 1 %g\n" % (len(rows), band.wl)
    body = "".join("%.10g %.10g %.10g %.10g\n" % r for r in rows)

    argv = ([binary] + scene.gortt_geometry() + gortt_spectral_args(band) + list(extra))
    out = _gortt(argv, header + body)

    brf = []
    for line in out.splitlines()[1:]:      # first line echoes the header
        vals = line.split()
        if len(vals) < 5:
            continue
        brf.append(float(vals[4]))
    if len(brf) != len(rows):
        sys.exit("scene: asked gortt for %d directions, parsed %d" % (len(rows), len(brf)))
    return brf


def gortt_probabilities(scene, binary=None):
    """
    gortt -W: the gap probability tables, written instead of a BRF.  Returns the
    file verbatim; see the format comment on gortt_write_prob_file for how to
    read it, and validate.gortt_nadir_gap for the row this is usually wanted for.
    """
    binary = binary or config.GORTT_BIN
    out = _gortt([binary, "-W"] + scene.gortt_geometry(), "")
    return out


def gortt_lidar(scene, band, waveform=True, binary=None, extra=()):
    """
    gortt -lidar.  Returns (header fields, rows).  With -lidar_waveform the
    rows are (z, w/binsize) on a uniform height grid; without it they are the
    per-layer profile (z, P, dP/dz, dlogP/dz, Rv).

    -lidar does not read stdin at all and exits before the per-angle loop.
    """
    binary = binary or config.GORTT_BIN
    argv = ([binary, "-lidar",
             "-lidar_wl", "%g" % band.wl,
             "-nlayers", "%d" % lidar_nlayers(scene)]
            + scene.gortt_geometry() + gortt_spectral_args(band) + list(extra))
    if waveform:
        argv += ["-lidar_waveform",
                 "-lidar_binsize", "%.10g" % lidar_height_binsize(),
                 "-lidar_pulse", "%.10g" % lidar_pulse_height_sigma()]
    out = _gortt(argv, "")
    lines = [l for l in out.splitlines() if l.strip()]
    head = [float(v) for v in lines[0].split()]
    rows = np.array([[float(v) for v in l.split()] for l in lines[1:]])
    return head, rows


def lidar_nlayers(scene):
    """
    gortt layers giving a thickness of about config.LIDAR_DZ.

    gortt lays nlayers layers across z1..z2, so dz = (z2 - z1)/(nlayers - 1);
    inverting that keeps the thickness, not the count, fixed across scenes.
    See the note by LIDAR_DZ for why that matters.
    """
    return int(round((scene.z2 - scene.z1) / config.LIDAR_DZ)) + 1


def lidar_height_binsize():
    """
    FLIGHT bins the round-trip path; one height metre is two path metres.
    g_path_limit = 2*(REC_HIGH - REC_LOW) over LIDAR_NBINS bins.
    """
    path_limit = 2.0 * (config.LIDAR_REC_HIGH - config.LIDAR_REC_LOW)
    return 0.5 * path_limit / config.LIDAR_NBINS


def lidar_pulse_height_sigma():
    """FLIGHT's pulse sigma is in ns; s = pulse_width*0.3 m of path (flight.c:6951)."""
    return 0.5 * config.LIDAR_PULSE_NS * 0.3
