"""
Shared helpers for the figure scripts: run gortt and parse what it prints.

Kept in one place so the scripts cannot drift apart on how they invoke the
model or read its output.
"""

import os
import subprocess
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Table I of Ni et al. (1999): input tree parameters for the SOBS site, the
# BOREAS Old Black Spruce stand that both Fig. 5 and Figs. 6-7 are drawn for.
#
# These are near enough gortt's own defaults -- gortt_set_defaults IS the SOBS
# site -- but are passed explicitly so the figures do not move if a default
# ever does.
#
# (The paper's Table I heading gives FAVD as 1/m^2, which is where gortt's
# usage message used to get it from.  It is an area per unit volume: 1/m.)
# ---------------------------------------------------------------------------
SOBS = [
    "-r", "0.76",       # horizontal crown radius (m)
    "-b", "2.7",        # vertical crown radius (m)
    "-lambda", "0.41",  # crown count density (1/m2)
    "-h1", "3.0",       # lower boundary of crown centres (m)
    "-h2", "8.5",       # upper boundary of crown centres (m)
    "-favd", "0.858",   # foliage area volume density (1/m)
]


def default_binary():
    """../gortt, relative to this file."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, os.pardir, "gortt")


def principal_plane(vza, sza):
    """
    (vza, vaa, sza, saa) rows for a signed principal-plane viewing zenith.

    The papers put the hotspot at *positive* viewing zenith, so positive is the
    backscatter direction (relative azimuth 0) and negative is forward scatter
    (relative azimuth 180).  gortt itself takes a non-negative viewing zenith
    plus azimuths, so a negative plotting angle becomes |vza| with the sun
    moved to the opposite azimuth.
    """
    rows = []
    for v in np.atleast_1d(vza):
        if v < 0.0:
            rows.append((abs(float(v)), 0.0, float(sza), 180.0))
        else:
            rows.append((float(v), 0.0, float(sza), 0.0))
    return rows


def run(binary, rows, wavelength, args=(), want="prop"):
    """
    Run gortt once over every angle in `rows` at a single wavelength.

    want="prop" returns the (n, 4) viewed proportions  Kc, Kg, Kt, Kz  -- the
    "[ ... ]" block that -prnprop prints.
    want="spec" returns the (n, 4) component signatures  C, G, T, Z  -- the
    "{ ... }" block that -prnspec prints.

    Either way the four values come back in the same C/G/T/Z order, which is
    the order gortt_rsurf composes the scene in.
    """
    flag, opener, closer = {
        "prop": ("-prnprop", "[", "]"),
        "spec": ("-prnspec", "{", "}"),
    }[want]

    header = "%d 1 %g\n" % (len(rows), wavelength)
    body = "".join("%g %g %g %g\n" % r for r in rows)

    # gortt prints %f -- six decimals -- unless GORTT_PREC asks for more.  At
    # the default the proportions only sum to one to about 1e-6.
    env = dict(os.environ, GORTT_PREC="17")

    argv = [binary, flag] + SOBS + list(args)

    try:
        proc = subprocess.run(
            argv, input=header + body, capture_output=True, text=True, env=env
        )
    except OSError as err:
        sys.exit("figures: cannot run %s: %s" % (binary, err))

    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit("figures: %s exited %d" % (binary, proc.returncode))

    out = []
    for line in proc.stdout.splitlines()[1:]:  # first line echoes the header
        if opener not in line:
            continue
        inside = line.split(opener, 1)[1].split(closer, 1)[0]
        vals = [float(x) for x in inside.split()]
        if len(vals) != 4:
            sys.exit("figures: expected 4 values, got %d: %s" % (len(vals), line))
        out.append(vals)

    if len(out) != len(rows):
        sys.exit("figures: asked for %d angles, parsed %d" % (len(rows), len(out)))

    return np.array(out)


def check_binary(path):
    if not os.path.exists(path):
        sys.exit("figures: %s not found -- run make first" % path)
