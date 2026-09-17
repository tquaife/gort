#!/usr/bin/env python3
"""
Build a lidar-capable FLIGHT binary without touching the FLIGHT repository.

FLIGHT ships with its waveform lidar mode compiled out:

    flight.c:233   #define LMODE_PERMITTED 0

and check_data() exits with "'l' and 'p' lidar modes not supported" if you ask
for it.  The macro also gates MAX_NO_WAVE_BINS, which is 1 when lidar is off, so
it genuinely is a compile-time decision and cannot be worked around at runtime.

This script therefore copies flight.c out of the FLIGHT tree, flips that one
character in the copy, and builds it with -O3.  /home/tquaife/src/flight is only
ever read, and nothing from it lands in the gort repository -- the copy and the
binary go to the scratchpad (or wherever --builddir points), which is outside
both.

The same binary is used for the BRDF runs as well as the lidar runs, so there is
one build and one set of compiler flags behind every FLIGHT number quoted here.
Nothing else in flight.c is altered: run `diff` on the two sources and the only
hit is the macro.
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys

FLIGHT_SRC = "/home/tquaife/src/flight"
GORT_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def default_builddir():
    """A scratch directory outside both repositories."""
    scratch = os.environ.get("CLAUDE_SCRATCHPAD")
    if scratch and os.path.isdir(scratch):
        return os.path.join(scratch, "flight_build")
    return os.path.join(os.path.expanduser("~"), ".cache", "gort_flight_build")


def refuse_if_inside_a_repo(path):
    """The whole point is that neither tree is written to.  Enforce it."""
    real = os.path.realpath(path)
    for repo, name in ((FLIGHT_SRC, "the FLIGHT tree"), (GORT_SRC, "the gort repository")):
        rr = os.path.realpath(repo)
        if real == rr or real.startswith(rr + os.sep):
            sys.exit("build_flight: refusing to build inside %s: %s" % (name, real))


POINTER = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".flight_binary")


def record_binary(path):
    """
    Leave the binary's location where config.flight_binary() can find it.

    The build lands in a session scratchpad whose name is not predictable from
    the repository, so the path is written down rather than recomputed.  The
    file is gitignored: it is a local build artefact, not a result.
    """
    with open(POINTER, "w") as fp:
        fp.write(path + "\n")


def patch(raw):
    """
    Turn the lidar macro on.  Returns (new_bytes, n_replacements).

    Byte-for-byte: flight.c has CRLF line endings, and reading it as text would
    silently rewrite every line.  The copy must differ by exactly this macro so
    that `diff` is a meaningful check on it.
    """
    old = b"#define LMODE_PERMITTED 0"
    new = b"#define LMODE_PERMITTED 1"
    n = raw.count(old)
    return raw.replace(old, new), n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builddir", default=default_builddir(),
                    help="where to copy and build (default: %(default)s)")
    ap.add_argument("--flight", default=FLIGHT_SRC,
                    help="the FLIGHT source tree, read only (default: %(default)s)")
    ap.add_argument("--cc", default="gcc", help="compiler (default: %(default)s)")
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if the binary is already up to date")
    args = ap.parse_args()

    src = os.path.join(args.flight, "flight.c")
    if not os.path.exists(src):
        sys.exit("build_flight: %s not found" % src)

    refuse_if_inside_a_repo(args.builddir)
    os.makedirs(args.builddir, exist_ok=True)

    with open(src, "rb") as fp:
        original = fp.read()
    patched, n = patch(original)
    if n != 1:
        sys.exit("build_flight: expected exactly one LMODE_PERMITTED definition, found %d.\n"
                 "                flight.c has changed; check the macro by hand." % n)

    out_c = os.path.join(args.builddir, "flight.c")
    out_bin = os.path.join(args.builddir, "flight")
    digest = hashlib.sha256(patched).hexdigest()
    stamp = os.path.join(args.builddir, "flight.c.sha256")

    up_to_date = (os.path.exists(out_bin) and os.path.exists(stamp)
                  and open(stamp).read().strip() == digest)
    if up_to_date and not args.force:
        record_binary(out_bin)
        sys.stderr.write("build_flight: %s is up to date\n" % out_bin)
        print(out_bin)
        return

    with open(out_c, "wb") as fp:
        fp.write(patched)

    # FLIGHT's own `make release`.  The build is noisy -- dozens of pre-existing
    # -Wformat-overflow and -Wunused warnings -- so only show them if it fails.
    cmd = [args.cc, out_c, "-lm", "-O3", "-o", out_bin]
    sys.stderr.write("build_flight: %s\n" % " ".join(cmd))
    proc = subprocess.run(cmd, cwd=args.builddir, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        sys.exit("build_flight: compilation failed")

    with open(stamp, "w") as fp:
        fp.write(digest + "\n")

    # A patched copy is only trustworthy if it differs by the one line.
    diff = subprocess.run(["diff", src, out_c], capture_output=True, text=True)
    changed = [l for l in diff.stdout.splitlines() if l[:1] in "<>"]
    if len(changed) != 2:
        sys.stderr.write("build_flight: WARNING -- expected a one-line diff, got:\n")
        sys.stderr.write(diff.stdout)

    record_binary(out_bin)
    sys.stderr.write("build_flight: built %s\n" % out_bin)
    print(out_bin)


if __name__ == "__main__":
    main()
