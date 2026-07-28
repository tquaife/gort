# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make            # builds ./gortt (C + Fortran, linked with gfortran)
make clean      # removes *.o *.mod etc (does NOT remove the ./gortt binary)
make gortt.o    # rebuild a single translation unit
```

`gcc` and `gfortran` are both required: the model core is C, PROSPECT-D is Fortran 90 and is
linked in directly (`makefile` uses `vpath` to find `PROSPECT-D/*.f90`). Headers live in
`include/` and are pulled in with `-I./include`, so sources use `#include<gortt.h>`, not
a relative path.

Expected build noise (pre-existing, not introduced by your change):
* many `-Wc-binding-type` warnings from the `bind(C)` interface in `PROSPECT-D/prospect_DB.f90`
* `gortt.c:118` implicit declaration of `gortt_gap_probabilities_Q08` — the function is defined
  in `gortt_pn_kopen.c` but missing from `include/gortt.h`
* `gortt.c:116` dangling-else warning

## Testing

```bash
make test       # 36 golden-file cases, ~30 s; or: cd tests && ./run.sh
make bless      # re-capture goldens -- only on a commit MEANT to change the numbers
```

`tests/` is a golden-file regression suite: `gortt` is a deterministic stdin→stdout filter, so a
byte comparison of its output is a complete check of the numerical path. Cases are registered in
`tests/cmds`; stdout, stderr and exit status are all captured. A failure prints the largest
absolute and relative deviation via `tests/numdiff.awk`, so a last-digit rounding artefact is
distinguishable from a broken model. See `tests/README.md`.

**`run.sh` sets `GORTT_PREC=17`.** Without it `gortt` prints `%f` — 4 significant figures — which
is too coarse to catch a small regression: a 1e-7 perturbation of `favd` was visible in 1 of 30
BRF cases at `%f`, and 25 of 36 at full precision. The `-W` cases already print `%0.40f` and are
the most sensitive probe in the suite.

Two cases deliberately pin *known-bad* behaviour (`useprobs-extreme`, `badanth`) so that fixing
those bugs shows up as a reviewable change — see `tests/README.md`.

## Running by hand

Verify changes by comparing model output before and after:

```bash
# BRF at 4 wavelengths for one view/sun geometry
echo -en "1 4 450 600 800 1000\n10 0 30 20\n" | ./gortt -LAI 4.0

# with per-component spectra {C G T Z}, viewed proportions [Kc Kg Kt Kz] and energy terms
echo -en "1 2 550 860\n0 0 30 0\n" | ./gortt -LAI 3.0 -prnspec -prnprop -energy

./gortt -u    # full option list (written to stderr)
```

Input on stdin is a header line `<n_angles> <n_wavelengths> <wl1> <wl2> ...` (nm, must be
400–2500) followed by `n_angles` lines of `vza vaa sza saa` in degrees. The counts on the
header line are checked against what is actually read and a mismatch is a fatal error. The
header line is echoed verbatim to stdout, then one output line per input angle.

Gap probabilities are the expensive part of a run and depend only on canopy/crown geometry,
so they can be cached and reused across spectral variations:

```bash
./gortt -LAI 4.0 -W < /dev/null > probs.dat            # compute probabilities and exit
echo -en "1 1 800\n10 0 30 20\n" | ./gortt -LAI 4.0 -P probs.dat -Cab 50
```
The geometry options must be identical between the `-W` and `-P` runs; only spectral options
may vary.

## Architecture

Single binary, one pass over stdin. `main()` in `gortt.c` is the whole driver.

**Data flow:** `gortt_cl_parser` → `gortt_init_params` → gap probabilities (computed, written,
or read from file) → parse stdin header → `gortt_price_soil` + `gortt_prospect_interface`
(spectra are computed once, they do not depend on angle) → per-angle loop
{ `gortt_set_zenith_dependant_probabilities` → `gortt_rsurf` → optionally `gortt_energy` }.

**Four state structs are threaded through nearly every function**, conventionally named
`p`, `g`, `s`, `c` (defined in `include/gortt.h`):
* `gortt_parameters` — crown/canopy geometry, all the precomputed probability arrays, and
  their allocations. This is the big one.
* `gortt_geometry` — one view/sun geometry plus the four viewed proportions `Kc/Kg/Kt/Kz`.
* `gortt_spectra` — wavelength grid, leaf/soil optics, PROSPECT and Price inputs, outputs
  (`rsurf`, `scomp`, `albedo`, `favegt`, `fasoil`).
* `gortt_control` — output/mode switches only.

**The scene model** (see the comment block on `gortt_rsurf`): reflectance is
`Kc*C + Kg*G + Kt*T + Kz*Z`, where `Kc/Kg/Kt/Kz` are the viewed areal proportions of sunlit
crown, sunlit background, shaded crown and shaded background, and `C/G/T/Z` are their
signatures. The `K` terms are wavelength-independent and computed once per geometry; the
signatures are computed inside a loop over wavebands. Each signature is a direct/diffuse
mixture weighted by `p->fd` (fraction of direct radiation), so most scattering functions come
in `_df` (direct) and `_ff` (diffuse) pairs.

**Prime ("_p", "_prime") variables are the ellipsoid→sphere transform.** GORT models crowns as
ellipsoids of horizontal radius `r` and vertical radius `b`; the analytic derivation is done for
spheres by scaling the vertical axis by `ellipticity = b/r`. Anything named `*_p`, `*_prime`,
or produced by `gortt_prime_theta` lives in that transformed space. Mixing primed and unprimed
quantities in one expression is almost always a bug.

**File responsibilities:**
* `gortt.c` — `main()`, CLI parsing and usage, `gortt_rsurf` (top-level scene composition),
  init/alloc/dealloc helpers, the PROSPECT and Price interfaces, soil-file reader.
* `gortt_brdf.c` — geometric-optical kernels (`gortt_kg`, `gortt_kc*`, `gortt_overlap`) and the
  radiative-transfer terms (`gortt_t_*`, `gortt_p_*`, `gortt_T_inf_*`, `gortt_R_inf_*`,
  `gortt_kuusk`, phase-function asymmetry).
* `gortt_pn_kopen.c` — the within/between-crown gap probability machinery: P(n=0), P(s), Vgamma,
  `k_open`/`k_openep`, and the geometry integrals behind them. Largely translated from
  Wenge Ni-Meister's original code; the `static` helpers at the bottom (`gortt_trisec`,
  `gortt_sector`, `gortt_triang`, `gortt_cylind`, ...) are internal to this file. Also holds
  `gortt_gap_probabilities_Q08`, the Lewis functional approximation used in Quaife et al. (2008),
  selected with `-q08_pn_kopen` — it is only valid for h=0.
* `gortt_albedo.c` — `gortt_albedo` (Gauss-Legendre hemispherical integral of the BRF, calls
  `gortt_rsurf` at each quadrature point) and `gortt_energy` (SW energy balance → albedo,
  fraction absorbed by vegetation, fraction absorbed by soil). `gauleg` lives here too.
* `gortt_lidar.c` — stub. `gortt_lidar` is called unconditionally from the per-angle loop and
  currently returns immediately; `-lidar` sets `cont.lidar_mode` but nothing reads it yet.
* `include/soil_rho.h` — **defines** (not just declares) the four Price soil EOF vectors, so it
  must only ever be included by `gortt.c`.

## Conventions and gotchas

**K&R-style function definitions are used throughout** (`void gortt_rsurf( p, g, s )` followed by
untyped-in-parens parameter declarations). Match the surrounding style when editing existing
files; prototypes in `include/gortt.h` carry the types.

**The CLI parser is prefix-matched and order-dependent.** `gortt_cl_parser` walks argv with
`strncasecmp`/`strncmp` on fixed lengths, so a shorter option would shadow a longer one — this
is why `-b` and `-r` are deliberately tested last ("leave these at the bottom!"). Adding a new
option that shares a prefix with an existing one requires care about both the length argument and
the position in the if-chain. Note `-HB`, `-BR`, `-PCC`, `-LAI` are case-*sensitive* (`strncmp`)
while most other options are not.

**Two ways to specify geometry, and the new style silently overrides the old.** Setting any of
`-HB`/`-BR`/`-PCC` makes the parser recompute `r` (fixed to 10 m), `b`, `h1`, `h2` and `lambda`
at the end of parsing, discarding any `-r`/`-b`/`-h1`/`-h2`/`-lambda` given. Similarly `-LAI`
overrides `-favd` and is converted using `lambda`, `r`, `b` — so it is applied after the
geometry conversion.

**Output ordering in `main()` matters.** `gortt_energy` → `gortt_albedo` re-runs `gortt_rsurf`
over the viewing hemisphere, mutating `g->vza`, `g->vaa` and overwriting `s->rsurf`/`s->scomp`.
The BRF and component values must therefore be printed *before* the `-energy` block, as they
currently are.

**`gortt_energy` assumes a Lambertian background.** If `G` or `Z` in `gortt_rsurf` are ever given
a BRDF, the energy balance has to change too (flagged in comments in both files).

**Known broken/incomplete paths:**
* `-soil_spectra <file>`: `gortt_read_soil_lut` (end of `gortt.c`) resamples the file correctly
  but then dumps the interpolated table to stdout and `exit(EXIT_FAILURE)`s — leftover debug
  code. The declared `gortt_get_rsoil_lut` is never defined, and `s->soil_spectra` is never
  consulted by `gortt_price_soil`, so this option does not work.
* `-diffuse <frac>` sets `p->fd = 1 - frac` and is undocumented in `-u`. Without it, `fd` is
  derived per-angle from the solar zenith following Ni et al. (1999). TODO.md flags the current
  default as "a bit dangerous".
* The leaf angle distribution is effectively hardwired: `params.lad = LAD_05` and
  `gortt_leaf_angle_distribution` returns 0.5 (and warns) for anything else.

`TODO.md` tracks the intended direction — notably breaking `main()` up and moving calculation
code out of `gortt.c`.

## Reference

Ni, W. et al. (1999), *An analytical hybrid GORT model for bidirectional reflectance over
discontinuous plant canopies*, IEEE TGRS 37(2), 987-999 — the equations the C code implements.
Quaife, T. et al. (2008), RSE 112(4), 1347-1364 — the "new style" parameterisation (HB/BR/PCC)
and the Q08 probability approximation. Price, J.C. (1990) — the soil EOFs.
