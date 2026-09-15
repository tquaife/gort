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

The C side builds clean under `-Wall`; a new warning is something you introduced. Expected
build noise (all from the vendored Fortran, pre-existing): `-Wc-binding-type` from the
`bind(C)` interface in `PROSPECT-D/prospect_DB.f90`, and `-Wtabs` from both PROSPECT sources.

`CFLAGS` is `-Wall -g`, i.e. `-O0`. That is deliberate: at `-O0` there is no FP contraction or
reassociation, so hoisting a repeated call to a pure function is guaranteed bit-identical and
the golden files can be used to prove a refactor changed nothing. Turning optimisation on
should be its own commit (with `-ffp-contract=off`).

## Testing

```bash
make test       # 41 golden-file cases, ~3 s; or: cd tests && ./run.sh
make bless      # re-capture goldens -- only on a commit MEANT to change the numbers
```

`tests/` is a golden-file regression suite: `gortt` is a deterministic stdin→stdout filter, so a
byte comparison of its output is a complete check of the numerical path. Cases are registered in
`tests/cmds`; stdout, stderr and exit status are all captured. A failure prints the largest
absolute and relative deviation via `tests/numdiff.awk`, so a last-digit rounding artefact is
distinguishable from a broken model. See `tests/README.md`.

**`run.sh` sets `GORTT_PREC=17`.** Without it `gortt` prints `%f` — 4 significant figures — which
is too coarse to catch a small regression: a 1e-7 perturbation of `favd` was visible in 1 of 30
BRF cases at `%f`, and 25 of 36 at full precision. The `-W` cases print `%.17g` and are the
most sensitive probe in the suite.

`badanth` deliberately pins *known-bad* behaviour so that fixing it shows up as a reviewable
change, and `useprobs-{extreme,grazing}` hold `-P` to reproducing the direct path exactly —
see `tests/README.md`.

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
{ `gortt_set_geometry` → `gortt_run_geometry` → `gortt_print_result` → optionally
`gortt_energy` }. `main()` is a ~110-line driver over eleven `static` helpers, all declared at
the top of `gortt.c`.

**Four state structs are threaded through nearly every function**, conventionally named
`p`, `g`, `s`, `c` (defined in `include/gortt.h`). They are separated by *lifetime*, and that
separation is the load-bearing invariant of the codebase — keep it when adding anything:
* `gortt_parameters` (`p`) — the canopy, fixed for the whole run. Three bannered groups:
  user inputs (written by `gortt_cl_parser`), derived constants (written by
  `gortt_init_params`), and the gap probability tables. **Const everywhere downstream**; if
  you find yourself wanting to write to `p` from the BRF, albedo or lidar path, the field
  belongs in `g` instead.
* `gortt_geometry` (`g`) — one view/sun geometry and everything derived from it: the primed
  zeniths, `k_vza`, the direct fraction `fd`, the interpolated gap probabilities, and the four
  viewed proportions `Kc/Kg/Kt/Kz`. Rewritten per input angle; nothing survives between them.
* `gortt_spectra` (`s`) — wavelength grid, leaf/soil optics, PROSPECT and Price inputs,
  outputs (`rsurf`, `scomp`, `albedo`, `favegt`, `fasoil`).
* `gortt_control` (`c`) — output/mode switches only.

A fifth, `gortt_rt`, is a local rather than a threaded parameter: `gortt_rt_init` fills it with
every radiative-transfer term for one (geometry, waveband) pair, in the order Ni et al. derive
them. `gortt_rsurf`'s waveband loop then reads `rt.t_prime_ff` and friends as plain values.
Add new RT terms there, not as another function that recomputes its own dependencies.

**The scene model** (see the comment block on `gortt_rsurf`): reflectance is
`Kc*C + Kg*G + Kt*T + Kz*Z`, where `Kc/Kg/Kt/Kz` are the viewed areal proportions of sunlit
crown, sunlit background, shaded crown and shaded background, and `C/G/T/Z` are their
signatures. The `K` terms are wavelength-independent and computed once per geometry; the
signatures are computed inside a loop over wavebands. Each signature is a direct/diffuse
mixture weighted by `g->fd` (fraction of direct radiation), so most RT terms come in `_df`
(direct) and `_ff` (diffuse) pairs.

**Prime ("_p", "_prime") variables are the ellipsoid→sphere transform.** GORT models crowns as
ellipsoids of horizontal radius `r` and vertical radius `b`; the analytic derivation is done for
spheres by scaling the vertical axis by `ellipticity = b/r`. Anything named `*_p`, `*_prime`,
or produced by `gortt_prime_theta` lives in that transformed space. Mixing primed and unprimed
quantities in one expression is almost always a bug.

**File responsibilities:**
* `gortt.c` — `main()`, CLI parsing and usage, `gortt_rsurf` (top-level scene composition),
  init/alloc/dealloc helpers, the PROSPECT and Price interfaces, soil-file reader.
* `gortt_brdf.c` — geometric-optical kernels (`gortt_kg`, `gortt_kc*`, `gortt_overlap`),
  `gortt_kuusk`, and `gortt_rt_init`, which is the whole radiative-transfer chain in one
  function.
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
  currently returns immediately. `-lidar` sets `cont.lidar_mode`, which nothing reads yet, and
  `p->nh_probs = p->nlayers`, which does have an effect: see below.
* `include/soil_rho.h` — **defines** (not just declares) the four Price soil EOF vectors, so it
  must only ever be included by `gortt.c`.

## Conventions and gotchas

**Layout is idiosyncratic and deliberately unnormalised** — tabs, spaces before `;`, spaces
inside parens (`gortt_rsurf( p, g, s )`). Match the surrounding file rather than reformatting;
`.editorconfig` covers only tabs, trailing whitespace and the final newline. Function
definitions are ANSI (they were K&R until they were converted wholesale).

**Only rows `h < p->nh_probs` of the height-indexed probability tables are filled.** It
defaults to 1, because the BRF only needs the gap probabilities at the ground, and computing
all fifteen layers cost about 95% of the runtime. `-lidar` raises it to `nlayers` and the full
table comes back — a waveform is a function of height. If you add something that indexes
`p_n0[h][t]`, `epgap[h][t]` or `pd_s[h][t][s]` above `h=0`, it must raise `nh_probs` too.

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
over the viewing hemisphere, overwriting `s->rsurf` and `s->scomp`. The BRF and component
values must therefore be printed *before* the `-energy` block, as they currently are.
(`gortt_albedo` integrates over a *copy* of `g`, so the caller's geometry is safe; it is only
the spectra that are clobbered.)

**`gortt_energy` assumes a Lambertian background.** If `G` or `Z` in `gortt_rsurf` are ever given
a BRDF, the energy balance has to change too (flagged in comments in both files) — and so does
its reading `G`/`Z` back out of `s->scomp` *after* `gortt_albedo` has overwritten it, which is
safe today only because neither depends on the view direction.

There are **three sources of soil reflectance**, selected in `gortt_price_soil` and mutually
exclusive in the parser: `-alb_soil` (a constant), `-soil_spectra <file>` (resampled onto a 1 nm
grid by `gortt_read_soil_lut`, then interpolated), and otherwise the Price EOFs with the `-rsl*`
weights.

**Known broken/incomplete paths:**
* `-diffuse <frac>` sets `p->user_fd = 1 - frac` and is undocumented in `-u`. Without it,
  `g->fd` is derived per-angle from the solar zenith following Ni et al. (1999). TODO.md flags
  the current default as "a bit dangerous".
* `-Anth` is advertised by `gortt -u` but the parser only matches `-canth`, so the documented
  spelling is rejected. Pinned by the `badanth` test case.
* The leaf angle distribution is effectively hardwired: `params.lad = LAD_05` and
  `gortt_leaf_angle_distribution` returns 0.5 (and warns) for anything else.

**One numerical bug is pinned as known-bad rather than fixed**, so that fixing it shows up as a
reviewable change in the goldens: `gortt_calc_epgap` fills `t < nth-1` while `gortt_calc_kopen`
integrates to `nth-1`, so the last trapezium panel is built from a zero. It needs checking
against the papers before being touched.

**The geometric-optical kernels were checked against a Monte Carlo** (`KC-INVESTIGATION.md`),
which led to three changes: F² dropped from `gortt_kc`, `Po = MIN(PvMv, PiMi)` in
`gortt_kc_fFbeta`, and the shadow overlap averaged over crown-centre heights in [h1, h2]
(`NH_OVERLAP` Gauss-Legendre points, set up in `gortt_init_params`) instead of using the canopy
depth. `figures/goms_kernels.py` is a numpy port of those kernels, checked against the binary
with `--selftest`, with the rejected alternatives kept as switches. `figures/kc_investigation.py`
scores any change there against the MC in seconds, because the MC is cached. Use it before
touching `gortt_kg`/`gortt_kc*`/`gortt_overlap`. **β is still wrong** (OPEN-QUESTIONS #9): the MC
wants it large in backscatter and small in forward scatter, and gortt's formula is symmetric.

**`-h1` equal to `-h2` hangs** in the gap-probability code, and nearly equal gives `-nan`
(OPEN-QUESTIONS #10).

**Both the direct path and `-P` return `-nan` at exactly 90° zenith**, because `theta[nth-1]`
gives `sec(90°) = inf`. `gortt_init_params` clamps `theta[90]` to 89° to keep the probability
tables finite, but the BRF path does not. Treating the endpoint analytically is a modelling
decision, not a tidy-up.

`TODO.md` tracks the intended direction and what is already done. `OPEN-QUESTIONS.md` collects
the things found during the tidy-up that were deliberately left alone — β in `gortt_kc`,
the 90° endpoint, the hardwired leaf angle distribution, and the `p->k`
evaluated at the wrong solar zenith that has to be fixed before any real LAD goes in. Read it
before starting the lidar work.

## Reference

Ni, W. et al. (1999), *An analytical hybrid GORT model for bidirectional reflectance over
discontinuous plant canopies*, IEEE TGRS 37(2), 987-999 — the equations the C code implements.
Quaife, T. et al. (2008), RSE 112(4), 1347-1364 — the "new style" parameterisation (HB/BR/PCC)
and the Q08 probability approximation. Price, J.C. (1990) — the soil EOFs.
