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
make test       # 107 golden-file cases, ~20 s; or: cd tests && ./run.sh
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

# -nlayers is the reporting grid and -nquad the internal quadrature orders;
# the BRF must not move with the first, and should not move with the second
./gortt -W -nlayers 47 < /dev/null | cmp - <( ./gortt -W < /dev/null )

# lidar: no stdin redirection, because lidar mode does not read stdin
./gortt -lidar -LAI 4.0                                  # per-layer profile
./gortt -lidar -lidar_waveform -lidar_pulse 1.0 -LAI 4.0 # binned waveform
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
`gortt_energy` }. `main()` is a ~110-line driver over twelve `static` helpers, all declared at
the top of `gortt.c`.

**Two options short-circuit that flow and return before stdin is read at all**: `-W`, which
writes the gap probabilities and exits, and `-lidar`, which builds a one-waveband spectral
grid from `-lidar_wl` (`gortt_lidar_spectra`, the lidar counterpart of `gortt_read_header`),
reuses `gortt_alloc_outputs`/`gortt_price_soil`/`gortt_prospect_interface` unchanged, calls
`gortt_lidar` and exits. Neither reaches the per-angle loop.

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

`gortt_lidar_params` is a fifth, but it is not threaded: `gortt_cl_parser` writes it and
`gortt_lidar` reads it, and nothing in between ever sees it, because lidar mode short-circuits
`main()` before the per-angle loop exists to thread it through. It holds the instrument
(wavelength, `J0`, pulse width, bin size, output layout), which belongs neither in `p` (the
canopy) nor in `c` (documented as output/mode switches only).

A sixth, `gortt_rt`, is a local rather than a threaded parameter: `gortt_rt_init` fills it with
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
* `gortt_pn_kopen.c` — the within/between-crown gap probability machinery: P(n=0), Vgamma,
  `epgap`, `k_open`/`k_openep`, and the geometry integrals behind them. Largely translated from
  Wenge Ni-Meister's original code, though the four quadratures have since been replaced (see
  "No quadrature in the model steps with `nlayers`" below); the `static` helpers at the bottom
  (`gortt_trisec`, `gortt_sector`, `gortt_triang`, `gortt_cylind`, ...) are internal to this
  file. Also holds
  `gortt_gap_probabilities_Q08`, the Lewis functional approximation used in Quaife et al. (2008),
  selected with `-q08_pn_kopen` — it is only valid for h=0.
* `gortt_albedo.c` — `gortt_albedo` (Gauss-Legendre hemispherical integral of the BRF, calls
  `gortt_rsurf` at each quadrature point) and `gortt_energy` (SW energy balance → albedo,
  fraction absorbed by vegetation, fraction absorbed by soil). `gauleg` lives here too.
* `gortt_lidar.c` — the lidar waveform of Ni-Meister et al. (2001). `gortt_lidar` is *not*
  called from the per-angle loop: `-lidar` short-circuits `main()` the way `-W` does, because
  the device is at nadir and no input geometry is meaningful, so stdin is never read. It takes
  no `gortt_geometry` — nadir is `theta[0] == 0` exactly, so nothing is interpolated. Two
  output layouts (per-layer profile, or binned and pulse-convolved waveform), chosen by
  `-lidar_waveform`.
* `include/soil_rho.h` — **defines** (not just declares) the four Price soil EOF vectors, so it
  must only ever be included by `gortt.c`.

## Conventions and gotchas

**Layout is idiosyncratic and deliberately unnormalised** — tabs, spaces before `;`, spaces
inside parens (`gortt_rsurf( p, g, s )`). Match the surrounding file rather than reformatting;
`.editorconfig` covers only tabs, trailing whitespace and the final newline. Function
definitions are ANSI (they were K&R until they were converted wholesale).

**Only rows `h < p->nh_probs` of `epgap` are filled.** It defaults to 1, because the BRF only
needs the gap probabilities at the ground, and computing all fifteen layers cost about 95% of
the runtime. `-lidar` raises it to `nlayers` and the full table comes back — a waveform is a
function of height. If you add something that indexes `epgap[h][t]` above `h=0`, it must raise
`nh_probs` too. (`p_n0` is the exception: it is cheap and is filled for every layer
regardless.)

**The raise happens after the argv loop, not in the `-lidar` branch.** `-nlayers` means the
two options are order-dependent, and `gortt -lidar -nlayers 30` would otherwise leave
`nh_probs` at 15 and read calloc zeros for the top half of the canopy. `gortt_lidar` has a
defensive check for exactly this.

**No quadrature in the model steps with `nlayers`.** `nlayers` is the *output* grid: how finely
the lidar profile is reported, and the `h` index of the probability tables. It is not an
accuracy knob and must not become one again. Four quadratures used to take their step from it,
and because the BRF reads only `h=0` it was silently inheriting whatever error the default
fifteen layers happened to give — `k_openep` was 26% above its converged value. They now each
have their own resolution:

| what | where | rule |
|---|---|---|
| Vgamma | `gortt_crown_proj_volume` | closed form on three of four pieces, `NQ_GL` Gauss-Legendre points across the band between the `sin t` breakpoints |
| first contact | `gortt_calc_epgap_ht` | `NQ_CONTACT` panels over `[z1_p, z2_p]`, each represented by its **midpoint** |
| path length `s` | `gortt_calc_epgap_ht` | none — `exp(-s tau')` is accumulated with the `s` that was computed |
| `es` | `gortt_get_es` | closed form outside `z_h ± r`, `NQ_GL` Gauss-Legendre points between |

All three `NQ_*` are in `include/gortt.h` beside `NH_OVERLAP`, and all are multiplied by
`p->nquad` (`-nquad`, default 1) so convergence can be shown without a recompile.
(`NQ_BRDF_PHI` sits with them but is not one of these: it is the azimuth rule for
`-rt_crown_brdf`'s hemisphere integral, it follows `-nquad` like the rest, and because it
never touches the gap probabilities the `-W` output cannot move with it.) **The
regression suite pins the decoupling**: `nlayers-noop-{a,b,c}` are the same BRF run at 15, 31
and 120 layers and their goldens must stay byte-identical to each other, and
`probs-nlayers.out` must stay byte-identical to `probs.out`.

**Two gap probabilities have closed forms, and both are worth checking against.** `z1 = h1 - b`
and `b` becomes `r` under the sphere transform, so the ground layer sits exactly one crown
radius below the lowest crown centre. Every crown is then entirely above it, which makes
`p_n0[0][t] = exp( -lv' pi r^2 sec(theta') ( h2_p - h1_p ) )` and `gortt_get_es` at `h=0` the
constant `4r/3`, both at *every* zenith angle. gortt matches the first to 2.8e-14 across
0–90°. `figures/quadrature_convergence.py` draws this, the `-nlayers` independence and the
`-nquad` convergence; run it after touching anything in `gortt_pn_kopen.c`.

**The CLI parser is prefix-matched and order-dependent.** `gortt_cl_parser` walks argv with
`strncasecmp`/`strncmp` on fixed lengths, so a shorter option would shadow a longer one — this
is why `-b` and `-r` are deliberately tested last ("leave these at the bottom!"). Adding a new
option that shares a prefix with an existing one requires care about both the length argument and
the position in the if-chain. Note `-HB`, `-BR`, `-PCC`, `-LAI` are case-*sensitive* (`strncmp`)
while most other options are not.

Two orderings are load-bearing and are commented as such: `-nlayers` and `-nquad` are tested
**first**, because `-N` is `strncasecmp`'d at two characters and would otherwise swallow both
as the PROSPECT structure parameter; and every `-lidar_*` option is tested **above** `-lidar`,
which is matched at six characters and would otherwise swallow all of them.

**An option that is inert is worse than one that errors, because `gortt` blocks on stdin.**
A `-lidar_*` option without `-lidar` used to set its field and change nothing else, so the run
fell through to the BRF path and waited on a terminal with no output and no message — which
reads as a hang, not as a usage error. The parser now refuses the combination. Anything new
that is only meaningful in one mode needs the same treatment.

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

**The `-rt_*` options are experimental and all default to off.** Most of them come
from Ni et al. (1999) section V, which says the model is validated over dense
forests and defers three things that matter as a canopy opens up. Two do not:
`-rt_crown_brdf` comes from the single-crown measurement
(`CROWN-BRDF-INVESTIGATION.md`) and `-rt_upscatter` from the closure underneath
eqs (8)-(11) rather than from anything the paper defers. They live in `p`
(group 1, beside `use_user_beta`), not in `c`, because `gortt_rt_init` and
`gortt_rsurf` never see a `gortt_control` and threading one through would touch
five signatures in the documented `const` chain. In the argv if-chain they sit
above `-r`, which is matched at two characters at the bottom and would otherwise
swallow them *and* eat the next word as a crown radius.

The invariant that matters: **with no `-rt_*` option the model is the paper as
published**, and `rt-path-scene` pins it — `-rt_crown_path -rt_path scene` asks
for the depth eqs (15)-(16) already use, so its golden must stay byte-identical
to `rt-crown-path-off`. `gortt_two_stream` was pulled out of `gortt_rt_init` so
the crown path could be evaluated at a second depth; that case is what catches an
algebra slip in the extraction. See `flight_compare/variants/` for what each one
is worth against FLIGHT.

`p->elai` is the *scene* LAI, `favd * lambda * V`, and that is correct: it is the
paper's `tau*h` with `h = lambda*V`, the vertical squeeze section III-C describes.
Verified equation by equation against `papers/Ni_etal_1999.pdf`, and gortt
reproduces the published Fig. 7 component signatures with it. Do not "fix" it.
`p->crown_lai` (`favd * 4b/3`, identically `elai/PCC`) and `p->squeeze_lai`
(`elai/(1-K_open)`) are the alternatives `-rt_crown_path` can select, and exist
only for that.

**`-rt_upscatter` is the only option that moves the red band, and `beta` means two
different things.** Eqs (8)-(11) fix the two-stream upscatter parameter at 1/2
implicitly — there is no `beta` in the published algebra, because `rt->gamma =
sqrt(1-w)` *is* that assumption: the closure has `a = 2k[1-w(1-beta)]` and
`b = 2k w beta`, whose eigenvalue collapses to `2k sqrt(1-w)` only at `beta = 1/2`.
`gortt_two_stream` therefore branches: with `-rt_upscatter` it computes `lambda`
directly from `rt->beta_up = 1/2 + (rl-tl)/(8w)` and takes `R_inf_ff` as
`b/(a+lambda)` rather than the algebraically equal `(a-lambda)/b`, which is 0/0 for
a black leaf. **Only the four semi-infinite quantities change** — the adding half of
the function is generic in them, which is what `latex_notes/two_stream_derivation.tex`
establishes. The beam upscatter is taken equal to the diffuse one; the proper one is
direction-dependent (Sellers 1985) and that is a known approximation, not a result.

`rt->beta_up` is **not** `p->beta`, which is Li & Strahler's weight between uniform
and random crown-centre heights in the geometric-optical mutual shadowing (`-beta`,
OPEN-QUESTIONS #9). Two different quantities from two different papers with the same
Greek letter, in two different halves of the model.

The option is identically inert when `rl == tl`, so it does nothing in the near
infrared of the FLIGHT comparison: `rt-upscatter-rteq{,-off}` pin that pair, which
agree to 4e-16 against 6e-2 when `rl != tl`. It cannot be a byte-identical control
because `sqrt((1-w/2)^2-(w/2)^2)` is not `sqrt(1-w)` in floating point.

**Known broken/incomplete paths:**
* `-diffuse <frac>` sets `p->user_fd = 1 - frac` and is undocumented in `-u`. Without it,
  `g->fd` is derived per-angle from the solar zenith following Ni et al. (1999). TODO.md flags
  the current default as "a bit dangerous".
* `-Anth` is advertised by `gortt -u` but the parser only matches `-canth`, so the documented
  spelling is rejected. Pinned by the `badanth` test case.
* The leaf angle distribution is effectively hardwired: `params.lad = LAD_05` and
  `gortt_leaf_angle_distribution` returns 0.5 (and warns) for anything else.
* **Two removable singularities in `gortt_two_stream` print `-nan`**, both in the published
  path and both reachable from the command line: `2*gamma*mu_i == 1` in `T_inf_df` (e.g.
  `-refl_leaf 0.45 -trans_leaf 0.30` with an overhead sun) and `omega == 1` in eqs (8)-(9)
  (`-refl_leaf 0.5 -trans_leaf 0.5`). Both have closed-form limits; OPEN-QUESTIONS #13 has
  them. Pinned as known-bad by `badtwostream-{gamma,omega1}` the way `badanth` is, so a fix
  shows up as a reviewable golden change. Not fixed in the crown BRDF work because the fix
  moves the published path and wants its own `make bless`.

**The geometric-optical kernels were checked against a Monte Carlo** (`KC-INVESTIGATION.md`),
which led to these changes:
* F² dropped from `gortt_kc`.
* `Po` per Schaaf et al. (1994) eq. 23A in `gortt_kc_fFbeta`: PiMi beyond the hotspot in
  backscatter, `MIN(PvMv, PiMi)` elsewhere. `MIN` everywhere goes negative with the sun near
  overhead; the `lowsun` test case pins that.
* The shadow overlap averaged over crown-centre heights in [h1, h2] (`NH_OVERLAP`
  Gauss-Legendre points, set up in `gortt_init_params`) instead of using the canopy depth.
* Γi in β is the illumination shadow area, and β = 1 (its limit) for an overhead sun.

The MC grid covers sza 0–60; check both ends of it. `figures/goms_kernels.py` is a numpy port of those kernels, checked against the binary
with `--selftest`, with the rejected alternatives kept as switches. `figures/kc_investigation.py`
scores any change there against the MC in seconds, because the MC is cached. Use it before
touching `gortt_kg`/`gortt_kc*`/`gortt_overlap`. β is Li & Strahler (1992, IGARSS) eq. 7 with Γi
the illumination shadow area. **β is still wrong** (OPEN-QUESTIONS #9): the MC wants it large in
backscatter and small in forward scatter, and the published formula does not depend on view
direction at all.

**`-h1` equal to `-h2` hangs** in the gap-probability code, and nearly equal gives `-nan`
(OPEN-QUESTIONS #10).

**Both the direct path and `-P` return `-nan` at exactly 90° zenith**, because `theta[nth-1]`
gives `sec(90°) = inf`. `gortt_init_params` clamps `theta[90]` to 89° to keep the probability
tables finite, but the BRF path does not. Treating the endpoint analytically is a modelling
decision, not a tidy-up.

`TODO.md` tracks the intended direction and what is already done. `OPEN-QUESTIONS.md` collects
the things found during the tidy-up that were deliberately left alone — β in `gortt_kc`,
the 90° endpoint, the hardwired leaf angle distribution, the `p->k`
evaluated at the wrong solar zenith that has to be fixed before any real LAD goes in, and
(#13) the two `-nan` singularities in `gortt_two_stream`. Read it before starting the lidar
work.

**`CROWN-BRDF-INVESTIGATION.md` is the record of the crown BRDF work and is the one to
read before touching `gortt_rt_init`, `gortt_crown_reflectance` or `gortt_kuusk`.** It
supersedes `CROWN-BRDF-NEXT.md`, which was the handover that set the work up. Section 3
of that handover — measure a single crown against FLIGHT before writing model code — was
done (`flight_compare/run_crown_image.py`) and is decisive: crowns are strongly
limb-darkened when optically thin and uniform when thick, which rules out putting a BRDF
on the crown envelope; and the sunlit crown is close to Lambertian *except* for a hotspot
spike GORT under-predicts by up to 1.8x and a forward-scatter limb it does not have at
all. Its recommended fix, `-rt_crown_brdf`, was implemented and **makes the FLIGHT
agreement worse** — it is kept off by default for the negative result, as
`-rt_recollision` is. The next target is the hotspot, not the phase function.

`CROWN-BRDF-NEXT.md` is still worth reading for how it framed the problem — `-rt_crown_path`
(v3.1) fixed the crown's optical *depth*, leaving `rho'_df` a directional-hemispherical
reflectance from a plane-parallel slab where a bidirectional one is needed — and for the
order it put the work in, which was the right order. Its two candidate routes are both now
closed: Option B by the measurement, Option A by having been built and scored.
`KC-INVESTIGATION.md` is the precedent both records follow.

**`latex_notes/` is untracked** (`.gitignore` excludes the whole directory) but is where the
physics is written up, and a later paper is meant to come out of it:
`crown_brdf_notes.tex` is this work with the equations and four figures,
`two_stream_derivation.tex` reconstructs eqs (8)-(11) — which Ni et al. (1999) defers to a
reference that was never published — and establishes that they are exact adding relations in
the four semi-infinite quantities, `crown_path_notes.tex` is v3.1, `lidar_notes.tex` the
lidar. `references.bib` and `paper_summaries.tex` index everything in `papers/`. Each `.tex`
has a `*_figs.py` beside it that regenerates its figures from the `flight_compare` CSVs.

## Reference

Ni, W. et al. (1999), *An analytical hybrid GORT model for bidirectional reflectance over
discontinuous plant canopies*, IEEE TGRS 37(2), 987-999 — the equations the C code implements.
Ni, W. et al. (1997), *Transmission of solar radiation in boreal conifer forests: measurements
and models*, JGR 102(D24), 29555-29566 — the source of the default direct beam fraction
`mu/(mu+0.09)` (its eq. 5, fitted to BOREAS clear-day measurements); `-rt_fd_paper` selects the
`mu/(mu+0.05)` of Ni et al. (1999) §IV instead, and `-diffuse` overrides both.
Quaife, T. et al. (2008), RSE 112(4), 1347-1364 — the "new style" parameterisation (HB/BR/PCC)
and the Q08 probability approximation. Price, J.C. (1990) — the soil EOFs.
