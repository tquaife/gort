# TODO list for GORT

Status as of the 2026 tidy-up (branch `tidy-up`; `v0-pre-refactor` tags the code before it).
See `OPEN-QUESTIONS.md` for things deliberately left alone and why.

## code

* ~~put most of the stuff in main() into seperate functions~~ — **done.** 413 lines to ~110,
  over eleven `static` helpers declared at the top of `gortt.c`.
* move any functions doing actual calculations out of gortt.c — **not done.** `gortt_rsurf`,
  `gortt_init_params`, `gortt_prime_theta`, `gortt_leaf_angle_distribution`,
  `gortt_set_zenith_dependant_probabilities`, the Price and PROSPECT interfaces and the
  alloc/dealloc helpers are all still there. The decomposition of `main()` was the prerequisite
  and is done, so this is now a pure-move exercise. `include/soil_rho.h` *defines* the Price
  vectors, so whatever ends up including it must be exactly one translation unit.
  * e.g. some to go in gortt_utils.c
* ~~consolidate style~~ — **mostly done.** All ~67 K&R definitions converted to ANSI, verified
  codegen-identical with `objdump`. `.editorconfig` added. Indentation was **deliberately not**
  normalised: the mixed tabs and spaces are original and reformatting would have buried the
  real changes.
* add functions to read in leaf spectra — **not done.** The soil equivalent now works
  (`-soil_spectra`, see below) and `gortt_read_soil_lut` is the obvious template; leaf optics
  are still PROSPECT-D or a pair of constants via `-refl_leaf`/`-trans_leaf`.
* ~~lidar routines~~ — **done for the nadir, single-layer case.** `gortt_lidar` implements
  Ni-Meister et al. (2001): total nadir gap probability `P(z)` from `p_n0+ epgap`, the canopy
  return `J0 rho_v dP/dz`, the ground return `J0 rho_g P(0)`, and the apparent foliage profile
  `dlogP/dz`, with an optional binned and pulse-convolved waveform. The two backscattering
  coefficients are those of Ni-Meister et al. (2018), `rho_v = 2*rleaf/3` and `rho_g = rsoil`,
  not the `2*tleaf/3` and `rsoil/2` of the 2001 paper — see the comment at the top of
  `gortt_lidar.c`. `-nlayers` went in with it
  and sets the reporting resolution of the whole model. Still open:
  * off-nadir. `gortt_set_zenith_dependant_probabilities` is hardwired to `[0]` in four places
    and would need a height-resolved variant; and the one-way attenuation of eq. 1 is only
    right at the hotspot, which nadir guarantees. Item 4 of `OPEN-QUESTIONS.md` bites here.
  * `-P`/`-W` for lidar. The file format holds the ground row only, so the combination is
    refused rather than silently returning zeros. Extending it would make cached probabilities
    usable for waveforms, which is worth doing — a lidar run is ~6x a BRF run at 15 layers and
    ~20x at 60.
  * multi-layer and multi-species canopies (paper section III-B-2).
* ~~remove compilation warnings~~ — **done for C.** Zero warnings under `-Wall`; a new one is
  something you introduced.
  * especially those from C binds added to PROSPECT — **not done.** `-Wc-binding-type` and
    `-Wtabs` still come from `PROSPECT-D/*.f90`, which is vendored code.
* Turn off diffuse radiation by default (a bit dangerous the way it's currently set up).
  — **not done.** Still derived per angle from the solar zenith unless `-diffuse` is given, and
  `-diffuse` is still undocumented in `-u`.

* ~~the quadratures behind the gap probabilities~~ — **done.** Four of them stepped with
  `nlayers`, which is the output grid, so the accuracy of the model depended on how finely you
  asked for it to be reported, and the BRF was inheriting whatever error the default fifteen
  layers gave (`k_openep` 26% high). Each now has its own resolution, scaled by `-nquad`:
  Vgamma and the mean within-crown path length are integrated in closed form except across one
  band each, the first-contact integral moved off the layer grid and onto its panel midpoints,
  and the path-length histogram was deleted outright. `p_n0` at the ground now matches its
  closed form to 2.8e-14 at every zenith angle, and the `-W` tables are byte-identical across
  any two values of `-nlayers`. See `figures/quadrature_convergence.py`. Still open:
  * `NQ_CONTACT` is the least converged of the three orders — 1.4e-5 relative on `k_openep`
    for the built-in canopy, 8e-4 for a sparse one. Raising it is linear in cost. Worth
    revisiting if anyone needs the gap probabilities to better than about 1e-4.
  * `gortt_get_s` divides the crown volume above the cut plane by the cross section *at* that
    plane, while its own limiting branch divides by the crown's full projected area. The two
    agree in the limit — the `sec(theta')` terms cancel — so this is not a bug, but the
    general branch is worth checking against the derivation.
  * `maxcrowns` overflows if raised much past 130, because `P_n` is evaluated from its closed
    form rather than recurrently. Harmless at its current 30, but a trap — item 12 of
    `OPEN-QUESTIONS.md`.

## added since

* `-O3` — deliberately **not** enabled. `CFLAGS` is `-Wall -g`, i.e. `-O0`, because at `-O0`
  there is no FP contraction or reassociation and the golden files can therefore prove a
  refactor changed nothing. Turning it on should be its own commit, with `-ffp-contract=off`,
  and it will move the numbers.
* ~~Promote the `if(TRUE)`/`if(FALSE)` alternative formulations in `gortt_overlap` and
  `gortt_kc_PP_only` to a selectable enum~~ — **done differently:** they were scored against a
  Monte Carlo, all lost, and they were deleted (5399f53). The alternatives are kept as switches
  in `figures/goms_kernels.py`. See `KC-INVESTIGATION.md`.
* ~~Resolve the `F` versus `F²` discrepancy in `gortt_kc`~~ — **done** (129ed80), together with
  `Po` (0204222) and the overlap height (59121b9). See `OPEN-QUESTIONS.md` item 1.
* Resolve β in `gortt_kc` — `OPEN-QUESTIONS.md` item 9. (Γi was corrected to the illumination
  area in c0ef2c3; the view-direction dependence is what remains.) It is now the largest error in the
  geometric-optical proportions, and fixing it will move most BRF goldens.
* Make `-h1` equal to `-h2` work, or refuse it — `OPEN-QUESTIONS.md` item 10.

## general

* ~~Add example data for testing~~ — **done.** `tests/` is a 51-case golden-file suite:
  `make test` (~3 s), `make bless` to re-capture. Covers stdout, stderr and exit status, so the
  error paths are tested too. `tests/README.md` explains the layout, the `GORTT_PREC=17`
  precision rationale, and the cases that pin known-bad behaviour.
* Documentation — **partly done.** `CLAUDE.md` covers build, test, architecture and the
  gotchas; `tests/README.md` covers the suite; `OPEN-QUESTIONS.md` covers what is still
  unresolved. There is still no user-facing documentation of the model or its options beyond
  `gortt -u`.

## fixed along the way

Memory safety (Phase 1) — a heap write overflow in `gortt_albedo` (accumulators sized by
`npoints` but indexed by waveband, so any `-energy` run with more than 32 wavebands corrupted
the heap silently), an unbounded write into the wavelength array, two out-of-range table reads
at 2500 nm, a NULL `argv[argc]` dereference, and two leaks.

Numerics (Phase 3) — the `-P` probability file was written in a *fixed* point format, so gap
probabilities below 1e-40 became zero and `-P` returned `-nan` above 89°; `gortt_calc_epgap`
and `gortt_calc_kopen` disagreed about the theta grid; the Q08 trapezium was seeded from
uninitialised array elements; `-soil_spectra` had never worked at all; and `-alb_soil 0` made
the energy balance divide 0 by 0.

Performance — `gortt -W` went from 1.41 s to 0.08 s (a provably dead subtree removed, plus
`nh_probs`), and the BRF/albedo path is about 5x faster (`gortt_rt_init` evaluates each
radiative-transfer term once per waveband instead of up to 39 times).
