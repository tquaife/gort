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
* lidar routines — **not done, but unblocked.** `gortt_lidar` is still a stub. The groundwork
  is in: `p->nh_probs` makes the full height-resolved probability tables available (`-lidar`
  sets it to `nlayers`), `gortt_parameters` is `const` throughout the BRF/albedo path so the
  compiler now enforces what a new routine may touch, and per-geometry state lives in
  `gortt_geometry` where it belongs. Read items 4 and 6 of `OPEN-QUESTIONS.md` first.
* ~~remove compilation warnings~~ — **done for C.** Zero warnings under `-Wall`; a new one is
  something you introduced.
  * especially those from C binds added to PROSPECT — **not done.** `-Wc-binding-type` and
    `-Wtabs` still come from `PROSPECT-D/*.f90`, which is vendored code.
* Turn off diffuse radiation by default (a bit dangerous the way it's currently set up).
  — **not done.** Still derived per angle from the solar zenith unless `-diffuse` is given, and
  `-diffuse` is still undocumented in `-u`.

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
* Resolve β in `gortt_kc` — `OPEN-QUESTIONS.md` item 9. It is now the largest error in the
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
