# flight_compare

`gortt` against **FLIGHT**, the Monte Carlo photon-tracing model of North (1996),
over canopies the two models are told about in the same terms. Bidirectional
reflectance at 660 and 860 nm, and the nadir lidar waveform at 1064 nm, for
twelve scenes spanning sparse/dense cover, low/high leaf area and oblate/
spherical/prolate crowns.

`gortt` had been checked against published figures and against a Monte Carlo of
its *own* geometric-optical kernels (`KC-INVESTIGATION.md`). FLIGHT is the first
independent radiative transfer model it has been put next to.

Open `index.html` for the results. Everything here is generated; nothing is
hand-edited.

```sh
./build_flight.py     # patched FLIGHT build, outside both repositories
./validate.py         # the two models really do have the same canopy
./run_brdf.py         # ~2 h on four cores
./run_lidar.py        # ~35 min
./plot_brdf.py && ./plot_lidar.py && ./make_html.py
```

`config.py` is the single source of truth for the scene matrix, the optical
properties and the angle grids. Change it there and rerun.

## The FLIGHT tree is never written to

`/home/tquaife/src/flight` is read-only throughout, and no FLIGHT source is
copied into this repository. FLIGHT's waveform lidar is disabled at compile
time — `#define LMODE_PERMITTED 0` at `flight.c:233`, which also pins
`MAX_NO_WAVE_BINS` to 1, so it is genuinely a build-time decision — and
`check_data()` exits if you ask for it. `build_flight.py` therefore copies
`flight.c` into a scratch directory outside both repositories, flips that one
macro, and builds with `-O3`. It refuses to build inside either tree, and checks
that the copy differs from the original by exactly one line. The same binary is
used for the BRDF runs, so there is one build behind every FLIGHT number here.

Each run gets its own working directory under `work/` (gitignored), because
FLIGHT's `FPATH` is `"./"` and every path it uses is relative to the current
directory. `DATA/` is symlinked to FLIGHT's own.

## How the two models were matched

| | gortt | FLIGHT |
|---|---|---|
| crowns | spheroids, semi-axes (r, r, b) | `crowns.data`, `type 1`, `Exy`, `Ez` |
| positions | Poisson in plan at `lambda`, centres uniform in [h1, h2] | the explicit list we generate |
| foliage | `-favd` (m<sup>-1</sup>) | per-crown column 7, with `TOTAL_LAI` negative |
| leaf angles | `LAD_05`, G = 0.5 | spherical LAD, `cos(θ_lo) − cos(θ_hi)` per 10° bin |
| leaf optics | `-refl_leaf`, `-trans_leaf` (PROSPECT off) | `SPEC/leaf.spec` |
| soil | `-alb_soil` (Price off), Lambertian | `SPEC/soil.spec`, `SOILROUGH 0` |
| illumination | `-diffuse 0` ⇒ `fd = 1` | `AER_OPT < 0` ⇒ `diffuse_est()` returns exactly 0 |
| hotspot | `raa = saa − vaa`, 0 = backscatter | `VIEW_AZIMUTH − SOLAR_AZIMUTH`, 0 = backscatter |

Crowns are written out explicitly rather than left to FLIGHT's `forest_gen3()`,
which places a hardcoded 499 of them and derives the scene size from the
requested cover. An explicit list gives GORT's canopy exactly, and lets a scene
realisation be reproducible: FLIGHT's own RNG is seeded fixed unless `argv[1]`
is `random`, and `random` is `srand(time(NULL))`, which collides between
processes started in the same second. Ensemble members are therefore different
*scenes* from our own seeded generator, with FLIGHT left on its fixed seed.

**`-favd`, not `-LAI`.** gortt's `-LAI` converts with
`favd = LAI / (lambda (4/3) pi r² b)` (`gortt.c:1435`), which counts crown
volume without correcting for overlap. Giving both models `favd` makes the
foliage density identical by construction, and leaves the overlap treatment as
something this comparison measures rather than assumes. Each scene's nominal
LAI is quoted as a label; FLIGHT reports its own overlap-corrected scene LAI,
which is lower.

### Three things that had to be worked around

**FLIGHT does not replicate crowns across its periodic walls.** Photons wrap
(`replace_if_out`), but the part of a crown that lies outside the box simply is
not there. With 2 m crowns in a 46 m half-scene that loses a couple of percent
of the foliage — a systematic deficit against GORT's infinite canopy that would
read as a model difference. `scene.tile_ghosts()` adds the periodic images.
Before it, `dense-low-sphere` gave a fractional cover of 0.465 against the 0.478
that `1 − exp(−λπr²)` predicts; with it, 0.472.

**FLIGHT's lidar mode double-frees.** `read_active_source()` closes the
`solar.spec` handle at `flight.c:6422` and closes the same `FILE*` again at
`flight.c:6444`. With afterpulsing off nothing reopens it in between, so every
lidar run aborts with `free(): double free detected in tcache 2` before tracing
a photon. Rather than patch more source, the harness switches afterpulsing *on*
and supplies an `afterpulse.spec`: the second `fclose()` then lands on a file we
control. This is presumably why the lidar mode ships disabled.

**A zero afterpulse energy fraction does not switch the afterpulse off.**
`flight.c:6954` draws the selector as `r = rand() * RAND_MAX`, an `int * int`
product that overflows, so `r` is negative about half the time and the test
`r < g_afterpulse_energy_frac[i]` on the next line passes about half the time
whatever the fraction says. Measured: set the pulse width to 1e-4 ns and leave
the afterpulse sd at 1 ns and 43% of the return picks up the afterpulse spread;
set the afterpulse sd to 1e-4 ns as well and the whole return collapses into a
single bin. The harness therefore writes the afterpulse to be *indistinguishable
from the real pulse* — sd equal to `LIDAR_PULSE_NS`, delay zero — rather than
writing it to be inert. Whichever branch fires, the photon gets the same
`gaussrand() * 0.3 * LIDAR_PULSE_NS` metres of path, and the ground-return
standard deviation is the 0.15 m the 1 ns pulse asks for. Anyone reusing
`scene.write_lidar()` with a different pulse width must keep the two in step.

**gortt's lidar layers have to be set from a thickness, not a count.** gortt
spreads `-nlayers` layers over `z1..z2` (`gortt.c:987`), and
`gortt_lidar_print_waveform` drops the whole of each layer-to-layer increment of
P into the *nearest* range bin (`gortt_lidar.c:197`). Once `dz` exceeds the bin
size the bins alternate between receiving one layer and receiving none, and the
waveform comes out as a staircase. The prolate crowns span 14 m against the
oblate crowns' 6 m, so at a fixed 120 layers they sat at 0.85 layers per 0.1 m
bin and were the only visibly stepped scenes. Measured raggedness (rms second
difference of the normalised canopy return, relative to its peak) against layer
count for `sparse-low-prolate`:

| nlayers | 120 | 180 | 240 | 360 | 480 | 1200 |
|---|---|---|---|---|---|---|
| layers per bin | 0.85 | 1.28 | 1.71 | 2.56 | 3.42 | 8.56 |
| raggedness | 0.037 | 0.016 | 0.011 | 0.007 | 0.006 | 0.003 |

Two things came out of that. **gortt itself was fixed** (gort commit `0ce47c1`):
`gortt_lidar_deposit` now splits each contribution between the two bins it falls
between instead of dropping it whole into the nearest one, which is exactly
conservative in the total energy and in its first moment. The prolate raggedness
at 120 layers fell from 0.037 to 0.007 — what 360 layers used to buy — and the
ground delta now lands on z = 0 exactly rather than up to half a bin away.

**And `config.LIDAR_DZ = 0.02` m fixes the layer *thickness* rather than the
count**, giving five layers per bin in every scene (301 layers oblate, 401
sphere, 701 prolate). That is no longer about aliasing but about comparability:
`z2 − z1` is 6 m for the oblate crowns and 14 m for the prolate ones, so a fixed
layer count resolves them differently and the per-scene RMSDs against FLIGHT
would not be measuring the same thing. The prolate waveforms went from the
roughest at 0.037 to the smoothest at 0.001–0.003; what is left in the oblate
scenes, 0.005–0.009, is the genuine curvature of a profile squeezed into 6 m.
It changed none of the reported statistics — `canopy_shape_rmsd` is identical to five significant
figures in all twelve scenes and the mean return heights moved by at most
0.004 m — because the staircase oscillated about the right curve rather than
displacing it. The cost is linear in the layer count and the whole gortt lidar
pass is about three minutes.

**`flightline.data`'s scattering-mode field does nothing.** It is read into
`g_lidar_mode` at `flight.c:6465` and never referenced again, so single
scattering is selected with `BOUNCE_THRESHOLD 1` in `reverse.data` instead.

## Validation

`validate.py` runs before the sweep and checks the models against each other on
quantities that do not depend on any of the radiative transfer:

* **Fractional cover** — our own Monte Carlo of the crown list, FLIGHT's
  `find_frac_cover()`, and `1 − exp(−λπr²)`, agreeing to better than 0.01.
* **Nadir gap probability** — black leaves (reflectance = transmittance = 0) and
  a white Lambertian soil, viewed and illuminated at nadir. FLIGHT's scene-mean
  BRF is then exactly the probability a vertical ray reaches the ground, and is
  directly comparable with `p_n0[0][0] + epgap[0][0]` from `gortt -W`. **The two
  agree to within 0.8%** across sparse, dense, spherical and prolate scenes.
  That is the decisive check that `favd`, `lambda`, `r`, `b`, `h1` and `h2` mean
  the same thing in both, and it settles a question the LAI diagnostics cannot:
  crown overlap attenuates twice in FLIGHT, as GORT's `elai` assumes.
* **Hotspot azimuth** — both models brighter at `raa = 0` than at `raa = 180`
  when `vza = sza`, so the convention mapping is pinned against the models
  rather than against a reading of the source.

## What it shows

Headline numbers, over all 12 scenes, 3 solar zeniths and 64 view directions,
as bias and RMSD relative to the mean FLIGHT BRF:

| | red (660 nm) | NIR (860 nm) |
|---|---|---|
| **all directions** | **+0.2% bias, 15% RMSD** | **−7.4% bias, 19% RMSD** |
| view zenith 0–20° | −4.6%, 6.6% | −9.7%, 14% |
| view zenith 20–40° | −3.0%, 7.7% | −7.6%, 15% |
| view zenith 40–60° | +0.6%, 16% | −6.6%, 19% |
| view zenith 60–75° | +9.5%, 27% | −6.6%, 23% |

The overall bias is small; it is the *shape* of the BRDF that differs.

**Agreement degrades with view zenith, in both bands.** Near nadir gortt and
FLIGHT are within a few percent in the red. Beyond 60° the red RMSD is 27% and
gortt reads high. gortt's polar plots have a markedly stronger bright limb than
FLIGHT's, which are comparatively flat.

**Crown shape is the strongest control.** In the red, RMSD is 7.6% for oblate
crowns, 11% for spheres and 24% for prolate. The worst single scene is
`dense-low-prolate` (b/r = 2.5, low foliage density in tall crowns), at 34–52%
RMSD in the red with the correlation across directions collapsing to −0.25 at
a solar zenith of 15°. gortt puts a strong hotspot and a bright limb on a scene
that FLIGHT renders almost Lambertian.

**In the red, sparse canopies are fine and dense ones are not**: 5.1% RMSD
sparse against 24.5% dense. In the NIR the split is by leaf area instead —
gortt reads 16% low at low LAI and is unbiased at high LAI — which points at
the multiple-scattering treatment rather than at the geometry.

**The forward-scatter half is the worse one in the NIR**: −10.8% bias and 24%
RMSD beyond 20° view zenith, against −4.6% and 16% in backscatter. The
principal-plane figures show why: gortt's forward limb turns back up into a
symmetric U, where FLIGHT falls away monotonically. That is the signature of
`OPEN-QUESTIONS.md` item 9 — β in `gortt_kc` has no view-direction dependence —
reached here independently of the Monte Carlo in `KC-INVESTIGATION.md`, which
is of gortt's own kernels rather than of a separate model.

### Lidar

**The waveform shape is the best-agreeing thing in this comparison.** With the
ground return removed and each profile renormalised, the canopy shapes differ
by an RMSD of 0.009–0.020 on a profile that integrates to one, and the mean
return height agrees to 0.5 m for the high-LAI scenes. gortt's mean return
height is biased low — by 1.2–1.4 m in the low-LAI prolate canopies, where it
puts the return too far down the crown.

**The ground-to-canopy energy split agrees — after this comparison turned up an
error in the published model.** The nadir gap probability was validated between
the two models to within 0.8% before any of this was run, so the partition can
only come from the backscatter coefficients. Backing the implied ratio out of
each model isolated it:

| | implied ρ_g/ρ_v | |
|---|---|---|
| gortt, as first written | 0.426 exactly, every scene | `= 0.5 r_soil / (2 t_leaf/3)` |
| gortt, Ni-Meister et al. (2018) | 0.750 exactly, every scene | `= r_soil / (2 r_leaf/3)` |
| FLIGHT | 0.79 mean (0.61–0.92) | |

gortt reproduces `f_ground/(1−f_ground) = (ρ_g/ρ_v)·P(0)/(1−P(0))` to the
digit, so the arithmetic was never in doubt and the whole difference was in that
constant. This page recorded it as a question for Ni-Meister et al. (2001)
rather than something to change on the evidence here. That was the right call
and the answer came back yes: **Ni-Meister et al. (2018), RSE 204, 229-243
corrects both coefficients** — `ρ_v = 2 r_leaf/3`, from the leaf *reflectance*,
because the hotspot is at scattering angle ξ = π and not 0; and `ρ_g = r_soil`,
not half of it. gortt was changed to match (see `gortt_lidar.c`).

The result is a genuine agreement rather than a disclaimer: the ground fraction
is now within 7% of FLIGHT in every scene, median ratio 0.97 against 0.84
before, and 0.750 sits comfortably inside FLIGHT's 0.61–0.92 scene-to-scene
spread. The canopy shape RMSDs are unchanged to four decimal places, as they
must be — the shape is `dP/dz` alone.

This is independent of the evidence in the 2018 paper itself, which recovers
ρ_v/ρ_g ≈ 1 from LVIS waveforms and 0.77–0.84 from RAMI-IV field data.

**Multiple scattering shows up as a tail below the ground, not as a broader
canopy.** The `orders 1-8` curve is FLIGHT's *cumulative* waveform over
scattering orders 1 to 8 — its lidar loop runs the order from 1 to
`BOUNCE_THRESHOLD` and accumulates a contribution for each (`flight.c:6879`) —
so it is the single-scattering curve with more physics in it, not order 8
alone. A multiply scattered photon takes a longer route to return from the same
height, and a longer route is recorded as a lower height, so the ground peak
grows a tail *downwards*: the standard deviation of the sub-1 m return goes from
0.18 m at single scattering to 0.8–1.8 m over orders 1–8, and 1.3–9.2% of the
total return arrives below z = −1 m, running to the floor of the recording
window at −10 m. It is largest in the dense, low-LAI canopies (9%), where there
is plenty of crown to bounce off and plenty of gap to escape through. The
canopy shape above z1 barely moves. gortt has no counterpart, Ni-Meister et al.
(2001) being a single-scattering waveform, which is the point of plotting it.

**FLIGHT's ground return sits at −0.198 m, not at 0**, in every scene, and that
is two separate things:

* −0.15 m of *range walk*, which is real. The source is a point 1000 m up, not a
  collimated beam, so a photon entering the 25 m footprint off-axis travels
  `sqrt(range² + ρ²)`. Over a uniform disc that is `E[ρ²]/(2·range)` = 0.156 m
  of extra one-way path. Shrinking the footprint to 1 m removes it (measured:
  −0.198 → −0.052 m) and so does moving the source to 10 km. Being one-sided it
  also skews the peak, which is why it has a longer tail below than above.
* −0.05 m of fixed offset in FLIGHT's range origin, which is not. With a
  pinhole footprint, a near-zero pulse and no crowns the entire return lands in
  one bin whose true height is −0.125 ± 0.025 m — about 0.25 m of two-way path —
  and it is invariant to canopy top (tested at 5, 20 and 50 m), footprint,
  sensor range, pulse width and bin count. So it is a bookkeeping constant, not
  physics, and it is the same size as the 0.205 m gap FLIGHT leaves between the
  photon launch plane (`canopy_top + 0.01`) and its sky plane
  (`canopy_top + 0.215`). Truncating rather than rounding the bin index cancels
  half a bin of it back.

Neither affects a shape comparison — both are inside one 0.1 m bin — but they
are why the two models' ground spikes do not coincide to the bin. gortt's own
ground return is now centred on zero exactly: its grid is anchored on `z2 + 3σ`
so the delta rarely falls on a bin centre, but the linear deposition added in
gort commit `0ce47c1` — which this comparison prompted — splits it between the
two bins either side in proportion and so preserves the centroid.

## Experimental scattering variants

The near-infrared bias above is the paper's approximation reaching the edge of
its range, not a bug: `gortt` was checked equation by equation against
`../papers/Ni_etal_1999.pdf` and reproduces the Fig. 7 component signatures of
the paper's own SOBS canopy. Section V says the model has been validated over
dense forests and defers three things that matter as a canopy opens up, and
`gort`'s `-rt_*` options switch them on one at a time. **All are off by default**
— with none of them given the model is the paper as published, and the original
82 regression goldens are byte-identical.

`variants/index.html` is the contents page; one page per option, each showing
the baseline, that option, and the same FLIGHT Monte Carlo. **FLIGHT was not
re-run for any of it.** Bias and RMSD against FLIGHT, as a percentage of
FLIGHT's mean:

| variant | red bias | red RMSD | NIR bias | NIR RMSD | NIR bias, sparse | NIR RMSD, sparse |
|---|---|---|---|---|---|---|
| baseline | +0.2% | 14.9% | −7.4% | 18.6% | −9.5% | 18.4% |
| **crown-path** | +1.0% | 15.2% | **+1.7%** | **16.7%** | **+2.5%** | 13.9% |
| upscatter | +1.2% | 15.1% | −7.4% | 18.6% | −9.5% | 18.4% |
| crown-path-up | +1.9% | 15.5% | +1.7% | 16.7% | +2.5% | 13.9% |
| crown-path-bg | +1.2% | 15.2% | +5.9% | 17.8% | +9.4% | 16.4% |
| squeeze-path | +1.5% | 15.9% | +6.6% | 19.3% | +1.6% | **13.3%** |
| recollision | +1.1% | 16.5% | +6.4% | 20.3% | +4.2% | 14.8% |
| crown-brdf | +2.7% | 19.6% | +3.1% | 21.9% | +4.0% | 18.6% |
| bg-coupling | +0.5% | 14.8% | −4.0% | 18.0% | −4.3% | 17.4% |
| crown-crown | +0.5% | 14.9% | +13.7% | 37.4% | −0.9% | 22.1% |
| all | +3.2% | 20.8% | +44.8% | 68.0% | +26.9% | 36.4% |

**`-rt_crown_path` is the one that works.** NIR bias −7.4% → +1.7% with the RMSD
down as well, which a rescaling could not do — it is correcting a shape error.

### The upscatter parameter: right band, right direction, wrong aggregate

`-rt_upscatter` replaces the `beta = 1/2` that eqs (8)-(11) assume implicitly
with `beta = 1/2 + (rl - tl)/(8w)`, the standard canopy value for the spherical
leaf angle distribution `gortt` is hardwired to. It is the only option here that
is aimed at the red, and it is **exactly inert in the near infrared** — `rl = tl`
there, so `beta` is 1/2 and the near-infrared columns above are identical to the
rows they sit under, to 4e-16.

On the 12-scene aggregate it is a small net loss, 14.9% to 15.1%. That number is
misleading, and the breakdown says why:

| red RMSD, subset | baseline | upscatter | crown-path | crown-path-up |
|---|---|---|---|---|
| all 12 | 14.87% | 15.12% | 15.23% | 15.55% |
| less `dense-low-prolate` | 9.66% | 9.63% | 9.85% | 9.90% |
| less all three `dense-low` | 7.04% | **6.54%** | 6.72% | **6.23%** |

**Nine of the twelve scenes improve** (median −0.26 points). The three that do not
are exactly the three `dense-low` canopies, and they carry the aggregate because
they are the darkest scenes in the red and so have the smallest denominator —
`dense-low-prolate` alone is at 41% RMSD before anything is switched on.

The reason is not obscure, and it joins up with the single-crown measurement in
`CROWN-BRDF-INVESTIGATION.md`. A larger `beta` sends more of what is scattered
backwards, so crowns get brighter. Measured against FLIGHT with an isolated crown
over a black soil, GORT's red sunlit-crown signature is **3.4x too bright** for
`dense-low-prolate` and 1.2x for `dense-low-sphere`, and **0.35 to 0.46 times too
dark** for every scene with an optically larger crown. Brightening helps the
second group and hurts the first, which is precisely the split in the table.

So the upscatter parameter is doing what the derivation says it should — the
bias moves toward zero nearly everywhere, −4.48% to −3.60% over the nine
well-behaved scenes, and −3.98% to −3.06% on top of `-rt_crown_path` — and the
aggregate is being set by a pre-existing red overestimate in optically thin
crowns that has nothing to do with it. Fixing that is a separate problem, and
until it is fixed this option cannot show a clean win on this scene matrix.

`crown-path-up` is the two together. They do stack, unlike `crown-path-bg` below:
they act on different bands, so the near infrared is `crown-path`'s number
unchanged and the red is `upscatter`'s.

### The two that work do not stack

`crown-path-bg` is `-rt_crown_path -rt_bg_coupling`, the only two options that
each improve on the baseline, run together. It is **worse than `-rt_crown_path`
alone** — NIR RMSD 16.7% → 17.8%, and 13.9% → 16.4% in the sparse scenes, which
is where both of them do their work.

The reason is not an interaction. In near-infrared bias the two corrections are
very nearly additive:

| subset | baseline | +crown_path | +bg_coupling | both | sum of the two | excess |
|---|---|---|---|---|---|---|
| sparse | −9.5 | +2.5 | −4.3 | +9.4 | +7.7 | +1.7 |
| dense | −5.7 | +1.0 | −3.8 | +3.1 | +2.9 | +0.2 |
| all | −7.4 | +1.7 | −4.0 | +5.9 | +5.1 | +0.9 |

They add to within 2 points of independence — and that is the problem. The
sparse deficit is 9.5 points and the two corrections between them supply 17, so
applying both overshoots by as much as the baseline undershot. They are not two
complementary fixes; they are two descriptions of the same missing brightness,
and `-rt_crown_path` is the better of the two accounts of it.

Worth knowing for anything built on `-rt_crown_path`: its success in the sparse
near infrared is not evidence that the deficit is purely a crown optical depth
problem, because a background–canopy coupling term with no crown geometry in it
at all recovers half the same deficit.

The variant is kept for that negative result, as `crown-brdf` and
`-rt_recollision` are.
In the sparse scenes, −9.5% → +2.5% and the RMSD 18.4% → 13.9%. Red moves by
0.8 points. And it is near-inert where the paper *was* validated: on SOBS it
moves the Fig. 7b sunlit-crown signature from 0.470 to 0.474.

**The recollision probability does not earn its place, and `run_flight_single.py
--black-soil` says why.** With a black soil the canopy alone satisfies
`BRF(all)/BRF(single) = 1/(1-omega*p)`, so p can be measured rather than
assumed. Against gortt's Wigner estimate `p = sigma*l/(1+sigma*l)` on the
crown's mean chord: 0.435 against a measured 0.37–0.46 for the thin crowns, but
0.909 against 0.53–0.57 for the thick ones. The first-flight escape probability
assumes scattering is spread uniformly through the crown; in a thick crown it
happens near the lit surface where escape is easier. Right in the thin limit,
wrong where it matters.

**`-rt_crown_crown` is not supported** (NIR RMSD 18.6% → 37.4%) and the options
do not compose — `-rt_all` overshoots by +43%, because `crown-path` and
`recollision` are two answers to the same question and `crown-crown` counts the
same photons a third time.

**The lidar is untouched by all of it.** `gortt_lidar` takes its backscatter
coefficients straight from the leaf reflectance and the soil albedo and never
calls `gortt_rt_init`, so no `-rt_*` option can reach it; `gortt -lidar -rt_all`
is byte-identical to `gortt -lidar`. The variant pages show the baseline lidar
figure and say so, and `run_lidar.py` defaults to the baseline alone.

### What this does not say

FLIGHT is a reference, not truth. It has its own approximations — leaves as a
turbid medium with a disc-based hotspot correction, an isotropic diffuse sky
(`g_sky_model` is forced to 0 on both branches of `flight.c:7409`, so the
anisotropic LUT is never read), and Monte Carlo noise of a few tenths of a
percent on the averaged BRFs here. The scenes are also idealised in both
models: no trunks, no understory, no senescent foliage, a Lambertian soil and
direct-beam illumination only.

## Imaging a single crown

`run_crown_image.py` is a different kind of comparison from the rest of this directory: it
does not compare BRFs at all, it measures what a *crown* looks like so that the crown
submodel has a target. One crown at the centre of a black scene, FLIGHT in reverse mode
with `OUTPUT_MODE 2` (the floating-point image), and the per-pixel radiance factor binned
against the impact parameter `u` of the view ray after the ellipsoid-to-sphere transform —
the same `u` as the chord average in `gortt_crown_reflectance`. `plot_crown_image.py` draws
it and `plot_crown_image.py --check` verifies the machinery.

Two traps, both guarded in the script:

* **FLIGHT is periodic**, so a ray that leaves a wall comes back on the opposite one at the
  same height and can meet the crown a second time, whenever `hw < r + 2b tan(vza)`. In a
  fixed 12 m box that put a clean second image of the crown out at `u = 2.5` for the prolate
  crowns at 60 degrees. The box is therefore sized per view zenith, which also keeps the
  cost down, since it goes as `hw^2`.
* **The mapping from pixel to `u` is analytic**, so it can be asserted on: every run checks
  that no lit pixel falls outside the predicted silhouette, with a 0.01 tolerance for rays
  that graze the limb. That check is what caught the wrapping.

What it found, and what came of it, is in `CROWN-BRDF-INVESTIGATION.md`.

## Scope

Deliberately **not** compared: the areal proportions of the scene components
(gortt's Kc/Kg/Kt/Kz against FLIGHT's per-component fractions).

**Lidar waveforms are compared as normalised shapes**, with the energy
partition shown separately, because the shape is `dP/dz` alone and is
independent of both backscatter coefficients. The partition is now a real
comparison and it agrees — see the lidar findings above.

It did not start that way. gortt originally built `rho_v` from the leaf
**transmittance** and halved the soil reflectance, following Ni-Meister et al.
(2001), and the resulting factor-of-two discrepancy in the energy split was
recorded here as a finding rather than worked around. It turned out to be an
error in that paper, corrected by Ni-Meister et al. (2018) and now fixed in
gortt. The change was made in the model, on the paper's authority, not tuned
against FLIGHT.

## Cost

FLIGHT does one view direction per process — `MULT_VIEWS` is compiled out and
marked "Not yet properly implemented" — so the BRDF sweep is
12 scenes × 3 realisations × 3 solar zeniths × 64 directions ≈ 6,900 FLIGHT
processes at about 4.5 s each. `run_brdf.py` checkpoints per
(scene, realisation) into `work/brdf/checkpoints/`, so an interrupted sweep
resumes. `--dry-run` prints the count and estimate; `--quick` is a smoke test.

`DATA/hotspot2.data` (239,700 lines) is re-parsed on every run, which is about
half a second of the 4.5. Avoiding it would mean editing FLIGHT.

## Files

| | |
|---|---|
| `config.py` | scene matrix, optical properties, angle grids — change things here |
| `build_flight.py` | the patched, lidar-enabled FLIGHT build |
| `scene.py` | FLIGHT input writers, crown generation, both models' runners and parsers |
| `validate.py` | the consistency checks above |
| `run_brdf.py`, `run_lidar.py` | the sweeps; write `results/*.csv` |
| `plot_brdf.py`, `plot_lidar.py` | the figures; write `figs/*.png` and the stats CSVs |
| `make_html.py` | `index.html` |
| `run_flight_single.py` | single-scattering and black-soil FLIGHT diagnostics for the `-rt_*` options |
| `run_crown_image.py`, `plot_crown_image.py` | the single-crown imaging measurement — see below |
| `variants/` | one comparison page per `-rt_*` option, plus an index |
| `results/` | tidy CSVs — **not tracked**, see below |
| `figs/` | the figures |
| `work/` | per-run scratch and the sweep checkpoints — not tracked |

### What is in git

The scripts, this file, `figs/` and `index.html`. **Nothing either model
produced is tracked**: `results/` and `work/` are gitignored, along with
`.flight_binary` (the path to the patched build, which lives outside both repos
and differs per machine).

Nothing is lost by that — every number this comparison reports is in the
figures and baked into `index.html`, which opens over `file://` with no server
and no CSVs. What it means is that the plotting scripts cannot redraw until the
models have been re-run:

```sh
./build_flight.py && ./validate.py && ./run_brdf.py && ./run_lidar.py
./plot_brdf.py && ./plot_lidar.py && ./make_html.py
```

`run_brdf.py` is the expensive part at about three hours, and it checkpoints per
(scene, realisation) so it can be interrupted and resumed. `run_lidar.py` is ten
minutes and `run_brdf.py --gortt-only`, which is all the contour panels need, is
seconds.
