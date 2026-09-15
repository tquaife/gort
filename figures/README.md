# figures

Scripts that reproduce published figures with `gortt`. Each takes `--help` and
writes its image into this directory. `gortt_run.py` is the shared helper that
invokes the model and parses its output; the rest are the figures.

The papers themselves live in `papers/`, which is gitignored; nothing here
reads them, so the scripts work without it.

The site is the same throughout: the BOREAS Old Black Spruce stand (SOBS),
Table I of Ni et al. (1999) — R=0.76, b=2.7, lambda=0.41, h1=3.0, h2=8.5,
FAVD=0.858. That is near enough `gortt`'s own defaults; `gortt_set_defaults`
*is* the SOBS site. The scripts pass them explicitly so the figures do not move
if a default ever does. (The paper's Table I heading gives FAVD as 1/m², which
is where `gortt`'s usage message used to get it from. It is 1/m.)

## fig5_ni1999.py

Figure 5 of Ni et al. (1999): the modelled areal proportions of the four scene
components (Kc sunlit crown, Kg sunlit background, Kt shaded crown, Kz shaded
background) in the principal plane at the BOREAS Old Black Spruce (SOBS) site,
for solar zenith angles of 30 and 60 degrees.

```sh
make                      # if ../gortt is not built yet
figures/fig5_ni1999.py                     # -> figures/fig5_ni1999.png
figures/fig5_ni1999.py --style colour      # coloured, same line styles
figures/fig5_ni1999.py --csv fig5.csv      # also dump the numbers
```

It takes about a second: one `gortt -prnprop` run per panel, since the four
proportions are exactly what `-prnprop` prints.

Positive viewing zenith is the backscatter direction (relative azimuth 0), so
the hotspot falls at +30 and +60 as in the paper; negative viewing zenith is
forward scatter (relative azimuth 180). The default range stops at +/-89
degrees because both the direct and `-P` paths return `-nan` at exactly 90.

### How good is the reproduction

> **Update.** The discussion below was written before `KC-INVESTIGATION.md`. Since then F² has
> been removed from `gortt_kc`, `Po` changed and the overlap height fixed. Kc at sza 30 now reads
> 0.817 at the hotspot, 0.811 at +44, 0.868 at +60, 0.242 at nadir, 0.143 at −40 and 0.191 at
> −60. A Monte Carlo of the same canopy gives 0.320 at nadir, twice the published ~0.15, so the
> published Kc curve is not a reliable target (KC-INVESTIGATION.md, section 7). What remains
> between gortt and the MC is mostly β (OPEN-QUESTIONS #9).

Good, with one exception.

* Kt, Kz and Kg reproduce the published shapes: Kt is a sine wave that is zero
  at the hotspot, Kz peaks at nadir and is zero at the hotspot, and Kg is
  close to zero throughout the szn=60 panel -- which is exactly what the
  published caption says of it.
* Kc matches at the two points that are easiest to read off the published
  figure: 0.817 against ~0.81 at the hotspot (szn=30), and 0.167 against ~0.15
  at nadir. The four proportions sum to one to machine precision.
* **Kc is too low away from the hotspot, in both directions.** The published
  curve is roughly flat at ~0.80 from the hotspot out to +60 and sits at ~0.13
  through the middle of the forward-scatter side. We give 0.59 at +44 and
  0.04-0.09 across the forward side.

That last point looks like the F-versus-F-squared question recorded in
`OPEN-QUESTIONS.md`. `gortt_kc` builds its principal-plane value as

```c
f = (1.-frac)*f0*F0 + frac*f180*F180 ;
```

and `gortt_kc_fFbeta` already returns `f` with a factor of `F` folded in, so
this carries `F` squared. `gortt_kc_PP_only`, which uses `F` to the first
power for the same equation, is defined but never called -- and the principal
plane is precisely where it would apply, which is what makes this figure a
test of it.

Dropping the extra factor in a scratch build (the model itself is untouched)
moves Kc toward the published curve on both sides:

| vza | Kc, as shipped | Kc, extra F dropped | Fig. 5 |
|----:|---------------:|--------------------:|-------:|
| +30 | 0.817 | 0.817 | ~0.81 |
| +44 | 0.588 | 0.732 | ~0.79 |
| +60 | 0.712 | 0.856 | ~0.80 |
| -40 | 0.059 | 0.143 | ~0.13 |
| -60 | 0.043 | 0.191 | ~0.15 |

Closer on every point off the hotspot, though it does not fully close the gap
at +44. This is evidence, not a diagnosis: `OPEN-QUESTIONS.md` asks for the
expression to be checked against Ni et al. and Li & Strahler before anything
is changed, and any change here moves the regression goldens.

### Reading the published figure

For comparison, the caption of the original describes: Kg dashed, peak at the
hotspot and a local peak at nadir; Kc solid, a bowl shape peaking at the
hotspot; Kt dash-dot, a sine wave with a local minimum at nadir and zero at
the hotspot; Kz dotted, peak at nadir and zero at the hotspot. This script
uses the same four line styles, and adds a legend.

## fig67_ni1999.py

Figures 6 and 7 of Ni et al. (1999): the component signatures C, G, T and Z for
the SOBS site, in the visible (Fig. 6) and the near infrared (Fig. 7). Each
figure has two panels — (a) the signatures against solar zenith angle, with
vertical bars showing the range of C over viewing direction, and (b) C at a
solar zenith of 45 degrees against viewing zenith in the principal plane.

```sh
figures/fig67_ni1999.py                    # -> fig6_ni1999.png, fig7_ni1999.png
figures/fig67_ni1999.py --band red         # just Fig. 6
figures/fig67_ni1999.py --soil 0.29         # restore Table II's red background
```

About five seconds; the vertical bars need a full principal-plane sweep at each
solar zenith they are drawn at.

### Bypassing PROSPECT and the Price soil

The figures use the fixed optical properties of the paper's Table II, so both
sub-models are switched off from the command line — `-refl_leaf`/`-trans_leaf`
for the leaf (which turn PROSPECT off) and `-alb_soil` for the background
(which turns the Price EOFs off):

| band | `-refl_leaf` | `-trans_leaf` | `-alb_soil` |
|------|-------------:|--------------:|------------:|
| 540 nm | 0.2   | 0.06  | 0.069 |
| 670 nm | 0.09  | 0.002 | 0.117 (see below; Table II says 0.29) |
| 865 nm | 0.515 | 0.305 | 0.29 |

Reflectance and transmittance go in separately because the radiative transfer
needs both: the single-scattering albedo is their sum, but `gortt_brdf.c:425`
builds the scattering asymmetry from their difference,

```c
rt->g_asym = -(4.0/9.0)*( *(s->rleaf+iw) - *(s->tleaf+iw) )/( rt->w ) ;
```

Table II's leaves are strongly asymmetric — `g_asym = -0.43` at 670 nm and
`-0.11` at 865 nm — so this matters. The older `-alb_leaf` set `r = t = arg/2`
and therefore always gave `g_asym = 0`; it was removed when the pair went in,
because a single albedo cannot describe a real leaf.

Using the true r and t raises C at the hotspot appreciably and barely moves it
at nadir:

| | nadir | hotspot |
|---|---:|---:|
| 670 nm | 0.0395 → 0.0409 | 0.098 → 0.124 |
| 865 nm | 0.4791 → 0.4797 | 0.584 → 0.596 |

### How good is the reproduction

The shapes come out right and the magnitudes are close but not exact.

* Panel (a) reproduces well in the red: C starts at 0.037 against a published
  ~0.035 and peaks at 0.149 against a published ~0.14. Our peak sits at ~55
  degrees where the paper's is nearer 68.
* Panel (b) has the right shape in both bands — a flat forward-scatter side
  with a sharp hotspot spike at +45 — but overestimates the spike. We give
  0.041 flat / 0.124 peak in the red against a published ~0.032 / ~0.09, and
  0.48 / 0.596 in the NIR against ~0.365 / ~0.45.
* The two panels are mutually consistent, and so is the paper: interpolating
  the published panel (a) to a solar zenith of 45 gives ~0.09, which is what
  its panel (b) shows at the hotspot. Ours reads ~0.12 both ways.
* Giving Table II's true asymmetric leaf (rather than the symmetric one the
  old `-alb_leaf` forced) pulled panel (a) toward the published curve — the
  red peak went 0.117 → 0.149 against ~0.14 — and pushed panel (b)'s hotspot
  away, 0.098 → 0.124 against ~0.09. It is the right input either way; the
  residual disagreement is the model's, not the option's.
* **The published solid curve in panel (a) is C at the hotspot, not at a fixed
  view.** The paper says C "has a maximum value at the hotspot viewing
  direction", and C at a fixed nadir view *falls* slightly with solar zenith
  instead of rising. G, T and Z are drawn at nadir, where the paper notes they
  are almost constant with view direction.
* **The red background albedo is 0.117, not Table II's 0.29** — see below.

### The red background albedo

G is not modelled: `gortt.c:697` sets `Gd = Gf = rsoil`, so the direct/diffuse
mixture `G = fd*Gd + (1-fd)*Gf` collapses to `rsoil` exactly, for any `fd`.
Verified: G tracks `-alb_soil` to floating-point representation across 0.069,
0.12, 0.29 and 0.5. So the flat dashed line in the published Fig. 6(a) *is*
the background albedo that was used to draw it, read directly off the page.

Measured against the panel's own axis calibration (top frame 0.4, bottom axis
0.0), that line sits at **0.117** — not the 0.29 that Table II gives for
670 nm. Two things support 0.117 being the value actually used:

* Table II lists 0.29 for both 670 nm and 865 nm. An identical red and NIR
  background is not what a boreal forest floor does, and it has the look of
  the NIR row duplicated.
* 0.117 brings panel (b) closer to the published curve at both ends than 0.29
  does, and panel (a) closer still: with it, C starts at 0.037 against a
  published ~0.035.

So the red band here defaults to `soil=0.117`. `--soil 0.29` restores the
tabulated value if you want to see it.

Fig. 6 is 670 nm rather than 540 nm: panel (b) at nadir is ~0.032 in the paper,
against 0.041 at 670 nm and ~0.10 at 540 nm from the model.

## Kernel investigation: goms_kernels.py, goms_mc.py, kc_investigation.py, kc_fsquared_sensitivity.py

These are the tools behind `KC-INVESTIGATION.md` at the top of the repo. They are not figure
reproductions:

* `goms_kernels.py` is a numpy port of gortt's geometric-optical kernels, with every contested
  choice (F², Po, the P-term forms, the overlap height and azimuthal form, β) behind a switch.
  `--selftest` checks that the defaults reproduce `gortt -prnprop` to about 1e-13.
* `goms_mc.py` is a Monte Carlo of the GOMS scene: solid spheroids at random positions and
  heights. It gives Kc/Kt/Kg/Kz and the L&S mutual-shadowing terms for one geometry.
* `kc_investigation.py` runs the MC over about 1,300 geometries (about 4 minutes on 4 cores,
  cached) and writes the tables and plots into `figures/kc_investigation/`, which is gitignored.
* `kc_fsquared_sensitivity.py --old <binary>` compares two gortt binaries across a sweep of
  canopies and sun angles. It was written for the F² change.

Note that the Fig. 5 discussion above predates that work. The Monte Carlo puts nadir Kc for SOBS
at about 0.32, not the ~0.15 of the published figure, so Fig. 5 is not a good target for Kc
(KC-INVESTIGATION.md, section 7).

