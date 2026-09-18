# figures

Scripts that reproduce published figures with `gortt`. Each takes `--help` and
writes its image into this directory. `gortt_run.py` is the shared helper that
invokes the model and parses its output; the rest are the figures.

The papers themselves live in `papers/`, which is gitignored; nothing here
reads them, so the scripts work without it.

The images are gitignored too, with one exception: `fig5_ni1999.png` is tracked,
because the top-level `README.md` shows it on the GitHub page. Regenerate it
with `fig5_ni1999.py` after any change that moves Kc, Kg, Kt or Kz, and commit
the new image along with that change. It is byte-identical when regenerated
from the same code with the same matplotlib, so it only shows up in `git status`
when something really has changed.

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
> been removed from `gortt_kc`, `Po` changed, the overlap height fixed and Γi in β corrected
> (c0ef2c3). Kc at sza 30 now reads 0.817 at the hotspot, 0.798 at +44, 0.824 at +60, 0.252 at
> nadir, 0.136 at −40 and 0.148 at −60. That is within a few hundredths of the published curve
> everywhere except nadir (~0.15). A Monte Carlo of the same canopy gives 0.320 at nadir, so the
> published Kc curve is still not a reliable target (KC-INVESTIGATION.md, section 7). What
> remains between gortt and the MC is mostly β (OPEN-QUESTIONS #9).
>
> Kc and Kt go flat below about −83° at sza 30. There P_vM_v overtakes P_iM_i, Po switches to
> P_iM_i, and the uniform-height part of Kc equals its bound 1 − Kg − Kz (all visible crown is
> sunlit). What is left varies only with F under a β that is fixed for the solar zenith. This is
> expected behaviour; the MC gives ~0.30 at −85.

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
at nadir. The figures below were measured before the gap probability
quadratures were replaced, so they are right about the *effect of the option*
but no longer match what gortt prints:

| | nadir | hotspot |
|---|---:|---:|
| 670 nm | 0.0395 → 0.0409 | 0.098 → 0.124 |
| 865 nm | 0.4791 → 0.4797 | 0.584 → 0.596 |

For reference, gortt now gives 0.0380 nadir / 0.0986 hotspot at 670 nm and
0.4720 / 0.5573 at 865 nm, with the true r and t.

### How good is the reproduction

The shapes come out right and the magnitudes are close but not exact.

* Panel (a) in the red: C starts at 0.036 against a published ~0.035, and
  peaks at 0.109 against a published ~0.14. Our peak sits at ~55 degrees where
  the paper's is nearer 68, and that offset is unchanged by anything done so
  far.
* Panel (b) has the right shape in both bands — a flat forward-scatter side
  with a sharp hotspot spike at +45. We give 0.035 flat / 0.099 peak in the
  red against a published ~0.032 / ~0.09, and 0.470 / 0.557 in the NIR against
  ~0.365 / ~0.45. The red is close; the NIR is still high on both.
* The two panels are mutually consistent, and so is the paper: interpolating
  the published panel (a) to a solar zenith of 45 gives ~0.09, which is what
  its panel (b) shows at the hotspot. Ours reads ~0.099 both ways.

Replacing the gap probability quadratures moved all of these, because C is
built from `k_open`/`k_openep`. It pulled panel (b) toward the published curves
and pushed panel (a)'s red peak past them:

| | before | now | published |
|---|---:|---:|---:|
| (a) red, C at sza 0 | 0.037 | 0.036 | ~0.035 |
| (a) red, peak C | 0.149 | 0.109 | ~0.14 |
| (b) red, forward-scatter | 0.037 | 0.035 | ~0.032 |
| (b) red, hotspot | 0.124 | 0.099 | ~0.09 |
| (b) NIR, forward-scatter | 0.477 | 0.470 | ~0.365 |
| (b) NIR, hotspot | 0.596 | 0.557 | ~0.45 |

Four of the six got closer and the red peak in panel (a) got further away,
going from a slight overestimate to an underestimate. That is the expected
shape of this kind of change: the old numbers were unconverged, so any
agreement they showed was partly luck. See `quadrature_convergence.py` below
for why the new ones are the right ones.
* Giving Table II's true asymmetric leaf (rather than the symmetric one the
  old `-alb_leaf` forced) pulled panel (a) toward the published curve and
  pushed panel (b)'s hotspot away. It is the right input either way; the
  residual disagreement is the model's, not the option's. (The numbers that
  used to be quoted here, 0.117 → 0.149 and 0.098 → 0.124, were measured
  before the quadrature work and are not directly comparable with the table
  above.)
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

## fig4_nimeister2001.py

Figure 4, case B, of Ni-Meister et al. (2001): the three lidar profiles against
height — the total nadir gap probability `P(z)`, the waveform `dP/dz`, and the
apparent foliage profile `dlogP/dz`. These are the first four columns
`gortt -lidar` prints, so the script is one model run per curve.

```sh
figures/fig4_nimeister2001.py                  # -> figures/fig4_nimeister2001.png
figures/fig4_nimeister2001.py --h2 both        # both candidate geometries (see below)
figures/fig4_nimeister2001.py --nlayers 120    # finer vertical grid, slower
figures/fig4_nimeister2001.py --csv fig4.csv   # also dump the model numbers
```

It prints rms and max differences against the published curves when it finishes.

Unlike the other scripts here this one needs a reference dataset, because the
published curves are not something `gortt` can produce: they are checked in as
`fig4_nimeister2001_digitised.csv`, traced off a 400 dpi render of the scanned
page. (`papers/` is gitignored, so the paper itself is not available to the
script.) That file's header documents how, and how accurately.

### Table I says h2 = 12; the figure was drawn with h2 = 10

Table I gives case B as R=1.2, b=3.6, lambda=0.24, h1=6, **h2=12**, Fa=0.44.
It does not reproduce the figure, and the figure itself says why. All three
panels carry dotted reference lines labelled z1, h1, h2 and z2; traced off the
page they sit at

    z1 = 2.40    h1 = 6.00    h2 = 10.00    z2 = 13.60   (m)

agreeing across the three panels to better than 0.03 m. GORT has z1 = h1 − b
and z2 = h2 + b, so z1 = 2.40 with h1 = 6 forces b = 3.60 — Table I's value
exactly — and z2 = 13.60 with that b then forces **h2 = 10.00**. No (h1, b)
pair gives both z1 = 2.4 and z2 = 13.6 with h2 = 12, and h2 = 12 would put z2 at
15.6 m, off the top of the figure's own 0–15 m axis.

So the two disagree, and this script follows the figure. `--h2 both` overlays
them, which is the quickest way to see which one `gortt` reproduces:

| geometry | P rms / max | dP/dz rms / max | dlogP/dz rms / max |
|---|---|---|---|
| h2 = 10 (figure)  | 0.007 / 0.014 | 0.0027 / 0.0068 | 0.0034 / 0.0118 |
| h2 = 12 (Table I) | 0.053 / 0.081 | 0.0130 / 0.0269 | 0.0147 / 0.0293 |

At h2 = 10 the rms difference in `P` is 0.007, which is *below* the digitisation
error itself (the two independent tracings of `P` in the csv differ by 0.01–0.03).
In other words the agreement is better than this comparison can resolve, and the
figure above is the limit of the tracing rather than of the model. The sigmoid
`P`, the single-peaked `dP/dz` and its position in the upper canopy all come out
right.

These numbers were 0.016 / 0.039 for `P` and 0.0093 / 0.0336 for `dP/dz` before
the gap probability quadratures were replaced. The fit to the published figure
improved by a factor of two to three, which is worth noting because nothing in
that work was tuned against this figure — it is an independent check that fell
out the right way.

### The wobble in the derivative panels

`gortt`'s `dP/dz` and `dlogP/dz` used to be visibly ragged where the published
curves are smooth. It was never the differencing — a centred difference has zero
gain at the Nyquist frequency, so it cannot manufacture a sawtooth. The roughness
was in the gap probability tables and differencing only exposed it, and it came
from four quadratures that stepped with the layer thickness. They have been
replaced; see `quadrature_convergence.py` below and the "No quadrature in the
model steps with `nlayers`" section of `CLAUDE.md`.

Measured as the rms second difference of `dP/dz` over its peak, for this canopy:

| `-nlayers` | dz (m) | before | now |
|---|---|---|---|
| 20  | 0.59 | 16% | 4.3% |
| 40  | 0.29 | 9%  | 1.1% |
| 80  | 0.14 | 6%  | 0.3% |
| 120 | 0.09 | 4%  | 0.1% |

It now falls as roughly `dz^2` rather than `sqrt(dz)`, so raising `-nlayers`
actually buys something. For a waveform the physical answer is still
`-lidar_pulse`, which convolves the remainder away along with everything else
below the pulse width.

## quadrature_convergence.py

Evidence that the gap probabilities no longer depend on `-nlayers`, and that the
quadratures behind them are converged at the default `-nquad`.

```sh
./quadrature_convergence.py
./quadrature_convergence.py --ref /path/to/an/older/gortt
```

Four panels, for three canopies (gortt's default, a dense one and a sparse one):

* **(a)** `k_openep` against `-nlayers`. Flat is the whole point. `--ref` overlays a
  second binary, so a build from before the change shows the old O(1/n) drift
  next to the new flat line.
* **(b)** `k_openep` against `-nquad`, i.e. whether the default orders are converged.
* **(c)** `p_n0[0][t]` against its closed form. GORT puts the lowest foliage at
  `z1 = h1 - b`, and `b` becomes `r` under the sphere transform, so the ground
  layer sits exactly one crown radius below the lowest crown centre. Every crown
  is then entirely above it, the projection cross section is the constant
  `pi r^2 sec(theta')`, and `p_n0[0][t] = exp(-lv' pi r^2 sec(theta') (h2_p - h1_p))`
  with no quadrature at all. gortt matches this to 2.8e-14 across 0–90°, including
  at 89° where the probability itself is 1e-65.
* **(d)** the raggedness of the lidar `dP/dz` profile against `-nlayers`.

Run it after anything that touches `gortt_pn_kopen.c`. It is the counterpart of
`kc_investigation.py` for the gap probabilities: that one scores the
geometric-optical kernels against a Monte Carlo, this one scores the
probabilities against closed forms and against their own convergence.

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


## Gap probability geometry: capsule_volume_check.c

The odd one out — C rather than python, because what it checks are internal functions with no
command-line path to them. It links against `gortt_pn_kopen.o`:

```sh
gcc -Wall -g -I../include -o capsule_volume_check capsule_volume_check.c ../gortt_pn_kopen.o -lm
```

It checks `gortt_capsule_cross_section` against a pixel count of the region it describes, and
`gortt_vol` against a Monte Carlo of the volume *it* computes. The first passes; the second
does not, by up to 6% where the crown centre plane cuts the region, which is
`OPEN-QUESTIONS.md` #14. Run it before touching `gortt_vol`, `gortt_trisec` or `gortt_cylind`,
and run it after, because it is the only independent check those functions have.
