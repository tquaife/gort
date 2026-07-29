# figures

Scripts that reproduce published figures with `gortt`. Each one is standalone,
takes `--help`, and writes its image into this directory.

The papers themselves live in `papers/`, which is gitignored; nothing here
reads them, so the scripts work without it.

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

The tree parameters are Table I of the paper (R=0.76, b=2.7, lambda=0.41,
h1=3.0, h2=8.5, FAVD=0.858), which are near enough `gortt`'s own defaults --
`gortt_set_defaults` *is* the SOBS site. They are passed explicitly anyway so
the figure does not move if the defaults ever do.

Positive viewing zenith is the backscatter direction (relative azimuth 0), so
the hotspot falls at +30 and +60 as in the paper; negative viewing zenith is
forward scatter (relative azimuth 180). The default range stops at +/-89
degrees because both the direct and `-P` paths return `-nan` at exactly 90.

### How good is the reproduction

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
