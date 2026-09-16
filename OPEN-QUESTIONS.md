# Open questions

Things found during the 2026 tidy-up that were **deliberately not changed**, because fixing
them is a modelling decision rather than a refactor, or because they need checking against the
papers first. Nothing here is a regression: all of it predates the tidy-up.

The tidy-up itself is on the `tidy-up` branch; `v0-pre-refactor` tags the code as it was.
Every commit outside the "behaviour change" ones is byte-identical on all 51 test cases, so
`make test` is the instrument for measuring anything below: change one thing, run it, and the
harness prints the largest absolute and relative deviation per case.

**Needs resolving first: item 9, β.** Since items 1 and 7 were fixed, β is the largest
remaining error in the geometric-optical proportions.

---

## 1. `gortt_kc` carries F squared — affects every BRF, principal plane included

> **Status: resolved** on branch `kc-investigation`.
> - F² removed in 129ed80: L&S (1992) eq. 25 interpolates f, and eq. 17 already contains F.
> - That exposed an unbounded forward-scatter Kc, which came from `Po` in `gortt_kc_fFbeta`.
>   Fixed in 0204222 with Po = min(PvMv, PiMi). That rule in turn made Kc negative with the
>   sun near overhead, so it was replaced by Schaaf, Li & Strahler (1994) eq. 23A in ed690df:
>   Po = PiMi beyond the hotspot in backscatter, min(PvMv, PiMi) elsewhere.
> - The overlap was also using the canopy depth where L&S use the crown-centre height. Fixed in
>   59121b9.
>
> All three were decided against a Monte Carlo reference; see **KC-INVESTIGATION.md**. The `Po`
> change has since been confirmed in the source: Li & Strahler (1992, IGARSS), conclusion 1 and
> Fig. 1, model Po = PiMi beyond the hotspot and Po = PvMv inside it "until the large zenith in
> the forward direction is reached where PvMv > PiMi again", i.e. the smaller of the two. gortt
> had dropped that last clause. β for an overhead
> sun was also corrected to its limit of 1 (b59ce42). What is left is β (item 9) and a small
> residual Kt clamp (item 11). The text below is the original
> question, kept for the record.

**Where:** `gortt_brdf.c`, `gortt_kc` and `gortt_kc_fFbeta`.

`gortt_kc_fFbeta` returns `f` with `F` already folded in:

```c
*f = *F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;
```

`gortt_kc` then interpolates between the two principal-plane values:

```c
f = (1.-frac)*f0*F0 + frac*f180*F180 ;
```

so each product is `F·(…)·F` — **F squared**. `gortt_kc_PP_only`, which is the same equation
written out directly, uses `F` to the first power:

```c
f1 = F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;
f  = beta*f1 + (1.0-beta)*F ;
```

Both cannot be right. Either `gortt_kc_fFbeta` should not multiply by `*F`, or `gortt_kc`
should not multiply again.

**Also:** the first `gortt_kc_fFbeta` call — the one at the actual geometry — has its `f`
discarded. Only its `F` and `beta` are used. That may be intentional (the whole point of the
routine is to interpolate the principal-plane `f`s) but it is worth confirming at the same
time, because it means the `f` expression is evaluated three times and used twice.

**To check:** Li & Strahler IGARSS '92 and RSE '94 for the `f`/`F`/`beta` decomposition, and
Ni et al. (1999) for how `Kc` is assembled from them.

**Blast radius:** every geometry, *including* the principal plane — which is wider than the
heading above suggests. At `raa = 0` the interpolation collapses to `frac = 0` and
`f = f0*F0`, and at `raa = 180` to `f = f180*F180`; both still carry the extra `F`. The
interpolation weight is what varies off the principal plane, not whether the extra factor is
there. `default`, `beta`, `hbrpcc` and `oldgeom` in the test suite all contain such angles, so
the delta will be immediately visible. There is a comment in `gortt_kc` recording all of this.

**Evidence (added when reproducing Fig. 5 of Ni et al., see `figures/`):** that figure is the
areal proportions in the principal plane for the SOBS site, and SOBS is our own default
canopy, so it is close to a direct test. `Kt`, `Kz` and `Kg` reproduce the published curves,
and `Kc` matches where the published figure can be read off precisely — 0.817 against ~0.81 at
the hotspot, 0.167 against ~0.15 at nadir. But `Kc` is too low everywhere else in both
directions, and dropping the extra `F` in a scratch build moves it toward the published curve
at every point off the hotspot:

| vza (szn=30) | as shipped | extra `F` dropped | Fig. 5 |
|-------------:|-----------:|------------------:|-------:|
| +30 (hotspot) | 0.817 | 0.817 | ~0.81 |
| +44 | 0.588 | 0.732 | ~0.79 |
| +60 | 0.712 | 0.856 | ~0.80 |
| -40 | 0.059 | 0.143 | ~0.13 |
| -60 | 0.043 | 0.191 | ~0.15 |

Closer at every one of them, though it does not fully close the gap at +44, so this does not
by itself prove the single-`F` form is the right one — it is a reason to do the paper check
above rather than a substitute for it. Note also that `gortt_kc_PP_only` is never called from
anywhere, so the principal plane does not currently take the single-`F` branch that exists for
it.

---

## 2. `-nan` at exactly 90 degrees zenith

**Where:** the `SEC()` calls throughout `gortt_brdf.c` and `gortt_rsurf`.

At vza or sza of exactly 90°, `cos` is zero and `SEC` is infinite, so the BRF is `-nan`. Both
the direct path and `-P` do this, identically — it is not a round-trip problem.

`gortt_init_params` already clamps the *probability table's* last angle:

```c
p->theta[i] = p->dth * (double) i;
if( p->theta[i] >= M_PI/2.0 ) p->theta[i] = M_PI/2.0 - 1.0*M_PI/180.0;
```

but nothing clamps the geometry the user asks for.

**Options:** clamp the input geometry the same way, return zero (a horizon view sees no
canopy), or handle the limit analytically. The first is what the tables already do and is
therefore at least self-consistent; the third is the honest one.

**Pinned by:** the last row of `extreme` and `useprobs-extreme`.

---

## 3. The top of the theta grid is a duplicated abscissa

**Where:** `gortt_init_params` (the clamp above) and `gortt_calc_kopen`.

Because everything at or above 90° is clamped to 89°, `theta[89]` and `theta[90]` hold the
*same* angle. `gortt_calc_kopen` then integrates with the trapezium rule at the full `dth`
width, so the final panel spans zero actual angle but is given a full panel's weight — in
effect the 89°–90° interval is filled with the 89° value.

Since Phase 3.2 filled `epgap` at `t = nth-1`, that panel is now `E*dth` rather than
`E*dth/2`. Either way it is worth about 1e-51 against a `k_openep` of 0.100, so there is no
numerical consequence **at the current one-degree resolution**. It would start to matter if
`dth` were coarsened, or if item 2 above were fixed by extending rather than clamping.

---

## 4. `p->k` is evaluated at the wrong solar zenith — a trap for the lidar work

**Where:** `gortt_init_params`.

```c
p->k = gortt_leaf_angle_distribution( p, g->sza ) ;
```

This runs before the per-angle loop, and `gortt_cl_parser` does not set the geometry, so
`g->sza` is still the default zero. `p->k` is therefore the leaf projection at nadir sun, used
for every geometry.

Harmless **only** because the leaf angle distribution is hardwired to return 0.5 regardless of
angle (item 5). The moment a real LAD goes in, `p->k` has to move out of `gortt_parameters` and
into the per-geometry path alongside `g->k_vza`, which is already there and already computed
per angle. There is a comment at the assignment saying so.

---

## 5. The leaf angle distribution is hardwired

`gortt_set_defaults` sets `p->lad = LAD_05` and no command-line option changes it.
`gortt_leaf_angle_distribution` warns and returns 0.5 for anything else. The `LAD_PLANOPHILE`
… `LAD_UNIFORM` constants and their branches exist and are unreachable.

Fixing this needs item 4 fixed first.

---

## 6. `nh_probs` bounds the height dimension

Not a bug — a constraint to remember. Only rows `h < p->nh_probs` of `epgap` are filled; the
rest are left at their `calloc` zero. It defaults to 1 because the BRF only needs the ground,
and computing all fifteen layers was about 95% of the runtime.

`-lidar` already raises it to `nlayers`. **Anything else that indexes those tables above
`h = 0` must raise it too**, or it will silently read zeros.

Three refinements since `-nlayers` and the lidar waveform went in:

* `p_n0` is the exception — it is cheap and is filled for every layer at any `nh_probs`. It is
  `epgap`, `k_open` and `k_openep` that stop at `nh_probs`.
* the raise now happens **after** the argv loop, not inside the `-lidar` branch. With
  `-nlayers` in the parser the two options are order-dependent, and `gortt -lidar -nlayers 30`
  would otherwise leave `nh_probs` at the default 15. `gortt_lidar` checks
  `nh_probs >= nlayers` on entry and refuses rather than reading zeros.
* `pd_s` is gone. It was a histogram of within-crown path lengths, written in one place and
  read in one, and the round trip through it cost a quantisation of every path length before
  it was exponentiated. `gortt_calc_epgap_ht` now accumulates `exp(-s tau')` directly.

---

## 7. Alternative published formulations are behind `if(TRUE)` / `if(FALSE)`

> **Status: resolved by deletion** in 5399f53. Every one of these was scored against a Monte
> Carlo reference in KC-INVESTIGATION.md (sections 3 and 4). None beats the branch that was
> taken, and `gortt_kc_PP_only` and Wenge's overlap were broken outright. Their provenance is
> kept as switches in `figures/goms_kernels.py`, where A/B-ing them is a keyword argument. The
> text below is the original item.

**Where:** `gortt_overlap` (Ambrals vs Li & Strahler '92; and Wenge's off-principal-plane
variant vs Li & Strahler '92) and `gortt_kc_PP_only` (IGARSS '92 vs RSE '94, and TGARSS '92).

These were kept, deliberately: they encode real provenance and deleting them would destroy
information. But `if(TRUE)` is a trap — it reads as dead code and the compiler does not stop
anyone from "cleaning it up".

**Suggested:** promote each to an enum in group 1 of `gortt_parameters`, defaulted to the
branch currently taken and selectable from the command line. That is bit-identical by
construction and turns A/B-ing the published formulations into a command-line operation, which
is the thing you would actually want when writing this up.

---

## 8. Smaller usability items

* **`-Anth` is documented but not accepted.** `gortt_usage` advertises it; `gortt_cl_parser`
  matches `-canth` at three characters (i.e. `-ca`), so the documented spelling is rejected as
  an unknown option. Pinned by the `badanth` case. Fix either the usage text or the parser —
  and note the parser is prefix-matched and order-dependent, so adding a new spelling needs
  care about both the length argument and its position in the if-chain.
* **`-diffuse` is undocumented in `-u`**, and TODO.md flags the default as "a bit dangerous".
  It sets `p->user_fd = 1 - frac`; without it `g->fd` is derived per angle from the solar
  zenith following Ni et al. (1999).
* **The stdin header line is capped at `MAX_LINE_LEN` (1000 characters)**, so roughly 200
  wavebands. This is not silent corruption — a truncation mid-number yields a wavelength below
  400 and trips the range check, and a truncation between numbers trips the
  declared-versus-found count check — but the resulting message is misleading. Reading the
  header with a growable buffer, or saying so explicitly, would both be better.
* **`-soil_spectra` files have no comment syntax.** Two numeric columns only; anything else is
  a fatal error. Worth documenting, or worth allowing `#`.

---

## 9. β, the weight between uniform-height and random mutual shadowing — **needs resolving**

**Where:** `gortt_kc_fFbeta` (the `*beta = …` expression) and `gortt_kc`
(`f = beta*f + (1-beta)*F`).

Kc = β·Kc(uniform height) + (1 − β)·Kc(random). The uniform-height extreme (β = 1) is now right:
for same-height crowns Kc matches the Monte Carlo with RMSE 0.074. How gortt mixes the two
extremes for real height distributions is not:

* **The MC's effective β depends strongly on view direction.** It is 0.8–1.0 in backscatter,
  where you see sunlit crown tops whatever their heights. It is 0.2–0.6 in forward scatter,
  where the height spread decides whether lower crowns' shadows land on the crowns you see.
  Near nadir it is 0.5–0.8. Table in KC-INVESTIGATION.md, section 6.
* **gortt's β** is IGARSS '92 eq. 7, with D = r cot(i′/2) and, since c0ef2c3, Γi = πr² sec(i′),
  the illumination shadow area (it used to be Γv). It depends only on the canopy and sza, so it
  cannot tell backscatter from forward scatter. It is too low for backscatter, which leaves a
  notch in Kc just beyond the hotspot where the (1 − β)F(1 − Kg) random term falls sharply
  (`figures/kc_investigation/kc_variants.png`). Current Kc RMSE against the MC: 0.064 overall;
  backscatter 0.095, near nadir 0.044, forward 0.060.
* **L&S (1992) eq. 23**, β = (1 − (h2 − h1)/4b)², is constant in angle: 0.24 for SOBS, and 0 once
  h2 − h1 > 4b. It is better in forward scatter and worse in backscatter; overall Kc RMSE 0.101
  against gortt's 0.064. β = 1 everywhere is worse again (0.141, and it restores 88 Kt clamps
  off the principal plane).
* **L&S apply β to the P terms** (eq. 24, P = βP1 + (1 − β)P2). gortt applies it to f. These are
  not the same thing, because F, M and the Γs differ between the extremes.
* **The "RSE '94" variant**, deleted with `gortt_kc_PP_only` and marked "totally wrong as far as I
  can make out", is **Abuelgasim & Strahler (1994), RSE 47:261–275, eq. 8**
  (`papers/Abuelgasim_and_Strahler_1994_RSE.pdf`, p. 265), transcribed exactly. As printed, eq. 8
  is IGARSS eq. 7 with both exponents' signs flipped: `1 − e^(−λΓi + (h2−h1)/D)` over
  `1 − e^(λΓi)`. It is a misprint. The paper's own text says β "will decrease from 1 to 0 with
  increasing canopy depth" and "increases … to 1 with increasing coverage". IGARSS eq. 7 does
  that, while eq. 8 as printed gives β = −0.61 at zero depth (λΓi = 0.5) and 12.5 at
  (h2 − h1)/D = 5. The same passage calls λΓi the "coverage and sun position", which is a second
  source for Γi being the illumination-direction area rather than gortt's Γv. Nothing further
  is needed on this variant.

**What Li & Strahler (1992, IGARSS) settles** (`papers/Li_and_Strahler_1992_IGARSS.pdf`,
pp. 766–768):
* **gortt's β is its eq. 7, except for Γi** (reprinted, with sign misprints, as Abuelgasim &
  Strahler 1994 eq. 8; see above). The paper defines Γi(h) as "the average shadowing
  areas at level h, of all the crowns which may cast shadows on that level". It uses the
  subscript i for illumination throughout (θi, Mi, Pi, and Γ(h,z) = Γi + Γv − O). gortt sets
  `Gamma_i = Gamma_v`, and that alone is what makes its β depend on vza. It reads like a
  transcription slip, though it could have been deliberate. Abuelgasim & Strahler (1994), who
  describe λΓi as the "coverage and sun position", point the same way.
* **Adopted in c0ef2c3.** Against the MC, Γi = πr² sec(i′) took overall Kc RMSE from 0.080 to
  0.064 and the worst error from 0.58 to 0.33. It halved the forward-scatter error (0.130 → 0.060
  in the principal plane) and made backscatter worse (0.065 → 0.095). SOBS 0.084 → 0.065, dense
  0.107 → 0.071, same-height unchanged. Better on balance, but not the fix. The old form is
  `beta="code"` in `figures/goms_kernels.py`.
* **The published β cannot depend on view direction, by construction.** Eq. 6, the "probability
  weighted overlap function" it comes from, is defined at nadir viewing. Option 3 below is
  therefore done, and the derivation cannot supply the forward/backscatter asymmetry the MC
  shows. The paper's own conclusion 4, "the Random Case hardly exists for spheroid shaped primary
  grains … realistic structures can be approximated as UC with good accuracy", matches the MC in
  backscatter but not in forward scatter.

**What Schaaf, Li & Strahler (1994) adds** (`papers/Schaaf_et_al_1994_TGRS.pdf`, appendix): eq. 24A
states Γi = πR² sec θi′, "the illumination shadow projected on the ground", which confirms
c0ef2c3. Eq. 23A applies β to f exactly as gortt does, so option 2 below is not how the model's
authors implemented it. Option 2 is still a legitimate experiment, but it is not a correction.
With the sun overhead, β now takes its limit of 1 (b59ce42). Against the MC at sza 0 that helps
SOBS and the dense canopy and overshoots the sparse one at large view zeniths, which is this
item again (KC-INVESTIGATION.md, section 9).

**Options, roughly in order of how much they change:**
1. A β that is larger in backscatter than in forward scatter, e.g. a function of cos(raa) or of
   the signed zenith difference. It could be fitted to the MC, or derived as the probability that
   a viewed crown's shadowing neighbours lie in the same height layer.
2. Apply β to ΔP as eq. 24 does, rather than to f, and re-score it.
3. ~~Find the IGARSS '92 derivation and check the formula in use against it.~~ Done, see above.
   The one discrepancy, Γi, is fixed in c0ef2c3.
4. The paper's eqs 1–2 give PiMi for a "fuzzy" Mi boundary as a numerical integral over the
   crown surface; eq. 5, which gortt uses, is its closed-form approximation. It is not a β fix,
   but it is the more physical form of the uniform-height extreme that β weights.

**To test any of them:** add a `beta` option to `figures/goms_kernels.py` and rerun
`figures/kc_investigation.py`. The MC is cached, so rescoring every variant takes seconds.
Compare with the "Effective beta" and "Kc" tables in `tables.md`. The change will move most BRF
goldens.

---

## 10. `-h1` equal to `-h2` hangs; nearly equal gives `-nan`

`gortt -h1 5.75 -h2 5.75` does not return: the gap-probability calculation never finishes.
`-h1 5.75 -h2 5.76` returns `-nan`. Both predate the kc-investigation changes; master
(e109ff0) does the same. It has not been traced. `p->lv = lambda/(h2 - h1)` in
`gortt_init_params` is the obvious suspect for the first.

Same-height crowns are a legitimate canopy (L&S Fig. 3a, and the cleanest test of the
geometric-optical kernels), so gortt should either handle them or refuse `h2 <= h1` with a
message. The overlap average added in 59121b9 already handles h1 == h2.

---

## 11. A small residual Kt clamp off the principal plane, near the hotspot

`gortt_kc` interpolates f in azimuth between raa = 0 and π (L&S eq. 25) and multiplies by
(1 − Kg) at the actual azimuth. Just off the hotspot, f0 ≈ 1/(1 − Kg0) is paired with a Kg that
has already fallen and a Kz that has already grown. Kc can then exceed 1 − Kg − Kz, and
`gortt_rsurf` clamps Kt.

On the investigation's grid this happens in 19 of 1,235 geometries, all on the principal cone
vza = sza, 10–90° of azimuth from the hotspot, mostly at sza 15°, and only for same-height and
sparse canopies. The excess is at most 0.10.

The alternatives tried (interpolating ΔP, or the P terms at the actual azimuth) remove the clamp
but double to triple the off-plane error (KC-INVESTIGATION.md, section 5). It may go away with
item 9, since β = 1 makes it much worse.

Li & Strahler (1992, IGARSS), conclusion 3, says a cosine in azimuth fits the sunlit crown area
Ac better than the linear interpolation their model used. Against the MC it does not help here:
* interpolating f with weight (1 − cos raa)/2 (`az_interp="fcos"`) gives off-plane RMSE 0.061
  against 0.057 for linear, with 28 clamps instead of 19;
* interpolating Kc itself removes the clamps but is worse again, 0.079 linear and 0.092 cosine.

---

## What is *not* on this list

The four bugs fixed in Phase 3 (`-P` round trip, the `epgap` grid, the Q08 seeding,
`-soil_spectra` and the black soil), and the six memory-safety bugs fixed in Phase 1. Those are
done and covered by the suite. See `git log` on the `tidy-up` branch — each commit records the
measured delta.
