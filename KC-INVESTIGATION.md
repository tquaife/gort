# The geometric-optical kernels: P terms, Po, the overlap, and beta

This investigation followed from OPEN-QUESTIONS #1. Dropping the F² from `gortt_kc` (Li &
Strahler 1992 eq. 25, commit 129ed80) matches the paper, but it exposed two problems:
- Kc > 1 − Kg − Kz at large forward-scatter view zeniths, where `gortt_rsurf` clamps Kt to zero.
  Kc reaches 3.8 at sza 15° for the SOBS stand.
- Nadir Kc in the Fig. 5 reproduction moved away from the published value.

This document checks the P terms, the `if(FALSE)` branches and the overlap function against the
paper, and against an independent Monte Carlo (MC) of the scene GOMS describes. It was written
against `gortt_brdf.c` as of 129ed80, and every "this branch" or "gortt" comparison below refers
to that version.

> **Status:** changes 1–3 below are made: `Po` in 0204222, dead code deleted in 5399f53, and the
> overlap height in 59121b9. `figures/goms_kernels.py` reproduces the resulting binary to
> 3.7e-13, and the MC scores it at Kc RMSE 0.080 (the "changes 1 + 2" row). After the papers in
> section 8, Γi in β was also changed to the illumination shadow area (c0ef2c3). That takes the
> MC score to 0.064, with a worst error of 0.33. β itself (section 6) is still open as
> OPEN-QUESTIONS #9, and the residual clamp (section 5) as #11.
>
> **Correction (section 9).** Change 1's rule, Po = min(PvMv, PiMi) everywhere, made Kc negative
> with the sun near overhead, which the sza 15–60 grid never sampled. It is replaced by Schaaf,
> Li & Strahler (1994) eq. 23A (ed690df), which is identical at sza ≥ 15. β for an overhead sun
> was also corrected to its limit of 1 (b59ce42). Where sections 2 and "Proposed C changes"
> describe Po = min(PvMv, PiMi) as the rule, read eq. 23A.

"L&S" is Li & Strahler (1992), IEEE TGRS 30(2) 276–292 (`papers/Li_and_Strahler_1992.pdf`).
Primes are the ellipsoid→sphere transform, as everywhere else in gortt.

---

## Summary

| # | finding | evidence | proposal |
|---|---|---|---|
| 1 | **Removing F² was right.** | Same-height crowns, where β plays no part: Kc RMSE against MC is 0.074 without F², 0.23 with it. | Keep 129ed80. |
| 2 | **`Po` is the cause of the unbounded Kc.** In forward scatter gortt sets `Po = PvMv`, which removes view-direction mutual shadowing entirely. | With `Po = min(PvMv, PiMi)`, Kc can provably never exceed 1 − Kg − Kz in the principal plane. The MC shows Po ≈ min(PvMv, PiMi) everywhere. | **Change 1.** Kc RMSE over all canopies 0.292 → 0.086, max error 6.6 → 0.58. |
| 3 | **The overlap uses the wrong height.** `p->h = 2b + h2 − h1` is the canopy depth; eqs 4–6 want the height of the crown centre. | Kg RMSE is 6–10× lower with the overlap averaged over centre heights. Depth produces a spurious dip on the inner side of the hotspot and a too-narrow hotspot in azimuth. | **Change 2.** |
| 4 | **The IGARSS P terms gortt uses are good.** Eqs 21–22 as printed, and the `gortt_kc_PP_only` TGARSS branch, are not. | ΔP = PvMv + PiMi − Po from gortt tracks the MC at L&S's own Fig. 3 parameters. The printed forms give Kc from −1.2 to 10.8. | Delete `gortt_kc_PP_only` (change 3). |
| 5 | **Interpolating f in azimuth (L&S eq. 25) is the best of the options tried.** | Off-PP RMSE 0.057, against 0.158 interpolating the P terms and 0.180 computing them at the actual azimuth. | No change. There is a small residual clamp near the hotspot off the PP. |
| 6 | **β is the main remaining error, and it is not a one-line fix.** | The MC's effective β is 0.8–1.0 in backscatter and 0.2–0.6 in forward scatter. gortt's β is symmetric in vza, and eq. 23 is a constant; neither can do that. | Open question (below). |
| 7 | **Ni et al. (1999) Fig. 5 is not a usable target for Kc.** | The MC gives nadir Kc 0.32 for SOBS against a published ~0.15, and forward-scatter Kc 0.18 against ~0.13. | Stop tuning to it. |

The recommended combination (changes 1 + 2, F² removed) takes Kc RMSE against the MC from
**0.174 on master** to **0.080**. The Kt clamp goes from 14 geometries to 19, all of them on the
principal cone near the hotspot, and the worst error drops from 6.6 on this branch (0.60 on
master) to 0.58. The full table is in section 2.

---

## Method

**`figures/goms_kernels.py`** is a numpy port of `gortt_overlap`, `gortt_kg`,
`gortt_kc_fFbeta`, `gortt_kc` and the Kz/Kt lines of `gortt_rsurf`, with every contested choice
behind a switch.
- With the switches at their defaults it reproduces `gortt -prnprop` to 3.5e-13 over 11,100
  geometries: four canopies, three β settings, principal plane and full azimuth, and both the
  F² and non-F² binaries.
- Run `figures/goms_kernels.py --selftest [--old gortt_master]` to repeat that check.

**`figures/goms_mc.py`** is the Monte Carlo. It builds the scene GOMS assumes:
- solid spheroids of radius r and half-axis b;
- centres uniformly random at density λ on periodic tiles of about 1,500 crowns;
- centre heights uniform in [h1, h2].

The MC is computed in sphere space, where a vertical scaling by r/b preserves projected areas,
occlusion and the sign of n·d. It uses two estimators:
- **rays:** parallel view rays, first hit, then a sun ray from the hit point. This gives Kc, Kt,
  Kg, Kz directly.
- **elements:** samples over crowns' view-projected discs. These tally the eq. 16 terms Gc/Gv,
  Mv, PvMv, PiMi, Po, plus Mi from the sun-projected disc.

By construction Kc = λπr²sec v′ (Gc/Gv − PvMv − PiMi + Po), and the two estimators agree.

The MC was validated as follows:
- The four proportions sum to 1, and Kt = Kz = 0 at the hotspot.
- In the sparse limit, Kc matches λπr²sec v′ (1 + cos ξ′)/2 and Kg matches the zero-overlap
  Boolean model, both within 1.3σ.
- MC Mi and Mv reproduce L&S eqs 10–11 to within 0.004.
- Changing the seed or quadrupling the tile moves the results by about 2σ or less.
- A Poisson crown count per tile made scene-to-scene λ wander by ±2.6%, so the count is fixed
  and results are averaged over 8 tiles. Standard errors come from the spread between tiles
  (typically 0.002–0.005 on Kc).

**`figures/kc_investigation.py`** runs about 1,300 MC geometries and writes the tables and plots
below to `figures/kc_investigation/`. The run takes about 4 minutes on 4 cores; the MC results
are cached in a CSV. The canopies are:
- SOBS (gortt's default)
- same-height: SOBS with h1 = h2 = 5.75, so β = 1
- sparse: λ = 0.1
- dense: λ = 1.0
- wide-spread: h1 = 3, h2 = 14, so h2 − h1 > 4b

Each canopy is run at sza 15/30/45/60° over:
- the principal plane, ±85° in 5° steps;
- the principal cone vza = sza, raa 0–180°;
- vza = 45° across azimuth.

L&S's Fig. 3/4 canopy is run as well.

```sh
figures/goms_kernels.py --selftest
figures/kc_investigation.py            # tables.md + five PNGs in figures/kc_investigation/
figures/goms_mc.py --vza 0 --sza 30    # one geometry, MC against gortt
```

---

## 1. F squared

The same-height canopy isolates it: β is 1 there under every β formula, so only the P terms and
F remain. RMSE of Kc against the MC:

| | same-height |
|---|---:|
| F² (master) | 0.230 |
| F² + overlap fix + Po fix | 0.240 |
| no F² + overlap fix + Po fix | **0.074** |

Across all five canopies the F² form is worse in every comparison once the Po fix is in
(0.170 against 0.080). The old F² suppressed the unbounded forward-scatter Kc, but it did so by
depressing Kc everywhere off the hotspot. See `kc_variants.png`: master's curve sits at about half
the MC across the whole forward side.

## 2. Po, and why Kc was unbounded

L&S eq. 17: `f = F [1 − Γv (PvMv + PiMi − Po)/Γc] / (1 − M)`. Since F = Γc/Γ and
1 − M = (1 − Kg)/(λΓ), this rearranges to

    f (1 − Kg) = λ (Γc − Γv ΔP),    ΔP = PvMv + PiMi − Po

so with β = 1, **Kc = λΓv (Γc/Γv − ΔP)**, independent of Kg.

The IGARSS PvMv gortt uses, `Mv − (1 − cos(v′cosφ − i′))/2`, satisfies an identity in the
principal plane, where Γc/Γv = (1 + cos(v′cosφ − i′))/2:

    Γc/Γv − PvMv = (1 + cos θMv)/2 = 1 − Mv

and λΓv(1 − Mv) = 1 − e^(−λΓv) = 1 − Kg − Kz. **So Kc ≤ 1 − Kg − Kz holds exactly whenever
ΔP ≥ PvMv, i.e. whenever Po ≤ PiMi.**

gortt's Po logic:
- **Backscatter beyond the hotspot:** `Po = PiMi`, so ΔP = PvMv. The bound holds with equality:
  all visible crown is sunlit. This is L&S eq. 20.
- **Between nadir and the hotspot:** `Po = PvMv`, so ΔP = PiMi. In every geometry tested here
  PiMi ≥ PvMv, so the bound holds, but that is observed, not proven.
- **Forward scatter:** `Po = PvMv`, so ΔP = PiMi. This follows L&S p. 282 ("PvMv … equal to Po
  in fact"). At large forward vza, PvMv grows past PiMi and the bound fails. Kc then grows like
  sec v′, reaching 3.8 at sza 15° and 7.5 for the dense canopy.

**The MC shows the smaller mutually-shadowed set is almost entirely inside the larger one**:
Po ≈ min(PvMv, PiMi).
- Over 1,235 geometries the median of Po / min is 0.94, and the RMSE of Po − min is 0.028.
- For same-height crowns at L&S's parameters the two agree to within 0.001.
- See `nesting.png`.

That is also exactly what gortt's case logic encodes wherever it works. `Po = min(PvMv, PiMi)`
agrees with the case logic everywhere except the grazing forward views where the case logic
breaks, and it restores the bound there.

**Kc against the MC, all variants.** RMSE per canopy, with sza 15–60° and every geometry pooled.
`clamp` counts geometries where Kc > 1 − Kg − Kz; the MC never has any.

| variant | SOBS | same-height | sparse | dense | wide-spread | all | max abs | clamp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| gortt, F² (master) | 0.162 | 0.230 | 0.067 | 0.227 | 0.129 | 0.174 | 0.603 | 14 |
| gortt (this branch, 129ed80) | 0.232 | 0.229 | 0.043 | 0.517 | 0.221 | 0.292 | 6.639 | 39 |
| + Po = min(PvMv, PiMi) (change 1) | 0.093 | 0.074 | 0.039 | 0.113 | 0.095 | 0.086 | 0.580 | 19 |
| + overlap h averaged (changes 1 + 2) | **0.084** | **0.074** | **0.033** | **0.107** | **0.084** | **0.080** | 0.580 | 19 |
| ... with β = eq. 23 | 0.093 | 0.074 | 0.030 | 0.149 | 0.118 | 0.101 | 0.499 | 16 |
| ... with β = 1 | 0.142 | 0.074 | 0.097 | 0.174 | 0.183 | 0.141 | 0.776 | 88 |
| overlap h averaged only, gortt Po | 0.229 | 0.229 | 0.038 | 0.516 | 0.217 | 0.290 | 6.639 | 39 |
| lune model (h averaged) | 0.094 | 0.118 | 0.034 | 0.162 | 0.092 | 0.108 | 0.749 | 23 |
| eqs 21–22 as printed | 0.500 | 0.582 | 0.145 | 1.119 | 0.444 | 0.642 | 10.797 | 63 |
| `gortt_kc_PP_only` forms | 0.682 | 0.856 | 0.339 | 1.021 | 0.586 | 0.734 | 10.393 | 40 |
| F² + changes 1 + 2 | 0.149 | 0.240 | 0.058 | 0.219 | 0.116 | 0.170 | 0.822 | 12 |

By direction, for SOBS: the Po fix removes almost all of the forward-scatter error at low sun
(RMSE at sza 15°: 0.934 → 0.181; at sza 30°: 0.241 → 0.189). At sza 45–60° it changes nothing
there, because PvMv never overtakes PiMi.

> **A correction to the plan this was investigated under.** I had derived a "lune model" of
> the P terms on the primed sphere. It predicted Po ≈ 0 in forward scatter and blamed the
> blow-up on Po being too large. The MC contradicts that: Po is as large as the smaller of the
> two terms. The lune model is kept in `goms_kernels.py` (`lune=True`) as a scored alternative
> (RMSE 0.108, worse than the fix), not as a recommendation.

## 3. The P terms and the `if(FALSE)` branch

`fig3_P_terms.png` reproduces L&S Figs 3–4 (λ = 30/900 m⁻², R = 3, b = 4.5, h = 25, θi = 30°) with
three columns: Pv, Pi, and ΔP, the only combination f depends on.

- **ΔP from gortt's IGARSS forms tracks the MC closely** for uniform heights (panel a), and
  within about 0.1 as the height spread grows (b, c), where it falls low in forward scatter. Since Po = min(PvMv, PiMi) over this whole
  range, gortt and the proposed fix coincide here.
- **L&S eqs 19–22 taken by case** (Pv = 1 between nadir and hotspot, Pi = 1 beyond) reproduce the
  shape of the published Fig. 3a. Its ΔP, however, is far too small everywhere except beyond the
  hotspot. Eq. 21's clamp to zero ("if the Mi boundary is too low") removes PiMi in forward
  scatter, where the MC has 0.15–0.25.
- **Pv and Pi individually are convention-dependent** and do not match the published MC (Fig. 4).
  L&S's simulation classified elements inside an intersecting crown as facing away from the sun.
  Ours counts them as both hidden and shadowed, which adds equally to PvMv, PiMi and Po, so ΔP is
  unaffected.

**`gortt_kc_PP_only`, the TGARSS '92 branch.** It is never called, and it is broken in four ways:
- `Po` is used uninitialised; it is set only in the `else` branch.
- `Pi` has no zero clamp.
- `Pv`'s cosine argument is `θMv − v′ + v′cosφ` while its condition is `θMv − i′ + v′cosφ`.
- There is no per-case logic, and β is never applied.

With Po supplied from the case logic so that it can be scored at all, it gives Kc RMSE 0.73 and
values up to 10.4. Eqs 21–22 exactly as printed give 0.64, from −1.2 to 10.8. The comment "This is
really not correct" is right in the sense suspected: both the coding bugs and the printed forms
themselves produce unphysical outputs.

## 4. The overlap function

`kg_overlap.png`. RMSE of Kg against the MC. The bracketed figure is the principal cone
vza = sza, where the azimuthal shape of the hotspot lives. The largest MC standard error on Kg is
0.005.

| variant | SOBS | same-height | sparse | dense | wide-spread |
|---|---:|---:|---:|---:|---:|
| **h = depth, ambrals (gortt)** | 0.0126 (0.0146) | 0.0018 (0.0020) | 0.0213 (0.0176) | 0.0029 (0.0038) | 0.0107 (0.0123) |
| h = depth, eq. 7 | 0.0119 (0.0124) | 0.0036 (0.0062) | 0.0208 (0.0154) | 0.0026 (0.0030) | 0.0103 (0.0109) |
| h = centre, ambrals | 0.0023 (0.0033) | 0.0012 (0.0013) | 0.0047 (0.0054) | 0.0006 (0.0009) | 0.0027 (0.0037) |
| **h = avg, ambrals** | **0.0013 (0.0017)** | **0.0012 (0.0013)** | **0.0028 (0.0031)** | **0.0005 (0.0008)** | **0.0012 (0.0014)** |
| h = avg, eq. 7 | 0.0024 (0.0040) | 0.0027 (0.0046) | 0.0033 (0.0045) | 0.0005 (0.0009) | 0.0020 (0.0033) |
| h = avg, pp (gortt's `else`) | 0.0247 (0.0406) | 0.0255 (0.0422) | 0.0531 (0.0828) | 0.0054 (0.0097) | 0.0235 (0.0393) |
| h = avg, wenge | 0.0049 (0.0088)* | 0.0049 (0.0087)* | 0.0083 (0.0144)* | 0.0013 (0.0022)* | 0.0037 (0.0064)* |

\* NaN at the hotspot itself, where `wenge` divides by D = 0.

**Height.**
- `p->h = 2b + h2 − h1` (`gortt.c:857`, "total canopy depth") is the only use of `p->h`, and it
  enters the overlap as h/b.
- L&S eqs 4–6 define h as the height of the spheroid's *centre*, and eq. 3 asks for the overlap
  *averaged* over the centre-height distribution.
- For SOBS, depth gives h/b = 4.04; the mean centre height gives 2.13. The effect is a hotspot
  about half as wide as it should be, with a spurious dip in Kg on its inner side (sza 30°,
  vza +15: 0.066 against MC 0.094 ± 0.001; the averaged overlap gives 0.094).
- For same-height crowns depth = 2b and centre = 5.75 are close, which is why that column barely
  changes.
- Averaging over [h1, h2] (24-point Gauss-Legendre) beats the mean centre height wherever the
  heights are spread.

**Azimuthal form.**
- In the principal plane, Ambrals `t2 = √(D² + (tan i′ tan v′ sin φ)²)` reduces exactly to
  L&S eq. 6, so the `if(TRUE)` choice matters only off the PP.
- There Ambrals beats both L&S's own eq. 7 (a linear interpolation to a cutoff azimuth) and
  gortt's `else` branch, which is the PP-only formula applied at all azimuths and is by far the
  worst.
- Wenge's `t1` variant is singular at the hotspot and worse than Ambrals everywhere else.

## 5. Azimuthal interpolation of f

`gortt_kc` evaluates f only at raa = 0 and π and interpolates it linearly (L&S eq. 25). RMSE of
Kc off the principal plane, all with the Po fix:

| | SOBS | same-height | sparse | dense | wide-spread | all | clamp |
|---|---:|---:|---:|---:|---:|---:|---:|
| **interpolate f (gortt)** | 0.047 | 0.087 | 0.016 | 0.065 | 0.042 | **0.057** | 19 |
| interpolate ΔP | 0.080 | 0.202 | 0.016 | 0.274 | 0.055 | 0.158 | 44 |
| ΔP at actual azimuth | 0.102 | 0.167 | 0.017 | 0.343 | 0.073 | 0.180 | 0 |

The IGARSS P terms are written with `cos φ` but extend poorly off the PP, so interpolating f is
right. It is the source of the **19 remaining clamps**. All of them are on the principal cone
vza = sza, 10–90° of azimuth from the hotspot and mostly at sza 15°. They occur only for same-height
and sparse crowns, with an excess of at most 0.10. f0 at the
hotspot is 1/(1 − Kg0), but it is multiplied by (1 − Kg) at the actual azimuth, where Kz has
already grown.

## 6. β, the main remaining error

Kc = β·Kc(uniform height) + (1 − β)·Kc(random). The MC therefore implies
β_eff = (Kc_MC − Kc_random)/(Kc_uniform − Kc_random). Medians over the principal plane:

| | back | nadir | fwd 25–60 | fwd 70–85 | gortt β (range over vza) | eq. 23 |
|---|---:|---:|---:|---:|---:|---:|
| SOBS, sza 15 | 0.88 | 0.73 | 0.47 | 0.45 | 0.38–0.91 | 0.24 |
| SOBS, sza 30 | 0.93 | 0.79 | 0.59 | 0.32 | 0.27–0.87 | 0.24 |
| SOBS, sza 60 | 0.99 | 0.81 | 0.82 | 0.46 | 0.21–0.83 | 0.24 |
| dense, sza 30 | 0.99 | 0.72 | 0.53 | 0.33 | 0.34–0.94 | 0.24 |
| wide-spread, sza 30 | 0.79 | 0.49 | 0.34 | 0.20 | 0.14–0.77 | 0.00 |
| same-height, sza 30 | 1.00 | 1.19 | 1.08 | 0.75 | 1.00 | 1.00 |

The full table is in `figures/kc_investigation/tables.md`.

- **The height spread matters in forward scatter and hardly at all in backscatter.** Beyond the
  hotspot you see sunlit crown tops whatever their heights. Looking toward the sun, the spread
  decides whether lower crowns' shadows fall on the crowns you see.
- **gortt's β** (the IGARSS formula, with Γi = Γv) grows with vza, but symmetrically, so it
  cannot separate the two sides. At nadir view it is 0.2–0.4, too low. That produces the notch
  just beyond the hotspot in `kc_variants.png`: the (1 − β)·F·(1 − Kg) random-case term falls
  sharply there, where the MC wants β ≈ 1.
- **Eq. 23** (β = 0.24 for SOBS, 0 for the wide-spread canopy) is better in forward scatter and
  worse in backscatter. Overall it is slightly worse (RMSE 0.101 against 0.080).
- **β = 1 everywhere** is exactly right for same-height crowns and wrong for every spread canopy
  (0.141, with 88 clamps off the PP).
- **Even for same-height crowns** the uniform-height model underestimates forward and near-nadir
  Kc at high sun (β_eff > 1 at sza 60°). It also saturates at Kc = 1 at grazing forward view with
  sza 15°–30°, where the MC gives 0.67–0.98. Some error is in the P terms themselves, not only in
  β.

A direction-dependent β fitted to the MC would remove most of what is left. It would be an
empirical change to a published model, though, so it is left as a decision.

## 7. Ni et al. (1999) Fig. 5

SOBS, sza 30°, principal plane:

| vza | MC | master (F²) | this branch | changes 1 + 2 | Fig. 5 |
|---:|---:|---:|---:|---:|---:|
| +60 | 0.973 | 0.712 | 0.856 | 0.868 | ~0.80 |
| +45 | 0.920 | 0.593 | 0.739 | 0.812 | ~0.79 |
| +30 (hotspot) | 0.814 | 0.817 | 0.817 | 0.817 | ~0.81 |
| 0 | 0.320 | 0.167 | 0.241 | 0.242 | ~0.15 |
| −40 | 0.183 | 0.059 | 0.143 | 0.143 | ~0.13 |
| −60 | 0.184 | 0.043 | 0.191 | 0.191 | ~0.15 |

Every model agrees with the MC at the hotspot, but the published figure disagrees with the MC by
2× at nadir. Master's nadir value was close to the figure only because F² depresses Kc everywhere.
The figure was presumably produced with Wenge's code, which carries its own choices. It tells us
about that code, not about the GO scene, and should not be used to choose between variants.

---

## 8. Addendum: what Li & Strahler (1992, IGARSS) and Abuelgasim & Strahler (1994) confirm

`papers/Li_and_Strahler_1992_IGARSS.pdf`, "Mutual shadowing and directional reflectance of a
rough surface", pp. 766–768, arrived after the sections above were written.

- **Change 1 is in the source.** Conclusion 1 and Fig. 1 (their own MC, uniform heights,
  principal plane): "it is quite accurate to model Po = PiMi at θv > θi; and Po = PvMv at θv < θi,
  until the large zenith in the forward direction is reached where PvMv > PiMi again". That is
  Po = min(PvMv, PiMi). gortt's old logic had kept Po = PvMv across the whole forward side.
- **The P terms gortt uses are the paper's.** Eq. 4 is gortt's PvMv. Eq. 5 is gortt's PiMi: a
  "fuzzy" Mi boundary replacing eq. 3, the TGRS eq. 21 form, which "guarantees the right result
  at the hotspot" and vanishes only when θv′ − θi′cosφ reaches π. That is why it has no zero
  clamp, and why it beat eq. 21 against the MC in section 3.
- **gortt's β is the paper's eq. 7, except that gortt sets Γi = Γv.** The paper defines Γi as
  the shadowing area cast on a level, i.e. the illumination direction. With Γi = πr² sec(i′)
  (`beta="illum"`, adopted in c0ef2c3), Kc RMSE over the five canopies goes 0.080 → 0.064. Forward scatter improves
  (0.130 → 0.060) and backscatter gets worse (0.065 → 0.095). Eq. 6, from which eq. 7 follows, is
  defined at nadir viewing, so the published β is view-independent by construction and cannot
  give the asymmetry of section 6. OPEN-QUESTIONS #9 stays open.
- **Abuelgasim & Strahler (1994), RSE 47:261–275, eq. 8** (`papers/Abuelgasim_and_Strahler_1994_RSE.pdf`)
  reprints the same β. It is the "RSE '94" variant that was deleted with `gortt_kc_PP_only`, and
  as printed it has both exponent signs flipped. The paper's text (β falls from 1 to 0 with
  canopy depth and rises to 1 with coverage) describes IGARSS eq. 7, not the printed form, which
  is negative at zero depth. It also calls λΓi the "coverage and sun position", which supports
  Γi being the illumination area. The paper does not say how β enters Kc.
- **The cosine-in-azimuth suggestion (conclusion 3) is not supported by our MC.** Using
  (1 − cos raa)/2 as the interpolation weight gives off-plane RMSE 0.061 against 0.057 for linear,
  with more clamps. Interpolating Kc rather than f removes the clamps but costs more: 0.079
  linear, 0.092 cosine.

---

## 9. Addendum: Schaaf, Li & Strahler (1994), and the sun near overhead

`papers/Schaaf_et_al_1994_TGRS.pdf`, IEEE TGRS 32(6) 1186–1193. The appendix restates the
Li–Strahler model in full.

- **Eq. 24A confirms Γi:** "Γi = πR² sec θi′ is the illumination shadow projected on the ground"
  (c0ef2c3). The β formula carries the same exponent-sign misprint as Abuelgasim & Strahler
  (1994).
- **Eq. 23A confirms that β weights f**, as gortt does, not the P terms:
  f = βF(1 − ΓvPvMv/Γc)/(1 − M) + (1 − β)F.
- **Eqs 5A–6A** define h as the height to the crown centre (change 2).
- **Eq. 23A also gives the Po rule, and it corrects change 1.** It keeps the decrement at PvMv
  "while θv′ > θi′; otherwise, the PvMv can be replaced with PiMi if it is a larger value". Read
  as backscatter beyond the hotspot, that is Po = PiMi there and Po = min(PvMv, PiMi) elsewhere.
  It is the old backscatter logic, plus the fix for large forward zeniths.
  - It differs from min(PvMv, PiMi) only when PiMi > PvMv beyond the hotspot, which happens only
    for sza below about 15°. There the min rule subtracts PiMi and **Kc goes negative**: −4.2,
    with BRF −1.9, for λ = 1 at sza 5, vza 80, and down to −59 on a 1° scan.
  - The main MC grid (sza 15–60) never sampled it. Neither did the test suite, which had no
    low-sun backscatter row; `lowsun` now pins it.
  - Read literally on both sides of nadir, the rule is much worse (Kc RMSE 0.153 on the main
    grid, 0.291 in forward scatter).
- **The overhead-sun guard in β was wrong.** `gortt_kc_fFbeta` returned β = 0 for sza′ < 1e-9.
  Eq. 7's limit there is 1, since D = r cot(sza′/2) → ∞, and β is already 0.99 at sza 0.1°.
  Fixed in b59ce42.
- **Eqs 7A–9A**, an exact overlap on the principal cone, differ from the Ambrals form near the
  hotspot by up to 1.85 (in units of πr²) at steep sun. Against the MC on the principal cone the
  two are indistinguishable: Kg RMSE 0.0015 against 0.0017 for SOBS, with MC standard errors of
  0.003–0.005. The off-plane approximation, eqs 11A–18A, was not tested.

The MC grid now also covers sza 0, 5 and 10: principal plane, plus vza = 45 across azimuth. Kc
RMSE against the MC on that grid (`tables.md`, "Low sun"; `kc_low_sun.png`):

| variant | SOBS | same-height | sparse | dense | wide-spread | all | max abs | Kc<0 | clamp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gortt at 129ed80 | 0.971 | 1.021 | 0.179 | 2.163 | 0.932 | 1.230 | 17.085 | 0 | 106 |
| Po = min, β(0) = 0 (8d5fe7c) | 0.168 | 0.204 | 0.073 | 1.235 | 0.135 | 0.569 | 10.602 | 11 | 0 |
| Po = eq. 23A, β(0) = 0 | 0.167 | 0.205 | 0.074 | 0.225 | 0.134 | 0.170 | 0.500 | 0 | 0 |
| **gortt now: eq. 23A, β(0) = 1** | **0.095** | **0.098** | **0.088** | **0.085** | **0.127** | **0.100** | **0.299** | **0** | **0** |

By solar zenith for gortt now: RMSE 0.138 at sza 0, 0.078 at 5, 0.084 at 10. The main grid is
unchanged at 0.064. At sza 0 the sparse canopy overshoots the MC at large view zeniths: β = 1 is
the formula's limit, but it is too large for a sparse canopy seen obliquely, which is the β
problem of section 6.

---

## Proposed C changes, in order

*All three are now made; the measured golden deltas are in their commit messages. Change 1
(0204222) moved only the six `-energy` cases, by up to 0.0027. Change 2 (59121b9) moved 37 of 55,
by up to 0.044. Change 3 (5399f53) was byte-identical.*

Each change is measured with `goms_kernels.py` over the sensitivity-sweep grid (sza 15–60,
principal plane ±80°, and vza 20/45/70 across azimuth; 480 geometries per canopy).

1. **`Po = MIN( PvMv, PiMi )` in `gortt_kc_fFbeta`** (superseded by Schaaf et al. 1994 eq. 23A
   in ed690df; see section 9), replacing the `raa`/`fabs(vza) > fabs(sza)`
   case logic.
   - Kc RMSE against the MC: 0.292 → 0.086. The clamp count drops from 39 to 19, and the largest
     error from 6.6 to 0.58.
   - Changes only the geometries where the case logic breaks: 19 of 480 for SOBS, 72 for dense,
     78 for the HB/BR/PCC canopy, none for sparse. Where it does change, the change is large
     (up to 2.9).
   - Expected goldens: forward-scatter, large-vza rows only.

2. **The overlap averaged over centre heights in [h1, h2]**, replacing `h = 2b + h2 − h1`.
   The Gauss-Legendre points already used elsewhere in gortt would do. The simplest version,
   h = (h1 + h2)/2, recovers most of the gain (Kg RMSE 0.0023 against 0.0013 for SOBS).
   - Kg RMSE improves 6–10× on every spread canopy.
   - Changes most geometries near the hotspot: max |ΔKg| 0.07 and max |ΔKc| 0.09 for SOBS.
     Expect nearly every BRF golden to move, and the gap-probability goldens not at all.
   - `p->h` is used nowhere else.

3. **Delete `gortt_kc_PP_only`, and the `wenge` and `pp` branches of `gortt_overlap`.**
   - None is reachable.
   - Each is scored above and is worse than, or broken relative to, what is used.
   - Their provenance is recorded here and in `goms_kernels.py`.
   - This is the alternative to OPEN-QUESTIONS #7's suggestion of promoting them to enums. The
     numbers say none is worth selecting.
   - It should be byte-identical on the goldens.

4. **Not proposed, needs a decision:** β (section 6), and the 19 residual off-PP clamps near the
   hotspot (section 5).

## Caveats

- The MC is the GOMS scene: solid crowns, a random Boolean field, uniform centre heights. It
  arbitrates the geometric-optical approximations, not how well GOMS describes a real forest.
  Within-crown gaps belong to the signatures, not to these proportions.
- MC "facing the sun" is n·s > 0. Elements inside an intersecting crown count as hidden and
  shadowed, which leaves ΔP, and therefore Kc, unaffected. It does make Pv and Pi individually
  differ from L&S's published Fig. 4.
- Sections 1–7 use sza 15–60. The sun near overhead (sza 0–10) was added later and has its own
  table in section 9.
- The azimuth experiments cover the principal cone vza = sza and vza = 45°. Grazing azimuthal
  geometries are sampled only in the principal plane.
- L&S's IGARSS '92 paper was not available when this was written. It has since been checked;
  see section 8.
