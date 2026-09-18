# The crown BRDF: what a crown actually looks like, and why Option A failed

Written 18 September 2026, following `CROWN-BRDF-NEXT.md`, which set out what to do and
in what order. This is the record of doing it. `KC-INVESTIGATION.md` is the model for the
form.

**Summary.** The measurement in section 3 of the handover was made and it is decisive, but
it does not say what the handover expected. The recommended fix, Option A — single
scattering exactly, multiple scattering isotropically, normalised to `rho_df` — was
implemented as `-rt_crown_brdf` and **makes the agreement with FLIGHT worse**, in every
phase-angle bin and in both bands. It is kept, off by default, for the negative result, the
way `-rt_recollision` is. The measurement says why it fails and points somewhere else.

---

## 1. What was measured

`flight_compare/run_crown_image.py`. One crown at the centre of a black scene, FLIGHT in
reverse mode with `OUTPUT_MODE 2`, and the per-pixel radiance factor binned against the
impact parameter `u` of the view ray after the ellipsoid-to-sphere transform — the same `u`
as the chord average in `gortt_crown_reflectance`. Five canopies spanning the crown optical
radius `sigma_r = k favd r` from 0.23 to 15, two bands, nineteen directions in the principal
plane at a 30 degree solar zenith. `./plot_crown_image.py` draws it;
`./plot_crown_image.py --check` verifies the machinery the analysis leans on.

Two things had to be got right before any of it meant anything, and both are pinned in the
script:

* **The box has to be sized per view zenith.** FLIGHT is periodic (`replace_if_out`,
  `flight.c:943`), and a ray that leaves a wall comes back on the opposite one at the same
  height and can meet the crown a second time. It can do so whenever
  `hw < r + 2b tan(vza)`. In a fixed 12 m box that put 299 lit pixels out to `u = 2.5` at
  `vza` 60 for the prolate crowns — a clean second image of the crown, which would have
  read as a broad bright halo rather than as an artefact.
* **The geometry is analytic, not fitted.** Each pixel is a known ray, so `u` is computed,
  not registered. Every run asserts that no lit pixel lies outside the predicted silhouette;
  the tolerance is 0.01 in `u`, for rays that graze the limb, against the 2.5 a real fault
  produced.

---

## 2. Finding 1: crowns are strongly limb-darkened, and it scales with optical radius

Ratio of the outermost ring of the projected disc to the innermost, at zero phase:

| scene | `sigma_r` | red | nir |
|---|---|---|---|
| dense-low-prolate | 0.23 | 0.250 | 0.240 |
| dense-low-sphere | 0.58 | 0.270 | 0.273 |
| sparse-low-sphere | 1.88 | 0.483 | 0.429 |
| dense-high-sphere | 2.31 | 0.579 | 0.484 |
| sparse-high-oblate | 15.0 | 0.939 | 0.720 |

A uniform disc — what the published model assumes — is 1.000. A Lambert sphere is 0.222.
So a thin crown is almost exactly Lambert-sphere limb-darkened and a thick one is almost
exactly uniform.

**The reason is not a surface law.** For a thin crown the brightness tracks the chord
length, and the chord through a sphere at impact parameter `u` is `2r sqrt(1-u^2)` — the
same shape a Lambert sphere gets from `cos(incidence)`, for an unrelated reason. As the
crown thickens the medium saturates and the chord stops mattering.

**This rules out Option B of the handover.** Its section 5 said to revisit the
tangent-plane envelope BRDF "only if section 3's measurement says the crowns behave like
surfaces ... rather than flattening out as the crowns get thin". The measurement says the
opposite: they flatten as the crowns get *thick*, and the surface-like profile appears
exactly where a surface law is least defensible. Do not go back to Option B.

**It also says the chord average was more right than the handover credited.** Weighted RMS
of the normalised profile against FLIGHT, over all five scenes and both bands:

| law | RMS |
|---|---|
| uniform disc (published) | 0.192 |
| `rho_df(L(u))`, i.e. `-rt_crown_path` | 0.097 |
| `R_ss(L(u))` alone | 0.122 |
| Option A | 0.088 |

The chord average halves the error on its own, and has the right behaviour in both limits
by construction. It does not go far enough — it predicts limb ratios of 0.40 to 0.61 where
FLIGHT says 0.24 to 0.27 for the thin crowns — but the residual is geometric, not angular:
near the limb the slab picture is wrong because photons leave the side of a laterally finite
sliver, and no amount of angular structure in the source function fixes that.

---

## 3. Finding 2: the crown's phase function is asymmetric, and it is not the leaf's

The disc mean varies by a factor of 2.7 with view angle, but most of that is the crown's own
shadow: for a sphere the sunlit fraction of the projected disc is exactly
`(1 + cos g)/2`. Dividing it out leaves the mean radiance factor of the *sunlit* crown,
normalised to nadir (positive `vza` is backscatter, sun at 30 degrees):

```
dense-low-sphere  nir    vza    -70    -40    -25      0    +30    +40    +70
  FLIGHT sunlit crown           1.65   1.09   1.00   1.00   1.30   1.04   0.99
  GORT C, baseline              1.05   1.00   1.00   1.00   1.04   1.03   1.07
  GORT C, -rt_crown_brdf        1.39   1.09   1.03   1.00   1.08   1.11   1.42

dense-high-sphere red   vza    -70    -40    -25      0    +30    +40    +70
  FLIGHT sunlit crown           0.70   0.75   0.83   1.00   1.90   1.20   0.92
  GORT C, baseline              0.99   0.99   0.99   1.00   1.08   1.05   1.03
  GORT C, -rt_crown_brdf        0.94   0.86   0.89   1.00   1.25   1.28   1.49
```

Three things in there:

1. **A narrow hotspot spike of 30 to 90 per cent** at zero phase, gone within 5 degrees of
   it. GORT has a hotspot term (Kuusk) and it is far too weak: 1.08 against a measured 1.90
   in the dense-high red, 1.00 against 1.47 in its near infrared.
2. **A strong forward/backward asymmetry away from the hotspot.** The sunlit crown brightens
   toward the forward-scatter limb (1.65, 1.28) and does not brighten at all toward the
   backscatter limb (0.99, 0.97) — and in the red it *darkens* toward both.
3. **Away from those two, the sunlit crown really is close to Lambertian** in the near
   infrared, within about 10 per cent over the whole mid-range. The published model's one
   number for `C` is not a bad description of the part of the hemisphere that is neither
   hotspot nor limb.

**Point 2 is what kills Option A.** In the near infrared `r_L = t_L = 0.45`, so
`g_asym = -(4/9)(r_L - t_L)/omega` is exactly zero, the Henyey-Greenstein phase function is
identically 1, and `R_ss` reduces to pure Lommel-Seeliger. Lommel-Seeliger brightens the
limb *symmetrically*: 1.39 forward and 1.42 back, against a measurement of 1.65 and 0.99.
It gets the forward limb roughly right and invents the backward one.

So the crown's forward/backward asymmetry cannot come from the leaf phase function, because
in the band where the asymmetry is largest the leaf phase function is isotropic. It has to be
geometric — the correlation between which volume elements are lit and which are seen, which
is what a hotspot function is for, and which in GORT is the Kuusk term that finding 1 above
says is already too weak.

---

## 4. What Option A does to the FLIGHT comparison

`-rt_crown_brdf`, over 12 canopies x 3 solar zeniths x 64 directions. Relative bias and
RMSD, per cent:

| variant | nir bias | nir RMSD | red bias | red RMSD | nir sparse RMSD |
|---|---|---|---|---|---|
| baseline | -7.4 | 18.6 | +0.2 | 14.9 | 18.4 |
| crown-path | +1.7 | **16.7** | +1.0 | 15.2 | **13.9** |
| crown-brdf | +3.1 | 21.9 | +2.7 | 19.6 | 18.6 |

Worse than the baseline in the red, and worse than `-rt_crown_path` everywhere. By
phase-angle bin, near infrared:

| variant | 0-10 | 10-25 | 25-45 | 45-70 | 70-120 |
|---|---|---|---|---|---|
| baseline | 17.1 | 14.0 | 14.8 | 17.2 | 27.7 |
| crown-path | 13.8 | 12.5 | 13.9 | 16.0 | 24.4 |
| crown-brdf | 17.1 | 18.9 | 21.9 | 22.8 | 25.0 |

The damage is spread across the mid-range, not concentrated at the hotspot — so it is **not**
the double-counting the handover warned about. It is the Lommel-Seeliger shape being wrong,
exactly as section 3 predicts. The one bin it helps is the extreme forward scatter, which is
the one place the measurement says the crown really does brighten.

---

## 5. What was implemented, and what it is worth keeping for

`-rt_crown_brdf`, requiring `-rt_crown_path`, off by default. Inside the chord average it
replaces `rho_df(L)` with

```
R = R_ss + R_ms,    R_ms isotropic and chosen so that (1/pi) int R mu_v dOmega = rho_df(L)
R_ss = omega P(Theta) / (4 (mu_i + mu_v)) [ 1 - exp( -k L (1/mu_i + 1/mu_v) ) ]
```

with `P` Henyey-Greenstein at `rt->g_asym`. Energy is conserved by construction, so it can
only move radiation between directions. `rho'_ff` is untouched — it is
hemispherical-hemispherical and has no direction to get right.

It is kept because the negative result is worth having in a form somebody can re-run, and
because the machinery is right even though the law is not: if a better angular law is tried,
it goes in the same place, under the same normalisation, on the same quadrature.

**Verifications, all reproducible:**

* `rt-crown-brdf-black` is byte-identical to `rt-crown-path-black`. With black leaves `R_ss`
  is proportional to `omega` and `rho_df` is zero, so the option must be exactly inert; this
  is its control, the way `rt-path-scene` is `-rt_crown_path`'s.
* The 82 published-path goldens did not move. Only `usage`, `rt-all` and `rt-energy` did, and
  the last two only because `-rt_all` now includes the new option.
* **Energy conservation holds to 1e-11** for the spherical crowns: the hemispherical integral
  of `C` is unchanged to that precision by switching the option on, which is what the
  normalisation promises. `./plot_crown_image.py --check` runs it.
* `-nquad 1` to `-nquad 8` moves the BRF by 3.2e-4 with the option on against 1.8e-4 with
  only `-rt_crown_path`, so the added azimuth quadrature is not the binding one.

**And one thing the verification turned up.** For the *non*-spherical crowns the energy
check is off by 2 to 17 per cent, and that is not a fault in this change: the normalisation
is exact in the primed frame, which is where `gortt_two_stream` and the rest of the RT chain
are evaluated, while the BRF is reported in the unprimed one. Before this change `rho'_df`
had no view dependence at all, so the frame it lived in could never show up in an integral
over view angle. **This is the first quantitative handle on the prime-transform question in
`OPEN-QUESTIONS.md`**, which until now had "nothing here to say which way is right".

---

## 6. Where to pick it up

In order of what the measurement supports:

1. **The hotspot, not the phase function.** Finding 2 says the sunlit crown is close to
   Lambertian except for a narrow spike that GORT under-predicts by a factor of up to 1.8,
   and a forward limb it does not have at all. That is a much better-posed target than a
   general BRDF, and `gortt_kuusk` is where it lives. Note this collides with
   `OPEN-QUESTIONS.md` #9 (beta) and #4 (`p->k` at the wrong solar zenith), which are in the
   same kernels.
2. **The prime frame**, now that section 5 gives a number to it.
3. **The limb residual**, if anybody wants it: the chord average is half the error and the
   rest is photons escaping the side of a thin sliver, which is a geometric problem about
   laterally finite slabs and not a source-function problem. It is worth the least of the
   three — the profile across the disc enters the BRF only through its mean, and the mean is
   what `-rt_crown_path` already gets right.

Do **not** go back to Option B; see section 2.

---

## 7. Files

| | |
|---|---|
| `flight_compare/run_crown_image.py` | the measurement; `--dry-run`, `--quick`, about 35 min in full |
| `flight_compare/plot_crown_image.py` | the figures and the summary table; `--check` verifies the two-stream port and the energy normalisation |
| `flight_compare/figs/crown-profile-phase0.png` | finding 1, every scene on one pair of axes |
| `flight_compare/figs/crown-phase.png` | finding 2 |
| `flight_compare/variants/crown-brdf.html` | section 4, against FLIGHT |
| `gortt_brdf.c` | `gortt_hg_phase`, `gortt_cos_scatter`, `gortt_r_ss`, `gortt_rho_ss`, `gortt_crown_brdf` |
| `latex_notes/two_stream_derivation.tex` | where eqs (8)-(11) come from, which is the layer this all sits on |
