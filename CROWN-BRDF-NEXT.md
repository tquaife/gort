# The crown BRDF: where to pick this up

> **SUPERSEDED, 18 September 2026.** The work this document sets out has been done.
> Read **`CROWN-BRDF-INVESTIGATION.md`** instead; this file is kept for what it asked
> for and in what order, which is still the right order, and because the investigation
> refers back to its section numbers.
>
> In short: section 3's measurement was made and is decisive, but it does not say what
> this document expected. Option B (section 5) is ruled out — crowns look *less* like
> surfaces as they get thinner, not more. Option A (section 4) was implemented as
> `-rt_crown_brdf` and makes the FLIGHT agreement worse, because the crown's
> forward/backward asymmetry is geometric and cannot come from the leaf phase function:
> in the near infrared `r_L = t_L`, so that phase function is isotropic and `R_ss`
> reduces to a symmetric Lommel-Seeliger that invents a backscatter limb the crown does
> not have. The section 4 warning about double-counting the hotspot turned out not to be
> the problem; the hotspot being too *weak* is.

Written 17 September 2026, at the end of the session that produced v3.1, as a handover
to a later one. It is instructions, not a record: `KC-INVESTIGATION.md` is the model for
how the record should look once something has actually been done.

**Read this before touching `gortt_rt_init` or `gortt_crown_reflectance`.**

---

## 1. Where things stand

`v3.1` on `master`. `-rt_crown_path` takes the crown path reflectances from a crown
rather than from a layer holding the whole scene's foliage, chord-averaged over the
spheroid. Off by default. Against FLIGHT over 12 canopies × 3 solar zeniths × 64
directions it moves the near-infrared bias from −7.4% to +1.7% and the RMSD from 18.6%
to 16.7%; in the sparse scenes, −9.5% to +2.5% and 18.4% to 13.9%. Red moves by 0.8
points. `flight_compare/variants/index.html` ranks all the options; `latex_notes/`
carries the working note with the equations.

Three facts that were expensive to establish and must not be re-litigated:

* **`p->elai` is correct.** It is the paper's `τh` with `h = λV`, the vertical squeeze
  §III-C describes. Verified equation by equation against `papers/Ni_etal_1999.pdf` and
  by reproducing the paper's own Fig. 7 component signatures. Do not "fix" it.
* **The two canopies are matched.** `flight_compare/validate.py` puts black leaves over
  a white soil and agrees with `gortt -W` on the nadir gap probability to within 0.8%.
  Whatever else is wrong, the structure is not.
* **The default path must not move.** With no `-rt_*` option the model is the paper as
  published and the 82 original goldens are byte-identical. `rt-path-scene` pins it:
  `-rt_crown_path -rt_path scene` asks for the published depth, so its golden must stay
  byte-identical to `rt-crown-path-off`.

---

## 2. The problem that is left

`ρ'_df` is a **directional–hemispherical** reflectance being used where a
**bidirectional** one is needed, and it comes from a **plane-parallel** slab.

Three approximations are stacked in it, and they have different fixes:

1. **Angular** — one incoming direction, all outgoing directions, no dependence on
   where the sensor is. `-rt_crown_path` did not touch this.
2. **Geometric** — a slab with a flat horizontal top, infinite extent, no sides, and a
   well-defined "down". Chord-averaging made each chord its own slab; they are still
   slabs.
3. **Positional** — brightness varies across a crown and so does the envelope normal.
   Chord-averaging is a *partial* fix: it is the paper's own acknowledged "less
   brightness near the edges of the crown and increased brightness near the center"
   (§III-B, where it is explicitly neglected) put back in, but only the depth part.

Size of what is left: the chord average gives `ρ_df = 0.398` for the sparse sphere in
the NIR where FLIGHT says ≈0.32, a ratio of 0.80. A nadir BRF 20% below the DHR is an
unremarkable amount of bowl shape for a canopy, so the arithmetic is consistent with
this being the remaining term. **That is a consistency check, not evidence.** Section 3
is how to turn it into one.

### Why "just impose a BRDF on the crown envelope" is the wrong instinct

A BRDF is defined against a surface normal with a half-space of medium beneath it. The
crown envelope is not an interface — photons do not scatter *at* it, they scatter in the
volume and cross it. So imposing a BRDF means the **tangent-plane approximation**:
assume the escaping radiance equals the plane-parallel BRDF for the local normal and
local depth. That needs local optical depth ≫ 1 *and* radius of curvature ≫ mean free
path. Both fail at the limb, and the limb is where the area is:

| | fraction of projected crown area |
|---|---|
| `u > 0.7` | 51% |
| `u > 0.8` | 36% (chord < 0.60 of max) |
| `u > 0.9` | 19% (chord < 0.44 of max) |

(`u = s/r`; the projected-area weight is `2u du`, largest at the limb.)

And these crowns are not optically large. Crown optical radius `σr = k·FAVD·r`, i.e.
mean free paths across a crown, over the scene matrix:

| scene | σr | σb |
|---|---|---|
| dense-low-prolate | **0.23** | 0.58 |
| dense-low-sphere | 0.58 | 0.58 |
| sparse-low-sphere | 1.88 | 1.88 |
| dense-high-sphere | 2.31 | 2.31 |
| sparse-high-oblate | **15.0** | 7.5 |

A 60-fold range, and at the bottom the crown is a quarter of a mean free path across —
a wisp, with no surface to put a BRDF on. Note the prolate crowns are thin across
(σr = 0.23) and thicker along (σb = 0.58); one BRDF cannot carry that.

---

## 3. Do this first: measure the thing before modelling it

**Do not write any model code until this is done.** It is minutes, not hours, and it
turns the argument above into a target.

**A single crown, big black scene, FLIGHT image output.**

* One crown at the centre of a scene large enough that nothing else is in the footprint.
* Black soil (`write_spec(..., rsoil=0.0)`), so the only thing in the image is the crown.
* `OUTPUT_MODE` in `in_flight.data` — line 23, currently hardwired to `0` in
  `flight_compare/scene.py::write_in_flight`. Setting it non-zero makes FLIGHT write an
  image rather than just the scene-average spectrum. **Check what the modes are in
  `flight.c` before assuming**; the writer will need a parameter.
* Sweep view zenith in the principal plane at one or two solar zeniths.

What to get out of it:

1. **The brightness profile across the crown** — i.e. limb darkening, and how it changes
   with phase angle. This is the function the model has to reproduce.
2. **The crown's own phase function** — the integrated brightness against view angle,
   with no soil and no neighbours.

Why this is the right target: GORT currently treats the sunlit crown as a **disc of
uniform radiance**, no limb darkening at all. The candidate laws disagree sharply about
exactly that — a Lambertian sphere is limb-darkened as `√(1−u²)` at zero phase, while
Lommel–Seeliger is perfectly flat. So the measurement discriminates between them
immediately, and it is a far sharper target than the whole-scene BRDF, which convolves
the crown law with the Li–Strahler geometry and the background.

Do it over a range of `σr` — at minimum `dense-low-sphere` (0.58) and `sparse-high-oblate`
(15.0) — because the whole question is how the answer degrades as the crown stops being
optically large.

---

## 4. Option A (recommended): keep it a volume problem

**Do not put a BRDF on the envelope. Give the volume source function the angular
structure and integrate along the chord in the view direction.** This sidesteps the
"edge of a spheroid" question entirely, because there is no surface — and it degrades
gracefully, since at the limb the chord is simply short.

Single scattering exactly, multiple scattering isotropically (Hapke's IMSA construction):

```
R(L, θi, θv, φ) = R_ss + R_ms

R_ss = ω·P(g) / (4(μi + μv)) · [ 1 − exp( −kL (1/μi + 1/μv) ) ]

R_ms   chosen so that   ∫ R dΩ  ≡  ρ_df(L)
```

Why this one:

* **`ρ_df` becomes the normalisation, not the answer.** Energy is conserved by
  construction and the two-stream keeps doing what it is good at (how much is
  scattered) while the single-scattering term does what it is good at (where it goes).
* **It is well matched to the two bands.** In the NIR (ω = 0.9) multiple scattering
  dominates the magnitude but is genuinely closer to isotropic; in the red (ω = 0.12)
  single scattering dominates and is treated exactly. The band that needs the anisotropy
  most is the band where the approximation is best.
* **It drops into the existing quadrature.** In `gortt_crown_reflectance`
  (`gortt_brdf.c`), replace `ρ_df(L(u))` with `R(L(u), θi, θv, φ)` under the same
  `2u du` chord integral on the same Gauss–Legendre rule. Structurally one for one.
* **The inputs already exist.** The asymmetry `g = −(4/9)(rL − tL)/ω` is already computed
  in `gortt_rt_init` as `rt->g_asym`, so a Henyey–Greenstein `P(g)` or the existing
  bi-Lambertian form needs no new user input.

Watch for: `rt->g_asym` is guarded against ω = 0 (added in v3.1) — keep that. And the
hotspot term in eq. (30) already carries a `p(n>0)` bidirectional gap probability, so
check you are not double-counting the single-scattering peak when `R_ss` goes in.

---

## 5. Option B: tangent-plane BRDF over the envelope

Worth knowing about, but not where to start.

For a **sphere** with a local reflectance law, the integral over the illuminated-and-
viewed cap has classical closed forms — the Lambert sphere and Lommel–Seeliger sphere
phase functions are textbook. GORT's prime transform already maps the spheroid to a
sphere, so it would come out in primed space for free, and Li & Strahler's geometry
already knows which part of the envelope is sunlit and viewed (that is what `Kc` is —
the *area* of that region; what is missing is the *angular* weighting over it).

Why not first: it inherits the limb problem in §2 wholesale, and it forces a choice of
local law that Option A makes for you. Revisit only if §3's measurement says the crowns
behave like surfaces — i.e. if the limb-darkening profile matches a surface law across
the σr range rather than flattening out as the crowns get thin.

---

## 6. Things not to do

* **Do not re-run FLIGHT for the BRDF or lidar comparison.** `run_brdf.py --gortt-only`
  and `run_lidar.py --gortt-only` write the gortt CSV and return before the FLIGHT
  binary is even checked. `results/lidar_flight.csv` has **no checkpoints** and is the
  one irreplaceable artefact; `work/baseline_results/` holds a copy, keep it that way.
* **Do not revisit the recollision probability without new physics.** It was tried and
  measured. `run_flight_single.py --black-soil` gives p from FLIGHT directly, and the
  Wigner estimate `σl̄/(1+σl̄)` is right for thin crowns (0.435 vs a measured 0.37–0.46)
  and far too high for thick ones (0.909 vs 0.53–0.57), because first-flight escape
  assumes scattering is spread uniformly through the body when in a thick crown it
  happens near the lit surface. `-rt_recollision` is kept *for that negative result*.
* **Do not turn anything on by default.** Whatever comes of this goes behind a flag
  alongside the other `-rt_*` options, for the same reason: the goldens are the guard.
* **Do not touch the prime-transform question in the same change.** `gortt_rt_init`
  evaluates the RT chain at the primed zenith where the paper writes `μi`
  (the `05/06/06` comments; see `OPEN-QUESTIONS.md`). It interacts with all of the
  above and wants isolating, not bundling.

---

## 7. Practical

```sh
make && make test                 # 97 cases; the 82 originals must not move
python3 figures/fig67_ni1999.py --gortt-arg=-rt_crown_path --suffix=-crown-path
cd flight_compare
cp results/*.csv work/baseline_results/          # before anything
./run_brdf.py --gortt-only && ./plot_brdf.py && ./make_html.py --variants
```

Adding a variant: `flight_compare/config.py::GORTT_VARIANTS` is the single source of
truth; the runners, plotters and page builder all loop over it.

Reading:

* `latex_notes/crown_path_notes.tex` — the working note: equations, the three
  approximations, §7 "what I am still unhappy about", which is this document's §2.
* `latex_notes/crown_path_figs.py` — its three figures, regenerable.
* `flight_compare/variants/index.html` — every option ranked against FLIGHT.
* `flight_compare/README.md` — how the two models were matched and what was worked
  around in FLIGHT.
* `papers/Ni_etal_1999.pdf` — §III-B for the neglected non-uniformity, §III-C for the
  squeezing passage, §V for the admission that only dense forests were validated.
