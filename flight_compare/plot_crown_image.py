#!/usr/bin/env python3
"""
Draw what run_crown_image.py measured: how a single crown's brightness varies
across its projected disc, and how the disc mean varies with view angle.

Three figures:

  crown-profile-phase0.png   the profile across the disc at zero phase, every
                             scene on one pair of axes (one per band), with the
                             three candidate laws drawn over it.  This is the
                             discriminator: at zero phase a Lambert sphere is
                             limb-darkened as sqrt(1 - u^2), Lommel-Seeliger is
                             perfectly flat, and GORT's uniform disc is also
                             flat, so the measurement separates the first from
                             the other two immediately.

  crown-profile-<scene>.png  the same profile at every phase angle sampled, one
                             panel per band, coloured by phase.  This is how the
                             limb darkening changes as the sun moves off the
                             view direction, which is what a model of it has to
                             reproduce as well as the zero-phase shape.

  crown-phase.png            the disc mean against signed view zenith, i.e. the
                             crown's own phase function with no soil and no
                             neighbours.  Positive view zenith is backscatter,
                             following plot_brdf.py and the GORT papers.

Everything is normalised by the disc mean where a shape is wanted and left in
radiance-factor units where a magnitude is wanted; each panel says which.

    ./plot_crown_image.py
"""

import argparse
import collections
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import scene as sc

BAND_COLOUR = {"red": "#c0392b", "nir": "#2a78d6"}
LAW_STYLE = {"Lambert sphere": ("#111111", "-"),
             "Lommel-Seeliger / GORT disc": ("#888888", "--")}


def read_profile(path):
    """(scene, band, sza, vza, raa) -> (u, brf, sem, frac_lit, phase)"""
    out = collections.defaultdict(list)
    with open(path) as fp:
        for r in csv.DictReader(fp):
            key = (r["scene"], r["band"], float(r["sza"]), float(r["vza"]), float(r["raa"]))
            out[key].append((float(r["u"]), float(r["brf"]), float(r["sem"]),
                             float(r["frac_lit"]), float(r["phase"]),
                             float(r["n"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def read_phase(path):
    """(scene, band, sza) -> list of (vza_signed, phase, brf, sem)"""
    out = collections.defaultdict(list)
    with open(path) as fp:
        for r in csv.DictReader(fp):
            vza, raa = float(r["vza"]), float(r["raa"])
            signed = vza if raa == 0.0 else -vza
            out[(r["scene"], r["band"], float(r["sza"]))].append(
                (signed, float(r["phase"]), float(r["brf"]), float(r["sem"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def disc_mean(rows):
    """
    The projected-disc mean rebuilt from the binned profile.

    Weighted by the pixel count of each ring, which *is* the projected-area
    weight 2u du: the pixels are a regular grid in the projection plane, so
    counting them is measuring area.  Using the count rather than u itself also
    gets the partly-filled outermost bin right.
    """
    n, brf = rows[:, 5], rows[:, 1]
    return float(np.sum(n * brf) / np.sum(n))


def two_stream(w, k, L, mu):
    """
    A numpy port of gortt_two_stream (gortt_brdf.c), returning rho_df.

    Same role as figures/goms_kernels.py is to the geometric-optical kernels:
    a readable copy of the C, so the implied *shape* of a model choice can be
    drawn without instrumenting the binary.  Checked against gortt itself by
    check_two_stream(), which --check runs.
    """
    w = np.asarray(w, dtype=float)
    L = np.asarray(L, dtype=float)
    gamma = np.sqrt(1.0 - w)

    t_0 = np.exp(-k * L / mu)
    R_inf_ff = (1.0 - gamma) / (1.0 + gamma)
    R_inf_df = (1.0 - gamma) / (1.0 + 2.0 * mu * gamma)
    T_inf_ff = np.exp(-2.0 * gamma * k * L)
    T_inf_df = ((w / 2.0) * (1.0 + 2.0 * mu)
                / (1.0 - (2.0 * gamma * mu) ** 2) * (T_inf_ff - t_0))

    denom = 1.0 - (T_inf_ff * R_inf_ff) ** 2
    t_ff = T_inf_ff * (1.0 - R_inf_ff ** 2) / denom
    bracket = t_0 * R_inf_df + T_inf_df * R_inf_ff
    return R_inf_df - t_ff * bracket


def prime_theta(scn, za_deg):
    """
    gortt_prime_theta: atan((b/r) tan(theta)), the ellipsoid -> sphere
    transform.  gortt_two_stream is always called at the *primed* solar zenith
    (gortt_brdf.c), so anything comparing against it has to be as well -- using
    the unprimed angle here was wrong by 25 degrees for the prolate crowns, and
    check_two_stream caught it.
    """
    return np.degrees(np.arctan((scn.b / scn.r) * np.tan(np.radians(za_deg))))


def chord_profile(scn, band, u, sza):
    """
    What GORT's chord average implies the profile across the disc should be:
    rho_df evaluated at the depth of the chord at each u.

    This is the shape -rt_crown_path already carries, because the chord at
    impact parameter u is 2b sqrt(1 - u^2) and gortt_crown_reflectance
    integrates rho_df over exactly that with the area weight 2u du.  Nothing
    new is being proposed here -- the point is to see whether the shape the
    v3.1 option already has is the shape FLIGHT shows.
    """
    w = band.rleaf + band.tleaf
    L = scn.favd * 2.0 * scn.b * np.sqrt(np.maximum(0.0, 1.0 - u * u))
    return two_stream(w, 0.5, L, np.cos(np.radians(prime_theta(scn, sza))))


def gortt_crown(scn, band, sza, signed_vzas):
    """
    What GORT says the same crown looks like: the sunlit-crown signature C and
    the viewed proportions Kc and Kt, over a black soil.

    With rsoil = 0 every background term drops out -- G, Z, CdG, CfG, CdCG and
    CfCG are all identically zero, and so is the shaded-crown signature T, since
    T is CdCG/CfCG -- which leaves C = rho'_df + the hotspot term and a scene
    reflectance of exactly Kc*C.  That is what makes this comparable with the
    measurement at all.

    At zero phase the comparison is exact: Kt is ~1e-13 there, so every bit of
    crown in view is sunlit and GORT's C *is* the disc mean it predicts.  Off
    the hotspot GORT sees a shaded fraction Kt/(Kc+Kt) which it paints black,
    while FLIGHT's shaded crown is lit by scattering within the crown; the
    dashed curve is that mixture and the solid one is C alone.  Note also that
    Kc and Kt are the proportions of a *stand* of these crowns, with the mutual
    shadowing that goes with it, where the measurement has one crown on its
    own -- another reason the hotspot is the place to read this.
    """
    rows = [(abs(v), 0.0, sza, 0.0 if v >= 0 else 180.0) for v in signed_vzas]
    header = "%d 1 %g\n" % (len(rows), band.wl)
    body = "".join("%.10g %.10g %.10g %.10g\n" % r for r in rows)
    argv = ([config.GORTT_BIN] + scn.gortt_geometry()
            + sc.gortt_spectral_args(band, rsoil=0.0) + ["-prnspec", "-prnprop"])
    out = sc._gortt(argv, header + body)

    c, kc, kt = [], [], []
    for line in out.splitlines()[1:]:
        if "{" not in line:
            continue
        spec = line.split("{")[1].split("}")[0].split()
        prop = line.split("[")[1].split("]")[0].split()
        c.append(float(spec[0]))
        kc.append(float(prop[0]))
        kt.append(float(prop[2]))
    return np.array(c), np.array(kc), np.array(kt)


def check_two_stream(scenes, sza):
    """
    Check the numpy port against the binary.

    gortt does not print rho'_df on its own, but over a black soil it prints
    C = rho'_df + (1-w) kuusk w (1-g) / (2 mu_i' mu_v'), and the Kuusk hotspot
    factor is a product of probabilities times H, bounded by 0 and its value at
    the hotspot.  So the residual C - rho_df(port) has to sit between zero and
    that coefficient times the largest H can be.  A port with the wrong algebra
    puts the residual outside the interval: sign errors and a wrong R_inf_df
    were both caught this way while this was being written.

    The bound is loose by construction -- it is a bracket, not an equality --
    but it is over two bands and every scene in the matrix at once, and the
    two-stream terms are the large ones in the near infrared, where rho_df is
    0.2 to 0.4 and the whole hotspot coefficient is 0.06.
    """
    bad = 0
    for scn in scenes:
        for band in config.BRDF_BANDS:
            w = band.rleaf + band.tleaf
            g = -(4.0 / 9.0) * (band.rleaf - band.tleaf) / w if w > 0 else 0.0
            # p->elai is favd * lambda * V, the depth the published model uses.
            elai = scn.favd * scn.lam * scn.crown_volume
            mu = np.cos(np.radians(prime_theta(scn, sza)))
            rho = float(two_stream(w, 0.5, elai, mu))
            c = gortt_crown(scn, band, sza, [sza])[0][0]
            hot_max = (1.0 - w) * w * (1.0 - g) / (2.0 * mu * mu)
            resid = c - rho
            ok = -1e-9 <= resid <= hot_max + 1e-9
            bad += not ok
            print("  %-22s %4s  rho_df(port) %7.4f  C %7.4f  residual %8.5f  "
                  "in [0, %.5f]  %s" % (scn.name, band.name, rho, c, resid,
                                        hot_max, "ok" if ok else "FAIL"))
    print("  %s\n" % ("two_stream port agrees with the binary" if not bad
                      else "%d FAILURES -- the port does not match gortt" % bad))
    return bad == 0


def check_energy(scenes, sza, nmu=40, nphi=48):
    """
    Check that -rt_crown_brdf conserves energy, which is the one thing its
    construction promises: R_ms is set so that the hemispherical integral of R
    is exactly rho_df, so switching the option on must move radiation between
    directions and never create any.

    The observable is (1/pi) int C mu_v dOmega_v over the view hemisphere with
    a black soil.  C is rho'_df plus the Kuusk hotspot term, and the hotspot
    term is untouched by the option, so that integral has to come out
    identical with it on and off.

    It does, to 1e-11, for the spherical crowns.  For the others it is off by
    a few per cent, and that is not this option's doing: the normalisation is
    performed in the primed frame, where gortt_two_stream and the whole RT
    chain are evaluated, while the BRF is reported in the unprimed one.  The
    two frames coincide only when b = r.  Before -rt_crown_brdf, rho'_df had
    no view dependence at all, so which frame it lived in could not show up in
    an integral over view angle -- see OPEN-QUESTIONS.md on the prime
    transform, and CROWN-BRDF-INVESTIGATION.md section 5.
    """
    x, wq = np.polynomial.legendre.leggauss(nmu)
    muv, wmu = 0.5*(x + 1), 0.5*wq
    vzas = np.degrees(np.arccos(muv))
    phis = np.degrees((np.arange(nphi) + 0.5) * 2.0*np.pi/nphi)
    dphi = 2.0*np.pi/nphi

    def hemi(scn, band, extra):
        rows = [(v, 0.0, sza, a) for v in vzas for a in phis]
        out = sc._gortt([config.GORTT_BIN] + scn.gortt_geometry()
                        + sc.gortt_spectral_args(band, rsoil=0.0)
                        + ["-prnspec"] + list(extra),
                        "%d 1 %g\n" % (len(rows), band.wl)
                        + "".join("%.10g %.10g %.10g %.10g\n" % r for r in rows))
        c = np.array([float(l.split("{")[1].split("}")[0].split()[0])
                      for l in out.splitlines() if "{" in l]).reshape(nmu, nphi)
        return float(np.sum(wmu * muv * np.sum(c, axis=1) * dphi) / np.pi)

    print("  %-22s %4s %12s %12s %10s %s"
          % ("scene", "band", "crown-path", "+crown_brdf", "rel diff", ""))
    bad = 0
    for scn in scenes:
        spherical = abs(scn.b - scn.r) < 1e-12
        for band in config.BRDF_BANDS:
            a = hemi(scn, band, ["-rt_crown_path"])
            b = hemi(scn, band, ["-rt_crown_path", "-rt_crown_brdf"])
            d = abs(b - a) / a
            # Only the spherical crowns are a test.  For the others the two
            # frames genuinely differ and the difference is the open question,
            # not a failure, so it is reported and not asserted on.
            if spherical:
                note = "ok" if d < 1e-8 else "FAIL"
                bad += d >= 1e-8
            else:
                note = "b != r: primed frame, see docstring"
            print("  %-22s %4s %12.8f %12.8f %10.1e  %s"
                  % (scn.name, band.name, a, b, d, note))
    print("  %s\n" % ("energy is conserved where the two frames coincide" if not bad
                      else "%d FAILURES -- the normalisation is broken" % bad))
    return bad == 0


def zero_phase_key(profiles, scene, band, sza):
    """The sampled direction closest to zero phase for this scene and band."""
    cands = [(v[0, 4], k) for k, v in profiles.items()
             if k[0] == scene and k[1] == band and k[2] == sza]
    if not cands:
        return None
    return min(cands)[1]


def plot_phase0(profiles, scenes, sza, outdir):
    """The zero-phase profile, every scene, with the candidate laws over it."""
    bands = [b.name for b in config.BRDF_BANDS]
    fig, axes = plt.subplots(1, len(bands), figsize=(11, 4.6), sharey=True)

    u = np.linspace(0.0, 1.0, 200)
    for ax, band in zip(np.atleast_1d(axes), bands):
        for scn, colour in zip(scenes, plt.cm.viridis(np.linspace(0, 0.85, len(scenes)))):
            key = zero_phase_key(profiles, scn.name, band, sza)
            if key is None:
                continue
            rows = profiles[key]
            m = disc_mean(rows)
            sigma_r = 0.5 * scn.favd * scn.r
            ax.errorbar(rows[:, 0], rows[:, 1] / m, yerr=rows[:, 2] / m,
                        color=colour, marker="o", ms=3, lw=1.2,
                        label=r"%s  $\sigma_r$=%.2f" % (scn.name, sigma_r))
        # What the chord average already implies, for the scene whose optical
        # radius sits in the middle of the set -- one curve, or the panel is
        # unreadable.  It is a shape, so it is normalised the same way.
        mid = scenes[len(scenes) // 2]
        band_obj = [b for b in config.BRDF_BANDS if b.name == band][0]
        pu = chord_profile(mid, band_obj, u, sza)
        wgt = 2.0 * u
        ax.plot(u, pu / (np.trapz(wgt * pu, u) / np.trapz(wgt, u)),
                color="#e67e22", ls="-.", lw=1.6,
                label="GORT chord average, %s" % mid.name)
        # (3/2) sqrt(1-u^2) has disc mean 1 under the 2u du weight.
        c, ls = LAW_STYLE["Lambert sphere"]
        ax.plot(u, 1.5 * np.sqrt(np.maximum(0.0, 1.0 - u * u)), color=c, ls=ls,
                lw=1.6, label=r"Lambert sphere, $\frac{3}{2}\sqrt{1-u^2}$")
        c, ls = LAW_STYLE["Lommel-Seeliger / GORT disc"]
        ax.axhline(1.0, color=c, ls=ls, lw=1.6,
                   label="Lommel-Seeliger, and GORT's uniform disc")
        ax.set_xlabel("$u$ = impact parameter / crown radius")
        ax.set_title("%s, phase $\\approx 0$" % band)
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.25)
    np.atleast_1d(axes)[0].set_ylabel("radiance factor / disc mean")
    np.atleast_1d(axes)[-1].legend(fontsize=7, loc="lower left", framealpha=0.9)
    fig.suptitle("Brightness across the projected disc of a single crown, "
                 "black soil, sza %g" % sza, fontsize=11)
    fig.tight_layout()
    path = os.path.join(outdir, "crown-profile-phase0.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_profiles(profiles, scn, sza, outdir):
    """Every sampled phase angle for one scene, coloured by phase."""
    bands = [b.name for b in config.BRDF_BANDS]
    fig, axes = plt.subplots(1, len(bands), figsize=(11, 4.6), sharey=True)

    keys = [k for k in profiles if k[0] == scn.name and k[2] == sza]
    phases = sorted(set(profiles[k][0, 4] for k in keys))
    norm = plt.Normalize(min(phases), max(phases))

    for ax, band in zip(np.atleast_1d(axes), bands):
        for k in sorted(keys, key=lambda k: profiles[k][0, 4]):
            if k[1] != band:
                continue
            rows = profiles[k]
            ph = rows[0, 4]
            ax.plot(rows[:, 0], rows[:, 1] / disc_mean(rows),
                    color=plt.cm.plasma(norm(ph)), lw=1.3)
        ax.plot(np.linspace(0, 1, 200),
                1.5 * np.sqrt(np.maximum(0.0, 1.0 - np.linspace(0, 1, 200) ** 2)),
                color="#111111", ls="-", lw=1.4, label="Lambert sphere")
        ax.axhline(1.0, color="#888888", ls="--", lw=1.4, label="uniform disc")
        ax.set_xlabel("$u$")
        ax.set_title(band)
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.25)
    np.atleast_1d(axes)[0].set_ylabel("radiance factor / disc mean")
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="lower left")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.plasma)
    fig.colorbar(sm, ax=np.atleast_1d(axes).tolist(), label="phase angle (deg)")
    fig.suptitle("%s: profile across the disc at every phase angle sampled "
                 "(sza %g)" % (scn.name, sza), fontsize=11)
    path = os.path.join(outdir, "crown-profile-%s.png" % scn.name)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_phase(phase, scenes, sza, outdir):
    """The disc mean against signed view zenith: the crown's phase function."""
    bands = [b.name for b in config.BRDF_BANDS]
    fig, axes = plt.subplots(1, len(bands), figsize=(11, 4.6))

    for ax, bname in zip(np.atleast_1d(axes), bands):
        band = [b for b in config.BRDF_BANDS if b.name == bname][0]
        for scn, colour in zip(scenes, plt.cm.viridis(np.linspace(0, 0.85, len(scenes)))):
            rows = phase.get((scn.name, bname, sza))
            if rows is None:
                continue
            ax.errorbar(rows[:, 0], rows[:, 2], yerr=rows[:, 3],
                        color=colour, marker="o", ms=3, lw=1.2, ls="none",
                        label="%s  FLIGHT" % scn.name)
            c, kc, kt = gortt_crown(scn, band, sza, rows[:, 0])
            ax.plot(rows[:, 0], c, color=colour, lw=1.2,
                    label="%s  GORT $C$" % scn.name)
            ax.plot(rows[:, 0], c * kc / (kc + kt), color=colour, lw=1.0, ls="--",
                    label="%s  GORT $K_cC/(K_c+K_t)$" % scn.name)
        ax.axvline(sza, color="#888888", lw=0.8, ls=":")
        ax.text(sza, ax.get_ylim()[1], " hotspot", fontsize=7, va="top", color="#666")
        ax.set_xlabel("view zenith (deg), positive = backscatter")
        ax.set_title(bname)
        ax.set_yscale("log")
        ax.grid(alpha=0.25)
    np.atleast_1d(axes)[0].set_ylabel("disc-mean radiance factor")
    np.atleast_1d(axes)[-1].legend(fontsize=5.5, ncol=1)
    fig.suptitle("The crown's own phase function: disc-mean radiance factor, "
                 "no soil, no neighbours (sza %g).  Points are FLIGHT, lines GORT."
                 % sza, fontsize=11)
    fig.tight_layout()
    path = os.path.join(outdir, "crown-phase.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def summarise(profiles, phase, scenes, sza):
    """
    The numbers the figures are there to support, as a table on stdout.

    limb ratio is the outermost u bin over the innermost, at zero phase: 1 for
    a uniform disc, 0.17 for a Lambert sphere with these 20 bins.
    """
    print("\nAt zero phase (Kt ~ 0, so every bit of crown in view is sunlit and "
          "GORT's C is\ndirectly the disc mean it predicts).\n")
    print("%-22s %5s %4s  %7s %7s %6s   %7s %7s %6s"
          % ("scene", "sig_r", "band", "u=.025", "u=.975", "limb", "FLIGHT", "GORT C", "G/F"))
    for scn in scenes:
        for bname in [b.name for b in config.BRDF_BANDS]:
            key = zero_phase_key(profiles, scn.name, bname, sza)
            if key is None:
                continue
            rows = profiles[key]
            m = disc_mean(rows)
            band = [b for b in config.BRDF_BANDS if b.name == bname][0]
            signed = key[3] if key[4] == 0.0 else -key[3]
            c = gortt_crown(scn, band, sza, [signed])[0][0]
            print("%-22s %5.2f %4s  %7.4f %7.4f %6.3f   %7.4f %7.4f %6.2f"
                  % (scn.name, 0.5 * scn.favd * scn.r, bname,
                     rows[0, 1], rows[-1, 1], rows[-1, 1] / rows[0, 1],
                     m, c, c / m if m else float("nan")))
    lam = 1.5 * np.sqrt(1.0 - 0.975 ** 2) / (1.5 * np.sqrt(1.0 - 0.025 ** 2))
    print("\nlimb = outermost u bin / innermost.  Uniform disc (what GORT assumes) "
          "= 1.000,\nLambert sphere = %.3f." % lam)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=config.RESULTS_DIR)
    ap.add_argument("--outdir", default=config.FIGS_DIR)
    ap.add_argument("--sza", type=float, default=config.CROWN_IMAGE_SZAS[0])
    ap.add_argument("--check", action="store_true",
                    help="check the numpy two_stream port against gortt and exit")
    args = ap.parse_args()

    profiles = read_profile(os.path.join(args.results, "crown_profile.csv"))
    phase = read_phase(os.path.join(args.results, "crown_phase.csv"))
    names = sorted(set(k[0] for k in profiles))
    scenes = [config.scene_by_name(n) for n in config.CROWN_IMAGE_SCENES
              if n in names]

    if args.check:
        print("\nchecking the numpy two_stream port against gortt:")
        ok = check_two_stream(scenes, args.sza)
        print("checking that -rt_crown_brdf conserves energy:")
        ok = check_energy(scenes, args.sza) and ok
        sys.exit(0 if ok else 1)

    os.makedirs(args.outdir, exist_ok=True)
    print(plot_phase0(profiles, scenes, args.sza, args.outdir))
    for scn in scenes:
        print(plot_profiles(profiles, scn, args.sza, args.outdir))
    print(plot_phase(phase, scenes, args.sza, args.outdir))
    summarise(profiles, phase, scenes, args.sza)


if __name__ == "__main__":
    main()
