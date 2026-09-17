#!/usr/bin/env python3
"""
Build index.html: a static, file://-openable page for browsing the comparison.

No external assets and no CDN -- the CSS and JS are inline and the images are
referenced relatively out of figs/, so the page works from a checkout with no
network and no server.  The scene table, the agreement statistics and the
figure list are all baked in at build time from the CSVs in results/, so the
page is a snapshot of a particular run rather than something that reads the
data live.
"""

import argparse
import csv
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fp:
        return list(csv.DictReader(fp))


def scene_records(figdir, outdir):
    out = []
    for s in config.scenes():
        figs = {}
        for band in config.BRDF_BANDS:
            p = os.path.join(figdir, "%s-polar-%s.png" % (s.name, band.name))
            if os.path.exists(p):
                figs["polar_" + band.name] = os.path.relpath(p, outdir)
        for key, fname in (("pp", "%s-pp.png" % s.name),
                           ("lidar", "%s-lidar.png" % s.name)):
            p = os.path.join(figdir, fname)
            if os.path.exists(p):
                figs[key] = os.path.relpath(p, outdir)
        out.append({
            "name": s.name,
            "cover": s.cover_name, "lai": s.lai_name, "shape": s.shape_name,
            "r": s.r, "b": s.b, "br": s.br, "lam": s.lam, "favd": s.favd,
            "pcc": s.pcc, "nadir_cover": s.nadir_cover,
            "h1": s.h1, "h2": s.h2, "z1": s.z1, "z2": s.z2,
            "target_lai": s.target_lai, "half_width": s.half_width,
            "ncrowns": s.ncrowns, "figs": figs,
        })
    return out


CSS = """
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { margin: 0; font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI",
       Roboto, Helvetica, Arial, sans-serif; color: #1a1a1a; background: #fafaf8; }
header { padding: 24px 28px 16px; border-bottom: 1px solid #e2e2dd; background: #fff; }
h1 { margin: 0 0 6px; font-size: 20px; font-weight: 600; }
h2 { font-size: 15px; font-weight: 600; margin: 28px 0 10px; }
.sub { color: #5a5a55; max-width: 80ch; margin: 0; }
main { padding: 20px 28px 60px; }
.controls { display: flex; flex-wrap: wrap; gap: 18px; padding: 14px 28px;
            background: #fff; border-bottom: 1px solid #e2e2dd;
            position: sticky; top: 0; z-index: 5; }
.group { display: flex; align-items: center; gap: 6px; }
.group .label { font-size: 11px; text-transform: uppercase; letter-spacing: .06em;
                color: #77776f; margin-right: 2px; }
button { font: inherit; font-size: 13px; padding: 3px 11px; border-radius: 999px;
         border: 1px solid #d3d3cc; background: #fff; color: #333; cursor: pointer; }
button.on { background: #1f6feb; border-color: #1f6feb; color: #fff; }
table { border-collapse: collapse; font-size: 12.5px; background: #fff; }
th, td { padding: 4px 9px; border-bottom: 1px solid #ececE6; text-align: right;
         white-space: nowrap; }
th { font-weight: 600; text-align: right; color: #55554f; border-bottom: 1px solid #d3d3cc; }
td.name, th.name { text-align: left; font-family: ui-monospace, Menlo, Consolas, monospace; }
tr.flag td { background: #fff8e6; }
.scene { background: #fff; border: 1px solid #e2e2dd; border-radius: 8px;
         margin: 0 0 22px; overflow: hidden; }
.scene > .head { padding: 10px 16px; border-bottom: 1px solid #ececE6;
                 display: flex; align-items: baseline; gap: 14px; flex-wrap: wrap; }
.scene > .head .n { font-family: ui-monospace, Menlo, Consolas, monospace;
                    font-weight: 600; font-size: 14px; }
.scene > .head .meta { color: #6a6a63; font-size: 12.5px; }
.figs { padding: 8px 16px 16px; }
.figs figure { margin: 12px 0 0; }
.figs figcaption { font-size: 12px; color: #6a6a63; margin-bottom: 4px; }
img { max-width: 100%; height: auto; display: block; border: 1px solid #ececE6;
      border-radius: 4px; }
.note { background: #fff; border-left: 3px solid #1f6feb; padding: 10px 14px;
        margin: 14px 0; max-width: 90ch; font-size: 13px; }
.warn { border-left-color: #d98014; background: #fffaf2; }
.hidden { display: none; }
.finding { background: #fff; border: 1px solid #e2e2dd; border-radius: 8px;
           padding: 4px 20px 12px; max-width: 95ch; }
.finding ul { padding-left: 20px; }
.finding li { margin: 7px 0; }
.scroll { overflow-x: auto; max-width: 100%; }
code { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12.5px;
       background: #f1f1ec; padding: 1px 4px; border-radius: 3px; }
"""

JS = """
const state = {cover: null, lai: null, shape: null, view: 'all'};
function apply() {
  document.querySelectorAll('.scene').forEach(el => {
    const ok = (!state.cover || el.dataset.cover === state.cover)
            && (!state.lai   || el.dataset.lai   === state.lai)
            && (!state.shape || el.dataset.shape === state.shape);
    el.classList.toggle('hidden', !ok);
  });
  document.querySelectorAll('figure[data-view]').forEach(el => {
    el.classList.toggle('hidden', state.view !== 'all' && el.dataset.view !== state.view);
  });
  document.querySelectorAll('tr[data-scene]').forEach(el => {
    const ok = (!state.cover || el.dataset.cover === state.cover)
            && (!state.lai   || el.dataset.lai   === state.lai)
            && (!state.shape || el.dataset.shape === state.shape);
    el.classList.toggle('hidden', !ok);
  });
}
document.addEventListener('click', ev => {
  const b = ev.target.closest('button[data-key]');
  if (!b) return;
  const key = b.dataset.key, val = b.dataset.val || null;
  state[key] = (state[key] === val) ? null : val;
  if (key === 'view' && state[key] === null) state[key] = 'all';
  document.querySelectorAll(`button[data-key="${key}"]`).forEach(o =>
    o.classList.toggle('on', o.dataset.val === (state[key] || '')));
  apply();
});
apply();
"""


def buttons(key, values, current=""):
    out = []
    for label, val in values:
        on = " on" if val == current else ""
        out.append('<button class="%s" data-key="%s" data-val="%s">%s</button>'
                   % (on.strip(), key, html.escape(val), html.escape(label)))
    return "".join(out)


def overview_block(figdir, outdir):
    p = os.path.join(figdir, "overview.png")
    if not os.path.exists(p):
        return ""
    return ("<h2>Overview</h2>"
            "<figure style='margin:0'>"
            "<figcaption>Every view direction in every scene, gortt against "
            "FLIGHT, and the relative disagreement scene by scene. The 1:1 line "
            "is drawn; points above it are gortt reading high.</figcaption>"
            "<img src='%s' alt='overview'></figure>"
            % html.escape(os.path.relpath(p, outdir)))


def scene_table(records):
    head = ("<tr><th class='name'>scene</th><th>r</th><th>b</th><th>b/r</th>"
            "<th>&lambda; (m<sup>-2</sup>)</th><th>favd (m<sup>-1</sup>)</th>"
            "<th>PCC</th><th>cover</th><th>h1</th><th>h2</th><th>z1</th><th>z2</th>"
            "<th>LAI</th><th>scene (m)</th><th>crowns</th></tr>")
    rows = []
    for r in records:
        flag = " class='flag'" if r["favd"] > 8.0 else ""
        rows.append(
            "<tr data-scene='%s' data-cover='%s' data-lai='%s' data-shape='%s'%s>"
            "<td class='name'>%s</td><td>%.1f</td><td>%.1f</td><td>%.1f</td>"
            "<td>%.4f</td><td>%.2f</td><td>%.2f</td><td>%.3f</td>"
            "<td>%.1f</td><td>%.1f</td><td>%.1f</td><td>%.1f</td>"
            "<td>%.1f</td><td>%.0f</td><td>%d</td></tr>"
            % (r["name"], r["cover"], r["lai"], r["shape"], flag, r["name"],
               r["r"], r["b"], r["br"], r["lam"], r["favd"], r["pcc"],
               r["nadir_cover"], r["h1"], r["h2"], r["z1"], r["z2"],
               r["target_lai"], 2 * r["half_width"], r["ncrowns"]))
    return "<div class='scroll'><table>%s%s</table></div>" % (head, "".join(rows))


def stats_table(rows, columns, lookup):
    if not rows:
        return "<p class='sub'>Not available &mdash; the run has not been done yet.</p>"
    by_scene = {s.name: s for s in config.scenes()}
    head = "<tr>" + "".join("<th class='%s'>%s</th>"
                            % ("name" if c == columns[0][0] else "", c[1])
                            for c in columns) + "</tr>"
    body = []
    for r in rows:
        scn = by_scene.get(r.get("scene"))
        if scn is None:
            continue
        cells = []
        for key, _ in columns:
            v = r.get(key, "")
            cells.append("<td class='%s'>%s</td>"
                         % ("name" if key == columns[0][0] else "", html.escape(str(v))))
        body.append("<tr data-scene='%s' data-cover='%s' data-lai='%s' data-shape='%s'>%s</tr>"
                    % (scn.name, scn.cover_name, scn.lai_name, scn.shape_name,
                       "".join(cells)))
    return "<div class='scroll'><table>%s%s</table></div>" % (head, "".join(body))


def build(outdir, figdir):
    records = scene_records(figdir, outdir)
    brdf_stats = read_csv(os.path.join(config.RESULTS_DIR, "brdf_stats.csv"))
    lidar_stats = read_csv(os.path.join(config.RESULTS_DIR, "lidar_stats.csv"))

    blocks = []
    for r in records:
        figs = []
        for band in config.BRDF_BANDS:
            key = "polar_" + band.name
            if key in r["figs"]:
                figs.append(
                    "<figure data-view='polar'><figcaption>Polar BRDF, %s "
                    "(%d nm) &mdash; gortt, FLIGHT, and their difference, at "
                    "three solar zeniths. Hotspot at the top.</figcaption>"
                    "<img loading='lazy' src='%s' alt='polar BRDF %s %s'></figure>"
                    % (band.name.upper(), band.wl, html.escape(r["figs"][key]),
                       r["name"], band.name))
        if "pp" in r["figs"]:
            figs.append(
                "<figure data-view='pp'><figcaption>Principal plane. Positive "
                "view zenith is backscatter; FLIGHT error bars are the standard "
                "error over scene realisations.</figcaption>"
                "<img loading='lazy' src='%s' alt='principal plane %s'></figure>"
                % (html.escape(r["figs"]["pp"]), r["name"]))
        if "lidar" in r["figs"]:
            figs.append(
                "<figure data-view='lidar'><figcaption>Nadir lidar waveform at "
                "%d nm, normalised to unit total energy, plus the energy "
                "partition. <em>Orders 1&ndash;8</em> is FLIGHT's cumulative "
                "waveform over scattering orders 1 to 8, single scattering "
                "included; its effect is the tail below the ground, where late "
                "arrivals are recorded as negative heights.</figcaption>"
                "<img loading='lazy' src='%s' alt='lidar waveform %s'></figure>"
                % (config.LIDAR_BAND.wl, html.escape(r["figs"]["lidar"]), r["name"]))
        if not figs:
            figs.append("<p class='sub'>No figures for this scene yet.</p>")

        blocks.append(
            "<section class='scene' data-cover='%s' data-lai='%s' data-shape='%s'>"
            "<div class='head'><span class='n'>%s</span>"
            "<span class='meta'>%s cover &middot; %s LAI (target %.1f) &middot; %s "
            "crowns, b/r = %.1f &middot; favd %.2f m<sup>-1</sup> &middot; "
            "%.0f &times; %.0f m scene</span></div>"
            "<div class='figs'>%s</div></section>"
            % (r["cover"], r["lai"], r["shape"], r["name"], r["cover"], r["lai"],
               r["target_lai"], r["shape"], r["br"], r["favd"],
               2 * r["half_width"], 2 * r["half_width"], "".join(figs)))

    page = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GORT vs FLIGHT</title>
<style>%(css)s</style></head>
<body>
<header>
  <h1>GORT vs FLIGHT</h1>
  <p class="sub">The analytic GORT model (<code>gortt</code>) against the FLIGHT
  Monte Carlo photon-tracing model, over matched canopies: the same spheroidal
  crowns at the same positions, the same in-crown foliage density, the same
  spherical leaf angle distribution, the same leaf and soil optical properties,
  and direct-beam illumination only. Bidirectional reflectance at %(bands)s, and
  the nadir lidar waveform at %(lidarwl)d&nbsp;nm.</p>
</header>

<div class="controls">
  <div class="group"><span class="label">Cover</span>%(cover)s</div>
  <div class="group"><span class="label">LAI</span>%(lai)s</div>
  <div class="group"><span class="label">Crown</span>%(shape)s</div>
  <div class="group"><span class="label">Show</span>%(view)s</div>
</div>

<main>
  <div class="note">
    <strong>How the two models were matched.</strong> Crowns are spheroids of
    semi-axes (r,&nbsp;r,&nbsp;b) with centres uniform in [h1,&nbsp;h2] and
    uniform in plan at density &lambda;. FLIGHT is given them as an explicit
    crown list rather than asked to generate its own, and <code>TOTAL_LAI</code>
    is set negative so it takes the in-crown foliage area density from the
    per-crown column &mdash; the same <code>favd</code> that
    <code>gortt&nbsp;-favd</code> gets. PROSPECT and the Price soil model are
    switched off in gortt (<code>-refl_leaf</code>, <code>-trans_leaf</code>,
    <code>-alb_soil</code>) so both models use identical constants. Before any
    of this was run, the nadir gap probability was checked between the two
    models and agrees to better than 1%%.
  </div>
  <div class="note warn">
    <strong>Deliberately not compared:</strong> the areal proportions of the
    scene components (gortt's Kc/Kg/Kt/Kz against FLIGHT's per-component
    fractions). <strong>Absolute lidar magnitudes</strong> are not compared
    either: gortt builds its volume backscatter coefficient from the leaf
    <em>transmittance</em>, which is not the same quantity FLIGHT returns, so
    the waveforms are compared as normalised shapes with the energy partition
    shown separately.
  </div>

  %(overview)s

  <h2>What it shows</h2>
  <div class="finding">
    <p><strong>The overall bias is small; it is the shape of the BRDF that
    differs.</strong> Over all scenes, solar zeniths and view directions, gortt
    sits at +0.2%% bias and 15%% RMSD in the red, and &minus;7.4%% bias and 19%%
    RMSD in the NIR, relative to the mean FLIGHT BRF.</p>
    <ul>
      <li><strong>Agreement degrades with view zenith.</strong> Red RMSD runs
      6.6%% within 20&deg; of nadir, 16%% at 40&ndash;60&deg;, and 27%% beyond
      60&deg;, where gortt reads high. Its polar plots carry a much stronger
      bright limb than FLIGHT's, which are comparatively flat.</li>
      <li><strong>Crown shape is the strongest control.</strong> Red RMSD is
      7.6%% for oblate crowns, 11%% for spheres, 24%% for prolate. The worst
      scene is <code>dense-low-prolate</code>, at 34&ndash;52%% in the red, with
      the correlation across directions collapsing to &minus;0.25 at a solar
      zenith of 15&deg;.</li>
      <li><strong>Red splits by cover</strong> (5.1%% sparse against 24.5%%
      dense); <strong>NIR splits by leaf area</strong> (&minus;16%% bias at low
      LAI, unbiased at high LAI), which points at the multiple-scattering
      treatment rather than the geometry.</li>
      <li><strong>The forward-scatter half is worse in the NIR</strong>
      (&minus;10.8%% bias, 24%% RMSD, against &minus;4.6%% and 16%% in
      backscatter). The principal-plane panels show gortt's forward limb turning
      back up into a symmetric U where FLIGHT falls away monotonically &mdash;
      the signature of <code>OPEN-QUESTIONS.md</code> item&nbsp;9,
      &beta; in <code>gortt_kc</code> having no view-direction dependence,
      reached here independently of gortt's own kernel Monte Carlo.</li>
      <li><strong>The lidar waveform shape is the best-agreeing thing here.</strong>
      With the ground return removed, the canopy profiles differ by an RMSD of
      0.009&ndash;0.020 on a profile integrating to one, and mean return heights
      agree to 0.5&nbsp;m in the high-LAI scenes. gortt's mean return height is
      biased low, by 1.2&ndash;1.4&nbsp;m in the low-LAI prolate canopies.</li>
      <li><strong>The lidar energy split differs by a factor of two, and it is
      not the canopy structure.</strong> The nadir gap probability was validated
      to within 0.8%% beforehand, so the difference can only be in the
      backscatter coefficients &mdash; and it is: the implied &rho;<sub>g</sub>/
      &rho;<sub>v</sub> is exactly 0.426 in gortt (every scene) against 0.79 in
      FLIGHT. Numerically the factor of 0.5 on the soil reflectance is all of
      it. It does not affect the waveform shape.</li>
      <li><strong>Multiple scattering shows up below the ground, not in the
      canopy.</strong> The <em>orders 1&ndash;8</em> curve is FLIGHT's
      cumulative waveform over scattering orders 1 to 8, single scattering
      included. A multiply scattered photon takes a longer route back from the
      same height, and a longer route is recorded as a lower height, so the
      ground peak grows a tail <em>downwards</em>: the standard deviation of the
      sub-1&nbsp;m return goes from 0.18&nbsp;m to 0.8&ndash;1.8&nbsp;m, and
      1.3&ndash;9.2%% of the return arrives below &minus;1&nbsp;m, running to the
      floor of the recording window at &minus;10&nbsp;m. Largest in the dense,
      low-LAI canopies. The canopy shape barely moves, and gortt has no
      counterpart &mdash; Ni-Meister et al. (2001) is single scattering.</li>
      <li><strong>FLIGHT's ground return sits at &minus;0.198&nbsp;m, not
      zero.</strong> About &minus;0.15&nbsp;m of that is real range walk: the
      source is a point 1&nbsp;km up, so a photon entering the 25&nbsp;m
      footprint off-axis travels further, by
      E[&rho;<sup>2</sup>]/(2&middot;range) = 0.156&nbsp;m on average. Shrinking
      the footprint to 1&nbsp;m removes it. The remaining &minus;0.05&nbsp;m is
      a fixed offset in FLIGHT's range origin, invariant to canopy top,
      footprint, range, pulse width and bin count. Both are inside one
      0.1&nbsp;m bin and neither affects a shape comparison.</li>
    </ul>
    <p><strong>FLIGHT is a reference, not truth.</strong> It has its own
    approximations &mdash; leaves as a turbid medium with a disc hotspot
    correction, an isotropic diffuse sky (<code>g_sky_model</code> is forced to
    0 on both branches of <code>flight.c:7409</code>, so the anisotropic lookup
    table is never read), and Monte Carlo noise of a few tenths of a percent on
    these averages. The scenes are idealised in both models: no trunks, no
    understory, Lambertian soil, direct beam only.</p>
  </div>

  <h2>Scenes</h2>
  <p class="sub">Highlighted rows have an in-crown foliage density above
  8&nbsp;m<sup>-1</sup>: that is what a high leaf area inside a shallow, sparse
  canopy implies, and it is a deliberate stress case rather than a typical
  stand.</p>
  %(scenetable)s

  <h2>BRDF agreement</h2>
  <p class="sub">gortt minus FLIGHT over the directions FLIGHT was actually run
  at, averaged over scene realisations.</p>
  %(brdfstats)s

  <h2>Lidar agreement</h2>
  <p class="sub">Shape statistics over the canopy only, with the ground return
  excluded and each profile renormalised.</p>
  %(lidarstats)s

  <h2>Results by scene</h2>
  %(blocks)s
</main>
<script>%(js)s</script>
</body></html>
""" % {
        "css": CSS,
        "js": JS,
        "bands": " and ".join("%d&nbsp;nm" % b.wl for b in config.BRDF_BANDS),
        "lidarwl": config.LIDAR_BAND.wl,
        "cover": buttons("cover", [(c, c) for c, _ in config.COVERS]),
        "lai": buttons("lai", [(n, n) for n, _ in config.LAIS]),
        "shape": buttons("shape", [(n, n) for n, _ in config.SHAPES]),
        "view": buttons("view", [("all", "all"), ("polar", "polar"),
                                 ("principal plane", "pp"), ("lidar", "lidar")], "all"),
        "overview": overview_block(figdir, outdir),
        "scenetable": scene_table(records),
        "brdfstats": stats_table(
            brdf_stats,
            [("scene", "scene"), ("band", "band"), ("sza", "sza"),
             ("mean_gortt", "mean gortt"), ("mean_flight", "mean FLIGHT"),
             ("bias", "bias"), ("rmsd", "RMSD"), ("rel_rmsd", "RMSD/mean"),
             ("r", "r")], None),
        "lidarstats": stats_table(
            lidar_stats,
            [("scene", "scene"), ("ground_frac_gortt", "ground frac gortt"),
             ("ground_frac_flight", "ground frac FLIGHT"),
             ("mean_height_gortt", "mean height gortt (m)"),
             ("mean_height_flight", "mean height FLIGHT (m)"),
             ("peak_z_gortt", "peak z gortt (m)"),
             ("peak_z_flight", "peak z FLIGHT (m)"),
             ("canopy_shape_rmsd", "shape RMSD")], None),
        "blocks": "".join(blocks),
    }

    out = os.path.join(outdir, "index.html")
    with open(out, "w") as fp:
        fp.write(page)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=config.HERE)
    ap.add_argument("--figs", default=config.FIGS_DIR)
    args = ap.parse_args()
    sys.stderr.write("make_html: wrote %s\n" % build(args.outdir, args.figs))


if __name__ == "__main__":
    main()
