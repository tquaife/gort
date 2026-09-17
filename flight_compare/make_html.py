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


def scene_records(figdir, outdir, suffix=""):
    """
    One record per scene, with whichever of its four figures exist.

    suffix picks a variant's set: plot_brdf.py and plot_lidar.py give the
    baseline the unsuffixed names and every other variant a -<name> set beside
    them, so a variant page is the same builder pointed at different files.
    """
    out = []
    for s in config.scenes():
        figs = {}
        for band in config.BRDF_BANDS:
            p = os.path.join(figdir, "%s-polar-%s%s.png" % (s.name, band.name, suffix))
            if os.path.exists(p):
                figs["polar_" + band.name] = os.path.relpath(p, outdir)
        for key, fname in (("pp", "%s-pp%s.png" % (s.name, suffix)),
                           ("lidar", "%s-lidar%s.png" % (s.name, suffix))):
            p = os.path.join(figdir, fname)
            if os.path.exists(p):
                figs[key] = os.path.relpath(p, outdir)
        # No -rt_* option reaches the lidar path, so a variant has no lidar
        # figure of its own; fall back to the baseline one and say so in the
        # caption rather than leaving a gap where a panel should be.
        if "lidar" not in figs:
            p = os.path.join(figdir, "%s-lidar.png" % s.name)
            if os.path.exists(p):
                figs["lidar"] = os.path.relpath(p, outdir)
                figs["lidar_is_baseline"] = True
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
.variant { background: #fff8e6; border-left: 3px solid #c7911f; padding: 10px 14px;
        margin: 12px 0 4px; max-width: 90ch; font-size: 13px; }
.variant p { margin: 4px 0; }
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


def overview_block(figdir, outdir, suffix=""):
    """
    The scatter and the per-scene bar chart.

    suffix picks the variant's own pair, the same way scene_records does. It is
    a parameter and not a constant because it was a constant once: every variant
    page showed the baseline overview while its own overview-<name>.png sat
    unreferenced beside it, which is exactly the kind of wrong that looks right.
    """
    p = os.path.join(figdir, "overview%s.png" % suffix)
    if not os.path.exists(p):
        p = os.path.join(figdir, "overview.png")
    if not os.path.exists(p):
        return ""
    return ("<h2>Overview</h2>"
            "<figure style='margin:0'>"
            "<figcaption>Every view direction in every scene, %s against "
            "FLIGHT, and the relative disagreement scene by scene. The 1:1 line "
            "is drawn; points above it are gortt reading high.</figcaption>"
            % ("gortt" if not suffix else "gortt with <code>%s</code>"
               % html.escape(suffix.lstrip("-"))) +
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



# The block at the top of a variant page. Kept out of the main template so that
# the baseline page -- the one published at tquaife.github.io -- renders exactly
# as it did before any of this.
VARIANT_BANNER = """
  <div class="variant">
    <p><strong>Experimental variant: <code>%(name)s</code></strong>
    &mdash; %(blurb)s.</p>
    <p class="sub">Every gortt curve labelled <code>%(name)s</code> was produced
    with <code>%(flags)s</code>; every curve labelled <code>baseline</code> is
    Ni et al. (1999) as published, which is what <code>gortt</code> does with no
    <code>-rt_*</code> option. <strong>FLIGHT was not re-run</strong>: it is the
    same Monte Carlo on both pages, so any change here is entirely gortt's.
    The polar panels show one formulation at a time &mdash; the third column,
    the difference against FLIGHT, is the one to read against the baseline
    page. The principal-plane and lidar panels show both at once.</p>
    <p class="sub"><a href="index.html">&larr; all variants</a> &nbsp;
    <a href="../index.html">baseline comparison</a></p>
  </div>
"""


def variant_index(outdir, written):
    """A contents page for variants/, linking each page with what it changes."""
    rows = "".join(
        "<tr><td><a href='%s'>%s</a></td><td><code>%s</code></td><td>%s</td></tr>"
        % (html.escape(fname), html.escape(name),
           html.escape(" ".join(flags)), html.escape(blurb))
        for name, flags, blurb, fname in written)
    page = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GORT vs FLIGHT &mdash; variants</title>
<style>%(css)s</style></head>
<body>
<header>
  <h1>Experimental scattering variants</h1>
  <p class="sub">Ni et al. (1999) section V states that the model has been
  validated over dense forests, and defers three things that matter as a canopy
  opens up: the bidirectional effect of the background, the exchange of
  radiation between background and canopy, and scattering between neighbouring
  crowns. The <code>-rt_*</code> options switch them on one at a time, and each
  page below is that option measured against the same FLIGHT Monte Carlo.</p>
  <p class="sub"><strong>None of them is on by default.</strong> With no
  <code>-rt_*</code> option gortt is the paper as published, and all 82 original
  regression goldens are byte-identical &mdash; that is the guard. FLIGHT was
  not re-run for any of this.</p>
  <p class="sub"><a href="../index.html">&larr; baseline comparison</a></p>
</header>
<main>
  <table><thead><tr><th>variant</th><th>options</th><th>what it changes</th>
  </tr></thead><tbody>%(rows)s</tbody></table>

  <div class="finding">
  <h2>What they are worth</h2>
  <p>Bias and RMSD against FLIGHT as a percentage of FLIGHT's own mean, over all
  12 scenes &times; 3 solar zeniths &times; 64 directions. The red band is the
  guard: it already agrees, and a formulation that improves the near infrared at
  its expense is not an improvement.</p>
  <table><thead><tr><th>variant</th><th>red bias</th><th>red RMSD</th>
  <th>NIR bias</th><th>NIR RMSD</th><th>NIR bias, sparse only</th>
  <th>NIR RMSD, sparse only</th></tr></thead><tbody>
  <tr><td>baseline</td><td>+0.2%%</td><td>14.9%%</td><td>&minus;7.4%%</td><td>18.6%%</td><td>&minus;9.5%%</td><td>18.4%%</td></tr>
  <tr><td><strong>crown-path</strong></td><td>+1.0%%</td><td>15.2%%</td><td><strong>+1.7%%</strong></td><td><strong>16.7%%</strong></td><td><strong>+2.5%%</strong></td><td>13.9%%</td></tr>
  <tr><td>squeeze-path</td><td>+1.5%%</td><td>15.9%%</td><td>+6.6%%</td><td>19.3%%</td><td>+1.6%%</td><td><strong>13.3%%</strong></td></tr>
  <tr><td>recollision</td><td>+1.1%%</td><td>16.5%%</td><td>+6.4%%</td><td>20.3%%</td><td>+4.2%%</td><td>14.8%%</td></tr>
  <tr><td>bg-coupling</td><td>+0.5%%</td><td>14.8%%</td><td>&minus;4.0%%</td><td>18.0%%</td><td>&minus;4.3%%</td><td>17.4%%</td></tr>
  <tr><td>crown-crown</td><td>+0.5%%</td><td>14.9%%</td><td>+13.7%%</td><td>37.4%%</td><td>&minus;0.9%%</td><td>22.1%%</td></tr>
  <tr><td>all</td><td>+1.7%%</td><td>16.6%%</td><td>+43.4%%</td><td>63.5%%</td><td>+25.5%%</td><td>31.9%%</td></tr>
  </tbody></table>

  <ul>
  <li><strong><code>crown-path</code> is the one that works.</strong> It takes
  the near-infrared bias from &minus;7.4%% to +1.7%% and lowers the RMSD as well,
  which a pure rescaling would not do &mdash; it is fixing a shape error, not
  adding a constant. In the sparse scenes, where the published approximation is
  furthest outside the range section&nbsp;V claims, it goes from &minus;9.5%% to
  +2.5%% and the RMSD falls from 18.4%% to 13.9%%. The red band moves by
  0.8&nbsp;points, i.e. not at all.</li>

  <li><strong>It is also near-inert where the paper was validated.</strong> On
  the paper's own SOBS canopy it moves the Fig.&nbsp;7b sunlit-crown signature
  from 0.470 to 0.474. That is the property to want: it should change the sparse
  canopies and leave the dense ones alone.</li>

  <li><strong><code>squeeze-path</code> is the better answer for sparse canopies
  and the worse one overall</strong> (sparse RMSD 13.3%% against 13.9%%, but
  +6.6%% bias overall against +1.7%%). The two bracket the truth rather than
  settling it, which is honest: one spreads the foliage over the crowns' own
  footprint, the other over the part of the hemisphere that is not open, and
  nothing here says which a crown really is.</li>

  <li><strong>The recollision probability does not earn its place, and the
  reason is measurable.</strong> <code>run_flight_single.py --black-soil</code>
  measures p from FLIGHT directly: with a black soil the canopy alone satisfies
  BRF(all&nbsp;orders)/BRF(single) = 1/(1&minus;&omega;p). Against gortt's
  Wigner estimate p = &sigma;l/(1+&sigma;l) on the crown's mean chord:
  0.435 against 0.37&ndash;0.46 for the thin crowns (good), but 0.909 against
  0.53&ndash;0.57 for the thick ones. The first-flight escape probability
  assumes the scattering is uniformly distributed through the crown; in an
  optically thick crown it happens near the lit surface, where escape is easier,
  so p comes out too high. Correct in the thin limit, wrong where it matters.</li>

  <li><strong><code>bg-coupling</code> is small, safe and in the right
  direction</strong> &mdash; &minus;7.4%% to &minus;4.0%% with the RMSD very
  slightly better in both bands. It is the least arguable of the five: as it
  stands the sunlit background cannot see the canopy at all.</li>

  <li><strong><code>crown-crown</code> is not supported.</strong> It improves the
  sparse bias and wrecks the RMSD (18.6%% to 37.4%% overall), i.e. it adds
  brightness in the wrong places. The paper calls it "usually small and very
  complicated to model"; on this evidence the first-order form here is not the
  way to put it back.</li>

  <li><strong>They do not compose.</strong> <code>all</code> overshoots by
  +43.4%%. <code>crown-path</code> and <code>recollision</code> are two answers
  to the same question, and adding <code>crown-crown</code> on top counts the
  same photons a third time.</li>
  </ul>

  <p class="sub"><strong>The lidar is untouched by all of this</strong> and is
  not repeated per variant: <code>gortt_lidar</code> takes its backscatter
  coefficients straight from the leaf transmittance and the soil albedo and
  never calls <code>gortt_rt_init</code>. Checked, not assumed &mdash;
  <code>gortt -lidar -rt_all</code> is byte-identical to
  <code>gortt -lidar</code>.</p>
  </div>
</main>
</body></html>
""" % {"css": CSS, "rows": rows}
    out = os.path.join(outdir, "index.html")
    with open(out, "w") as fp:
        fp.write(page)
    return out


def build(outdir, figdir, variant=None, name="index.html", extra_head=""):
    """
    Assemble one page.

    variant is the gortt formulation this page is about; None means the
    baseline, which is the main page and keeps every filename it had. For a
    variant the statistics tables are filtered to the baseline and that variant
    so the two can be read against each other, and the figures come from the
    suffixed set.
    """
    baseline = config.BASELINE
    suffix = "" if variant in (None, baseline) else "-" + variant
    keep = [baseline] if variant in (None, baseline) else [baseline, variant]

    records = scene_records(figdir, outdir, suffix)
    brdf_stats = [r for r in read_csv(os.path.join(config.RESULTS_DIR, "brdf_stats.csv"))
                  if r.get("variant", baseline) in keep]
    lidar_stats = [r for r in read_csv(os.path.join(config.RESULTS_DIR, "lidar_stats.csv"))
                   if r.get("variant", baseline) in keep]

    blocks = []
    for r in records:
        figs = []
        for band in config.BRDF_BANDS:
            key = "polar_" + band.name
            # baseline only -- see the note in plot_brdf.main
            if key in r["figs"]:
                caption = ("Polar BRDF, %s (%d nm) &mdash; gortt, FLIGHT, and "
                           "their difference, at three solar zeniths. Hotspot "
                           "at the top." % (band.name.upper(), band.wl))
                if suffix:
                    # the variant figure is the difference pair, not the
                    # absolute fields -- see plot_brdf.plot_polar_diff
                    caption = ("Polar BRDF, %s (%d nm) &mdash; how the "
                               "disagreement with FLIGHT changes, baseline "
                               "beside <code>%s</code>, at three solar zeniths "
                               "on one shared scale. Hotspot at the top; blue "
                               "is gortt reading low. The absolute fields are "
                               "on the <a href='../index.html'>baseline "
                               "page</a> and FLIGHT is the same Monte Carlo "
                               "either way."
                               % (band.name.upper(), band.wl,
                                  html.escape(variant or "")))
                figs.append(
                    "<figure data-view='polar'><figcaption>%s</figcaption>"
                    "<img loading='lazy' src='%s' alt='polar BRDF %s %s'></figure>"
                    % (caption, html.escape(r["figs"][key]), r["name"], band.name))
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
                "arrivals are recorded as negative heights.%s</figcaption>"
                "<img loading='lazy' src='%s' alt='lidar waveform %s'></figure>"
                % (config.LIDAR_BAND.wl,
                   (" <strong>Unchanged by this variant</strong> &mdash; "
                    "<code>gortt_lidar</code> takes its backscatter "
                    "coefficients straight from the leaf transmittance and the "
                    "soil albedo and never calls <code>gortt_rt_init</code>, so "
                    "no <code>-rt_*</code> option can reach it. This is the "
                    "baseline figure."
                    if r["figs"].get("lidar_is_baseline") else ""),
                   html.escape(r["figs"]["lidar"]), r["name"]))
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
<title>%(title)s</title>
<style>%(css)s</style></head>
<body>
<header>
  <h1>%(title)s</h1>
  %(banner)s
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
        "overview": overview_block(figdir, outdir, suffix),
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
        "title": ("GORT vs FLIGHT" if variant in (None, baseline)
                  else "GORT vs FLIGHT &mdash; %s" % html.escape(variant)),
        "banner": ("" if variant in (None, baseline) else VARIANT_BANNER % {
            "name": html.escape(variant),
            "flags": html.escape(" ".join(config.variant_args(variant))),
            "blurb": html.escape(config.variant_blurb(variant)),
        }),
    }

    out = os.path.join(outdir, name)
    with open(out, "w") as fp:
        fp.write(page)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=config.HERE)
    ap.add_argument("--figs", default=config.FIGS_DIR)
    ap.add_argument("--variants", action="store_true",
                    help="also write one page per gortt variant into "
                         "variants/, plus an index linking them. The main "
                         "index.html is rebuilt either way and is always the "
                         "baseline -- it is the published page.")
    args = ap.parse_args()
    sys.stderr.write("make_html: wrote %s\n" % build(args.outdir, args.figs))

    if args.variants:
        vdir = os.path.join(args.outdir, "variants")
        os.makedirs(vdir, exist_ok=True)
        written = []
        for name, flags, blurb in config.GORTT_VARIANTS:
            if name == config.BASELINE:
                continue
            page = build(vdir, args.figs, variant=name, name="%s.html" % name)
            written.append((name, flags, blurb, os.path.basename(page)))
            sys.stderr.write("make_html: wrote %s\n" % page)
        idx = variant_index(vdir, written)
        sys.stderr.write("make_html: wrote %s\n" % idx)


if __name__ == "__main__":
    main()
