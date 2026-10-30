"""Interactive HTML map (vis-network) with mechanism toggles and per-target path views."""
import json

VIS_JS = "https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.9/standalone/umd/vis-network.min.js"


def _clean(v):
    return v if v == v else ""  # NaN -> ""


def build_data(G, ranking, overview_nodes, views, index=None, best=None, title=""):
    """Everything the page needs, as plain JSON-able data.

    views: {target: rows from mechanisms.target_paths} (rows may lack mechanism data).
    """
    index = index or {}
    info = ranking.set_index("protein")
    shown = set(overview_nodes)
    vdata = {}
    for t, rows in views.items():
        nodes, pairs = {t}, set()
        for r in rows:
            chain = [t, r["via"], r["input"]] if r["via"] else [t, r["input"]]
            nodes.update(chain)
            pairs.update(tuple(sorted(p)) for p in zip(chain, chain[1:]))
        shown |= nodes
        vdata[t] = {"nodes": sorted(nodes), "pairs": [list(p) for p in sorted(pairs)], "rows": [
            {k: r[k] for k in ("input", "hops", "via", "effect_on_input", "how",
                               "effect_on_target", "how_reverse")} for r in rows]}

    nodes = []
    for n in sorted(shown):
        r = info.loc[n]
        nodes.append({"id": n, "input": bool(G.nodes[n]["is_input"]), "overview": n in overview_nodes,
                      "score": float(r.score), "direct": int(r.n_input_partners),
                      "twohop": int(r.n_two_hop), "partners": _clean(r.input_partners),
                      "twohop_text": _clean(r.two_hop_partners)})
    pairs = [{"a": a, "b": b, "sources": sorted(d["sources"])}
             for a, b, d in G.subgraph(shown).edges(data=True)]
    mechs = [{"a": a, "b": b, "sign": m["sign"], "effect": m["effect"], "kinds": m["kinds"],
              "details": m["details"], "sources": m["sources"], "refs": m["n_refs"],
              "ppi": G.has_edge(a, b)}
             for (a, b), m in index.items() if a in shown and b in shown]
    kinds = sorted({k for m in mechs for k in m["kinds"]})
    return {"title": title, "best": best, "has_mech": bool(index), "kinds": kinds,
            "nodes": nodes, "pairs": pairs, "mechs": mechs, "views": vdata,
            "order": list(views)}


def draw_html(G, ranking, path, overview_nodes, views, index=None, best=None, title=""):
    data = build_data(G, ranking, overview_nodes, views, index, best, title)
    blob = json.dumps(data).replace("</", "<\\/")
    html = TEMPLATE.replace("__DATA__", blob).replace("__VIS__", VIS_JS)
    path.write_text(html, encoding="utf-8")


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PPI hub map</title>
<script src="__VIS__"></script>
<style>
  :root { --bg:#fff; --fg:#1d1d1f; --muted:#6b6b70; --line:#e3e3e6; --panel:#fafafb;
          --act:#2E9E44; --inh:#D1342F; --mix:#8E5CC4; --unk:#5A5A5A; --none:#C4C4C4; }
  * { box-sizing:border-box; }
  body { margin:0; font:13px/1.4 -apple-system, system-ui, "Segoe UI", sans-serif;
         color:var(--fg); background:var(--bg); display:flex; height:100vh; }
  #side { width:340px; min-width:280px; border-right:1px solid var(--line); background:var(--panel);
          overflow:auto; padding:14px 16px; }
  #net { flex:1; position:relative; }
  h1 { font-size:15px; margin:0 0 2px; } .sub { color:var(--muted); margin-bottom:12px; }
  h2 { font-size:11px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted);
       margin:16px 0 6px; }
  label { display:flex; align-items:center; gap:6px; padding:2px 0; cursor:pointer; }
  select, input[type=text] { width:100%; padding:5px 6px; font:inherit; border:1px solid var(--line);
       border-radius:6px; background:#fff; }
  .sw { width:22px; height:0; border-top:3px solid; display:inline-block; position:relative; }
  .sw.bar::after { content:""; position:absolute; right:-1px; top:-7px; height:11px; border-right:3px solid; border-color:inherit; }
  .sw.arr::after { content:""; position:absolute; right:-3px; top:-6px; border:4.5px solid transparent;
                   border-left:7px solid; border-left-color:inherit; border-right:0; }
  .sw.dash { border-top-style:dashed; }
  .dot { width:11px; height:11px; border-radius:50%; display:inline-block; }
  .dia { width:10px; height:10px; transform:rotate(45deg); display:inline-block; margin:0 1px; }
  .muted { color:var(--muted); } .small { font-size:12px; }
  table { border-collapse:collapse; width:100%; font-size:12px; }
  td, th { text-align:left; padding:3px 4px; border-bottom:1px solid var(--line); vertical-align:top; }
  th { color:var(--muted); font-weight:500; }
  .eff { font-weight:600; white-space:nowrap; }
  .counts span { display:inline-block; margin:0 8px 4px 0; }
  .row { display:flex; gap:6px; } .row > * { flex:1; }
  button { font:inherit; padding:4px 8px; border:1px solid var(--line); border-radius:6px;
           background:#fff; cursor:pointer; }
  .tip { white-space:pre-line; max-width:420px; font:12px/1.4 -apple-system, system-ui, sans-serif; }
  #nomech { display:none; }
  @media (max-width: 700px) { body { flex-direction:column; } #side { width:auto; height:45vh;
       border-right:0; border-bottom:1px solid var(--line); } }
</style>
</head>
<body>
<div id="side">
  <h1 id="ttl">PPI hub map</h1>
  <div class="sub" id="subttl"></div>

  <h2>View</h2>
  <select id="view"></select>
  <div class="small muted" style="margin-top:4px">Pick a target to see only its paths to your list.
    Double-click a protein to open its view.</div>
  <div style="margin-top:8px"><input type="text" id="find" placeholder="Find protein + Enter" list="names"></div>
  <datalist id="names"></datalist>

  <div id="mechbox">
    <h2>Effect (OmniPath)</h2>
    <div id="effects"></div>
    <h2>Mechanism</h2>
    <div class="row" style="margin-bottom:4px"><button id="allk">All</button><button id="nok">None</button></div>
    <div id="kinds"></div>
    <h2>Links</h2>
    <label><input type="checkbox" id="omni"> <span class="sw dash arr" style="border-color:var(--unk)"></span>
      Links only in OmniPath (not in STRING/BioGRID)</label>
  </div>
  <div id="nomech" class="small muted" style="margin-top:12px">Run with <b>--mechanism</b> to colour links
    by activation / inhibition and toggle mechanisms.</div>

  <h2>Display</h2>
  <label><input type="checkbox" id="oo"> Overview: links between two outside proteins</label>
  <label><input type="checkbox" id="iso" checked> Hide outside proteins with no visible link</label>
  <label><input type="checkbox" id="phys" checked> Physics (auto layout)</label>

  <h2>Legend</h2>
  <div class="small">
    <div><span class="dot" style="background:#4C72B0"></span> Your input protein &nbsp;
         <span class="dia" style="background:#DD8452"></span> Outside protein</div>
    <div class="muted">Thick dark border = selected target. Node size = score.
      Thicker line = found in both STRING and BioGRID.</div>
  </div>

  <div id="detail"></div>
</div>
<div id="net"></div>

<script>
const D = __DATA__;
const C = {activates:"#2E9E44", inhibits:"#D1342F", mixed:"#8E5CC4", "sign unknown":"#5A5A5A",
           "no mechanism data":"#C4C4C4"};
const EFFECTS = ["activates", "inhibits", "mixed", "sign unknown", "no mechanism data"];
const EFFECT_TEXT = {activates:"Activates", inhibits:"Inhibits", mixed:"Mixed / context-dependent",
  "sign unknown":"Acts on (sign unknown)", "no mechanism data":"Binding / association only"};
const INPUT_C = "#4C72B0", OUT_C = "#DD8452";
const $ = id => document.getElementById(id);
const key = (a, b) => a < b ? a + "|" + b : b + "|" + a;
const byId = Object.fromEntries(D.nodes.map(n => [n.id, n]));

function tip(text) { const d = document.createElement("div"); d.className = "tip"; d.innerText = text; return d; }
function esc(s) { return String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c])); }

// ---------- datasets ----------
const mechByPair = {};
for (const m of D.mechs) (mechByPair[key(m.a, m.b)] ||= []).push(m);

const edges = [];
for (const p of D.pairs) {
  const k = key(p.a, p.b), ms = (mechByPair[k] || []);
  const ppi = "PPI: " + p.sources.join(" + ");
  if (ms.length) {
    for (const m of ms) edges.push(mechEdge(m, p.sources.length, ppi));
  } else {
    edges.push({id: "p:" + k, from: p.a, to: p.b, pair: k, effect: "no mechanism data", kinds: [],
      ppi: true, width: 1 + 2 * (p.sources.length - 1), color: {color: C["no mechanism data"]},
      title: tip(p.a + " - " + p.b + "\n" + ppi + "\nNo directed mechanism in OmniPath")});
  }
}
for (const m of D.mechs) if (!m.ppi) edges.push(mechEdge(m, 1, "Not in STRING/BioGRID (OmniPath only)"));

function mechEdge(m, nsrc, ppiText) {
  const col = C[m.effect];
  return {id: "m:" + m.a + ">" + m.b, from: m.a, to: m.b, pair: key(m.a, m.b), effect: m.effect,
    kinds: m.kinds, ppi: m.ppi, dashes: !m.ppi, width: 1.5 + 2 * (nsrc - 1),
    color: {color: col, highlight: col, hover: col},
    arrows: {to: {enabled: true, type: m.effect === "inhibits" ? "bar" : "arrow", scaleFactor: 0.7}},
    smooth: {type: "curvedCW", roundness: 0.15},
    title: tip(m.a + " → " + m.b + ": " + m.effect.toUpperCase() + "\n" + m.details.join("\n") +
      "\n\nDatabases: " + m.sources + "\nReferences: " + m.refs + "\n" + ppiText)};
}

function nodeTitle(n) {
  let t = n.id + (n.input ? " (your input)" : " (outside your list)") + "\nScore " + n.score +
    "\nDirect links to " + n.direct + " inputs" + (n.partners ? ": " + n.partners : "");
  if (n.twohop) t += "\nVia one intermediate: " + n.twohop + (n.twohop_text ? "\n  " + n.twohop_text.split("), ").join(")\n  ") : "");
  return tip(t);
}
const nodes = D.nodes.map(n => ({id: n.id, label: n.id, shape: n.input ? "dot" : "diamond",
  size: 8 + 1.6 * n.score, color: {background: n.input ? INPUT_C : OUT_C, border: "#ffffff"},
  borderWidth: 1.5, font: {size: 14, face: "-apple-system, system-ui, sans-serif"}, title: nodeTitle(n)}));

const nodesDS = new vis.DataSet(nodes), edgesDS = new vis.DataSet(edges);
let visN = new Set(), visE = new Set();
const nodesView = new vis.DataView(nodesDS, {filter: n => visN.has(n.id)});
const edgesView = new vis.DataView(edgesDS, {filter: e => visE.has(e.id)});

// ---------- controls ----------
$("ttl").textContent = D.title || "PPI hub map";
$("subttl").textContent = D.nodes.filter(n => n.input).length + " input proteins" +
  (D.has_mech ? " | mechanisms from OmniPath" : "");
const sel = $("view");
sel.add(new Option("Overview (all top proteins)", ""));
for (const t of D.order) {
  const n = byId[t];
  sel.add(new Option(t + " — " + n.direct + " direct" + (n.twohop ? " + " + n.twohop + " via one" : ""), t));
}
for (const n of D.nodes) $("names").appendChild(new Option(n.id));

const effOn = new Set(EFFECTS), kindOn = new Set(D.kinds);
for (const e of EFFECTS) {
  const sw = e === "inhibits" ? "bar" : e === "no mechanism data" ? "" : "arr";
  $("effects").insertAdjacentHTML("beforeend", `<label><input type="checkbox" data-eff="${esc(e)}" checked>
    <span class="sw ${sw}" style="border-color:${C[e]}"></span> ${esc(EFFECT_TEXT[e])}</label>`);
}
for (const k of D.kinds) {
  const n = D.mechs.filter(m => m.kinds.includes(k)).length;
  $("kinds").insertAdjacentHTML("beforeend", `<label><input type="checkbox" data-kind="${esc(k)}" checked>
    ${esc(k)} <span class="muted">(${n})</span></label>`);
}
if (!D.has_mech) { $("mechbox").style.display = "none"; $("nomech").style.display = "block"; }

document.querySelectorAll("[data-eff]").forEach(cb => cb.onchange = () => {
  cb.checked ? effOn.add(cb.dataset.eff) : effOn.delete(cb.dataset.eff); refresh(); });
document.querySelectorAll("[data-kind]").forEach(cb => cb.onchange = () => {
  cb.checked ? kindOn.add(cb.dataset.kind) : kindOn.delete(cb.dataset.kind); refresh(); });
$("allk").onclick = () => { document.querySelectorAll("[data-kind]").forEach(cb => { cb.checked = true; kindOn.add(cb.dataset.kind); }); refresh(); };
$("nok").onclick = () => { document.querySelectorAll("[data-kind]").forEach(cb => { cb.checked = false; kindOn.delete(cb.dataset.kind); }); refresh(); };
$("omni").onchange = $("iso").onchange = $("oo").onchange = refresh;
function refit() {
  network.once("stabilizationIterationsDone", () => network.fit({animation: {duration: 400}}));
  network.stabilize(400);
}
sel.onchange = () => { refresh(); refit(); };
$("phys").onchange = () => network.setOptions({physics: {enabled: $("phys").checked}});
$("find").onkeydown = ev => {
  if (ev.key !== "Enter") return;
  const id = Object.keys(byId).find(n => n.toLowerCase() === ev.target.value.trim().toLowerCase());
  if (!id) return;
  if (!visN.has(id) && D.views[id]) { sel.value = id; refresh(); }
  if (visN.has(id)) { network.selectNodes([id]); network.focus(id, {scale: 1.3, animation: true}); }
};

// ---------- filtering ----------
function edgePasses(e) {
  if (!effOn.has(e.effect)) return false;
  if (e.effect === "no mechanism data") return true;
  if (!e.ppi && !$("omni").checked) return false;
  return e.kinds.some(k => kindOn.has(k));
}

function refresh() {
  const t = sel.value, view = D.views[t];
  let candN, candPairs = null;
  if (view) { candN = new Set([...view.nodes, ...D.nodes.filter(n => n.input).map(n => n.id)]);
              candPairs = new Set(view.pairs.map(p => key(p[0], p[1]))); }
  else candN = new Set(D.nodes.filter(n => n.overview).map(n => n.id));

  visE = new Set();
  for (const e of edges) {
    if (!candN.has(e.from) || !candN.has(e.to)) continue;
    if (candPairs && !candPairs.has(e.pair)) continue;
    if (!candPairs && !$("oo").checked && !byId[e.from].input && !byId[e.to].input) continue;
    if (edgePasses(e)) visE.add(e.id);
  }
  const touched = new Set();
  for (const id of visE) { const e = edgesDS.get(id); touched.add(e.from); touched.add(e.to); }
  visN = new Set([...candN].filter(n => byId[n].input || n === t || !$("iso").checked || touched.has(n)));
  for (const id of [...visE]) { const e = edgesDS.get(id); if (!visN.has(e.from) || !visN.has(e.to)) visE.delete(id); }

  const sel_ = t || D.best;
  nodesDS.update(D.nodes.map(n => ({id: n.id, borderWidth: n.id === sel_ ? 5 : 1.5,
    color: {background: n.input ? INPUT_C : OUT_C, border: n.id === sel_ ? "#111111" : "#ffffff"},
    font: {size: n.id === sel_ ? 18 : 14, bold: n.id === sel_}})));
  nodesView.refresh(); edgesView.refresh();
  showDetail(t);
}

// ---------- side panel ----------
function showDetail(t) {
  const box = $("detail");
  if (!D.views[t]) {
    const top = D.order.slice(0, 12).map(x => `<tr><td><a href="#" data-go="${esc(x)}">${esc(x)}</a></td>
      <td>${byId[x].direct}</td><td>${byId[x].twohop}</td><td>${byId[x].score}</td></tr>`).join("");
    box.innerHTML = `<h2>Top targets</h2><table><tr><th>Protein</th><th>Direct</th><th>Via one</th><th>Score</th></tr>${top}</table>`;
  } else {
    const rows = D.views[t].rows, per = {};
    for (const r of rows) (per[r.input] ||= []).push(r);
    const overall = {};
    for (const [y, rs] of Object.entries(per)) {
      const effs = new Set(rs.map(r => r.effect_on_input).filter(e => e !== "no mechanism data"));
      const signed = [...effs].filter(e => e !== "sign unknown");
      overall[y] = signed.length > 1 ? "mixed" : signed.length ? signed[0] : effs.size ? "sign unknown" : "no mechanism data";
    }
    const counts = EFFECTS.map(e => [e, Object.values(overall).filter(v => v === e).length]).filter(x => x[1]);
    const head = D.has_mech ? `<div class="counts">${counts.map(([e, n]) =>
      `<span style="color:${C[e] === C["no mechanism data"] ? "var(--muted)" : C[e]}"><b>${n}</b> ${esc(e)}</span>`).join("")}</div>` : "";
    const body = Object.entries(per).sort((a, b) => a[1][0].hops - b[1][0].hops || a[0].localeCompare(b[0])).map(([y, rs]) => {
      const via = rs.map(r => {
        const name = r.via || "direct", e = r.effect_on_input;
        const sym = {activates: "\u2192", inhibits: "\u22a3", mixed: "\u00b1", "sign unknown": "?"}[e] || "";
        return D.has_mech && sym ? `<span style="color:${C[e]}">${esc(name)} ${sym}</span>` : `<span class="muted">${esc(name)}</span>`;
      }).join(", ");
      const how = rs.map(r => (r.via ? "via " + r.via + ": " : "") + (r.how || "no directed data") +
        (r.how_reverse ? "\nREVERSE (" + y + " acts on " + t + "): " + r.how_reverse : "")).join("\n\n");
      const e = overall[y];
      return `<tr title="${esc(how)}"><td><b>${esc(y)}</b></td><td>${via}</td>
        <td class="eff" style="color:${C[e] === C["no mechanism data"] ? "var(--muted)" : C[e]}">${D.has_mech ? esc(e) : ""}</td></tr>`;
    }).join("");
    box.innerHTML = `<h2>How ${esc(t)} acts on your list</h2>${head}
      <table><tr><th>Input</th><th>Route (\u2192 activates, \u22a3 inhibits)</th><th>${D.has_mech ? "Overall" : ""}</th></tr>${body}</table>
      <div class="small muted" style="margin-top:6px">Hover a row for each step's mechanism.
      Effects of 2-step paths are the product of both steps (inhibit × inhibit = activate).</div>`;
  }
  box.querySelectorAll("[data-go]").forEach(a => a.onclick = ev => {
    ev.preventDefault(); sel.value = a.dataset.go; refresh(); refit(); });
}

// ---------- network ----------
const network = new vis.Network($("net"), {nodes: nodesView, edges: edgesView}, {
  interaction: {hover: true, tooltipDelay: 120, navigationButtons: false},
  physics: {solver: "barnesHut", barnesHut: {gravitationalConstant: -12000, springLength: 220,
            springConstant: 0.02, avoidOverlap: 0.4}, stabilization: {iterations: 400}},
  edges: {smooth: {type: "continuous"}, selectionWidth: 2, hoverWidth: 1.5},
  nodes: {borderWidth: 1.5}});
network.on("doubleClick", p => { if (p.nodes.length && D.views[p.nodes[0]]) { sel.value = p.nodes[0]; refresh(); refit(); } });
refresh(); refit();
</script>
</body>
</html>
"""
