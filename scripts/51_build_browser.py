"""
Step 24: Build the interactive PEN-ASSEMBLE design browser.

Generates a self-contained HTML file (zero CDN dependencies) that loads the
design catalog inline and provides:
  - Sortable, filterable table across all 1029 triaged designs
  - Strategy filter buttons (A / B / C / D / All)
  - PenScore slider and search box
  - Expandable per-design panel: full axis breakdown, provenance, sequence
  - IS621 beater highlighting (pen_score > 0.929)
  - Download CSV button

Output: catalog/release_v0.5.0/browser/index.html

Usage:
  py 51_build_browser.py
"""
from __future__ import annotations
import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
BASE        = SCRIPTS_DIR.parent / "pipeline_results_local_test"
RELEASE     = SCRIPTS_DIR.parent / "catalog" / "release_v0.5.0"
BROWSER_DIR = RELEASE / "browser"
BROWSER_DIR.mkdir(parents=True, exist_ok=True)

AXES = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature"]

def build_catalog_json(tri: pd.DataFrame) -> str:
    records = []
    for _, row in tri.iterrows():
        rec = {
            "id":       row["design_id"],
            "strategy": row["strategy"],
            "pen":      round(float(row["pen_score"]), 4),
            "beats":    bool(row["beats_is621"]),
            "len":      int(row["protein_length_aa"]),
            "plddt":    round(float(row["final_mean_plddt"]), 1),
            "as_plddt": round(float(row["active_site_plddt"]), 1),
            "ptm":      round(float(row["ptm"]), 3),
            "axes": {a: round(float(row[a]), 4) for a in AXES},
            "tier_a":   str(row["tier_a"]),
            "composite": bool(row["composite"]),
            "comp_prob": round(float(row["composite_prob"]), 4),
            "organism": str(row["organism"]) if pd.notna(row.get("organism")) else "",
            "organism_short": (str(row["organism"])[:30] if pd.notna(row.get("organism")) else ""),
            "seq":      str(row["protein_sequence"]),
            "stab":     str(row["stability_gate_status"]),
        }
        records.append(rec)
    return json.dumps(records, separators=(",", ":"))

def build_html(catalog_json: str, n_total: int, n_beaters: int) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PEN-ASSEMBLE Design Browser v0.5.0</title>
<style>
:root {{
  --c-bg:#f8f9fa;--c-card:#fff;--c-border:#dee2e6;--c-primary:#2c7be5;
  --c-pass:#1a936f;--c-fail:#c0392b;--c-warn:#e67e22;
  --c-a:#8e44ad;--c-b:#2980b9;--c-c:#1a936f;--c-d:#d35400;
  --radius:6px;--shadow:0 1px 4px rgba(0,0,0,.1);
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:system-ui,sans-serif;background:var(--c-bg);color:#212529;font-size:14px}}
header{{background:#1a2035;color:#fff;padding:18px 24px;display:flex;align-items:center;justify-content:space-between}}
header h1{{font-size:1.3em;font-weight:700;letter-spacing:.5px}}
header .meta{{font-size:.8em;opacity:.7}}
.controls{{background:var(--c-card);border-bottom:1px solid var(--c-border);padding:12px 24px;display:flex;gap:10px;flex-wrap:wrap;align-items:center}}
.strat-btn{{border:1px solid var(--c-border);background:#fff;padding:4px 12px;border-radius:20px;cursor:pointer;font-size:.85em;font-weight:600}}
.strat-btn.active{{background:var(--c-primary);color:#fff;border-color:var(--c-primary)}}
.strat-btn[data-s="A"].active{{background:var(--c-a);border-color:var(--c-a)}}
.strat-btn[data-s="B"].active{{background:var(--c-b);border-color:var(--c-b)}}
.strat-btn[data-s="C"].active{{background:var(--c-c);border-color:var(--c-c)}}
.strat-btn[data-s="D"].active{{background:var(--c-d);border-color:var(--c-d)}}
input[type=text]{{border:1px solid var(--c-border);border-radius:var(--radius);padding:5px 10px;font-size:.9em;width:220px}}
.slider-wrap{{display:flex;align-items:center;gap:6px;font-size:.85em}}
input[type=range]{{width:100px}}
#count-badge{{margin-left:auto;font-size:.82em;color:#666}}
.table-wrap{{overflow-x:auto;padding:0 24px 24px}}
table{{width:100%;border-collapse:collapse;margin-top:12px;font-size:.85em}}
th{{background:#1a2035;color:#fff;padding:8px 10px;text-align:left;cursor:pointer;white-space:nowrap;user-select:none}}
th:hover{{background:#2c3a5a}}
th.asc::after{{content:" ▲"}}th.desc::after{{content:" ▼"}}
tr:nth-child(even){{background:#f2f4f8}}
tr.beater{{background:#eafaf1}}
tr.beater td:first-child{{border-left:3px solid var(--c-pass)}}
tr.expanded{{background:#fff7e6}}
td{{padding:7px 10px;border-bottom:1px solid var(--c-border);vertical-align:middle}}
.strat-chip{{display:inline-block;padding:1px 7px;border-radius:10px;font-size:.78em;font-weight:700;color:#fff}}
.s-A{{background:var(--c-a)}}.s-B{{background:var(--c-b)}}.s-C{{background:var(--c-c)}}.s-D{{background:var(--c-d)}}
.score-bar{{display:inline-block;height:8px;border-radius:4px;background:var(--c-primary);min-width:2px}}
.detail-row td{{padding:0;background:#fffdf5}}
.detail-panel{{padding:14px 16px;border-left:3px solid var(--c-warn)}}
.detail-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:10px;margin-top:8px}}
.detail-card{{background:var(--c-bg);border-radius:var(--radius);padding:10px 12px;border:1px solid var(--c-border)}}
.detail-card h4{{font-size:.78em;text-transform:uppercase;color:#888;margin-bottom:6px}}
.axis-row{{display:flex;justify-content:space-between;margin-bottom:3px;font-size:.85em}}
.axis-bar-wrap{{height:6px;background:#e9ecef;border-radius:3px;margin-top:2px}}
.axis-bar{{height:6px;border-radius:3px;background:var(--c-primary)}}
.seq-box{{font-family:monospace;font-size:.75em;word-break:break-all;background:#f1f3f5;padding:8px;border-radius:4px;max-height:80px;overflow-y:auto;margin-top:6px}}
.badge{{font-size:.78em;padding:2px 8px;border-radius:10px;font-weight:600}}
.badge-pass{{background:#d4edda;color:#155724}}.badge-fail{{background:#f8d7da;color:#721c24}}
.badge-warn{{background:#fff3cd;color:#856404}}
.dl-btn{{display:inline-block;margin-top:4px;font-size:.8em;padding:3px 10px;border:1px solid var(--c-primary);border-radius:4px;color:var(--c-primary);cursor:pointer;background:#fff;text-decoration:none}}
.dl-btn:hover{{background:var(--c-primary);color:#fff}}
</style>
</head>
<body>
<header>
  <div>
    <h1>PEN-ASSEMBLE Design Browser</h1>
    <div style="font-size:.82em;opacity:.8;margin-top:3px">
      v0.5.0 &nbsp;&bull;&nbsp; {n_total:,} triaged designs &nbsp;&bull;&nbsp;
      <span style="color:#7fffc4">{n_beaters} beat IS621 (PenScore &gt; 0.929)</span>
      &nbsp;&bull;&nbsp; PEN-ASSEMBLE pre-registered: 5/5 PASS
    </div>
  </div>
  <div class="meta" style="text-align:right">
    IS621 reference: 0.9290<br>
    Calibrated lockpoint: 0.9255
  </div>
</header>
<div class="controls">
  <button class="strat-btn active" data-s="ALL">All strategies</button>
  <button class="strat-btn" data-s="A">A &mdash; Chimeras (3)</button>
  <button class="strat-btn" data-s="B">B &mdash; Orthologs (992)</button>
  <button class="strat-btn" data-s="C">C &mdash; Deimmunized (2)</button>
  <button class="strat-btn" data-s="D">D &mdash; ProtMPNN (32)</button>
  <input type="text" id="search" placeholder="Search design ID or organism&hellip;">
  <div class="slider-wrap">
    PenScore &ge;
    <input type="range" id="pen-slider" min="0.90" max="0.97" step="0.001" value="0.90">
    <span id="pen-val">0.900</span>
  </div>
  <label style="font-size:.85em;display:flex;align-items:center;gap:5px">
    <input type="checkbox" id="beaters-only"> IS621 beaters only
  </label>
  <span id="count-badge"></span>
  <button class="dl-btn" id="dl-csv">&#11015; Download CSV</button>
</div>
<div class="table-wrap">
<table id="main-table">
<thead><tr>
  <th data-col="id">Design ID</th>
  <th data-col="strategy">Strat</th>
  <th data-col="pen">PenScore</th>
  <th data-col="axes.S_DSB">S_DSB</th>
  <th data-col="axes.S_Spec">S_Spec</th>
  <th data-col="axes.S_Cargo">S_Cargo</th>
  <th data-col="axes.S_Deliv">S_Deliv</th>
  <th data-col="axes.S_Immuno">S_Immuno</th>
  <th data-col="axes.S_Prog">S_Prog</th>
  <th data-col="axes.S_Mature">S_Mature</th>
  <th data-col="len">Length</th>
  <th data-col="plddt">pLDDT</th>
</tr></thead>
<tbody id="tbody"></tbody>
</table>
</div>
<script>
const RAW = {catalog_json};
const AXES = ["S_DSB","S_Spec","S_Cargo","S_Deliv","S_Immuno","S_Prog","S_Mature"];
const W = {{S_DSB:0.25,S_Spec:0.10,S_Cargo:0.20,S_Deliv:0.15,S_Immuno:0.10,S_Prog:0.15,S_Mature:0.05}};
let sortCol="pen", sortAsc=false, expanded=null;

function getVal(d,col){{
  if(col.startsWith("axes."))return d.axes[col.slice(5)];
  return d[col];
}}

function fmt(v,col){{
  if(col==="pen")return v.toFixed(4);
  if(col.startsWith("axes."))return v.toFixed(4);
  if(col==="len")return v+"aa";
  if(col==="plddt")return v.toFixed(1);
  return v;
}}

function bar(v,max=1){{
  const pct=Math.round(v/max*100);
  return `<div class="score-bar" style="width:${{pct}}px" title="${{v.toFixed(4)}}"></div>`;
}}

function renderRow(d){{
  const cls=[d.beats?"beater":"",expanded===d.id?"expanded":""].filter(Boolean).join(" ");
  const shortId = d.id.length>52 ? d.id.slice(0,49)+"..." : d.id;
  return `<tr class="${{cls}}" data-id="${{d.id}}" onclick="toggleDetail('${{d.id}}')">
    <td title="${{d.id}}">${{shortId}}</td>
    <td><span class="strat-chip s-${{d.strategy}}">${{d.strategy}}</span></td>
    <td><b>${{d.pen.toFixed(4)}}</b> ${{bar(d.pen,0.97)}} ${{d.beats?'<span class="badge badge-pass">&#10003; beats IS621</span>':''}}</td>
    ${{AXES.map(a=>`<td>${{d.axes[a].toFixed(4)}} ${{bar(d.axes[a])}}</td>`).join("")}}
    <td>${{d.len}}</td>
    <td>${{d.plddt.toFixed(1)}}</td>
  </tr>`;
}}

function renderDetail(d){{
  const axHtml = AXES.map(a=>{{
    const v=d.axes[a], pct=Math.round(v*100);
    return `<div class="axis-row"><span>${{a}}</span><b>${{v.toFixed(4)}}</b></div>
            <div class="axis-bar-wrap"><div class="axis-bar" style="width:${{pct}}%"></div></div>`;
  }}).join("");
  const seqFasta = ">"+d.id+"\\n"+d.seq.match(/.{{1,60}}/g).join("\\n");
  const strat_note={{
    A:"Domain-swap chimera (IS621 RuvC + IS621 bRNA module)",
    B:"IS110-family bridge recombinase ortholog",
    C:"Computationally deimmunized IS621 variant (surface redesign)",
    D:"ProtMPNN backbone-conditioned redesign of IS621 scaffold"
  }}[d.strategy]||"";
  return `<tr class="detail-row"><td colspan="12"><div class="detail-panel">
    <b>${{d.id}}</b>
    <div style="color:#666;margin:4px 0;font-size:.85em">${{strat_note}}</div>
    <div class="detail-grid">
      <div class="detail-card"><h4>PenScore Axes</h4>${{axHtml}}</div>
      <div class="detail-card"><h4>Structure Quality</h4>
        <div class="axis-row"><span>Global pLDDT</span><b>${{d.plddt.toFixed(1)}}</b></div>
        <div class="axis-row"><span>Active-site pLDDT</span><b>${{d.as_plddt.toFixed(1)}}</b></div>
        <div class="axis-row"><span>pTM</span><b>${{d.ptm.toFixed(3)}}</b></div>
        <div class="axis-row"><span>Length</span><b>${{d.len}} aa</b></div>
        <div style="margin-top:6px;font-size:.8em;color:#888">${{d.stab}}</div>
      </div>
      <div class="detail-card"><h4>Classification</h4>
        <div class="axis-row"><span>Tier A</span><b style="font-size:.78em">${{d.tier_a.replace("DSB_FREE_TRANSEST_RECOMBINASE","DSB-FREE RECOMBINASE")}}</b></div>
        <div class="axis-row"><span>Composite arch.</span><b>${{d.composite?'YES':'NO'}}</b></div>
        <div class="axis-row"><span>Composite prob.</span><b>${{d.comp_prob.toFixed(3)}}</b></div>
        ${{d.organism?`<div class="axis-row"><span>Organism</span><b style="font-size:.82em;font-style:italic">${{d.organism_short}}</b></div>`:""}}
      </div>
      <div class="detail-card"><h4>Sequence (${{d.len}} aa)</h4>
        <div class="seq-box">${{d.seq}}</div>
        <a class="dl-btn" href="data:text/plain,${{encodeURIComponent(seqFasta)}}"
           download="${{d.id}}.fasta" style="margin-top:6px">&#11015; FASTA</a>
        <a class="dl-btn" href="../designs/${{d.id}}.json" download="${{d.id}}.json">&#11015; JSON</a>
      </div>
    </div>
  </div></td></tr>`;
}}

function currentFiltered(){{
  const strat=document.querySelector(".strat-btn.active").dataset.s;
  const q=document.getElementById("search").value.toLowerCase();
  const minPen=parseFloat(document.getElementById("pen-slider").value);
  const beatersOnly=document.getElementById("beaters-only").checked;
  return RAW.filter(d=>
    (strat==="ALL"||d.strategy===strat) &&
    (d.pen>=minPen) &&
    (!beatersOnly||d.beats) &&
    (!q||d.id.toLowerCase().includes(q)||d.organism.toLowerCase().includes(q))
  );
}}

function render(){{
  const data=currentFiltered();
  data.sort((a,b)=>{{
    const va=getVal(a,sortCol), vb=getVal(b,sortCol);
    return sortAsc?(va>vb?1:va<vb?-1:0):(va<vb?1:va>vb?-1:0);
  }});
  let html="";
  for(const d of data){{
    html+=renderRow(d);
    if(expanded===d.id) html+=renderDetail(d);
  }}
  document.getElementById("tbody").innerHTML=html;
  document.getElementById("count-badge").textContent=
    data.length+" / "+RAW.length+" designs";
  document.querySelectorAll("th[data-col]").forEach(th=>{{
    th.classList.remove("asc","desc");
    if(th.dataset.col===sortCol) th.classList.add(sortAsc?"asc":"desc");
  }});
}}

function toggleDetail(id){{
  expanded=expanded===id?null:id;
  render();
}}

document.querySelectorAll("th[data-col]").forEach(th=>{{
  th.addEventListener("click",()=>{{
    if(sortCol===th.dataset.col) sortAsc=!sortAsc;
    else {{ sortCol=th.dataset.col; sortAsc=false; }}
    render();
  }});
}});

document.querySelectorAll(".strat-btn").forEach(btn=>{{
  btn.addEventListener("click",()=>{{
    document.querySelectorAll(".strat-btn").forEach(b=>b.classList.remove("active"));
    btn.classList.add("active");
    render();
  }});
}});

document.getElementById("search").addEventListener("input",render);
document.getElementById("pen-slider").addEventListener("input",e=>{{
  document.getElementById("pen-val").textContent=parseFloat(e.target.value).toFixed(3);
  render();
}});
document.getElementById("beaters-only").addEventListener("change",render);

document.getElementById("dl-csv").addEventListener("click",()=>{{
  const data=currentFiltered();
  const cols=["id","strategy","pen","len","plddt","beats",
    ...AXES.map(a=>"axes."+a)];
  const hdr=["design_id","strategy","pen_score","length_aa","plddt","beats_is621",
    ...AXES];
  const rows=data.map(d=>cols.map(c=>getVal(d,c)));
  const csv=[hdr.join(","),...rows.map(r=>r.join(","))].join("\\n");
  const a=document.createElement("a");
  a.href="data:text/csv,"+encodeURIComponent(csv);
  a.download="pen_assemble_filtered.csv";
  a.click();
}});

render();
</script>
</body>
</html>"""

def main() -> None:
    tri = pd.read_parquet(BASE / "part_d" / "triaged_designs.parquet")
    n_beaters = int((tri["pen_score"] > 0.929).sum())
    print(f"Building browser for {len(tri)} designs ({n_beaters} IS621 beaters)...")
    cat_json = build_catalog_json(tri)
    html = build_html(cat_json, len(tri), n_beaters)
    out = BROWSER_DIR / "index.html"
    out.write_text(html, encoding="utf-8")
    sz = out.stat().st_size
    print(f"Browser written: {out}")
    print(f"  Size: {sz:,} bytes ({sz/1024:.0f} KB)")
    print(f"  Self-contained: yes (no CDN dependencies)")
    print(f"  Features: sort, filter by strategy/pen_score/search, expandable rows, CSV download")

if __name__ == "__main__":
    main()
