#!/usr/bin/env python3
"""Structure evidence map: paint a predicted fold by PSM support.

Colours every residue of a predicted structure by the BEST peptide that covers it,
where "best" is the 2-D category (spectral grade x genome-uniqueness) from the manual
PSM-verification work:  grade in {excellent, good, caution, low_information}  x
uniqueness in {unique, multiple loci}.  Uncovered residues stay grey.

Inputs are decoupled so the structure source can change without touching the evidence:
  * evidence   -> PSM_Verification/genome_uniqueness/combined_peptide_table.tsv
                  (best-of-all-PSMs grade + genome_search uniqueness verdict, per peptide)
  * structure  -> any .pdb or .cif (ESMFold now; AlphaFold3 later — Biopython reads both,
                  3Dmol renders both).  Residue numbering is VALIDATED against the ORF
                  sequence, so a full-length AF3 model drops in with no edits.

Swap in AlphaFold3 outputs with one command, no code change, e.g.:
    build_evidence_map.py SEMR=af3/semr_model.cif SYNIR=af3/synir_model.cif

Outputs (in this directory):
  * evidence_map.html            self-contained interactive viewer (3Dmol inlined)
  * evidence_map.artifact.html   same body, no <html>/<head> wrapper (for Artifact publish)
  * <KEY>_residue_categories.tsv per-residue audit table
  * <KEY>_color.pml              PyMOL script for a ray-traced publication render later

Env: the ms2rescore python (has Biopython).  3Dmol-min.js + staged models live in ./models/.
"""
import os, sys, csv, json, io
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
# Default evidence = the vendored snapshot in data/ (keeps this dir self-contained); override with
# --evidence to point at the live PSM_Verification/genome_uniqueness/combined_peptide_table.tsv.
EVIDENCE_TSV = os.path.join(HERE, "data", "combined_peptide_table.tsv")
TDMOL_JS = os.path.join(HERE, "models", "3Dmol-min.js")

# ---- the two ORFs (exact sequences the combined-table positions are numbered against) ----
SEMR_SEQ = ("MQALPLVLYLSFFSLYSADLGPHQPRNLTWQVISGAGDTIWSLSKVTTPSTWWPDLFPDICKLAIGAPPGWDLEGYSDIQRAPSTPP"
            "PYVEKHPRDPWGGCSNQRDRSMLRIHPFYVCPGPHRSQSLNPKCGGKTDFFCKSWGCETSGEVYWKPSSSWDYIKVRANYSLESYVP"
            "GGFDLDECTDWCHPLRVTFTEPGKRALEWTKGYTWGLRLYKERYDEGLLFTIRLKIETPYNPLGPNKVLAPDEINPRPTPPSSSRPTP"
            "PPNSRGTVAPSKRITETVSL")
SYNIR_SEQ = ("MKGDEVSITRREMLPSAPDGIKTPMTTIRFPKAKTLQILTLLISIGVVQTSPMPHEPHNLTWVITNTATGTVINSTSHLAPIGTWFP"
             "DLVFDLCALADGTWDPLGLVWGCQHPAQHSELRRTPFYVCPGEGRSPQQTATCGGPESYFCAYWGCESTGSISWTPPIKDDLIKVQRS"
             "PGSESQGWGYGLQPKLDANDLGGPCSSNCNPITVQFTPKGKESVGWEKGKTWGLRLYVSGYDYGVMFTIQLSARQGMPLGPNPALPPV"
             "GLTVTQSQITVPISQELRPRSGPWDLLKATYLVLNNSKPELTSKCWLCLDAQPPYYEGIAVPGNYTTSTDNKDCRWQYPGKGRLTLEL"
             "IQGKGLCFGTIPHTHRHLCNSIQEAPDRDAYIIPPINTWWACTTGLTSCIHTQALNLTNGFCVLVQLAPRILRYSDEMGTRLSPFPAA"
             "RAKRATVAEILGASAFDFQDENYKVLSSAIDADLMELGKSLSKSKTSLTSLREATLRNQREQDFQSLQQDGLCKPLEKRCCTFVDNLK"
             "HARELLAKVRKRLEDRKREMEQGQPVWITVIATLSGFLGLVILIVICKLYVVNPLMDSSRTRSPPSDTY")

# per-locus config: sequence, furin SU|TM boundary (None = SU-only), TM landmarks, default model
LOCI = {
    "SEMR": dict(
        title="SEMR / Brip1os", gene="ENSMUST00000128960 · chr11 · 282 aa",
        note="Retroviral-Env-derived microprotein; SU-like, no furin/TM split — wholly secreted.",
        seq=SEMR_SEQ, boundary=None, landmarks={},
        model=os.path.join(HERE, "models", "SEMR_esmfold.pdb"), method="ESMFold v1"),
    "SYNIR": dict(
        title="SYNIR / Env", gene="ENSMUSG00000120992 · chr17 · 596 aa",
        note="Env-like ORF, furin-cleaved into SU (1–443, secreted) and TM (444–596, membrane-anchored).",
        seq=SYNIR_SEQ, boundary=443,
        landmarks={"furin RAKR": (440, 443), "CX6CC ISD": (508, 516), "MSD (TM helix)": (553, 574)},
        model=None, method="AlphaFold3"),   # provide via CLI:  SYNIR=/path/to/af3.cif
}

# ---- the evidence categories:  grade (hue) x uniqueness (dark=unique / light=multiple) ----
GRADE_RANK = {"excellent": 0, "good": 1, "caution": 2, "low_information": 3}
GRADE_LABEL = {"excellent": "Excellent", "good": "Good", "caution": "Caution", "low_information": "Low-info"}
COLORS = {
    "excellent|unique":   "#124e8c", "excellent|multiple":   "#7fb0da",
    "good|unique":        "#1c7d47", "good|multiple":        "#8fca9f",
    "caution|unique":     "#c0561a", "caution|multiple":     "#eaa06a",
    "low_information|unique": "#6a4b9c", "low_information|multiple": "#b4a3d6",
}
UNCOVERED = "#9aa4a0"


def cat_key(grade, unique):
    return f"{grade}|{'unique' if unique else 'multiple'}"


def cat_priority(grade, unique):          # lower wins: grade-major, unique-minor (matches spec order)
    return GRADE_RANK.get(grade, 9) * 2 + (0 if unique else 1)


# ---------------------------------------------------------------- evidence
def load_evidence(path=EVIDENCE_TSV):
    """peptides per locus from the combined table."""
    per = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["pos"] = int(r["pos"]); r["len"] = int(r["len"])
            r["unique"] = (r["genome_unique"].strip().lower() == "true")
            r["grade"] = r["best_grade"].strip()
            per[r["locus"]].append(r)
    return per


# ---------------------------------------------------------------- structure
THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U",
}


def parse_structure(path):
    """Return (chain_id, [ (resid:int, one_letter, plddt:float) ...], raw_text, fmt)."""
    from Bio.PDB import PDBParser, MMCIFParser
    fmt = "cif" if path.lower().endswith((".cif", ".mmcif")) else "pdb"
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)
    struct = parser.get_structure("m", path)
    model = next(iter(struct))
    chain = max(model, key=lambda c: sum(1 for _ in c))     # longest chain
    residues = []
    for res in chain:
        if "CA" not in res:
            continue
        aa = THREE2ONE.get(res.resname.strip().upper())
        if aa is None:
            continue
        residues.append((res.id[1], aa, float(res["CA"].get_bfactor())))
    raw = open(path).read()
    return chain.id, residues, raw, fmt


def align_offset(model_seq, orf_seq, key):
    """Find k so model_seq == orf_seq[k:k+len].  Full-length exact -> 0.  Warn on mismatch."""
    if model_seq == orf_seq:
        return 0
    k = orf_seq.find(model_seq)
    if k >= 0:
        print(f"  [{key}] model is residues {k+1}-{k+len(model_seq)} of the {len(orf_seq)}-aa ORF (fragment).")
        return k
    # last resort: best single-offset match (handles a few ragged residues)
    best_k, best_score = 0, -1
    for k in range(0, max(1, len(orf_seq) - len(model_seq) + 1)):
        seg = orf_seq[k:k + len(model_seq)]
        score = sum(a == b for a, b in zip(model_seq, seg))
        if score > best_score:
            best_k, best_score = k, score
    frac = best_score / max(len(model_seq), 1)
    print(f"  [{key}] WARNING: model sequence does not exactly match the ORF; "
          f"best offset {best_k} matches {frac:.0%}. Check the construct/numbering.")
    return best_k


# ---------------------------------------------------------------- per-residue categories
def residue_categories(peptides, seq, chain_residues, offset):
    """Map ORF position -> winning category + all covering peptides. Also a pdb_resid<->orf lookup."""
    orf_to_pdb = {}
    for (resid, aa, plddt) in chain_residues:
        orf_pos = resid + offset          # full-length: offset 0 -> orf_pos == resid
        orf_to_pdb[orf_pos] = resid
    win = {}                              # orf_pos -> (priority, cat, grade, unique)
    covering = defaultdict(list)          # orf_pos -> [peptide dicts]
    for p in peptides:
        pr = cat_priority(p["grade"], p["unique"])
        ck = cat_key(p["grade"], p["unique"])
        for pos in range(p["pos"], p["pos"] + p["len"]):
            covering[pos].append(p)
            if pos not in win or pr < win[pos][0]:
                win[pos] = (pr, ck, p["grade"], p["unique"])
    resinfo = {}
    plddt_by_orf = {}
    for (resid, aa, plddt) in chain_residues:
        orf_pos = resid + offset
        plddt_by_orf[orf_pos] = plddt
        ck = win[orf_pos][1] if orf_pos in win else None
        peps = [dict(peptide=q["peptide"], pos=q["pos"], grade=q["grade"], unique=q["unique"],
                     cov=q.get("best_cov", ""), ppm=q.get("best_ppm", ""), pep=q.get("best_PEP", ""),
                     dataset=q.get("best_dataset", ""), source=q.get("best_source", ""))
                for q in covering.get(orf_pos, [])]
        resinfo[orf_pos] = dict(resid=resid, aa=aa, cat=ck, plddt=round(plddt, 2), peps=peps)
    return resinfo, orf_to_pdb


# ---------------------------------------------------------------- emitters
def write_residue_tsv(key, seq, resinfo):
    path = os.path.join(HERE, f"{key}_residue_categories.tsv")
    with open(path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["orf_pos", "aa", "pdb_resid", "plddt", "category", "n_peptides", "covering_peptides"])
        for pos in sorted(resinfo):
            ri = resinfo[pos]
            w.writerow([pos, ri["aa"], ri["resid"], ri["plddt"], ri["cat"] or "uncovered",
                        len(ri["peps"]), ";".join(f"{q['peptide']}({q['grade']}/"
                        f"{'U' if q['unique'] else 'M'})" for q in ri["peps"])])
    return path


def contiguous_runs(items):
    """[(resid, cat)] sorted -> [(cat, start, end)] merged runs (same cat, consecutive resid)."""
    runs = []
    for resid, cat in items:
        if runs and runs[-1][0] == cat and resid == runs[-1][2] + 1:
            runs[-1][2] = resid
        else:
            runs.append([cat, resid, resid])
    return runs


def write_pml(key, cfg, resinfo):
    """A portable PyMOL script for a ray-traced publication figure (run wherever PyMOL exists)."""
    model_base = os.path.basename(cfg["model"])
    by_resid = sorted((ri["resid"], ri["cat"]) for ri in resinfo.values())
    runs = contiguous_runs([(r, c) for r, c in by_resid if c])
    lines = [f"# PyMOL evidence map for {key} — run:  pymol -cq {key}_color.pml",
             f"load models/{model_base}, {key}", "hide everything", "show cartoon",
             "set cartoon_transparency, 0", "bg_color white", "set ray_shadows, 0",
             "set ray_opaque_background, 0", f"color grey70, {key}"]
    for cat, colors_hex in [(c, COLORS[c]) for c in COLORS]:
        sel = "+".join(f"{s}-{e}" for cc, s, e in runs if cc == cat)
        if sel:
            lines.append(f'set_color c_{cat.replace("|","_")}, [{int(colors_hex[1:3],16)/255:.3f},'
                         f'{int(colors_hex[3:5],16)/255:.3f},{int(colors_hex[5:7],16)/255:.3f}]')
            lines.append(f'color c_{cat.replace("|","_")}, {key} and resi {sel}')
    if cfg.get("boundary"):
        lines.append(f'# SU|TM furin boundary at {cfg["boundary"]}')
    lines += ["orient", "set ray_trace_mode, 1", f"ray 1800, 1350", f"png {key}_evidence.png, dpi=300"]
    path = os.path.join(HERE, f"{key}_color.pml")
    open(path, "w").write("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- HTML build
def struct_payload(key, cfg, peptides, resinfo, raw, fmt, offset):
    """JSON-serialisable data for one structure card."""
    # colour lists for the viewer, keyed by pdb resid
    cat_resids = defaultdict(list)
    for pos, ri in resinfo.items():
        if ri["cat"]:
            cat_resids[ri["cat"]].append(ri["resid"])
    # resinfo keyed by pdb resid for hover
    hover = {ri["resid"]: dict(orf=pos, aa=ri["aa"], cat=ri["cat"], plddt=ri["plddt"],
                               peps=[f"{q['peptide']} · {GRADE_LABEL[q['grade']]} · "
                                     f"{'unique' if q['unique'] else 'multiple loci'}" for q in ri["peps"]])
             for pos, ri in resinfo.items()}
    # sequence ribbon by ORF position
    ribbon = [dict(pos=pos, aa=resinfo[pos]["aa"] if pos in resinfo else "·",
                   cat=resinfo[pos]["cat"] if pos in resinfo else None,
                   resid=resinfo[pos]["resid"] if pos in resinfo else None)
              for pos in range(1, len(cfg["seq"]) + 1)]
    covered = sum(1 for ri in resinfo.values() if ri["cat"])
    counts = defaultdict(int)
    for ri in resinfo.values():
        if ri["cat"]:
            counts[ri["cat"]] += 1
    return dict(key=key, title=cfg["title"], gene=cfg["gene"], note=cfg["note"],
                method=cfg["method"], model_name=os.path.basename(cfg["model"]),
                boundary=cfg.get("boundary"), landmarks=cfg.get("landmarks", {}),
                length=len(cfg["seq"]), covered=covered,
                pdb=raw, fmt=fmt, cat_resids=cat_resids, hover=hover, ribbon=ribbon,
                counts=counts, offset=offset,
                peptides=sorted(peptides, key=lambda p: p["pos"]))


CSS = r"""
:root{
  --bg:#eef2f0; --panel:#f8faf9; --panel2:#ffffff; --ink:#12191c; --muted:#5c6b6b;
  --hair:#cdd7d3; --accent:#2f8f83; --viewer:#f4f7f6; --chip:#e7ecea; --shadow:rgba(20,40,38,.08);
}
@media (prefers-color-scheme:dark){:root{
  --bg:#0d1214; --panel:#141b1e; --panel2:#0f1619; --ink:#e7eeeb; --muted:#8fa19d;
  --hair:#26312f; --accent:#4fb9ab; --viewer:#0b1113; --chip:#1b2426; --shadow:rgba(0,0,0,.4);}}
:root[data-theme="light"]{
  --bg:#eef2f0; --panel:#f8faf9; --panel2:#ffffff; --ink:#12191c; --muted:#5c6b6b;
  --hair:#cdd7d3; --accent:#2f8f83; --viewer:#f4f7f6; --chip:#e7ecea; --shadow:rgba(20,40,38,.08);}
:root[data-theme="dark"]{
  --bg:#0d1214; --panel:#141b1e; --panel2:#0f1619; --ink:#e7eeeb; --muted:#8fa19d;
  --hair:#26312f; --accent:#4fb9ab; --viewer:#0b1113; --chip:#1b2426; --shadow:rgba(0,0,0,.4);}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;line-height:1.5;
  -webkit-font-smoothing:antialiased}
.mono{font-family:ui-monospace,"SF Mono","Cascadia Code",Menlo,monospace;font-variant-numeric:tabular-nums}
.wrap{max-width:1160px;margin:0 auto;padding:28px 20px 64px}
header.top{border-bottom:1px solid var(--hair);padding-bottom:18px;margin-bottom:22px}
.eyebrow{font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:var(--accent);font-weight:600}
h1{font-size:26px;line-height:1.15;margin:6px 0 6px;text-wrap:balance;letter-spacing:-.01em}
.sub{color:var(--muted);max-width:70ch;font-size:14.5px}
.card{background:var(--panel);border:1px solid var(--hair);border-radius:14px;
  box-shadow:0 1px 2px var(--shadow);margin:20px 0;overflow:hidden}
.card-h{display:flex;flex-wrap:wrap;align-items:baseline;gap:8px 12px;padding:16px 18px 12px}
.card-h h2{font-size:18px;margin:0;letter-spacing:-.01em}
.gene{color:var(--muted);font-size:12.5px}
.badge{margin-left:auto;font-size:11px;font-weight:600;letter-spacing:.03em;color:var(--accent);
  border:1px solid var(--hair);background:var(--panel2);padding:3px 9px;border-radius:999px}
.note{padding:0 18px 6px;color:var(--muted);font-size:13px;max-width:80ch}
.viewrow{display:grid;grid-template-columns:1fr;gap:0}
.viewer{position:relative;width:100%;height:480px;background:var(--viewer);
  border-top:1px solid var(--hair);border-bottom:1px solid var(--hair)}
.viewer canvas{border-radius:0}
.controls{display:flex;flex-wrap:wrap;gap:6px;padding:10px 18px;align-items:center}
.btn{font:inherit;font-size:12.5px;color:var(--ink);background:var(--panel2);border:1px solid var(--hair);
  border-radius:8px;padding:5px 11px;cursor:pointer}
.btn:hover{border-color:var(--accent)}
.btn:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.btn[aria-pressed="true"]{background:var(--accent);color:#fff;border-color:var(--accent)}
.stat{display:inline-flex;gap:6px;align-items:baseline;margin-left:auto;color:var(--muted);font-size:12.5px}
.stat b{color:var(--ink);font-size:14px}
.ribbon-wrap{padding:8px 18px 16px;overflow-x:auto}
.ribbon{display:flex;flex-wrap:wrap;gap:1px;font-size:0}
.res{width:11px;height:16px;display:inline-block;border-radius:2px;cursor:pointer;position:relative}
.res:hover{outline:2px solid var(--ink);outline-offset:-1px;z-index:2}
.tick{font-size:10px;color:var(--muted);width:100%;margin:2px 0 -2px}
.legend{display:flex;flex-wrap:wrap;gap:14px 26px;padding:16px 18px;align-items:flex-start}
.lgrid{display:grid;grid-template-columns:auto auto auto;gap:6px 10px;align-items:center;font-size:12.5px}
.lgrid .hcol{font-size:10.5px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);text-align:center}
.lgrid .rlab{color:var(--muted);text-align:right;padding-right:2px}
.sw{width:26px;height:14px;border-radius:3px;border:1px solid rgba(0,0,0,.15)}
.legend .aside{font-size:12px;color:var(--muted);max-width:34ch}
table{border-collapse:collapse;width:100%;font-size:12.5px}
.tbl-wrap{overflow-x:auto;padding:4px 18px 18px}
th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--hair);white-space:nowrap}
th{font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}
td.catcell{border-left:4px solid transparent}
.pill{font-size:11px;padding:1px 7px;border-radius:999px;color:#fff;font-weight:600}
footer{color:var(--muted);font-size:12px;border-top:1px solid var(--hair);margin-top:26px;padding-top:16px;max-width:82ch}
footer code{background:var(--chip);padding:1px 5px;border-radius:4px;font-size:11.5px}
a{color:var(--accent)}
"""

JS_APP = r"""
(function(){
  const D = window.__EVIDENCE__;
  const COLORS = D.colors, UNCOVERED = D.uncovered;
  const viewers = {};
  function hx(c){return '0x'+c.slice(1);}

  function themeBg(){
    const r=document.documentElement, t=r.getAttribute('data-theme');
    const dark = t? t==='dark' : matchMedia('(prefers-color-scheme:dark)').matches;
    return dark? '#0b1113' : '#f4f7f6';
  }

  function build(struct){
    const host = document.getElementById('viewer-'+struct.key);
    if(!window.$3Dmol){ host.innerHTML='<div style="padding:20px;color:#c0561a">3D viewer failed to load.</div>'; return; }
    const v = $3Dmol.createViewer(host, {backgroundColor: themeBg()});
    v.addModel(struct.pdb, struct.fmt);
    v.setStyle({}, {cartoon:{color: hx(UNCOVERED)}});
    for(const cat in struct.cat_resids){
      v.setStyle({resi: struct.cat_resids[cat]}, {cartoon:{color: hx(COLORS[cat])}});
    }
    // furin boundary marker (SYNIR): a thin dark ring at the SU|TM cut
    if(struct.boundary){
      v.addStyle({resi:[struct.boundary]}, {cartoon:{color:'0x12191c'}});
    }
    v.setHoverable({}, true,
      function(atom,vwr){
        if(atom._lab) return;
        const info = struct.hover[atom.resi];
        if(!info) return;
        let t = info.aa+info.orf+'   '+(info.cat? catName(info.cat):'no PSM support');
        t += '   pLDDT '+info.plddt;
        if(info.peps && info.peps.length){ t += '\n'+info.peps.join('\n'); }
        atom._lab = vwr.addLabel(t, {position:atom, backgroundColor:'0x12191c', backgroundOpacity:0.92,
          fontColor:'white', fontSize:12, alignment:'bottomLeft'});
        vwr.render();
      },
      function(atom,vwr){ if(atom._lab){ vwr.removeLabel(atom._lab); delete atom._lab; vwr.render(); } });
    v.zoomTo(); v.render();
    viewers[struct.key] = {v, struct, spinning:false, hideUnc:false};
  }

  function catName(c){
    const [g,u]=c.split('|');
    const G={excellent:'Excellent',good:'Good',caution:'Caution',low_information:'Low-info'};
    return G[g]+' · '+(u==='unique'?'unique':'multiple loci');
  }

  function ribbon(struct){
    const wrap = document.getElementById('ribbon-'+struct.key);
    let html='';
    struct.ribbon.forEach(r=>{
      const col = r.cat? COLORS[r.cat] : UNCOVERED;
      const tip = r.aa+r.pos+(r.cat? ' — '+catName(r.cat):' — no support');
      html += '<span class="res" title="'+tip+'" data-resid="'+(r.resid||'')+'" '+
              'style="background:'+col+'"></span>';
    });
    wrap.innerHTML = html;
    wrap.querySelectorAll('.res').forEach(el=>{
      el.addEventListener('click',()=>{
        const rid = parseInt(el.getAttribute('data-resid'));
        const V = viewers[struct.key]; if(!V||!rid) return;
        V.v.zoomTo({resi:rid}); V.v.render();
      });
    });
  }

  function controls(struct){
    const spin=document.getElementById('spin-'+struct.key);
    const reset=document.getElementById('reset-'+struct.key);
    const hide=document.getElementById('hide-'+struct.key);
    spin.addEventListener('click',()=>{const V=viewers[struct.key]; V.spinning=!V.spinning;
      V.v.spin(V.spinning?'y':false); spin.setAttribute('aria-pressed',V.spinning);});
    reset.addEventListener('click',()=>{const V=viewers[struct.key]; V.v.zoomTo(); V.v.render();});
    hide.addEventListener('click',()=>{const V=viewers[struct.key]; V.hideUnc=!V.hideUnc;
      hide.setAttribute('aria-pressed',V.hideUnc);
      V.v.setStyle({}, {cartoon:{color:hx(UNCOVERED), hidden:V.hideUnc}});
      for(const cat in struct.cat_resids){ V.v.setStyle({resi:struct.cat_resids[cat]},{cartoon:{color:hx(COLORS[cat])}}); }
      if(struct.boundary) V.v.addStyle({resi:[struct.boundary]},{cartoon:{color:'0x12191c'}});
      V.v.render();});
  }

  function refreshBg(){ for(const k in viewers){ viewers[k].v.setBackgroundColor(themeBg()); } }
  new MutationObserver(refreshBg).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
  matchMedia('(prefers-color-scheme:dark)').addEventListener('change',refreshBg);

  D.structures.forEach(s=>{ build(s); ribbon(s); controls(s); });
})();
"""


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def html_legend():
    cols = ["Unique", "Multiple loci"]
    grades = ["excellent", "good", "caution", "low_information"]
    cells = ['<div class="hcol"></div><div class="hcol">Unique locus</div><div class="hcol">Multiple loci</div>']
    for g in grades:
        cells.append(f'<div class="rlab">{GRADE_LABEL[g]}</div>')
        for u in ("unique", "multiple"):
            cells.append(f'<div class="sw" style="background:{COLORS[g+"|"+u]}" '
                         f'title="{GRADE_LABEL[g]} · {u}"></div>')
    grid = '<div class="lgrid">' + "".join(cells) + '</div>'
    aside = ('<div class="aside"><b>Hue = spectral grade</b> (best of all PSMs). '
             '<b>Shade = genome uniqueness</b>: dark if the peptide maps only to this locus, '
             'light if it also occurs at other genomic loci (assigned here by ribo-seq translation '
             'specificity). Isobaric-only off-targets (I/L) count as unique. '
             f'<span class="sw" style="display:inline-block;vertical-align:-2px;background:{UNCOVERED}"></span> '
             'grey = no PSM support.</div>')
    return f'<div class="legend">{grid}{aside}</div>'


def html_card(s):
    lm = ""
    if s["landmarks"]:
        lm = " · ".join(f'{k} {v[0]}–{v[1]}' for k, v in s["landmarks"].items())
        lm = f'<div class="note mono" style="font-size:12px">TM landmarks: {esc(lm)}</div>'
    # peptide table
    rows = ""
    for p in s["peptides"]:
        ck = cat_key(p["grade"], p["unique"])
        col = COLORS[ck]
        sub = ""
        if s["boundary"]:
            end = p["pos"] + p["len"] - 1
            sub = "SU" if end <= s["boundary"] else ("TM" if p["pos"] > s["boundary"] else "SU/TM")
        rows += (f'<tr><td class="catcell mono" style="border-left-color:{col}">{p["pos"]}–{p["pos"]+p["len"]-1}</td>'
                 f'<td class="mono">{esc(p["peptide"])}</td><td>{sub}</td>'
                 f'<td><span class="pill" style="background:{col}">{GRADE_LABEL[p["grade"]]}</span></td>'
                 f'<td>{"unique" if p["unique"] else esc(p.get("genome_status",""))}</td>'
                 f'<td class="mono">{esc(p.get("best_cov",""))}</td>'
                 f'<td class="mono">{esc(p.get("best_ppm",""))}</td>'
                 f'<td class="mono">{esc(p.get("best_PEP",""))}</td>'
                 f'<td>{esc(p.get("best_source",""))}</td></tr>')
    table = ('<div class="tbl-wrap"><table><thead><tr>'
             '<th>Residues</th><th>Peptide</th><th>Sub</th><th>Grade</th><th>Genome</th>'
             '<th>cov</th><th>ppm</th><th>best PEP</th><th>Dataset</th></tr></thead>'
             f'<tbody>{rows}</tbody></table></div>')
    pct = round(100 * s["covered"] / s["length"])
    return f'''<section class="card">
  <div class="card-h"><h2>{esc(s["title"])}</h2><span class="gene">{esc(s["gene"])}</span>
    <span class="badge">{esc(s["method"])} · {esc(s["model_name"])}</span></div>
  <div class="note">{esc(s["note"])}</div>{lm}
  <div class="viewrow"><div class="viewer" id="viewer-{s["key"]}"></div></div>
  <div class="controls">
    <button class="btn" id="spin-{s["key"]}" aria-pressed="false">Spin</button>
    <button class="btn" id="reset-{s["key"]}">Reset view</button>
    <button class="btn" id="hide-{s["key"]}" aria-pressed="false">Hide unsupported</button>
    <span class="stat"><b>{s["covered"]}</b>/{s["length"]} residues supported ({pct}%) · <b>{len(s["peptides"])}</b> peptides</span>
  </div>
  <div class="ribbon-wrap"><div class="tick mono">sequence 1 → {s["length"]} (click a residue to focus)</div>
    <div class="ribbon" id="ribbon-{s["key"]}"></div></div>
  {table}
</section>'''


def build_html(structures):
    data = dict(colors=COLORS, uncovered=UNCOVERED, structures=structures)
    tdmol = open(TDMOL_JS).read()
    cards = "\n".join(html_card(s) for s in structures)
    methods = ", ".join(sorted({s["method"] for s in structures}))
    body = f'''<style>{CSS}</style>
<div class="wrap">
  <header class="top">
    <div class="eyebrow">PSM verification · structural evidence</div>
    <h1>Where the mass-spec evidence lands on the fold</h1>
    <p class="sub">Each residue of the predicted structure is coloured by the strongest peptide-spectrum
    match that covers it — <b>hue for spectral grade</b>, <b>shade for genome uniqueness</b>. Grey means the
    fold is predicted but no confident PSM covers it. Hover any residue for its peptide and metrics; click the
    sequence ribbon to focus the 3-D view. Structures: {esc(methods)}.</p>
  </header>
  {html_legend()}
  {cards}
  <footer>
    Evidence: best-of-all-PSMs grade (q&lt;0.01 &amp; PEP&lt;0.1; Deutsch/TransCODE rubric) ×
    genome-uniqueness verdict (6-frame exact translate-match, I/L folded), from
    <code>combined_peptide_table.tsv</code>. Overlapping peptides: the residue takes the best category
    (grade first, then unique over multiple). Numbering is validated against the ORF sequence, so an
    AlphaFold3 model swaps in unchanged. Regenerate with <code>build_evidence_map.py</code>.
  </footer>
</div>
<script>{tdmol}</script>
<script>window.__EVIDENCE__ = {json.dumps(data)};</script>
<script>{JS_APP}</script>'''
    standalone = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                  '<meta name="viewport" content="width=device-width,initial-scale=1">'
                  '<title>PSM structural evidence map</title></head><body>' + body + '</body></html>')
    open(os.path.join(HERE, "evidence_map.html"), "w").write(standalone)
    open(os.path.join(HERE, "evidence_map.artifact.html"), "w").write(body)


# ---------------------------------------------------------------- main
def main():
    overrides = {}
    for a in sys.argv[1:]:
        if "=" in a:
            k, p = a.split("=", 1)
            overrides[k.strip().upper()] = p.strip()
    evidence = load_evidence()
    structures = []
    for key, cfg in LOCI.items():
        if key in overrides:
            cfg["model"] = overrides[key]
            cfg["method"] = "AlphaFold3" if overrides[key].lower().endswith((".cif", ".mmcif")) else cfg["method"]
        if not cfg["model"] or not os.path.exists(cfg["model"]):
            print(f"[{key}] no model yet — skipping (provide with {key}=/path/to/model.cif)")
            continue
        if key not in evidence:
            print(f"[{key}] no peptides in evidence table — skipping"); continue
        print(f"[{key}] {cfg['method']}  {os.path.basename(cfg['model'])}")
        chain_id, chain_res, raw, fmt = parse_structure(cfg["model"])
        model_seq = "".join(a for _, a, _ in chain_res)
        offset = align_offset(model_seq, cfg["seq"], key)
        resinfo, orf_to_pdb = residue_categories(evidence[key], cfg["seq"], chain_res, offset)
        tsv = write_residue_tsv(key, cfg["seq"], resinfo)
        pml = write_pml(key, cfg, resinfo)
        payload = struct_payload(key, cfg, evidence[key], resinfo, raw, fmt, offset)
        structures.append(payload)
        cov = payload["covered"]
        print(f"       chain {chain_id}, {len(chain_res)} res, offset {offset}; "
              f"{cov}/{payload['length']} supported ({round(100*cov/payload['length'])}%)")
        print(f"       -> {os.path.basename(tsv)}, {os.path.basename(pml)}")
    if not structures:
        print("No structures rendered."); return
    build_html(structures)
    print(f"\n-> evidence_map.html + evidence_map.artifact.html  ({len(structures)} structure/s)")


if __name__ == "__main__":
    main()
