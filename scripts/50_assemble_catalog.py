"""
Step 23: Assemble the final PEN-ASSEMBLE public design catalog.

Outputs (all in catalog/release_v0.5.0/):
  pen_assemble_catalog.csv / .parquet  -- full 1029-design scorecard
  p1_beaters_catalog.csv               -- 16 designs > IS621 0.929
  p5_top5_catalog.csv                  -- 5 P5-compliant top designs
  all_designs.fasta                    -- all 1029 sequences
  p1_beaters.fasta                     -- 16 P1-beating sequences
  top5_compliant.fasta                 -- 5 P5-compliant sequences
  designs/                             -- per-design JSON (1029 files)
  validation/                          -- P1-P5 result JSONs (copied)
  README_CATALOG.md                    -- catalog documentation

Note on PDB structures:
  ESMFold PDB files reside on VM at /home/anees_22phd0670/esm_tier1_output/pdbs/.
  VM paths are stored in the final_pdb column. Copy to release via SFTP when needed.

Usage:
  py 50_assemble_catalog.py
"""
from __future__ import annotations
import json, hashlib, shutil, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
BASE        = SCRIPTS_DIR.parent / "pipeline_results_local_test"
RELEASE     = SCRIPTS_DIR.parent / "catalog" / "release_v0.5.0"
RELEASE.mkdir(parents=True, exist_ok=True)
(RELEASE / "designs").mkdir(exist_ok=True)
(RELEASE / "validation").mkdir(exist_ok=True)

# Export columns (adapted from execution plan — no tier_b/composite_flag/ddg in actual data)
CATALOG_COLS = [
    "design_id", "strategy", "protein_length_aa",
    "tier_a", "composite", "composite_prob",
    "S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature",
    "pen_score", "beats_is621",
    "final_mean_plddt", "active_site_plddt",
    "ddg_kcal_mol", "ddg_method", "stability_gate_status",
    "organism", "genus", "protein_name",
    "gate_7_pf01548", "gate_7_pf02371", "gate_7_pass",
]

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def write_fasta(designs: pd.DataFrame, path: Path, label: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for _, row in designs.iterrows():
            seq = row["protein_sequence"]
            if not seq or (isinstance(seq, float)):
                continue
            hdr = (f">{row['design_id']} "
                   f"[strategy={row['strategy']}] "
                   f"[pen_score={row['pen_score']:.4f}] "
                   f"[length={row['protein_length_aa']}aa] "
                   f"[S_DSB={row['S_DSB']:.3f}] "
                   f"[S_Immuno={row['S_Immuno']:.4f}] "
                   f"[beats_is621={'True' if row['beats_is621'] else 'False'}]")
            f.write(hdr + "\n")
            # 60-char wrapped sequence
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    print(f"  {label}: {path.name}  ({len(designs)} sequences)")

def main() -> None:
    tri   = pd.read_parquet(BASE / "part_d" / "triaged_designs.parquet")
    top5  = pd.read_parquet(BASE / "p5_compliant_top5.parquet")
    p1_b  = tri[tri["pen_score"] > 0.929].sort_values("pen_score", ascending=False).copy()
    print(f"Triaged catalog: {len(tri)} designs")
    print(f"P1 beaters (>0.929): {len(p1_b)}")
    print(f"P5-compliant top-5: {len(top5)}")

    # 1. Scorecard CSV + Parquet
    avail_cols = [c for c in CATALOG_COLS if c in tri.columns]
    export = tri[avail_cols + ["protein_sequence"]].copy()
    export.to_csv(RELEASE / "pen_assemble_catalog.csv", index=False)
    export.to_parquet(RELEASE / "pen_assemble_catalog.parquet", index=False)
    print(f"\n  Full catalog: {len(export)} rows, {len(avail_cols)+1} cols")

    # Subset catalogs
    p1_cols = [c for c in avail_cols if c in p1_b.columns]
    p1_b[p1_cols + ["protein_sequence"]].to_csv(RELEASE / "p1_beaters_catalog.csv", index=False)
    top5_cols = [c for c in avail_cols if c in top5.columns]
    top5[top5_cols + [c for c in ["protein_sequence"] if c in top5.columns]].to_csv(
        RELEASE / "p5_top5_catalog.csv", index=False)

    # 2. FASTA files
    write_fasta(tri,  RELEASE / "all_designs.fasta",    "all designs")
    write_fasta(p1_b, RELEASE / "p1_beaters.fasta",     "P1 beaters")
    write_fasta(top5, RELEASE / "top5_compliant.fasta", "P5-compliant top-5")

    # 3. Per-design JSON metadata
    n_json = 0
    for _, row in tri.iterrows():
        meta = {}
        for col in avail_cols + ["protein_sequence", "final_pdb",
                                  "mech_class_source", "pen_score_method"]:
            if col in row.index:
                v = row[col]
                if hasattr(v, "item"):
                    v = v.item()
                meta[col] = v if not (isinstance(v, float) and v != v) else None
        (RELEASE / "designs" / f"{row['design_id']}.json").write_text(
            json.dumps(meta, indent=2, default=str))
        n_json += 1
    print(f"  Per-design JSONs: {n_json} files in designs/")

    # 4. Copy validation JSONs
    val_src = BASE / "validation"
    for jf in val_src.glob("*.json"):
        shutil.copy(jf, RELEASE / "validation" / jf.name)
    print(f"  Validation JSONs: {len(list((RELEASE/'validation').glob('*.json')))} files copied")

    # 5. Structure note
    pdb_note = (
        "ESMFold PDB structures reside on VM at "
        "/home/anees_22phd0670/esm_tier1_output/pdbs/. "
        "VM path per design is stored in the final_pdb column of "
        "pen_assemble_catalog.csv. Transfer via SFTP when preparing the Zenodo deposit."
    )
    (RELEASE / "STRUCTURES_NOTE.txt").write_text(pdb_note)
    print(f"  STRUCTURES_NOTE.txt written")

    # 6. Checksums
    checksums = {}
    for fpath in sorted(RELEASE.rglob("*")):
        if fpath.is_file() and fpath.name != "checksums.sha256":
            rel = fpath.relative_to(RELEASE)
            checksums[str(rel)] = sha256_file(fpath)
    ck_path = RELEASE / "checksums.sha256"
    with open(ck_path, "w") as f:
        for rel, h in sorted(checksums.items()):
            f.write(f"{h}  {rel}\n")
    print(f"  checksums.sha256: {len(checksums)} entries")

    print(f"\nRelease directory: {RELEASE}")
    for item in sorted(RELEASE.iterdir()):
        if item.is_dir():
            n = len(list(item.rglob("*")))
            print(f"  {item.name}/  ({n} files)")
        else:
            sz = item.stat().st_size
            print(f"  {item.name}  ({sz:,} bytes)")

if __name__ == "__main__":
    main()
