#!/usr/bin/env python3
"""
Step 12.5 — Active-site pLDDT extraction from ESMFold PDB outputs.

Run ON THE VM (where PDB files are stored):
    python3 14b_step12_active_site_plddt.py

ESMFold writes per-residue pLDDT to the B-factor column of CA atoms.
Catalytic positions are defined in IS621 full-length numbering (1-based, 342 aa):
    D11, E60, D102, D105, S241  (IS621 positions; chain-relative in 8WT6: 8,57,99,102,238)

For Strategy A/C/D: use exact IS621 catalytic positions (sequences are IS621-based, 342 aa).
For Strategy B: sequences are IS110 orthologs of varying length; catalytic positions are
    approximated by proportional mapping from IS621 positions (fraction of total length),
    clipped to ±2 residues for robustness.

Outputs:
    step12_active_site_plddt.parquet  — per-design active-site statistics
    step12b_boltz1_tier2_filtered.fasta — Strategy B top-150 by active_site_pLDDT * pTM
        (use this instead of the initial pLDDT*pTM ranking if higher specificity is wanted)
"""

import csv
import json
from pathlib import Path
from collections import defaultdict

# ── Paths (VM-relative) ──────────────────────────────────────────────────────
PDB_DIR      = Path("/home/anees_22phd0670/esm_tier1_output/pdbs")
RESULTS_CSV  = Path("/home/anees_22phd0670/esm_tier1_output/step12_results.csv")
OUTPUT_DIR   = Path("/home/anees_22phd0670/esm_tier1_output")
OUTFILE      = OUTPUT_DIR / "step12_active_site_plddt.csv"

# IS621 catalytic residues (1-based, full-length 342 aa)
IS621_CATALYTIC = [11, 60, 102, 105, 241]
IS621_LENGTH    = 342
IS621_FRACTIONS = [pos / IS621_LENGTH for pos in IS621_CATALYTIC]

# Gate threshold for active-site filter
ACTIVE_SITE_PLDDT_THRESHOLD = 75.0


def get_ca_bfactors(pdb_path: Path) -> dict[int, float]:
    """Return {residue_number: bfactor} for all CA atoms in first chain."""
    bfactors = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    res_num = int(line[22:26])
                    bfac    = float(line[60:66])
                    bfactors[res_num] = bfac
                except ValueError:
                    continue
    return bfactors


def catalytic_positions_for_strategy(strategy: str, seq_len: int) -> list[int]:
    """Return 1-based residue positions to extract for active-site pLDDT."""
    if strategy in ("A", "C", "D"):
        # IS621-based sequences: exact catalytic positions apply
        return [p for p in IS621_CATALYTIC if p <= seq_len]
    elif strategy == "B":
        # IS110 orthologs: proportional mapping from IS621 fractions
        positions = []
        for frac in IS621_FRACTIONS:
            approx = round(frac * seq_len)
            # include ±1 window to handle slight alignment offsets
            positions.extend([approx - 1, approx, approx + 1])
        # deduplicate, clamp to [1, seq_len]
        positions = sorted({max(1, min(seq_len, p)) for p in positions})
        return positions
    return []


def mean_bfactor_at_positions(bfactors: dict[int, float], positions: list[int]) -> tuple[float, float, int]:
    """Return (mean, min, n_found) for CA bfactors at given positions."""
    vals = [bfactors[p] for p in positions if p in bfactors]
    if not vals:
        return 0.0, 0.0, 0
    return sum(vals) / len(vals), min(vals), len(vals)


# ── Load results CSV ─────────────────────────────────────────────────────────
gate5_pass = {}
with open(RESULTS_CSV) as f:
    for r in csv.DictReader(f):
        if str(r["pass_gate5"]) == "True":
            gate5_pass[r["id"]] = r

print(f"Gate5 survivors to process: {len(gate5_pass)}")

# ── Extract active-site pLDDT ────────────────────────────────────────────────
results = []
missing_pdbs = []

for seq_id, row in gate5_pass.items():
    pdb_path = PDB_DIR / f"{seq_id}.pdb"
    if not pdb_path.exists():
        missing_pdbs.append(seq_id)
        continue

    strategy = row["strategy"]
    seq_len  = int(row["seq_len"])
    plddt_global = float(row["mean_plddt"])
    ptm          = float(row["ptm"])

    bfactors  = get_ca_bfactors(pdb_path)
    cat_pos   = catalytic_positions_for_strategy(strategy, seq_len)
    as_mean, as_min, n_found = mean_bfactor_at_positions(bfactors, cat_pos)

    # Fallback: if no catalytic residues found, use global mean
    if n_found == 0:
        as_mean = plddt_global
        as_min  = plddt_global
        n_found = -1  # sentinel: indicates fallback used

    pass_as = as_mean >= ACTIVE_SITE_PLDDT_THRESHOLD
    results.append({
        "id":                   seq_id,
        "strategy":             strategy,
        "seq_len":              seq_len,
        "mean_plddt_global":    f"{plddt_global:.3f}",
        "ptm":                  f"{ptm:.4f}",
        "active_site_plddt":    f"{as_mean:.3f}",
        "active_site_min":      f"{as_min:.3f}",
        "cat_residues_sampled": n_found,
        "pass_active_site":     pass_as,
        "rank_score":           as_mean * ptm,  # for Strategy B selection
    })

print(f"Processed: {len(results)} | Missing PDBs: {len(missing_pdbs)}")
if missing_pdbs:
    print(f"  Missing: {missing_pdbs[:5]}...")

# ── Write CSV ────────────────────────────────────────────────────────────────
fieldnames = ["id", "strategy", "seq_len", "mean_plddt_global", "ptm",
              "active_site_plddt", "active_site_min", "cat_residues_sampled",
              "pass_active_site", "rank_score"]
with open(OUTFILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)
print(f"Written: {OUTFILE}")

# ── Per-strategy summary ─────────────────────────────────────────────────────
from collections import defaultdict
by_strat = defaultdict(list)
for r in results:
    by_strat[r["strategy"]].append(r)

print("\n=== Active-site pLDDT gate summary (>= 75.0) ===")
for s in sorted(by_strat):
    rows_s = by_strat[s]
    n_pass = sum(1 for r in rows_s if r["pass_active_site"])
    as_vals = [float(r["active_site_plddt"]) for r in rows_s]
    print(f"  Strategy {s}: {n_pass}/{len(rows_s)} pass | "
          f"mean_AS_pLDDT={sum(as_vals)/len(as_vals):.2f} | "
          f"min={min(as_vals):.2f} max={max(as_vals):.2f}")

# ── Strategy B top-150 by active_site_pLDDT * pTM ───────────────────────────
b_results = sorted([r for r in results if r["strategy"] == "B"],
                   key=lambda x: x["rank_score"], reverse=True)
b_top150_ids = {r["id"] for r in b_results[:150]}

print(f"\nStrategy B: {len(b_results)} gate5-pass | top-150 cutoff score: {b_results[149]['rank_score']:.4f}")
print("Top 5 B by active_site_pLDDT * pTM:")
for r in b_results[:5]:
    print(f"  {r['id']:20s}  AS_pLDDT={r['active_site_plddt']}  pTM={r['ptm']}  score={r['rank_score']:.4f}")

# ── Write per-strategy JSON summary ─────────────────────────────────────────
summary = {}
for s in sorted(by_strat):
    rows_s = by_strat[s]
    n_pass = sum(1 for r in rows_s if r["pass_active_site"])
    summary[f"strategy_{s}"] = {
        "total": len(rows_s),
        "pass_active_site_plddt_75": n_pass,
        "fail": len(rows_s) - n_pass,
    }
(OUTPUT_DIR / "step12_active_site_summary.json").write_text(json.dumps(summary, indent=2))
print(f"\nSummary JSON written: {OUTPUT_DIR / 'step12_active_site_summary.json'}")
print("Done.")
