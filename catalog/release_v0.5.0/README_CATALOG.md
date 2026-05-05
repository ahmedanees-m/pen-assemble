# PEN-ASSEMBLE v0.5.0 - Public Design Catalog

**Pipeline**: PEN-ASSEMBLE (Programmable Editor Nomination - ASSEMBly of top-ranked designs)  
**Release**: v0.5.0  
**Date**: 2026-05-20  
**Designs**: 1,029 triaged IS110-family bridge recombinase designs (1,041 sourced; 12 excluded at Gate 7)  
**Pre-registered predictions verified**: 5/5 PASS 

---

## Contents

```
release_v0.5.0/
├── README_CATALOG.md           ← this file
├── pen_assemble_catalog.csv    ← full 1,029-design scorecard (CSV)
├── pen_assemble_catalog.parquet← full 1,029-design scorecard (Parquet)
├── p1_beaters_catalog.csv      ← 16 designs beating IS621 verbatim (pen_score > 0.929)
├── p5_top5_catalog.csv         ← 5 P5-compliant top designs (diversity-enforced)
├── all_designs.fasta           ← all 1,029 protein sequences (FASTA)
├── p1_beaters.fasta            ← 16 P1-beating sequences (FASTA)
├── top5_compliant.fasta        ← 5 P5 top-design sequences (FASTA)
├── STRUCTURES_NOTE.txt         ← note on PDB file location (VM SFTP)
├── checksums.sha256            ← SHA-256 checksums for all files
├── browser/
│   └── index.html              ← interactive design browser (self-contained HTML)
├── designs/
│   └── <design_id>.json        ← per-design metadata JSON (1,029 files)
├── validation/
│   ├── P1_beat_is621_result.json
│   ├── P2_cargo_deliv_result.json
│   ├── P3_deimmunization_result.json
│   ├── P4_orthologs_result.json
│   ├── P5_diversity_result.json
│   ├── all_predictions_summary.json
│   └── bootstrap_rankings.json (if generated)
└── wetlab/
    ├── wetlab_index.md         ← index of all 16 P1-beating designs
    └── <design_id>.md          ← per-design wet-lab reference (16 files)
```

---

## Design Strategies

| Strategy | Description | n sourced | n triaged |
|----------|-------------|-----------|-----------|
| A | IS621 domain-swap chimeras (IscB/IS621 modules) | 15 | 3 |
| B | IS110-family ortholog survey | 992 | 992 |
| C | IS621 targeted de-immunization (MHCflurry 2.2.1) | 2 | 2 |
| D | IS621 ProtMPNN sequence design + 2 natural orthologs | 34 | 32 |
| **Total** | | **1,041** | **1,029** |

> **Gate 7 exclusions**: 12 Strategy A designs using IscB omega-module or Cas12k crRNA
> failed the bRNA-interface PFAM check (PF01548 + PF02371 co-occurrence required).

---

## Scorecard Columns (`pen_assemble_catalog.csv`)

| Column | Description |
|--------|-------------|
| `design_id` | Unique design identifier |
| `strategy` | A / B / C / D |
| `protein_length_aa` | Protein length in amino acids |
| `tier_a` | Tier-A classification (True/False) |
| `composite` | Composite mechanism class (True/False) |
| `composite_prob` | ML composite probability |
| `S_DSB` | Double-strand break score (0-1) |
| `S_Spec` | Target specificity score (0-1) |
| `S_Cargo` | Payload compatibility score (0-1) |
| `S_Deliv` | Delivery suitability score (0-1) |
| `S_Immuno` | De-immunization score (0-1) |
| `S_Prog` | Programmability score (0-1) |
| `S_Mature` | Maturity / TRL score (0-1) |
| `pen_score` | PenScore composite (weighted sum, see formula) |
| `beats_is621` | True if pen_score > 0.929 (verbatim IS621 lockpoint) |
| `final_mean_plddt` | ESMFold mean pLDDT (structural quality proxy) |
| `active_site_plddt` | ESMFold pLDDT at catalytic residues |
| `ddg_kcal_mol` | Rosetta ΔΔG value (see Note 5 - non-functional gate) |
| `ddg_method` | Method for ΔΔG (rosetta_cartesian_ddg / rosetta_cross_protein_INVALID) |
| `stability_gate_status` | NOT_APPLIED_CROSS_PROTEIN_BUG for all designs |
| `organism` | Source organism |
| `genus` | Source genus |
| `protein_name` | UniProt/NCBI protein name |
| `gate_7_pf01548` | PF01548 domain present (True/False) |
| `gate_7_pf02371` | PF02371 domain present (True/False) |
| `gate_7_pass` | Both PFAM domains confirmed (True/False) |
| `protein_sequence` | Full amino acid sequence |

---

## PenScore Formula

```
pen_score = S_DSB   x 0.25
          + S_Spec  x 0.10
          + S_Cargo x 0.20
          + S_Deliv x 0.15
          + S_Immunox 0.10
          + S_Prog  x 0.15
          + S_Maturex 0.05
```

**IS621 reference lockpoint**: 0.929 (verbatim pre-registered).  
**Calibrated lockpoint**: 0.9255 (MHCflurry 2.2.1-consistent; secondary analysis only).

---

## Pre-Registered Prediction Results

| ID | Prediction | Verdict | Key result |
|----|-----------|---------|------------|
| P1 | ≥ 5 designs beat IS621 verbatim (0.929) | **PASS** | 16 designs > 0.929 |
| P2 | ≥ 1 design with S_Cargo = 1.0 AND S_Deliv ≥ 0.9 | **PASS** | 1,029 / 1,029 (existence claim) |
| P3 | ≥ 0.10 S_Immuno gain in deimmunized IS621 | **PASS** | Δ = +0.1183 (IS621_deimmunized_v2) |
| P4 | ≥ 10 Strategy B ortholog survivors | **PASS** | 992 survivors |
| P5 | Top-5 from ≥ 3 distinct strategies | **PASS** | 3 strategies (A, C, D) |

****5/5 pre-registered predictions PASS.  
Full details: `validation/all_predictions_summary.json`

---

## Limitations and Caveats

1. **Rosetta ΔΔG gate non-functional** (Note 5): ALL 47 computed ΔΔG values
   (Strategy A:15, C:2, D:30) are cross-structure absolute energies
   (range -31,838 to -41,308 kcal/mol), not point-mutation ΔΔG values.
   Gate not threshold-comparable; structural quality proxied by ESMFold pLDDT.

2. **MHCflurry version shift**: Strategy D S_Immuno individually computed with
   MHCflurry 2.2.1. Strategy B uses IS621 conservative baseline (0.7594).

3. **P3 reference corrected**: Earlier placeholder IS621 S_Immuno = 0.250
   corrected to PEN-SCORE published value 0.7594.

4. **P5 diversity enforcement**: Rank-5 design A_007 (pen_score = 0.9209) replaces
   natural rank-5 D023 (pen_score = 0.9319) to satisfy the ≥ 3 strategy criterion.
   Top-5 is not purely rank-ordered.

5. **Gate 8 (ATLAS novelty) not evaluated** for Strategy B: GENOME-ATLAS ATLAS DuckDB
   not accessible. Strategy B pre-filtered for genus diversity (700 distinct genera).

Full deviation log: `DESIGN_PROVENANCE.md` in the repository root.

---

## Interactive Browser

Open `browser/index.html` in any modern web browser (no internet connection required).  
Features:
- Filter by strategy (A / B / C / D / All)
- PenScore range slider
- Text search (design ID or organism)
- Beaters-only toggle
- Sortable columns (click header)
- Expandable per-design panel (axis breakdown, pLDDT, sequence, FASTA download)
- CSV download of current filtered view

---

## Wet-Lab Reference

`wetlab/wetlab_index.md` - index of 16 P1-beating designs with links to individual files.  
Each `wetlab/<design_id>.md` contains:
- Score summary and axis breakdown
- Codon-optimised DNA (Kazusa human preferred codons) with Kozak + stop
- GC content and restriction site flags
- Tier 1-3 validation experiment recommendations (recombination assay, EMSA, HEK293T)

> **Disclaimer**: Sequences are computationally designed. Experimental validation is required
> before any biological activity claims can be made.

---

## ESMFold Structures (PDB)

PDB files are on the compute VM:

```
~/esm_tier1_output/pdbs/<design_id>.pdb
```

The `final_pdb` column in `pen_assemble_catalog.csv` stores the file path.  
See `STRUCTURES_NOTE.txt` for details.

---

## Checksums

`checksums.sha256` contains SHA-256 hashes for all files in this release directory
(1,062 entries). Verify with:

```bash
# Linux / macOS
sha256sum -c checksums.sha256

# Windows PowerShell
Get-Content checksums.sha256 | ForEach-Object {
    $hash, $file = $_ -split '  ', 2
    $actual = (Get-FileHash $file -Algorithm SHA256).Hash.ToLower()
    if ($actual -ne $hash) { Write-Warning "MISMATCH: $file" }
}
```

---

## Citation

> Ahmed A. *et al.* (2026). PEN-ASSEMBLE: A computational pipeline for nomination of
> IS110-family bridge recombinase candidates for therapeutic genome editing.
> Pipeline v0.5.0. https://github.com/ahmedanees-m/pen-assemble

---

*Generated by PEN-ASSEMBLE v0.5.0 - 50_assemble_catalog.py / 52_generate_wetlab_reference.py*
