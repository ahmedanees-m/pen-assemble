# Changelog

All notable changes to PEN-ASSEMBLE are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.5.2] - 2026-05-25

### Added

- **Catalog schema v0.5.2** (`data/catalog_v0.5.2_current.parquet`): two new boolean columns
  required by PEN-COMPARE v3.2:
  - `intrinsic_cargo_mechanism: bool` - `True` for all 1,029 designs (all passed IS110-family
    PFAM triage: PF01548 and PF02371). IS110 bridge recombinases insert cargo as part of their
    catalytic mechanism; no external HDR donor template required.
  - `cell_based_evidence: bool` - `False` for all 1,029 designs. All designs are computational
    predictions with no peer-reviewed mammalian cell data. Enforces PEN-COMPARE v3.2
    pre-registered prediction P2: "Zero designs classified TRUE_WRITER (all cap at
    PROBABLE_WRITER due to missing cell_based_evidence)."
- `parent_editor` column now populated with canonical editor names:
  - Strategy A, C, D (ProteinMPNN): `"IS621"` (IS621-derived designs)
  - Strategy B, D (natural orthologs): UniProt accession (the design IS the editor)
- `scripts/upgrade_catalog_to_v052.py` - reproducible catalog upgrade script
- `tests/test_catalog_v052.py` - 11 new tests verifying v3.2 compatibility invariants

### Changed

- Dependency pins updated to v3.2-compatible upstream packages:
  - `genome-atlas>=0.7.1,<0.8.0` -> `>=0.7.2,<0.8.0` (ISCro4 canonical naming)
  - `mech-class>=0.5.3,<0.6.0` -> `>=0.5.4,<0.6.0` (ISCro4 holdout probe renamed)
  - `pen-score>=0.1.2,<0.2.0` -> `>=0.1.3,<0.2.0` (new `get_editor_metadata()` API)
- `parent_editor` column: no IS622 values (0 rows affected; no IS622-derived designs
  were in the catalog; column was blank in v0.5.1 and is now populated)
- Inline documentation: IS622 -> ISCro4 (deprecated alias) in markdown files where
  contextually appropriate

### Pre-registration integrity

All v0.5.0 pre-registration results unchanged. The v0.5.2 upgrade adds two columns
without modifying any PenScore, strategy, or design_id values.

- P1: **16 designs beat IS621 verbatim lockpoint 0.929** - FINAL, unchanged
- Total designs: **1,029** - unchanged

### Compatibility

- Requires pen-score v0.1.3+ for full `get_editor_metadata()` API support
- Required by PEN-COMPARE v3.2 (Gates 3 and TRUE_WRITER tier)
- Backward compatible: all v0.5.1 catalog columns preserved unchanged at same positions

---

## [0.5.1] - 2026-05-24

### Changed

**Dependency updates**
- `pen-score>=0.1.0,<0.2.0` -> `>=0.1.2,<0.2.0`: 8-axis PenScore with S_Energy axis.
  All IS110-family designs gain `S_Energy = 1.0` (no Walker A/B motifs in PF01548/PF02371).
  Weight redistribution is nearly neutral for IS110 designs: ΔPenScore ~ -0.002 to +0.001.
- `mech-class>=0.5.2,<0.6.0` -> `>=0.5.3,<0.6.0`: ISCro4 (D2TGM5) added as 6th OOD holdout
  probe; atlas pin bumped to >=0.7.1.
- `genome-atlas>=0.6.0,<0.7.0` -> `>=0.7.1,<0.8.0`: restores SIMILAR_TO/HAS_RNA/PART_OF
  edges via `graph_view='full'`; ISCro4 (formerly IS622 in Perry 2025 bioRxiv) added as System node.

**Re-scoring catalog** (v0.5.1 current best estimate)
- `data/catalog_v0.5.1_current.parquet` produced from `scripts/rescore_v012.py` (8-axis, pen-score v0.1.2).
- IS621 reference lockpoint: 0.929 -> **0.957** (S_Energy axis + weight redistribution).
- Designs beating new lockpoint (0.957): **2** (IS621_deimmunized_v2 ~ 0.966, C_targeted_001 ~ 0.959).
- `catalog/pen_assemble_catalog.parquet` (v0.5.0, 7-axis) is the **frozen pre-registration record** - NOT modified.

**Repository**
- Version bumped to `0.5.1`.
- Repository made **public** at https://github.com/ahmedanees-m/pen-assemble.
- `RESCORING_v0.1.2.md` - full accounting of v0.1.2 re-scoring (weight table, key numbers, pre-reg integrity statement).
- `scripts/rescore_v012.py` - reproducible re-scoring script.

### Pre-registration integrity

The P1 prediction (">=5 designs beat IS621 verbatim lockpoint 0.929") was pre-registered and tested
against v0.1.0 scores - **PASS, 16 designs, FINAL**. v0.5.1 reports the v0.1.2 lockpoint (0.957)
as a secondary current-best-estimate analysis only.

---

## [0.5.0] - 2026-05-20

### Added

**Design generation (Strategies A-D)**
- Strategy A: 15 domain-swap chimera designs (IS110 scaffold x guide-module combinations)
- Strategy B: 992 IS110-family ortholog candidates from NCBI, passing 7-gate triage
- Strategy C: 2 Monte Carlo deimmunized IS621 variants (best: IS621_deimmunized_v2, PenScore = 0.9673)
- Strategy D: 32 ProteinMPNN backbone-redesigned sequences conditioned on IS621 ESMFold structure (PDB: 8WT6)

**Scoring & triage**
- `pen_assemble.pen_score` - seven-axis PenScore formula with verbatim IS621 lockpoint (0.929)
- `pen_assemble.catalog` - catalog loader for the 1,029-design scorecard (CSV, Parquet)
- `pen_assemble.codon` - human-preferred codon optimisation (Kazusa table), restriction-site scanner, expression ORF builder
- `pen_assemble.api` - high-level `Designer` API covering all four strategies
- `pen_assemble.triage` - multi-gate triage and diversity analysis
- `pen_assemble.verification` - axis evaluation modules (stability, mechclass, penscore, structure, active_site)

**Catalog release artifacts** (`catalog/release_v0.5.0/`)
- `pen_assemble_catalog.{csv,parquet}` - 1,029 designs with all PenScore columns
- `p1_beaters_catalog.{csv,parquet}` - 16 verbatim IS621-beating designs
- `p5_top5_catalog.{csv,parquet}` - diversity-enforced top-5
- FASTAs: full catalog, P1 beaters, P5 top-5
- 1,029 individual design JSON files
- `browser/index.html` - interactive design browser
- `wetlab/` - 16 wet-lab reference Markdown sheets (codon-optimised ORF, restriction sites, expression notes)
- `validation/` - pre-registered prediction result JSONs (P1-P5) and combined summary
- `checksums.sha256` - SHA-256 checksums for all catalog files

**CI/CD**
- GitHub Actions CI: lint (ruff + mypy) -> test (Python 3.10/3.11/3.12) -> build -> GitHub Release
- PyPI trusted publishing workflow (`publish.yml`)
- GitHub Pages documentation workflow (`docs.yml`)

**Tests** - 63 pytest tests covering PenScore formula, catalog integrity, codon optimisation

### Changed

- Rosetta ΔΔG stability gate: non-functional for all 47 designs (cross-protein absolute energies);
  gate not threshold-comparable, so ESMFold pLDDT (>= 90 global, >= 95 active-site) used as proxy
- IS621 S_Immuno reference corrected from an earlier placeholder (0.250) to the computed value (0.7594)
- P5 rank-5 diversity-enforced: A_007 (0.9209) replaces D023 (0.9319) to ensure >= 3 strategies in top-5
- MHCflurry 2.2.1 recalibration: IS621 S_Immuno = 0.7243 under current tool version;
  calibrated lockpoint = 0.9255 (secondary analysis, 32 beaters)
- MECH-CLASS ML model (v0.5.1.dev2) misclassifies IS110-family as DSB_NUCLEASE;
  corrected via PFAM gate (PF01548 + PF02371) with `is110_reclassified = True`

### Pre-registered prediction results

| ID | Result |
|----|--------|
| P1 | **PASS** - 16 designs (PenScore > 0.929) |
| P2 | **PASS** - 1,029 designs satisfy joint constraint |
| P3 | **PASS** - Δ S_Immuno = +0.118 (threshold 0.10) |
| P4 | **PASS** - 992 PFAM-verified IS110-family candidates |
| P5 | **PASS** - A, C, D represented in diversity-enforced top-5 |

**5/5 pre-registered predictions PASS**

---

[0.5.0]: https://github.com/ahmedanees-m/pen-assemble/releases/tag/v0.5.0
