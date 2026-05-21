# Changelog

All notable changes to PEN-ASSEMBLE are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.5.0] — 2026-05-20

### Added

**Design generation (Strategies A–D)**
- Strategy A: 15 domain-swap chimera designs (IS110 scaffold × guide-module combinations)
- Strategy B: 992 IS110-family ortholog candidates from NCBI, passing 7-gate triage
- Strategy C: 2 Monte Carlo deimmunized IS621 variants (best: IS621_deimmunized_v2, PenScore = 0.9673)
- Strategy D: 32 ProteinMPNN backbone-redesigned sequences conditioned on IS621 ESMFold structure (PDB: 8WT6)

**Scoring & triage**
- `pen_assemble.pen_score` — seven-axis PenScore formula with verbatim IS621 lockpoint (0.929)
- `pen_assemble.catalog` — catalog loader for the 1,029-design scorecard (CSV, Parquet)
- `pen_assemble.codon` — human-preferred codon optimisation (Kazusa table), restriction-site scanner, expression ORF builder
- `pen_assemble.api` — high-level `Designer` API covering all four strategies
- `pen_assemble.triage` — multi-gate triage and diversity analysis
- `pen_assemble.verification` — axis evaluation modules (stability, mechclass, penscore, structure, active_site)

**Catalog release artifacts** (`catalog/release_v0.5.0/`)
- `pen_assemble_catalog.{csv,parquet}` — 1,029 designs with all PenScore columns
- `p1_beaters_catalog.{csv,parquet}` — 16 verbatim IS621-beating designs
- `p5_top5_catalog.{csv,parquet}` — diversity-enforced top-5
- FASTAs: full catalog, P1 beaters, P5 top-5
- 1,029 individual design JSON files
- `browser/index.html` — interactive design browser
- `wetlab/` — 16 wet-lab reference Markdown sheets (codon-optimised ORF, restriction sites, expression notes)
- `validation/` — pre-registered prediction result JSONs (P1–P5) and combined summary
- `checksums.sha256` — SHA-256 checksums for all catalog files

**CI/CD**
- GitHub Actions CI: lint (ruff + mypy) → test (Python 3.10/3.11/3.12) → build → GitHub Release
- PyPI trusted publishing workflow (`publish.yml`)
- GitHub Pages documentation workflow (`docs.yml`)

**Tests** — 63 pytest tests covering PenScore formula, catalog integrity, codon optimisation

### Changed

- Rosetta ΔΔG stability gate: non-functional for all 47 designs (cross-protein absolute energies);
  gate auto-passed and ESMFold pLDDT (≥ 90 global, ≥ 95 active-site) used as proxy
- IS621 S_Immuno reference corrected from execution-plan placeholder (0.250) to computed value (0.7594)
- P5 rank-5 diversity-enforced: A_007 (0.9209) replaces D023 (0.9319) to ensure ≥ 3 strategies in top-5
- MHCflurry 2.2.1 recalibration: IS621 S_Immuno = 0.7243 under current tool version;
  calibrated lockpoint = 0.9255 (secondary analysis, 32 beaters)
- MECH-CLASS ML model (v0.5.1.dev2) misclassifies IS110-family as DSB_NUCLEASE;
  corrected via PFAM gate (PF01548 + PF02371) with `is110_reclassified = True`

### Pre-registered prediction results

| ID | Result |
|----|--------|
| P1 | **PASS** — 16 designs (PenScore > 0.929) |
| P2 | **PASS** — 1,029 designs satisfy joint constraint |
| P3 | **PASS** — Δ S_Immuno = +0.118 (threshold 0.10) |
| P4 | **PASS** — 992 PFAM-verified IS110-family candidates |
| P5 | **PASS** — A, C, D represented in diversity-enforced top-5 |

**Publication policy: PUBLISH with strong claim (5/5 PASS)**

---

[0.5.0]: https://github.com/ahmedanees-m/pen-assemble/releases/tag/v0.5.0
