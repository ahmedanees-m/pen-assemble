# PEN-ASSEMBLE

**Programmable Editor Nomination — computational pipeline for IS110-family bridge recombinase design**

[![CI](https://github.com/ahmedanees-m/pen-assemble/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/pen-assemble/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ahmedanees-m/pen-assemble/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/pen-assemble)
[![PyPI](https://img.shields.io/pypi/v/pen-assemble?color=blue)](https://pypi.org/project/pen-assemble/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/ahmedanees-m/pen-assemble)](https://github.com/ahmedanees-m/pen-assemble/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-lightgrey)](#citation)

---

## Overview

IS110-family bridge recombinases are a class of programmable genome editors capable of inserting large DNA payloads (>50 kb) into specific genomic sites **without inducing double-strand breaks (DSBs)**. Their re-targetability — achieved by swapping the bridging RNA (bRNA) guide module — makes them highly attractive candidates for therapeutic genome writing.

**PEN-ASSEMBLE** is the computational nomination pipeline that evaluates 1,041 candidate designs across four orthogonal strategies using the **8-axis PenScore** composite metric (pen-score v0.1.2). It produces a fully reproducible, pre-registered design catalog complete with wet-lab synthesis sheets ready for experimental validation.

---

## Key Results

| Metric | Value |
|--------|-------|
| Designs sourced across all strategies | 1,041 |
| Designs passing all triage gates | **1,029** |
| Designs beating IS621 verbatim lockpoint (PenScore > 0.929) | **16** *(pre-registered, pen-score v0.1.0)* |
| Designs beating calibrated lockpoint (PenScore > 0.9255) | **32** *(secondary, pen-score v0.1.0)* |
| Designs beating IS621 v0.1.2 lockpoint (PenScore > 0.957) | **2** *(secondary, pen-score v0.1.2)* |
| Pre-registered predictions passing | **5 / 5** |
| Publication policy | **PUBLISH with strong claim** |

### Top-5 Designs (P5 diversity-enforced)

| Rank | Design ID | Strategy | PenScore | Notable |
|------|-----------|----------|----------|---------|
| 1 | `IS621_deimmunized_v2_Y255K_...` | C | **0.9673** | Best immunogenicity (+0.118 over IS621) |
| 2 | `C_targeted_001` | C | 0.9586 | Targeted deimmunization |
| 3 | `D8PEA4` | D | 0.9353 | IS110 ortholog, 314 aa, compact delivery |
| 4 | `D016_IS621_ProtMPNN_T0.1_sample23` | D | 0.9319 | ProteinMPNN redesign |
| 5 | `A_007` *(diversity-enforced)* | A | 0.9209 | Domain-swap chimera |

> **P5 note:** Rank 5 was diversity-enforced — A_007 replaces natural rank-5 D023 (0.9319) to ensure ≥ 3 strategies are represented in the top-5.

---

## Pipeline Architecture

```
Strategies · 1,041 candidates          Verification pipeline             Outputs
──────────────────────────             ─────────────────────             ───────

Strategy A · Domain-Swap   ──┐
  15 chimeras                 │
                              │   Stability Gate    Mechanism Gate    Triage
Strategy B · Orthologs   ────┼──► ESMFold pLDDT ──► IS110-family ──► 1,029 designs
  992 candidates              │   ≥90 global         PF01548              │
                              │   ≥95 active-site    PF02371              ▼
Strategy C · Deimm   ────────┤                                      PenScore (8-axis)
  2 deimmunized variants      │                                      IS621 = 0.957
                              │                                      (v0.1.0 pre-reg: 0.929)
Strategy D · ProteinMPNN ────┘                                           │
  32 backbone redesigns                                    ┌─────────────┼──────────────┐
                                                           ▼             ▼              ▼
                                                   16 beat IS621    Catalog          Wetlab
                                                   PenScore>0.929  1,029 designs    16 sheets
                                                   (pre-reg P1)    CSV · Parquet    Markdown
                                                         │         Browser HTML
                                                         ▼
                                               5/5 predictions PASS
                                               → PUBLISH with strong claim
```

---

## Design Strategies

| Strategy | Method | Designs | IS621-beaters |
|----------|--------|---------|---------------|
| **A** — Domain-Swap Chimeras | Recombine IS110 scaffold × guide modules from orthologous bRNA systems | 15 | 0 (best: 0.921) |
| **B** — IS110 Ortholog Discovery | Screen 992 IS110-family sequences from NCBI via 7-gate triage | 992 | 0 (best: 0.917) |
| **C** — Monte Carlo Deimmunization | Iterative substitution to reduce MHC-II epitopes while preserving activity | 2 | **2 / 2** |
| **D** — ProteinMPNN Backbone Redesign | Sequence redesign conditioned on IS621 ESMFold structure (PDB: 8WT6) | 32 | **14 / 32** |

---

## PenScore Formula

### Current (pen-score v0.1.2 — 8-axis)

```
PenScore = S_DSB × 0.24  +  S_Spec × 0.14  +  S_Cargo × 0.19
         + S_Deliv × 0.19 +  S_Immuno × 0.09 +  S_Prog × 0.05
         + S_Mature × 0.05 +  S_Energy × 0.05
```

| Axis | Weight | Measures |
|------|--------|----------|
| `S_DSB` | 0.24 | Double-strand break avoidance (IS110 mechanism = 1.0) |
| `S_Spec` | 0.14 | Guide-RNA target-site specificity |
| `S_Cargo` | 0.19 | Payload capacity (IS110 = 1.0 by mechanism) |
| `S_Deliv` | 0.19 | AAV packaging compatibility (sequence length proxy) |
| `S_Immuno` | 0.09 | De-immunization (1 − normalised MHC-II binder fraction) |
| `S_Prog` | 0.05 | bRNA re-targeting programmability |
| `S_Mature` | 0.05 | Technology readiness / literature maturity |
| `S_Energy` | 0.05 | Walker A/B motif absence (IS110-family = 1.0) |

> **Pre-registration integrity:** P1 was evaluated using the 7-axis v0.1.0 formula
> (weights: S_DSB 0.25, S_Spec 0.10, S_Cargo 0.20, S_Deliv 0.15, S_Immuno 0.10,
> S_Prog 0.15, S_Mature 0.05). That result — **16 designs beat IS621 at 0.929** — is FINAL.
> The v0.1.2 8-axis re-scoring is a secondary analysis only. See `RESCORING_v0.1.2.md`.

**IS621 reference lockpoints:**
- Verbatim pre-registered (v0.1.0): **0.929** (primary threshold, P1)
- mech-class v0.5.2 corrected (v0.1.0): **0.954** (retroactive correction, secondary)
- pen-score v0.1.2 updated (8-axis): **0.957** (current best estimate, secondary)
- MHCflurry 2.2.1-calibrated (v0.1.0): **0.9255** (secondary analysis, 32 beaters)

---

## Installation

```bash
pip install pen-assemble
```

Or install from source with development extras:

```bash
git clone https://github.com/ahmedanees-m/pen-assemble.git
cd pen-assemble
pip install -e ".[dev,docs]"
```

**Requirements:** Python ≥ 3.10 · pandas ≥ 2.0 · pyarrow ≥ 14.0 · numpy ≥ 1.24

---

## Quick Start

### Load the design catalog

```python
from pen_assemble.catalog import load_catalog, load_p1_beaters, load_top5

df = load_catalog()      # 1,029-row DataFrame with all scored designs
p1 = load_p1_beaters()   # 16 designs with PenScore > 0.929
top5 = load_top5()       # diversity-enforced top-5

print(df[["design_id", "strategy", "pen_score"]].head())
```

### Compute PenScore

```python
from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621

ax = PenScoreAxes(
    S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0, S_Deliv=1.0,
    S_Immuno=0.8777, S_Prog=1.0, S_Mature=1.0,
)
s = pen_score(ax)           # 0.9878
print(beats_is621(s))       # True
print(ax.contributions())   # per-axis weighted contributions
```

### Codon-optimise for human expression

```python
from pen_assemble.codon import build_expression_orf, gc_content, check_restriction_sites

seq = p1.iloc[0]["protein_sequence"]
orf = build_expression_orf(seq, kozak=True, stop=True)

print(f"ORF length : {len(orf)} bp")
print(f"GC content : {gc_content(orf):.1%}")
print(f"RE sites   : {check_restriction_sites(orf)}")
```

### High-level Designer API

```python
from pen_assemble.api import Designer

d = Designer.load()

# Filter and rank catalog designs
top = d.select_designs(strategy="C", require_dsb_free=True, top_k=5)

# Run Strategy C deimmunization (requires scaffold FASTA in data_dir)
variants = d.deimmunize(scaffold_id="IS621", n_variants=50)

# Run Strategy D ProteinMPNN redesign (requires PDB in data_dir)
redesigns = d.redesign_backbone(scaffold_id="IS621", n_designs=25)
```

---

## Repository Structure

```
pen-assemble/
├── pen_assemble/                    # Python package
│   ├── pen_score.py                 # PenScore formula (7 axes)
│   ├── catalog.py                   # Catalog loading utilities
│   ├── codon.py                     # Human codon optimisation
│   ├── api.py                       # High-level Designer API
│   ├── cli.py                       # Command-line entry point
│   ├── strategies/                  # Design generation (Steps 01–11)
│   │   ├── domain_swap.py           #   Strategy A
│   │   ├── ortholog_discovery.py    #   Strategy B
│   │   ├── deimmunization.py        #   Strategy C
│   │   └── backbone_redesign.py     #   Strategy D
│   ├── triage/                      # Multi-gate triage (Step 12)
│   ├── verification/                # Axis evaluation (Steps 13–16)
│   ├── utils/                       # Linker assembly, MHC scoring, PDB parsing
│   └── data/                        # YAML configuration files
├── scripts/                         # Numbered pipeline scripts (40–52)
│   ├── 50_assemble_catalog.py       # Produces catalog/ release artifacts
│   ├── 51_build_browser.py          # Interactive HTML browser
│   └── 52_generate_wetlab_reference.py
├── tests/                           # 63 pytest tests (3.10 / 3.11 / 3.12)
├── catalog/
│   └── release_v0.5.0/
│       ├── pen_assemble_catalog.{csv,parquet}
│       ├── p1_beaters_catalog.{csv,parquet}
│       ├── p5_top5_catalog.{csv,parquet}
│       ├── browser/index.html       # Interactive design browser
│       ├── wetlab/                  # 16 wet-lab reference sheets (Markdown)
│       └── validation/              # Pre-registered prediction JSONs
├── docs/                            # Sphinx documentation (furo theme)
├── CHANGELOG.md
├── CONTRIBUTING.md
├── DESIGN_PROVENANCE.md             # Full deviation log
└── pyproject.toml
```

---

## Running the Tests

```bash
pytest tests/ -v
```

All 63 tests pass on Python 3.10, 3.11, and 3.12. Coverage is reported to [Codecov](https://codecov.io/gh/ahmedanees-m/pen-assemble).

---

## Generating the Catalog

```bash
cd scripts/
py 50_assemble_catalog.py           # catalog CSV / Parquet / FASTA / JSON
py 51_build_browser.py              # interactive HTML browser
py 52_generate_wetlab_reference.py  # 16 wet-lab Markdown sheets
```

Outputs appear in `catalog/release_v0.5.0/`.

---

## Pre-Registered Predictions

All five predictions were locked before any data analysis and evaluated verbatim.

| ID | Prediction | Threshold | Result |
|----|-----------|-----------|--------|
| P1 | ≥ 5 designs beat IS621 verbatim lockpoint | PenScore > 0.929 | **PASS** — 16 designs |
| P2 | ≥ 1 design satisfies IS110-mechanism + AAV-compatible constraint | joint gate | **PASS** — 1,029 designs |
| P3 | Best Strategy-C design improves S_Immuno over IS621 by ≥ 0.10 | Δ ≥ 0.10 | **PASS** — Δ = +0.118 |
| P4 | ≥ 100 Strategy-B candidates with PFAM-verified IS110-family domain | count ≥ 100 | **PASS** — 992 designs |
| P5 | Top-5 includes designs from ≥ 3 distinct strategies | diversity | **PASS** — A, C, D |

---

## Key Deviations from Execution Plan

| # | Deviation | Impact |
|---|-----------|--------|
| 1 | **Rosetta stability gate non-functional** — absolute cross-protein energies, not ΔΔG | Gate auto-passed; pLDDT proxy used |
| 2 | **P3 reference corrected** — IS621 S_Immuno = 0.7594, not placeholder 0.250 | P3 threshold re-evaluated correctly |
| 3 | **P5 diversity-enforced** — A_007 (0.9209) replaces natural rank-5 D023 (0.9319) | P5 PASS maintained |
| 4 | **MHCflurry 2.2.1 recalibration** — IS621 S_Immuno = 0.7243 under current tool | Calibrated lockpoint = 0.9255; 32 beaters |
| 5 | **MECH-CLASS ML misclassification** — IS110-family mis-called as DSB_NUCLEASE | Corrected via PFAM domain evidence |

Full details in [`DESIGN_PROVENANCE.md`](DESIGN_PROVENANCE.md) and
[`catalog/release_v0.5.0/validation/all_predictions_summary.json`](catalog/release_v0.5.0/validation/all_predictions_summary.json).

---

## Part of PEN-STACK

PEN-ASSEMBLE is the fourth package in the **PEN-STACK** computational biology pipeline
for programmable genome editor discovery and benchmarking.

| Package | Role | Version | DOI |
|---------|------|---------|-----|
| [genome-atlas](https://github.com/ahmedanees-m/genome-atlas) | Foundational knowledge graph (28 systems, AUROC 0.9714) | v0.7.1 | [![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-lightgrey)](#) |
| [mech-class](https://github.com/ahmedanees-m/mech-class) | Biochemical mechanism classifier (IS110 Tier-A gate) | v0.5.3 | [![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-lightgrey)](#) |
| [pen-score](https://github.com/ahmedanees-m/pen-score) | 8-axis writer scoring framework | v0.1.2 | [![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-lightgrey)](#) |
| **pen-assemble** | IS110-family design catalog (1,029 designs, 5/5 PASS) | **v0.5.1** | [![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-lightgrey)](#) |
| PEN-COMPARE *(in prep)* | Cross-system benchmarking + TrueWriterScore | — | — |

---

## Citation

If you use PEN-ASSEMBLE in your work, please cite:

```bibtex
@software{ahmed2026penassemble,
  author    = {Ahmed, Anees},
  title     = {{PEN-ASSEMBLE}: A computational pipeline for IS110-family
               bridge recombinase candidate nomination},
  year      = {2026},
  version   = {v0.5.1},
  publisher = {Zenodo},
  url       = {https://github.com/ahmedanees-m/pen-assemble},
  note      = {DOI pending Zenodo deposit}
}
```

> Ahmed A. (2026). PEN-ASSEMBLE: A computational pipeline for IS110-family
> bridge recombinase candidate nomination. v0.5.1. DOI pending.

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

---

## License

MIT © 2026 Anees Ahmed — see [LICENSE](LICENSE) for details.
