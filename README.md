# PEN-ASSEMBLE v0.5.0

**Programmable Editor Nomination — Assembly pipeline for IS110-family bridge recombinase design**

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

PEN-ASSEMBLE is a computational pipeline for nominating IS110-family bridge recombinase
candidates for therapeutic genome editing applications. It evaluates 1,041 designs across
four strategies using the seven-axis **PenScore** composite metric and five pre-registered
formal predictions.

## Results Summary

| Metric | Value |
|--------|-------|
| Total designs sourced | 1,041 |
| Designs passing all triage gates | 1,029 |
| Designs beating IS621 (pen_score > 0.929) | **16** |
| Pre-registered predictions passing | **5 / 5** |
| Publication policy | **PUBLISH with strong claim** |

### Top-5 Designs (P5-compliant)

| Rank | Design | Strategy | PenScore |
|------|--------|----------|----------|
| 1 | IS621_deimmunized_v2_Y255K_... | C | 0.9673 |
| 2 | C_targeted_001 | C | 0.9586 |
| 3 | D8PEA4 | D | 0.9353 |
| 4 | D016_IS621_ProtMPNN_T0.1_sample23 | D | 0.9319 |
| 5 | A_007 *(diversity-enforced)* | A | 0.9209 |

## Installation

```bash
pip install -e ".[dev,docs]"
```

Requirements: Python ≥ 3.10, pandas ≥ 2.0, pyarrow ≥ 14.0, numpy ≥ 1.24.

## Quick Start

```python
from pen_assemble.catalog import load_catalog, load_p1_beaters
from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621
from pen_assemble.codon import build_expression_orf

# Load catalog
df = load_catalog()
p1 = load_p1_beaters()  # 16 IS621-beating designs

# Compute PenScore
ax = PenScoreAxes(S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0, S_Deliv=1.0,
                  S_Immuno=0.8777, S_Prog=1.0, S_Mature=1.0)
s = pen_score(ax)  # 0.9678
print(beats_is621(s))  # True

# Codon-optimise for human expression
orf = build_expression_orf(p1.iloc[0]["protein_sequence"], kozak=True, stop=True)
```

## PenScore Formula

```
pen_score = S_DSB×0.25 + S_Spec×0.10 + S_Cargo×0.20 + S_Deliv×0.15
          + S_Immuno×0.10 + S_Prog×0.15 + S_Mature×0.05
```

IS621 reference lockpoint: **0.929** (verbatim pre-registered).

## Running the Tests

```bash
py -m pytest tests/ -v
```

## Generating the Catalog

```bash
cd scripts/
py 50_assemble_catalog.py          # catalog CSV/Parquet/FASTA/JSON
py 51_build_browser.py             # interactive HTML browser
py 52_generate_wetlab_reference.py # wet-lab Markdown files
```

Outputs appear in `catalog/release_v0.5.0/`.

## Repository Structure

```
pen-assemble/
├── pen_assemble/          # Python package
│   ├── pen_score.py       # PenScore formula
│   ├── catalog.py         # catalog loading
│   └── codon.py           # codon optimisation
├── scripts/               # pipeline scripts (40–52)
├── tests/                 # pytest test suite
├── docs/                  # Sphinx documentation
├── catalog/
│   └── release_v0.5.0/    # public design catalog
└── DESIGN_PROVENANCE.md   # full deviation log
```

## Key Deviations

- **Rosetta stability gate non-functional** for all 47 designs (A:15, C:2, D:30) —
  values are cross-protein absolute energies, not ΔΔG. Gate auto-passed; pLDDT used as proxy.
- **P3 reference corrected**: IS621 S_Immuno = 0.7594 (not execution-plan placeholder 0.250).
- **P5 diversity-enforced**: rank-5 is A_007 (0.9209), not natural rank-5 D023 (0.9319).

Full details in `DESIGN_PROVENANCE.md` and `catalog/release_v0.5.0/validation/all_predictions_summary.json`.

## License

MIT © 2026 Anees Ahmed

## Citation

> Ahmed A. *et al.* (2026). PEN-ASSEMBLE: A computational pipeline for IS110-family
> bridge recombinase candidate nomination. *[journal pending]*. v0.5.0.
> DOI: [pending Zenodo deposit].
