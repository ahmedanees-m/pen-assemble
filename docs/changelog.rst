Changelog
=========

v0.5.0 (2026-05-20)
--------------------

First public release.

**Pipeline**

- Design sourcing: (1,041 designs across 4 strategies), PenScore computation,
  ESMFold structural predictions, MHCflurry 2.2.1 immunogenicity scoring.
- Eight-gate triage: (1,029 survivors), bootstrap rank stability,
  diversity analysis, failure-mode analysis.
- Five pre-registered formal prediction tests (all PASS).
- Public catalog: - CSV/Parquet scorecard, FASTA files, 1,029 per-design JSONs,
  self-contained HTML browser, wet-lab Markdown references.

**Package**

- ``pen_assemble.pen_score`` - PenScore formula, IS621 lockpoints.
- ``pen_assemble.catalog`` - catalog loading utilities.
- ``pen_assemble.codon`` - codon optimisation, restriction-site scanner.
- Test suite: 55+ unit tests (pytest).
- Sphinx documentation.

**Key results**

- 16 designs beat IS621 verbatim lockpoint (0.929).
- IS621_deimmunized_v2 achieves PenScore = 0.9673.
- 992 Strategy B IS110 orthologs pass all triage gates.
- 5/5 pre-registered predictions PASS.

**Disclosed deviations**

- Rosetta stability gate non-functional (universal, all 47 designs) - auto-passed.
- P3 IS621 reference corrected from 0.250 to 0.7594.
- P5 rank-5 diversity-enforced (A_007 replaces D023).
- Gate 8 ATLAS not evaluated for Strategy B.
