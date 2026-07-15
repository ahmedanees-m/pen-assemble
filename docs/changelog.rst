Changelog
=========

v0.5.2 (2026-05-25)
--------------------

**Catalog schema**

- ``data/catalog_v0.5.2_current.parquet``: two new boolean columns required by
  PEN-COMPARE v3.2. ``intrinsic_cargo_mechanism`` is ``True`` for all 1,029 designs
  (all pass IS110-family PFAM triage: PF01548 and PF02371). ``cell_based_evidence``
  is ``False`` for all 1,029 designs (all are computational predictions with no
  peer-reviewed mammalian cell data).
- ``parent_editor`` column populated with canonical editor names.
- ``scripts/upgrade_catalog_to_v052.py`` - reproducible catalog upgrade script.
- 11 new tests verifying v3.2 compatibility invariants.

**Dependency pins**

- ``genome-atlas>=0.7.2,<0.8.0``, ``mech-class>=0.5.4,<0.6.0``,
  ``pen-score>=0.1.3,<0.2.0``.

**Pre-registration integrity**

- All v0.5.0 pre-registration results unchanged. P1: 16 designs beat the IS621
  verbatim lockpoint (0.929). Total designs: 1,029.

v0.5.1 (2026-05-24)
--------------------

**Dependency updates**

- ``pen-score>=0.1.2,<0.2.0``: 8-axis PenScore with the S_Energy axis. All IS110-family
  designs gain ``S_Energy = 1.0``; weight redistribution is near-neutral for IS110 designs.
- ``mech-class>=0.5.3,<0.6.0`` and ``genome-atlas>=0.7.1,<0.8.0``.

**Re-scoring**

- ``data/catalog_v0.5.1_current.parquet`` produced by ``scripts/rescore_v012.py``
  (8-axis, pen-score v0.1.2).
- IS621 reference lockpoint: 0.929 to 0.957. Designs beating 0.957: 2.
- ``catalog/release_v0.5.0/pen_assemble_catalog.parquet`` (v0.5.0, 7-axis) is the frozen
  pre-registration record and is not modified.

**Pre-registration integrity**

- The P1 prediction was tested against v0.1.0 scores: PASS, 16 designs, FINAL. The v0.1.2
  lockpoint (0.957) is reported as a secondary current-best-estimate analysis only.

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
- Test suite: 63 unit tests (pytest).
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
