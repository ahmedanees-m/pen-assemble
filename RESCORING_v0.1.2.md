# RESCORING_v0.1.2.md - pen-score v0.1.2 Re-scoring of pen-assemble Catalog

**Date:** 2026-05-24
**pen-assemble version:** v0.5.1 (current best estimate)
**pen-score version:** v0.1.2 (8 axes)
**Profile:** `human_therapeutic_aav_insertion`

---

## Summary

pen-assemble catalog v0.5.0 was scored against **pen-score v0.1.0** (7 axes, IS621 lockpoint = 0.929).
pen-score v0.1.2 introduces **two breaking changes**:

1. **IS110 Tier-A gate (mech-class v0.5.2 / v0.5.3):** IS110 proteins get S_DSB = 1.0 via
   hard gate (PF01548 and PF02371) - but pen-assemble designs already had S_DSB = 1.0 via PFAM
   evidence at confidence=0.99 (see DESIGN_PROVENANCE.md Step 15). No change for designs.

2. **S_Energy axis (8th axis, pen-score v0.1.1+):** Walker A/B motif scan on the primary
   sequence detects ATP-dependent subunits. All IS110-family designs score S_Energy = 1.0
   (no Walker A/B motifs in PF01548/PF02371 domain families). Weight redistribution shifts
   S_Prog from 0.15->0.05 and increases S_Spec (0.10->0.14) and S_Deliv (0.15->0.19).

---

## Pre-Registration Integrity

> **The pre-registered prediction P1** (">=5 designs beat IS621 verbatim lockpoint 0.929")
> **was tested against v0.1.0 scores and is FINAL as tested.**
>
> This document does **NOT** re-test P1. The v0.5.0 catalog (`catalog/pen_assemble_catalog.parquet`)
> is the pre-registration record. v0.5.1 (`data/catalog_v0.5.1_current.parquet`) is the
> current best-estimate update.

P1 primary result (v0.1.0 scores, verbatim pre-registered lockpoint 0.929):
**PASS - 16 designs exceed IS621 lockpoint (2 C + 14 D)** <- final, not re-tested.

---

## Weight Changes: v0.1.0 -> v0.1.2 (human_therapeutic_aav_insertion)

| Axis | v0.1.0 weight | v0.1.2 weight | Δ weight | IS110 axis value |
|------|--------------|--------------|---------|-----------------|
| S_DSB | 0.25 | 0.24 | -0.01 | 1.0 |
| S_Spec | 0.10 | 0.14 | **+0.04** | ~0.9891 |
| S_Cargo | 0.20 | 0.19 | -0.01 | 1.0 |
| S_Deliv | 0.15 | 0.19 | **+0.04** | 0.94-0.99 |
| S_Immuno | 0.10 | 0.09 | -0.01 | 0.76-0.88 |
| S_Prog | 0.15 | 0.05 | **-0.10** | 1.0 |
| S_Mature | 0.05 | 0.05 | 0 | varies |
| S_Energy | - | 0.05 | **+0.05** | 1.0 (IS110) |
| **Sum** | **1.00** | **1.00** | - | - |

For IS110-family designs (S_DSB = S_Cargo = S_Prog = S_Energy = 1.0):

```
ΔPenScore = -0.01x1.0 + 0.04xS_Spec - 0.01x1.0 + 0.04xS_Deliv
            - 0.01xS_Immuno - 0.10x1.0 + 0.05x1.0
          = -0.07 + 0.04xS_Spec + 0.04xS_Deliv - 0.01xS_Immuno
          ~ -0.07 + 0.04x0.989 + 0.04x0.955 - 0.01x0.810
          ~ -0.07 + 0.0396 + 0.0382 - 0.0081
          ~ -0.0003
```

**The weight redistribution is essentially neutral for IS110-family designs.**
PenScore changes across the 1,029 designs range from approximately -0.002 to +0.001.

---

## Key Numbers

| Metric | v0.1.0 (pre-reg, frozen) | v0.1.2 (current best estimate) |
|--------|--------------------------|-------------------------------|
| IS621 lockpoint | **0.929** (verbatim, pre-registered) | **0.957** (updated) |
| Designs > lockpoint | **16** (2 C + 14 D) -> **PASS** (P1) | **2** (2 C only) |
| Median PenScore delta | - | ~ -0.0003 |
| S_Energy for all IS110 | - | **1.0** |
| Axes scored | 7 | 8 |

---

## Top-5 Designs (v0.1.2)

| Rank | Design | Strategy | PenScore v0.1.0 | PenScore v0.1.2 | Δ | Beats v0.1.2 lockpoint (0.957)? |
|------|--------|---------|----------------|----------------|---|---|
| 1 | IS621_deimmunized_v2 | C | 0.9673 | ~0.9658 | -0.002 | YES |
| 2 | C_targeted_001 | C | 0.9586 | ~0.9589 | +0.0003 | YES |
| 3 | D8PEA4 (314 aa) | D | 0.9353 | ~0.9354 | +0.0001 | No |
| 4 | D016 | D | 0.9319 | ~0.9318 | -0.0001 | No |
| 5 | A_007 (281 aa) | A | 0.9209 | ~0.9195 | -0.001 | No |
| - | IS621 WT reference | - | 0.957 | **0.957** | - | (reference) |

*Approximate values; run `scripts/rescore_v012.py` for exact computed values.*

---

## Strategy B Designs (992 orthologs)

Strategy B designs have S_Mature ~ 0.0 (novel orthologs, no publications). Maximum possible
PenScore in v0.1.2 for a B design (S_Deliv = 1.0):

```
Max_B = 0.24x1.0 + 0.14x0.989 + 0.19x1.0 + 0.19x1.0
      + 0.09x0.759 + 0.05x1.0 + 0.05x0.0 + 0.05x1.0
      = 0.927
```

**No Strategy B design can exceed 0.927** in v0.1.2 (same conclusion as v0.1.0: 0.924 max).
992 ortholog designs remain below both lockpoints. P4 (>=10 B designs pass triage) unaffected.

---

## Reporting of Secondary Analysis

Against the v0.1.2 IS621 lockpoint (0.957):

| Result | Count | Details |
|--------|-------|---------|
| Designs beating 0.957 | **2** | IS621_deimmunized_v2 (~0.966), C_targeted_001 (~0.959) |
| Designs beating 0.929 (pre-reg) | **16** | Unchanged; this is the P1 primary result |

**P1 conclusion (primary):** 16 designs beat the pre-registered verbatim lockpoint 0.929.
**P1 conclusion (secondary):** Against the updated IS621 reference (0.957), only 2 C designs
exceed the lockpoint. This is reported for transparency; it does NOT invalidate P1.

---

## Catalog Files

| File | Description | Version | Do NOT modify |
|------|-------------|---------|---------------|
| `catalog/pen_assemble_catalog.parquet` | v0.5.0 frozen, pre-registration record | pen-score v0.1.0, 7-axis |  YES - pre-registration |
| `data/catalog_v0.5.1_current.parquet` | v0.5.1 current best estimate | pen-score v0.1.2, 8-axis | No |
| `results/rescore_comparison_v010_v012.csv` | Side-by-side v0.1.0 vs v0.1.2 | - | No |

To reproduce: `python scripts/rescore_v012.py --frozen catalog/pen_assemble_catalog.parquet`

---

## Usage Guidance

- **For P1 verification:** Use `catalog/pen_assemble_catalog.parquet` (frozen v0.5.0).
- **For PEN-COMPARE v3.2 analysis:** Use `data/catalog_v0.5.2_current.parquet` (adds
  `intrinsic_cargo_mechanism` and `cell_based_evidence` columns required by v3.2 Gates).
- **For PEN-COMPARE TrueWriterScore analysis:** Use `data/catalog_v0.5.1_current.parquet` or
  `data/catalog_v0.5.2_current.parquet` (v0.5.2 is a superset).
- **For fair comparison with ISCro4** (canonical; formerly "IS622" in Perry 2025 bioRxiv):
  Score ISCro4 using `pen_score.score('ISCro4')` with pen-score v0.1.3, optionally with
  `exclude_axes=['S_Mature']` for biophysical-only comparison (ISCro4 has S_Mature ~ 0.0
  as a brand-new editor; pen-assemble C/D designs inherit IS621 S_Mature ~ 0.8 from
  their parent scaffold). Note: `pen_score.score('IS622')` still works via alias
  resolution but emits a DeprecationWarning.

---

*Generated 2026-05-24. pen-assemble v0.5.1 / pen-score v0.1.2.*
*Updated 2026-05-25. pen-assemble v0.5.2 / pen-score v0.1.3 (ISCro4 canonical naming).*
*Pre-registration record: catalog/pen_assemble_catalog.parquet (DO NOT MODIFY).*
