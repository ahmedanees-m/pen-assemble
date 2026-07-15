# Design Provenance

PEN-ASSEMBLE nominates IS110-family bridge recombinase designs for the human
therapeutic AAV insertion use case. This document records how the catalog of
1,029 scored designs was generated and verified.

## Design strategies

Designs come from four orthogonal strategies:

- **Strategy A (mechanism-guided chimeras).** Composite designs that pair an
  IS110-family catalytic scaffold with alternative targeting and guide modules,
  guided by the recombinase mechanism.
- **Strategy B (triage-catalog recovery).** Natural IS110-family orthologs
  recovered from a candidate catalog and retained through the multi-gate triage.
- **Strategy C (deimmunization).** IS621 variants in which surface residues are
  mutated to lower predicted MHC binding while retaining at least 80% sequence
  identity to IS621.
- **Strategy D (backbone redesign).** ProteinMPNN inverse-folding of the IS621
  backbone (structure 8WT6), with catalytic and guide-contact residues held fixed.

## Verification gates

Each candidate is evaluated through a fixed sequence of gates:

- **Mechanism (mech-class tier check).** Confirms the design classifies as an
  IS110-family, double-strand-break-free mechanism.
- **PenScore.** Composite score over the mechanistic axes, referenced to the
  IS621 verbatim lockpoint (PenScore = 0.929).
- **Structure confidence (pLDDT filter).** Predicted structures must meet a
  minimum pLDDT threshold, both globally and at the active site.
- **Stability (PyRosetta ddG gate).** Change in folding free energy relative to
  the reference backbone.
- **Immunogenicity.** Predicted MHC class I and class II binding load.

## Result

Sixteen designs exceed the IS621 verbatim lockpoint (PenScore > 0.929). The
frozen scorecard, per-design records, and wet-lab reference sheets are in
`catalog/release_v0.5.0/`.

## Limitation

All designs in this catalog are computational predictions. None has been
experimentally validated. The catalog nominates candidates for laboratory
testing; it does not report measured recombinase activity, specificity, or
immunogenicity.
