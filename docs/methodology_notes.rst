Methodology Notes
=================

The following methodology notes are documented for transparency.
Full detail is in ``DESIGN_PROVENANCE.md`` and in the
``methodology_notes`` field of ``validation/all_predictions_summary.json``.

Note 1 - Rosetta ΔΔG for D8PEA4 / D7BKC8 (cross-protein invalid)
------------------------------------------------------------------------

*Note 4 in DESIGN_PROVENANCE.md.*

D8PEA4 and D7BKC8 are natural IS110 orthologs included in Strategy D.
Their Rosetta ``cartesian_ddg`` values were computed using IS621 as the
reference structure - this produces cross-protein absolute energies, not
thermodynamic ΔΔG values. Values set to ``None``; ``ddg_method`` =
``rosetta_cross_protein_INVALID``.

Note 2 - Rosetta stability gate universally non-functional
---------------------------------------------------------------

*Note 5 in DESIGN_PROVENANCE.md.*

ALL 47 Rosetta ``cartesian_ddg`` values (Strategy A: 15, C: 2, D: 30) are
cross-structure absolute energies (range -31,838 to -41,308 kcal/mol),
not point-mutation ΔΔG values. The gate was not applied to any design.
Structural quality is instead proxied by ESMFold pLDDT (global mean > 70;
active-site > 70 for all triaged designs).

**Disclosure (mandatory)**: *"The Rosetta stability gate was
non-functional for all designs due to cross-protein absolute energy
computation. No designs were excluded on stability grounds. ESMFold pLDDT
is used as a structural quality proxy."*

Note 3 - MHCflurry version shift
--------------------------------------

Strategy D S_Immuno scores individually recomputed with MHCflurry 2.2.1.
Strategy B uses the IS621 conservative baseline (0.7594) rather than
per-design computation. The calibrated IS621 lockpoint (0.9255) is reported
as a secondary analysis alongside the verbatim 0.929.

Note 4 - P3 reference corrected
--------------------------------------

The earlier placeholder ``WT_IS621_S_IMMUNO = 0.250`` is incorrect.
The correct value is **0.7594**.
Delta for IS621_deimmunized_v2 = +0.1183 >= 0.10 -> PASS.

Note 5 - P5 diversity enforcement
----------------------------------------

Rank-5 was diversity-enforced: A_007 (pen_score = 0.9209, Strategy A)
replaces the natural rank-5 D023 (pen_score = 0.9319, Strategy D) to satisfy
the >= 3 strategy criterion. The top-5 is **not purely rank-ordered**.

Note 6 - Gate 8 ATLAS novelty not evaluated for Strategy B
---------------------------------------------------------------

The GENOME-ATLAS DuckDB embedding-distance gate was not evaluated in this
run (database access pending). Strategy B designs were pre-filtered for genus
diversity (700 distinct genera) at the sourcing step, providing de facto
novelty evidence. Disclosed in ``validation/P4_orthologs_result.json``.
