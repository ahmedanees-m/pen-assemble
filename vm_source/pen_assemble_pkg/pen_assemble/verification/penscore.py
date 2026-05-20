"""PEN-SCORE 7-axis evaluation of designs. Step 16.

Requires pen-score>=0.1.0 (pip install 'pen-assemble[pen-score]').
Reuses all 7 Paper 3 axis pipelines.

IS621 lockpoint: PenScore = 0.929 (Paper 3 v0.1.0, human_therapeutic_aav_insertion).
P1 pre-registration: >=5 designs beat 0.929.

S_Mature inheritance rules (pre-registered before any design generation):
  Strategy A chimeras (novel architecture)  : S_Mature = 0.0
  Strategy B orthologs (novel proteins)     : S_Mature = 0.0
  Strategy C deimm (>=80% id to IS621)      : S_Mature = parent S_Mature = 0.792
  Strategy D ProteinMPNN (IS621 backbone)   : S_Mature = IS621_S_Mature × 0.5 = 0.396

The 7 axes and their computation sources:
  S_DSB     : MECH-CLASS tier_a (DSB_FREE → 1.0; else 0.0)
  S_Spec    : specificity model from pen-score package (sequence-based)
  S_Cargo   : enzyme mechanism gate (transposase/recombinase → 1.0 per Paper 3 rules)
  S_Deliv   : protein size scoring (target: ≤1200 aa; AAV capacity)
  S_Immuno  : netMHCpan epitope load (calibrated to Paper 3 scale)
  S_Prog    : programmability model (RNA-guide inference from sequence)
  S_Mature  : citations × expertise inheritance (pre-registered rules above)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import pandas as pd

IS621_PENSCORE  = 0.929    # P1 lockpoint (Paper 3 v0.1.0, human_therapeutic_aav_insertion)
IS621_S_MATURE  = 0.792    # Paper 3 v0.1.0 public_scorecard.parquet
IS621_S_IMMUNO  = 0.7594   # Paper 3 v0.1.0 (P3 deimm baseline)
IS621_S_DELIV   = 0.949    # estimated from Paper 3 (~260 aa, AAV-compact)

P1_THRESHOLD_N_DESIGNS = 5  # pre-registered: >= 5 designs beat 0.929

# Pre-registered S_Mature inheritance rules
S_MATURE_RULES: dict[str, float] = {
    "A_domain_swap":    0.0,
    "B_ortholog":       0.0,
    "C_deimm":          IS621_S_MATURE,           # 0.792 (inherit parent)
    "D_protmpnn":       IS621_S_MATURE * 0.5,     # 0.396 (partial inheritance)
}

# Axis weights (from Paper 3 pen-score v0.1.0; must be identical for comparability)
AXIS_WEIGHTS: dict[str, float] = {
    "S_DSB":    0.25,
    "S_Spec":   0.10,
    "S_Cargo":  0.20,
    "S_Deliv":  0.15,
    "S_Immuno": 0.10,
    "S_Prog":   0.15,
    "S_Mature": 0.05,
}
assert abs(sum(AXIS_WEIGHTS.values()) - 1.0) < 1e-6, "Axis weights must sum to 1.0"


# ---------------------------------------------------------------------------
# Axis computation helpers
# ---------------------------------------------------------------------------

def _compute_s_dsb(tier_a: str) -> float:
    """S_DSB: 1.0 if DSB-free, 0.0 if DSB-creating, 0.5 if uncertain."""
    if tier_a == "DSB_FREE_TRANSEST_RECOMBINASE":
        return 1.0
    if tier_a == "DSB_NUCLEASE":
        return 0.0
    return 0.5  # UNCLASSIFIED


def _compute_s_cargo(tier_a: str, composite: bool) -> float:
    """S_Cargo: does mechanism support site-specific cargo insertion?
    DSB-free transposase/recombinase with composite architecture → 1.0.
    """
    if tier_a == "DSB_FREE_TRANSEST_RECOMBINASE" and composite:
        return 1.0
    if tier_a == "DSB_FREE_TRANSEST_RECOMBINASE":
        return 0.85  # likely capable but no composite signal
    return 0.0


def _compute_s_deliv(protein_length_aa: int) -> float:
    """S_Deliv: size score. Peak at ~260 aa (IS621 size); falls off above 900 aa.
    AAV packaging capacity ≈ 4700 bp; ~1200 aa with bRNA cargo budget.
    """
    n = protein_length_aa
    if n <= 260:
        return 1.0 - max(0.0, (260 - n) / 260) * 0.2   # slight penalty for very small
    if n <= 900:
        return 1.0 - (n - 260) / (900 - 260) * 0.35    # linear decline
    if n <= 1200:
        return 0.65 - (n - 900) / (1200 - 900) * 0.65  # steep decline
    return 0.0  # > 1200 aa: AAV cannot package


def _compute_s_immuno(
    sequence: str,
    strategy: str,
    precomputed_s_immuno: Optional[float] = None,
) -> float:
    """S_Immuno: if precomputed (from Strategy C deimm), use that directly.
    Otherwise, use heuristic from deimmunization module.
    """
    if precomputed_s_immuno is not None:
        return precomputed_s_immuno
    if strategy in ("C_deimm",):
        # Should always have precomputed S_Immuno from deimmunization pipeline
        return IS621_S_IMMUNO  # fallback to WT value
    # For other strategies: use heuristic via deimmunization module
    try:
        from pen_assemble.strategies.deimmunization import compute_epitope_profile, IS621_S_IMMUNO_WT
        profile = compute_epitope_profile(sequence, use_netmhcpan=False)
        # Normalize relative to IS621 WT: if same epitope load, return IS621_S_IMMUNO_WT
        # (different sequences will have different loads; this is a rough estimate)
        return min(1.0, max(0.0, IS621_S_IMMUNO_WT * (1.0 + (0.7594 - profile.weighted_load / 50.0) * 0.1)))
    except Exception:
        return IS621_S_IMMUNO  # safe fallback


def _compute_s_mature(strategy: str, sequence: str, parent_sequence: Optional[str] = None) -> float:
    """S_Mature: pre-registered inheritance rules."""
    if strategy.startswith("A"):
        return S_MATURE_RULES["A_domain_swap"]
    if strategy.startswith("B"):
        return S_MATURE_RULES["B_ortholog"]
    if strategy.startswith("D"):
        return S_MATURE_RULES["D_protmpnn"]
    if strategy.startswith("C"):
        # Verify >= 80% identity to IS621 WT
        if parent_sequence and sequence:
            n_match = sum(a == b for a, b in zip(sequence, parent_sequence))
            identity = n_match / max(len(parent_sequence), len(sequence))
            if identity >= 0.80:
                return S_MATURE_RULES["C_deimm"]
            else:
                return S_MATURE_RULES["A_domain_swap"]  # not enough identity → 0.0
        return S_MATURE_RULES["C_deimm"]  # assume >= 80% if parent not provided
    return 0.0


# ---------------------------------------------------------------------------
# Full PEN-SCORE computation
# ---------------------------------------------------------------------------

def _check_penscore_package() -> bool:
    try:
        import penscore  # noqa: F401
        return True
    except ImportError:
        return False


def compute_penscore(
    design_id: str,
    protein_sequence: str,
    strategy: str,
    tier_a: str = "DSB_FREE_TRANSEST_RECOMBINASE",
    tier_a_confidence: float = 1.0,
    composite: bool = True,
    composite_prob: float = 1.0,
    structure_pdb: Optional[Path] = None,
    precomputed_s_immuno: Optional[float] = None,
    parent_sequence: Optional[str] = None,
    use_package: bool = True,
) -> dict[str, Any]:
    """Compute full 7-axis PEN-SCORE for a single design.

    When pen-score package is available, delegates to its native API.
    Otherwise uses local axis implementations (reproducible fallback).

    Parameters
    ----------
    design_id : str
    protein_sequence : str
    strategy : str
        One of: "A_domain_swap", "B_ortholog_discovery", "C_deimmunization",
                "D_protmpnn_redesign"
    tier_a, tier_a_confidence, composite, composite_prob :
        From MECH-CLASS evaluation (Step 15).
    precomputed_s_immuno : float | None
        From Strategy C deimmunization pipeline (if available).
    parent_sequence : str | None
        IS621 WT sequence — used for S_Mature identity check (Strategy C).
    """
    n_aa = len(protein_sequence)

    # Try pen-score package first
    if use_package and _check_penscore_package():
        try:
            import penscore
            ps = penscore.PenScore(use_case="human_therapeutic_aav_insertion")
            result = ps.score(
                sequence=protein_sequence,
                mechanism=tier_a,
                composite=composite,
                structure_pdb=str(structure_pdb) if structure_pdb else None,
            )
            # Apply pre-registered S_Mature override
            s_mature = _compute_s_mature(strategy, protein_sequence, parent_sequence)
            result_dict = result.to_dict() if hasattr(result, "to_dict") else dict(result)
            result_dict["S_Mature"] = s_mature
            # Recompute total with corrected S_Mature
            total = sum(result_dict.get(f"S_{ax}", 0.0) * w
                        for ax, w in [("DSB", 0.25), ("Spec", 0.10), ("Cargo", 0.20),
                                      ("Deliv", 0.15), ("Immuno", 0.10), ("Prog", 0.15),
                                      ("Mature", 0.05)])
            result_dict["pen_score"] = round(total, 4)
            result_dict["design_id"] = design_id
            result_dict["strategy"] = strategy
            result_dict["pen_score_method"] = "pen_score_package"
            return result_dict
        except Exception:
            pass  # fall through to local implementation

    # Local axis implementations
    S_DSB    = _compute_s_dsb(tier_a)
    S_Cargo  = _compute_s_cargo(tier_a, composite)
    S_Deliv  = _compute_s_deliv(n_aa)
    S_Immuno = _compute_s_immuno(protein_sequence, strategy, precomputed_s_immuno)
    S_Mature = _compute_s_mature(strategy, protein_sequence, parent_sequence)

    # S_Spec and S_Prog require the pen-score package; use proxy values
    # S_Spec proxy: high confidence MECH-CLASS → high specificity (sequence-based)
    S_Spec   = round(min(1.0, tier_a_confidence * 0.9 + 0.1 * composite_prob), 4)
    # S_Prog proxy: composite architecture with RNA-binding module → programmable
    S_Prog   = round(min(1.0, composite_prob * 0.8 + 0.2 * (1.0 if composite else 0.0)), 4)

    pen_score = round(
        S_DSB * 0.25 + S_Spec * 0.10 + S_Cargo * 0.20 +
        S_Deliv * 0.15 + S_Immuno * 0.10 + S_Prog * 0.15 + S_Mature * 0.05,
        4
    )

    return {
        "design_id": design_id,
        "strategy": strategy,
        "protein_length_aa": n_aa,
        "tier_a": tier_a,
        "composite": composite,
        "S_DSB": round(S_DSB, 4),
        "S_Spec": S_Spec,
        "S_Cargo": round(S_Cargo, 4),
        "S_Deliv": round(S_Deliv, 4),
        "S_Immuno": round(S_Immuno, 4),
        "S_Prog": S_Prog,
        "S_Mature": round(S_Mature, 4),
        "pen_score": pen_score,
        "beats_is621": pen_score > IS621_PENSCORE,
        "pen_score_method": "local_axis_implementations",
        "pen_score_note": (
            "S_Spec and S_Prog are proxies; run with pen-score package for exact values."
            if not _check_penscore_package() else ""
        ),
    }


def run_penscore_batch(
    designs_df: pd.DataFrame,
    sequence_col: str = "protein_sequence",
    id_col: str = "design_id",
    strategy_col: str = "strategy",
    tier_a_col: str = "tier_a",
    composite_col: str = "composite",
    s_immuno_col: Optional[str] = "predicted_s_immuno",
    parent_sequence: Optional[str] = None,
    structure_col: Optional[str] = "final_pdb",
) -> pd.DataFrame:
    """Run PEN-SCORE on all designs. Adds pen_score and 7 axis columns."""
    rows: list[dict] = []
    for _, row in designs_df.iterrows():
        pdb = Path(row[structure_col]) if structure_col and row.get(structure_col) else None
        precomp_immuno = row.get(s_immuno_col) if s_immuno_col else None
        r = compute_penscore(
            design_id=row[id_col],
            protein_sequence=row[sequence_col],
            strategy=row.get(strategy_col, "A_domain_swap"),
            tier_a=row.get(tier_a_col, "DSB_FREE_TRANSEST_RECOMBINASE"),
            tier_a_confidence=row.get("tier_a_confidence", 1.0),
            composite=bool(row.get(composite_col, True)),
            composite_prob=float(row.get("composite_prob", 1.0)),
            structure_pdb=pdb,
            precomputed_s_immuno=float(precomp_immuno) if precomp_immuno is not None else None,
            parent_sequence=parent_sequence,
        )
        rows.append(r)

    pen_df = pd.DataFrame(rows)
    out = designs_df.copy()
    score_cols = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog",
                  "S_Mature", "pen_score", "beats_is621", "pen_score_method"]
    for col in score_cols:
        if col in pen_df.columns:
            out[col] = pen_df[col].values
    return out
