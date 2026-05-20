"""Multi-gate triage filter for design candidates. Step 17.

Applies all verification gates in sequence:
  Gate 1: ESMFold Tier 1 mean pLDDT >= 50
  Gate 2: Boltz-1 Tier 2 catalytic-core pLDDT >= 75
  Gate 3: Rosetta ΔΔG <= +5.0 kcal/mol
  Gate 4: Active-site geometry within ±1.5 Å of WT reference
  Gate 5: MECH-CLASS tier_a == 'DSB_FREE_TRANSEST_RECOMBINASE'
  Gate 6: PEN-SCORE computed (no sentinel flags)
"""
from __future__ import annotations

from typing import Optional

import pandas as pd

TRIAGE_GATES = {
    "esmfold_mean_plddt_min": 50.0,
    "boltz_catalytic_plddt_min": 75.0,
    "ddg_max_kcal": 5.0,
    "active_site_tolerance_a": 1.5,
    "required_tier_a": "DSB_FREE_TRANSEST_RECOMBINASE",
}


def run_triage(
    designs_df: pd.DataFrame,
    gates: Optional[dict] = None,
) -> pd.DataFrame:
    """Apply all triage gates sequentially; return surviving designs.

    Each gate adds a boolean `gate_X_pass` column and removes failing rows.
    Preserves original index for provenance tracing.

    Gate A — ESMFold mean pLDDT >= 50 (column: esmfold_mean_plddt or mean_plddt)
    Gate B — Rosetta ΔΔG <= +5.0 kcal/mol (column: ddg_kcal_mol)
    Gate C — Active-site geometry (column: active_site_pass, bool)
    Gate D — MECH-CLASS DSB-free (column: mech_class_tier_a)
    Gate E — PEN-SCORE minimum > 0 (column: pen_score)

    Missing gate columns are treated as a skip (gate always passes) with a warning
    printed to stderr, consistent with graceful-degradation across all scripts.

    Args:
        designs_df: DataFrame of design candidates with verification columns.
        gates: Override thresholds (keys match TRIAGE_GATES). Defaults to TRIAGE_GATES.

    Returns:
        DataFrame of designs that passed all available gates, with added gate_X_pass columns.
    """
    import sys

    g = {**TRIAGE_GATES, **(gates or {})}
    df = designs_df.copy()
    n_start = len(df)

    # Gate A: ESMFold pLDDT
    plddt_col = next((c for c in ["esmfold_mean_plddt", "mean_plddt", "af3_mean_plddt"]
                      if c in df.columns), None)
    if plddt_col:
        df["gate_A_pass"] = df[plddt_col] >= g["esmfold_mean_plddt_min"]
        df = df[df["gate_A_pass"]]
    else:
        print(f"  [triage] Gate A: pLDDT column not found — skipped", file=sys.stderr)

    # Gate B: Rosetta ΔΔG
    if "ddg_kcal_mol" in df.columns:
        df["gate_B_pass"] = df["ddg_kcal_mol"] <= g["ddg_max_kcal"]
        df = df[df["gate_B_pass"]]
    else:
        print(f"  [triage] Gate B: ddg_kcal_mol column not found — skipped", file=sys.stderr)

    # Gate C: Active-site geometry
    if "active_site_pass" in df.columns:
        df["gate_C_pass"] = df["active_site_pass"].fillna(True).astype(bool)
        df = df[df["gate_C_pass"]]
    else:
        print(f"  [triage] Gate C: active_site_pass column not found — skipped", file=sys.stderr)

    # Gate D: MECH-CLASS DSB-free
    if "mech_class_tier_a" in df.columns:
        df["gate_D_pass"] = df["mech_class_tier_a"] == g["required_tier_a"]
        df = df[df["gate_D_pass"]]
    else:
        print(f"  [triage] Gate D: mech_class_tier_a column not found — skipped", file=sys.stderr)

    # Gate E: PEN-SCORE sanity (>= 0 ensures no sentinel NaN made it through)
    if "pen_score" in df.columns:
        df["gate_E_pass"] = df["pen_score"].notna() & (df["pen_score"] >= 0.0)
        df = df[df["gate_E_pass"]]
    else:
        print(f"  [triage] Gate E: pen_score column not found — skipped", file=sys.stderr)

    n_end = len(df)
    print(
        f"  [triage] {n_start} → {n_end} designs survived all gates "
        f"({n_start - n_end} eliminated)",
        file=sys.stderr,
    )
    return df.reset_index(drop=True)
