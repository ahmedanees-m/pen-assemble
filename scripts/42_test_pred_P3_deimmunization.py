"""
Pre-registered P3: Best Strategy C deimmunized IS621 variant achieves
S_Immuno > WT IS621 S_Immuno + 0.10 on human_therapeutic_aav_insertion.

IMPORTANT — execution plan had placeholder WT_IS621_S_IMMUNO = 0.250 (incorrect).
Corrected to 0.7594, which is IS621's actual S_Immuno from Paper 3 v0.1.0.

IS621 S_Immuno computation:
  9-mer sliding window, HLA-I alleles [HLA-A*02:01, A*01:01, B*07:02, B*44:02],
  n_binders = sum(presentation_score > 0.5), S_Immuno = clip(1 - n_binders/L / 0.35, 0, 1).
  IS621 published: 33 binders -> S_Immuno = 0.7594 (older MHCflurry version).
  IS621 with MHCflurry 2.2.1: 33 binders -> S_Immuno = 0.7243 (model update).

Both reference values are reported for transparency. Primary criterion uses
published 0.7594 (Paper 3 lockpoint, consistent with pre-registration).

Usage:
  py 42_test_pred_P3_deimmunization.py
"""
from __future__ import annotations
import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd

BASE    = Path(__file__).resolve().parent.parent / "pipeline_results_local_test"
OUT_DIR = BASE / "validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WT_IS621_S_IMMUNO_PUBLISHED = 0.7594   # Paper 3 v0.1.0 (pre-registered reference)
WT_IS621_S_IMMUNO_MHCFL221  = 0.7243   # MHCflurry 2.2.1 re-measurement
DELTA_THRESHOLD = 0.10


def main() -> None:
    triaged = pd.read_parquet(BASE / "part_d" / "triaged_designs.parquet")
    strat_c = triaged[triaged["strategy"] == "C"].copy()
    print(f"Strategy C designs in triaged catalog: {len(strat_c)}")

    if strat_c.empty:
        result = {
            "test_id": "P3_deimmunization",
            "verdict": "NOT_EVALUABLE",
            "reason": "No Strategy C designs in triaged catalog.",
        }
        (OUT_DIR / "P3_deimmunization_result.json").write_text(json.dumps(result, indent=2))
        print("P3 Verdict: NOT_EVALUABLE")
        return

    best_idx      = strat_c["S_Immuno"].idxmax()
    best_s_immuno = float(strat_c.loc[best_idx, "S_Immuno"])
    best_design   = strat_c.loc[best_idx, "design_id"]
    delta_primary = best_s_immuno - WT_IS621_S_IMMUNO_PUBLISHED
    delta_mhcfl   = best_s_immuno - WT_IS621_S_IMMUNO_MHCFL221
    pass_p3 = delta_primary >= DELTA_THRESHOLD

    print(f"\n=== P3 Test Results ===")
    print(f"IS621 S_Immuno (published, pre-registered reference): {WT_IS621_S_IMMUNO_PUBLISHED}")
    print(f"IS621 S_Immuno (MHCflurry 2.2.1 re-measurement):      {WT_IS621_S_IMMUNO_MHCFL221}")
    print(f"Delta threshold: >= {DELTA_THRESHOLD}")
    print(f"\nAll Strategy C designs:")
    for _, row in strat_c.sort_values("S_Immuno", ascending=False).iterrows():
        d = row["S_Immuno"] - WT_IS621_S_IMMUNO_PUBLISHED
        flag = "  <- PASS" if d >= DELTA_THRESHOLD else ""
        print(f"  {row['design_id'][:62]}: S_Immuno={row['S_Immuno']:.4f}  "
              f"delta={d:+.4f}{flag}")
    print(f"\nBest deimmunized design: {best_design}")
    print(f"  S_Immuno = {best_s_immuno:.4f}")
    print(f"  delta vs published IS621: {delta_primary:+.4f} (threshold >= {DELTA_THRESHOLD})")
    print(f"  delta vs MHCflurry 2.2.1 IS621: {delta_mhcfl:+.4f}")
    print(f"\nP3 Verdict: {'PASS' if pass_p3 else 'FAIL'}")

    result = {
        "test_id": "P3_deimmunization",
        "verdict": "PASS" if pass_p3 else "FAIL",
        "criterion": f"best Strategy C S_Immuno > IS621 S_Immuno + {DELTA_THRESHOLD}",
        "wt_is621_s_immuno_published": WT_IS621_S_IMMUNO_PUBLISHED,
        "wt_is621_s_immuno_mhcfl221":  WT_IS621_S_IMMUNO_MHCFL221,
        "best_deimm_design":  best_design,
        "best_deimm_s_immuno": best_s_immuno,
        "delta_vs_published": round(float(delta_primary), 4),
        "delta_vs_mhcfl221":  round(float(delta_mhcfl), 4),
        "threshold_delta": DELTA_THRESHOLD,
        "correction_note": (
            "Execution plan placeholder WT_IS621_S_IMMUNO=0.250 was incorrect. "
            "Corrected to Paper 3 published value 0.7594. "
            "IS621_deimmunized_v2 uses a 14-mutation greedy MHC-I+II anchor reduction strategy "
            "with cumDDG budget 15.0 kcal/mol."
        ),
        "all_strategy_c": strat_c[["design_id", "S_Immuno", "pen_score",
                                    "S_DSB", "S_Deliv"]].to_dict(orient="records"),
    }
    out = OUT_DIR / "P3_deimmunization_result.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"\nResult saved -> {out}")


if __name__ == "__main__":
    main()
