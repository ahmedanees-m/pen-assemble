"""
Pre-registered P2: >= 1 design achieves S_Cargo = 1.0 AND S_Deliv >= 0.9
on human_therapeutic_aav_insertion (megabase-cargo + single-AAV joint constraint).

IMPORTANT - this is an EXISTENCE CLAIM, not a per-design gate.
The pre-registration asks: does at least one design satisfy both criteria simultaneously?
Reporting "N/N pass" would misrepresent this as a gate test.

Context:
  S_Cargo = 1.0 for all IS110-family designs by mechanism (transesterification supports
  arbitrarily large cargo without size constraint). S_Deliv >= 0.9 is satisfied by
  compact IS110 orthologs (200-400 aa fit within AAV packaging limit).

Usage:
  py 41_test_pred_P2_cargo_deliv.py
"""
from __future__ import annotations
import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd

BASE    = Path(__file__).resolve().parent.parent / "pipeline_results_local_test"
OUT_DIR = BASE / "validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    triaged = pd.read_parquet(BASE / "part_d" / "triaged_designs.parquet")
    print(f"Triaged catalog: {len(triaged)} designs")

    mask = (triaged["S_Cargo"] == 1.0) & (triaged["S_Deliv"] >= 0.9)
    satisfying = triaged[mask].sort_values("pen_score", ascending=False)
    pass_p2 = len(satisfying) >= 1

    print(f"\n=== P2 Test Results (EXISTENCE CLAIM) ===")
    print(f"Criterion: S_Cargo == 1.0 AND S_Deliv >= 0.9")
    print(f"Designs satisfying both: {len(satisfying)}")
    print(f"\nS_Cargo distribution: {triaged['S_Cargo'].value_counts().to_dict()}")
    print(f"S_Deliv >= 0.9: {(triaged['S_Deliv'] >= 0.9).sum()} of {len(triaged)}")
    print(f"\nTop 10 by pen_score:")
    for _, row in satisfying.head(10).iterrows():
        print(f"  [{row['strategy']}] {row['design_id'][:52]}: "
              f"pen={row['pen_score']:.4f}  S_Cargo={row['S_Cargo']:.3f}  "
              f"S_Deliv={row['S_Deliv']:.4f}  ({row.get('protein_length_aa', '?')} aa)")
    print(f"\nP2 Verdict: {'PASS' if pass_p2 else 'FAIL'} (existence claim)")

    result = {
        "test_id": "P2_megabase_aav",
        "verdict": "PASS" if pass_p2 else "FAIL",
        "criterion": "S_Cargo == 1.0 AND S_Deliv >= 0.9 (joint constraint)",
        "threshold": ">= 1 design satisfies both",
        "n_designs_satisfying_joint_constraint": int(len(satisfying)),
        "interpretation": (
            "EXISTENCE CLAIM - not a per-design gate. "
            "IS110-family mechanism gives S_Cargo=1.0 to all designs by definition "
            "(transesterification imposes no cargo-size limit). "
            "S_Deliv >= 0.9 satisfied by compact IS110 orthologs (200-400 aa). "
            "Future experimental validation needed for delivery efficiency in vivo."
        ),
        "top5_designs": satisfying.head(5)[
            ["design_id", "strategy", "S_Cargo", "S_Deliv", "pen_score"]
        ].to_dict(orient="records"),
    }
    out = OUT_DIR / "P2_cargo_deliv_result.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"\nResult saved -> {out}")


if __name__ == "__main__":
    main()
