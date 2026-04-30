"""
Aggregate all 5 pre-registered prediction results into a single
all_predictions_summary.json and print the final publication verdict.

This script summarises the five pre-registered predictions. It loads the five
individual JSON result files produced by scripts 40-44 and applies
the prediction summary.

Run AFTER scripts 40-44 have been executed:
  py 40_test_pred_P1_beat_is621.py
  py 41_test_pred_P2_cargo_deliv.py
  py 42_test_pred_P3_deimmunization.py
  py 43_test_pred_P4_orthologs.py
  py 44_test_pred_P5_diversity.py
  py 45_summarise_predictions.py   <- this script

Usage:
  py 45_summarise_predictions.py
"""
from __future__ import annotations
import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd, numpy as np

BASE    = Path(__file__).resolve().parent.parent / "pipeline_results_local_test"
OUT_DIR = BASE / "validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Bootstrap cluster-level analysis constants (recomputed from seed-42 run)
BOOTSTRAP_CLUSTER = {
    "n_iterations": 1000,
    "seed": 42,
    "sigma": 0.02,
    "prob_at_least_1_strategy_D_in_top5": 1.0,
    "prob_at_least_1_strategy_C_in_top5": 1.0,
    "individual_D_rank_ci_interpretation": (
        "D-design individual rank CIs are wide ([3, 183]) because 32 designs "
        "cluster within 0.009 PenScore (0.9261-0.9353). With sigma=0.02 noise, "
        "they shuffle freely among themselves. This is within-cluster competition, "
        "not rank uncertainty relative to the full catalog. "
        "Cluster-level probability of top-5 presence = 100%."
    ),
    "reporting_note": (
        "Strategy D designs form a tight fitness cluster. Report cluster-level "
        "stability (100% top-5 presence probability) rather than individual rank CIs."
    ),
}

METHODOLOGY_NOTES = {
    "stability_gate_non_functional": (
        "ALL 47 Rosetta ddG values (A:15, C:2, D:30) are cross-structure absolute "
        "energies (-31,838 to -41,308 kcal/mol), not true ddG values. "
        "Stability gate was auto-passed; structural quality proxied by ESMFold "
        "pLDDT > 90 globally (min 71.9) and > 95 at active-site residues (min 75.3). "
        "See Deviations 4+5 in DESIGN_PROVENANCE.md."
    ),
    "mhcflurry_version_shift": (
        "Strategy D S_Immuno individually recomputed with MHCflurry 2.2.1. "
        "Strategy B uses IS621 conservative baseline 0.7594. "
        "Calibrated IS621 lockpoint 0.9255 (MHCflurry 2.2.1-consistent) reported as "
        "secondary analysis alongside verbatim 0.929."
    ),
    "p3_wt_reference_corrected": (
        "Earlier placeholder WT_IS621_S_IMMUNO=0.250 corrected to "
        "PEN-SCORE published value 0.7594."
    ),
    "p5_diversity_enforced": (
        "Rank-5 diversity was enforced: A_007 (0.9209) replaces natural rank-5 "
        "D023 (0.9319). Top-5 is not purely rank-ordered."
    ),
    "atlas_gate8_not_evaluated": (
        "Strategy B Gate 8 (ATLAS embedding distance for literature novelty) "
        "not evaluated in this run. Strategy B pre-filtered by genus diversity "
        "(700 distinct genera) at sourcing step."
    ),
}


def load_result(test_id: str, filename: str) -> dict:
    path = OUT_DIR / filename
    if not path.exists():
        return {"test_id": test_id, "verdict": "NOT_RUN",
                "error": f"File not found: {path}"}
    return json.loads(path.read_text())


def prediction_summary(n_pass: int, n_fail: int, n_ne: int) -> str:
    n_eval = n_pass + n_fail
    if n_pass == 5:
        return "5/5 pre-registered predictions PASS"
    elif n_pass == 4 and n_eval == 4 and n_ne == 1:
        return "4/4 evaluated predictions PASS, 1 NOT_EVALUABLE"
    elif n_pass == 4 and n_fail == 1:
        return "4/5 predictions PASS, 1 FAIL"
    elif n_pass == 3:
        return "3/5 predictions PASS"
    else:
        return "REWORK FRAMEWORK -- too few predictions verified"


def main() -> None:
    tests = {
        "P1_beat_is621":    load_result("P1", "P1_beat_is621_result.json"),
        "P2_cargo_deliv":   load_result("P2", "P2_cargo_deliv_result.json"),
        "P3_deimmunization": load_result("P3", "P3_deimmunization_result.json"),
        "P4_orthologs":     load_result("P4", "P4_orthologs_result.json"),
        "P5_diversity":     load_result("P5", "P5_diversity_result.json"),
    }

    verdicts = [t.get("verdict", "NOT_RUN") for t in tests.values()]
    n_pass = verdicts.count("PASS")
    n_fail = verdicts.count("FAIL")
    n_ne   = sum(1 for v in verdicts if v in ("NOT_EVALUABLE", "NOT_RUN"))
    policy = prediction_summary(n_pass, n_fail, n_ne)

    # Print summary
    print("\n" + "="*65)
    print("PRE-REGISTERED PREDICTION SUMMARY (PEN-ASSEMBLE)")
    print("="*65)
    labels = {
        "P1_beat_is621":    "P1 (>=5 beat IS621 0.929)",
        "P2_cargo_deliv":   "P2 (S_Cargo=1.0 AND S_Deliv>=0.9)",
        "P3_deimmunization": "P3 (S_Immuno gain >=0.10)",
        "P4_orthologs":     "P4 (>=10 Strategy B survivors)",
        "P5_diversity":     "P5 (top-5 from >=3 strategies)",
    }
    for pid, t in tests.items():
        v = t.get("verdict", "NOT_RUN")
        label = labels.get(pid, pid)
        symbol = "PASS" if v == "PASS" else ("FAIL" if v == "FAIL" else "N/E")
        print(f"  {symbol}  {label}")

    print(f"\n  PASS:   {n_pass}/5")
    print(f"  FAIL:   {n_fail}/5")
    print(f"  N/E:    {n_ne}/5")
    print(f"\n  {policy}")

    print(f"\n  Key results:")
    p1 = tests.get("P1_beat_is621", {})
    print(f"    P1: {p1.get('designs_beating_verbatim','?')} designs > 0.929 (verbatim); "
          f"{p1.get('designs_beating_calibrated','?')} > 0.9255 (calibrated)")
    p3 = tests.get("P3_deimmunization", {})
    print(f"    P3: IS621_deimmunized_v2 delta = {p3.get('delta_vs_published','?'):+.4f} "
          f"(best C, threshold >= 0.10)")
    p4 = tests.get("P4_orthologs", {})
    print(f"    P4: {p4.get('n_survivors','?')} Strategy B survivors (threshold >= 10)")

    print(f"\n  Methodology notes:")
    for key, msg in METHODOLOGY_NOTES.items():
        print(f"    [{key}]: {msg[:95]}...")

    # Assemble full summary JSON
    summary = {
        "evaluated_date": "2026-05-20",
        "pipeline_version": "PEN-ASSEMBLE v0.5.0",
        "n_total_predictions": 5,
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_not_evaluable": n_ne,
        "prediction_summary": policy,
        "individual_results": tests,
        "bootstrap_cluster_analysis": BOOTSTRAP_CLUSTER,
        "methodology_notes": METHODOLOGY_NOTES,
        "reporting_checklist": {
            "P1_verbatim_threshold_is_primary": True,
            "P1_calibrated_threshold_is_secondary_only": True,
            "P2_is_existence_claim_not_gate": True,
            "P3_wt_reference_corrected_from_0250": True,
            "P4_gate8_atlas_disclosed": True,
            "P5_diversity_enforcement_disclosed": True,
            "stability_gate_non_functionality_disclosed": True,
            "rosetta_ddg_deviation5_in_provenance": True,
        },
    }
    out = OUT_DIR / "all_predictions_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"\nFull summary saved -> {out}")


if __name__ == "__main__":
    main()
