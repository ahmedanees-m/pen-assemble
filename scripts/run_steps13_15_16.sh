#!/bin/bash
# ============================================================
# PEN-ASSEMBLE: Steps 13 → 15 → 16 (ESMFold-based, no Docker)
# Run directly on VM after pipeline_data/ is uploaded.
#
# Usage:
#   bash run_steps13_15_16.sh [DATA_DIR] [RESULTS_DIR]
#
# Defaults:
#   DATA_DIR    = /home/anees_22phd0670/pen_pipeline_data
#   RESULTS_DIR = /home/anees_22phd0670/pen_pipeline_results
# ============================================================
set -euo pipefail

DATA_DIR="${1:-/home/anees_22phd0670/pen_pipeline_data}"
RESULTS_DIR="${2:-/home/anees_22phd0670/pen_pipeline_results}"
REPO="/home/anees_22phd0670/pen-stack/code/repos/pen-assemble"
PARENT_PDB="/home/anees_22phd0670/8WT6.pdb"
LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

export PYTHONPATH="${REPO}:${PYTHONPATH:-}"
echo "[pipeline] PYTHONPATH = ${PYTHONPATH}"
echo "[pipeline] DATA_DIR   = ${DATA_DIR}"
echo "[pipeline] RESULTS_DIR= ${RESULTS_DIR}"

# ── Verify parent PDB ────────────────────────────────────────────────────────
if [ ! -f "${PARENT_PDB}" ]; then
    echo "[pipeline] WARNING: Parent PDB not found at ${PARENT_PDB} — Rosetta will use Grantham proxy"
fi

# ── Check PyRosetta ──────────────────────────────────────────────────────────
echo "[pipeline] Checking PyRosetta..."
python3 -c "import pyrosetta; print('[pipeline] PyRosetta:', pyrosetta.__version__)" 2>&1 || \
    echo "[pipeline] WARNING: PyRosetta not available — will use RaSP/Grantham fallback"

# ── Check MECH-CLASS ─────────────────────────────────────────────────────────
echo "[pipeline] Checking MECH-CLASS..."
python3 -c "import mech_class; print('[pipeline] mech_class:', mech_class.__version__)" 2>&1 || {
    echo "[pipeline] ERROR: mech_class not found. Install: pip install mech-class>=0.5.1"
    exit 1
}

# ── Check pen_score ──────────────────────────────────────────────────────────
echo "[pipeline] Checking pen_score..."
python3 -c "import pen_score; print('[pipeline] pen_score:', pen_score.__version__)" 2>&1 || \
    echo "[pipeline] WARNING: pen_score not available — will use local axis implementations"

echo ""
echo "==========================================================="
echo " STEP 13: Rosetta ΔΔG Stability Filtering"
echo "==========================================================="

python3 "${REPO}/scripts/21_run_rosetta_ddg.py" \
    --structures_parquet "${DATA_DIR}/structures/all_designs_structures.parquet" \
    --designs_dir        "${DATA_DIR}/designs" \
    --parent_pdb         "${PARENT_PDB}" \
    --output_dir         "${RESULTS_DIR}/step13_stability" \
    --method             auto \
    2>&1 | tee "${LOG_DIR}/step13_$(date +%Y%m%d_%H%M%S).log"

echo "[pipeline] Step 13 done."

echo ""
echo "==========================================================="
echo " STEP 15: MECH-CLASS Re-evaluation"
echo "==========================================================="

python3 "${REPO}/scripts/23_run_mech_class.py" \
    --input_parquet  "${RESULTS_DIR}/step13_stability/stability_passed.parquet" \
    --designs_dir    "${DATA_DIR}/designs" \
    --output_dir     "${RESULTS_DIR}/step15_mechclass" \
    --min_confidence 0.80 \
    2>&1 | tee "${LOG_DIR}/step15_$(date +%Y%m%d_%H%M%S).log"

echo "[pipeline] Step 15 done."

echo ""
echo "==========================================================="
echo " STEP 16: PEN-SCORE 7-Axis Scoring"
echo "==========================================================="

python3 "${REPO}/scripts/24_run_pen_score.py" \
    --input_parquet "${RESULTS_DIR}/step15_mechclass/mech_class_passed.parquet" \
    --designs_dir   "${DATA_DIR}/designs" \
    --structures_parquet "${DATA_DIR}/structures/all_designs_structures.parquet" \
    --output_dir    "${RESULTS_DIR}/step16_penscore" \
    2>&1 | tee "${LOG_DIR}/step16_$(date +%Y%m%d_%H%M%S).log"

echo "[pipeline] Step 16 done."

echo ""
echo "==========================================================="
echo " PIPELINE COMPLETE"
echo "==========================================================="
echo "Results at: ${RESULTS_DIR}"
echo ""
echo "Key outputs:"
echo "  Step 13: ${RESULTS_DIR}/step13_stability/stability_passed.parquet"
echo "  Step 15: ${RESULTS_DIR}/step15_mechclass/mech_class_passed.parquet"
echo "  Step 16: ${RESULTS_DIR}/step16_penscore/p1_candidates.parquet"
echo "  Step 16: ${RESULTS_DIR}/step16_penscore/pen_score_summary.json"
