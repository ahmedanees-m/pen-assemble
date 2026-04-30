#!/bin/bash
# Boltz-1 Tier 2 wrapper - runs inside pen-stack/design:0.1.0
# Installs any missing deps, verifies imports, then launches runner.
set -e

echo "[boltz-wrapper] Start $(date)"
echo "[boltz-wrapper] CUDA: $(python3 -c 'import torch; print(torch.cuda.is_available(), torch.version.cuda)')"

# 1. Verify / install Python deps
echo "[boltz-wrapper] Checking deps..."
pip install -q biopython pandas pyarrow 2>&1 | tail -3

# 2. Verify Boltz-1 import
echo "[boltz-wrapper] Verifying Boltz-1..."
python3 -c "
import boltz
print('[boltz-wrapper] boltz version:', boltz.__version__)
import torch
print('[boltz-wrapper] CUDA devices:', torch.cuda.device_count())
# Verify boltz predict CLI is callable
import subprocess
r = subprocess.run(['boltz', 'predict', '--help'], capture_output=True)
print('[boltz-wrapper] boltz predict CLI rc:', r.returncode)
"

# 3. Verify PyRosetta (non-fatal - only needed for Step 13)
echo "[boltz-wrapper] Checking PyRosetta (non-fatal)..."
python3 -c "
import pyrosetta
print('[boltz-wrapper] PyRosetta OK:', pyrosetta.__version__)
" 2>&1 || echo "[boltz-wrapper] WARNING: PyRosetta not available - OK for Boltz-1 step, fix before Step 13"

# 4. Launch Boltz-1 runner
echo "[boltz-wrapper] Launching Boltz-1 runner..."
python3 /input/step12b_boltz1_runner.py

echo "[boltz-wrapper] All done $(date)"
