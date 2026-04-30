#!/bin/bash
# ESMFold Tier 1 wrapper - runs inside esm-fold:v1 with --gpus all
# Uses the exact openfold commit pinned by fair-esm README:
#   https://github.com/aqlaboratory/openfold@4b41059694619831a7db195b7e0988fc4ff3a307
set -e

echo "[wrapper] Start $(date)"
echo "[wrapper] torch CUDA: $(python3 -c 'import torch; print(torch.cuda.is_available(), torch.cuda.device_count())')"

OPENFOLD_DIR="/opt/openfold"
OPENFOLD_COMMIT="4b41059694619831a7db195b7e0988fc4ff3a307"
SO_PATH="$OPENFOLD_DIR/attn_core_inplace_cuda.cpython-310-x86_64-linux-gnu.so"

# 1. Python deps
echo "[wrapper] Installing Python deps..."
pip install -q \
    omegaconf einops biotite \
    ml-collections dm-tree \
    biopython pandas scipy networkx \
    modelcif
echo "[wrapper] Deps installed."

# 2. Clone openfold at the pinned fair-esm compatible commit
# Remove any previously cloned version that is on the wrong commit
if [ -d "$OPENFOLD_DIR/openfold" ]; then
    CURRENT=$(git -C "$OPENFOLD_DIR" rev-parse HEAD 2>/dev/null || echo "unknown")
    if [ "$CURRENT" != "$OPENFOLD_COMMIT" ]; then
        echo "[wrapper] Wrong openfold commit ($CURRENT), removing..."
        rm -rf "$OPENFOLD_DIR"
    else
        echo "[wrapper] Openfold at correct commit $OPENFOLD_COMMIT."
    fi
fi

if [ ! -d "$OPENFOLD_DIR/openfold" ]; then
    echo "[wrapper] Cloning openfold at pinned commit $OPENFOLD_COMMIT..."
    git clone -q https://github.com/aqlaboratory/openfold.git "$OPENFOLD_DIR"
    git -C "$OPENFOLD_DIR" checkout -q "$OPENFOLD_COMMIT"
    echo "[wrapper] Cloned and checked out."
fi

# 3. Compile CUDA extension if not yet built
if [ ! -f "$SO_PATH" ]; then
    echo "[wrapper] Compiling openfold CUDA extensions (~10-15 min)..."
    export PATH="/usr/local/cuda/bin:$PATH"
    cd "$OPENFOLD_DIR"
    python3 setup.py build_ext --inplace 2>&1
    echo "[wrapper] CUDA extensions compiled."
    cd /
else
    echo "[wrapper] CUDA extension already built - skipping."
fi

export PYTHONPATH="$OPENFOLD_DIR:$PYTHONPATH"

# 4. Verify full ESMFold import
echo "[wrapper] Verifying ESMFold imports..."
python3 -c "
import sys
sys.path.insert(0, '/opt/openfold')
import openfold.data.data_transforms;       print('  data_transforms OK')
import openfold.model.primitives;           print('  primitives OK')
import openfold.model.structure_module;     print('  structure_module OK')
import openfold.model.triangular_attention; print('  triangular_attention OK')
import openfold.utils.loss;                print('  loss OK')
import esm.esmfold.v1.pretrained;          print('  ESMFold v1 pretrained OK')
print('[wrapper] All imports verified.')
"

# 5. ESMFold batch inference
echo "[wrapper] Launching ESMFold runner..."
PYTHONPATH="$OPENFOLD_DIR" python3 /input/step12_esmfold_runner.py

echo "[wrapper] All done $(date)"
