"""Strategy C VM runner — deploy this to the VM and execute directly.

Self-contained: no 'rich' dependency, hardcoded IS621 sequence as fallback.
Uses real netMHCpan-4.1 + NetMHCIIpan-4.0 already installed at ~/netmhc/.

Designed to be SFTPed to ~/strategy_C_runner.py and run as:
    nohup python3 ~/strategy_C_runner.py > ~/strategy_C.log 2>&1 &

Output directory: ~/strategy_C_results/
  deimmunized_variants.parquet       — zstd-compressed, all metadata
  deimmunized_variants.fasta         — sequences for structure prediction
  deimmunized_variants_manifest.json — summary + P3 pass status
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Set up sys.path to find pen_assemble (uploaded to ~/pen_assemble_pkg/)
# ---------------------------------------------------------------------------
_pkg_dir = Path.home() / "pen_assemble_pkg"
if _pkg_dir.exists() and str(_pkg_dir) not in sys.path:
    sys.path.insert(0, str(_pkg_dir))

# ---------------------------------------------------------------------------
# IS621 sequence — UniProt A0A2X3M8B0, 342 aa (Paper 3 baseline, locked)
# Cross-check: n_I=15, n_II=92, S_Immuno=0.7594 with Paper 3 MHC panel
# ---------------------------------------------------------------------------
IS621_SEQUENCE = (
    "MDRFFPVIRICKVGFTMEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAH"
    "ICIEATGTYMEPVAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERA"
    "LRALVVRHQALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGI"
    "GEKTSAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
    "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA"
)

assert len(IS621_SEQUENCE) == 342, f"IS621 sequence length mismatch: {len(IS621_SEQUENCE)}"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR = Path.home() / "strategy_C_results"
N_VARIANTS = 10
MAX_MUTATIONS = 15
MAX_DDG_TOTAL = 8.0
MAX_PER_MUT_DDG = 3.0
N_MC_STEPS = 2000
BASE_SEED = 42

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    print("=" * 70)
    print("Strategy C: Computational Deimmunization of IS621")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"IS621 length: {len(IS621_SEQUENCE)} aa")
    print(f"Output dir:   {OUTPUT_DIR}")
    print("=" * 70)

    # Import pen_assemble
    try:
        from pen_assemble.strategies.deimmunization import (
            run_deimmunization,
            _check_mhcpan_available,
            _check_mhciipan_available,
            IS621_S_IMMUNO_WT, P3_TARGET_MINIMUM,
        )
    except ImportError as e:
        print(f"ERROR importing pen_assemble: {e}")
        print(f"sys.path = {sys.path}")
        sys.exit(1)

    # Verify MHC tools
    mhc_i_ok = _check_mhcpan_available()
    mhc_ii_ok = _check_mhciipan_available()
    print(f"\nTool check:")
    print(f"  netMHCpan-4.1:    {'FOUND' if mhc_i_ok else 'MISSING'}")
    print(f"  NetMHCIIpan-4.0:  {'FOUND' if mhc_ii_ok else 'MISSING'}")
    if not mhc_i_ok:
        print("ERROR: netMHCpan-4.1 not found. Aborting.")
        sys.exit(1)

    print(f"\nWT IS621 baseline (Paper 3): S_Immuno = {IS621_S_IMMUNO_WT}")
    print(f"P3 target: S_Immuno >= {P3_TARGET_MINIMUM} (delta >= 0.10)")
    print(f"\nRunning {N_VARIANTS}x3 MC trajectories ({N_MC_STEPS} steps each)...")
    print("MC inner loop: HEURISTIC (fast)")
    print("WT + final variant eval: REAL netMHCpan/netMHCIIpan")
    print()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    variants = run_deimmunization(
        scaffold_id="IS621",
        scaffold_sequence=IS621_SEQUENCE,
        scaffold_pdb=None,           # no PDB on VM; uses sequence heuristic for surface
        max_mutations=MAX_MUTATIONS,
        max_ddg_total=MAX_DDG_TOTAL,
        max_per_mutation_ddg=MAX_PER_MUT_DDG,
        n_variants=N_VARIANTS,
        n_mc_steps=N_MC_STEPS,
        seed=BASE_SEED,
        output_dir=OUTPUT_DIR,
        use_netmhcpan=True,
    )

    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"DONE. Elapsed: {elapsed/60:.1f} min")
    print(f"Variants produced: {len(variants)}")
    if variants:
        n_p3 = sum(1 for v in variants if v.passes_p3())
        best = variants[0]
        print(f"P3 passing (delta >= 0.10): {n_p3}/{len(variants)}")
        print(f"Best: {best.variant_id}, S_Immuno={best.predicted_s_immuno:.4f}, "
              f"delta={best.predicted_s_immuno_delta:+.4f}, "
              f"mutations={best.total_mutations}")
        print(f"\nResults written to: {OUTPUT_DIR}")
        for f in sorted(OUTPUT_DIR.iterdir()):
            print(f"  {f.name}  ({f.stat().st_size/1024:.0f} KB)")

    # Write a completion marker so deploy script can poll for it
    (OUTPUT_DIR / "DONE.txt").write_text(
        json.dumps({
            "status": "complete",
            "elapsed_min": round(elapsed/60, 1),
            "n_variants": len(variants),
            "n_p3_pass": sum(1 for v in variants if v.passes_p3()) if variants else 0,
            "best_delta": variants[0].predicted_s_immuno_delta if variants else None,
        }, indent=2)
    )
    print("\n[DONE.txt written — deploy script will detect completion]")


if __name__ == "__main__":
    main()
