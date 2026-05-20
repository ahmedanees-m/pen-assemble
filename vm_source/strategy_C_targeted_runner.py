"""Strategy C targeted deimmunization — anchor-position mutagenesis on the VM.

This replaces the heuristic-MC approach (which gave best delta=+0.0236) with a
proper targeted anchor-position scan that directly tests mutations at residues
known to anchor MHC-I and MHC-II binding.

Algorithm:
  1. Use known WT epitope positions (from our real netMHCpan run: n_I=15, n_II=92)
  2. For each MHC-I binder (15 positions):
     - Test conservative mutations at anchor residues (positions p+1 and p+8 in 9-mer)
     - Run netMHCpan on each modified sequence (one call per mutation candidate)
     - Score: n_binders_eliminated / ddg_cost
  3. Greedy selection: pick highest-scoring mutation, apply to current sequence,
     update remaining binder set, repeat for max_mutations rounds
  4. After MHC-I targeted phase: use remaining mutation budget for MHC-II reduction
  5. Final full netMHCpan + netMHCIIpan evaluation on assembled variant

Estimated runtime: ~200-400 netMHCpan calls x ~3s each = ~10-20 min on VM
Expected improvement: delta >= 0.06-0.12 vs WT (targeting P3 threshold of 0.10)

Deploy and run:
    python scripts/deploy_strategy_C_targeted_vm.py
OR directly on VM:
    python3 ~/strategy_C_targeted_runner.py

Output: ~/strategy_C_targeted_results/
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import subprocess
import time
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# sys.path setup
# ---------------------------------------------------------------------------
_pkg_dir = Path.home() / "pen_assemble_pkg"
if _pkg_dir.exists() and str(_pkg_dir) not in sys.path:
    sys.path.insert(0, str(_pkg_dir))

# ---------------------------------------------------------------------------
# Constants — IS621 (locked)
# ---------------------------------------------------------------------------
IS621_SEQUENCE = (
    "MDRFFPVIRICKVGFTMEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAH"
    "ICIEATGTYMEPVAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERA"
    "LRALVVRHQALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGI"
    "GEKTSAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
    "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA"
)
assert len(IS621_SEQUENCE) == 342

# WT epitope positions from real netMHCpan run (2026-05-16, Paper 3 panel)
# MHC-I: %Rank_EL < 0.5, alleles HLA-A02:01, HLA-A01:01, HLA-B07:02, HLA-B44:02
WT_MHC_I_POSITIONS = [15, 81, 86, 116, 139, 157, 176, 182, 200, 212, 221, 223, 267, 313, 317]
# MHC-II: %Rank_EL <= 10.0, alleles DRB1_0101, DRB1_0301, DRB1_0401
# We know n_II=92; positions will be recomputed in WT scan if needed
WT_N_II = 92

IS621_S_IMMUNO_WT = 0.7594
P3_TARGET_DELTA = 0.10
S_IMMUNO_MAX_TOTAL = 253.5

# Frozen positions (active site + bRNA interface — must not mutate)
IS621_FROZEN: frozenset[int] = frozenset([
    9, 10, 11, 12, 13,
    58, 59, 60, 61, 62,
    100, 101, 102, 103, 104, 105, 106,
    239, 240, 241, 242, 243,
    111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
    121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
])

# Conservative amino acid substitution set (avoid G/P)
_CONSERVATIVE_AA = list("ACDEFHIKLMNQRSTVWY")

# VM paths
HOME = Path.home()
BINARY_I  = str(HOME / "netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan")
PERL_II   = str(HOME / "netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl")
ALLELES_I  = "HLA-A02:01,HLA-A01:01,HLA-B07:02,HLA-B44:02"
ALLELES_II = "DRB1_0101,DRB1_0301,DRB1_0401"

OUTPUT_DIR = HOME / "strategy_C_targeted_results"
MAX_MUTATIONS = 15
MAX_DDG_TOTAL = 8.0
MAX_PER_MUT_DDG = 3.0


# ---------------------------------------------------------------------------
# MHC calling helpers (Paper 3 convention)
# ---------------------------------------------------------------------------

def _make_env() -> dict:
    env = os.environ.copy()
    env["NETMHCpan"] = str(HOME / "netmhc/netMHCpan-4.1/Linux_x86_64")
    env["NETMHCIIpan"] = str(HOME / "netmhc/netMHCIIpan-4.0")
    env["NetMHCIIpanPLAT"] = str(HOME / "netmhc/netMHCIIpan-4.0/Linux_x86_64")
    env["TMPDIR"] = "/tmp"
    return env


def run_netmhcpan_i(sequence: str) -> tuple[int, list[int]]:
    """Run Class I. Returns (n_unique_positions, [sorted_positions])."""
    env = _make_env()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, dir="/tmp") as f:
        f.write(f">query\n{sequence}\n")
        fpath = f.name
    try:
        r = subprocess.run(
            [BINARY_I, "-f", fpath, "-a", ALLELES_I, "-l", "9", "-BA"],
            capture_output=True, text=True, timeout=120, env=env,
        )
    finally:
        Path(fpath).unlink(missing_ok=True)
    positions: set[int] = set()
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"): continue
        parts = line.split()
        if len(parts) < 13 or not parts[0].lstrip("-").isdigit(): continue
        try:
            if float(parts[12]) < 0.5:
                positions.add(int(parts[0]))
        except (ValueError, IndexError): continue
    return len(positions), sorted(positions)


def run_netmhcpan_ii(sequence: str) -> tuple[int, list[int]]:
    """Run Class II. Returns (n_unique_positions, [sorted_positions])."""
    env = _make_env()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, dir="/tmp") as f:
        f.write(f">query\n{sequence}\n")
        fpath = f.name
    try:
        r = subprocess.run(
            ["perl", PERL_II, "-f", fpath, "-a", ALLELES_II],
            capture_output=True, text=True, timeout=180, env=env,
        )
    finally:
        Path(fpath).unlink(missing_ok=True)
    positions: set[int] = set()
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"): continue
        parts = line.split()
        if len(parts) < 9 or not parts[0].lstrip("-").isdigit(): continue
        try:
            if float(parts[8]) <= 10.0:
                positions.add(int(parts[0]))
        except (ValueError, IndexError): continue
    return len(positions), sorted(positions)


def compute_s_immuno(n_i: int, n_ii: int) -> float:
    return round(min(1.0, max(0.0, 1.0 - (n_i + 0.5 * n_ii) / S_IMMUNO_MAX_TOTAL)), 4)


# ---------------------------------------------------------------------------
# Grantham ddG proxy
# ---------------------------------------------------------------------------
_AA_MW = {
    "G": 57, "A": 71, "V": 99, "L": 113, "I": 113, "P": 97, "F": 147,
    "W": 186, "M": 131, "S": 87, "T": 101, "C": 103, "Y": 163, "H": 137,
    "D": 115, "E": 129, "N": 114, "Q": 128, "K": 128, "R": 156,
}
_POLAR = set("RKDENQHST")

def grantham_ddg(wt: str, mut: str) -> float:
    cross = (wt in _POLAR) != (mut in _POLAR)
    size_change = abs(_AA_MW.get(wt, 115) - _AA_MW.get(mut, 115)) / 30.0
    return round(min(3.0, 0.3 + (0.8 if cross else 0.0) + size_change * 0.5), 2)


# ---------------------------------------------------------------------------
# Targeted mutagenesis scan
# ---------------------------------------------------------------------------

def get_epitope_neighborhood(binder_positions: list[int], peptide_len: int = 9) -> set[int]:
    """All residue positions (1-indexed) within any binder's peptide window."""
    positions = set()
    for start in binder_positions:
        for offset in range(peptide_len):
            p = start + offset
            if 1 <= p <= len(IS621_SEQUENCE):
                positions.add(p)
    return positions


def get_mutable_positions(binder_positions: list[int], peptide_len: int = 9,
                          frozen: frozenset[int] = IS621_FROZEN,
                          anchor_only: bool = False) -> list[int]:
    """Get mutable surface positions within epitope windows.

    anchor_only=True: only anchor positions (2nd and 9th residue of each 9-mer)
    anchor_only=False: all positions within epitope windows
    """
    if anchor_only:
        positions = set()
        for start in binder_positions:
            anchor2 = start + 1  # 2nd residue (0-indexed offset 1)
            anchor9 = start + 8  # 9th residue (0-indexed offset 8)
            for p in (anchor2, anchor9):
                if 1 <= p <= len(IS621_SEQUENCE) and p not in frozen:
                    positions.add(p)
    else:
        positions = get_epitope_neighborhood(binder_positions, peptide_len) - frozen
    return sorted(p for p in positions if 1 <= p <= len(IS621_SEQUENCE))


def scan_single_mutations(
    current_seq: str,
    candidate_positions: list[int],
    current_i_positions: list[int],
    n_candidates_per_pos: int = 5,
) -> list[dict]:
    """For each candidate position, test top-N conservative substitutions with netMHCpan.

    Returns list of {pos, wt_aa, mut_aa, ddg, new_n_i, new_i_positions, binders_removed}.
    Sorted by binders_removed descending (best first).
    """
    seq_list = list(current_seq)
    results = []
    n_wt = len(current_i_positions)

    for pos in candidate_positions:
        wt_aa = seq_list[pos - 1]
        # Generate candidate substitutions ordered by estimated utility:
        # Prefer polar/charged replacements for MHC-anchor disruption
        # (anchor positions typically hydrophobic — swap to charged breaks binding)
        candidates = [
            aa for aa in _CONSERVATIVE_AA
            if aa != wt_aa and grantham_ddg(wt_aa, aa) <= MAX_PER_MUT_DDG
        ]
        # Sort: cross-group substitutions (more likely to break anchors) first
        candidates.sort(key=lambda aa: (
            0 if (wt_aa in _POLAR) != (aa in _POLAR) else 1,   # cross-group first
            grantham_ddg(wt_aa, aa)                             # then by ddg
        ))
        candidates = candidates[:n_candidates_per_pos]

        for mut_aa in candidates:
            ddg = grantham_ddg(wt_aa, mut_aa)
            proposed = seq_list.copy()
            proposed[pos - 1] = mut_aa
            proposed_seq = "".join(proposed)

            new_n_i, new_i_pos = run_netmhcpan_i(proposed_seq)
            binders_removed = n_wt - new_n_i

            results.append({
                "pos": pos,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "ddg": ddg,
                "new_n_i": new_n_i,
                "new_i_positions": new_i_pos,
                "binders_removed": binders_removed,
                # Score: epitopes removed per unit ddG (efficiency)
                "score": (binders_removed + 0.01) / (ddg + 0.1),
            })
            print(f"    Pos {pos} {wt_aa}->{mut_aa}: MHC-I={new_n_i} "
                  f"(removed={binders_removed}), ddG={ddg:.2f}")

    results.sort(key=lambda r: (-r["binders_removed"], r["ddg"]))
    return results


def greedy_deimmunize(
    max_mutations: int = MAX_MUTATIONS,
    max_ddg: float = MAX_DDG_TOTAL,
    anchor_only: bool = False,
    n_candidates_per_pos: int = 4,
) -> dict:
    """Greedy targeted deimmunization.

    Phase 1: Target MHC-I binders (known 15 positions)
    Phase 2: Use remaining mutation budget for MHC-II reduction
    """
    seq = IS621_SEQUENCE
    current_i_pos = list(WT_MHC_I_POSITIONS)
    mutations_applied: list[dict] = []
    cumulative_ddg = 0.0

    print(f"\n=== Phase 1: MHC-I targeted mutagenesis ===")
    print(f"  Starting positions: {len(current_i_pos)} binders, {len(mutations_applied)} mutations")
    print(f"  Anchor-only mode: {anchor_only}")

    for iteration in range(max_mutations):
        if not current_i_pos:
            print(f"  All MHC-I binders eliminated after {iteration} mutations.")
            break

        # Get mutable positions within current epitope windows
        candidate_positions = get_mutable_positions(
            current_i_pos, peptide_len=9, anchor_only=anchor_only
        )
        # Exclude already-mutated positions
        already_mutated = {m["position"] for m in mutations_applied}
        candidate_positions = [p for p in candidate_positions if p not in already_mutated]

        if not candidate_positions:
            print(f"  No more mutable positions in epitope windows. Stopping Phase 1.")
            break

        print(f"\n  Iteration {iteration+1}: {len(current_i_pos)} MHC-I binders, "
              f"testing {len(candidate_positions)} positions × {n_candidates_per_pos} substitutions "
              f"= {len(candidate_positions)*n_candidates_per_pos} netMHCpan calls...")

        scan_results = scan_single_mutations(
            seq, candidate_positions, current_i_pos,
            n_candidates_per_pos=n_candidates_per_pos,
        )

        if not scan_results:
            print(f"  No candidates found. Stopping Phase 1.")
            break

        # Pick best mutation that passes ddG budget
        best = None
        for r in scan_results:
            if cumulative_ddg + r["ddg"] <= max_ddg:
                best = r
                break

        if best is None or best["binders_removed"] <= 0:
            print(f"  Best mutation removes 0 binders. Stopping Phase 1.")
            break

        # Apply best mutation
        seq_list = list(seq)
        seq_list[best["pos"] - 1] = best["mut_aa"]
        seq = "".join(seq_list)
        current_i_pos = best["new_i_positions"]
        cumulative_ddg += best["ddg"]
        mutations_applied.append({
            "position": best["pos"],
            "wt_aa": best["wt_aa"],
            "mut_aa": best["mut_aa"],
            "ddg_pred": best["ddg"],
            "binders_removed": best["binders_removed"],
            "new_n_i": best["new_n_i"],
            "phase": "MHC-I",
        })

        print(f"  APPLIED: Pos {best['pos']} {best['wt_aa']}->{best['mut_aa']} "
              f"| MHC-I: {len(WT_MHC_I_POSITIONS)} -> {best['new_n_i']} "
              f"(-{best['binders_removed']}) | cumDDG={cumulative_ddg:.2f}")

        if len(mutations_applied) >= max_mutations:
            break

    # Phase 2: MHC-II targeted reduction with remaining budget
    remaining_mutations = max_mutations - len(mutations_applied)
    print(f"\n=== Phase 2: MHC-II targeted mutagenesis ===")
    print(f"  Remaining mutation budget: {remaining_mutations}")

    if remaining_mutations > 0 and cumulative_ddg < max_ddg:
        # Get current MHC-II binder positions
        print("  Running netMHCpan-II on current sequence to get binder positions...")
        n_ii_current, ii_positions = run_netmhcpan_ii(seq)
        print(f"  Current MHC-II: {n_ii_current} binders")

        if ii_positions:
            # Get mutable positions within MHC-II windows
            already_mutated = {m["position"] for m in mutations_applied}
            ii_candidates = get_mutable_positions(
                ii_positions[:30],  # Focus on top-30 highest-rank binders
                peptide_len=15, anchor_only=False
            )
            ii_candidates = [p for p in ii_candidates if p not in already_mutated]

            print(f"  Testing {min(len(ii_candidates), 20)} MHC-II epitope positions "
                  f"× {n_candidates_per_pos} substitutions...")

            for iteration in range(remaining_mutations):
                if not ii_positions or cumulative_ddg >= max_ddg:
                    break

                candidate_positions = get_mutable_positions(
                    ii_positions[:20], peptide_len=15, anchor_only=False
                )
                already_mutated = {m["position"] for m in mutations_applied}
                candidate_positions = [p for p in candidate_positions if p not in already_mutated]

                if not candidate_positions:
                    break

                # For MHC-II, use a proxy scan (test MHC-II directly)
                # Limit to top-10 positions to control call count
                candidate_positions = candidate_positions[:10]
                print(f"\n  MHC-II iteration {iteration+1}: testing {len(candidate_positions)} positions...")

                seq_list_current = list(seq)
                n_ii_wt = len(ii_positions)
                best_ii = None

                for pos in candidate_positions:
                    wt_aa = seq_list_current[pos - 1]
                    candidates = [aa for aa in _CONSERVATIVE_AA
                                  if aa != wt_aa and grantham_ddg(wt_aa, aa) <= MAX_PER_MUT_DDG]
                    candidates.sort(key=lambda aa: (
                        0 if (wt_aa in _POLAR) != (aa in _POLAR) else 1,
                        grantham_ddg(wt_aa, aa)
                    ))
                    candidates = candidates[:3]  # fewer for phase 2

                    for mut_aa in candidates:
                        ddg = grantham_ddg(wt_aa, mut_aa)
                        if cumulative_ddg + ddg > max_ddg:
                            continue
                        proposed = seq_list_current.copy()
                        proposed[pos - 1] = mut_aa
                        proposed_seq = "".join(proposed)
                        new_n_ii, new_ii_pos = run_netmhcpan_ii(proposed_seq)
                        removed = n_ii_wt - new_n_ii
                        print(f"    Pos {pos} {wt_aa}->{mut_aa}: MHC-II={new_n_ii} (removed={removed}), ddG={ddg:.2f}")
                        if best_ii is None or removed > best_ii["binders_removed"]:
                            best_ii = {"pos": pos, "wt_aa": wt_aa, "mut_aa": mut_aa,
                                       "ddg": ddg, "new_n_ii": new_n_ii,
                                       "new_ii_positions": new_ii_pos, "binders_removed": removed}

                if best_ii is None or best_ii["binders_removed"] <= 0:
                    print("  No MHC-II improvement found. Stopping Phase 2.")
                    break

                seq_list = list(seq)
                seq_list[best_ii["pos"] - 1] = best_ii["mut_aa"]
                seq = "".join(seq_list)
                ii_positions = best_ii["new_ii_positions"]
                cumulative_ddg += best_ii["ddg"]
                mutations_applied.append({
                    "position": best_ii["pos"],
                    "wt_aa": best_ii["wt_aa"],
                    "mut_aa": best_ii["mut_aa"],
                    "ddg_pred": best_ii["ddg"],
                    "binders_removed": best_ii["binders_removed"],
                    "new_n_ii": best_ii["new_n_ii"],
                    "phase": "MHC-II",
                })
                print(f"  APPLIED: Pos {best_ii['pos']} {best_ii['wt_aa']}->{best_ii['mut_aa']} "
                      f"| MHC-II: {n_ii_wt} -> {best_ii['new_n_ii']} "
                      f"(-{best_ii['binders_removed']}) | cumDDG={cumulative_ddg:.2f}")

    # Final evaluation
    print(f"\n=== Final evaluation (real netMHCpan) ===")
    print(f"  Sequence with {len(mutations_applied)} mutations, cumDDG={cumulative_ddg:.2f}")
    final_n_i, final_i_pos = run_netmhcpan_i(seq)
    final_n_ii, final_ii_pos = run_netmhcpan_ii(seq)
    final_combined = final_n_i + 0.5 * final_n_ii
    final_s_immuno = compute_s_immuno(final_n_i, final_n_ii)
    final_delta = round(final_s_immuno - IS621_S_IMMUNO_WT, 4)

    print(f"  Final MHC-I:  {len(WT_MHC_I_POSITIONS)} -> {final_n_i} "
          f"(-{len(WT_MHC_I_POSITIONS)-final_n_i})")
    print(f"  Final MHC-II: {WT_N_II} -> {final_n_ii} (-{WT_N_II-final_n_ii})")
    print(f"  Final combined: {61.0:.1f} -> {final_combined:.1f}")
    print(f"  Final S_Immuno: {IS621_S_IMMUNO_WT} -> {final_s_immuno:.4f} "
          f"(delta={final_delta:+.4f})")
    print(f"  P3 pass (delta>=0.10): {'YES' if final_delta >= 0.10 else 'NO'}")

    return {
        "variant_sequence": seq,
        "mutations": mutations_applied,
        "wt_n_i": len(WT_MHC_I_POSITIONS),
        "wt_n_ii": WT_N_II,
        "final_n_i": final_n_i,
        "final_n_ii": final_n_ii,
        "final_s_immuno": final_s_immuno,
        "final_delta": final_delta,
        "cumulative_ddg": round(cumulative_ddg, 2),
        "passes_p3": final_delta >= P3_TARGET_DELTA,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 70)
    print("Strategy C: Targeted Anchor-Position Deimmunization of IS621")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"WT: n_I={len(WT_MHC_I_POSITIONS)}, n_II={WT_N_II}, S_Immuno={IS621_S_IMMUNO_WT}")
    print(f"P3 target: delta >= {P3_TARGET_DELTA} (S_Immuno >= {IS621_S_IMMUNO_WT+P3_TARGET_DELTA:.4f})")
    print("=" * 70)

    # Verify tools
    if not Path(BINARY_I).exists():
        print(f"ERROR: netMHCpan not found at {BINARY_I}")
        sys.exit(1)
    if not Path(PERL_II).exists():
        print(f"ERROR: NetMHCIIpan not found at {PERL_II}")
        sys.exit(1)
    print(f"netMHCpan-4.1: FOUND")
    print(f"NetMHCIIpan-4.0: FOUND")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Run targeted greedy optimization.
    # anchor_only=True: only tests positions p+1 (anchor 2) and p+8 (anchor 9) of each 9-mer.
    # This reduces calls from ~320/iteration to ~30/iteration → ~750 total (~40 min on VM).
    result = greedy_deimmunize(
        max_mutations=MAX_MUTATIONS,
        max_ddg=MAX_DDG_TOTAL,
        anchor_only=True,     # anchor positions only (fast: ~30 positions per iteration)
        n_candidates_per_pos=3,  # top-3 cross-group substitutions per anchor position
    )

    # Save outputs
    import pandas as pd

    rows = [{
        "variant_id": "C_targeted_001",
        "protein_sequence": result["variant_sequence"],
        "mutations_introduced": result["mutations"],
        "total_mutations": len(result["mutations"]),
        "n_mhc_i_epitopes_removed": result["wt_n_i"] - result["final_n_i"],
        "n_mhc_ii_epitopes_removed": result["wt_n_ii"] - result["final_n_ii"],
        "predicted_s_immuno": result["final_s_immuno"],
        "predicted_s_immuno_delta": result["final_delta"],
        "predicted_ddg_total": result["cumulative_ddg"],
        "passes_p3": result["passes_p3"],
        "strategy": "C_targeted_deimmunization",
        "scaffold_id": "IS621",
        "s_immuno_wt_baseline": IS621_S_IMMUNO_WT,
        "method": "targeted_anchor_position_greedy",
    }]
    df = pd.DataFrame(rows)

    parquet_path = OUTPUT_DIR / "deimmunized_targeted.parquet"
    df.to_parquet(str(parquet_path), index=False, compression="zstd")

    fasta_path = OUTPUT_DIR / "deimmunized_targeted.fasta"
    with fasta_path.open("w") as f:
        f.write(
            f">C_targeted_001 IS621_deimm_targeted "
            f"delta_s_immuno={result['final_delta']:+.4f} "
            f"n_mut={len(result['mutations'])} "
            f"ddg={result['cumulative_ddg']:.2f}\n"
            f"{result['variant_sequence']}\n"
        )

    manifest = {
        "method": "targeted_anchor_position_greedy",
        "n_variants": 1,
        "n_p3_pass": 1 if result["passes_p3"] else 0,
        "wt_n_i": result["wt_n_i"],
        "wt_n_ii": result["wt_n_ii"],
        "final_n_i": result["final_n_i"],
        "final_n_ii": result["final_n_ii"],
        "final_s_immuno": result["final_s_immuno"],
        "final_delta": result["final_delta"],
        "cumulative_ddg": result["cumulative_ddg"],
        "total_mutations": len(result["mutations"]),
        "passes_p3": result["passes_p3"],
        "elapsed_min": round((time.time() - t0) / 60, 1),
        "netmhcpan_used": True,
        "netmhciipan_used": True,
        "mhc_i_rank_threshold": 0.5,
        "mhc_ii_rank_threshold": 10.0,
        "s_immuno_max_total": S_IMMUNO_MAX_TOTAL,
        "mutations": result["mutations"],
    }
    (OUTPUT_DIR / "deimmunized_targeted_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )

    (OUTPUT_DIR / "DONE.txt").write_text(json.dumps({
        "status": "complete",
        "elapsed_min": manifest["elapsed_min"],
        "final_delta": result["final_delta"],
        "passes_p3": result["passes_p3"],
        "total_mutations": len(result["mutations"]),
    }, indent=2))

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"DONE. Elapsed: {elapsed/60:.1f} min")
    print(f"P3 pass: {'YES' if result['passes_p3'] else 'NO'} (delta={result['final_delta']:+.4f})")
    print(f"Results in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
