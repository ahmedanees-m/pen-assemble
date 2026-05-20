"""Strategy C v2 — MHC-II-first anchor-targeted deimmunization.

Addresses the two failures of 12c (delta=+0.0315, P3 not achieved):

  FAILURE 1: DDG budget exhausted before Phase 2.
    12c used Phase 1 (MHC-I) first, consuming 6.91/8.0 DDG, leaving only
    1.09 for Phase 2 — so ALL cross-group mutations in Phase 2 were skipped
    (minimum cross-group ddg = 1.1). Fix: run Phase 2 (MHC-II) FIRST with
    dedicated budget of 10.0, then Phase 1 (MHC-I) with remaining budget.
    Total budget raised to 15.0.

  FAILURE 2: Wrong positions tested for MHC-II.
    12c's Phase 2 took all positions within 15-mer windows and tested the
    first 10 in sorted (N-terminal) order — NOT the actual MHC-II anchor
    positions. Fix: parse NetMHCIIpan output 'Of' field to compute actual
    core anchors P1/P4/P6/P9 (DR allele anchor map: offset 0, 3, 5, 8 from
    core start). Target those specific residues with cross-group substitutions.

Algorithm:
  Phase 2 (MHC-II, FIRST):
    1. Run NetMHCIIpan on WT IS621 → parse (start, Of, core, rank_el)
    2. For each binder, compute P1=start+Of, P4=P1+3, P6=P1+5, P9=P1+8
    3. Greedy scan: test top-5 cross-group substitutions at each anchor
    4. Apply best mutation (most binders removed / lowest ddg)
    5. Repeat for min(9, max_mutations) rounds with budget 10.0

  Phase 1 (MHC-I, SECOND):
    1. Apply anchor-position scan (p+1, p+8 of each 9-mer) — same as 12c
    2. Budget = max_ddg - phase2_cumddg
    3. Repeat for remaining mutation count

  Final: real netMHCpan + netMHCIIpan full evaluation

Mathematical target:
  WT: n_I=15, n_II=92, combined=61.0, S_Immuno=0.7594
  P3: combined <= 35.65 (delta >= 0.10)
  Need: reduce combined by 25.35
  Scenario: n_I=3 (reduce by 12, 8 mutations), n_II=65 (reduce by 27, 7 mutations)
            combined = 3 + 32.5 = 35.5 -> delta = 0.1001 -> P3 PASS

Deploy:
    python scripts/deploy_strategy_C_v2.py
OR directly on VM:
    python3 ~/strategy_C_v2_runner.py

Output: ~/strategy_C_v2_results/
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

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

WT_MHC_I_POSITIONS = [15, 81, 86, 116, 139, 157, 176, 182, 200, 212, 221, 223, 267, 313, 317]
WT_N_I = 15
WT_N_II = 92
IS621_S_IMMUNO_WT = 0.7594
P3_TARGET_DELTA = 0.10
S_IMMUNO_MAX_TOTAL = 253.5

# Frozen positions (active site + bRNA interface)
IS621_FROZEN: frozenset[int] = frozenset([
    9, 10, 11, 12, 13,
    58, 59, 60, 61, 62,
    100, 101, 102, 103, 104, 105, 106,
    239, 240, 241, 242, 243,
    111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
    121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
    145, 146, 147, 148, 149, 150,
])

# Conservative AA set (no G, P which destabilize backbone)
_CONSERVATIVE_AA = list("ACDEFHIKLMNQRSTVWY")
_POLAR = set("RKDENQHST")

# VM paths
HOME = Path.home()
BINARY_I  = str(HOME / "netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan")
PERL_II   = str(HOME / "netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl")
ALLELES_I  = "HLA-A02:01,HLA-A01:01,HLA-B07:02,HLA-B44:02"
ALLELES_II = "DRB1_0101,DRB1_0301,DRB1_0401"

OUTPUT_DIR = HOME / "strategy_C_v2_results"

# Budget parameters (raised from 8.0 to address Phase 2 budget exhaustion)
MAX_MUTATIONS     = 15
MAX_DDG_TOTAL     = 15.0   # raised from 8.0
MAX_PER_MUT_DDG   = 4.0    # raised from 3.0
PHASE2_DDG_BUDGET = 10.0   # MHC-II gets priority budget
# Phase 1 (MHC-I) gets remainder: MAX_DDG_TOTAL - cumddg_after_phase2

# MHC-II DR anchor offsets from core P1 position (0-based within core)
# P1=0, P4=3, P6=5, P9=8 — key anchor positions for DRB1 alleles
MHC_II_ANCHOR_OFFSETS = [0, 3, 5, 8]


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
def _make_env() -> dict:
    env = os.environ.copy()
    env["NETMHCpan"]       = str(HOME / "netmhc/netMHCpan-4.1/Linux_x86_64")
    env["NETMHCIIpan"]     = str(HOME / "netmhc/netMHCIIpan-4.0")
    env["NetMHCIIpanPLAT"] = str(HOME / "netmhc/netMHCIIpan-4.0/Linux_x86_64")
    env["TMPDIR"]          = "/tmp"
    return env


# ---------------------------------------------------------------------------
# MHC calling
# ---------------------------------------------------------------------------
def run_mhcpan_i(sequence: str) -> tuple[int, list[int]]:
    """Class I — returns (n_unique_positions, sorted_positions)."""
    env = _make_env()
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, dir="/tmp") as f:
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
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        if len(parts) < 13 or not parts[0].lstrip("-").isdigit():
            continue
        try:
            if float(parts[12]) < 0.5:
                positions.add(int(parts[0]))
        except (ValueError, IndexError):
            continue
    return len(positions), sorted(positions)


def run_mhcpan_ii_full(sequence: str) -> tuple[int, list[dict]]:
    """Class II — returns (n_unique_positions, binder_records).

    binder_records: list of {pos, offset, core, allele, rank_el}
    offset = 'Of' field from NetMHCIIpan output (core start offset in 15-mer)
    P1_abs = pos + offset  (1-indexed absolute position in sequence)
    """
    env = _make_env()
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, dir="/tmp") as f:
        f.write(f">query\n{sequence}\n")
        fpath = f.name
    try:
        r = subprocess.run(
            ["perl", PERL_II, "-f", fpath, "-a", ALLELES_II],
            capture_output=True, text=True, timeout=300, env=env,
        )
    finally:
        Path(fpath).unlink(missing_ok=True)

    seen_pos: set[int] = set()
    binders: list[dict] = []

    for line in r.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        # NetMHCIIpan-4.0 columns (verified):
        # [0]=Pos [1]=MHC [2]=Peptide(15mer) [3]=Of [4]=Core(9mer)
        # [5]=Core_Rel [6]=Identity [7]=Score_EL [8]=%Rank_EL [9]=Exp_Bind
        if len(parts) < 9 or not parts[0].lstrip("-").isdigit():
            continue
        try:
            pos      = int(parts[0])
            allele   = parts[1]
            offset   = int(parts[3])   # Of: offset of 9-mer core within 15-mer
            core     = parts[4]        # 9-mer core sequence
            rank_el  = float(parts[8])
        except (ValueError, IndexError):
            continue
        if rank_el <= 10.0:
            seen_pos.add(pos)
            binders.append({
                "pos":    pos,
                "offset": offset,
                "core":   core,
                "allele": allele,
                "rank_el": rank_el,
                # Absolute 1-indexed position of core P1 in the full sequence
                "p1_abs": pos + offset,
            })

    # Sort by rank_el ascending (best binders first — most important to target)
    binders.sort(key=lambda b: b["rank_el"])
    return len(seen_pos), binders


def run_mhcpan_ii_count(sequence: str) -> tuple[int, list[int]]:
    """Class II — returns (n_unique_positions, sorted_positions). Lightweight."""
    env = _make_env()
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, dir="/tmp") as f:
        f.write(f">query\n{sequence}\n")
        fpath = f.name
    try:
        r = subprocess.run(
            ["perl", PERL_II, "-f", fpath, "-a", ALLELES_II],
            capture_output=True, text=True, timeout=300, env=env,
        )
    finally:
        Path(fpath).unlink(missing_ok=True)
    positions: set[int] = set()
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        if len(parts) < 9 or not parts[0].lstrip("-").isdigit():
            continue
        try:
            if float(parts[8]) <= 10.0:
                positions.add(int(parts[0]))
        except (ValueError, IndexError):
            continue
    return len(positions), sorted(positions)


# ---------------------------------------------------------------------------
# Grantham ddG proxy
# ---------------------------------------------------------------------------
_AA_MW = {
    "G": 57,  "A": 71,  "V": 99,  "L": 113, "I": 113, "P": 97,  "F": 147,
    "W": 186, "M": 131, "S": 87,  "T": 101, "C": 103, "Y": 163, "H": 137,
    "D": 115, "E": 129, "N": 114, "Q": 128, "K": 128, "R": 156,
}

def grantham_ddg(wt: str, mut: str) -> float:
    cross = (wt in _POLAR) != (mut in _POLAR)
    size_change = abs(_AA_MW.get(wt, 115) - _AA_MW.get(mut, 115)) / 30.0
    return round(min(MAX_PER_MUT_DDG, 0.3 + (0.8 if cross else 0.0) + size_change * 0.5), 2)


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------
def get_mhc_ii_anchor_positions(binders: list[dict],
                                  frozen: frozenset[int] = IS621_FROZEN,
                                  max_binders: int = 30) -> list[int]:
    """Get P1/P4/P6/P9 anchor positions from NetMHCIIpan binder records.

    For DRB1 alleles, the 9-mer core anchors are at offsets 0,3,5,8 within
    the core. Given core starts at absolute position p1_abs (= pos + offset):
      P1 = p1_abs + 0
      P4 = p1_abs + 3
      P6 = p1_abs + 5
      P9 = p1_abs + 8

    Returns sorted unique list of non-frozen mutable anchor positions.
    """
    anchor_positions: set[int] = set()
    seen_binders = set()  # deduplicate by (pos, allele)
    n_used = 0
    for b in binders:
        key = (b["pos"], b["allele"])
        if key in seen_binders:
            continue
        seen_binders.add(key)
        n_used += 1
        if n_used > max_binders:
            break
        p1 = b["p1_abs"]
        for off in MHC_II_ANCHOR_OFFSETS:
            p = p1 + off
            if 1 <= p <= len(IS621_SEQUENCE) and p not in frozen:
                anchor_positions.add(p)
    return sorted(anchor_positions)


def get_mhc_i_anchor_positions(binder_positions: list[int],
                                 frozen: frozenset[int] = IS621_FROZEN) -> list[int]:
    """MHC-I anchor positions: p+1 (anchor 2) and p+8 (anchor 9) of each 9-mer."""
    positions: set[int] = set()
    for start in binder_positions:
        for p in (start + 1, start + 8):
            if 1 <= p <= len(IS621_SEQUENCE) and p not in frozen:
                positions.add(p)
    return sorted(positions)


def substitution_candidates(wt_aa: str, n: int = 5,
                              budget: float = MAX_PER_MUT_DDG) -> list[tuple[str, float]]:
    """Top-N substitutions ordered by: cross-group first, then by ddg."""
    candidates = [
        (aa, grantham_ddg(wt_aa, aa))
        for aa in _CONSERVATIVE_AA
        if aa != wt_aa and grantham_ddg(wt_aa, aa) <= budget
    ]
    candidates.sort(key=lambda x: (
        0 if (wt_aa in _POLAR) != (x[0] in _POLAR) else 1,
        x[1]
    ))
    return candidates[:n]


# ---------------------------------------------------------------------------
# Phase 2: MHC-II core-anchor scan
# ---------------------------------------------------------------------------
def phase2_mhcii(
    seq: str,
    n_rounds: int,
    ddg_budget: float,
    mutations_applied: list[dict],
    n_candidates: int = 5,
) -> tuple[str, float, list[dict]]:
    """Greedy MHC-II deimmunization targeting core anchor positions P1/P4/P6/P9.

    Returns (updated_seq, cumulative_ddg, mutations_applied).
    """
    cumddg = sum(m["ddg_pred"] for m in mutations_applied)
    already_mutated = {m["position"] for m in mutations_applied}

    print(f"\n=== Phase 2: MHC-II core-anchor targeted mutagenesis ===")
    print(f"  Budget: {ddg_budget:.1f} DDG, {n_rounds} rounds max")
    print(f"  Getting WT MHC-II binder positions with core offsets...")

    n_ii_wt, binders_wt = run_mhcpan_ii_full(seq)
    print(f"  MHC-II binders: {n_ii_wt} (WT)")

    current_n_ii = n_ii_wt
    binders_current = binders_wt

    for rnd in range(n_rounds):
        if not binders_current:
            print(f"  All MHC-II binders eliminated after {rnd} rounds.")
            break
        if cumddg >= ddg_budget:
            print(f"  DDG budget exhausted ({cumddg:.2f}/{ddg_budget:.1f}). Stopping Phase 2.")
            break

        # Get anchor positions from current binder set
        anchor_pos = get_mhc_ii_anchor_positions(binders_current, max_binders=40)
        anchor_pos = [p for p in anchor_pos if p not in already_mutated]

        if not anchor_pos:
            print(f"  No more mutable MHC-II anchor positions. Stopping Phase 2.")
            break

        remaining_budget_per_mut = ddg_budget - cumddg
        seq_list = list(seq)
        best = None

        print(f"\n  Round {rnd+1}: {current_n_ii} binders, "
              f"{len(anchor_pos)} anchor positions, "
              f"budget_remaining={remaining_budget_per_mut:.2f}")

        for pos in anchor_pos:
            wt_aa = seq_list[pos - 1]
            candidates = substitution_candidates(
                wt_aa, n=n_candidates, budget=min(MAX_PER_MUT_DDG, remaining_budget_per_mut)
            )
            for mut_aa, ddg in candidates:
                proposed = seq_list.copy()
                proposed[pos - 1] = mut_aa
                proposed_seq = "".join(proposed)

                new_n_ii, new_pos = run_mhcpan_ii_count(proposed_seq)
                removed = current_n_ii - new_n_ii
                score = (removed + 0.01) / (ddg + 0.1)

                print(f"    Pos {pos} {wt_aa}->{mut_aa}: MHC-II={new_n_ii} "
                      f"(removed={removed:+d}), ddG={ddg:.2f}, score={score:.3f}")

                if best is None or removed > best["removed"] or (
                        removed == best["removed"] and ddg < best["ddg"]):
                    best = {
                        "pos": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "ddg": ddg,
                        "removed": removed, "new_n_ii": new_n_ii,
                        "new_ii_positions": new_pos,
                    }

        if best is None or best["removed"] <= 0:
            print(f"  No MHC-II improvement found in round {rnd+1}. Stopping Phase 2.")
            break

        # Apply best mutation
        seq_list[best["pos"] - 1] = best["mut_aa"]
        seq = "".join(seq_list)
        cumddg += best["ddg"]
        already_mutated.add(best["pos"])
        current_n_ii = best["new_n_ii"]

        # Refresh binder records for next round
        _, binders_current = run_mhcpan_ii_full(seq)

        mutations_applied.append({
            "position":       best["pos"],
            "wt_aa":          best["wt_aa"],
            "mut_aa":         best["mut_aa"],
            "ddg_pred":       best["ddg"],
            "binders_removed": best["removed"],
            "new_n_ii":       best["new_n_ii"],
            "phase":          "MHC-II",
        })
        print(f"  APPLIED: Pos {best['pos']} {best['wt_aa']}->{best['mut_aa']} "
              f"| MHC-II: {n_ii_wt} -> {best['new_n_ii']} (-{best['removed']}) "
              f"| cumDDG={cumddg:.2f}")

        if len(mutations_applied) >= MAX_MUTATIONS:
            break

    return seq, cumddg, mutations_applied


# ---------------------------------------------------------------------------
# Phase 1: MHC-I anchor scan
# ---------------------------------------------------------------------------
def phase1_mhci(
    seq: str,
    n_rounds: int,
    ddg_budget: float,
    cumddg_start: float,
    mutations_applied: list[dict],
    n_candidates: int = 5,
) -> tuple[str, float, list[dict]]:
    """Greedy MHC-I deimmunization targeting anchor positions p+1 and p+8.

    Returns (updated_seq, cumulative_ddg, mutations_applied).
    """
    cumddg = cumddg_start
    already_mutated = {m["position"] for m in mutations_applied}

    print(f"\n=== Phase 1: MHC-I anchor-position mutagenesis ===")
    print(f"  Budget remaining: {ddg_budget - cumddg:.2f} DDG, {n_rounds} rounds max")

    current_n_i, current_i_pos = run_mhcpan_i(seq)
    print(f"  Current MHC-I: {current_n_i} binders (WT: {WT_N_I})")

    for rnd in range(n_rounds):
        if not current_i_pos:
            print(f"  All MHC-I binders eliminated after {rnd} rounds.")
            break
        if cumddg >= ddg_budget:
            print(f"  DDG budget exhausted. Stopping Phase 1.")
            break

        anchor_pos = get_mhc_i_anchor_positions(current_i_pos)
        anchor_pos = [p for p in anchor_pos if p not in already_mutated]

        if not anchor_pos:
            print(f"  No more mutable MHC-I anchor positions. Stopping Phase 1.")
            break

        remaining_budget_per_mut = ddg_budget - cumddg
        seq_list = list(seq)
        best = None

        print(f"\n  Round {rnd+1}: {current_n_i} binders, "
              f"{len(anchor_pos)} anchor positions, "
              f"budget_remaining={remaining_budget_per_mut:.2f}")

        for pos in anchor_pos:
            wt_aa = seq_list[pos - 1]
            candidates = substitution_candidates(
                wt_aa, n=n_candidates, budget=min(MAX_PER_MUT_DDG, remaining_budget_per_mut)
            )
            for mut_aa, ddg in candidates:
                proposed = seq_list.copy()
                proposed[pos - 1] = mut_aa
                proposed_seq = "".join(proposed)

                new_n_i, new_i_pos = run_mhcpan_i(proposed_seq)
                removed = current_n_i - new_n_i

                print(f"    Pos {pos} {wt_aa}->{mut_aa}: MHC-I={new_n_i} "
                      f"(removed={removed:+d}), ddG={ddg:.2f}")

                if best is None or removed > best["removed"] or (
                        removed == best["removed"] and ddg < best["ddg"]):
                    best = {
                        "pos": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "ddg": ddg,
                        "removed": removed, "new_n_i": new_n_i,
                        "new_i_positions": new_i_pos,
                    }

        if best is None or best["removed"] <= 0:
            print(f"  No MHC-I improvement found in round {rnd+1}. Stopping Phase 1.")
            break

        seq_list[best["pos"] - 1] = best["mut_aa"]
        seq = "".join(seq_list)
        cumddg += best["ddg"]
        already_mutated.add(best["pos"])
        current_n_i = best["new_n_i"]
        current_i_pos = best["new_i_positions"]

        mutations_applied.append({
            "position":       best["pos"],
            "wt_aa":          best["wt_aa"],
            "mut_aa":         best["mut_aa"],
            "ddg_pred":       best["ddg"],
            "binders_removed": best["removed"],
            "new_n_i":        best["new_n_i"],
            "phase":          "MHC-I",
        })
        print(f"  APPLIED: Pos {best['pos']} {best['wt_aa']}->{best['mut_aa']} "
              f"| MHC-I: {WT_N_I} -> {best['new_n_i']} (-{best['removed']}) "
              f"| cumDDG={cumddg:.2f}")

        if len(mutations_applied) >= MAX_MUTATIONS:
            break

    return seq, cumddg, mutations_applied


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 70)
    print("Strategy C v2: MHC-II-first Core-Anchor Deimmunization of IS621")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"WT: n_I={WT_N_I}, n_II={WT_N_II}, S_Immuno={IS621_S_IMMUNO_WT}")
    print(f"P3 target: delta >= {P3_TARGET_DELTA}")
    print(f"Budget: Phase2_MHC-II={PHASE2_DDG_BUDGET}, Total={MAX_DDG_TOTAL}, "
          f"Per-mut={MAX_PER_MUT_DDG}, Max_muts={MAX_MUTATIONS}")
    print("=" * 70)

    if not Path(BINARY_I).exists():
        print(f"ERROR: netMHCpan not found at {BINARY_I}")
        sys.exit(1)
    if not Path(PERL_II).exists():
        print(f"ERROR: NetMHCIIpan not found at {PERL_II}")
        sys.exit(1)
    print("netMHCpan-4.1:  FOUND")
    print("NetMHCIIpan-4.0: FOUND")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mutations_applied: list[dict] = []
    seq = IS621_SEQUENCE

    # --- Phase 2: MHC-II first (priority budget = 10.0, up to 9 rounds) ---
    phase2_rounds = min(9, MAX_MUTATIONS)
    seq, cumddg, mutations_applied = phase2_mhcii(
        seq=seq,
        n_rounds=phase2_rounds,
        ddg_budget=PHASE2_DDG_BUDGET,
        mutations_applied=mutations_applied,
        n_candidates=5,
    )

    # --- Phase 1: MHC-I with remaining budget ---
    phase1_rounds = MAX_MUTATIONS - len(mutations_applied)
    if phase1_rounds > 0 and cumddg < MAX_DDG_TOTAL:
        seq, cumddg, mutations_applied = phase1_mhci(
            seq=seq,
            n_rounds=phase1_rounds,
            ddg_budget=MAX_DDG_TOTAL,
            cumddg_start=cumddg,
            mutations_applied=mutations_applied,
            n_candidates=5,
        )

    # --- Final evaluation ---
    print(f"\n=== Final evaluation (real netMHCpan) ===")
    print(f"  Variant: {len(mutations_applied)} mutations, cumDDG={cumddg:.2f}")
    final_n_i, _ = run_mhcpan_i(seq)
    final_n_ii, _ = run_mhcpan_ii_count(seq)
    final_combined = final_n_i + 0.5 * final_n_ii
    final_s_immuno = round(min(1.0, max(0.0, 1.0 - final_combined / S_IMMUNO_MAX_TOTAL)), 4)
    final_delta = round(final_s_immuno - IS621_S_IMMUNO_WT, 4)

    print(f"  MHC-I:  {WT_N_I} -> {final_n_i} (-{WT_N_I - final_n_i})")
    print(f"  MHC-II: {WT_N_II} -> {final_n_ii} (-{WT_N_II - final_n_ii})")
    print(f"  Combined: {61.0:.1f} -> {final_combined:.1f}")
    print(f"  S_Immuno: {IS621_S_IMMUNO_WT} -> {final_s_immuno:.4f} (delta={final_delta:+.4f})")
    print(f"  P3 pass: {'YES ✓' if final_delta >= P3_TARGET_DELTA else 'NO'}")

    # Build mutation summary with sequence context
    mut_summary = []
    for m in mutations_applied:
        pos = m["position"]
        ctx_start = max(0, pos - 4)
        ctx_end   = min(len(seq), pos + 4)
        ctx = seq[ctx_start:ctx_end]
        mut_summary.append({**m, "sequence_context": ctx})

    # Variant sequence with mutation markers
    print(f"\n  Mutations applied:")
    for m in mutations_applied:
        print(f"    [{m['phase']}] Pos {m['position']}: "
              f"{m['wt_aa']}->{m['mut_aa']} (ddG={m['ddg_pred']:.2f})")

    elapsed = time.time() - t0
    print(f"\nElapsed: {elapsed/60:.1f} min")

    # --- Save outputs ---
    manifest = {
        "method": "mhcii_first_core_anchor_greedy_v2",
        "n_variants": 1,
        "n_p3_pass": 1 if final_delta >= P3_TARGET_DELTA else 0,
        "wt_n_i": WT_N_I,
        "wt_n_ii": WT_N_II,
        "final_n_i": final_n_i,
        "final_n_ii": final_n_ii,
        "final_combined": final_combined,
        "final_s_immuno": final_s_immuno,
        "final_delta": final_delta,
        "cumulative_ddg": round(cumddg, 2),
        "total_mutations": len(mutations_applied),
        "passes_p3": final_delta >= P3_TARGET_DELTA,
        "elapsed_min": round(elapsed / 60, 1),
        "netmhcpan_used": True,
        "netmhciipan_used": True,
        "mhc_i_rank_threshold": 0.5,
        "mhc_ii_rank_threshold": 10.0,
        "s_immuno_max_total": S_IMMUNO_MAX_TOTAL,
        "budget_phase2_ddg": PHASE2_DDG_BUDGET,
        "budget_total_ddg": MAX_DDG_TOTAL,
        "budget_per_mut_ddg": MAX_PER_MUT_DDG,
        "max_mutations": MAX_MUTATIONS,
        "mutations": mut_summary,
    }

    manifest_path = OUTPUT_DIR / "deimmunized_v2_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nManifest: {manifest_path}")

    # FASTA
    fasta_path = OUTPUT_DIR / "deimmunized_v2.fasta"
    mut_str = "_".join(f"{m['wt_aa']}{m['position']}{m['mut_aa']}" for m in mutations_applied)
    fasta_path.write_text(
        f">IS621_deimmunized_v2_{mut_str}\n"
        f"{seq}\n"
    )
    print(f"FASTA:    {fasta_path}")

    # Parquet
    try:
        import pandas as pd
        df = pd.DataFrame([{
            "variant_id": f"IS621_deimmunized_v2",
            "scaffold_id": "IS621",
            "variant_sequence": seq,
            "wt_sequence": IS621_SEQUENCE,
            "total_mutations": len(mutations_applied),
            "cumulative_ddg": round(cumddg, 2),
            "final_n_i": final_n_i,
            "final_n_ii": final_n_ii,
            "final_combined": final_combined,
            "predicted_s_immuno": final_s_immuno,
            "predicted_s_immuno_delta": final_delta,
            "passes_p3": final_delta >= P3_TARGET_DELTA,
            "method": "mhcii_first_core_anchor_v2",
            "mutations_json": json.dumps(mut_summary),
        }])
        parquet_path = OUTPUT_DIR / "deimmunized_v2.parquet"
        df.to_parquet(str(parquet_path), compression="zstd", index=False)
        print(f"Parquet:  {parquet_path}")
    except ImportError:
        print("pandas not available — skipping parquet")

    # DONE marker
    done_info = {
        "status": "complete",
        "elapsed_min": round(elapsed / 60, 1),
        "final_delta": final_delta,
        "passes_p3": final_delta >= P3_TARGET_DELTA,
        "total_mutations": len(mutations_applied),
        "final_n_i": final_n_i,
        "final_n_ii": final_n_ii,
    }
    (OUTPUT_DIR / "DONE.txt").write_text(json.dumps(done_info, indent=2))

    print(f"\n{'=' * 70}")
    print(f"DONE. Elapsed: {elapsed/60:.1f} min")
    print(f"P3 pass: {'YES ✓' if final_delta >= P3_TARGET_DELTA else 'NO'} "
          f"(delta={final_delta:+.4f})")
    print(f"Results: {OUTPUT_DIR}")
    print("[DONE.txt written]")


if __name__ == "__main__":
    main()
