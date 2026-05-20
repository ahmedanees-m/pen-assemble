"""Strategy D VM runner — ProteinMPNN backbone-conditioned IS621 redesign.

Self-contained: deploys to VM and runs ProteinMPNN CLI on 8WT6.pdb (chain A).
Pins 118 catalytic + bRNA contact positions (computed from 8WT6 structure).
Generates 30 designs, filters by sequence identity 50–75% vs IS621 WT.

Deploy:
    # 1. Upload PDB and this script:
    #    sftp -> put 8WT6.pdb ~/8WT6.pdb
    #    sftp -> put 13b_strategy_D_vm_runner.py ~/strategy_D_runner.py
    # 2. Run:
    #    nohup python3 -u ~/strategy_D_runner.py > ~/strategy_D.log 2>&1 &

Output: ~/strategy_D_results/
    strategy_D_designs.fasta
    strategy_D_designs.parquet
    strategy_D_manifest.json
    DONE.txt

Pre-requisites on VM:
    pip install git+https://github.com/dauparas/ProteinMPNN.git@8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
    (torch already available)

IS621 WT sequence (342 aa, UniProt A0A2X3M8B0):
    Catalytic residues (Hiraizumi 2024 Nature 630:994–1002):
      D11, E60, D102, D105, S241  — MUST pin
    bRNA contacts (computed from 8WT6 chain A vs RNA chains E-J, <5A):
      117 positions  — MUST pin
    Total pinned: 118 (5 catalytic + 117 bRNA, with overlap)
    Redesignable: 224 positions (65% of IS621)
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
# Constants
# ---------------------------------------------------------------------------
HOME = Path.home()

IS621_SEQUENCE = (
    "MDRFFPVIRICKVGFTMEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAH"
    "ICIEATGTYMEPVAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERA"
    "LRALVVRHQALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGI"
    "GEKTSAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
    "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA"
)
assert len(IS621_SEQUENCE) == 342

# Catalytic residues (1-indexed, Hiraizumi 2024)
CATALYTIC = [11, 60, 102, 105, 241]

# bRNA contact residues (1-indexed, from 8WT6 chain A vs RNA E-J, <5A)
BRNA_CONTACTS = [
    11, 12, 13, 14, 15, 26, 27, 28, 29, 30, 46, 60, 61, 62, 63, 64, 65, 66, 67,
    83, 84, 85, 86, 89, 98, 99, 100, 101, 102, 103, 104, 105, 107, 132, 136, 138,
    142, 143, 145, 146, 147, 149, 150, 151, 153, 154, 156, 199, 201, 202, 203, 204,
    205, 206, 207, 208, 209, 221, 222, 223, 224, 225, 226, 227, 228, 230, 231, 234,
    235, 236, 237, 250, 251, 252, 253, 254, 255, 256, 257, 260, 261, 264, 265, 266,
    268, 269, 272, 273, 274, 280, 283, 284, 287, 289, 290, 291, 292, 293, 295, 296,
    297, 299, 300, 301, 304, 305, 307, 308, 309, 311, 312, 314, 315, 316, 317,
    320, 321,
]

# Total pinned (merged)
PINNED_POSITIONS = sorted(set(CATALYTIC) | set(BRNA_CONTACTS))

OUTPUT_DIR   = HOME / "strategy_D_results"
PDB_PATH     = HOME / "8WT6.pdb"
N_DESIGNS    = 30
SEED         = 42
SAMPLING_TEMP = 0.1   # lower = more conservative (closer to WT at free positions)
CHAIN_ID     = "A"    # IS621 monomer

# Identity filter: keep designs with 50–75% similarity to IS621 WT
IDENTITY_LO = 0.50
IDENTITY_HI = 0.75


# ---------------------------------------------------------------------------
# Find ProteinMPNN script
# ---------------------------------------------------------------------------
def find_proteinmpnn_script() -> Path:
    """Locate protein_mpnn_run.py from the pip-installed package or git clone."""
    import importlib.util

    # Try pip-installed package first
    candidates = [
        HOME / ".local/lib/python3.10/site-packages/ProteinMPNN/protein_mpnn_run.py",
        HOME / ".local/lib/python3.10/site-packages/protein_mpnn_run.py",
        HOME / "ProteinMPNN/protein_mpnn_run.py",   # git clone fallback
        HOME / "proteinmpnn/protein_mpnn_run.py",
    ]
    # Also check site-packages subdirs
    import site
    for sp in site.getsitepackages() + [site.getusersitepackages()]:
        candidates.append(Path(sp) / "ProteinMPNN" / "protein_mpnn_run.py")
        candidates.append(Path(sp) / "protein_mpnn_run.py")

    for c in candidates:
        if c.exists():
            return c

    raise FileNotFoundError(
        "protein_mpnn_run.py not found. Install with:\n"
        "  pip install --user git+https://github.com/dauparas/ProteinMPNN.git"
        "@8907e6671bfbfc92303b5f79c4b5e6ce47cdef57\n"
        "OR clone: git clone https://github.com/dauparas/ProteinMPNN.git ~/ProteinMPNN"
    )


# ---------------------------------------------------------------------------
# Prepare ProteinMPNN input files
# ---------------------------------------------------------------------------
def write_chain_jsonl(out_dir: Path) -> Path:
    """Write chain_id.jsonl specifying which chain to redesign."""
    data = {"8WT6": [CHAIN_ID]}
    p = out_dir / "chain_ids.jsonl"
    p.write_text(json.dumps(data) + "\n")
    return p


def write_fixed_positions_jsonl(out_dir: Path) -> Path:
    """Write fixed_positions.jsonl to pin catalytic + bRNA contact positions."""
    data = {"8WT6": {CHAIN_ID: PINNED_POSITIONS}}
    p = out_dir / "fixed_positions.jsonl"
    p.write_text(json.dumps(data) + "\n")
    return p


# ---------------------------------------------------------------------------
# Sequence identity
# ---------------------------------------------------------------------------
def sequence_identity(seq1: str, seq2: str) -> float:
    """Fraction of identical positions (aligned at common length)."""
    n = min(len(seq1), len(seq2))
    if n == 0:
        return 0.0
    return sum(a == b for a, b in zip(seq1[:n], seq2[:n])) / n


# ---------------------------------------------------------------------------
# Parse ProteinMPNN output FASTA
# ---------------------------------------------------------------------------
def parse_proteinmpnn_fasta(fasta_path: Path) -> list[dict]:
    """Parse ProteinMPNN output FASTA.

    ProteinMPNN FASTA headers look like:
      >8WT6, score=0.5432, seq_recovery=0.9123, T=0.1, sample=1
    Returns list of {header, sequence, score, recovery, sample_id}.
    """
    text = fasta_path.read_text()
    entries = []
    cur_h, cur_s = None, ""

    for line in text.splitlines():
        if line.startswith(">"):
            if cur_h is not None and cur_s:
                entries.append(_parse_entry(cur_h, cur_s))
            cur_h, cur_s = line, ""
        elif line.strip():
            cur_s += line.strip()
    if cur_h and cur_s:
        entries.append(_parse_entry(cur_h, cur_s))

    return entries


def _parse_entry(header: str, sequence: str) -> dict:
    score, recovery, sample = 0.0, 0.0, 0
    for part in header.split(","):
        part = part.strip()
        if "score=" in part:
            try:
                score = float(part.split("=")[1])
            except ValueError:
                pass
        if "seq_recovery=" in part:
            try:
                recovery = float(part.split("=")[1])
            except ValueError:
                pass
        if "sample=" in part:
            try:
                sample = int(part.split("=")[1])
            except ValueError:
                pass
    return {
        "header": header,
        "sequence": sequence,
        "protmpnn_score": score,
        "seq_recovery": recovery,
        "sample_id": sample,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 70)
    print("Strategy D: ProteinMPNN Backbone-Conditioned IS621 Redesign")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"PDB: {PDB_PATH}")
    print(f"Chain: {CHAIN_ID}, N_designs: {N_DESIGNS}, T: {SAMPLING_TEMP}, seed: {SEED}")
    print(f"Pinned positions: {len(PINNED_POSITIONS)} (catalytic + bRNA contacts)")
    print(f"Redesignable: {342 - len(PINNED_POSITIONS)} positions")
    print("=" * 70)

    if not PDB_PATH.exists():
        print(f"ERROR: PDB not found at {PDB_PATH}")
        sys.exit(1)

    # Find ProteinMPNN script
    try:
        mpnn_script = find_proteinmpnn_script()
        print(f"\nProteinMPNN script: {mpnn_script}")
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        # Try git clone as fallback
        print("Attempting git clone of ProteinMPNN...")
        rc = subprocess.run(
            ["git", "clone", "https://github.com/dauparas/ProteinMPNN.git",
             str(HOME / "ProteinMPNN")],
            capture_output=True, text=True, timeout=120,
        ).returncode
        if rc == 0:
            mpnn_script = HOME / "ProteinMPNN/protein_mpnn_run.py"
            print(f"Cloned. Script at: {mpnn_script}")
        else:
            print("Git clone also failed. Exiting.")
            sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tmp_input  = OUTPUT_DIR / "mpnn_input"
    tmp_input.mkdir(exist_ok=True)
    tmp_output = OUTPUT_DIR / "mpnn_output"
    tmp_output.mkdir(exist_ok=True)

    # Prepare PDB symlink/copy in tmp_input
    pdb_input = tmp_input / "8WT6.pdb"
    if not pdb_input.exists():
        import shutil
        shutil.copy(str(PDB_PATH), str(pdb_input))

    # Write JSON input files
    chain_jsonl  = write_chain_jsonl(tmp_input)
    fixed_jsonl  = write_fixed_positions_jsonl(tmp_input)

    print(f"\nRunning ProteinMPNN...")
    print(f"  {N_DESIGNS} designs × T={SAMPLING_TEMP} × chain {CHAIN_ID}")
    print(f"  Pinned: {len(PINNED_POSITIONS)} positions")
    print(f"  Output: {tmp_output}")

    cmd = [
        sys.executable, str(mpnn_script),
        "--pdb_path",              str(pdb_input),
        "--chain_id_jsonl",        str(chain_jsonl),
        "--fixed_positions_jsonl", str(fixed_jsonl),
        "--out_folder",            str(tmp_output),
        "--num_seq_per_target",    str(N_DESIGNS),
        "--sampling_temp",         str(SAMPLING_TEMP),
        "--seed",                  str(SEED),
        "--batch_size",            "1",   # CPU-friendly
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if result.returncode != 0:
        print(f"ERROR: ProteinMPNN failed (rc={result.returncode})")
        print(f"STDOUT: {result.stdout[:500]}")
        print(f"STDERR: {result.stderr[:500]}")
        sys.exit(1)

    print(f"ProteinMPNN completed in {(time.time()-t0)/60:.1f} min")
    print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)

    # Find output FASTA
    seqs_dir = tmp_output / "seqs"
    fasta_files = list(seqs_dir.glob("*.fa")) if seqs_dir.exists() else []
    if not fasta_files:
        fasta_files = list(tmp_output.glob("**/*.fa")) + list(tmp_output.glob("**/*.fasta"))

    if not fasta_files:
        print(f"ERROR: No output FASTA found in {tmp_output}")
        sys.exit(1)

    fasta_file = fasta_files[0]
    print(f"\nParsing output: {fasta_file}")
    all_entries = parse_proteinmpnn_fasta(fasta_file)
    print(f"  Total sequences generated: {len(all_entries)}")

    # Filter: skip WT (first entry), apply identity filter
    designs = []
    for i, e in enumerate(all_entries):
        if i == 0 and e["sample_id"] == 0:
            print(f"  Skipping WT entry (sample_id=0)")
            continue
        seq = e["sequence"]
        if len(seq) < 200:
            continue   # likely a fragment
        identity = sequence_identity(seq, IS621_SEQUENCE)
        e["sequence_identity_to_wt"] = round(identity, 4)
        e["passes_identity_filter"] = IDENTITY_LO <= identity <= IDENTITY_HI

        # Verify pinned positions are preserved
        pinned_ok = all(
            pos <= len(seq) and seq[pos - 1] == IS621_SEQUENCE[pos - 1]
            for pos in PINNED_POSITIONS
            if pos <= min(len(seq), len(IS621_SEQUENCE))
        )
        e["pinned_positions_preserved"] = pinned_ok

        designs.append(e)

    print(f"  After length filter: {len(designs)} sequences")
    passes = [d for d in designs if d["passes_identity_filter"]]
    print(f"  Pass identity filter ({IDENTITY_LO:.0%}–{IDENTITY_HI:.0%}): {len(passes)}")

    # Final selection: prefer designs in identity range, sorted by ProteinMPNN score
    final = sorted(passes, key=lambda x: x["protmpnn_score"]) if passes else designs
    final = final[:N_DESIGNS]  # cap at requested count

    print(f"\n  Final designs: {len(final)}")
    for d in final[:5]:
        print(f"    sample={d['sample_id']}, score={d['protmpnn_score']:.4f}, "
              f"identity={d['sequence_identity_to_wt']:.3f}, "
              f"pinned_ok={d['pinned_positions_preserved']}")

    # Write outputs
    out_fasta = OUTPUT_DIR / "strategy_D_designs.fasta"
    with out_fasta.open("w") as f:
        for i, d in enumerate(final, 1):
            f.write(f">D{i:03d}_IS621_ProtMPNN_T{SAMPLING_TEMP}_sample{d['sample_id']} "
                    f"score={d['protmpnn_score']:.4f} "
                    f"identity={d['sequence_identity_to_wt']:.3f} "
                    f"pinned_ok={d['pinned_positions_preserved']}\n")
            seq = d["sequence"]
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + "\n")
    print(f"\nFASTA: {out_fasta}")

    # Parquet
    try:
        import pandas as pd
        rows = []
        for i, d in enumerate(final, 1):
            rows.append({
                "design_id": f"D{i:03d}",
                "strategy": "D",
                "variant_sequence": d["sequence"],
                "wt_sequence": IS621_SEQUENCE,
                "sequence_length": len(d["sequence"]),
                "sequence_identity_to_wt": d["sequence_identity_to_wt"],
                "protmpnn_score": d["protmpnn_score"],
                "protmpnn_seq_recovery": d["seq_recovery"],
                "protmpnn_sample_id": d["sample_id"],
                "protmpnn_temp": SAMPLING_TEMP,
                "n_pinned": len(PINNED_POSITIONS),
                "n_redesignable": 342 - len(PINNED_POSITIONS),
                "pinned_positions_preserved": d["pinned_positions_preserved"],
                "passes_identity_filter": d["passes_identity_filter"],
            })
        df = pd.DataFrame(rows)
        parquet_path = OUTPUT_DIR / "strategy_D_designs.parquet"
        df.to_parquet(str(parquet_path), compression="zstd", index=False)
        print(f"Parquet: {parquet_path}")
    except ImportError:
        print("pandas not available — skipping parquet")

    elapsed = time.time() - t0
    manifest = {
        "strategy": "D",
        "method": "ProteinMPNN",
        "pdb": "8WT6",
        "chain": CHAIN_ID,
        "n_designs": len(final),
        "n_pinned": len(PINNED_POSITIONS),
        "n_redesignable": 342 - len(PINNED_POSITIONS),
        "sampling_temp": SAMPLING_TEMP,
        "seed": SEED,
        "identity_filter": [IDENTITY_LO, IDENTITY_HI],
        "n_pass_identity": len(passes),
        "elapsed_min": round(elapsed / 60, 1),
        "proteinmpnn_commit": "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57",
    }
    (OUTPUT_DIR / "strategy_D_manifest.json").write_text(json.dumps(manifest, indent=2))
    (OUTPUT_DIR / "DONE.txt").write_text(json.dumps({
        "status": "complete", "n_designs": len(final),
        "elapsed_min": round(elapsed / 60, 1),
    }))
    print(f"\nDONE. Elapsed: {elapsed/60:.1f} min. {len(final)} designs ready.")


if __name__ == "__main__":
    main()
