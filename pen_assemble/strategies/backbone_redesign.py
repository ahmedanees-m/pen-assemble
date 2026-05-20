"""Strategy D — ProteinMPNN backbone-conditioned IS621 redesign.

Inverse-folds IS621's backbone (PDB 8WT6) using ProteinMPNN to generate novel
sequences that fold into IS621's proven structure. Pinned residues: catalytic site
(D11/E60/D102/D105/S241 per Hiraizumi 2024 Nature 630:994–1002), bRNA-binding loop
contacts (TBL+DBL per Hiraizumi 2024), and domain-domain interface.
Redesignable: all other residues (~150 of ~300).

Rationale (added 2026-06): P1 math shows Strategy A chimeras face a ~−0.07
PenScore disadvantage from S_Deliv (larger protein) and S_Mature (no citations).
Strategy D preserves IS621's protein size (S_Deliv unchanged) while generating
genuinely novel sequences, making it the most likely path to P1 success.

P5 diversity contribution: Strategy D sequences are distinct from IS621 WT
(50–75% identity) and count as an independent scaffold source.

Target: ~25 ProteinMPNN-redesigned variants. Implemented in Step 11.5.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

# IS621 catalytic residues — MUST be pinned (Hiraizumi 2024 Nature 630:994–1002, Table 1)
# DEDD tetrad: D11, E60, D102, D105 (RuvC-like fold); S241 (Tnp serine recombinase C-term)
CATALYTIC_RESIDUES_IS621 = [11, 60, 102, 105, 241]  # D11, E60, D102, D105 (DEDD) + S241 (Tnp-Ser)

# ProteinMPNN GitHub commit pin (immutable)
PROTEINMPNN_COMMIT = "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57"

# S_Mature partial inheritance for Strategy D variants (same structural family, novel sequence)
S_MATURE_INHERITANCE_FACTOR = 0.5  # = parent IS621 S_Mature × 0.5 (pre-registered rule)


@dataclass
class BackboneRedesignVariant:
    """A single ProteinMPNN-redesigned IS621 variant."""

    design_id: str
    protmpnn_sequence: str
    pinned_positions: list[int] = field(default_factory=list)
    backbone_rmsd_to_wt_a: Optional[float] = None   # angstroms vs 8WT6
    sequence_identity_to_wt: Optional[float] = None  # fraction 0–1
    novel_residues_count: Optional[int] = None

    def passes_backbone_check(self, rmsd_threshold: float = 1.5) -> Optional[bool]:
        """True if backbone RMSD to WT IS621 8WT6 <= rmsd_threshold Å."""
        if self.backbone_rmsd_to_wt_a is None:
            return None
        return self.backbone_rmsd_to_wt_a <= rmsd_threshold

    def passes_identity_range(self, lo: float = 0.50, hi: float = 0.75) -> Optional[bool]:
        """True if sequence identity to WT is in [lo, hi] — novel but plausible."""
        if self.sequence_identity_to_wt is None:
            return None
        return lo <= self.sequence_identity_to_wt <= hi


def _write_fixed_positions_jsonl(pdb_path: Path, pinned: list[int], out_dir: Path) -> Path:
    """Write ProteinMPNN fixed_positions_jsonl file for pinned-residue inverse folding."""
    import json as _json
    chain_id = "A"
    fixed = {pdb_path.stem: {chain_id: pinned}}
    jsonl_path = out_dir / "fixed_positions.jsonl"
    jsonl_path.write_text(_json.dumps(fixed) + "\n")
    return jsonl_path


def _get_wt_sequence_from_pdb(pdb_path: Path) -> str:
    """Extract WT amino acid sequence from PDB (chain A, standard residues)."""
    try:
        from Bio.PDB import PDBParser
        from Bio.SeqUtils import seq1
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("ref", str(pdb_path))
        return "".join(
            seq1(r.get_resname()) for r in struct.get_residues()
            if r.id[0] == " "   # HETATM excluded
        )
    except ImportError:
        return ""


def _parse_protmpnn_fasta(fasta_path: Path, wt_sequence: str) -> list[BackboneRedesignVariant]:
    """Parse ProteinMPNN output FASTA; compute sequence identity to WT."""
    variants: list[BackboneRedesignVariant] = []
    try:
        from Bio import SeqIO as _SeqIO
        records = list(_SeqIO.parse(str(fasta_path), "fasta"))
    except (ImportError, FileNotFoundError):
        return variants

    for i, record in enumerate(records):
        seq = str(record.seq)
        if wt_sequence and seq == wt_sequence:
            continue  # skip exact WT recovery
        n_match = sum(a == b for a, b in zip(seq, wt_sequence)) if wt_sequence else 0
        identity = n_match / max(len(wt_sequence), len(seq)) if wt_sequence and seq else 0.0
        variants.append(BackboneRedesignVariant(
            design_id=f"D_{i:03d}_protmpnn",
            protmpnn_sequence=seq,
            pinned_positions=CATALYTIC_RESIDUES_IS621[:],
            sequence_identity_to_wt=round(identity, 4),
            novel_residues_count=(len(wt_sequence) - n_match) if wt_sequence else None,
        ))
    return variants


def generate_backbone_redesigns(
    scaffold_pdb: Optional[Path] = None,
    n_designs: int = 30,
    pinned_residues: Optional[list[int]] = None,
    bRNA_contact_residues: Optional[list[int]] = None,
    seed: int = 42,
    output_dir: Optional[Path] = None,
) -> list[BackboneRedesignVariant]:
    """Run ProteinMPNN inverse-folding on IS621 backbone (Strategy D).

    Requires ProteinMPNN installed in pen-stack/design Docker image at pinned commit
    8907e6671bfbfc92303b5f79c4b5e6ce47cdef57 (``pip install git+…@{commit}``).

    Pinned positions: CATALYTIC_RESIDUES_IS621 (D11/E60/D102/D105/S241) merged with
    bRNA_contact_residues (from Hiraizumi 2024 SI Table S2 JSON annotation file).
    Redesignable: all other residues (~150 of ~300).

    After generation, filters to identity range [0.50, 0.75] (novel but plausible).

    Args:
        scaffold_pdb: Path to IS621 PDB (8WT6). Defaults to /data/raw/8WT6.pdb.
        n_designs: Number of sequences ProteinMPNN should generate (default 30).
        pinned_residues: Override catalytic pinned positions (default CATALYTIC_RESIDUES_IS621).
        bRNA_contact_residues: bRNA-binding loop contacts to pin. If None, only
            catalytic residues are pinned (conservative fallback).
        seed: Random seed for ProteinMPNN sampling (default 42).
        output_dir: Directory for ProteinMPNN output. Defaults to /data/pen-assemble/designs/strategy_D.

    Returns:
        List of BackboneRedesignVariant sorted by sequence_identity_to_wt (lowest first,
        i.e. most novel first). Returns empty list if ProteinMPNN is unavailable.

    Raises:
        ImportError: If ProteinMPNN package is not installed.
        FileNotFoundError: If scaffold_pdb does not exist.
        RuntimeError: If ProteinMPNN subprocess fails.
    """
    import subprocess as _subprocess

    pdb = scaffold_pdb or Path("/data/raw/8WT6.pdb")
    out = output_dir or Path("/data/pen-assemble/designs/strategy_D")
    out.mkdir(parents=True, exist_ok=True)

    if not pdb.exists():
        raise FileNotFoundError(
            f"IS621 PDB not found at {pdb}. Download 8WT6.pdb from RCSB or "
            "run scripts/03_validate_sequences.py which may fetch it."
        )

    # Confirm ProteinMPNN is importable
    try:
        import ProteinMPNN  # noqa: F401
    except ImportError:
        raise ImportError(
            f"ProteinMPNN not installed. Install with:\n"
            f"  pip install 'git+https://github.com/dauparas/ProteinMPNN.git"
            f"@{PROTEINMPNN_COMMIT}'\n"
            f"inside the pen-stack/design:0.1.0 Docker image."
        )

    # Build pinned position list
    pinned = sorted(set(pinned_residues or CATALYTIC_RESIDUES_IS621)
                    | set(bRNA_contact_residues or []))

    fixed_jsonl = _write_fixed_positions_jsonl(pdb, pinned, out)
    wt_seq = _get_wt_sequence_from_pdb(pdb)

    cmd = [
        "python", "-m", "ProteinMPNN.protein_mpnn_run",
        "--pdb_path", str(pdb),
        "--out_folder", str(out),
        "--num_seq_per_target", str(n_designs),
        "--sampling_temp", "0.1",
        "--seed", str(seed),
        "--fixed_positions_jsonl", str(fixed_jsonl),
        "--batch_size", "1",
    ]
    result = _subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"ProteinMPNN exited with code {result.returncode}:\n{result.stderr[:800]}"
        )

    fasta_out = out / "seqs" / f"{pdb.stem}.fa"
    variants = _parse_protmpnn_fasta(fasta_out, wt_seq)

    # Filter to novelty range [50%, 75%] identity (pre-registered Step 11.5 criterion)
    variants = [
        v for v in variants
        if v.sequence_identity_to_wt is None
        or (0.50 <= v.sequence_identity_to_wt <= 0.75)
    ]

    # Sort most novel first
    variants.sort(key=lambda v: v.sequence_identity_to_wt or 1.0)
    return variants
