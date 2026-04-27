"""PDB / structure file helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def parse_plddt_from_pdb(pdb_path: Path) -> np.ndarray:
    """Extract per-residue pLDDT from B-factor column of a PDB file."""
    plddts = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    plddts.append(float(line[60:66].strip()))
                except ValueError:
                    pass
    return np.array(plddts, dtype=float)


def compute_backbone_rmsd(pdb_a: Path, pdb_b: Path) -> float | None:
    """Compute Cα RMSD between two PDB structures (requires BioPython). Step 12 helper."""
    try:
        from Bio.PDB import PDBParser, Superimposer

        parser = PDBParser(QUIET=True)
        struct_a = parser.get_structure("A", str(pdb_a))
        struct_b = parser.get_structure("B", str(pdb_b))
        atoms_a = [a for a in struct_a.get_atoms() if a.get_name() == "CA"]
        atoms_b = [a for a in struct_b.get_atoms() if a.get_name() == "CA"]
        n = min(len(atoms_a), len(atoms_b))
        if n == 0:
            return None
        sup = Superimposer()
        sup.set_atoms(atoms_a[:n], atoms_b[:n])
        return float(sup.rms)
    except ImportError as exc:
        raise ImportError("BioPython required: pip install biopython>=1.83") from exc
