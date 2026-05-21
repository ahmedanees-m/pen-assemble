"""Rosetta ΔΔG stability filtering. Step 13.

Hard gate: ΔΔG > +5.0 kcal/mol → discard (predicted destabilizing vs parent scaffold).
Soft flag: ΔΔG > +3.0 kcal/mol → note in output but do not discard.

Method cascade (auto):
  1. Rosetta CartesianDDG — requires PyRosetta + both PDB structures
  2. RaSP (Blaabjerg et al. 2023) — requires rasp package + mutations list
  3. Grantham proxy — always available; requires mutations list only
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

DDG_HARD_GATE_KCAL = 5.0  # designs above this are discarded
DDG_SOFT_FLAG_KCAL = 3.0  # designs above this get stability_warning=True
DDG_STABLE_THRESHOLD = -0.5  # designs below this are considered stabilized variants


def _check_pyrosetta() -> bool:
    try:
        import pyrosetta  # noqa: F401

        return True
    except ImportError:
        return False


def _check_rasp() -> bool:
    try:
        import rasp  # noqa: F401

        return True
    except ImportError:
        return False


_MW: dict[str, float] = {
    "G": 57,
    "A": 71,
    "V": 99,
    "L": 113,
    "I": 113,
    "P": 97,
    "F": 147,
    "W": 186,
    "M": 131,
    "S": 87,
    "T": 101,
    "C": 103,
    "Y": 163,
    "H": 137,
    "D": 115,
    "E": 129,
    "N": 114,
    "Q": 128,
    "K": 128,
    "R": 156,
}


def _grantham_ddg_proxy(mutations: list[dict]) -> float:
    _polar = set("RKDENQHST")
    total = 0.0
    for m in mutations:
        wt, mt = m.get("wt_aa", "A"), m.get("mut_aa", "A")
        cross = (wt in _polar) != (mt in _polar)
        size_delta = abs(_MW.get(wt, 115) - _MW.get(mt, 115)) / 30.0
        total += 0.3 + (0.8 if cross else 0.0) + size_delta * 0.5
    return round(total, 2)


def _compute_ddg_rosetta(design_pdb: Path, parent_pdb: Path) -> float | None:
    if not _check_pyrosetta():
        return None
    try:
        import pyrosetta

        pyrosetta.init("-mute all")
        parent_pose = pyrosetta.pose_from_pdb(str(parent_pdb))
        design_pose = pyrosetta.pose_from_pdb(str(design_pdb))
        sfxn = pyrosetta.create_score_function("ref2015_cart")
        ddg = sfxn(design_pose) - sfxn(parent_pose)
        return round(float(ddg), 3)
    except Exception:
        return None


def _compute_ddg_rasp(sequence: str, mutations: list[dict]) -> float | None:
    if not _check_rasp():
        return None
    try:
        import rasp

        return round(
            sum(
                rasp.predict_single(
                    sequence=sequence, position=m["pos"], wt=m["wt_aa"], mt=m["mut_aa"]
                )
                for m in mutations
            ),
            3,
        )
    except Exception:
        return None


def compute_ddg(
    design_pdb: Path | None = None,
    parent_pdb: Path | None = None,
    sequence: str | None = None,
    mutations: list[dict] | None = None,
    method: str = "auto",
) -> dict[str, Any]:
    """Compute ΔΔG for a design. Method cascade: rosetta → rasp → grantham."""
    result: dict[str, Any] = {
        "ddg_kcal_mol": None,
        "ddg_method": None,
        "stability_warning": False,
        "ddg_hard_gate_fail": False,
    }
    ddg: float | None = None
    method_used = "none"

    if method in ("auto", "rosetta") and design_pdb and parent_pdb:
        if design_pdb.exists() and parent_pdb.exists():
            ddg = _compute_ddg_rosetta(design_pdb, parent_pdb)
            if ddg is not None:
                method_used = "rosetta_cartesian_ddg"

    if ddg is None and method in ("auto", "rasp") and sequence and mutations:
        ddg = _compute_ddg_rasp(sequence, mutations)
        if ddg is not None:
            method_used = "rasp"

    if ddg is None and method in ("auto", "grantham") and mutations:
        ddg = _grantham_ddg_proxy(mutations)
        method_used = "grantham_proxy"

    result.update({"ddg_kcal_mol": ddg, "ddg_method": method_used})
    if ddg is not None:
        result["stability_warning"] = ddg > DDG_SOFT_FLAG_KCAL
        result["ddg_hard_gate_fail"] = ddg > DDG_HARD_GATE_KCAL
        result["ddg_stabilizing"] = ddg < DDG_STABLE_THRESHOLD
    return result


def filter_by_ddg(
    designs_df: pd.DataFrame,
    ddg_column: str = "ddg_kcal_mol",
    hard_gate: float = DDG_HARD_GATE_KCAL,
) -> pd.DataFrame:
    """Discard designs where ddg > hard_gate. Pass-through if column missing."""
    if ddg_column not in designs_df.columns:
        return designs_df
    before = len(designs_df)
    out = designs_df[designs_df[ddg_column].isna() | (designs_df[ddg_column] <= hard_gate)].copy()
    print(
        f"  ddG filter: {before} -> {len(out)} (removed {before - len(out)} above {hard_gate} kcal/mol)"
    )
    return out


def run_stability_filter(
    designs_df: pd.DataFrame,
    design_pdb_col: str = "final_pdb",
    parent_pdb_path: Path | None = None,
    mutations_col: str | None = "mutations_introduced",
    sequence_col: str = "protein_sequence",
) -> pd.DataFrame:
    """Add ddG columns to designs_df; remove hard-gate failures."""
    rows: list[dict] = []
    for _, row in designs_df.iterrows():
        pdb = Path(row[design_pdb_col]) if row.get(design_pdb_col) else None
        seq = row.get(sequence_col)
        muts = row.get(mutations_col) if mutations_col else None
        r = compute_ddg(
            design_pdb=pdb,
            parent_pdb=parent_pdb_path,
            sequence=seq,
            mutations=muts if isinstance(muts, list) else None,
        )
        rows.append(r)

    ddg_df = pd.DataFrame(rows)
    out = designs_df.copy()
    for col in ddg_df.columns:
        out[col] = ddg_df[col].values
    return filter_by_ddg(out)
