"""Structure prediction wrappers: Boltz-1 (primary), ESMFold (pre-filter), ColabFold (fallback).

Two-tier funnel (Step 12):
  Tier 1 - ESMFold on all ~500 candidates (~30 sec each on V100); drops mean pLDDT < 50.
  Tier 2 - Boltz-1 AF3-class on ~200 survivors (~5-10 min each on A100).

Boltz-1 fallback policy:
  If Boltz-1 RNA complex prediction fails (boltz==0.4.1 RNA multimer support maturing),
  fallback is:
    1. Rerun monomer-only Boltz-1 for protein chain alone.
    2. Run RNAfold (ViennaRNA) for bRNA secondary structure separately.
    3. Use ColabFold v1.5.5 (AF2-multimer) for protein+bRNA complex.
  Logged to structure_fallback_log.parquet. Designs with fallback get
  boltz_rna_fallback=True - reported transparently in paper.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

# ---------------------------------------------------------------------------
# Tool availability checks
# ---------------------------------------------------------------------------


def _check_tool(cmd: str, args: list[str] | None = None) -> bool:
    try:
        r = subprocess.run([cmd] + (args or ["--help"]), capture_output=True, timeout=10)
        return r.returncode in (0, 1)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _check_boltz() -> bool:
    return _check_tool("boltz", ["--help"])


def _check_esmfold() -> bool:
    try:
        import esm  # noqa: F401

        return True
    except ImportError:
        return False


def _check_colabfold() -> bool:
    return _check_tool("colabfold_batch", ["--help"])


# ---------------------------------------------------------------------------
# PDB mean pLDDT extraction
# ---------------------------------------------------------------------------


def _extract_plddt_from_pdb(pdb_path: Path) -> dict[str, float]:
    """Parse B-factor column (pLDDT) from PDB. Returns {mean, min, max, per_residue:[...]}."""
    plddts: list[float] = []
    seen_residues: set[int] = set()
    try:
        for line in pdb_path.read_text().splitlines():
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":  # Cα only, one per residue
                continue
            res_seq = int(line[22:26].strip())
            if res_seq in seen_residues:
                continue
            seen_residues.add(res_seq)
            try:
                bfactor = float(line[60:66].strip())
                plddts.append(bfactor)
            except ValueError:
                pass
    except Exception:
        return {"mean": 0.0, "min": 0.0, "max": 0.0, "per_residue": []}
    if not plddts:
        return {"mean": 0.0, "min": 0.0, "max": 0.0, "per_residue": []}
    import statistics

    return {
        "mean": round(statistics.mean(plddts), 2),
        "min": round(min(plddts), 2),
        "max": round(max(plddts), 2),
        "per_residue": [round(v, 2) for v in plddts],
    }


# ---------------------------------------------------------------------------
# Tier 1: ESMFold pre-filter
# ---------------------------------------------------------------------------


def predict_structure_esmfold(
    sequence: str,
    design_id: str,
    output_dir: Path,
) -> dict[str, Any]:
    """Run ESMFold on a single sequence. Returns {pdb_path, mean_plddt, success}.

    Requires: pip install fair-esm (pen-stack/plm image).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = output_dir / f"{design_id}_esmfold.pdb"
    result: dict[str, Any] = {
        "design_id": design_id,
        "tool": "ESMFold",
        "pdb_path": None,
        "mean_plddt": None,
        "success": False,
    }

    if not _check_esmfold():
        result["error"] = "ESMFold (fair-esm) not installed"
        return result

    try:
        import esm
        import torch

        model = esm.pretrained.esmfold_v1()
        model = model.eval().cuda() if torch.cuda.is_available() else model.eval()

        with torch.no_grad():
            output = model.infer_pdb(sequence)

        pdb_path.write_text(output)
        plddt_info = _extract_plddt_from_pdb(pdb_path)
        result.update(
            {
                "pdb_path": str(pdb_path),
                "mean_plddt": plddt_info["mean"],
                "success": True,
            }
        )
    except Exception as e:
        result["error"] = str(e)

    return result


# ---------------------------------------------------------------------------
# Tier 2: Boltz-1
# ---------------------------------------------------------------------------


def _write_boltz_input(
    sequence: str,
    design_id: str,
    tmpdir: Path,
    brna_sequence: str | None = None,
) -> Path:
    """Write Boltz-1 input YAML (or FASTA for monomer mode)."""
    if brna_sequence:
        # Multimer mode: protein + RNA
        content = {
            "sequences": [
                {"protein": {"id": "A", "sequence": sequence}},
                {"rna": {"id": "B", "sequence": brna_sequence}},
            ]
        }
        yaml_path = tmpdir / f"{design_id}.yaml"
        import yaml as _yaml

        yaml_path.write_text(_yaml.dump(content, default_flow_style=False))
        return yaml_path
    else:
        fasta_path = tmpdir / f"{design_id}.fasta"
        fasta_path.write_text(f">{design_id}\n{sequence}\n")
        return fasta_path


def predict_structure_boltz(
    sequence: str,
    design_id: str,
    output_dir: Path,
    num_models: int = 1,
    num_recycle: int = 3,
    brna_sequence: str | None = None,
) -> dict[str, Any]:
    """Tier 2: Boltz-1 AF3-class prediction. Returns {pdb_path, mean_plddt, success, tool}.

    If brna_sequence provided (RNA multimer mode) and Boltz-1 RNA complex fails,
    falls back to ColabFold v1.5.5 (see module docstring).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "design_id": design_id,
        "tool": "Boltz-1",
        "pdb_path": None,
        "mean_plddt": None,
        "success": False,
        "boltz_rna_fallback": False,
    }

    if not _check_boltz():
        result["error"] = "boltz not installed (pip install boltz)"
        return result

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_p = Path(tmpdir)
        input_path = _write_boltz_input(sequence, design_id, tmpdir_p, brna_sequence)
        out_sub = tmpdir_p / "out"

        cmd = [
            "boltz",
            "predict",
            str(input_path),
            "--out_dir",
            str(out_sub),
            "--num_models",
            str(num_models),
            "--num_recycle",
            str(num_recycle),
            "--override",
        ]
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=1800)  # 30 min max
            boltz_ok = r.returncode == 0
        except subprocess.TimeoutExpired:
            boltz_ok = False
            result["error"] = "Boltz-1 timed out (30 min)"

        if boltz_ok:
            # Find output PDB
            pdb_candidates = list(out_sub.rglob("*.pdb"))
            if pdb_candidates:
                src_pdb = pdb_candidates[0]
                dst_pdb = output_dir / f"{design_id}_boltz.pdb"
                dst_pdb.write_bytes(src_pdb.read_bytes())
                plddt_info = _extract_plddt_from_pdb(dst_pdb)
                result.update(
                    {
                        "pdb_path": str(dst_pdb),
                        "mean_plddt": plddt_info["mean"],
                        "success": True,
                    }
                )
                return result

        # Boltz-1 failed - if RNA complex mode, try ColabFold fallback
        if brna_sequence:
            result["boltz_rna_fallback"] = True
            fallback = predict_structure_colabfold(
                sequence=sequence,
                design_id=design_id,
                output_dir=output_dir,
                brna_sequence=brna_sequence,
            )
            result.update(
                {
                    "tool": "Boltz-1+ColabFold_fallback",
                    "pdb_path": fallback.get("pdb_path"),
                    "mean_plddt": fallback.get("mean_plddt"),
                    "success": fallback.get("success", False),
                    "error": fallback.get("error"),
                }
            )

    return result


# ---------------------------------------------------------------------------
# Boltz-1 fallback: ColabFold
# ---------------------------------------------------------------------------


def predict_structure_colabfold(
    sequence: str,
    design_id: str,
    output_dir: Path,
    brna_sequence: str | None = None,
    use_templates: bool = False,
) -> dict[str, Any]:
    """ColabFold v1.5.5 AF2-multimer fallback. Returns {pdb_path, mean_plddt, success}."""
    output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "design_id": design_id,
        "tool": "ColabFold",
        "pdb_path": None,
        "mean_plddt": None,
        "success": False,
    }

    if not _check_colabfold():
        result["error"] = "colabfold_batch not on PATH"
        return result

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_p = Path(tmpdir)
        fasta_path = tmpdir_p / f"{design_id}.fasta"
        if brna_sequence:
            fasta_path.write_text(
                f">{design_id}_protein\n{sequence}\n>{design_id}_rna\n{brna_sequence}\n"
            )
        else:
            fasta_path.write_text(f">{design_id}\n{sequence}\n")
        out_sub = tmpdir_p / "out"

        cmd = [
            "colabfold_batch",
            str(fasta_path),
            str(out_sub),
            "--num-models",
            "1",
            "--num-recycle",
            "3",
        ]
        if not use_templates:
            cmd.append("--templates")
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=3600)
            cf_ok = r.returncode == 0
        except subprocess.TimeoutExpired:
            result["error"] = "ColabFold timed out (60 min)"
            return result

        if cf_ok:
            pdbs = sorted(out_sub.rglob("*relaxed*.pdb")) or sorted(out_sub.rglob("*.pdb"))
            if pdbs:
                dst_pdb = output_dir / f"{design_id}_colabfold.pdb"
                dst_pdb.write_bytes(pdbs[0].read_bytes())
                plddt_info = _extract_plddt_from_pdb(dst_pdb)
                result.update(
                    {
                        "pdb_path": str(dst_pdb),
                        "mean_plddt": plddt_info["mean"],
                        "success": True,
                    }
                )
        else:
            result["error"] = "ColabFold exited non-zero"

    return result


# ---------------------------------------------------------------------------
# Unified entry point
# ---------------------------------------------------------------------------


def predict_structure(
    sequence: str,
    design_id: str,
    output_dir: Path,
    brna_sequence: str | None = None,
    esmfold_plddt_cutoff: float = 50.0,
    skip_esmfold: bool = False,
) -> dict[str, Any]:
    """Unified two-tier structure prediction.

    1. ESMFold (fast pre-filter) - if mean_pLDDT < esmfold_plddt_cutoff, discard.
    2. Boltz-1 (AF3-class) - with ColabFold fallback for RNA complex.

    Returns merged result dict with all keys from both tiers.
    """
    result: dict[str, Any] = {
        "design_id": design_id,
        "esmfold_plddt": None,
        "esmfold_pass": None,
        "boltz_pdb": None,
        "boltz_plddt": None,
        "boltz_rna_fallback": False,
        "final_tool": None,
        "final_pdb": None,
        "final_mean_plddt": None,
        "success": False,
    }

    # Tier 1: ESMFold pre-filter
    if not skip_esmfold and _check_esmfold():
        esm_result = predict_structure_esmfold(sequence, design_id, output_dir / "esmfold")
        result["esmfold_plddt"] = esm_result.get("mean_plddt")
        if esm_result.get("success") and esm_result["mean_plddt"] is not None:
            if esm_result["mean_plddt"] < esmfold_plddt_cutoff:
                result["esmfold_pass"] = False
                result["esmfold_discard_reason"] = (
                    f"ESMFold mean pLDDT {esm_result['mean_plddt']:.1f} < {esmfold_plddt_cutoff}"
                )
                return result  # discard before expensive Boltz-1 run
            result["esmfold_pass"] = True
    else:
        result["esmfold_pass"] = None  # not evaluated (tool absent)

    # Tier 2: Boltz-1 (with ColabFold fallback for RNA complex)
    boltz_result = predict_structure_boltz(
        sequence=sequence,
        design_id=design_id,
        output_dir=output_dir / "boltz",
        brna_sequence=brna_sequence,
    )
    result["boltz_pdb"] = boltz_result.get("pdb_path")
    result["boltz_plddt"] = boltz_result.get("mean_plddt")
    result["boltz_rna_fallback"] = boltz_result.get("boltz_rna_fallback", False)
    result["final_tool"] = boltz_result.get("tool")
    result["final_pdb"] = boltz_result.get("pdb_path")
    result["final_mean_plddt"] = boltz_result.get("mean_plddt")
    result["success"] = boltz_result.get("success", False)

    return result


# ---------------------------------------------------------------------------
# Batch prediction (used by driver script)
# ---------------------------------------------------------------------------


def run_batch_prediction(
    designs_df: pd.DataFrame,
    output_dir: Path,
    esmfold_plddt_cutoff: float = 50.0,
    brna_col: str | None = None,
    sequence_col: str = "protein_sequence",
    id_col: str = "design_id",
    fallback_log_path: Path | None = None,
) -> pd.DataFrame:
    """Predict structures for all designs in DataFrame. Adds structure columns in-place.

    Columns added: esmfold_plddt, esmfold_pass, boltz_pdb, boltz_plddt,
    boltz_rna_fallback, final_pdb, final_mean_plddt, structure_success.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []

    for i, row in designs_df.iterrows():
        design_id = row[id_col]
        sequence = row[sequence_col]
        brna_seq = row.get(brna_col) if brna_col else None

        print(f"  [{i + 1}/{len(designs_df)}] {design_id} ({len(sequence)} aa) ...")
        r = predict_structure(
            sequence=sequence,
            design_id=design_id,
            output_dir=output_dir,
            brna_sequence=brna_seq,
            esmfold_plddt_cutoff=esmfold_plddt_cutoff,
        )
        results.append(r)

    results_df = pd.DataFrame(results)

    # Log fallback cases
    fallback_df = results_df[results_df["boltz_rna_fallback"] == True]  # noqa: E712
    if len(fallback_df) > 0:
        log_path = fallback_log_path or (output_dir / "structure_fallback_log.parquet")
        fallback_df.to_parquet(log_path, index=False)
        print(f"\n  Boltz-1 RNA fallback applied: {len(fallback_df)} designs -> {log_path}")

    # Merge back onto input DataFrame
    out = designs_df.copy()
    for col in [
        "esmfold_plddt",
        "esmfold_pass",
        "boltz_pdb",
        "boltz_plddt",
        "boltz_rna_fallback",
        "final_pdb",
        "final_mean_plddt",
        "success",
    ]:
        if col in results_df.columns:
            out[col] = results_df[col].values

    return out
