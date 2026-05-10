"""
Tests for pen-assemble v0.5.2 catalog - PEN-COMPARE v3.2 compatibility.

Verifies that catalog_v0.5.2_current.parquet contains the two new boolean
columns required by PEN-COMPARE v3.2:
  • intrinsic_cargo_mechanism - True for all IS110-family designs (by PFAM triage)
  • cell_based_evidence       - False for all designs (computational; no wet-lab)

Also verifies ISCro4 canonical naming (no IS622 in parent_editor).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

CATALOG_V052 = Path(__file__).resolve().parents[1] / "data" / "catalog_v0.5.2_current.parquet"

pytestmark = pytest.mark.skipif(
    not CATALOG_V052.exists(),
    reason=(
        "catalog_v0.5.2_current.parquet not found. Run: python scripts/upgrade_catalog_to_v052.py"
    ),
)


@pytest.fixture(scope="module")
def cat() -> pd.DataFrame:
    return pd.read_parquet(CATALOG_V052)


# Existence and row count


def test_catalog_loads(cat: pd.DataFrame) -> None:
    """Catalog must load and have exactly 1 029 designs."""
    assert len(cat) == 1029, f"Expected 1 029 rows, got {len(cat)}"


# New v3.2 columns present


def test_intrinsic_cargo_mechanism_column_present(cat: pd.DataFrame) -> None:
    assert "intrinsic_cargo_mechanism" in cat.columns, (
        "Column 'intrinsic_cargo_mechanism' missing - run upgrade_catalog_to_v052.py"
    )


def test_cell_based_evidence_column_present(cat: pd.DataFrame) -> None:
    assert "cell_based_evidence" in cat.columns, (
        "Column 'cell_based_evidence' missing - run upgrade_catalog_to_v052.py"
    )


# PEN-COMPARE v3.2 invariants


def test_all_designs_intrinsic_cargo_true(cat: pd.DataFrame) -> None:
    """All IS110-family designs have intrinsic_cargo_mechanism = True.

    Invariant: every design passed the IS110-family PFAM triage gate
    (PF01548 AND PF02371 present). IS110 bridge recombinases insert cargo
    via their catalytic mechanism - no external HDR donor template needed.
    """
    n_false = int((~cat["intrinsic_cargo_mechanism"]).sum())
    assert n_false == 0, (
        f"{n_false} designs have intrinsic_cargo_mechanism=False. "
        "All catalog designs are IS110-family (passed PFAM triage) - "
        "all must have intrinsic_cargo_mechanism=True."
    )


def test_all_designs_no_cell_based_evidence(cat: pd.DataFrame) -> None:
    """PEN-COMPARE v3.2 P2 invariant: zero designs have cell_based_evidence = True.

    All designs are computational predictions. None have peer-reviewed mammalian
    cell data at > 1% efficiency. This ensures no pen-assemble design can be
    classified TRUE_WRITER by PEN-COMPARE (all cap at PROBABLE_WRITER).
    """
    n_true = int(cat["cell_based_evidence"].sum())
    assert n_true == 0, (
        f"{n_true} designs have cell_based_evidence=True. "
        "All designs are computational - none have wet-lab cell data. "
        "This violates the PEN-COMPARE v3.2 P2 pre-registered prediction."
    )


def test_no_is622_in_parent_editor(cat: pd.DataFrame) -> None:
    """All IS622 parent references must be migrated to ISCro4.

    IS622 is the deprecated preprint label (Perry et al. 2025 bioRxiv).
    The canonical name is ISCro4 (Pelea 2026 Science + UniProt D2TGM5).
    """
    n_is622 = int((cat["parent_editor"] == "IS622").sum())
    assert n_is622 == 0, (
        f"{n_is622} designs still have parent_editor='IS622'. "
        "Run upgrade_catalog_to_v052.py to migrate to 'ISCro4'."
    )


# IS621-derived design parent assignment


def test_strategy_c_parent_is_is621(cat: pd.DataFrame) -> None:
    """Strategy C designs (deimmunised IS621 variants) have parent_editor='IS621'."""
    strategy_c = cat[cat["strategy"] == "C"]
    if len(strategy_c) > 0:
        assert (strategy_c["parent_editor"] == "IS621").all(), (
            f"Strategy C parent_editor values: {strategy_c['parent_editor'].unique()}"
        )


def test_strategy_d_protmpnn_parent_is_is621(cat: pd.DataFrame) -> None:
    """Strategy D ProteinMPNN redesigns (D0XX_IS621_ProtMPNN_...) have parent='IS621'."""
    strategy_d = cat[cat["strategy"] == "D"]
    protmpnn_mask = strategy_d["design_id"].str.match(r"^D\d+_IS621_ProtMPNN_", na=False)
    protmpnn_designs = strategy_d[protmpnn_mask]
    if len(protmpnn_designs) > 0:
        assert (protmpnn_designs["parent_editor"] == "IS621").all(), (
            f"ProteinMPNN redesign parent_editor values: {protmpnn_designs['parent_editor'].unique()}"
        )


# Backward-compatibility: pre-registration results unchanged


def test_p1_beaters_count_unchanged(cat: pd.DataFrame) -> None:
    """P1 pre-registration result must not change: exactly 16 designs > 0.929."""
    n = int(cat["beats_v010_lockpoint_0929"].sum())
    assert n == 16, (
        f"P1 beater count changed: expected 16, got {n}. "
        "The v0.5.2 upgrade must not alter any PenScore columns."
    )


def test_total_row_count_preserved(cat: pd.DataFrame) -> None:
    """All 1 029 designs must be preserved - no rows added or removed."""
    assert len(cat) == 1029


def test_column_count_increased_by_two(cat: pd.DataFrame) -> None:
    """v0.5.2 adds exactly 2 columns to the v0.5.1 catalog (25 -> 27)."""
    # v0.5.1 had 25 columns; v0.5.2 adds intrinsic_cargo_mechanism + cell_based_evidence
    assert len(cat.columns) == 27, (
        f"Expected 27 columns (25 from v0.5.1 + 2 new), got {len(cat.columns)}: "
        f"{cat.columns.tolist()}"
    )
