"""
Tests for pen_assemble.catalog
"""
import pytest
import pandas as pd
from pen_assemble.catalog import (
    load_catalog,
    load_p1_beaters,
    load_top5,
    filter_strategy,
    IS621_LOCKPOINT,
)

# These tests require the release catalog to exist.
# Mark the whole module to skip gracefully if catalog is absent.
pytestmark = pytest.mark.skipif(
    not __import__("pathlib").Path(
        __import__("pen_assemble.catalog", fromlist=["RELEASE_DIR"]).RELEASE_DIR
    ).exists(),
    reason="Release catalog not built (run 50_assemble_catalog.py first)",
)


class TestLoadCatalog:
    def test_row_count(self):
        df = load_catalog()
        assert len(df) == 1029

    def test_required_columns(self):
        df = load_catalog()
        required = {"design_id", "strategy", "pen_score", "protein_sequence",
                    "S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno",
                    "S_Prog", "S_Mature", "beats_is621"}
        assert required.issubset(set(df.columns))

    def test_pen_score_range(self):
        df = load_catalog()
        assert df["pen_score"].between(0.0, 1.0).all()

    def test_strategies_present(self):
        df = load_catalog()
        assert set(df["strategy"].unique()) == {"A", "B", "C", "D"}

    def test_no_null_sequences(self):
        df = load_catalog()
        assert df["protein_sequence"].notna().all()

    def test_csv_format(self):
        df = load_catalog(fmt="csv")
        assert len(df) == 1029

    def test_beats_is621_consistent(self):
        # beats_is621 column was computed against the MHCflurry 2.2.1-calibrated
        # lockpoint (0.9255) in the upstream pipeline, not the verbatim 0.929.
        # All 32 D-designs (pen_score 0.9261-0.9353) are marked True.
        # Verify: every verbatim beater (pen_score > 0.929) has beats_is621 = True.
        df = load_catalog()
        verbatim_beaters = df["pen_score"] > IS621_LOCKPOINT
        assert df.loc[verbatim_beaters, "beats_is621"].astype(bool).all(), (
            "All designs with pen_score > 0.929 must have beats_is621 = True"
        )


class TestLoadP1Beaters:
    def test_count(self):
        p1 = load_p1_beaters()
        assert len(p1) == 16

    def test_all_beat_is621(self):
        p1 = load_p1_beaters()
        assert (p1["pen_score"] > IS621_LOCKPOINT).all()

    def test_sorted_descending(self):
        p1 = load_p1_beaters()
        assert list(p1["pen_score"]) == sorted(p1["pen_score"], reverse=True)


class TestLoadTop5:
    def test_count(self):
        t5 = load_top5()
        assert len(t5) == 5

    def test_strategy_diversity(self):
        t5 = load_top5()
        assert t5["strategy"].nunique() >= 3

    def test_sorted_descending(self):
        t5 = load_top5()
        assert list(t5["pen_score"]) == sorted(t5["pen_score"], reverse=True)

    def test_top_design_is_strategy_c(self):
        """Top-ranked design should be IS621_deimmunized_v2 (Strategy C)."""
        t5 = load_top5()
        assert t5.iloc[0]["strategy"] == "C"


class TestFilterStrategy:
    def test_strategy_b_count(self):
        df = load_catalog()
        sub = filter_strategy(df, "B")
        assert len(sub) == 992

    def test_strategy_c_count(self):
        df = load_catalog()
        sub = filter_strategy(df, "C")
        assert len(sub) == 2

    def test_strategy_a_count(self):
        df = load_catalog()
        sub = filter_strategy(df, "A")
        # 15 sourced; 12 excluded at gate_7; 3 triaged
        assert len(sub) == 3

    def test_strategy_d_count(self):
        df = load_catalog()
        sub = filter_strategy(df, "D")
        # 32 = 30 ProtMPNN + 2 natural orthologs
        assert len(sub) == 32

    def test_sorted_descending(self):
        df = load_catalog()
        sub = filter_strategy(df, "D")
        assert list(sub["pen_score"]) == sorted(sub["pen_score"], reverse=True)

    def test_invalid_strategy_raises(self):
        df = load_catalog()
        with pytest.raises(ValueError, match="strategy must be one of"):
            filter_strategy(df, "Z")

    def test_returns_dataframe(self):
        df = load_catalog()
        sub = filter_strategy(df, "C")
        assert isinstance(sub, pd.DataFrame)
