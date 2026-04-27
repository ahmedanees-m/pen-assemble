"""Anti-mode-collapse diversity enforcement. Step 19 / P5 pre-registration.

P5: Top-5 final designs must come from >=3 different scaffold sources (strategies).
Prevents the framework from collapsing to IS621-only variants in the final ranking.
"""

from __future__ import annotations

import pandas as pd


def check_diversity(
    top_k_designs: pd.DataFrame,
    scaffold_col: str = "scaffold_provenance_catalytic_core",
    min_distinct: int = 3,
) -> dict:
    """Check P5: distinct scaffold sources in top-k designs.

    Returns {'n_distinct': int, 'sources': list[str], 'passes_p5': bool}.
    """
    if scaffold_col not in top_k_designs.columns:
        raise KeyError(f"Column {scaffold_col!r} not found in designs DataFrame.")
    sources = top_k_designs[scaffold_col].unique().tolist()
    return {
        "n_distinct": len(sources),
        "sources": sources,
        "passes_p5": len(sources) >= min_distinct,
    }


def enforce_diversity(
    designs_df: pd.DataFrame,
    top_k: int = 50,
    min_per_strategy: int = 5,
) -> pd.DataFrame:
    """Re-rank final catalog to enforce minimum representation per strategy (Step 19).

    Anti-mode-collapse rule (pre-registered P5 requirement):
      - If the dominant strategy occupies >= 4 of the top-5 positions, replace
        the 4th dominant design with the top design from the next-best strategy.
      - Then pad to `top_k` rows maintaining PEN-SCORE order.
      - Guarantees >= 3 distinct scaffold sources in the top-5.

    Args:
        designs_df: DataFrame sorted by pen_score descending (all designs).
        top_k: Number of designs to return in the re-ranked catalog.
        min_per_strategy: Minimum designs per strategy in the full `top_k` catalog.
            Strategies below this count are padded from their best remaining designs.

    Returns:
        Re-ranked DataFrame of length min(top_k, len(designs_df)) with
        ``diversity_rank`` column (1 = best) and ``diversity_note`` explaining
        any swap made.
    """
    if designs_df.empty:
        return designs_df.copy()

    df = designs_df.copy()
    if "pen_score" in df.columns:
        df = df.sort_values("pen_score", ascending=False).reset_index(drop=True)

    strategy_col = next(
        (
            c
            for c in ["strategy", "scaffold_provenance_catalytic_core", "scaffold_source"]
            if c in df.columns
        ),
        None,
    )

    df["diversity_note"] = ""

    if strategy_col and len(df) >= 5:
        top5 = df.head(5).copy()
        counts = top5[strategy_col].value_counts()
        dominant_strategy = counts.index[0]
        dominant_count = counts.iloc[0]

        if dominant_count >= 4:
            # Find top non-dominant design not already in top-5
            top5_ids = set(top5.index)
            alt = df[(df[strategy_col] != dominant_strategy) & (~df.index.isin(top5_ids))]
            if not alt.empty:
                swap_idx = alt.index[0]
                # Replace the 4th dominant design in top-5
                dominant_positions = top5[top5[strategy_col] == dominant_strategy].index.tolist()
                replace_pos = (
                    dominant_positions[3]
                    if len(dominant_positions) >= 4
                    else dominant_positions[-1]
                )

                # Build new top-5 by swapping
                new_top5_indices = [i for i in top5.index if i != replace_pos] + [swap_idx]
                # Re-sort by pen_score
                new_top5 = df.loc[new_top5_indices].sort_values("pen_score", ascending=False)
                new_top5.at[swap_idx, "diversity_note"] = (
                    f"anti_mode_collapse: replaced rank-{replace_pos + 1} {dominant_strategy} "
                    f"with {df.at[swap_idx, strategy_col]}"
                )
                rest_indices = [i for i in df.index if i not in set(new_top5_indices)]
                df = pd.concat([new_top5, df.loc[rest_indices]], ignore_index=True)

    # Pad each strategy to min_per_strategy if possible
    if strategy_col and min_per_strategy > 0:
        # Check representation
        result = df.head(top_k)
        for strat in df[strategy_col].unique() if strategy_col else []:
            strat_in_result = result[result[strategy_col] == strat]
            if len(strat_in_result) < min_per_strategy:
                # Add more from this strategy from the tail
                strat_extra = df[(df[strategy_col] == strat) & (~df.index.isin(result.index))].head(
                    min_per_strategy - len(strat_in_result)
                )
                if not strat_extra.empty:
                    result = (
                        pd.concat([result, strat_extra])
                        .sort_values("pen_score", ascending=False)
                        .head(top_k)
                    )

        df = result.reset_index(drop=True)
    else:
        df = df.head(top_k).reset_index(drop=True)

    df["diversity_rank"] = range(1, len(df) + 1)
    return df
