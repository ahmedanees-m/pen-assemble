"""
Catalog loading utilities for PEN-ASSEMBLE release data.

Functions load the release CSV/Parquet files into pandas DataFrames and
provide convenience views (P1 beaters, P5 top-5, strategy subsets).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

__all__ = [
    "load_catalog",
    "load_p1_beaters",
    "load_top5",
    "filter_strategy",
    "RELEASE_DIR",
]

# Default release directory (relative to package root)
_PKG_ROOT = Path(__file__).resolve().parent.parent
RELEASE_DIR: Path = _PKG_ROOT / "catalog" / "release_v0.5.0"

# The frozen v0.5.0 release catalog ships with the repository, not with the wheel:
# it is the pre-registration record, so it is not duplicated into the package.
_RELEASE_HINT = (
    "The installed wheel ships the library, not the frozen v0.5.0 release catalog. "
    "Clone the repository, or pass release_dir= pointing at a downloaded "
    "catalog/release_v0.5.0 directory."
)

#: IS621 verbatim lockpoint used for P1 beater definition.
IS621_LOCKPOINT = 0.929


def load_catalog(
    release_dir: str | Path | None = None,
    fmt: str = "parquet",
) -> pd.DataFrame:
    """Load the full 1,029-design PEN-ASSEMBLE scorecard.

    Parameters
    ----------
    release_dir:
        Path to the release directory. Defaults to
        ``catalog/release_v0.5.0/`` relative to the repository root.
    fmt:
        ``"parquet"`` (default) or ``"csv"``.

    Returns
    -------
    pd.DataFrame
        Scorecard with all catalog columns plus ``protein_sequence``.

    Raises
    ------
    FileNotFoundError
        If the release directory or catalog file does not exist.

    Examples
    --------
    >>> from pen_assemble.catalog import load_catalog
    >>> df = load_catalog()
    >>> len(df)
    1029
    >>> "pen_score" in df.columns
    True
    """
    rdir = Path(release_dir) if release_dir else RELEASE_DIR
    if not rdir.exists():
        raise FileNotFoundError(
            f"Release catalog not found: {rdir}\n{_RELEASE_HINT}\n"
            "From a source checkout you can regenerate it with scripts/50_assemble_catalog.py."
        )
    if fmt == "parquet":
        fpath = rdir / "pen_assemble_catalog.parquet"
    else:
        fpath = rdir / "pen_assemble_catalog.csv"

    if not fpath.exists():
        raise FileNotFoundError(f"Catalog file not found: {fpath}\n{_RELEASE_HINT}")

    return pd.read_parquet(fpath) if fmt == "parquet" else pd.read_csv(fpath)


def load_p1_beaters(
    release_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Load the 16 designs that beat the IS621 verbatim lockpoint (pen_score > 0.929).

    Parameters
    ----------
    release_dir:
        Path to the release directory.

    Returns
    -------
    pd.DataFrame
        Sorted by pen_score descending.

    Examples
    --------
    >>> from pen_assemble.catalog import load_p1_beaters
    >>> p1 = load_p1_beaters()
    >>> len(p1)
    16
    >>> (p1["pen_score"] > 0.929).all()
    True
    """
    rdir = Path(release_dir) if release_dir else RELEASE_DIR
    fpath = rdir / "p1_beaters_catalog.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"P1 beaters file not found: {fpath}\n{_RELEASE_HINT}")
    df = pd.read_csv(fpath)
    return df.sort_values("pen_score", ascending=False).reset_index(drop=True)


def load_top5(
    release_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Load the five P5-compliant top designs (diversity-enforced top-5).

    Note: the top-5 is not purely rank-ordered - rank-5 was diversity-enforced
    (A_007 replaces natural rank-5 D023). See validation/P5_diversity_result.json.

    Parameters
    ----------
    release_dir:
        Path to the release directory.

    Returns
    -------
    pd.DataFrame
        Five designs sorted by pen_score descending.

    Examples
    --------
    >>> from pen_assemble.catalog import load_top5
    >>> t5 = load_top5()
    >>> len(t5)
    5
    >>> t5["strategy"].nunique() >= 3
    True
    """
    rdir = Path(release_dir) if release_dir else RELEASE_DIR
    fpath = rdir / "p5_top5_catalog.csv"
    if not fpath.exists():
        raise FileNotFoundError(f"Top-5 file not found: {fpath}\n{_RELEASE_HINT}")
    df = pd.read_csv(fpath)
    return df.sort_values("pen_score", ascending=False).reset_index(drop=True)


def filter_strategy(
    df: pd.DataFrame,
    strategy: str,
) -> pd.DataFrame:
    """Filter a catalog DataFrame to a single strategy.

    Parameters
    ----------
    df:
        Catalog DataFrame (from :func:`load_catalog` or similar).
    strategy:
        One of ``"A"``, ``"B"``, ``"C"``, ``"D"``.

    Returns
    -------
    pd.DataFrame
        Filtered view, sorted by pen_score descending.

    Raises
    ------
    ValueError
        If strategy is not A/B/C/D or not present in ``df``.

    Examples
    --------
    >>> from pen_assemble.catalog import load_catalog, filter_strategy
    >>> df = load_catalog()
    >>> filter_strategy(df, "C")["design_id"].tolist()
    ['IS621_deimmunized_v2_Y255K_D41C_G65K_E187M_D27C_D203C_V285C_V152C_T224C_L318K_E87C_L193I_L275I_P177V', 'C_targeted_001']
    """
    valid = {"A", "B", "C", "D"}
    if strategy not in valid:
        raise ValueError(f"strategy must be one of {valid}; got {strategy!r}")
    sub = df[df["strategy"] == strategy].copy()
    if sub.empty:
        raise ValueError(f"No designs found for strategy {strategy!r}")
    return sub.sort_values("pen_score", ascending=False).reset_index(drop=True)
