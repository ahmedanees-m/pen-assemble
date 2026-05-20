"""MECH-CLASS re-evaluation of all designs. Step 15.

Requires mech-class>=0.5.1 (pip install 'mech-class>=0.5.1').

Critical v0.5.1 gate:
  - SpCas9 must NOT trigger composite_flag (was bug in v0.5.0)
  - Any chimera with active Cas9 nuclease fires DSB_NUCLEASE Tier A → discard
  - Strategy A dCas9 modules MUST use D10A+H840A double mutation

All designs must pass: tier_a == 'DSB_FREE_TRANSEST_RECOMBINASE'.
Designs with tier_a == 'UNCLASSIFIED' are flagged for manual review.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# MECH-CLASS result dataclass
# ---------------------------------------------------------------------------

@dataclass
class MechClassResult:
    """Result from a single MECH-CLASS evaluation."""
    design_id: str
    tier_a: str                         # DSB_FREE_TRANSEST_RECOMBINASE | DSB_NUCLEASE | UNCLASSIFIED
    tier_a_confidence: float            # 0–1
    composite: bool                     # IS110-class composite architecture detected
    composite_prob: float               # 0–1
    mechanism_label: str                # human-readable label
    pfam_hits: list[str]               # Pfam domains detected
    mech_class_version: str = "0.5.1"
    error: Optional[str] = None

    def is_dsb_free(self) -> bool:
        return self.tier_a == "DSB_FREE_TRANSEST_RECOMBINASE"

    def to_dict(self) -> dict[str, Any]:
        return {
            "design_id": self.design_id,
            "tier_a": self.tier_a,
            "tier_a_confidence": self.tier_a_confidence,
            "composite": self.composite,
            "composite_prob": self.composite_prob,
            "mechanism_label": self.mechanism_label,
            "pfam_hits": "|".join(self.pfam_hits),
            "mech_class_version": self.mech_class_version,
            "dsb_free_pass": self.is_dsb_free(),
            "mech_class_error": self.error,
        }


# ---------------------------------------------------------------------------
# MECH-CLASS invocation
# ---------------------------------------------------------------------------

def _check_mechclass() -> bool:
    try:
        import mechclass  # noqa: F401
        return True
    except ImportError:
        return False


def _check_mechclass_version() -> Optional[str]:
    try:
        import importlib.metadata
        return importlib.metadata.version("mech-class")
    except Exception:
        return None


def run_mech_class_on_sequence(
    sequence: str,
    design_id: str = "query",
    pfam_hmm_path: Optional[Path] = None,
) -> MechClassResult:
    """Run mech-class v0.5.1 Predictor on a single sequence.

    Requires: pip install 'mech-class>=0.5.1'
    """
    if not _check_mechclass():
        raise ImportError(
            "mech-class not installed. Run: pip install 'mech-class>=0.5.1'\n"
            "Or use the pen-stack/design Docker image which includes mech-class."
        )

    version = _check_mechclass_version() or "unknown"

    try:
        from mechclass import Predictor

        predictor = Predictor(pfam_db=str(pfam_hmm_path) if pfam_hmm_path else None)
        prediction = predictor.predict(sequence)

        return MechClassResult(
            design_id=design_id,
            tier_a=prediction.tier_a,
            tier_a_confidence=prediction.tier_a_confidence,
            composite=prediction.composite,
            composite_prob=prediction.composite_prob,
            mechanism_label=getattr(prediction, "mechanism_label", prediction.tier_a),
            pfam_hits=getattr(prediction, "pfam_hits", []),
            mech_class_version=version,
        )
    except AttributeError:
        # Older API: try mech_class.Classifier
        try:
            from mech_class import Classifier
            clf = Classifier()
            result = clf.classify(sequence)
            return MechClassResult(
                design_id=design_id,
                tier_a=result.get("tier_a", "UNCLASSIFIED"),
                tier_a_confidence=result.get("confidence", 0.0),
                composite=result.get("composite", False),
                composite_prob=result.get("composite_prob", 0.0),
                mechanism_label=result.get("label", "UNCLASSIFIED"),
                pfam_hits=result.get("pfam_hits", []),
                mech_class_version=version,
            )
        except Exception as e2:
            return MechClassResult(
                design_id=design_id,
                tier_a="UNCLASSIFIED",
                tier_a_confidence=0.0,
                composite=False,
                composite_prob=0.0,
                mechanism_label="UNCLASSIFIED",
                pfam_hits=[],
                mech_class_version=version,
                error=str(e2),
            )
    except Exception as e:
        return MechClassResult(
            design_id=design_id,
            tier_a="UNCLASSIFIED",
            tier_a_confidence=0.0,
            composite=False,
            composite_prob=0.0,
            mechanism_label="UNCLASSIFIED",
            pfam_hits=[],
            mech_class_version=version,
            error=str(e),
        )


def evaluate_design_mechanism(
    accession: str,
    sequence: str,
    pfam_hints: Optional[list[str]] = None,
    pfam_hmm_path: Optional[Path] = None,
) -> dict[str, Any]:
    """Run mech-class on a single design. Returns dict with all MECH-CLASS fields.

    Parameters
    ----------
    accession : str
        Design ID (e.g. "A_001", "B_003_ACC12345").
    sequence : str
        Protein sequence (1-letter code).
    pfam_hints : list[str] | None
        Pre-computed Pfam hits to pass as prior (skips internal HMMER scan if provided).
    pfam_hmm_path : Path | None
        Path to Pfam-A.hmm for internal HMMER scan (optional).
    """
    result = run_mech_class_on_sequence(sequence, design_id=accession, pfam_hmm_path=pfam_hmm_path)
    return result.to_dict()


def filter_dsb_free(
    designs_df: pd.DataFrame,
    tier_a_col: str = "tier_a",
    confidence_col: str = "tier_a_confidence",
    min_confidence: float = 0.80,
) -> pd.DataFrame:
    """Keep only designs with tier_a == 'DSB_FREE_TRANSEST_RECOMBINASE' and confidence >= min_confidence.

    Designs with tier_a == 'UNCLASSIFIED' are flagged in a separate column for manual review.
    """
    n_before = len(designs_df)

    unclassified_mask = designs_df[tier_a_col] == "UNCLASSIFIED"
    dsb_nuclease_mask = designs_df[tier_a_col] == "DSB_NUCLEASE"
    pass_mask = (
        (designs_df[tier_a_col] == "DSB_FREE_TRANSEST_RECOMBINASE")
        & (designs_df[confidence_col] >= min_confidence)
    )

    designs_df = designs_df.copy()
    designs_df["mech_class_flag"] = "PASS"
    designs_df.loc[unclassified_mask, "mech_class_flag"] = "UNCLASSIFIED_MANUAL_REVIEW"
    designs_df.loc[dsb_nuclease_mask, "mech_class_flag"] = "DSB_NUCLEASE_DISCARD"

    result = designs_df[pass_mask].copy()
    n_discarded_nuclease = dsb_nuclease_mask.sum()
    n_manual = unclassified_mask.sum()
    print(
        f"  MECH-CLASS filter: {n_before} -> {len(result)} DSB-free "
        f"(discarded {n_discarded_nuclease} DSB_NUCLEASE, "
        f"{n_manual} UNCLASSIFIED for manual review)"
    )
    return result


def run_mechclass_batch(
    designs_df: pd.DataFrame,
    sequence_col: str = "protein_sequence",
    id_col: str = "design_id",
    pfam_hmm_path: Optional[Path] = None,
) -> pd.DataFrame:
    """Run MECH-CLASS on all designs in DataFrame. Adds MECH-CLASS columns in-place."""
    rows: list[dict] = []
    for _, row in designs_df.iterrows():
        r = run_mech_class_on_sequence(
            row[sequence_col], design_id=row[id_col], pfam_hmm_path=pfam_hmm_path
        )
        rows.append(r.to_dict())

    mc_df = pd.DataFrame(rows)
    out = designs_df.copy()
    for col in mc_df.columns:
        if col != id_col:
            out[col] = mc_df[col].values
    return out
