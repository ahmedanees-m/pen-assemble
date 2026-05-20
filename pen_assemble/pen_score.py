"""
PenScore: composite scoring formula for IS110-family bridge recombinase designs.

The PenScore (Programmable Editor Nomination score) is a weighted linear combination
of seven mechanistic axes calibrated for the human_therapeutic_aav_insertion use-case.

References
----------
- PEN-ASSEMBLE v0.5.0 execution plan (Paper 4)
- IS621 published lockpoint: pen_score = 0.929 (verbatim pre-registered)
- MHCflurry 2.2.1-calibrated lockpoint: 0.9255 (secondary analysis)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    "WEIGHTS",
    "IS621_LOCKPOINT",
    "IS621_LOCKPOINT_CALIBRATED",
    "PenScoreAxes",
    "pen_score",
    "beats_is621",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Axis weights for human_therapeutic_aav_insertion use-case.
WEIGHTS: dict[str, float] = {
    "S_DSB":    0.25,
    "S_Spec":   0.10,
    "S_Cargo":  0.20,
    "S_Deliv":  0.15,
    "S_Immuno": 0.10,
    "S_Prog":   0.15,
    "S_Mature": 0.05,
}

#: Verbatim pre-registered IS621 lockpoint (primary threshold).
IS621_LOCKPOINT: float = 0.929

#: MHCflurry 2.2.1-calibrated lockpoint (secondary analysis only).
IS621_LOCKPOINT_CALIBRATED: float = 0.9255


# ---------------------------------------------------------------------------
# Data class for axis scores
# ---------------------------------------------------------------------------

@dataclass
class PenScoreAxes:
    """Container for the seven PenScore axis scores.

    All axes are in [0, 1]. Missing values (``None``) are treated as 0.0
    when computing the composite score.

    Parameters
    ----------
    S_DSB:
        Double-strand break induction score.
    S_Spec:
        Target-site specificity score.
    S_Cargo:
        Payload compatibility score (IS110-family = 1.0 by mechanism).
    S_Deliv:
        Delivery suitability score (AAV packaging compatibility).
    S_Immuno:
        De-immunization score (1 − normalised predicted immunogenicity).
    S_Prog:
        Programmability score (bRNA re-targeting tractability).
    S_Mature:
        Maturity / technology-readiness score.
    """
    S_DSB:    float = 0.0
    S_Spec:   float = 0.0
    S_Cargo:  float = 0.0
    S_Deliv:  float = 0.0
    S_Immuno: float = 0.0
    S_Prog:   float = 0.0
    S_Mature: float = 0.0

    def as_dict(self) -> dict[str, float]:
        """Return axes as a plain dict."""
        return {
            "S_DSB":    self.S_DSB,
            "S_Spec":   self.S_Spec,
            "S_Cargo":  self.S_Cargo,
            "S_Deliv":  self.S_Deliv,
            "S_Immuno": self.S_Immuno,
            "S_Prog":   self.S_Prog,
            "S_Mature": self.S_Mature,
        }

    def contributions(self) -> dict[str, float]:
        """Return per-axis weighted contributions."""
        return {ax: val * WEIGHTS[ax] for ax, val in self.as_dict().items()}


# ---------------------------------------------------------------------------
# Core formula
# ---------------------------------------------------------------------------

def pen_score(axes: PenScoreAxes) -> float:
    """Compute the composite PenScore for a single design.

    Parameters
    ----------
    axes:
        Seven axis scores packaged as a :class:`PenScoreAxes` instance.

    Returns
    -------
    float
        Composite PenScore in [0, 1].

    Examples
    --------
    >>> from pen_assemble.pen_score import pen_score, PenScoreAxes
    >>> ax = PenScoreAxes(S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0,
    ...                   S_Deliv=1.0, S_Immuno=0.8777, S_Prog=1.0, S_Mature=0.5)
    >>> round(pen_score(ax), 4)
    0.9628
    """
    total = sum(getattr(axes, ax) * w for ax, w in WEIGHTS.items())
    return float(total)


def beats_is621(score: float, calibrated: bool = False) -> bool:
    """Return True if *score* exceeds the IS621 reference lockpoint.

    Parameters
    ----------
    score:
        A PenScore value.
    calibrated:
        If ``True``, compare against the MHCflurry 2.2.1-calibrated lockpoint
        (0.9255). Default is ``False`` (verbatim pre-registered 0.929).

    Returns
    -------
    bool

    Examples
    --------
    >>> from pen_assemble.pen_score import beats_is621
    >>> beats_is621(0.935)
    True
    >>> beats_is621(0.928)
    False
    >>> beats_is621(0.928, calibrated=True)
    True
    """
    threshold = IS621_LOCKPOINT_CALIBRATED if calibrated else IS621_LOCKPOINT
    return score > threshold
