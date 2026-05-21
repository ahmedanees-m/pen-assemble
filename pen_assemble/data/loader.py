"""Load scaffold universe and pre-registration configs with Pydantic validation."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel

DATA_DIR = Path(__file__).parent


class ScaffoldEntry(BaseModel):
    """A single scaffold editor entry from scaffold_universe.yaml."""

    id: str
    canonical_accession: str
    canonical_pdb: str | None = None
    paper3_pen_score: float | None = None  # null for carve-out sentinels (MmeFz2)
    tier_b: str
    is_composite: bool = False
    contributes_modules: list[str] = []
    contributes_modules_carveout_only: list[str] = []
    canonical_protospacer: str | None = None
    cargo_capacity_bp: int | None = None
    paper3_axes: dict = {}
    paper3_pen_score_breakdown: dict = {}
    total_aa: int | None = None


class PreRegistrationPrediction(BaseModel):
    """A single pre-registered prediction entry."""

    id: str
    statement: str
    operationalization: str
    threshold: str
    test_script: str
    rationale: str = ""
    failure_mode_analysis: str = ""


def load_scaffold_universe() -> list[ScaffoldEntry]:
    """Load and validate scaffold_universe.yaml."""
    path = DATA_DIR / "scaffold_universe.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"scaffold_universe.yaml not found at {path}. "
            "Run scripts/01_curate_scaffold_universe.py first (Step 4)."
        )
    data = yaml.safe_load(path.read_text())
    return [ScaffoldEntry(**e) for e in data["scaffolds"]]


def load_pre_registration() -> list[PreRegistrationPrediction]:
    """Load and validate pre_registration.yaml."""
    path = DATA_DIR / "pre_registration.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"pre_registration.yaml not found at {path}. "
            "Run scripts/02_compute_prereg_hash.py first (Step 7)."
        )
    data = yaml.safe_load(path.read_text())
    return [PreRegistrationPrediction(**p) for p in data["predictions"]]


def load_design_strategies() -> dict:
    """Load design_strategies.yaml as a raw dict."""
    path = DATA_DIR / "design_strategies.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"design_strategies.yaml not found at {path}. Complete Step 5 first."
        )
    return yaml.safe_load(path.read_text())
