"""Linker design helpers for Strategy A domain-swap chimeras."""
from __future__ import annotations

from typing import Literal

LinkerType = Literal["GS", "EAAAK", "natural"]

LINKER_PRESETS: dict[str, dict] = {
    "GS_4": {"sequence": "GGGGS", "type": "GS", "flexible": True},
    "GS_8": {"sequence": "GGGGSGGGGS", "type": "GS", "flexible": True},
    "EAAAK_5": {"sequence": "EAAAK", "type": "EAAAK", "rigid": True},
    "EAAAK_10": {"sequence": "EAAAKEAAAK", "type": "EAAAK", "rigid": True},
}


def design_linker(
    n_term_domain_last_aa: int,
    c_term_domain_first_aa: int,
    linker_type: LinkerType = "GS",
    target_length_aa: int = 5,
) -> str:
    """Return a linker sequence connecting two domains. Strategy A Step 9 helper."""
    if linker_type == "GS":
        repeats = max(1, target_length_aa // 4)
        return "GGGGS" * repeats
    if linker_type == "EAAAK":
        repeats = max(1, target_length_aa // 5)
        return "EAAAK" * repeats
    raise ValueError(f"Unknown linker_type: {linker_type!r}. Use 'GS' or 'EAAAK'.")


def estimate_linker_plddt_expectation(linker_type: LinkerType) -> str:
    """Return expected pLDDT range for linker type (for Step 12 region annotation)."""
    if linker_type == "GS":
        return "30–60 (flexible; low pLDDT expected and acceptable)"
    if linker_type == "EAAAK":
        return "50–75 (rigid helix; moderate pLDDT)"
    return "unknown"
