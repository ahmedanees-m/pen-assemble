"""PEN-ASSEMBLE: Computational design framework for programmable genome-writing editors.

Part of the PEN-STACK trilogy (Paper 4). Builds on GENOME-ATLAS (Paper 1),
MECH-CLASS (Paper 2), and PEN-SCORE (Paper 3).
"""
from __future__ import annotations

try:
    from pen_assemble._version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ["__version__"]
