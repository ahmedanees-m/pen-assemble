"""
PEN-ASSEMBLE: Programmable Editor Nomination pipeline for IS110-family bridge
recombinase design and ranking.

Public API
----------
pen_score
    PenScore composite formula (seven mechanistic axes, weighted sum).
catalog
    Catalog loading utilities (load_catalog, load_p1_beaters, load_top5).
codon
    Codon optimisation for human expression (Kazusa table, restriction-site scan).

Quick start
-----------
>>> from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621
>>> ax = PenScoreAxes(S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0,
...                   S_Deliv=1.0, S_Immuno=0.8777, S_Prog=1.0, S_Mature=0.5)
>>> s = pen_score(ax)
>>> beats_is621(s)
True

>>> from pen_assemble.catalog import load_catalog
>>> df = load_catalog()
>>> df.shape[0]
1029
"""
from pen_assemble._version import __version__, __version_info__

__all__ = ["__version__", "__version_info__"]
