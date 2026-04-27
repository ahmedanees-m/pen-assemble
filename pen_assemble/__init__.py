"""
PEN-ASSEMBLE v0.5.0: Programmable Editor Nomination pipeline for IS110-family
bridge recombinase design and ranking.

Subpackages
-----------
pen_score
    PenScore composite formula (7 axes, weighted sum), IS621 lockpoints.
catalog
    Catalog loading utilities (load_catalog, load_p1_beaters, load_top5).
codon
    Human codon optimisation (Kazusa table, restriction-site scanner).
api
    High-level Designer API (compose_chimera, discover_orthologs, deimmunize,
    redesign_backbone, get_catalog, select_designs).
strategies
    Strategy A-D design generators (domain_swap, ortholog_discovery,
    deimmunization, backbone_redesign).
triage
    Multi-gate triage (multi_gate) and diversity analysis (diversity).
verification
    Axis evaluation: mechclass, penscore, stability, structure, active_site.
utils
    Low-level helpers: linker assembly, MHC epitope scoring, PDB parsing.
data
    YAML configuration data (scaffold_universe, design_strategies,
    pre_registration).

Quick start
-----------
>>> from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621
>>> ax = PenScoreAxes(S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0,
...                   S_Deliv=1.0, S_Immuno=0.8777, S_Prog=1.0, S_Mature=1.0)
>>> s = pen_score(ax)
>>> beats_is621(s)
True

>>> from pen_assemble.catalog import load_catalog
>>> df = load_catalog()
>>> df.shape[0]
1029

>>> from pen_assemble.api import Designer
>>> d = Designer.load()
"""

from pen_assemble._version import __version__, __version_info__

__all__ = ["__version__", "__version_info__"]
