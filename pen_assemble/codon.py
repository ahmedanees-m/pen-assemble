"""
Codon optimisation utilities for human expression.

Implements rule-based codon optimisation using the Kazusa Homo sapiens
high-expression preferred-codon table. Provides restriction-site scanning
and ORF assembly helpers.

Note
----
This is a rule-based optimiser (one preferred codon per amino acid). For
synthesis orders, verify the codon-adaptation index (CAI) with a commercial
tool (IDT CodonOpt, Twist, GeneArt) before submitting.
"""

from __future__ import annotations

__all__ = [
    "CODON_TABLE_HUMAN",
    "RESTRICTION_SITES",
    "codon_optimise",
    "gc_content",
    "check_restriction_sites",
    "build_expression_orf",
]

# ---------------------------------------------------------------------------
# Kazusa Homo sapiens preferred-codon table
# (most frequent codon per amino acid in human high-expression genes)
# ---------------------------------------------------------------------------
CODON_TABLE_HUMAN: dict[str, str] = {
    "A": "GCC",  # Ala
    "R": "AGG",  # Arg
    "N": "AAC",  # Asn
    "D": "GAC",  # Asp
    "C": "TGC",  # Cys
    "Q": "CAG",  # Gln
    "E": "GAG",  # Glu
    "G": "GGC",  # Gly
    "H": "CAC",  # His
    "I": "ATC",  # Ile
    "L": "CTG",  # Leu
    "K": "AAG",  # Lys
    "M": "ATG",  # Met
    "F": "TTC",  # Phe
    "P": "CCC",  # Pro
    "S": "AGC",  # Ser
    "T": "ACC",  # Thr
    "W": "TGG",  # Trp
    "Y": "TAC",  # Tyr
    "V": "GTG",  # Val
    "*": "TGA",  # Stop (preferred)
}

# ---------------------------------------------------------------------------
# Common restriction enzyme recognition sites
# ---------------------------------------------------------------------------
RESTRICTION_SITES: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "NheI": "GCTAGC",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def codon_optimise(aa_sequence: str) -> str:
    """Translate an amino acid sequence to human-preferred-codon DNA.

    Unknown amino acid symbols are replaced with ``NNN``.

    Parameters
    ----------
    aa_sequence:
        Single-letter amino acid sequence (case-insensitive).
        Stop codons (``*``) are included if present.

    Returns
    -------
    str
        DNA sequence (uppercase, no spaces).

    Examples
    --------
    >>> from pen_assemble.codon import codon_optimise
    >>> codon_optimise("MA")
    'ATGGCC'
    >>> codon_optimise("M*")
    'ATGTGA'
    """
    return "".join(CODON_TABLE_HUMAN.get(aa.upper(), "NNN") for aa in aa_sequence)


def gc_content(dna: str) -> float:
    """Return the GC fraction of a DNA string.

    Parameters
    ----------
    dna:
        DNA sequence (any case).

    Returns
    -------
    float
        GC fraction in [0, 1], or 0.0 for empty strings.

    Examples
    --------
    >>> from pen_assemble.codon import gc_content
    >>> round(gc_content("GCGCAT"), 4)
    0.6667
    """
    if not dna:
        return 0.0
    uc = dna.upper()
    return sum(1 for b in uc if b in "GC") / len(uc)


def check_restriction_sites(dna: str) -> list[str]:
    """Return names of restriction enzymes whose sites appear in *dna*.

    Checks both strands (forward only; palindromic sites are self-complementary
    so forward-strand search is sufficient for all sites in :data:`RESTRICTION_SITES`).

    Parameters
    ----------
    dna:
        DNA sequence (any case).

    Returns
    -------
    list[str]
        Sorted list of enzyme names with sites found in *dna*. Empty if none.

    Examples
    --------
    >>> from pen_assemble.codon import check_restriction_sites
    >>> check_restriction_sites("AAAGAATTCAAA")
    ['EcoRI']
    >>> check_restriction_sites("ACGTACGT")
    []
    """
    uc = dna.upper()
    return sorted(name for name, site in RESTRICTION_SITES.items() if site in uc)


def build_expression_orf(
    aa_sequence: str,
    kozak: bool = True,
    stop: bool = True,
) -> str:
    """Build a full expression-ready ORF from an amino acid sequence.

    Applies codon optimisation then optionally prepends a Kozak consensus
    (``GCCACC``) and appends a preferred stop codon (``TGA``).

    Parameters
    ----------
    aa_sequence:
        Amino acid sequence (case-insensitive). Must begin with ``M`` (Met).
    kozak:
        If ``True`` (default), prepend Kozak context ``GCCACC``.
    stop:
        If ``True`` (default), append stop codon ``TGA``.

    Returns
    -------
    str
        DNA sequence ready for gene synthesis.

    Raises
    ------
    ValueError
        If *aa_sequence* is empty or does not start with Met.

    Examples
    --------
    >>> from pen_assemble.codon import build_expression_orf
    >>> build_expression_orf("MA", kozak=False, stop=False)
    'ATGGCC'
    >>> build_expression_orf("MA", kozak=True, stop=True)
    'GCCACCATGGCCTGA'
    """
    if not aa_sequence:
        raise ValueError("aa_sequence must not be empty")
    if aa_sequence[0].upper() != "M":
        raise ValueError(f"aa_sequence must start with Met (M); got {aa_sequence[0]!r}")
    orf = codon_optimise(aa_sequence)
    if kozak:
        orf = "GCCACC" + orf
    if stop:
        orf = orf + "TGA"
    return orf
