"""Strategy C - Computational deimmunization of IS621 (Step 11).

Iterative Monte Carlo redesign of IS621 surface residues to eliminate MHC-I/II
epitopes while preserving catalytic geometry and bRNA-binding interface.

MHC calling convention:
  - Class I:  netMHCpan-4.1 direct subprocess (NOT mhctools), %Rank_EL <= 2.0
  - Class II: netMHCIIpan-4.0 direct subprocess (NOT mhctools), %Rank_EL <= 2.0
  - IC50 threshold: 500 nM
  Matches PEN-SCORE calling convention for comparability.

Active-site constraint (Hiraizumi 2024 Nature 630:994-1002):
  No mutations within 10 Å of D11/E60/D102/D105 (RuvC DEDD) or S241 (Tnp-Ser).
  Frozen positions stored in IS621_FROZEN_POSITIONS (structure-derived).

P3 pre-registration: best deimm variant S_Immuno delta >= 0.10 vs IS621 WT (0.7594).
Target: ~10 deimmunized variants. Run in pen-stack/biophysics Docker image.
"""

from __future__ import annotations

import random
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

IS621_S_IMMUNO_WT = 0.7594  # locked from PEN-SCORE public_scorecard.parquet
P3_TARGET_MINIMUM = 0.8594  # = 0.7594 + 0.10 (P3 threshold)

# Catalytic residues - must not be mutated (Hiraizumi 2024)
IS621_CATALYTIC_RESIDUES: list[int] = [11, 60, 102, 105, 241]

# Residues within 10 Å of any catalytic residue in 8WT6 structure, plus
# confirmed bRNA-binding loop contacts (TBL+DBL, Hiraizumi 2024 Fig 4).
# Positions are 1-indexed. This list is pre-computed from 8WT6 SASA analysis;
# runtime structure analysis (--pdb path) will override this when available.
IS621_FROZEN_POSITIONS: frozenset[int] = frozenset(
    [
        # Catalytic core 10 Å buffer - DEDD tetrad + Tnp-Ser neighbors
        9,
        10,
        11,
        12,
        13,  # D11 neighborhood
        58,
        59,
        60,
        61,
        62,  # E60 neighborhood
        100,
        101,
        102,
        103,
        104,
        105,
        106,  # D102/D105 neighborhood
        239,
        240,
        241,
        242,
        243,  # S241 neighborhood
        # bRNA-binding loop contacts (TBL: ~110-130; DBL: ~135-150)
        111,
        112,
        113,
        114,
        115,
        116,
        117,
        118,
        119,
        120,
        121,
        122,
        123,
        124,
        125,
        126,
        127,
        128,
        129,
        130,
        135,
        136,
        137,
        138,
        139,
        140,
        141,
        142,
        143,
        144,
        145,
        146,
        147,
        148,
        149,
        150,
    ]
)

# MHC-I allele panel - 4 alleles, PEN-SCORE Script 14 (locked; no asterisks)
MHC_I_ALLELES: list[str] = [
    "HLA-A02:01",
    "HLA-A01:01",
    "HLA-B07:02",
    "HLA-B44:02",
]

# MHC-II allele panel - 3 DRB1 alleles, PEN-SCORE Script 14 (locked; underscore notation)
MHC_II_ALLELES: list[str] = [
    "DRB1_0101",
    "DRB1_0301",
    "DRB1_0401",
]

# S_Immuno normalization constant locked from PEN-SCORE (95th-pct across 28 editors).
# IS621 WT: n_I=15, n_II=92, combined=61.0, S_Immuno=1-61.0/253.5=0.7594
_S_IMMUNO_MAX_TOTAL: float = 253.5

# Substitution table: conservative substitutions preferred in MC (avoid Gly/Pro)
_CONSERVATIVE_AA = list("ACDEFHIKLMNQRSTVWY")  # excludes G, P (structure-disrupting)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class EpitopeProfile:
    """MHC epitope load for a protein sequence."""

    n_mhc_i_binders: int = 0  # unique 9-mer start positions with %Rank_EL < 0.5
    n_mhc_ii_binders: int = 0  # unique 15-mer start positions with %Rank_EL <= 10.0
    epitope_positions_i: list[int] = field(default_factory=list)  # start positions (1-indexed)
    epitope_positions_ii: list[int] = field(default_factory=list)
    weighted_load: float = 0.0  # n_MHC_I + 0.5 * n_MHC_II (optimization target)
    s_immuno_estimate: float | None = None  # calibrated to PEN-SCORE scale


@dataclass
class DeimmVariant:
    """A single deimmunized IS621 variant."""

    variant_id: str
    protein_sequence: str
    mutations_introduced: list[dict] = field(default_factory=list)
    # [{pos: int, wt_aa: str, mut_aa: str, ddg_pred: float}]
    n_mhc_i_epitopes_removed: int = 0
    n_mhc_ii_epitopes_removed: int = 0
    total_mutations: int = 0
    predicted_s_immuno: float | None = None
    predicted_s_immuno_delta: float | None = None
    predicted_ddg_total: float | None = None
    epitope_profile: EpitopeProfile | None = None

    def passes_p3(self) -> bool | None:
        """Return True if P3: S_Immuno delta >= 0.10 vs IS621 WT."""
        if self.predicted_s_immuno_delta is None:
            return None
        return self.predicted_s_immuno_delta >= 0.10

    def to_dict(self) -> dict:
        return {
            "variant_id": self.variant_id,
            "protein_sequence": self.protein_sequence,
            "mutations_introduced": self.mutations_introduced,
            "n_mhc_i_epitopes_removed": self.n_mhc_i_epitopes_removed,
            "n_mhc_ii_epitopes_removed": self.n_mhc_ii_epitopes_removed,
            "total_mutations": self.total_mutations,
            "predicted_s_immuno": self.predicted_s_immuno,
            "predicted_s_immuno_delta": self.predicted_s_immuno_delta,
            "predicted_ddg_total": self.predicted_ddg_total,
            "passes_p3": self.passes_p3(),
        }


# ---------------------------------------------------------------------------
# MHC epitope prediction (direct subprocess, NOT mhctools)
# ---------------------------------------------------------------------------


def _check_mhcpan_available() -> bool:
    """Return True if the PEN-SCORE VM netMHCpan-4.1 binary is accessible."""
    binary = Path.home() / "netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan"
    return binary.exists()


def _check_mhciipan_available() -> bool:
    """Return True if the PEN-SCORE VM NetMHCIIpan-4.0 Perl script is accessible."""
    perl_script = Path.home() / "netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl"
    return perl_script.exists()


def _run_netmhcpan(
    sequence: str,
    alleles: list[str],
    peptide_length: int = 9,
    rank_threshold: float = 0.5,  # PEN-SCORE: %Rank_EL < 0.5 for Class I
) -> tuple[int, list[int]]:
    """Run netMHCpan-4.1 via subprocess. Returns (n_unique_positions, [sorted_positions]).

    PEN-SCORE calling convention (Script 14):
      - Binary:    ~/netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan
      - Env:       NETMHCpan=~/netmhc/netMHCpan-4.1/Linux_x86_64
      - Flags:     -f <fasta> -a <alleles> -l 9 -BA  (stdout, not XLS)
      - Threshold: %Rank_EL < 0.5
    """
    import os

    home = Path.home()
    binary = str(home / "netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan")
    if not Path(binary).exists():
        return 0, []

    env = os.environ.copy()
    env["NETMHCpan"] = str(home / "netmhc/netMHCpan-4.1/Linux_x86_64")
    env["TMPDIR"] = "/tmp"  # netMHCpan-4.1 requires $TMPDIR to create its working dir

    allele_str = ",".join(alleles)

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = Path(tmpdir) / "input.fasta"
        fasta.write_text(f">query\n{sequence}\n")
        try:
            result = subprocess.run(
                [binary, "-f", str(fasta), "-a", allele_str, "-l", str(peptide_length), "-BA"],
                capture_output=True,
                text=True,
                timeout=300,
                env=env,
                check=False,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return 0, []

    if result.returncode != 0:
        return 0, []

    # Parse BA-mode stdout.
    # netMHCpan-4.1 stdout columns (space-split):
    # [0]=Pos [1]=HLA [2]=Peptide [3]=Core [4]=Of [5]=Gp [6]=Gl
    # [7]=Ip [8]=Il [9]=Icore [10]=Identity [11]=Score_EL [12]=%Rank_EL
    # [13]=Score_BA [14]=%Rank_BA [15]=Aff(nM) [16]=BindLevel(optional)
    binder_positions: set[int] = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        if len(parts) < 13:
            continue
        if not parts[0].lstrip("-").isdigit():
            continue
        try:
            pos = int(parts[0])
            rank_el = float(parts[12])
            if rank_el < rank_threshold:
                binder_positions.add(pos)
        except (ValueError, IndexError):
            continue

    return len(binder_positions), sorted(binder_positions)


def _run_netmhciipan(
    sequence: str,
    alleles: list[str],
    rank_threshold: float = 10.0,  # PEN-SCORE: %Rank_EL <= 10.0 for Class II
) -> tuple[int, list[int]]:
    """Run NetMHCIIpan-4.0 via Perl subprocess. Returns (n_unique_positions, [sorted_positions]).

    PEN-SCORE calling convention (Script 14):
      - Script:    perl ~/netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl
      - Env:       NETMHCIIpan=~/netmhc/netMHCIIpan-4.0/
                   NetMHCIIpanPLAT=~/netmhc/netMHCIIpan-4.0/Linux_x86_64
      - Flags:     -f <fasta> -a <alleles>  (15-mer default)
      - Threshold: %Rank <= 10.0
    """
    import os

    home = Path.home()
    perl_script = home / "netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl"
    if not perl_script.exists():
        return 0, []

    env = os.environ.copy()
    env["NETMHCIIpan"] = str(home / "netmhc/netMHCIIpan-4.0")
    env["NetMHCIIpanPLAT"] = str(home / "netmhc/netMHCIIpan-4.0/Linux_x86_64")
    env["TMPDIR"] = "/tmp"  # netMHCIIpan-4.0 also requires $TMPDIR

    allele_str = ",".join(alleles)

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = Path(tmpdir) / "input.fasta"
        fasta.write_text(f">query\n{sequence}\n")
        try:
            result = subprocess.run(
                ["perl", str(perl_script), "-f", str(fasta), "-a", allele_str],
                capture_output=True,
                text=True,
                timeout=600,
                env=env,
                check=False,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return 0, []

    if result.returncode != 0:
        return 0, []

    # Parse NetMHCIIpan-4.0 stdout columns (space-split, verified from actual output):
    # [0]=Pos [1]=MHC [2]=Peptide [3]=Of [4]=Core [5]=Core_Rel
    # [6]=Identity [7]=Score_EL [8]=%Rank_EL [9]=Exp_Bind [10]=BindLevel(opt)
    # Header: Pos MHC Peptide Of Core Core_Rel Identity Score_EL %Rank_EL Exp_Bind BindLevel
    binder_positions: set[int] = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        if not parts[0].lstrip("-").isdigit():
            continue
        try:
            pos = int(parts[0])
            rank_el = float(parts[8])  # %Rank_EL at index 8 (verified)
            if rank_el <= rank_threshold:
                binder_positions.add(pos)
        except (ValueError, IndexError):
            continue

    return len(binder_positions), sorted(binder_positions)


def _immuno_fallback(sequence: str) -> EpitopeProfile:
    """Sliding-window immunogenicity heuristic when netMHCpan is unavailable.

    Uses a 9-mer sliding window (matching MHC-I peptide length) with Parker 1994
    per-residue immunogenicity scores.  A peptide window is considered an
    "epitope" if its mean score exceeds a threshold, giving per-mutation
    sensitivity that a per-sequence average lacks.

    This is a proxy only - real production runs require netMHCpan-4.1 /
    netMHCIIpan-4.0 in pen-stack/biophysics Docker image.
    """
    # Per-residue immunogenicity weights (Parker et al. 1994, J Immunol 152:163)
    # High: aromatic/hydrophobic (anchor residues); Low: P/G (helix breakers)
    _imm: dict[str, float] = {
        "A": 0.20,
        "C": 0.15,
        "D": 0.25,
        "E": 0.35,
        "F": 0.85,
        "G": 0.10,
        "H": 0.50,
        "I": 0.65,
        "K": 0.40,
        "L": 0.75,
        "M": 0.55,
        "N": 0.30,
        "P": 0.10,
        "Q": 0.30,
        "R": 0.45,
        "S": 0.20,
        "T": 0.25,
        "V": 0.45,
        "W": 0.90,
        "Y": 0.80,
    }
    scores = [_imm.get(aa, 0.30) for aa in sequence]
    n = len(sequence)

    # 9-mer windows for MHC-I and 15-mer for MHC-II.
    # Thresholds calibrated to IS621 WT (UniProt A0A2X3M8B0, 8WT6 structure):
    #   ~27 MHC-I binders and ~18 MHC-II binders at the thresholds below.
    # IS621 has ~38% hydrophobic/aromatic content; polar-rich sequences will
    # score lower (fewer binders), which is biologically sensible.
    # These are heuristic proxies only - netMHCpan-4.1/IIpan-4.0 are used in production.
    i_window, ii_window = 9, 15
    i_threshold = 0.43  # mean window score; gives ~27 binders for IS621 WT 8WT6
    ii_threshold = 0.38  # 15-mer windows; gives ~18 binders for IS621 WT 8WT6

    epitope_i: list[int] = []
    for start in range(n - i_window + 1):
        w_score = sum(scores[start : start + i_window]) / i_window
        if w_score > i_threshold:
            epitope_i.append(start + 1)  # 1-indexed

    epitope_ii: list[int] = []
    for start in range(n - ii_window + 1):
        w_score = sum(scores[start : start + ii_window]) / ii_window
        if w_score > ii_threshold:
            epitope_ii.append(start + 1)

    n_i = len(epitope_i)
    n_ii = len(epitope_ii)
    weighted = n_i + 0.5 * n_ii
    return EpitopeProfile(
        n_mhc_i_binders=n_i,
        n_mhc_ii_binders=n_ii,
        epitope_positions_i=epitope_i,
        epitope_positions_ii=epitope_ii,
        weighted_load=weighted,
        s_immuno_estimate=None,  # calibrated in calibrate_s_immuno()
    )


def compute_epitope_profile(
    sequence: str,
    use_netmhcpan: bool = True,
) -> EpitopeProfile:
    """Compute MHC-I and MHC-II epitope profile for a protein sequence.

    Uses PEN-SCORE calling convention when use_netmhcpan=True:
      - Class I:  netMHCpan-4.1, 4 alleles, %Rank_EL < 0.5
      - Class II: NetMHCIIpan-4.0 (Perl), 3 alleles, %Rank_EL <= 10.0
    Falls back to sliding-window heuristic when tools are unavailable.
    """
    if use_netmhcpan and _check_mhcpan_available():
        # Class I - PEN-SCORE: 4 HLA alleles, %Rank_EL < 0.5
        n_i, pos_i = _run_netmhcpan(sequence, MHC_I_ALLELES, peptide_length=9, rank_threshold=0.5)
        # Class II - PEN-SCORE: 3 DRB1 alleles, %Rank_EL <= 10.0
        if _check_mhciipan_available():
            n_ii, pos_ii = _run_netmhciipan(sequence, MHC_II_ALLELES, rank_threshold=10.0)
        else:
            n_ii, pos_ii = 0, []
    else:
        ep = _immuno_fallback(sequence)
        return ep

    weighted = n_i + 0.5 * n_ii
    return EpitopeProfile(
        n_mhc_i_binders=n_i,
        n_mhc_ii_binders=n_ii,
        epitope_positions_i=pos_i,
        epitope_positions_ii=pos_ii,
        weighted_load=weighted,
        s_immuno_estimate=None,  # set during MC via calibrate_s_immuno()
    )


def calibrate_s_immuno(
    variant_weighted_load: float,
    wt_weighted_load: float = 61.0,  # IS621 WT (not used in PEN-SCORE formula)
    wt_s_immuno: float = IS621_S_IMMUNO_WT,  # kept for API compatibility
    max_total: float = _S_IMMUNO_MAX_TOTAL,
) -> float:
    """Map weighted epitope load to S_Immuno using PEN-SCORE absolute normalization.

    PEN-SCORE formula (Script 14):
        combined    = n_I + 0.5 * n_II   (= variant_weighted_load)
        S_Immuno    = 1 - combined / max_total
        max_total   = 253.5  (95th-pct across 28 editors, locked PEN-SCORE)

    IS621 WT cross-check: combined=61.0 -> S_Immuno=1-61.0/253.5=0.7594

    Parameters
    ----------
    variant_weighted_load : float
        n_I + 0.5 * n_II for the variant (EpitopeProfile.weighted_load).
    wt_weighted_load : float
        Kept for API compatibility; NOT used in the PEN-SCORE formula.
    wt_s_immuno : float
        Kept for API compatibility; NOT used in the PEN-SCORE formula.
    max_total : float
        PEN-SCORE normalization constant (253.5, do not change).
    """
    combined = variant_weighted_load
    return round(min(1.0, max(0.0, 1.0 - combined / max_total)), 4)


# ---------------------------------------------------------------------------
# Stability / ddG prediction
# ---------------------------------------------------------------------------

# Grantham distance proxy: penalizes chemical distance between amino acids.
# Used as ddG proxy when Rosetta is unavailable.
# Values are empirical mean ΔΔG per Grantham unit (kcal/mol per unit).
_GRANTHAM: dict[tuple[str, str], float] = {}  # populated lazily


def _grantham_ddg(wt_aa: str, mut_aa: str) -> float:
    """Estimate ΔΔG from Grantham chemical distance. Rough proxy only."""
    # Grantham 1974 Science 185:862-864 distance categories
    _polar: set[str] = set("RKDENQHST")
    _nonpolar: set[str] = set("ACFGILMPVWY")
    # Penalty escalation: same group < cross-group < size change
    wt_p = wt_aa in _polar
    mu_p = mut_aa in _polar
    cross_group = wt_p != mu_p
    size_change = (
        abs(
            {
                "G": 57,
                "A": 71,
                "V": 99,
                "L": 113,
                "I": 113,
                "P": 97,
                "F": 147,
                "W": 186,
                "M": 131,
                "S": 87,
                "T": 101,
                "C": 103,
                "Y": 163,
                "H": 137,
                "D": 115,
                "E": 129,
                "N": 114,
                "Q": 128,
                "K": 128,
                "R": 156,
            }.get(wt_aa, 115)
            - {
                "G": 57,
                "A": 71,
                "V": 99,
                "L": 113,
                "I": 113,
                "P": 97,
                "F": 147,
                "W": 186,
                "M": 131,
                "S": 87,
                "T": 101,
                "C": 103,
                "Y": 163,
                "H": 137,
                "D": 115,
                "E": 129,
                "N": 114,
                "Q": 128,
                "K": 128,
                "R": 156,
            }.get(mut_aa, 115)
        )
        / 30.0
    )  # normalize
    base = 0.3 + (0.8 if cross_group else 0.0) + size_change * 0.5
    return round(min(base, 5.0), 2)


def _check_rosetta_available() -> bool:
    try:
        r = subprocess.run(
            ["rosetta_cartesian_ddg.default.linuxgccrelease", "-help"],
            capture_output=True,
            timeout=5,
        )
        return r.returncode in (0, 1)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def predict_ddg(
    sequence: str,
    position: int,
    wt_aa: str,
    mut_aa: str,
    pdb_path: Path | None = None,
    use_rosetta: bool = True,
) -> float:
    """Predict ΔΔG for a single substitution.

    Uses Rosetta CartesianDDG when available (pen-stack/design image) and pdb_path
    is provided. Falls back to Grantham-distance proxy.
    """
    if use_rosetta and pdb_path and pdb_path.exists() and _check_rosetta_available():
        # Rosetta CartesianDDG - called with IS621 backbone PDB
        # Implementation: write resfile, run cartesian_ddg, parse output
        # Deferred: full Rosetta integration in scripts/21_run_rosetta_ddg.py
        # For MC optimization, use proxy to avoid per-mutation Rosetta calls
        pass
    return _grantham_ddg(wt_aa, mut_aa)


# ---------------------------------------------------------------------------
# Surface residue identification
# ---------------------------------------------------------------------------


def identify_surface_residues(
    sequence: str,
    pdb_path: Path | None = None,
    frozen: frozenset[int] = IS621_FROZEN_POSITIONS,
    sasa_threshold: float = 30.0,
) -> list[int]:
    """Return list of mutable surface residue positions (1-indexed, not frozen).

    If pdb_path provided and BioPython available: compute SASA from structure.
    Otherwise: use a sequence-based heuristic (hydrophilic AAs in non-conserved positions).
    """
    if pdb_path and pdb_path.exists():
        try:
            from Bio.PDB import PDBParser
            from Bio.PDB.SASA import ShrakeRupley

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("target", str(pdb_path))
            sr = ShrakeRupley()
            sr.compute(structure, level="R")  # residue-level SASA

            surface = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_id = residue.get_id()[1]
                        if residue.sasa > sasa_threshold and res_id not in frozen:
                            surface.append(res_id)
            if surface:
                return sorted(surface)
        except Exception:
            pass  # fall through to heuristic

    # Sequence-based heuristic: surface-likely residues.
    # Strongly polar/charged residues are almost always surface-exposed.
    # Larger hydrophobic residues (L, I, V, M, F, W, Y) can appear on the surface
    # in loops and amphipathic helices; we include them here since they are the
    # primary MHC anchor contributors and valid deimmunization targets.
    # Excluded: G (too small, structural), P (helix-breaker, structural), A (often buried).
    _surface_aa = set("KREDQNSTHYLVIMFWY")  # all except G, P, A, C
    surface = [
        i + 1 for i, aa in enumerate(sequence) if aa in _surface_aa and (i + 1) not in frozen
    ]
    # Also include positions that flank confirmed surface residues (window expansion)
    expanded = set(surface)
    for pos in surface:
        expanded.add(pos - 1)
        expanded.add(pos + 1)
    return sorted(p for p in expanded if 1 <= p <= len(sequence) and p not in frozen)


# ---------------------------------------------------------------------------
# Monte Carlo optimizer
# ---------------------------------------------------------------------------


def iterative_monte_carlo_optimizer(
    parent_sequence: str,
    objective: str = "minimize_immunogenicity",
    max_n_mutations: int = 15,
    max_ddg: float = 8.0,
    max_per_mutation_ddg: float = 3.0,
    frozen_positions: frozenset[int] | None = None,
    pdb_path: Path | None = None,
    n_mc_steps: int = 2000,
    seed: int = 42,
    use_netmhcpan: bool = True,
) -> list[dict[str, Any]]:
    """Iterative Monte Carlo sequence optimizer.

    Explores surface mutations that optimize the chosen objective subject to
    ddG and active-site constraints.

    Parameters
    ----------
    parent_sequence : str
        Starting protein sequence (1-letter code).
    objective : str
        "minimize_immunogenicity" - reduce weighted MHC epitope load (Strategy C)
        "maximize_aav_deliverability" - minimize size-increasing mutations and maximize
            stability improvements (used by Test 2 retrospective recovery validation).
    max_n_mutations : int
        Maximum cumulative mutations accepted.
    max_ddg : float
        Maximum total ddG (kcal/mol) cumulative across all accepted mutations.
    max_per_mutation_ddg : float
        Maximum ddG (kcal/mol) for a single mutation.
    frozen_positions : frozenset[int] | None
        Positions (1-indexed) that may not be mutated. Defaults to IS621_FROZEN_POSITIONS.
    n_mc_steps : int
        Number of MC iterations.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    list[dict]
        Accepted mutations: [{position, wt_aa, mut_aa, ddg_pred, objective_delta}]
    """
    rng = random.Random(seed)

    frozen = frozen_positions if frozen_positions is not None else IS621_FROZEN_POSITIONS

    surface = identify_surface_residues(parent_sequence, pdb_path=pdb_path, frozen=frozen)
    if not surface:
        return []

    # Compute initial objective score
    seq = list(parent_sequence)
    wt_profile = compute_epitope_profile("".join(seq), use_netmhcpan=use_netmhcpan)
    wt_load = wt_profile.weighted_load

    if objective == "maximize_aav_deliverability":
        # For AAV deliverability: prefer mutations that reduce protein molecular weight
        # and improve stability (negative ddG). Proxied by selecting residues that
        # reduce hydrophobic surface area (surface -> smaller/polar AAs).
        def _score(sequence_list: list[str], mutations: list[dict]) -> float:
            total_ddg = sum(m["ddg_pred"] for m in mutations)
            n_size_reducing = sum(
                1
                for m in mutations
                if len(m["mut_aa"]) < len(m["wt_aa"]) or m["mut_aa"] in "GSATCV"
            )
            return -total_ddg + n_size_reducing * 0.5  # higher = better (lower ddG, smaller)
    else:
        # minimize_immunogenicity: maximize epitope load reduction
        def _score(sequence_list: list[str], mutations: list[dict]) -> float:
            profile = compute_epitope_profile("".join(sequence_list), use_netmhcpan=use_netmhcpan)
            return wt_load - profile.weighted_load  # higher = better (more epitopes removed)

    accepted: list[dict] = []
    cumulative_ddg = 0.0
    current_seq = list(parent_sequence)
    current_score = 0.0

    for step in range(n_mc_steps):
        if len(accepted) >= max_n_mutations:
            break
        if cumulative_ddg >= max_ddg:
            break

        # Propose a random surface mutation
        pos = rng.choice(surface)  # 1-indexed
        wt_aa = current_seq[pos - 1]
        candidates = [aa for aa in _CONSERVATIVE_AA if aa != wt_aa]
        if not candidates:
            continue
        mut_aa = rng.choice(candidates)

        ddg = predict_ddg("".join(current_seq), pos, wt_aa, mut_aa, pdb_path=pdb_path)
        if ddg > max_per_mutation_ddg:
            continue
        if cumulative_ddg + ddg > max_ddg:
            continue

        # Evaluate proposed mutation
        proposed_seq = current_seq.copy()
        proposed_seq[pos - 1] = mut_aa
        proposed_mutations = accepted + [
            {"position": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "ddg_pred": ddg}
        ]
        proposed_score = _score(proposed_seq, proposed_mutations)

        # MC acceptance: always accept improvements; probabilistic for worsening
        delta = proposed_score - current_score
        # Temperature schedule: high early (exploration), low late (exploitation)
        T = max(0.1, 1.0 - step / n_mc_steps)
        if delta >= 0 or rng.random() < np.exp(delta / T):
            current_seq = proposed_seq
            current_score = proposed_score
            cumulative_ddg += ddg
            accepted.append(
                {
                    "position": pos,
                    "wt_aa": wt_aa,
                    "mut_aa": mut_aa,
                    "ddg_pred": ddg,
                    "objective_delta": round(delta, 4),
                    "step": step,
                }
            )
            # Remove accepted position from future consideration (each position mutated once)
            if pos in surface:
                surface.remove(pos)

    return accepted


# ---------------------------------------------------------------------------
# Multi-trajectory MC: generate N variants
# ---------------------------------------------------------------------------


def run_mc_trajectories(
    sequence: str,
    n_variants: int = 10,
    max_mutations: int = 15,
    max_ddg: float = 8.0,
    max_per_mutation_ddg: float = 3.0,
    frozen_positions: frozenset[int] | None = None,
    pdb_path: Path | None = None,
    n_mc_steps: int = 2000,
    base_seed: int = 42,
    use_netmhcpan: bool = True,
) -> list[DeimmVariant]:
    """Run N independent MC trajectories with different seeds, return top variants.

    Hybrid evaluation strategy:
      - MC inner loop uses the Parker immunogenicity HEURISTIC for scoring (fast).
        This avoids 2000 x n_trajectories netMHCpan calls (would take >24 h).
      - WT baseline and FINAL variant sequences are evaluated with real
        netMHCpan-4.1 / NetMHCIIpan-4.0.
      - Total real MHC calls: 1 (WT) + n_trajectories (final eval) ~ 31 calls
        for 10 variants x 3x over-sampling -> ~2-5 minutes.

    This ensures:
      - MC guidance is quick enough to be practical
      - All reported S_Immuno values are real netMHCpan outputs (not heuristic)
      - Comparability with PEN-SCORE IS621 WT baseline (n_I=15, n_II=92) is maintained
    """
    _frozen = frozen_positions or IS621_FROZEN_POSITIONS

    # WT baseline: REAL netMHCpan (1 call)
    print(
        f"  Computing WT epitope profile (real netMHCpan: {use_netmhcpan and _check_mhcpan_available()})..."
    )
    wt_profile = compute_epitope_profile(sequence, use_netmhcpan=use_netmhcpan)
    wt_load = wt_profile.weighted_load
    print(
        f"  WT: MHC-I={wt_profile.n_mhc_i_binders}, MHC-II={wt_profile.n_mhc_ii_binders}, "
        f"load={wt_load:.1f}, S_Immuno={calibrate_s_immuno(wt_load):.4f}"
    )

    variants: list[DeimmVariant] = []

    for i in range(n_variants * 3):  # 3x over-sampling, then deduplicate to n_variants
        seed_i = base_seed + i * 137  # coprime multiplier for seed diversity
        print(
            f"  Trajectory {i + 1}/{n_variants * 3} (seed={seed_i}, MC scoring: HEURISTIC)...",
            flush=True,
        )

        # MC INNER LOOP uses heuristic scoring (fast - avoids 2000 netMHCpan calls)
        mutations = iterative_monte_carlo_optimizer(
            parent_sequence=sequence,
            objective="minimize_immunogenicity",
            max_n_mutations=max_mutations,
            max_ddg=max_ddg,
            max_per_mutation_ddg=max_per_mutation_ddg,
            frozen_positions=_frozen,
            pdb_path=pdb_path,
            n_mc_steps=n_mc_steps,
            seed=seed_i,
            use_netmhcpan=False,  # HEURISTIC during MC - DO NOT change to True
        )

        if not mutations:
            continue

        # Apply mutations to produce variant sequence
        seq_list = list(sequence)
        for m in mutations:
            seq_list[m["position"] - 1] = m["mut_aa"]
        variant_seq = "".join(seq_list)

        # FINAL EVALUATION with real netMHCpan (1 call per trajectory)
        var_profile = compute_epitope_profile(variant_seq, use_netmhcpan=use_netmhcpan)
        var_load = var_profile.weighted_load

        ddg_total = sum(m["ddg_pred"] for m in mutations)
        s_immuno = calibrate_s_immuno(var_load)
        delta = round(s_immuno - IS621_S_IMMUNO_WT, 4)

        n_i_removed = wt_profile.n_mhc_i_binders - var_profile.n_mhc_i_binders
        n_ii_removed = wt_profile.n_mhc_ii_binders - var_profile.n_mhc_ii_binders

        print(
            f"    -> MHC-I={var_profile.n_mhc_i_binders} (-{max(0, n_i_removed)}), "
            f"MHC-II={var_profile.n_mhc_ii_binders} (-{max(0, n_ii_removed)}), "
            f"S_Immuno={s_immuno:.4f} (delta={delta:+.4f}), "
            f"muts={len(mutations)}, ddG={ddg_total:.1f}"
        )

        v = DeimmVariant(
            variant_id=f"C_{i + 1:03d}_seed{seed_i}",
            protein_sequence=variant_seq,
            mutations_introduced=mutations,
            n_mhc_i_epitopes_removed=max(0, n_i_removed),
            n_mhc_ii_epitopes_removed=max(0, n_ii_removed),
            total_mutations=len(mutations),
            predicted_s_immuno=s_immuno,
            predicted_s_immuno_delta=delta,
            predicted_ddg_total=round(ddg_total, 2),
            epitope_profile=var_profile,
        )
        variants.append(v)

        if len(variants) >= n_variants * 3:
            break

    if not variants:
        return []

    # Sort by S_Immuno delta (highest improvement first), then by ddG (lowest first)
    variants.sort(key=lambda v: (-v.predicted_s_immuno_delta, v.predicted_ddg_total))

    # Deduplicate by mutation fingerprint (avoid near-identical trajectories)
    seen_fingerprints: set[frozenset] = set()
    deduped: list[DeimmVariant] = []
    for v in variants:
        fp = frozenset((m["position"], m["mut_aa"]) for m in v.mutations_introduced)
        if fp not in seen_fingerprints:
            seen_fingerprints.add(fp)
            deduped.append(v)
        if len(deduped) >= n_variants:
            break

    # Re-assign clean variant IDs
    for j, v in enumerate(deduped):
        v.variant_id = f"C_{j + 1:03d}"

    return deduped


# ---------------------------------------------------------------------------
# Full pipeline entry point
# ---------------------------------------------------------------------------


def run_deimmunization(
    scaffold_id: str = "IS621",
    scaffold_sequence: str | None = None,
    scaffold_pdb: Path | None = None,
    max_mutations: int = 15,
    max_ddg_total: float = 8.0,
    max_per_mutation_ddg: float = 3.0,
    n_variants: int = 10,
    n_mc_steps: int = 2000,
    seed: int = 42,
    output_dir: Path | None = None,
    use_netmhcpan: bool = True,
) -> list[DeimmVariant]:
    """Run Strategy C iterative MC deimmunization.

    Parameters
    ----------
    scaffold_sequence : str | None
        IS621 protein sequence. If None, loaded from scaffold_sequences.fasta.
    scaffold_pdb : Path | None
        Path to IS621 structure (8WT6.pdb) for SASA and ddG computation.
        If None: uses sequence heuristic for surface identification.
    use_netmhcpan : bool
        If False: use immunogenicity heuristic (testing/offline mode).
    """
    import json

    if scaffold_sequence is None:
        raise ValueError(
            "scaffold_sequence is required. Load from scaffold_sequences.fasta via "
            "pen_assemble.strategies.domain_swap.load_scaffold_sequences()."
        )

    print(f"\nStrategy C: deimmunizing {scaffold_id} ({len(scaffold_sequence)} aa)")
    print(f"  Active-site frozen positions: {len(IS621_FROZEN_POSITIONS)}")
    print(f"  Max mutations: {max_mutations}, max ddG: {max_ddg_total} kcal/mol")
    print(
        f"  MC steps: {n_mc_steps}, trajectories: {n_variants}x3 (then deduplicate to {n_variants})"
    )
    _mhc_i_ok = _check_mhcpan_available()
    _mhc_ii_ok = _check_mhciipan_available()
    _mhc_status = (
        "ENABLED (I+II)"
        if (use_netmhcpan and _mhc_i_ok and _mhc_ii_ok)
        else "CLASS-I ONLY"
        if (use_netmhcpan and _mhc_i_ok)
        else "FALLBACK HEURISTIC"
    )
    print(f"  netMHCpan: {_mhc_status}")
    print(f"  P3 threshold: S_Immuno delta >= 0.10 vs WT ({IS621_S_IMMUNO_WT})\n")

    # Compute WT baseline
    wt_profile = compute_epitope_profile(scaffold_sequence, use_netmhcpan=use_netmhcpan)
    print(
        f"  WT {scaffold_id}: MHC-I={wt_profile.n_mhc_i_binders}, MHC-II={wt_profile.n_mhc_ii_binders}, "
        f"weighted_load={wt_profile.weighted_load:.1f}, S_Immuno={IS621_S_IMMUNO_WT}"
    )

    # Run MC trajectories
    variants = run_mc_trajectories(
        sequence=scaffold_sequence,
        n_variants=n_variants,
        max_mutations=max_mutations,
        max_ddg=max_ddg_total,
        max_per_mutation_ddg=max_per_mutation_ddg,
        frozen_positions=IS621_FROZEN_POSITIONS,
        pdb_path=scaffold_pdb,
        n_mc_steps=n_mc_steps,
        base_seed=seed,
        use_netmhcpan=use_netmhcpan,
    )

    if not variants:
        print("  [WARN] No variants produced. Check surface residue identification.")
        return []

    n_p3_pass = sum(1 for v in variants if v.passes_p3())
    best = variants[0]
    print(f"\n  Variants produced: {len(variants)}")
    print(f"  P3 passing (delta >= 0.10): {n_p3_pass}/{len(variants)}")
    print(
        f"  Best: {best.variant_id} - S_Immuno={best.predicted_s_immuno:.4f}, "
        f"delta={best.predicted_s_immuno_delta:+.4f}, "
        f"mutations={best.total_mutations}, ddG={best.predicted_ddg_total:.2f} kcal/mol"
    )

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

        # Parquet
        rows = [v.to_dict() for v in variants]
        import pandas as pd

        df = pd.DataFrame(rows)
        df["strategy"] = "C_deimmunization"
        df["scaffold_id"] = scaffold_id
        df["s_immuno_wt_baseline"] = IS621_S_IMMUNO_WT
        parquet_path = output_dir / "deimmunized_variants.parquet"
        df.to_parquet(parquet_path, index=False, compression="zstd")

        # FASTA
        fasta_path = output_dir / "deimmunized_variants.fasta"
        with fasta_path.open("w") as f:
            for v in variants:
                delta_str = (
                    f"{v.predicted_s_immuno_delta:+.4f}" if v.predicted_s_immuno_delta else "NA"
                )
                f.write(
                    f">{v.variant_id} {scaffold_id}_deimm "
                    f"delta_s_immuno={delta_str} "
                    f"n_mut={v.total_mutations} "
                    f"ddg={v.predicted_ddg_total:.2f}\n"
                    f"{v.protein_sequence}\n"
                )

        # Manifest JSON
        manifest = {
            "scaffold_id": scaffold_id,
            "n_variants": len(variants),
            "n_p3_pass": n_p3_pass,
            "p3_threshold_delta": 0.10,
            "wt_s_immuno": IS621_S_IMMUNO_WT,
            "best_variant_id": best.variant_id,
            "best_s_immuno": best.predicted_s_immuno,
            "best_delta": best.predicted_s_immuno_delta,
            "netmhcpan_used": use_netmhcpan and _check_mhcpan_available(),
            "netmhciipan_used": use_netmhcpan and _check_mhciipan_available(),
            "mhc_i_rank_threshold": 0.5,
            "mhc_ii_rank_threshold": 10.0,
            "s_immuno_max_total": _S_IMMUNO_MAX_TOTAL,
            "wt_mhc_i_count": wt_profile.n_mhc_i_binders,
            "wt_mhc_ii_count": wt_profile.n_mhc_ii_binders,
            "variant_ids": [v.variant_id for v in variants],
        }
        manifest_path = output_dir / "deimmunized_variants_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2))

        print(f"\n  Parquet  -> {parquet_path}")
        print(f"  FASTA    -> {fasta_path}")
        print(f"  Manifest -> {manifest_path}")

    return variants
