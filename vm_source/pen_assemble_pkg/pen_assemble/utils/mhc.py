"""MHC epitope prediction helpers for S_Immuno axis and Strategy C deimmunization.

CALLING CONVENTION (Item 4 — must be consistent with Paper 3 S_Immuno pipeline):
=================================================================================
Paper 3 (pen-score v0.1.0) computed S_Immuno via DIRECT subprocess calls to the
netMHCpan-4.1 and netMHCIIpan-4.0 binaries, NOT via the mhctools wrapper. The raw
output is parsed from netMHCpan's tab-separated stdout with column '%Rank_EL'.

Paper 4 MUST use the same calling convention and output parser to ensure epitope
counts are comparable across papers (the IS621 S_Immuno=0.7594 baseline was computed
with Paper 3's direct-subprocess approach).

The [immuno] optional extra installs mhctools>=1.9 as a Python wrapper around the
same netMHCpan binary, BUT mhctools uses a different calling convention:
  - Binary name: 'NetMHCpan' (uppercase N) vs Paper 3's 'netMHCpan-4.1'
  - Output format: mhctools returns BindingPrediction objects; Paper 3 parsed
    raw TSV with '%Rank_EL' column and IC50_nM threshold of 500 nM
  - IC50 thresholding: mhctools defaults differ from Paper 3's 500 nM strong-binder cutoff

DECISION: Paper 4 uses the same direct subprocess approach as Paper 3 (see
pen_score/axes/s_immuno.py for the reference implementation). The mhctools import
is kept for optional convenience but is NOT used in the production S_Immuno pipeline.
This is documented here to prevent a future maintainer from silently switching to
mhctools and breaking comparability with Paper 3 scores.

See also: Paper 3 execution summary Deviation D4 (MHCflurry 2.0 rejected as
class-I only; netMHCpan-4.1 used for class I + netMHCIIpan-4.0 for class II).
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

# Allele sets used in Paper 3 S_Immuno (locked — must match for comparability).
# Source: PAPER_3/pen-score/scripts/14_compute_S_Immuno.py
# Class I: 4 alleles, NO asterisks (netMHCpan-4.1 format)
# Class II: 3 DRB1 alleles, underscore notation (netMHCIIpan-4.0 Perl format)
MHC_I_ALLELES_PAPER3 = [
    "HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B44:02",
]
MHC_II_ALLELES_PAPER3 = [
    "DRB1_0101", "DRB1_0301", "DRB1_0401",
]

# %Rank_EL thresholds — Paper 3 calls (NOT IC50 nM).
RANK_EL_THRESHOLD_I = 0.5    # Class I: %Rank_EL < 0.5 → strong binder (Paper 3 Script 14)
RANK_EL_THRESHOLD_II = 10.0  # Class II: %Rank_EL ≤ 10.0 (Paper 3 Script 14)

# Paper 3 S_Immuno normalization constant (95th-pct across 28 editors, locked).
# IS621 WT: n_I=15, n_II=92 → combined=61.0 → S_Immuno=1-61.0/253.5=0.7594
S_IMMUNO_MAX_TOTAL = 253.5
IS621_S_IMMUNO_WT = 0.7594


def find_mhc_i_binders(
    sequence: str,
    alleles: Optional[list[str]] = None,
    peptide_length: int = 9,
    rank_el_threshold: float = RANK_EL_THRESHOLD_I,
    netmhcpan_binary: Optional[str] = None,
) -> list[dict]:
    """Find MHC-I strong binders via direct netMHCpan-4.1 subprocess call.

    Paper 3 calling convention (Script 14):
      - Binary: ~/netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan
      - Env: NETMHCpan=~/netmhc/netMHCpan-4.1/Linux_x86_64
      - Flags: -f <fasta> -a <alleles> -l 9 -BA
      - Threshold: %Rank_EL < 0.5 (NOT IC50 nM)
      - Returns: unique start positions of binding 9-mers

    Returns list of {peptide, start, allele, rank_el} — one entry per
    unique (position, allele) pair that passes the threshold.

    Raises RuntimeError if netMHCpan binary is not found.
    """
    import os
    import shutil
    import tempfile

    _alleles = alleles or MHC_I_ALLELES_PAPER3

    # Resolve binary: explicit path > PATH lookup > Paper 3 VM default
    home = Path.home()
    _vm_binary = str(home / "netmhc/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan")
    _binary = netmhcpan_binary or _vm_binary
    if not Path(_binary).exists() and not shutil.which(_binary):
        raise RuntimeError(
            f"netMHCpan binary not found at {_binary!r}. "
            "Run inside pen-stack VM or biophysics Docker image."
        )

    # Set required env vars (netMHCpan-4.1 needs NETMHCpan + TMPDIR)
    env = os.environ.copy()
    env["NETMHCpan"] = str(home / "netmhc/netMHCpan-4.1/Linux_x86_64")
    env["TMPDIR"] = "/tmp"   # netMHCpan-4.1 uses mkdtemp under $TMPDIR; must be set

    allele_str = ",".join(_alleles)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as fh:
        fh.write(f">query\n{sequence}\n")
        fasta_path = fh.name

    try:
        result = subprocess.run(
            [_binary, "-f", fasta_path, "-a", allele_str,
             "-l", str(peptide_length), "-BA"],
            capture_output=True, text=True, timeout=300, env=env,
        )
    finally:
        Path(fasta_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"netMHCpan failed (rc={result.returncode}): {result.stderr[:400]}"
        )

    binders: list[dict] = []
    seen: set[tuple] = set()  # deduplicate (pos, allele) pairs
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        # netMHCpan-4.1 BA-mode stdout columns (space-split):
        # [0]=Pos [1]=HLA [2]=Peptide [3]=Core [4]=Of [5]=Gp [6]=Gl
        # [7]=Ip [8]=Il [9]=Icore [10]=Identity [11]=Score_EL [12]=%Rank_EL
        # [13]=Score_BA [14]=%Rank_BA [15]=Aff(nM) [16]=BindLevel(optional)
        if len(parts) < 13:
            continue
        if not parts[0].lstrip("-").isdigit():
            continue  # skip header/separator rows
        try:
            pos = int(parts[0])
            allele = parts[1]
            peptide = parts[2]
            rank_el = float(parts[12])
        except (ValueError, IndexError):
            continue
        if rank_el < rank_el_threshold:
            key = (pos, allele)
            if key not in seen:
                seen.add(key)
                binders.append({
                    "peptide": peptide,
                    "start": pos,
                    "allele": allele,
                    "rank_el": rank_el,
                })

    return binders


def find_mhc_ii_binders(
    sequence: str,
    alleles: Optional[list[str]] = None,
    rank_el_threshold: float = RANK_EL_THRESHOLD_II,
    netmhciipan_perl: Optional[str] = None,
) -> list[dict]:
    """Find MHC-II strong binders via netMHCIIpan-4.0 Perl script call.

    Paper 3 calling convention (Script 14):
      - Script: perl ~/netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl
      - Env: NETMHCIIpan=~/netmhc/netMHCIIpan-4.0/
             NetMHCIIpanPLAT=~/netmhc/netMHCIIpan-4.0/Linux_x86_64
      - Flags: -f <fasta> -a <alleles>  (default 15-mer)
      - Threshold: %Rank <= 10.0 (NOT IC50 nM)
      - Returns: unique start positions of binding 15-mers

    Returns list of {peptide, start, allele, rank_el}.

    Raises RuntimeError if Perl script is not found.
    """
    import os
    import tempfile

    _alleles = alleles or MHC_II_ALLELES_PAPER3

    home = Path.home()
    _vm_perl = str(home / "netmhc/netMHCIIpan-4.0/NetMHCIIpan-4.0.pl")
    _perl_script = netmhciipan_perl or _vm_perl
    if not Path(_perl_script).exists():
        raise RuntimeError(
            f"NetMHCIIpan-4.0 Perl script not found at {_perl_script!r}. "
            "Run inside pen-stack VM or biophysics Docker image."
        )

    env = os.environ.copy()
    env["NETMHCIIpan"] = str(home / "netmhc/netMHCIIpan-4.0")
    env["NetMHCIIpanPLAT"] = str(home / "netmhc/netMHCIIpan-4.0/Linux_x86_64")
    env["TMPDIR"] = "/tmp"   # NetMHCIIpan-4.0 also requires $TMPDIR

    allele_str = ",".join(_alleles)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as fh:
        fh.write(f">query\n{sequence}\n")
        fasta_path = fh.name

    try:
        result = subprocess.run(
            ["perl", _perl_script, "-f", fasta_path, "-a", allele_str],
            capture_output=True, text=True, timeout=600, env=env,
        )
    finally:
        Path(fasta_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"NetMHCIIpan-4.0 failed (rc={result.returncode}): {result.stderr[:400]}"
        )

    binders: list[dict] = []
    seen: set[tuple] = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split()
        # NetMHCIIpan-4.0 output columns (space-split, verified from actual output):
        # [0]=Pos [1]=MHC [2]=Peptide [3]=Of [4]=Core [5]=Core_Rel
        # [6]=Identity [7]=Score_EL [8]=%Rank_EL [9]=Exp_Bind [10]=BindLevel(opt)
        if len(parts) < 9:
            continue
        if not parts[0].lstrip("-").isdigit():
            continue
        try:
            pos = int(parts[0])
            allele = parts[1]
            peptide = parts[2]
            rank_el = float(parts[8])   # %Rank_EL at index 8 (verified)
        except (ValueError, IndexError):
            continue
        if rank_el <= rank_el_threshold:
            key = (pos, allele)
            if key not in seen:
                seen.add(key)
                binders.append({
                    "peptide": peptide,
                    "start": pos,
                    "allele": allele,
                    "rank_el": rank_el,
                })

    return binders


def compute_weighted_epitope_load(
    mhc_i_binders: list[dict],
    mhc_ii_binders: list[dict],
    weight_mhc_ii: float = 0.5,
) -> float:
    """Weighted epitope count: n_MHC_I + weight_mhc_ii * n_MHC_II.

    Matches Paper 3 Strategy C optimization target (Script 14 formula).
    Unique positions are counted (deduplicated within each class).
    """
    # Count unique start positions (already deduplicated in find_mhc_*_binders)
    n_i = len({b["start"] for b in mhc_i_binders})
    n_ii = len({b["start"] for b in mhc_ii_binders})
    return n_i + weight_mhc_ii * n_ii


def compute_s_immuno_paper3(
    n_i: int,
    n_ii: int,
    max_total: float = S_IMMUNO_MAX_TOTAL,
) -> float:
    """Compute S_Immuno using Paper 3 absolute normalization formula.

    S_Immuno = 1 - (n_I + 0.5 * n_II) / max_total

    Parameters
    ----------
    n_i : int
        Number of unique MHC-I binding positions (%Rank_EL < 0.5).
    n_ii : int
        Number of unique MHC-II binding positions (%Rank_EL <= 10.0).
    max_total : float
        Normalization constant = 253.5 (95th-pct across 28 editors, locked Paper 3).

    Returns
    -------
    float
        S_Immuno in [0, 1]; higher = less immunogenic.
        IS621 WT: n_I=15, n_II=92 -> S_Immuno=0.7594.
    """
    combined = n_i + 0.5 * n_ii
    return round(min(1.0, max(0.0, 1.0 - combined / max_total)), 4)
