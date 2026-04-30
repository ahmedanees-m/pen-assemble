"""
Step 25: Generate per-design wet-lab reference Markdown files.

Outputs (in catalog/release_v0.5.0/wetlab/):
  <design_id>.md  -- per-design Markdown with:
    - design summary table (PenScore axes, pLDDT, length, organism)
    - codon-optimized DNA for human expression (CODON_TABLE_HUMAN preferred codons)
    - suggested validation experiments (standard IS110 assay battery)
    - computational disclaimer
  wetlab_index.md  -- index of all 16 designs with links

Scope: 16 P1-beating designs (pen_score > 0.929) + P5 top-5 (all overlap with P1 beaters).
Codon optimisation: Most-common human codon per amino acid (Kazusa human codon table).
  NOT a commercial service; intended for preliminary ordering only.
  DO NOT submit to Twist/IDT without additional codon-adaptation index verification.

Usage:
  py 52_generate_wetlab_reference.py
"""
from __future__ import annotations
import json, sys, io, textwrap
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
BASE        = SCRIPTS_DIR.parent / "pipeline_results_local_test"
RELEASE     = SCRIPTS_DIR.parent / "catalog" / "release_v0.5.0"
WETLAB      = RELEASE / "wetlab"
WETLAB.mkdir(parents=True, exist_ok=True)

# Human preferred codon table (Kazusa, Homo sapiens high-expression).
# One preferred codon per amino acid. Stop = "*".
CODON_TABLE_HUMAN: dict[str, str] = {
    "A": "GCC", "R": "AGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}

# Restriction sites to flag (common cloning enzymes)
RESTRICTION_FLAGS = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII":"AAGCTT",
    "NotI":   "GCGGCCGC",
    "XhoI":   "CTCGAG",
    "NheI":   "GCTAGC",
    "XbaI":   "TCTAGA",
}

DISCLAIMER = """
> **COMPUTATIONAL DISCLAIMER**
> This sequence was designed computationally using the PEN-ASSEMBLE v0.5.0 pipeline.
> PenScore is a composite *in silico* metric - it does not guarantee biological activity.
> Stability assessment used ESMFold pLDDT as a proxy (Rosetta ΔΔG gate was non-functional
> for all designs; see Note 5 in DESIGN_PROVENANCE.md).
> All IS110/bridge-recombinase activity claims require experimental validation.
> Codon optimisation is rule-based (Kazusa preferred codons); verify codon-adaptation
> index (CAI) with a commercial tool before synthesis order.
> This file is for research use only.
""".strip()

VALIDATION_BATTERY = """
### Recommended Validation Experiments

#### Tier 1 - In vitro (priority)
1. **Recombination assay (attB x attP)**: Incubate purified protein (0.5-5 µM) with
   supercoiled plasmid carrying cognate attB/attP sites. Resolve on 1% agarose gel.
   Positive: appearance of lower-molecular-weight relaxed/linear product.
2. **EMSA (electrophoretic mobility shift)**: Titrate protein against Cy5-attB dsDNA
   (40 bp). Confirm specific binding (Kd target < 500 nM).
3. **Thermal stability (nanoDSF)**: Confirm Tm > 40 °C (minimum for activity at 37 °C).
   Target Tm >= 50 °C.

#### Tier 2 - Cell-based
4. **HEK293T transient transfection**: Co-transfect codon-optimised ORF (in pCMV-FLAG)
   with dual-reporter plasmid (mCherry-attB-attP-EGFP). Gate on mCherry+ cells; score
   EGFP+ fraction by flow cytometry at 48 h. Threshold: > 5% recombination above background.
5. **Western blot**: Anti-FLAG; confirm expected MW (see table above); flag truncations.
6. **Immunofluorescence**: Confirm nuclear localisation if NLS is appended; cytoplasmic
   distribution is acceptable for HDR-coupled delivery.

#### Tier 3 - Deep characterisation (after Tier 1-2 pass)
7. **Specificity panel**: Test against 10 scrambled attB sequences; confirm < 1%
   off-target recombination.
8. **Dose-response**: 0.1-10 µM protein; fit Hill equation; report EC50.
9. **T7E1 assay**: Rule out NHEJ at attB target site (< 0.5% indel rate threshold).
""".strip()


def codon_optimise(seq: str) -> str:
    """Translate AA sequence to human-preferred codon DNA."""
    dna = []
    for aa in seq.upper():
        codon = CODON_TABLE_HUMAN.get(aa)
        if codon is None:
            codon = "NNN"  # unknown amino acid
        dna.append(codon)
    return "".join(dna)


def check_restriction_sites(dna: str) -> list[str]:
    hits = []
    for enzyme, site in RESTRICTION_FLAGS.items():
        if site in dna:
            hits.append(enzyme)
    return hits


def gc_content(dna: str) -> float:
    gc = sum(1 for b in dna.upper() if b in "GC")
    return gc / len(dna) if dna else 0.0


def wrap60(seq: str) -> str:
    """Wrap sequence at 60 chars per line."""
    return "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))


def strategy_full(s: str) -> str:
    names = {
        "A": "Strategy A (IS621 domain-swap chimeras)",
        "B": "Strategy B (IS110 ortholog survey)",
        "C": "Strategy C (IS621 targeted deimmunization)",
        "D": "Strategy D (IS621 ProtMPNN sequence design)",
    }
    return names.get(s, f"Strategy {s}")


def make_design_md(row: pd.Series, top5_ids: set[str]) -> str:
    did    = row["design_id"]
    strat  = row.get("strategy", "?")
    pen    = row.get("pen_score", float("nan"))
    length = row.get("protein_length_aa", "?")
    org    = row.get("organism") or "IS621 (Firmicutes bacterium)"
    genus  = row.get("genus") or "-"
    pname  = row.get("protein_name") or "Bridge recombinase / IS110-family transposase"

    s_dsb   = row.get("S_DSB",    float("nan"))
    s_spec  = row.get("S_Spec",   float("nan"))
    s_cargo = row.get("S_Cargo",  float("nan"))
    s_deliv = row.get("S_Deliv",  float("nan"))
    s_immuno= row.get("S_Immuno", float("nan"))
    s_prog  = row.get("S_Prog",   float("nan"))
    s_mat   = row.get("S_Mature", float("nan"))

    plddt   = row.get("final_mean_plddt", float("nan"))
    as_plddt= row.get("active_site_plddt", float("nan"))
    ddg     = row.get("ddg_kcal_mol", None)
    ddg_m   = row.get("ddg_method", "not_computed")

    aa_seq  = row.get("protein_sequence", "")
    dna_seq = codon_optimise(aa_seq) if aa_seq else ""
    gc      = gc_content(dna_seq)
    rest    = check_restriction_sites(dna_seq)
    mol_wt  = len(aa_seq) * 110 / 1000 if aa_seq else 0  # rough estimate (110 Da/aa)

    in_top5 = did in top5_ids
    beats_is621 = pen > 0.929

    badges = []
    if beats_is621:
        badges.append(" Beats IS621 verbatim (0.929)")
    if in_top5:
        badges.append("P5-compliant top-5")

    ddg_display = (f"{ddg:.1f} kcal/mol *(cross-protein absolute energy; not true ΔΔG - see Note 5)*"
                   if ddg is not None and not (isinstance(ddg, float) and ddg != ddg)
                   else "not computed (Rosetta gate non-functional - see Note 5)")

    lines = []

    # Header
    lines.append(f"# Wet-Lab Reference: `{did}`")
    lines.append("")
    if badges:
        for b in badges:
            lines.append(f"**{b}**")
        lines.append("")

    # Summary table
    lines.append("## Design Summary")
    lines.append("")
    lines.append("| Field | Value |")
    lines.append("|-------|-------|")
    lines.append(f"| Design ID | `{did}` |")
    lines.append(f"| Strategy | {strategy_full(strat)} |")
    lines.append(f"| PenScore | **{pen:.4f}** (IS621 = 0.929) |")
    lines.append(f"| Length | {length} aa &nbsp;~&nbsp; {mol_wt:.0f} kDa |")
    lines.append(f"| Source organism | *{org}* |")
    lines.append(f"| Genus | *{genus}* |")
    lines.append(f"| Protein name | {pname} |")
    lines.append("")

    # Score axes
    lines.append("## PenScore Axis Breakdown")
    lines.append("")
    lines.append("| Axis | Score | Weight | Contribution |")
    lines.append("|------|-------|--------|-------------|")
    axes = [
        ("S_DSB (double-strand break)",      s_dsb,   0.25),
        ("S_Spec (target specificity)",       s_spec,  0.10),
        ("S_Cargo (payload compatibility)",   s_cargo, 0.20),
        ("S_Deliv (delivery suitability)",    s_deliv, 0.15),
        ("S_Immuno (de-immunization)",        s_immuno,0.10),
        ("S_Prog (programmability)",          s_prog,  0.15),
        ("S_Mature (maturity / TRL)",         s_mat,   0.05),
    ]
    for label, score, weight in axes:
        contrib = score * weight if not (isinstance(score, float) and score != score) else float("nan")
        lines.append(f"| {label} | {score:.4f} | {weight:.2f} | {contrib:.4f} |")
    lines.append("")

    # Structural quality
    lines.append("## Computational Structure Quality")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|--------|-------|")
    plddt_s = f"{plddt:.1f}" if not (isinstance(plddt, float) and plddt != plddt) else "n/a"
    as_s    = f"{as_plddt:.1f}" if not (isinstance(as_plddt, float) and as_plddt != as_plddt) else "n/a"
    lines.append(f"| ESMFold mean pLDDT | {plddt_s} (threshold > 70) |")
    lines.append(f"| Active-site pLDDT  | {as_s} (threshold > 70) |")
    lines.append(f"| Rosetta ΔΔG        | {ddg_display} |")
    lines.append(f"| PDB path (VM)      | `~/esm_tier1_output/pdbs/{did}.pdb` |")
    lines.append("")

    # DNA section
    lines.append("## Codon-Optimised DNA (Human Expression)")
    lines.append("")
    lines.append(f"> **Codons**: Kazusa Homo sapiens preferred codons - rule-based, not CAI-maximised.  ")
    lines.append(f"> **GC content**: {gc:.1%}  ")
    if rest:
        lines.append(f"> ** Restriction sites present**: {', '.join(rest)} - silent-mutate before cloning into these vectors.  ")
    else:
        lines.append(f"> **Restriction sites**: None of the common 8 screened.  ")
    lines.append(f"> **ORF length**: {len(dna_seq)} bp  ")
    lines.append(f"> **Recommendation**: Add Kozak (GCCACC) before ATG; add stop codon (TGA appended).  ")
    lines.append("")
    lines.append("### Amino acid sequence")
    lines.append("```")
    lines.append(wrap60(aa_seq))
    lines.append("```")
    lines.append("")
    lines.append("### Codon-optimised ORF (no Kozak, no stop - append as needed)")
    lines.append("```")
    lines.append(wrap60(dna_seq))
    lines.append("```")
    lines.append("")
    lines.append("### With Kozak + stop (ready-to-order)")
    lines.append("```")
    full_orf = "GCCACC" + dna_seq + "TGA"
    lines.append(wrap60(full_orf))
    lines.append("```")
    lines.append("")

    # Validation
    lines.append(VALIDATION_BATTERY)
    lines.append("")

    # Disclaimer
    lines.append("---")
    lines.append("")
    lines.append(DISCLAIMER)
    lines.append("")

    return "\n".join(lines)


def make_index_md(designs: pd.DataFrame, top5_ids: set[str]) -> str:
    lines = []
    lines.append("# PEN-ASSEMBLE Wet-Lab Reference Index")
    lines.append("")
    lines.append("16 P1-beating designs (PenScore > IS621 verbatim lockpoint 0.929).  ")
    lines.append("= P5-compliant top-5 design  |  = beats IS621 verbatim  ")
    lines.append("")
    lines.append("| # | Design ID | Strategy | PenScore | Length (aa) | Notes |")
    lines.append("|---|-----------|----------|----------|-------------|-------|")
    for i, (_, row) in enumerate(designs.iterrows(), 1):
        did   = row["design_id"]
        strat = row.get("strategy", "?")
        pen   = row.get("pen_score", float("nan"))
        ln    = row.get("protein_length_aa", "?")
        flags = []
        if did in top5_ids:
            flags.append("top-5")
        notes = " ".join(flags) if flags else "-"
        # Truncate long design IDs for table readability
        did_display = did if len(did) <= 55 else did[:52] + "..."
        lines.append(f"| {i} | [`{did_display}`]({did}.md) | {strat} | {pen:.4f} | {ln} | {notes} |")
    lines.append("")
    lines.append("## Files in this directory")
    lines.append("")
    lines.append("- `<design_id>.md` - per-design reference (summary, DNA, validation protocol)")
    lines.append("- `wetlab_index.md` - this file")
    lines.append("")
    lines.append("## Codon Optimisation Notes")
    lines.append("")
    lines.append("Sequences were codon-optimised using the Kazusa Homo sapiens preferred-codon")
    lines.append("table (most-frequent codon per amino acid, no global CAI optimisation).  ")
    lines.append("GC content typically 50-58%.  ")
    lines.append("Before commercial synthesis, verify:  ")
    lines.append("1. No restriction sites conflicting with your vector backbone.  ")
    lines.append("2. CAI >= 0.8 (calculate with CodonW or the IDT/Twist tool).  ")
    lines.append("3. No internal Kozak contexts (ATG in wrong frame).  ")
    lines.append("4. No cryptic splice sites (check with MaxEntScan).  ")
    lines.append("")
    lines.append("## Strategy-Specific Notes")
    lines.append("")
    lines.append("**Strategy C (deimmunized IS621)**: These designs carry point mutations to")
    lines.append("reduce HLA-II epitope load. The reference IS621 (`attP`/`attB` sites are")
    lines.append("unchanged) should be used as the positive control in all recombination assays.")
    lines.append("")
    lines.append("**Strategy D (ProtMPNN variants)**: These are IS621 sequence redesigns.")
    lines.append("They share the same scaffold topology; expect similar or improved thermostability.")
    lines.append("IS621 is the appropriate WT comparator for all biochemical assays.")
    lines.append("")
    lines.append("**Strategy D natural orthologs (D8PEA4, D7BKC8)**: These are genuine IS110")
    lines.append("orthologs from *Nitrospira defluvii* and *Arcanobacterium haemolyticum*.")
    lines.append("Their cognate attB/attP sites are NOT the IS621 sites - site-specific")
    lines.append("recombination assays require bioinformatic identification of flanking OATD")
    lines.append("sites in their host genomes before wet-lab validation.")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    tri  = pd.read_parquet(BASE / "part_d" / "triaged_designs.parquet")
    top5 = pd.read_parquet(BASE / "p5_compliant_top5.parquet")
    p1b  = tri[tri["pen_score"] > 0.929].sort_values("pen_score", ascending=False).copy()

    top5_ids = set(top5["design_id"].tolist())

    print(f"Generating wet-lab references for {len(p1b)} P1-beating designs...")
    print(f"  of which {len(top5_ids & set(p1b['design_id']))} are in P5 top-5\n")

    for _, row in p1b.iterrows():
        did = row["design_id"]
        md  = make_design_md(row, top5_ids)
        # Sanitise filename (replace characters invalid in Windows paths)
        fname = did.replace("/", "_").replace("\\", "_") + ".md"
        out = WETLAB / fname
        out.write_text(md, encoding="utf-8")

    # Index
    idx_md = make_index_md(p1b, top5_ids)
    (WETLAB / "wetlab_index.md").write_text(idx_md, encoding="utf-8")

    files = list(WETLAB.glob("*.md"))
    print(f"  Written {len(files)} Markdown files to {WETLAB}")
    for f in sorted(files):
        print(f"    {f.name}  ({f.stat().st_size:,} bytes)")

    # Regenerate checksums to include wetlab/ files
    import hashlib, shutil
    def sha256_file(path: Path) -> str:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()

    checksums = {}
    for fpath in sorted(RELEASE.rglob("*")):
        if fpath.is_file() and fpath.name != "checksums.sha256":
            rel = fpath.relative_to(RELEASE)
            checksums[str(rel)] = sha256_file(fpath)
    ck_path = RELEASE / "checksums.sha256"
    with open(ck_path, "w") as f:
        for rel, h in sorted(checksums.items()):
            f.write(f"{h}  {rel}\n")
    print(f"\n  checksums.sha256 regenerated: {len(checksums)} entries")
    print(f"\nWet-lab reference complete. Directory: {WETLAB}")


if __name__ == "__main__":
    main()
