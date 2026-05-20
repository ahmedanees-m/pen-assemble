# Wet-Lab Reference: `D030_IS621_ProtMPNN_T0.1_sample7`

**✓ Beats IS621 verbatim (0.929)**

## Design Summary

| Field | Value |
|-------|-------|
| Design ID | `D030_IS621_ProtMPNN_T0.1_sample7` |
| Strategy | Strategy D (IS621 ProtMPNN sequence design) |
| PenScore | **0.9311** (IS621 = 0.929) |
| Length | 342 aa &nbsp;≈&nbsp; 38 kDa |
| Source organism | *IS621 (Firmicutes bacterium)* |
| Genus | *—* |
| Protein name | Bridge recombinase / IS110-family transposase |

## PenScore Axis Breakdown

| Axis | Score | Weight | Contribution |
|------|-------|--------|-------------|
| S_DSB (double-strand break) | 1.0000 | 0.25 | 0.2500 |
| S_Spec (target specificity) | 0.9531 | 0.10 | 0.0953 |
| S_Cargo (payload compatibility) | 1.0000 | 0.20 | 0.2000 |
| S_Deliv (delivery suitability) | 0.9552 | 0.15 | 0.1433 |
| S_Immuno (de-immunization) | 0.7494 | 0.10 | 0.0749 |
| S_Prog (programmability) | 0.9851 | 0.15 | 0.1478 |
| S_Mature (maturity / TRL) | 0.3960 | 0.05 | 0.0198 |

## Computational Structure Quality

| Metric | Value |
|--------|-------|
| ESMFold mean pLDDT | 92.0 (threshold > 70) |
| Active-site pLDDT  | 93.7 (threshold > 70) |
| Rosetta ΔΔG        | -40652.7 kcal/mol *(cross-protein absolute energy; not true ΔΔG — see Deviation 5)* |
| PDB path (VM)      | `/home/anees_22phd0670/esm_tier1_output/pdbs/D030_IS621_ProtMPNN_T0.1_sample7.pdb` |

## Codon-Optimised DNA (Human Expression)

> **Codons**: Kazusa Homo sapiens preferred codons — rule-based, not CAI-maximised.  
> **GC content**: 68.0%  
> **Restriction sites**: None of the common 8 screened.  
> **ORF length**: 1026 bp  
> **Recommendation**: Add Kozak (GCCACC) before ATG; add stop codon (TGA appended).  

### Amino acid sequence
```
MDRFFPVIRICKVGFTMEHTTVYIGIDTAKETLTVCVLLPDGRHRTATFANNAAGYEALV
AWLRELGITDADVVIEATGTYMEPVARTLIAAGYRVAVINPALFKAFRQSEGRRNKTDTV
DARALALFGRQKRPPAWTPPPPLEEELRALVVEHQRLTDMHTQVLNRLETARPEERPFLE
AHLLWLDARLADIAARIRELVAADPVLARKKALIESIPGIGEKTSAVALAFIGLEQRFAH
ARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVDFRRSLYMPAMVFVSKTKVGRKFAERL
RKNGKKGKVILGAMMRKLAQVIYGVLKSGVPFDASRHNPVAA
```

### Codon-optimised ORF (no Kozak, no stop — append as needed)
```
ATGGACAGGTTCTTCCCCGTGATCAGGATCTGCAAGGTGGGCTTCACCATGGAGCACACC
ACCGTGTACATCGGCATCGACACCGCCAAGGAGACCCTGACCGTGTGCGTGCTGCTGCCC
GACGGCAGGCACAGGACCGCCACCTTCGCCAACAACGCCGCCGGCTACGAGGCCCTGGTG
GCCTGGCTGAGGGAGCTGGGCATCACCGACGCCGACGTGGTGATCGAGGCCACCGGCACC
TACATGGAGCCCGTGGCCAGGACCCTGATCGCCGCCGGCTACAGGGTGGCCGTGATCAAC
CCCGCCCTGTTCAAGGCCTTCAGGCAGAGCGAGGGCAGGAGGAACAAGACCGACACCGTG
GACGCCAGGGCCCTGGCCCTGTTCGGCAGGCAGAAGAGGCCCCCCGCCTGGACCCCCCCC
CCCCCCCTGGAGGAGGAGCTGAGGGCCCTGGTGGTGGAGCACCAGAGGCTGACCGACATG
CACACCCAGGTGCTGAACAGGCTGGAGACCGCCAGGCCCGAGGAGAGGCCCTTCCTGGAG
GCCCACCTGCTGTGGCTGGACGCCAGGCTGGCCGACATCGCCGCCAGGATCAGGGAGCTG
GTGGCCGCCGACCCCGTGCTGGCCAGGAAGAAGGCCCTGATCGAGAGCATCCCCGGCATC
GGCGAGAAGACCAGCGCCGTGGCCCTGGCCTTCATCGGCCTGGAGCAGAGGTTCGCCCAC
GCCAGGCAGTTCGCCGCCTTCGCCGGCCTGACCCCCAGGAGGTACGAGAGCGGCAGCAGC
GTGAGGGGCGCCAGCAGGATGAGCAAGGCCGGCCACGTGGACTTCAGGAGGAGCCTGTAC
ATGCCCGCCATGGTGTTCGTGAGCAAGACCAAGGTGGGCAGGAAGTTCGCCGAGAGGCTG
AGGAAGAACGGCAAGAAGGGCAAGGTGATCCTGGGCGCCATGATGAGGAAGCTGGCCCAG
GTGATCTACGGCGTGCTGAAGAGCGGCGTGCCCTTCGACGCCAGCAGGCACAACCCCGTG
GCCGCC
```

### With Kozak + stop (ready-to-order)
```
GCCACCATGGACAGGTTCTTCCCCGTGATCAGGATCTGCAAGGTGGGCTTCACCATGGAG
CACACCACCGTGTACATCGGCATCGACACCGCCAAGGAGACCCTGACCGTGTGCGTGCTG
CTGCCCGACGGCAGGCACAGGACCGCCACCTTCGCCAACAACGCCGCCGGCTACGAGGCC
CTGGTGGCCTGGCTGAGGGAGCTGGGCATCACCGACGCCGACGTGGTGATCGAGGCCACC
GGCACCTACATGGAGCCCGTGGCCAGGACCCTGATCGCCGCCGGCTACAGGGTGGCCGTG
ATCAACCCCGCCCTGTTCAAGGCCTTCAGGCAGAGCGAGGGCAGGAGGAACAAGACCGAC
ACCGTGGACGCCAGGGCCCTGGCCCTGTTCGGCAGGCAGAAGAGGCCCCCCGCCTGGACC
CCCCCCCCCCCCCTGGAGGAGGAGCTGAGGGCCCTGGTGGTGGAGCACCAGAGGCTGACC
GACATGCACACCCAGGTGCTGAACAGGCTGGAGACCGCCAGGCCCGAGGAGAGGCCCTTC
CTGGAGGCCCACCTGCTGTGGCTGGACGCCAGGCTGGCCGACATCGCCGCCAGGATCAGG
GAGCTGGTGGCCGCCGACCCCGTGCTGGCCAGGAAGAAGGCCCTGATCGAGAGCATCCCC
GGCATCGGCGAGAAGACCAGCGCCGTGGCCCTGGCCTTCATCGGCCTGGAGCAGAGGTTC
GCCCACGCCAGGCAGTTCGCCGCCTTCGCCGGCCTGACCCCCAGGAGGTACGAGAGCGGC
AGCAGCGTGAGGGGCGCCAGCAGGATGAGCAAGGCCGGCCACGTGGACTTCAGGAGGAGC
CTGTACATGCCCGCCATGGTGTTCGTGAGCAAGACCAAGGTGGGCAGGAAGTTCGCCGAG
AGGCTGAGGAAGAACGGCAAGAAGGGCAAGGTGATCCTGGGCGCCATGATGAGGAAGCTG
GCCCAGGTGATCTACGGCGTGCTGAAGAGCGGCGTGCCCTTCGACGCCAGCAGGCACAAC
CCCGTGGCCGCCTGA
```

### Recommended Validation Experiments

#### Tier 1 — In vitro (priority)
1. **Recombination assay (attB × attP)**: Incubate purified protein (0.5–5 µM) with
   supercoiled plasmid carrying cognate attB/attP sites. Resolve on 1% agarose gel.
   Positive: appearance of lower-molecular-weight relaxed/linear product.
2. **EMSA (electrophoretic mobility shift)**: Titrate protein against Cy5-attB dsDNA
   (40 bp). Confirm specific binding (Kd target < 500 nM).
3. **Thermal stability (nanoDSF)**: Confirm Tm > 40 °C (minimum for activity at 37 °C).
   Target Tm ≥ 50 °C.

#### Tier 2 — Cell-based
4. **HEK293T transient transfection**: Co-transfect codon-optimised ORF (in pCMV-FLAG)
   with dual-reporter plasmid (mCherry-attB-attP-EGFP). Gate on mCherry+ cells; score
   EGFP+ fraction by flow cytometry at 48 h. Threshold: > 5% recombination above background.
5. **Western blot**: Anti-FLAG; confirm expected MW (see table above); flag truncations.
6. **Immunofluorescence**: Confirm nuclear localisation if NLS is appended; cytoplasmic
   distribution is acceptable for HDR-coupled delivery.

#### Tier 3 — Deep characterisation (after Tier 1–2 pass)
7. **Specificity panel**: Test against 10 scrambled attB sequences; confirm < 1%
   off-target recombination.
8. **Dose–response**: 0.1–10 µM protein; fit Hill equation; report EC50.
9. **T7E1 assay**: Rule out NHEJ at attB target site (< 0.5% indel rate threshold).

---

> **COMPUTATIONAL DISCLAIMER**
> This sequence was designed computationally using the PEN-ASSEMBLE v0.5.0 pipeline.
> PenScore is a composite *in silico* metric — it does not guarantee biological activity.
> Stability assessment used ESMFold pLDDT as a proxy (Rosetta ΔΔG gate was non-functional
> for all designs; see Deviation 5 in DESIGN_PROVENANCE.md).
> All IS110/bridge-recombinase activity claims require experimental validation.
> Codon optimisation is rule-based (Kazusa preferred codons); verify codon-adaptation
> index (CAI) with a commercial tool before synthesis order.
> This file is for research use only.
