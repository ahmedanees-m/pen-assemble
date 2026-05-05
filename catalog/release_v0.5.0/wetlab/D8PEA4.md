# Wet-Lab Reference: `D8PEA4`

** Beats IS621 verbatim (0.929)**
**P5-compliant top-5**

## Design Summary

| Field | Value |
|-------|-------|
| Design ID | `D8PEA4` |
| Strategy | Strategy D (IS621 ProtMPNN sequence design) |
| PenScore | **0.9353** (IS621 = 0.929) |
| Length | 314 aa &nbsp;~&nbsp; 35 kDa |
| Source organism | *Nitrospira defluvii* |
| Genus | *-* |
| Protein name | Transposase |

## PenScore Axis Breakdown

| Axis | Score | Weight | Contribution |
|------|-------|--------|-------------|
| S_DSB (double-strand break) | 1.0000 | 0.25 | 0.2500 |
| S_Spec (target specificity) | 0.9549 | 0.10 | 0.0955 |
| S_Cargo (payload compatibility) | 1.0000 | 0.20 | 0.2000 |
| S_Deliv (delivery suitability) | 0.9705 | 0.15 | 0.1456 |
| S_Immuno (de-immunization) | 0.7452 | 0.10 | 0.0745 |
| S_Prog (programmability) | 0.9993 | 0.15 | 0.1499 |
| S_Mature (maturity / TRL) | 0.3960 | 0.05 | 0.0198 |

## Computational Structure Quality

| Metric | Value |
|--------|-------|
| ESMFold mean pLDDT | 94.3 (threshold > 70) |
| Active-site pLDDT  | 97.0 (threshold > 70) |
| Rosetta ΔΔG        | not computed (Rosetta gate non-functional - see Note 5) |
| PDB path (VM)      | `~/esm_tier1_output/pdbs/D8PEA4.pdb` |

## Codon-Optimised DNA (Human Expression)

> **Codons**: Kazusa Homo sapiens preferred codons - rule-based, not CAI-maximised.  
> **GC content**: 71.1%  
> **Restriction sites**: None of the common 8 screened.  
> **ORF length**: 942 bp  
> **Recommendation**: Add Kozak (GCCACC) before ATG; add stop codon (TGA appended).  

### Amino acid sequence
```
MSAEIFVGMDVSQGGVDVAVQPGTAFQIAHNERGIAEAVARLQAVQPTLIVLEATGGLEV
PLTGALAAAGLPVVVINPRQVRDFARATGPLAKTDRLEAQILARFAEAIRPPVRPVPDEQ
TQALAALVARRRQLIEMLTAEKNRLRLAARPIQKRVQAHVTWLEKELASTNTDLTATIRE
SPVWRAKADVLRSVPGVGPVLTTTLFANLPELGTLTRKEVAALAGVAPFPRDSGTLKGRR
TIWGGRAHVRAALYMAALVATRRNSVIRAFYQRLCQAGKAKKLALTACMRKLLTILNAML
KHGTRWRVTASQPA
```

### Codon-optimised ORF (no Kozak, no stop - append as needed)
```
ATGAGCGCCGAGATCTTCGTGGGCATGGACGTGAGCCAGGGCGGCGTGGACGTGGCCGTG
CAGCCCGGCACCGCCTTCCAGATCGCCCACAACGAGAGGGGCATCGCCGAGGCCGTGGCC
AGGCTGCAGGCCGTGCAGCCCACCCTGATCGTGCTGGAGGCCACCGGCGGCCTGGAGGTG
CCCCTGACCGGCGCCCTGGCCGCCGCCGGCCTGCCCGTGGTGGTGATCAACCCCAGGCAG
GTGAGGGACTTCGCCAGGGCCACCGGCCCCCTGGCCAAGACCGACAGGCTGGAGGCCCAG
ATCCTGGCCAGGTTCGCCGAGGCCATCAGGCCCCCCGTGAGGCCCGTGCCCGACGAGCAG
ACCCAGGCCCTGGCCGCCCTGGTGGCCAGGAGGAGGCAGCTGATCGAGATGCTGACCGCC
GAGAAGAACAGGCTGAGGCTGGCCGCCAGGCCCATCCAGAAGAGGGTGCAGGCCCACGTG
ACCTGGCTGGAGAAGGAGCTGGCCAGCACCAACACCGACCTGACCGCCACCATCAGGGAG
AGCCCCGTGTGGAGGGCCAAGGCCGACGTGCTGAGGAGCGTGCCCGGCGTGGGCCCCGTG
CTGACCACCACCCTGTTCGCCAACCTGCCCGAGCTGGGCACCCTGACCAGGAAGGAGGTG
GCCGCCCTGGCCGGCGTGGCCCCCTTCCCCAGGGACAGCGGCACCCTGAAGGGCAGGAGG
ACCATCTGGGGCGGCAGGGCCCACGTGAGGGCCGCCCTGTACATGGCCGCCCTGGTGGCC
ACCAGGAGGAACAGCGTGATCAGGGCCTTCTACCAGAGGCTGTGCCAGGCCGGCAAGGCC
AAGAAGCTGGCCCTGACCGCCTGCATGAGGAAGCTGCTGACCATCCTGAACGCCATGCTG
AAGCACGGCACCAGGTGGAGGGTGACCGCCAGCCAGCCCGCC
```

### With Kozak + stop (ready-to-order)
```
GCCACCATGAGCGCCGAGATCTTCGTGGGCATGGACGTGAGCCAGGGCGGCGTGGACGTG
GCCGTGCAGCCCGGCACCGCCTTCCAGATCGCCCACAACGAGAGGGGCATCGCCGAGGCC
GTGGCCAGGCTGCAGGCCGTGCAGCCCACCCTGATCGTGCTGGAGGCCACCGGCGGCCTG
GAGGTGCCCCTGACCGGCGCCCTGGCCGCCGCCGGCCTGCCCGTGGTGGTGATCAACCCC
AGGCAGGTGAGGGACTTCGCCAGGGCCACCGGCCCCCTGGCCAAGACCGACAGGCTGGAG
GCCCAGATCCTGGCCAGGTTCGCCGAGGCCATCAGGCCCCCCGTGAGGCCCGTGCCCGAC
GAGCAGACCCAGGCCCTGGCCGCCCTGGTGGCCAGGAGGAGGCAGCTGATCGAGATGCTG
ACCGCCGAGAAGAACAGGCTGAGGCTGGCCGCCAGGCCCATCCAGAAGAGGGTGCAGGCC
CACGTGACCTGGCTGGAGAAGGAGCTGGCCAGCACCAACACCGACCTGACCGCCACCATC
AGGGAGAGCCCCGTGTGGAGGGCCAAGGCCGACGTGCTGAGGAGCGTGCCCGGCGTGGGC
CCCGTGCTGACCACCACCCTGTTCGCCAACCTGCCCGAGCTGGGCACCCTGACCAGGAAG
GAGGTGGCCGCCCTGGCCGGCGTGGCCCCCTTCCCCAGGGACAGCGGCACCCTGAAGGGC
AGGAGGACCATCTGGGGCGGCAGGGCCCACGTGAGGGCCGCCCTGTACATGGCCGCCCTG
GTGGCCACCAGGAGGAACAGCGTGATCAGGGCCTTCTACCAGAGGCTGTGCCAGGCCGGC
AAGGCCAAGAAGCTGGCCCTGACCGCCTGCATGAGGAAGCTGCTGACCATCCTGAACGCC
ATGCTGAAGCACGGCACCAGGTGGAGGGTGACCGCCAGCCAGCCCGCCTGA
```

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

---

> **COMPUTATIONAL DISCLAIMER**
> This sequence was designed computationally using the PEN-ASSEMBLE v0.5.0 pipeline.
> PenScore is a composite *in silico* metric - it does not guarantee biological activity.
> Stability assessment used ESMFold pLDDT as a proxy (Rosetta ΔΔG gate was non-functional
> for all designs; see Note 5 in DESIGN_PROVENANCE.md).
> All IS110/bridge-recombinase activity claims require experimental validation.
> Codon optimisation is rule-based (Kazusa preferred codons); verify codon-adaptation
> index (CAI) with a commercial tool before synthesis order.
> This file is for research use only.
