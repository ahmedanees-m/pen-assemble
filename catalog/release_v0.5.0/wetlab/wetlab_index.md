# PEN-ASSEMBLE Wet-Lab Reference Index

16 P1-beating designs (PenScore > IS621 verbatim lockpoint 0.929).  
★ = P5-compliant top-5 design  |  ✓ = beats IS621 verbatim  

| # | Design ID | Strategy | PenScore | Length (aa) | Notes |
|---|-----------|----------|----------|-------------|-------|
| 1 | [`IS621_deimmunized_v2_Y255K_D41C_G65K_E187M_D27C_D203...`](IS621_deimmunized_v2_Y255K_D41C_G65K_E187M_D27C_D203C_V285C_V152C_T224C_L318K_E87C_L193I_L275I_P177V.md) | C | 0.9673 | 342 | ★ top-5 |
| 2 | [`C_targeted_001`](C_targeted_001.md) | C | 0.9586 | 342 | ★ top-5 |
| 3 | [`D8PEA4`](D8PEA4.md) | D | 0.9353 | 314 | ★ top-5 |
| 4 | [`D016_IS621_ProtMPNN_T0.1_sample23`](D016_IS621_ProtMPNN_T0.1_sample23.md) | D | 0.9319 | 342 | ★ top-5 |
| 5 | [`D023_IS621_ProtMPNN_T0.1_sample8`](D023_IS621_ProtMPNN_T0.1_sample8.md) | D | 0.9319 | 342 | — |
| 6 | [`D022_IS621_ProtMPNN_T0.1_sample27`](D022_IS621_ProtMPNN_T0.1_sample27.md) | D | 0.9311 | 342 | — |
| 7 | [`D024_IS621_ProtMPNN_T0.1_sample15`](D024_IS621_ProtMPNN_T0.1_sample15.md) | D | 0.9311 | 342 | — |
| 8 | [`D030_IS621_ProtMPNN_T0.1_sample7`](D030_IS621_ProtMPNN_T0.1_sample7.md) | D | 0.9311 | 342 | — |
| 9 | [`D025_IS621_ProtMPNN_T0.1_sample19`](D025_IS621_ProtMPNN_T0.1_sample19.md) | D | 0.9303 | 342 | — |
| 10 | [`D7BKC8`](D7BKC8.md) | D | 0.9302 | 399 | — |
| 11 | [`D006_IS621_ProtMPNN_T0.1_sample10`](D006_IS621_ProtMPNN_T0.1_sample10.md) | D | 0.9294 | 342 | — |
| 12 | [`D008_IS621_ProtMPNN_T0.1_sample11`](D008_IS621_ProtMPNN_T0.1_sample11.md) | D | 0.9294 | 342 | — |
| 13 | [`D010_IS621_ProtMPNN_T0.1_sample17`](D010_IS621_ProtMPNN_T0.1_sample17.md) | D | 0.9294 | 342 | — |
| 14 | [`D011_IS621_ProtMPNN_T0.1_sample18`](D011_IS621_ProtMPNN_T0.1_sample18.md) | D | 0.9294 | 342 | — |
| 15 | [`D020_IS621_ProtMPNN_T0.1_sample6`](D020_IS621_ProtMPNN_T0.1_sample6.md) | D | 0.9294 | 342 | — |
| 16 | [`D026_IS621_ProtMPNN_T0.1_sample4`](D026_IS621_ProtMPNN_T0.1_sample4.md) | D | 0.9294 | 342 | — |

## Files in this directory

- `<design_id>.md` — per-design reference (summary, DNA, validation protocol)
- `wetlab_index.md` — this file

## Codon Optimisation Notes

Sequences were codon-optimised using the Kazusa Homo sapiens preferred-codon
table (most-frequent codon per amino acid, no global CAI optimisation).  
GC content typically 50–58%.  
Before commercial synthesis, verify:  
1. No restriction sites conflicting with your vector backbone.  
2. CAI ≥ 0.8 (calculate with CodonW or the IDT/Twist tool).  
3. No internal Kozak contexts (ATG in wrong frame).  
4. No cryptic splice sites (check with MaxEntScan).  

## Strategy-Specific Notes

**Strategy C (deimmunized IS621)**: These designs carry point mutations to
reduce HLA-II epitope load. The reference IS621 (`attP`/`attB` sites are
unchanged) should be used as the positive control in all recombination assays.

**Strategy D (ProtMPNN variants)**: These are IS621 sequence redesigns.
They share the same scaffold topology; expect similar or improved thermostability.
IS621 is the appropriate WT comparator for all biochemical assays.

**Strategy D natural orthologs (D8PEA4, D7BKC8)**: These are genuine IS110
orthologs from *Nitrospira defluvii* and *Arcanobacterium haemolyticum*.
Their cognate attB/attP sites are NOT the IS621 sites — site-specific
recombination assays require bioinformatic identification of flanking OATD
sites in their host genomes before wet-lab validation.
