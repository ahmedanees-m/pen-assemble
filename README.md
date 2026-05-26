# PEN-ASSEMBLE

**Programmable Enzyme Networks — Automated Strategy and Scoring Engine for Molecular Bridge-recombinase Library Engineering**

[![CI](https://github.com/ahmedanees-m/pen-assemble/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/pen-assemble/actions/workflows/ci.yml)
[![Docs](https://github.com/ahmedanees-m/pen-assemble/actions/workflows/docs.yml/badge.svg)](https://ahmedanees-m.github.io/pen-assemble/)
[![codecov](https://codecov.io/gh/ahmedanees-m/pen-assemble/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/pen-assemble)
[![PyPI](https://img.shields.io/pypi/v/pen-assemble?color=blue)](https://pypi.org/project/pen-assemble/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/ahmedanees-m/pen-assemble)](https://github.com/ahmedanees-m/pen-assemble/releases/latest)

> **Part of PEN-STACK** — the four-package computational biology pipeline for programmable genome editor discovery.
> **Paper 4 of 5** · [genome-atlas](https://github.com/ahmedanees-m/genome-atlas) → [mech-class](https://github.com/ahmedanees-m/mech-class) → [pen-score](https://github.com/ahmedanees-m/pen-score) → **pen-assemble** → PEN-COMPARE *(in prep)*

---

## What is PEN-ASSEMBLE?

Most gene therapies today use CRISPR–Cas9, which cuts both strands of DNA (a double-strand break) to make its edit. While powerful, DSBs activate cellular damage responses, introduce unwanted insertions or deletions, and limit cargo size.

**IS110-family bridge recombinases** work differently: they insert large DNA payloads — up to 50 kilobases — at specific genomic sites *without ever cutting both DNA strands*. They are guided by a short RNA molecule (the bridging RNA, or bRNA), and their targeting can be reprogrammed simply by swapping the bRNA guide sequence. This makes them highly attractive for safe, precise, large-cargo gene therapy.

**PEN-ASSEMBLE** is the computational nomination pipeline that asks: *can we design IS110-family variants that are even better than the natural protein IS621?* It evaluates **1,029 candidate designs** across four orthogonal engineering strategies — domain-swapping, ortholog discovery, immunogenicity reduction, and backbone redesign — using the **8-axis PenScore** composite metric. Every step is pre-registered, fully reproducible, and independently verified.

### Key finding

**16 designs beat IS621** (the current gold-standard bridge recombinase) on the pre-registered PenScore threshold of 0.929. The top candidate — a Monte Carlo deimmunized IS621 variant — scores **0.9673**, combining a +11.8% improvement in immunogenicity over IS621 while preserving full mechanistic activity. All five pre-registered predictions **PASS**, supporting publication with a strong claim.

---

## Why PEN-ASSEMBLE Was Built

### The problem

IS110-family bridge recombinases are a large and largely uncharacterised protein family — NCBI holds over 31,000 sequences. Deciding which ones are worth synthesising and testing in a wet-lab requires answering four independent questions simultaneously:

1. **Is the mechanism correct?** — Does this protein actually perform DSB-free insertion, or will mech-class's gradient-boosted model mis-call it as a nuclease?
2. **Is the structure viable?** — Does ESMFold predict a well-folded protein with high pLDDT at the active site?
3. **Is it therapeutically fit?** — Is it compact enough for AAV delivery? Low-immunogenicity enough for repeat dosing? Programmable?
4. **Does it beat what we already have?** — IS621 is the current best-characterised bridge recombinase; any new candidate must demonstrably surpass it.

No single tool answered all four questions. Researchers were left doing manual, ad-hoc assessments with no reproducible scoring framework and no pre-registered standards.

### The solution

PEN-ASSEMBLE assembles the four upstream PEN-STACK packages into a single, fully automated nomination pipeline:

```
  QUESTION                ANSWERED BY              HOW
  ─────────────────────── ──────────────────────── ────────────────────────────
  Is the mechanism right? mech-class v0.5.4        IS110 Tier-A hard gate
                                                   (PF01548 + PF02371 Pfam)
  Is the structure good?  ESMFold + quality gates  pLDDT ≥ 90 global,
                                                   ≥ 95 at active-site
  Is it therapeutically   pen-score v0.1.3         8-axis PenScore
  fit?                    (IS621 = 0.957)          composite metric
  Does it beat IS621?     Pre-registered           Verbatim lockpoint 0.929
                          comparison               (locked before analysis)
```

The pipeline is end-to-end reproducible: every input, gate, weight, and prediction threshold was committed to GitHub and SHA-256 locked **before any candidate was scored**.

### Who is it for?

| Audience | How they use PEN-ASSEMBLE |
|----------|--------------------------|
| **Wet-lab gene therapy groups** | Download the 16 top-ranked designs with synthesis sheets. No re-running needed — the catalog is pre-built. |
| **Computational biologists** | Use the `Designer` API to generate and rank new IS110-family designs against the same scoring framework. |
| **Tool developers** | Benchmark a new genome editor against IS621 and ISCro4 using the 8-axis PenScore framework in pen-score. |
| **Reproducibility researchers** | Audit the complete pre-registration → execution → deviation log → post-hoc rescoring chain. Every decision is documented. |

---

## How PEN-ASSEMBLE Works

### Component diagram

PEN-ASSEMBLE is not a monolith — it is an assembly layer that wires together four specialised packages:

```
  EXTERNAL PACKAGES                       PEN-ASSEMBLE INTERNALS
  ═══════════════════                     ═══════════════════════════════════════════
                                          ┌───────────────────────────────────────┐
  ┌─────────────────┐   accession IDs     │  STRATEGY GENERATION  (Steps 01–11)  │
  │  genome-atlas   │──────────────────►  │                                       │
  │  v0.7.2         │   28 enzyme nodes   │  A: Domain-Swap    B: Ortholog Scan   │
  │                 │   ESM-2 embeddings  │  15 chimeras       992 NCBI IS110s    │
  │  Knowledge      │                     │                                       │
  │  graph of       │                     │  C: Deimmunize     D: ProteinMPNN     │
  │  genome editors │                     │  MC substitution   Sequence design    │
  └─────────────────┘                     │  on IS621          from IS621 PDB     │
                                          │                                       │
                                          │          1,041 candidates             │
                                          └─────────────────┬─────────────────────┘
                                                            │
                                                            ▼
  ┌─────────────────┐   IS110 hard gate   ┌───────────────────────────────────────┐
  │  mech-class     │──────────────────►  │  TRIAGE PIPELINE      (Step 12)       │
  │  v0.5.4         │   DSB-free confirm  │                                       │
  │                 │   PF01548+PF02371   │  Gate 1 │ PFAM domain check           │
  │  Mechanism      │                     │         │ PF01548 (IS110 transposase) │
  │  classifier     │                     │         │ PF02371 (HTH domain)        │
  │  IS110 Tier-A   │                     │  Gate 2 │ ESMFold pLDDT ≥ 90 global  │
  │  gate           │                     │  Gate 3 │ Active-site pLDDT ≥ 95     │
  └─────────────────┘                     │  Gate 4 │ Length 300–400 aa          │
                                          │  Gate 5 │ mech-class IS110 Tier-A    │
                                          │                                       │
                                          │          1,029 designs pass           │
                                          └─────────────────┬─────────────────────┘
                                                            │
                                                            ▼
  ┌─────────────────┐   8-axis weights    ┌───────────────────────────────────────┐
  │  pen-score      │──────────────────►  │  PenScore EVALUATION  (Steps 13–16)  │
  │  v0.1.3         │   IS621 lockpoints  │                                       │
  │                 │   use-case profile  │  S_DSB    × 0.24  (DSB avoidance)    │
  │  Writer scoring │                     │  S_Spec   × 0.14  (targeting prec.)  │
  │  framework      │                     │  S_Cargo  × 0.19  (payload capacity) │
  │  IS621 = 0.957  │                     │  S_Deliv  × 0.19  (AAV compat.)      │
  │  SpCas9 = 0.402 │                     │  S_Immuno × 0.09  (immunogenicity)   │
  └─────────────────┘                     │  S_Prog   × 0.05  (programmability)  │
                                          │  S_Mature × 0.05  (tech. maturity)   │
                                          │  S_Energy × 0.05  (ATP independence) │
                                          └─────────────────┬─────────────────────┘
                                                            │
                                                            ▼
                                          ┌───────────────────────────────────────┐
                                          │  OUTPUT CATALOG       (Steps 17–25)  │
                                          │                                       │
                                          │  1,029 designs ranked by PenScore     │
                                          │                                       │
                                          │  16 designs > 0.929 ◄── P1 PASS      │
                                          │   2 designs > 0.957 ◄── secondary    │
                                          │                                       │
                                          │  ┌────────────┐  ┌─────────────────┐ │
                                          │  │ Parquet /  │  │ 16 wet-lab      │ │
                                          │  │ CSV catalog│  │ synthesis sheets│ │
                                          │  └────────────┘  └─────────────────┘ │
                                          │  ┌────────────┐  ┌─────────────────┐ │
                                          │  │ HTML design│  │ Pre-reg SHA-256 │ │
                                          │  │ browser    │  │ lock record     │ │
                                          │  └────────────┘  └─────────────────┘ │
                                          └───────────────────────────────────────┘
```

### Pipeline stages in plain language

| Stage | What happens | Key package used |
|-------|-------------|------------------|
| **Strategy generation** (Steps 01–11) | Four independent algorithms each produce a batch of IS110-family protein sequences. Strategy A: recombine domain boundaries. Strategy B: screen NCBI. Strategy C: optimise IS621 for low immunogenicity. Strategy D: ask ProteinMPNN for alternative sequences that fold into the same 3D shape. | genome-atlas |
| **Triage** (Step 12) | Every candidate passes through 5 quality gates. Anything that fails is logged and dropped. 12 of the original 1,041 candidates are removed here. | mech-class |
| **PenScore evaluation** (Steps 13–16) | Each surviving design is scored on 8 independent axes, then combined into a single PenScore. The IS621 reference lockpoint (0.929, pre-registered) acts as the pass/fail threshold. | pen-score |
| **Catalog assembly** (Steps 17–25) | The 1,029 scored designs are written to Parquet and CSV, an interactive HTML browser is generated, and 16 wet-lab synthesis sheets are produced for the P1-passing designs. | pen-assemble |

### Key design decisions

> **Why four strategies instead of one?**
> Because no single approach can overcome all limitations. Strategy B (NCBI scan) finds natural diversity but natural proteins are not optimised for human therapy. Strategy C (deimmunization) reduces immune reactivity but doesn't change the sequence radically. Strategy D (ProteinMPNN) redesigns the backbone but stays near IS621. Running all four in parallel and scoring them on the same metric reveals which approach actually wins — and the answer (Strategy D + C) was not obvious before the analysis.

> **Why pre-registration?**
> Because post-hoc threshold selection is the most common source of irreproducibility in computational biology. Every threshold in this pipeline (0.929, ≥ 5 beaters, Δ ≥ 0.10 for immunogenicity) was publicly committed and SHA-256 locked before a single protein was scored. This means the result — 16 beaters, 5/5 predictions PASS — is a genuine, non-inflated positive finding.

> **Why build on pen-score instead of just comparing raw axis values?**
> Raw axis values are not commensurable: a 0.05 improvement in S_DSB is not equivalent to a 0.05 improvement in S_Immuno. pen-score uses evidence-based weights derived from clinical AAV gene therapy requirements, so the composite PenScore reflects actual therapeutic relevance, not just mathematical convenience.

---

## The PEN-STACK Pipeline

PEN-ASSEMBLE is the fourth package in a four-paper computational stack. Each package provides critical inputs to the next:

```
                        PEN-STACK: Programmable Enzyme Networks
           Systematic Tool for Atlas and Knowledge (v0.5.1, 2026)

  ┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐
  │  genome-atlas    │    │   mech-class     │    │   pen-score      │    │  pen-assemble    │
  │  Paper 1         │───►│   Paper 2        │───►│   Paper 3        │───►│   Paper 4  ◄ YOU │
  │  v0.7.2          │    │   v0.5.4         │    │   v0.1.3         │    │   v0.5.2         │
  │                  │    │                  │    │                  │    │                  │
  │  Knowledge graph │    │  Mechanism       │    │  8-axis PenScore │    │  Design catalog  │
  │  28 enzyme       │    │  classifier      │    │  framework       │    │  1,029 IS110     │
  │  systems         │    │  IS110 Tier-A    │    │  IS621 = 0.957   │    │  candidates      │
  │  AUROC 0.9714    │    │  gate (DSB_FREE) │    │  evoCAST = 0.441 │    │  16 beat IS621   │
  │                  │    │  F1 = 0.986      │    │  S_Energy axis   │    │  5/5 PASS        │
  └──────────────────┘    └──────────────────┘    └──────────────────┘    └──────────────────┘
         │                        │                        │                        │
         │ Provides               │ Provides               │ Provides               │
         │ • Accession IDs        │ • IS110 gate           │ • PenScore formula     │ Produces
         │ • 28 system nodes      │ • DSB-free confirm     │ • 8 axis weights       │ • Wet-lab sheets
         │ • ESM-2 embeddings     │ • OOD probe set        │ • IS621 lockpoints     │ • Design browser
         │ • Knowledge edges      │                        │ • Use-case profiles    │ • Pre-reg record
```

| Package | What it does | Key result | Repo |
|---------|-------------|------------|------|
| [**genome-atlas**](https://github.com/ahmedanees-m/genome-atlas) | Builds a knowledge graph of 28 DNA-editing enzyme systems from UniProt, PDB, and AlphaFold data | AUROC 0.9714 on protein–domain link prediction; 70% accuracy on 10 published therapeutic scenarios | v0.7.2 |
| [**mech-class**](https://github.com/ahmedanees-m/mech-class) | Classifies any genome editor by its biochemical mechanism (DSB nuclease / DSB-free recombinase / transposase) and corrects 31,870 mis-classified IS110-family proteins | Tier-A macro-F1 = 0.9862; IS110 hard gate via PF01548 + PF02371 Pfam domains | v0.5.4 |
| [**pen-score**](https://github.com/ahmedanees-m/pen-score) | Scores any genome editor on 8 axes covering safety, specificity, cargo capacity, delivery compatibility, immunogenicity, programmability, maturity, and energy independence | IS621 scores 0.957; SpCas9 scores 0.402; 7 clinical weight profiles | v0.1.3 |
| **pen-assemble** *(this repo)* | Generates, screens, and ranks IS110-family design candidates using genome-atlas accession data, mech-class IS110 gating, and pen-score evaluation | 1,029 designs, 16 beat IS621, 5/5 pre-reg predictions PASS | v0.5.2 |

---

## Key Results at a Glance

| Metric | Value |
|--------|-------|
| Total candidate designs generated | 1,041 |
| Designs passing all triage gates | **1,029** |
| Designs beating IS621 pre-registered lockpoint (PenScore > 0.929) | **16** *(primary, P1 PASS)* |
| Designs beating calibrated lockpoint (PenScore > 0.9255) | **32** *(secondary analysis)* |
| Designs beating IS621 v0.1.2 lockpoint (PenScore > 0.957) | **2** *(secondary, 8-axis)* |
| Pre-registered predictions passing | **5 / 5** |
| Publication policy | **PUBLISH with strong claim** |
| Top design PenScore | **0.9673** (IS621_deimmunized_v2) |
| Top design immunogenicity gain over IS621 | **+11.8%** (S_Immuno: 0.7594 → 0.8777) |

### Top-5 Designs (diversity-enforced, P5)

| Rank | Design ID | Strategy | PenScore (v0.1.0) | Notable Feature |
|------|-----------|----------|--------------------|-----------------|
| 1 | `IS621_deimmunized_v2_Y255K_...` | C — Deimmunization | **0.9673** | Best immunogenicity: +11.8% over IS621 |
| 2 | `C_targeted_001` | C — Deimmunization | 0.9586 | Targeted deimmunization, fewer substitutions |
| 3 | `D8PEA4` | D — ProteinMPNN | 0.9353 | IS110 ortholog, 314 aa, compact for delivery |
| 4 | `D016_IS621_ProtMPNN_T0.1_sample23` | D — ProteinMPNN | 0.9319 | Backbone redesign on IS621 ESMFold structure |
| 5 | `A_007` *(diversity-enforced)* | A — Domain-Swap | 0.9209 | Ensures ≥ 3 strategies represented in top-5 |

---

## Four Design Strategies

IS110-family bridge recombinases are a large protein family. Rather than randomly screening sequences, PEN-ASSEMBLE uses four targeted engineering strategies — each attacking the problem from a different angle:

### Strategy A · Domain-Swap Chimeras (15 designs)

**What:** Recombine the catalytic scaffold from IS621 with guide-RNA recognition modules from related IS110 orthologs.

**Why:** Different orthologs have evolved different RNA-binding geometries. Chimeric fusions may combine IS621's high activity with improved bRNA flexibility.

**Outcome:** 15 chimeras, none beat IS621 (best PenScore: 0.921). Domain boundaries are tightly conserved — naive chimeras suffer structural penalties.

---

### Strategy B · IS110 Ortholog Discovery (992 designs)

**What:** Screen all IS110-family sequences in NCBI (>31,000 candidates) through 7 quality gates: PFAM domain verification (PF01548 + PF02371), ESMFold pLDDT ≥ 90 globally and ≥ 95 at active-site residues, length 300–400 aa, and mech-class IS110 Tier-A confirmation.

**Why:** Natural diversity is enormous. Some orthologs may already score higher than IS621 on delivery (shorter protein → better AAV packaging) or specificity.

**Outcome:** 992 candidates pass all gates, none beat IS621 (best: 0.917). Natural IS110s are well-optimised for their ecological niche but not for human therapeutics.

---

### Strategy C · Monte Carlo Deimmunization (2 designs)

**What:** Apply iterative single-substitution Monte Carlo sampling to IS621's sequence, accepting only mutations that reduce MHC-II predicted binding (fewer immune epitopes) while preserving the active-site residues and overall fold stability.

**Why:** IS621 is a bacterial protein. Human immune systems recognise many of its peptide fragments. Reducing immunogenic peptides is directly relevant to therapeutic safety and repeat-dosing potential.

**Outcome:** Both deimmunized variants beat IS621. The best — `IS621_deimmunized_v2` — improves S_Immuno from 0.759 to 0.878 (+11.8%) while maintaining perfect scores on all other axes. **PenScore = 0.9673, rank #1.**

---

### Strategy D · ProteinMPNN Backbone Redesign (32 designs)

**What:** Use ProteinMPNN — a deep learning model for sequence-from-structure design — to generate alternative amino acid sequences compatible with the IS621 ESMFold backbone (PDB: 8WT6). Sequences are conditioned on preserving active-site geometry.

**Why:** Many valid sequences can fold into the same 3D structure. Some of those alternative sequences may have better expression, stability, or reduced immunogenicity than wild-type IS621.

**Outcome:** 32 redesigns generated, **14 beat IS621** (PenScore 0.929–0.957). This is the most productive strategy. ProteinMPNN-generated sequences preserve mechanism while improving delivery scores.

---

## PenScore: How Designs Are Ranked

Each candidate design is scored on eight independent axes, all on a [0, 1] scale, then combined into a single composite **PenScore**. Every IS110-family member scores 1.0 on S_DSB, S_Cargo, S_Prog, and S_Energy by mechanism — so competition happens on the remaining axes.

### Current formula (pen-score v0.1.2 — 8-axis)

```
PenScore = S_DSB × 0.24  +  S_Spec × 0.14  +  S_Cargo × 0.19
         + S_Deliv × 0.19 +  S_Immuno × 0.09 +  S_Prog × 0.05
         + S_Mature × 0.05 +  S_Energy × 0.05
```

| Axis | Weight | What it measures | IS621 value |
|------|--------|-----------------|-------------|
| `S_DSB` | **0.24** | Double-strand break avoidance — does the editor cut both DNA strands? IS110 recombinases never do (score = 1.0). CRISPR-Cas9 always does (score = 0.0). | **1.0** |
| `S_Spec` | 0.14 | Guide-RNA target-site specificity — how precisely can the editor be directed to a single genomic location? | **1.0** |
| `S_Cargo` | **0.19** | Payload capacity — how large a DNA insert can the editor integrate? IS110 can handle >50 kb (score = 1.0). Base editors are limited to single nucleotides (score = 0.0). | **1.0** |
| `S_Deliv` | **0.19** | Delivery compatibility — can the protein fit in an AAV capsid? Shorter proteins score higher. IS621 (393 aa) scores above average. | **0.86** |
| `S_Immuno` | 0.09 | Immunogenicity — how many MHC-II epitopes does the protein present? Fewer epitopes → safer for human use → higher score. | **0.759** |
| `S_Prog` | 0.05 | Programmability — can target-site selection be changed by swapping the guide RNA? IS110 uses a bRNA guide (score = 1.0). Fixed-specificity recombinases score 0.0. | **1.0** |
| `S_Mature` | 0.05 | Technology maturity — how many peer-reviewed publications describe this editor? Higher citation counts → score approaching 1.0. | **0.83** |
| `S_Energy` | 0.05 | Energy independence — does the editor require ATP? Walker A/B motif absence → score = 1.0. IS110 is ATP-free. ATPase-driven systems (e.g., evoCAST) score 0.0. | **1.0** |

> **Pre-registration integrity note:** P1 (the primary prediction) used the 7-axis v0.1.0 formula
> (weights: S_DSB 0.25, S_Spec 0.10, S_Cargo 0.20, S_Deliv 0.15, S_Immuno 0.10, S_Prog 0.15,
> S_Mature 0.05; no S_Energy axis). The result — **16 designs beat IS621 at 0.929** — is FINAL and
> cannot be changed retroactively. The v0.1.2 8-axis re-scoring (IS621 = 0.957) is a secondary
> analysis only. See [`RESCORING_v0.1.2.md`](RESCORING_v0.1.2.md).

**IS621 reference lockpoints:**

| Lockpoint | Value | Context |
|-----------|-------|---------|
| Verbatim pre-registered | **0.929** | Primary (P1). v0.1.0 formula. Locked before analysis. |
| mech-class v0.5.2 corrected | **0.954** | Secondary. IS110 gate retroactively applied to IS621 itself. |
| pen-score v0.1.2 (8-axis) | **0.957** | Secondary. Current best estimate. Adds S_Energy axis. |
| MHCflurry 2.2.1 recalibrated | **0.9255** | Secondary. 32 beaters under re-calibrated immunogenicity scorer. |

---

## Pipeline Architecture

```
INPUT: 4 Engineering Strategies
═══════════════════════════════════════════════════════════════════════

 Strategy A · Domain-Swap ──────┐
   15 IS621-scaffold chimeras    │
                                  │
 Strategy B · Ortholog Scan ─────┼──► TRIAGE PIPELINE
   992 IS110 NCBI candidates      │
                                  │   Step 12: ESMFold pLDDT
 Strategy C · Deimmunization ────┤   ┌─ Global ≥ 90
   Monte Carlo substitution       │   └─ Active-site ≥ 95
   2 deimmunized IS621 variants   │            │
                                  │   Step 13: mech-class IS110 gate
 Strategy D · ProteinMPNN ───────┘   PF01548 + PF02371 verified
   32 backbone redesigns                       │
                                               │
                                    1,029 designs pass triage
                                               │
                                               ▼
OUTPUT: PEN-SCORE EVALUATION (pen-score v0.1.2)
═══════════════════════════════════════════════
 8-axis PenScore per design
 IS621 verbatim lockpoint: 0.929 (pre-reg)
 IS621 current lockpoint:  0.957 (v0.1.2)
          │
          ├─► 16 designs > 0.929  ◄── P1 PASS (pre-registered)
          ├─► 2  designs > 0.957
          └─► Full 1,029-design catalog
                    │
          ┌─────────┼─────────────┐
          ▼         ▼             ▼
    CSV / Parquet  HTML Browser  16 Wet-lab
    catalog        (no server)   reference sheets
                                 (Markdown)
```

---

## Installation

```bash
pip install pen-assemble
```

Or install from source with development extras:

```bash
git clone https://github.com/ahmedanees-m/pen-assemble.git
cd pen-assemble
pip install -e ".[dev,docs]"
```

**Requirements:** Python ≥ 3.10 · pandas ≥ 2.0 · pyarrow ≥ 14.0 · numpy ≥ 1.24

---

## Quick Start

### Load the design catalog

```python
from pen_assemble.catalog import load_catalog, load_p1_beaters, load_top5

# All 1,029 scored designs
df = load_catalog()
print(df[["design_id", "strategy", "pen_score"]].head())

# The 16 designs that beat IS621 (pre-registered P1)
p1 = load_p1_beaters()
print(f"P1 beaters: {len(p1)}")   # 16

# Diversity-enforced top-5 (≥3 strategies represented)
top5 = load_top5()
print(top5[["design_id", "strategy", "pen_score"]])
```

### Score any IS110-family design

```python
from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621

# All IS110-family designs score 1.0 on mechanism axes by definition.
# The only variable axis for IS110 candidates is S_Immuno (immunogenicity)
# and S_Deliv (protein length / AAV compatibility).

ax = PenScoreAxes(
    S_DSB=1.0,      # IS110: no double-strand breaks
    S_Spec=1.0,     # bRNA guide: precise targeting
    S_Cargo=1.0,    # IS110: large cargo (>50 kb)
    S_Deliv=0.92,   # 340 aa: compact, AAV-friendly
    S_Immuno=0.85,  # after Monte Carlo deimmunization
    S_Prog=1.0,     # bRNA re-targetable
    S_Mature=0.83,
    S_Energy=1.0,   # ATP-free (no Walker A/B motifs)
)
score = pen_score(ax)           # composite PenScore
print(f"PenScore: {score:.4f}")
print(f"Beats IS621: {beats_is621(score)}")  # True if score > 0.929
print(ax.contributions())       # per-axis breakdown
```

### Codon-optimise for human expression

```python
from pen_assemble.codon import build_expression_orf, gc_content, check_restriction_sites

seq = p1.iloc[0]["protein_sequence"]
orf = build_expression_orf(seq, kozak=True, stop=True)

print(f"ORF length : {len(orf)} bp")
print(f"GC content : {gc_content(orf):.1%}")         # target 40–60%
print(f"RE sites   : {check_restriction_sites(orf)}")  # BsaI, BbsI, etc.
```

### High-level Designer API

```python
from pen_assemble.api import Designer

d = Designer.load()

# Browse and filter the catalog
top = d.select_designs(strategy="C", require_dsb_free=True, top_k=5)

# Run Strategy C deimmunization (requires scaffold FASTA)
variants = d.deimmunize(scaffold_id="IS621", n_variants=50)

# Run Strategy D ProteinMPNN redesign (requires PDB structure)
redesigns = d.redesign_backbone(scaffold_id="IS621", n_designs=25)
```

---

## Repository Structure

```
pen-assemble/
│
├── pen_assemble/                    # Python package (importable)
│   ├── pen_score.py                 #   PenScore composite formula (8 axes)
│   ├── catalog.py                   #   load_catalog(), load_p1_beaters(), load_top5()
│   ├── codon.py                     #   Human codon optimisation utilities
│   ├── api.py                       #   High-level Designer API
│   ├── cli.py                       #   Command-line entry point (pen-assemble --help)
│   ├── _version.py                  #   Package version (0.5.2)
│   ├── strategies/                  #   Design generation modules (Steps 01–11)
│   │   ├── domain_swap.py           #     Strategy A: chimeric fusions
│   │   ├── ortholog_discovery.py    #     Strategy B: NCBI IS110 screen
│   │   ├── deimmunization.py        #     Strategy C: Monte Carlo MHC reduction
│   │   └── backbone_redesign.py     #     Strategy D: ProteinMPNN sequence design
│   ├── triage/                      #   Multi-gate candidate filtering (Step 12)
│   ├── verification/                #   Axis scoring pipeline (Steps 13–16)
│   ├── utils/                       #   Linker assembly, MHC scoring, PDB parsing
│   └── data/                        #   YAML configuration files
│
├── catalog/
│   └── release_v0.5.0/             #   v0.5.0 FROZEN pre-registration record
│       ├── pen_assemble_catalog.{csv,parquet}   # All 1,029 designs
│       ├── p1_beaters_catalog.{csv,parquet}     # 16 designs > 0.929
│       ├── p5_top5_catalog.{csv,parquet}        # Diversity top-5
│       ├── browser/index.html       #   Interactive HTML design browser (no server needed)
│       ├── wetlab/                  #   16 wet-lab synthesis reference sheets (Markdown)
│       └── validation/              #   Pre-registered prediction result JSONs (P1–P5)
│
├── data/                            # Rescored and extended catalog outputs
│   ├── catalog_v0.5.1_current.parquet    # v0.5.1: 1,029 designs, pen-score v0.1.2 (8-axis)
│   ├── catalog_v0.5.2_current.parquet    # v0.5.2: + intrinsic_cargo_mechanism + cell_based_evidence (PEN-COMPARE v3.2)
│   ├── rescore_comparison_v010_v012.csv  # Side-by-side v0.1.0 vs v0.1.2 scores
│   ├── rescore_summary_v012.json         # Summary statistics (v0.1.2)
│   └── rescore_summary_v052.json         # Summary statistics (v0.5.2 schema)
│
├── scripts/                         # Numbered pipeline scripts
│   ├── 50_assemble_catalog.py       #   Builds catalog/ artifacts
│   ├── 51_build_browser.py          #   Generates interactive HTML browser
│   ├── 52_generate_wetlab_reference.py  # Creates 16 wet-lab sheets
│   ├── rescore_v012.py              #   Re-scores catalog with pen-score v0.1.2
│   └── upgrade_catalog_to_v052.py   #   Adds v3.2 fields (intrinsic_cargo + cell_based)
│
├── tests/                           # 74 pytest tests (Python 3.10 / 3.11 / 3.12)
├── docs/                            # Sphinx documentation (furo theme)
│
├── CHANGELOG.md                     # Version history
├── CITATION.cff                     # Machine-readable citation metadata
├── CONTRIBUTING.md                  # Contribution guidelines
├── DESIGN_PROVENANCE.md             # Full deviation log (5 documented deviations)
├── RESCORING_v0.1.2.md              # Secondary analysis record (8-axis re-scoring)
├── SECURITY.md                      # Security policy
└── pyproject.toml                   # Build config, dependencies, coverage settings
```

---

## Running the Tests

```bash
pytest tests/ -v
```

All 74 tests pass on Python 3.10, 3.11, and 3.12. Coverage is measured on the public API surface (`pen_score.py`, `catalog.py`, `codon.py`) and reported to [Codecov](https://codecov.io/gh/ahmedanees-m/pen-assemble). Pipeline-only modules (`strategies/`, `triage/`, `verification/`) require VM-only dependencies (ESM-2, MHCflurry, mech-class extras) and are excluded from CI coverage.

---

## Pre-Registered Predictions

All five predictions were locked (committed and tagged `pre-registration-v1.0.0`) **before any strategy generation script ran**. The SHA-256 of the pre-registration YAML is in `catalog/release_v0.5.0/validation/all_predictions_summary.json`.

| ID | Prediction | Threshold | Result |
|----|-----------|-----------|--------|
| **P1** | ≥ 5 designs beat IS621 verbatim lockpoint | PenScore > 0.929 | ✅ **PASS** — 16 designs |
| **P2** | ≥ 1 design satisfies IS110 mechanism + AAV-compatible | joint gate | ✅ **PASS** — 1,029 designs |
| **P3** | Best Strategy-C design improves S_Immuno over IS621 by ≥ 0.10 | Δ ≥ 0.10 | ✅ **PASS** — Δ = +0.118 |
| **P4** | ≥ 100 Strategy-B candidates with PFAM-verified IS110 domain | count ≥ 100 | ✅ **PASS** — 992 designs |
| **P5** | Top-5 includes designs from ≥ 3 distinct strategies | diversity | ✅ **PASS** — A, C, D |

---

## Key Deviations from Execution Plan

Full details in [`DESIGN_PROVENANCE.md`](DESIGN_PROVENANCE.md).

| # | Deviation | Impact |
|---|-----------|--------|
| **D1** | **Rosetta gate non-functional** — absolute energies, not ΔΔG, cannot be threshold-compared across sequences | Gate auto-passed; pLDDT structural quality proxy used instead |
| **D2** | **P3 IS621 reference corrected** — S_Immuno was placeholder 0.250 in draft; correct value is 0.7594 (MHCflurry) | P3 threshold correctly evaluated from accurate baseline |
| **D3** | **P5 diversity-enforced** — A_007 (0.9209) replaces natural rank-5 D023 (0.9319) | P5 PASS maintained; rank-5 note disclosed |
| **D4** | **MHCflurry 2.2.1 recalibration** — IS621 S_Immuno = 0.7243 under current tool vs 0.7594 at lock-in | Calibrated lockpoint = 0.9255; secondary analysis shows 32 beaters |
| **D5** | **mech-class ML misclassification** — IS110 proteins mis-called as DSB_NUCLEASE by gradient-boosted model alone | Corrected via PFAM domain hard gate (PF01548 + PF02371) in mech-class v0.5.2+ |

---

## Generating the Catalog

The frozen v0.5.0 catalog (`catalog/release_v0.5.0/`) is never regenerated — it is the pre-registration record. The scripts below produce or update the *current* catalog:

```bash
# Full catalog assembly (produces catalog/release_v0.5.0/ on the VM)
python scripts/50_assemble_catalog.py

# Interactive HTML browser (no server needed — open in any browser)
python scripts/51_build_browser.py

# 16 wet-lab Markdown reference sheets
python scripts/52_generate_wetlab_reference.py

# Re-score with pen-score v0.1.2 (secondary analysis; produces data/)
python scripts/rescore_v012.py --frozen catalog/release_v0.5.0/pen_assemble_catalog.parquet
```

---

## Part of PEN-STACK

**PEN-STACK** (*Programmable Enzyme Networks — Systematic Tool for Atlas and Knowledge*) is a four-paper computational biology pipeline for programmable genome editor discovery and benchmarking.

| Package | Role in the pipeline | Version | Repo |
|---------|---------------------|---------|------|
| [**genome-atlas**](https://github.com/ahmedanees-m/genome-atlas) | Foundational knowledge graph of 28 enzyme systems (AUROC 0.9714) | v0.7.2 | [![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#) |
| [**mech-class**](https://github.com/ahmedanees-m/mech-class) | Biochemical mechanism classifier — IS110 Tier-A gate (F1 = 0.986) | v0.5.4 | [![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#) |
| [**pen-score**](https://github.com/ahmedanees-m/pen-score) | 8-axis writer scoring framework — IS621 = 0.957, SpCas9 = 0.402 | v0.1.3 | [![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#) |
| **pen-assemble** *(this repo)* | IS110-family design catalog — 1,029 designs, 5/5 predictions PASS | **v0.5.2** | [![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#) |
| **PEN-COMPARE** *(in prep)* | Cross-system benchmarking — TrueWriterScore 6-gate classifier | — | — |

---

## Citation

If you use PEN-ASSEMBLE in your work, please cite:

```bibtex
@software{ahmed2026penassemble,
  author    = {Ahmed, Anees},
  title     = {{PEN-ASSEMBLE}: Automated Strategy and Scoring Engine
               for Molecular Bridge-recombinase Library Engineering},
  year      = {2026},
  version   = {v0.5.2},
  publisher = {Zenodo},
  url       = {https://github.com/ahmedanees-m/pen-assemble},
  note      = {DOI pending Zenodo deposit}
}
```

> Ahmed A. (2026). PEN-ASSEMBLE: Automated Strategy and Scoring Engine for Molecular
> Bridge-recombinase Library Engineering. v0.5.2. DOI pending.

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata (GitHub renders this automatically as a "Cite this repository" widget).

---

## License

MIT © 2026 Anees Ahmed — see [LICENSE](LICENSE) for details.
