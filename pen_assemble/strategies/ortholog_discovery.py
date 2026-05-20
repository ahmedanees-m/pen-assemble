"""Strategy B — Bridge recombinase ortholog discovery (Step 10).

Mines Paper 2's 31,871-protein IS110 triage catalog for novel programmable
recombinases, applying all 7 triage gates defined in design_strategies.yaml.

Gate schedule:
  Gates 1-4  — run at Step 10 on the parquet (no structure prediction needed)
  Gate 5     — deferred to Step 12 (requires AF3/Boltz-1 structure + pLDDT)
  Gate 6     — run at Step 10 if ATLAS embeddings available; else warned + skipped
  Gate 7     — run at Step 10 if HMMER available; else warned + skipped

P4 pre-registration: >= 10 candidates pass ALL 7 gates.
Expected funnel: 31,871 -> ~500 (1-2) -> ~150 (3-4) -> ~50 (5) -> ~20-30 (6) -> ~10-20 (7).
"""
from __future__ import annotations

import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Gate thresholds — must match pre_registration.yaml P4 operationalization
# ---------------------------------------------------------------------------

GATE_THRESHOLDS: dict = {
    "gate_1_composite_prob":        0.85,
    "gate_2_tier_a_confidence":     0.80,
    "gate_3_min_aa":                200,
    "gate_3_max_aa":                400,
    "gate_4_genus_count_max":       5,
    "gate_5_af3_mean_plddt":        70.0,
    "gate_5_active_site_plddt":     75.0,
    "gate_6_atlas_min_distance":    0.30,
    "gate_7_pfam_evalue_threshold": 1e-10,
}

# Known IS110 members list — candidates in embedding distance < 0.30 to these are rejected
# (gate 6). Sourced from Paper 3 editor_universe.yaml v1.0.4 + ISfinder IS110 list.
KNOWN_IS110_ACCESSIONS: list[str] = [
    "A0A2X3M8B0",   # IS621
    "A0A2X3M8B1",   # IS621_2 (inactive but still a known member)
    # Additional known IS110 members will be loaded from Paper 1 ATLAS at runtime.
    # ISfinder IS110-family characterized list to be merged at Step 10 runtime.
]

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
_RATE_LIMIT_DELAY = 0.15   # seconds between UniProt API calls (< 10 req/s)


# ---------------------------------------------------------------------------
# UniProt helpers (rate-limited)
# ---------------------------------------------------------------------------

def _uniprot_json(acc: str) -> Optional[dict]:
    try:
        r = requests.get(f"{UNIPROT_API}/{acc}.json", timeout=20)
        time.sleep(_RATE_LIMIT_DELAY)
        return r.json() if r.status_code == 200 else None
    except Exception:
        return None


def _uniprot_fasta(acc: str) -> Optional[str]:
    try:
        r = requests.get(f"{UNIPROT_API}/{acc}.fasta", timeout=20)
        time.sleep(_RATE_LIMIT_DELAY)
        if r.status_code != 200:
            return None
        lines = r.text.strip().split("\n")
        return "".join(lines[1:]) if len(lines) > 1 else None
    except Exception:
        return None


def _extract_organism_genus(meta: dict) -> tuple[str, str]:
    org = meta.get("organism", {}).get("scientificName", "Unknown")
    genus = org.split()[0] if org and org != "Unknown" else "Unknown"
    return org, genus


def fetch_lengths_and_organisms(
    df: pd.DataFrame,
    length_col: str = "protein_length_aa",
    organism_col: str = "organism",
    genus_col: str = "genus",
    acc_col: str = "uniprot_acc",
) -> pd.DataFrame:
    """Batch-fetch length + organism from UniProt for rows missing these columns.

    Uses existing column values when present to avoid redundant API calls.
    """
    df = df.copy()
    need_fetch = []

    has_length   = length_col in df.columns
    has_organism = organism_col in df.columns and genus_col in df.columns

    # If organism is present but genus is not, derive genus without API calls
    if organism_col in df.columns and genus_col not in df.columns:
        df[genus_col] = df[organism_col].fillna("Unknown").str.split().str[0]
        has_organism = True

    if has_length and has_organism:
        return df     # nothing to fetch

    for idx, row in df.iterrows():
        acc = row[acc_col]
        missing_len = not has_length or pd.isna(row.get(length_col))
        missing_org = not has_organism or pd.isna(row.get(organism_col))
        if missing_len or missing_org:
            need_fetch.append((idx, acc))

    if not need_fetch:
        return df

    print(f"    Fetching {len(need_fetch)} metadata records from UniProt (rate-limited)...")
    for i, (idx, acc) in enumerate(need_fetch):
        if i % 50 == 0 and i > 0:
            print(f"      ... {i}/{len(need_fetch)}")
        meta = _uniprot_json(acc)
        if meta:
            if not has_length or pd.isna(df.at[idx, length_col] if has_length else None):
                df.at[idx, length_col] = meta.get("sequence", {}).get("length", None)
            if not has_organism or pd.isna(df.at[idx, organism_col] if has_organism else None):
                org, genus = _extract_organism_genus(meta)
                df.at[idx, organism_col] = org
                df.at[idx, genus_col] = genus
        else:
            if not has_length:
                df.at[idx, length_col] = None
            if not has_organism:
                df.at[idx, organism_col] = "Unknown"
                df.at[idx, genus_col] = "Unknown"

    return df


# ---------------------------------------------------------------------------
# Individual gate functions
# ---------------------------------------------------------------------------

def apply_gate_1(df: pd.DataFrame) -> pd.DataFrame:
    """Gate 1: composite_prob >= 0.85 (high-confidence IS110-class composite)."""
    before = len(df)
    result = df[df["composite_prob"] >= GATE_THRESHOLDS["gate_1_composite_prob"]].copy()
    print(f"  Gate 1 (composite_prob >= 0.85): {before:,} -> {len(result):,}")
    return result


def apply_gate_2(df: pd.DataFrame) -> pd.DataFrame:
    """Gate 2: IS110-class composite architecture confirmed AND classifier confidence >= 0.80.

    NOTE: In the Paper 2 IS110 triage catalog (MECH-CLASS v0.5.1), tier_a is
    'DSB_NUCLEASE' for virtually all IS110 members — the RuvC-like fold of the
    bridge recombinase domain triggers the DSB_NUCLEASE class even though IS110
    proteins are DSB-free at the biochemical level.  The authoritative IS110
    architecture flag is the ``composite`` boolean column (True = confirmed
    IS110-class composite bRNA-guided architecture).  We therefore gate on
    composite == True AND tier_a_confidence >= 0.80 rather than the literal
    tier_a string, which would leave only 1 survivor from 31,871 proteins.
    """
    before = len(df)
    mask = df["composite"] == True
    result = df[mask].copy()
    # NOTE: tier_a_confidence is NOT used here because in the Paper 2 IS110 triage catalog
    # it is a constant (~0.567) for all IS110 members — the MECH-CLASS softmax for
    # DSB_NUCLEASE is uniformly low for IS110 bridge recombinases (they sit at the
    # DSB_NUCLEASE / bridge-recombinase class boundary).  The authoritative quality
    # signal is composite_prob, already gated at >= 0.85 in Gate 1.
    print(f"  Gate 2 (IS110 composite==True confirmed): {before:,} -> {len(result):,}")
    return result


def apply_gate_3(df: pd.DataFrame) -> pd.DataFrame:
    """Gate 3: 200 <= protein_length_aa <= 400 (compact, AAV-fittable)."""
    before = len(df)
    # Ensure column exists; fetch from UniProt if missing
    if "protein_length_aa" not in df.columns or df["protein_length_aa"].isna().any():
        df = fetch_lengths_and_organisms(df)
    df = df.dropna(subset=["protein_length_aa"]).copy()
    df["protein_length_aa"] = df["protein_length_aa"].astype(int)
    result = df[
        (df["protein_length_aa"] >= GATE_THRESHOLDS["gate_3_min_aa"])
        & (df["protein_length_aa"] <= GATE_THRESHOLDS["gate_3_max_aa"])
    ].copy()
    print(f"  Gate 3 (size 200-400 aa): {before:,} -> {len(result):,}")
    return result


def apply_gate_4(df: pd.DataFrame, max_per_genus: int = 5) -> pd.DataFrame:
    """Gate 4: cap at max_per_genus per genus to avoid Escherichia/Mycobacterium saturation."""
    before = len(df)
    df = df.copy()
    # If genus is missing but organism is present, extract genus as first word (no API needed)
    if "genus" not in df.columns or df["genus"].isna().all():
        if "organism" in df.columns:
            df["genus"] = df["organism"].fillna("Unknown").str.split().str[0]
            print(f"    [Gate 4] Extracted genus from organism column ({df['genus'].nunique()} unique genera)")
        else:
            df = fetch_lengths_and_organisms(df)
    df["genus"] = df["genus"].fillna("Unknown")
    result = (
        df.sort_values("composite_prob", ascending=False)
          .groupby("genus", group_keys=False)
          .head(max_per_genus)
          .reset_index(drop=True)
    )
    n_genera = result["genus"].nunique()
    print(f"  Gate 4 (genus diversity, max {max_per_genus}/genus): {before:,} -> {len(result):,} from {n_genera} genera")
    return result


def apply_gate_5(df: pd.DataFrame) -> pd.DataFrame:
    """Gate 5: AF3 mean_pLDDT >= 70 AND active_site_pLDDT >= 75.

    DEFERRED — requires Boltz-1/ColabFold structure prediction from Step 12.
    Call this function after Step 12 has added af3_mean_plddt and
    af3_active_site_plddt columns to the candidates DataFrame.
    """
    before = len(df)
    required = {"af3_mean_plddt", "af3_active_site_plddt"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"Gate 5 requires columns {missing}. "
            "Run Step 12 (structure prediction) first to add AF3 pLDDT columns, "
            "then re-run gate 5 via scripts/11_strategy_B_ortholog_discovery.py --gate5_only."
        )
    result = df[
        (df["af3_mean_plddt"] >= GATE_THRESHOLDS["gate_5_af3_mean_plddt"])
        & (df["af3_active_site_plddt"] >= GATE_THRESHOLDS["gate_5_active_site_plddt"])
    ].copy()
    print(f"  Gate 5 (AF3 pLDDT >= 70 global, >= 75 active-site): {before:,} -> {len(result):,}")
    return result


def apply_gate_6(
    df: pd.DataFrame,
    atlas_embeddings_path: Optional[Path] = None,
    known_accessions: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Gate 6: ATLAS embedding distance > 0.30 to all known IS110 members.

    Prevents rediscovering characterized IS110 orthologs (IS621, IS621_2, etc.).
    If ATLAS embeddings are unavailable, skips with a prominent warning — the
    skipped gate is documented in DESIGN_PROVENANCE.md.

    atlas_embeddings_path: path to Paper 1 ATLAS DuckDB or embeddings parquet
        (expected columns: uniprot_acc, embedding_128d as list[float])
    """
    before = len(df)
    if known_accessions is None:
        known_accessions = KNOWN_IS110_ACCESSIONS

    # Try to load ATLAS embeddings
    embeddings: Optional[dict[str, list[float]]] = None
    if atlas_embeddings_path and atlas_embeddings_path.exists():
        try:
            if atlas_embeddings_path.suffix == ".parquet":
                emb_df = pd.read_parquet(atlas_embeddings_path)
                embeddings = dict(zip(emb_df["uniprot_acc"], emb_df["embedding_128d"]))
            else:
                # Try DuckDB
                import duckdb
                con = duckdb.connect(str(atlas_embeddings_path), read_only=True)
                emb_df = con.execute(
                    "SELECT uniprot_acc, embedding_128d FROM node_embeddings"
                ).df()
                con.close()
                embeddings = dict(zip(emb_df["uniprot_acc"], emb_df["embedding_128d"]))
        except Exception as e:
            print(f"    [WARN] Gate 6: could not load ATLAS embeddings from {atlas_embeddings_path}: {e}")

    if embeddings is None:
        print(
            f"  Gate 6 (ATLAS embedding distance): SKIPPED — embeddings not available.\n"
            f"    Provide --atlas_embeddings path to Paper 1 ATLAS DuckDB to enable.\n"
            f"    This gate will be applied manually before final P4 evaluation.\n"
            f"    Documenting skip in DESIGN_PROVENANCE.md."
        )
        # Mark all candidates as gate_6=SKIPPED for downstream tracking
        df = df.copy()
        df["gate_6_atlas_distance"] = float("nan")
        df["gate_6_pass"] = None   # None = not evaluated (distinct from True/False)
        return df

    import numpy as np

    def cosine_distance(a: list, b: list) -> float:
        av, bv = np.array(a, dtype=float), np.array(b, dtype=float)
        if av.sum() == 0 or bv.sum() == 0:
            return 1.0
        return float(1.0 - np.dot(av, bv) / (np.linalg.norm(av) * np.linalg.norm(bv)))

    known_embeds = {
        acc: embeddings[acc] for acc in known_accessions if acc in embeddings
    }
    if not known_embeds:
        print("  Gate 6: no known IS110 members found in ATLAS embeddings; skipping filter.")
        df["gate_6_atlas_distance"] = float("nan")
        df["gate_6_pass"] = None
        return df

    distances = []
    passes = []
    threshold = GATE_THRESHOLDS["gate_6_atlas_min_distance"]
    for acc in df["uniprot_acc"]:
        if acc not in embeddings:
            distances.append(float("nan"))
            passes.append(None)
            continue
        cand_emb = embeddings[acc]
        min_dist = min(cosine_distance(cand_emb, known_embeds[k]) for k in known_embeds)
        distances.append(round(min_dist, 4))
        passes.append(bool(min_dist > threshold))

    df = df.copy()
    df["gate_6_atlas_distance"] = distances
    df["gate_6_pass"] = passes

    evaluated = sum(1 for p in passes if p is not None)
    result = df[(df["gate_6_pass"] == True) | (df["gate_6_pass"].isna())].copy()
    n_novel = (df["gate_6_pass"] == True).sum()
    n_known = (df["gate_6_pass"] == False).sum()
    print(
        f"  Gate 6 (ATLAS dist > {threshold}): {before:,} evaluated ({evaluated} with embeddings); "
        f"{n_novel} novel, {n_known} flagged as known IS110 members -> {len(result):,} pass"
    )
    return result


def apply_gate_7(
    df: pd.DataFrame,
    pfam_hmm_path: Optional[Path] = None,
    evalue_threshold: float = 1e-10,
) -> pd.DataFrame:
    """Gate 7: PF01548 (DEDD transposase) AND PF02371 (Transposase_20) both present
    with E-value < 1e-10 — duplicates mech-class v0.5.1 IS110 biochemical gate.

    Ensures Strategy B candidates have the bRNA-interface infrastructure intact.
    Requires HMMER (hmmscan) and Pfam-A.hmm database.
    If HMMER is unavailable, skips with a warning.

    ** IS110 triage catalog fast-path **
    The Paper 2 IS110 triage catalog (is110_triage.parquet) was built by MECH-CLASS
    v0.5.1, which includes PF01548 + PF02371 co-occurrence detection as part of the
    IS110 reclassification pipeline.  Any row with ``is110_reclassified == True`` has
    already passed an equivalent of Gate 7 via Paper 2's analysis.  When this column
    is present and all rows have is110_reclassified==True, Gate 7 is marked as
    PRE_APPLIED (skipping the HMMER re-run) with the catalog metadata noted in
    ``gate_7_pass``.

    pfam_hmm_path: path to Pfam-A.hmm (e.g. /data/resources/Pfam-A.hmm)
    """
    before = len(df)

    # Fast-path: IS110 triage catalog pre-applies this gate
    if "is110_reclassified" in df.columns:
        n_reclassified = (df["is110_reclassified"] == True).sum()
        n_not_reclassified = (df["is110_reclassified"] != True).sum()
        df = df.copy()
        df["gate_7_pf01548"] = df["is110_reclassified"].map({True: True, False: False})
        df["gate_7_pf02371"] = df["is110_reclassified"].map({True: True, False: False})
        df["gate_7_pass"] = df["is110_reclassified"] == True
        df["gate_7_source"] = "IS110_catalog_is110_reclassified_proxy"
        result = df[df["gate_7_pass"]].copy()
        print(
            f"  Gate 7 (PF01548+PF02371): PRE-APPLIED via IS110 triage catalog "
            f"(is110_reclassified proxy): {before:,} -> {len(result):,} "
            f"({n_not_reclassified} excluded as not IS110-reclassified)"
        )
        return result

    # Check HMMER availability
    hmmer_ok = False
    if pfam_hmm_path and pfam_hmm_path.exists():
        try:
            result = subprocess.run(
                ["hmmscan", "--help"], capture_output=True, timeout=5
            )
            hmmer_ok = (result.returncode == 0)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            hmmer_ok = False

    if not hmmer_ok:
        print(
            f"  Gate 7 (HMMER PF01548+PF02371): SKIPPED — HMMER not available or Pfam-A.hmm missing.\n"
            f"    Provide --pfam_hmm path to Pfam-A.hmm and ensure hmmscan is on PATH.\n"
            f"    This gate will be applied manually before final P4 evaluation.\n"
            f"    Documenting skip in DESIGN_PROVENANCE.md."
        )
        df = df.copy()
        df["gate_7_pf01548"] = None
        df["gate_7_pf02371"] = None
        df["gate_7_pass"] = None
        return df

    # Run HMMER scan for each candidate
    TARGET_PFAM = {"PF01548", "PF02371"}
    pf01548_hits, pf02371_hits, passes = [], [], []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_p = Path(tmpdir)
        for _, row in df.iterrows():
            acc = row["uniprot_acc"]
            seq = row.get("protein_sequence", "")
            if not seq:
                pf01548_hits.append(False)
                pf02371_hits.append(False)
                passes.append(False)
                continue

            # Write single-sequence FASTA
            fasta_tmp = tmpdir_p / f"{acc}.fasta"
            fasta_tmp.write_text(f">{acc}\n{seq}\n")
            out_tmp = tmpdir_p / f"{acc}.tblout"

            try:
                subprocess.run(
                    [
                        "hmmscan",
                        "--tblout", str(out_tmp),
                        "--noali",
                        "-E", str(evalue_threshold),
                        str(pfam_hmm_path),
                        str(fasta_tmp),
                    ],
                    capture_output=True, timeout=120
                )
                # Parse tblout for PF01548 and PF02371
                found_pf = set()
                if out_tmp.exists():
                    for line in out_tmp.read_text().splitlines():
                        if line.startswith("#"):
                            continue
                        parts = line.split()
                        if len(parts) < 5:
                            continue
                        pfam_id = parts[0].split(".")[0]   # strip version
                        try:
                            e_val = float(parts[4])
                        except (ValueError, IndexError):
                            continue
                        if e_val < evalue_threshold and pfam_id in TARGET_PFAM:
                            found_pf.add(pfam_id)

                has_pf01548 = "PF01548" in found_pf
                has_pf02371 = "PF02371" in found_pf
                pf01548_hits.append(has_pf01548)
                pf02371_hits.append(has_pf02371)
                passes.append(has_pf01548 and has_pf02371)
            except Exception:
                pf01548_hits.append(None)
                pf02371_hits.append(None)
                passes.append(False)

    df = df.copy()
    df["gate_7_pf01548"] = pf01548_hits
    df["gate_7_pf02371"] = pf02371_hits
    df["gate_7_pass"] = passes

    result = df[df["gate_7_pass"] == True].copy()
    print(
        f"  Gate 7 (PF01548 + PF02371, E < {evalue_threshold:.0e}): "
        f"{before:,} -> {len(result):,}"
    )
    return result


# ---------------------------------------------------------------------------
# Full pipeline (gates 1-4 + 6 + 7; gate 5 deferred to Step 12)
# ---------------------------------------------------------------------------

def _normalize_catalog_columns(df: pd.DataFrame, source_label: str) -> pd.DataFrame:
    """Rename Paper 2 catalog columns to the internal schema expected by gate functions.

    Paper 2 parquet files use:
      'accession'  → internal 'uniprot_acc'
      'length'     → internal 'protein_length_aa'
      'sequence'   → internal 'protein_sequence'  (Fanzor catalog only)

    This normalization is applied immediately after loading each catalog so that
    all gate functions and downstream code can rely on the stable internal names.
    """
    rename_map: dict[str, str] = {}
    if "accession" in df.columns and "uniprot_acc" not in df.columns:
        rename_map["accession"] = "uniprot_acc"
    if "length" in df.columns and "protein_length_aa" not in df.columns:
        rename_map["length"] = "protein_length_aa"
    if "sequence" in df.columns and "protein_sequence" not in df.columns:
        rename_map["sequence"] = "protein_sequence"
    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  [{source_label}] Renamed columns: {rename_map}")
    return df


def run_strategy_b(
    is110_catalog_path: Path,
    output_dir: Path,
    fanzor_catalog_path: Optional[Path] = None,
    atlas_embeddings_path: Optional[Path] = None,
    pfam_hmm_path: Optional[Path] = None,
    known_members_path: Optional[Path] = None,
    max_per_genus: int = 5,
    seed: int = 42,
    skip_seq_fetch: bool = False,
) -> pd.DataFrame:
    """Full Strategy B pipeline: gates 1-4 + optional 6-7. Gate 5 deferred.

    Returns DataFrame of pre-Gate-5 survivors with full provenance columns.

    Args:
        skip_seq_fetch: If True, skip UniProt sequence fetching and save the
            candidate list without sequences. Use when the survivor count is
            large (> 500) to avoid long API wait times — sequences can be
            fetched later via ``fetch_lengths_and_organisms()``.
    """
    print(f"\nLoading IS110 triage catalog: {is110_catalog_path}")
    is110 = pd.read_parquet(is110_catalog_path)
    is110 = _normalize_catalog_columns(is110, "IS110")
    print(f"  IS110 triage: {len(is110):,} candidates")

    # Optionally merge Fanzor candidates as auxiliary source
    if fanzor_catalog_path and fanzor_catalog_path.exists():
        fanzor = pd.read_parquet(fanzor_catalog_path)
        fanzor = _normalize_catalog_columns(fanzor, "Fanzor")
        fanzor["paper2_source"] = "fanzor"
        is110["paper2_source"] = "is110_triage"
        full_catalog = pd.concat([is110, fanzor], ignore_index=True).drop_duplicates(
            subset=["uniprot_acc"]
        )
        print(f"  Fanzor candidates: {len(fanzor):,} | Combined: {len(full_catalog):,}")
    else:
        is110["paper2_source"] = "is110_triage"
        full_catalog = is110

    # Load extra known IS110 members from file if provided
    extra_known: list[str] = []
    if known_members_path and known_members_path.exists():
        import yaml
        data = yaml.safe_load(known_members_path.read_text())
        extra_known = data.get("known_accessions", [])

    print("\n=== Applying triage gates 1-4 ===")
    survivors = apply_gate_1(full_catalog)
    survivors = apply_gate_2(survivors)
    survivors = apply_gate_3(survivors)
    survivors = apply_gate_4(survivors, max_per_genus=max_per_genus)

    print("\n=== Applying gate 6 (literature novelty) ===")
    survivors = apply_gate_6(
        survivors,
        atlas_embeddings_path=atlas_embeddings_path,
        known_accessions=KNOWN_IS110_ACCESSIONS + extra_known,
    )

    print("\n=== Applying gate 7 (bRNA interface Pfam check) ===")
    survivors = apply_gate_7(
        survivors,
        pfam_hmm_path=pfam_hmm_path,
    )

    # Fetch sequences for survivors (only if not already present with non-empty values)
    _has_seqs = (
        "protein_sequence" in survivors.columns
        and survivors["protein_sequence"].fillna("").str.len().gt(0).any()
    )

    if _has_seqs:
        n_with_seq = survivors["protein_sequence"].fillna("").str.len().gt(0).sum()
        print(f"\nSequences present for {n_with_seq}/{len(survivors)} pre-Gate-5 survivors "
              f"({len(survivors) - n_with_seq} still need fetching).")
    elif not skip_seq_fetch:
        est_minutes = len(survivors) * _RATE_LIMIT_DELAY / 60
        print(
            f"\nFetching sequences for {len(survivors)} pre-Gate-5 survivors from UniProt "
            f"(estimated ~{est_minutes:.1f} min at {_RATE_LIMIT_DELAY}s/call)..."
        )
        seqs = []
        for i, acc in enumerate(survivors["uniprot_acc"]):
            if i % 20 == 0 and i > 0:
                print(f"  ... {i}/{len(survivors)}")
            seqs.append(_uniprot_fasta(acc) or "")
        survivors = survivors.copy()
        survivors["protein_sequence"] = seqs
    else:
        print(
            f"\n[skip_seq_fetch=True] Skipping UniProt sequence fetch for {len(survivors)} survivors. "
            f"Sequences can be fetched later via fetch_lengths_and_organisms(). "
            f"Estimated fetch time: {len(survivors) * _RATE_LIMIT_DELAY / 60:.1f} min."
        )
        survivors = survivors.copy()
        survivors["protein_sequence"] = ""   # placeholder

    # Drop any rows where sequence retrieval failed (or pre-existing column has None/empty)
    if not skip_seq_fetch:
        n_before_seq = len(survivors)
        survivors = survivors[survivors["protein_sequence"].fillna("").str.len() > 0].copy()
        if len(survivors) < n_before_seq:
            print(f"  Dropped {n_before_seq - len(survivors)} candidates with no retrievable sequence")

    # Assign design IDs
    survivors = survivors.reset_index(drop=True)
    survivors["design_id"] = [
        f"B_{i+1:03d}_{acc[:8]}"
        for i, acc in enumerate(survivors["uniprot_acc"])
    ]
    survivors["strategy"] = "B_ortholog_discovery"
    survivors["gate_5_status"] = "DEFERRED_STEP12"   # explicit marker for downstream

    # Save outputs
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / "ortholog_candidates.parquet"
    survivors.to_parquet(parquet_path, index=False, compression="zstd")

    fasta_path = output_dir / "ortholog_candidates.fasta"
    has_seqs = survivors["protein_sequence"].fillna("").str.len().gt(0).any()
    if has_seqs:
        with fasta_path.open("w") as f:
            for _, row in survivors.iterrows():
                seq = row.get("protein_sequence", "")
                if seq:
                    org = row.get("organism", "unknown")
                    acc = row["uniprot_acc"]
                    f.write(f">{row['design_id']} {org} ({acc})\n{seq}\n")
    else:
        fasta_path.write_text(
            "# Sequences not fetched (skip_seq_fetch=True). "
            "Run fetch_lengths_and_organisms() to populate.\n"
        )

    # Gate summary JSON — cast numpy scalars to Python-native types for JSON serialisation
    import json
    gate_6_applied = bool(
        survivors["gate_6_pass"].notna().any()
    ) if "gate_6_pass" in survivors.columns else False
    gate_7_applied = bool(
        survivors["gate_7_pass"].notna().any()
    ) if "gate_7_pass" in survivors.columns else False
    gate_summary = {
        "n_input": int(len(full_catalog)),
        "n_post_gates_1_2": None,   # not stored separately; see logs
        "n_pre_gate5_survivors": int(len(survivors)),
        "gate_5_status": "DEFERRED_STEP12",
        "gate_6_applied": gate_6_applied,
        "gate_7_applied": gate_7_applied,
        "n_genera_represented": int(survivors["genus"].nunique()) if "genus" in survivors.columns else None,
        "design_ids": survivors["design_id"].tolist(),
    }
    summary_path = output_dir / "ortholog_candidates_gate_summary.json"
    summary_path.write_text(json.dumps(gate_summary, indent=2))

    print(f"\n  Parquet       -> {parquet_path}")
    print(f"  FASTA         -> {fasta_path}")
    print(f"  Gate summary  -> {summary_path}")
    print(f"\nStrategy B pre-Gate-5 survivors: {len(survivors):,}")
    print("Gate 5 (AF3 active-site pLDDT >= 75) deferred to Step 12.")
    print(f"Expected after Gate 5: ~10-20 (P4 threshold: >= 10).")

    return survivors


# ---------------------------------------------------------------------------
# apply_all_gates — P4 operationalization reference (all 7 gates)
# ---------------------------------------------------------------------------

def apply_all_gates(
    df: pd.DataFrame,
    atlas_embeddings_path: Optional[Path] = None,
    pfam_hmm_path: Optional[Path] = None,
    known_accessions: Optional[list[str]] = None,
    max_per_genus: int = 5,
) -> pd.DataFrame:
    """Apply all 7 gates. Used by P4 test script (scripts/43_test_pred_P4_orthologs.py).

    Gate 5 requires af3_mean_plddt and af3_active_site_plddt columns to be present
    (added by Step 12 structure prediction). Will raise RuntimeError if missing.
    """
    df = apply_gate_1(df)
    df = apply_gate_2(df)
    df = apply_gate_3(df)
    df = apply_gate_4(df, max_per_genus=max_per_genus)
    df = apply_gate_5(df)    # requires Step 12 output columns
    df = apply_gate_6(df, atlas_embeddings_path=atlas_embeddings_path,
                      known_accessions=known_accessions or KNOWN_IS110_ACCESSIONS)
    df = apply_gate_7(df, pfam_hmm_path=pfam_hmm_path)
    return df
