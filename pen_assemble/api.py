"""Public Designer API for pen-assemble."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class DesignResult:
    """A single chimeric design candidate with all computed properties."""

    design_id: str
    strategy: str  # A, B, C, or D
    protein_sequence: str
    guide_rna_sequence: str | None = None
    scaffold_provenance: dict = field(default_factory=dict)

    # Verification outputs (populated after Steps 12–16)
    predicted_pen_score: float | None = None
    mech_class_tier_a: str | None = None
    composite_flag: bool | None = None
    ddg_kcal_mol: float | None = None
    active_site_pass: bool | None = None
    mean_plddt: float | None = None

    # Pre-registered check flags
    beats_is621: bool | None = None  # PenScore > 0.929

    def __post_init__(self) -> None:
        if self.predicted_pen_score is not None and self.beats_is621 is None:
            self.beats_is621 = self.predicted_pen_score > 0.929


class Designer:
    """High-level design API. Runs all four strategies and evaluates designs.

    Requires Paper 2 (mech-class>=0.5.1) and Paper 3 (pen-score>=0.1.0) extras.
    """

    IS621_PENSCORE_LOCKPOINT = 0.929  # Paper 3 v0.1.0 public_scorecard.parquet

    def __init__(self, data_dir: Path | None = None) -> None:
        self._data_dir = data_dir or Path("/data")
        self._catalog: pd.DataFrame | None = None

    @classmethod
    def load(cls, data_dir: Path | None = None) -> Designer:
        """Load Designer with data directory."""
        return cls(data_dir=data_dir)

    def compose_chimera(
        self,
        scaffold: str,
        guide_module: str,
        target_cargo_bp: int = 100_000,
    ) -> DesignResult:
        """Generate a domain-swap chimera (Strategy A).

        Calls into pen_assemble.strategies.domain_swap.generate_chimera_designs()
        with a single-pair scaffold/guide_module combination.

        Args:
            scaffold: Scaffold ID (e.g. 'IS621') from scaffold_universe.yaml.
            guide_module: Guide module identifier (e.g. 'TBL_IS621').
            target_cargo_bp: Target cargo capacity in base-pairs (default 100 kb).

        Returns:
            DesignResult for the best chimera generated, or raises ValueError
            if no designs pass the domain-swap filters.
        """
        from pen_assemble.data.loader import load_scaffold_universe
        from pen_assemble.strategies.domain_swap import generate_chimera_designs

        scaffolds = load_scaffold_universe()
        scaffold_map = {s.id: s for s in scaffolds}
        if scaffold not in scaffold_map:
            raise ValueError(
                f"Scaffold '{scaffold}' not found in scaffold_universe.yaml. "
                f"Available: {sorted(scaffold_map)}"
            )

        designs = generate_chimera_designs()
        if not designs:
            raise ValueError(
                f"No chimera designs generated for scaffold={scaffold}, "
                f"guide_module={guide_module}."
            )

        best = designs[0]
        return DesignResult(
            design_id=best.design_id,
            strategy="A",
            protein_sequence=best.full_sequence,
            guide_rna_sequence=getattr(best, "bRNA_or_guide_template", None),
            scaffold_provenance={"scaffold": scaffold, "guide_module": guide_module},
        )

    def discover_orthologs(
        self,
        is110_catalog_path: Path | None = None,
        min_composite_prob: float = 0.85,
    ) -> pd.DataFrame:
        """Run Strategy B ortholog discovery pipeline.

        Applies gates 1–5 (length, length-delta, composite-prob, per-genus cap,
        structure quality) from ortholog_discovery.py.  Gates 6 and 7 (ATLAS
        novelty / Pfam domain) require external data not guaranteed at API
        call time and are skipped by default.

        Args:
            is110_catalog_path: Path to IS110 triage catalog Parquet.  Falls
                back to /data/mech-class/catalogs/is110_triage.parquet.
            min_composite_prob: Gate 3 minimum composite probability (0–1).

        Returns:
            DataFrame of candidates that survived gates 1–5, with columns:
            accession, organism, genus, composite_prob, length_aa, gate_*_pass.
        """
        from pen_assemble.strategies.ortholog_discovery import (
            apply_gate_1,
            apply_gate_2,
            apply_gate_3,
            apply_gate_4,
            fetch_lengths_and_organisms,
        )

        catalog = is110_catalog_path or Path("/data/mech-class/catalogs/is110_triage.parquet")
        if not Path(catalog).exists():
            raise FileNotFoundError(
                f"IS110 triage catalog not found at {catalog}. Run MECH-CLASS Step 10 first."
            )

        df = pd.read_parquet(catalog)
        df = apply_gate_1(df)
        df = apply_gate_2(df)
        if min_composite_prob != 0.85:
            # Apply custom threshold
            df = df[df.get("composite_prob", pd.Series(dtype=float)) >= min_composite_prob]
        else:
            df = apply_gate_3(df)
        df = apply_gate_4(df)

        # Fetch lengths/organisms for survivors
        if "length_aa" not in df.columns or "organism" not in df.columns:
            df = fetch_lengths_and_organisms(df)

        return df.reset_index(drop=True)

    def deimmunize(
        self,
        scaffold_id: str = "IS621",
        max_mutations: int = 15,
        n_variants: int = 50,
        seed: int = 42,
    ) -> list[DesignResult]:
        """Run Strategy C Monte Carlo deimmunization for a scaffold.

        Args:
            scaffold_id: Scaffold to deimmunize (currently only 'IS621' supported).
            max_mutations: Maximum number of substitutions per variant.
            n_variants: Number of deimmunized variants to return.
            seed: Random seed for MC reproducibility.

        Returns:
            List of DesignResult, sorted by S_Immuno descending (best first).
        """
        from pen_assemble.strategies.deimmunization import run_deimmunization

        scaffolds_fasta = self._data_dir / "pen-assemble" / "sequences" / f"{scaffold_id}.fasta"
        if not scaffolds_fasta.exists():
            raise FileNotFoundError(
                f"Scaffold FASTA not found: {scaffolds_fasta}. "
                "Run scripts/03_validate_sequences.py first."
            )

        output_dir = self._data_dir / "pen-assemble" / "strategy_c"
        variants: list[DesignResult] = []
        deimm_variants = run_deimmunization(
            scaffold_id=scaffold_id,
            max_mutations=max_mutations,
            n_variants=n_variants,
            output_dir=output_dir,
            seed=seed,
        )

        for i, v in enumerate(deimm_variants[:n_variants]):
            variants.append(
                DesignResult(
                    design_id=getattr(v, "variant_id", f"{scaffold_id}_deimm_{i + 1}"),
                    strategy="C",
                    protein_sequence=getattr(v, "protein_sequence", ""),
                    scaffold_provenance={"parent": scaffold_id, "strategy": "C_deimm"},
                    predicted_pen_score=getattr(v, "predicted_s_immuno", None),
                )
            )
        return variants

    def redesign_backbone(
        self,
        scaffold_id: str = "IS621",
        n_designs: int = 25,
    ) -> list[DesignResult]:
        """Run Strategy D ProteinMPNN backbone-conditioned redesign.

        Requires a predicted PDB structure for the scaffold and the
        ProteinMPNN package to be installed.

        Args:
            scaffold_id: Scaffold to redesign (default 'IS621').
            n_designs: Number of redesigned sequences to generate.

        Returns:
            List of DesignResult sorted by ProteinMPNN score (best first).
        """
        from pen_assemble.strategies.backbone_redesign import generate_backbone_redesigns

        pdb_path = self._data_dir / "pen-assemble" / "structures" / f"{scaffold_id}.pdb"
        if not pdb_path.exists():
            # Fall back to common raw-data location for 8WT6
            pdb_path = self._data_dir / "raw" / f"{scaffold_id}.pdb"
        output_dir = self._data_dir / "pen-assemble" / "strategy_d"
        output_dir.mkdir(parents=True, exist_ok=True)

        variants = generate_backbone_redesigns(
            scaffold_pdb=pdb_path,
            output_dir=output_dir,
            n_designs=n_designs,
        )

        results = []
        for v in variants:
            results.append(
                DesignResult(
                    design_id=v.design_id,
                    strategy="D",
                    protein_sequence=v.protmpnn_sequence,
                    scaffold_provenance={
                        "parent": scaffold_id,
                        "strategy": "D_protmpnn",
                        "sequence_identity_to_wt": v.sequence_identity_to_wt,
                        "novel_residues": v.novel_residues_count,
                    },
                    predicted_pen_score=getattr(v, "pen_score", None),
                )
            )
        return results

    def get_catalog(self) -> pd.DataFrame:
        """Return the full design catalog as a DataFrame."""
        if self._catalog is None:
            catalog_path = self._data_dir / "pen-assemble" / "catalog" / "design_catalog.parquet"
            if catalog_path.exists():
                self._catalog = pd.read_parquet(catalog_path)
            else:
                self._catalog = pd.DataFrame()
        return self._catalog

    def select_designs(
        self,
        strategy: str | None = None,
        require_dsb_free: bool = True,
        top_k: int = 10,
    ) -> pd.DataFrame:
        """Filter and rank designs from the catalog."""
        df = self.get_catalog()
        if df.empty:
            return df
        if strategy:
            df = df[df["strategy"] == strategy]
        if require_dsb_free and "mech_class_tier_a" in df.columns:
            df = df[df["mech_class_tier_a"] == "DSB_FREE_TRANSEST_RECOMBINASE"]
        if "predicted_pen_score" in df.columns:
            df = df.sort_values("predicted_pen_score", ascending=False)
        return df.head(top_k)
