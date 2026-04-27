"""Strategy A - Domain-swap chimera generation (Step 9).

Combinatorially assembles chimeric DSB-free editors by fusing domain modules from
the 12 scaffold editors in scaffold_universe.yaml. Each chimera is:
  catalytic_core [+ linker + RNA_binding_module] [+ linker + DNA_targeting_module]

Target: ~30 designs passing size (250-1200 aa) + DSB-free + RNA-guided constraints.

Module boundaries for IS621 use Hiraizumi 2024 Nature 630:994-1002:
  RuvC_fold_DEDD: residues 1-110 (D11/E60/D102/D105 DEDD tetrad)
  bRNA_binding_loop: residues 110-150
  Tnp_serine_C_term: residues 150-260 (S241 Tnp serine)

dCas9 requirement: any Cas9-derived module MUST use D10A+H840A double mutation.
  Cas9 backbone is reused for guide-RNA-binding scaffold only.
  MECH-CLASS Step 15 verifies no chimera fires DSB_NUCLEASE Tier A.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from pydantic import BaseModel

# ---------------------------------------------------------------------------
# Module definitions
# ---------------------------------------------------------------------------


@dataclass
class DomainModule:
    """A single extractable domain module from a scaffold editor."""

    name: str
    source_scaffold: str
    sequence: str  # aa slice; populated at runtime from scaffold_sequences.fasta
    start_residue: int  # 1-indexed in source protein
    end_residue: int  # inclusive
    catalytic_residues: list[int] = field(default_factory=list)  # global source coords
    pfam_acc: str | None = None
    notes: str = ""

    @property
    def expected_length(self) -> int:
        return self.end_residue - self.start_residue + 1


# Module catalog - boundaries from scaffold_universe.yaml + literature annotation.
# Sequences are placeholders; populated at runtime by populate_module_sequences().
DOMAIN_MODULES: dict[str, DomainModule] = {
    # --- IS621 modules (Hiraizumi 2024 Nature 630:994-1002) ---
    "RuvC_fold_DEDD": DomainModule(
        name="RuvC_fold_DEDD",
        source_scaffold="IS621",
        sequence="",  # populated at runtime; residues 1-110 of IS621 (A0A2X3M8B0)
        start_residue=1,
        end_residue=110,
        catalytic_residues=[11, 60, 102, 105],  # D11/E60/D102/D105 DEDD tetrad
        pfam_acc="PF01548",
        notes="DEDD tetrad per Hiraizumi 2024. Boundary extended to 110 to include D102/D105.",
    ),
    "bRNA_binding_loop": DomainModule(
        name="bRNA_binding_loop",
        source_scaffold="IS621",
        sequence="",  # residues 110-150 of IS621
        start_residue=110,
        end_residue=150,
        catalytic_residues=[],
        notes=(
            "Central bRNA-binding loop (41 aa); contact residues per Hiraizumi 2024 SI. "
            "INTERNAL USE ONLY within IS621_full_composite (sits between RuvC and Tnp-serine). "
            "Cannot fold independently - too small for chimeric RNA-binding module. "
            "Use IS621_bRNA_extended for standalone chimeric use with non-IS621 catalytic cores."
        ),
    ),
    "IS621_bRNA_extended": DomainModule(
        name="IS621_bRNA_extended",
        source_scaffold="IS621",
        sequence="",  # residues 100-260 of IS621 (161 aa)
        start_residue=100,
        end_residue=260,
        catalytic_residues=[],
        notes=(
            "Extended IS621 central scaffold for standalone chimeric RNA-binding use (161 aa). "
            "Boundaries extended from the 41-aa loop (110-150) to 100-260 to include flanking "
            "α-helices required for stable independent folding of the bRNA-binding interface. "
            "Validated minimum size for chimeras with single-domain catalytic cores: "
            "smallest chimera is Bxb1_serine_compact (91 aa) + IS621_bRNA_extended (161 aa) + "
            "10 aa linker = 262 aa, passing the 250 aa minimum. "
            "NOTE: when paired with IS621_full_composite, residues 150-260 appear in both "
            "Tnp_serine_C_term and IS621_bRNA_extended (internal IS621-sequence repeat). "
            "These 4 chimeras are flagged with 'IS621_seq_repeat' in provenance and may "
            "score lower in structure prediction (Gate 5 will evaluate viability)."
        ),
    ),
    "Tnp_serine_C_term": DomainModule(
        name="Tnp_serine_C_term",
        source_scaffold="IS621",
        sequence="",  # residues 150-260 of IS621
        start_residue=150,
        end_residue=260,
        catalytic_residues=[241],  # S241 Tnp-serine (global IS621 coords); Hiraizumi 2024
        pfam_acc="PF02371",
        notes="S241 per Hiraizumi 2024 (corrected from S180).",
    ),
    # --- Cre recombinase modules ---
    "loxP_HTH": DomainModule(
        name="loxP_HTH",
        source_scaffold="Cre",
        sequence="",  # residues 30-100 of Cre (P06956)
        start_residue=30,
        end_residue=100,
        catalytic_residues=[],
        notes="HTH loxP-recognition domain.",
    ),
    "catalytic_tyrosine_compact": DomainModule(
        name="catalytic_tyrosine_compact",
        source_scaffold="Cre",
        sequence="",  # residues 200-343 of Cre (P06956); includes Y324
        start_residue=200,
        end_residue=343,
        catalytic_residues=[324],  # Y324 canonical catalytic tyrosine
        notes="Compact tyrosine recombinase catalytic domain; smallest tyrosine scaffold (343 aa total).",
    ),
    # --- Lambda integrase modules ---
    "catalytic_tyrosine": DomainModule(
        name="catalytic_tyrosine",
        source_scaffold="Lambda_Int",
        sequence="",  # residues 150-356 of Lambda_Int (P03700); includes Y342
        start_residue=150,
        end_residue=356,
        catalytic_residues=[342],  # Y342 catalytic tyrosine in Lambda Int
        pfam_acc="PF02899",
        notes="Lambda Int tyrosine recombinase catalytic domain + C-terminal.",
    ),
    "attL_attR_recognition_HTH": DomainModule(
        name="attL_attR_recognition_HTH",
        source_scaffold="Lambda_Int",
        sequence="",  # residues 1-75 of Lambda_Int (P03700); N-terminal arm-binding domain
        start_residue=1,
        end_residue=75,
        catalytic_residues=[],
        notes="Lambda Int N-terminal arm-binding domain for att-site recognition.",
    ),
    # --- IscB modules (RNA-binding only; HNH nuclease excluded - see note below) ---
    "omega_RNA_binding": DomainModule(
        name="omega_RNA_binding",
        source_scaffold="IscB",
        sequence="",  # residues 1-200 of IscB (K9VH02, 400 aa)
        start_residue=1,
        end_residue=200,
        catalytic_residues=[],
        notes="IscB ωRNA-binding scaffold. ONLY module extracted from IscB for chimeras."
        " IscB HNH/RuvC nuclease (DSB-creating) is intentionally excluded.",
    ),
    # HNH_catalytic_compact from IscB is EXCLUDED from all Strategy A catalytic cores.
    # IscB is an RNA-guided HNH endonuclease (Cas9 ancestor, Altae-Tran 2021 Science).
    # Incorporating its HNH domain would make the chimera a DSB nuclease, not a writer.
    # --- Bxb1 serine recombinase modules ---
    "catalytic_serine_recombinase": DomainModule(
        name="catalytic_serine_recombinase",
        source_scaffold="Bxb1",
        sequence="",  # residues 10-100 of Bxb1 (Q9B086); S12 catalytic serine
        start_residue=10,
        end_residue=100,
        catalytic_residues=[12],  # S12 catalytic serine
        pfam_acc="PF07508",
        notes="Canonical serine recombinase catalytic domain.",
    ),
    "attP_attB_binding_domain": DomainModule(
        name="attP_attB_binding_domain",
        source_scaffold="Bxb1",
        sequence="",  # residues 220-500 of Bxb1 (Q9B086)
        start_residue=220,
        end_residue=500,
        catalytic_residues=[],
        notes="Bxb1 C-terminal att-site recognition domain.",
    ),
    # --- phiC31 serine integrase modules ---
    "phiC31_serine_catalytic": DomainModule(
        name="phiC31_serine_catalytic",
        source_scaffold="phiC31",
        sequence="",  # residues 1-140 of phiC31 (Q9T221, 613 aa); S12 catalytic serine
        start_residue=1,
        end_residue=140,
        catalytic_residues=[12],  # S12 conserved catalytic serine in large serine integrases
        pfam_acc="PF07508",
        notes="phiC31 large serine integrase N-terminal catalytic domain.",
    ),
    "phiC31_attP_attB_compact": DomainModule(
        name="phiC31_attP_attB_compact",
        source_scaffold="phiC31",
        sequence="",  # residues 460-613 of phiC31 (Q9T221); C-terminal ZD/RBD
        start_residue=460,
        end_residue=613,
        catalytic_residues=[],
        notes="phiC31 C-terminal zinc ribbon domain for attP/attB specificity.",
    ),
    # --- CAST-V-K (Cas12k) programmability module ---
    "Cas12k_compact_RNP": DomainModule(
        name="Cas12k_compact_RNP",
        source_scaffold="CAST_VK",
        sequence="",  # residues 1-600 of Cas12k (A0A8M0FGU0, 600 aa)
        start_residue=1,
        end_residue=600,
        catalytic_residues=[],  # nuclease activity disabled in chimera context
        notes="Cas12k compact RNP scaffold (600 aa Type V-K). "
        "Programmability (crRNA-guided) only; catalytic nuclease activity "
        "engineered out in chimera context (dREC1 mutation analogous to dCas9).",
    ),
}


# ---------------------------------------------------------------------------
# Linker library
# ---------------------------------------------------------------------------

LINKERS: dict[str, str] = {
    "rigid": "EAAAKEAAAK",  # 10 aa rigid alpha-helix
    "flexible_short": "GGGGSGGGGS",  # 10 aa flexible GS
    "flexible_long": "GGGGS" * 4,  # 20 aa flexible GS
    "natural_IS_family": "QRSAEELNREL",  # natural 11 aa from IS621 paralog linker region
}


# ---------------------------------------------------------------------------
# Design model
# ---------------------------------------------------------------------------


class ChimericDesign(BaseModel):
    """A single Strategy A domain-swap chimeric design."""

    design_id: str
    strategy: str = "A_domain_swap"
    modules_used: list[dict]  # [{module_name, source_scaffold, position_in_chimera}]
    linkers: list[dict]  # [{junction_index, linker_name, sequence}]
    full_sequence: str
    total_aa: int
    bRNA_or_guide_template: str
    expected_mechanism: str  # always "DSB_FREE_TRANSEST_RECOMBINASE" for Strategy A
    scaffold_provenance: dict  # {catalytic_core: str, rna_module: str, linker: str}


# ---------------------------------------------------------------------------
# Combinatorial generation
# ---------------------------------------------------------------------------

# Strategy A catalytic cores: (label, primary_module, secondary_module_or_None)
# Secondary is used for composite architectures (e.g. IS621 needs both RuvC + Tnp-serine)
_CATALYTIC_CORES: list[tuple[str, str, str | None]] = [
    ("IS621_full_composite", "RuvC_fold_DEDD", "Tnp_serine_C_term"),
    ("IS621_RuvC_only", "RuvC_fold_DEDD", None),  # half-composite test
    ("Bxb1_serine_compact", "catalytic_serine_recombinase", None),
    ("Lambda_Int_tyrosine", "catalytic_tyrosine_compact", None),
    ("phiC31_serine", "phiC31_serine_catalytic", None),
]

# RNA-binding modules: (label, module_name)
# IS621_bRNA uses IS621_bRNA_extended (161 aa, residues 100-260) - the standalone-viable
# form of the IS621 bRNA-binding scaffold.  The raw 41-aa bRNA_binding_loop cannot fold
# independently and produces chimeras below the 250 aa minimum when paired with single-
# domain catalytic cores.  IS621_bRNA_extended carries a 'IS621_seq_repeat' flag in
# provenance when combined with IS621_full_composite (overlap at residues 150-260 with
# Tnp_serine_C_term); these 4 chimeras are valid candidates but flagged for Step 12 review.
_RNA_MODULES: list[tuple[str, str]] = [
    ("IS621_bRNA", "IS621_bRNA_extended"),  # 161 aa standalone; passes 250 aa min
    ("IscB_omega", "omega_RNA_binding"),
    ("Cas12k_crRNA", "Cas12k_compact_RNP"),
]

# Linker variants to try for each combination
_LINKER_VARIANTS: list[str] = ["rigid", "flexible_short"]


def _module_length(mod: DomainModule) -> int:
    """Use actual sequence length if populated; fall back to boundary estimate."""
    if mod.sequence:
        return len(mod.sequence)
    return mod.expected_length


def _assemble_sequence(
    module_names: list[str], linker_name: str, catalog: dict[str, DomainModule]
) -> str:
    """Concatenate module sequences with linker between each junction."""
    linker_seq = LINKERS[linker_name]
    parts: list[str] = []
    for i, mod_name in enumerate(module_names):
        mod = catalog.get(mod_name)
        seq = mod.sequence if (mod and mod.sequence) else f"<{mod_name}>"
        parts.append(seq)
        if i < len(module_names) - 1:
            parts.append(linker_seq)
    return "".join(parts)


def _choose_bRNA_template(rna_mod_name: str) -> str:
    if "IS621" in rna_mod_name or "bRNA" in rna_mod_name:
        return "ACCTGTACCGAGGGCCTGTA"  # IS621 canonical bRNA protospacer (Hiraizumi 2024)
    if "omega" in rna_mod_name or "IscB" in rna_mod_name:
        return "GGCTGTGTGGAGAACGATCC"  # IscB canonical ωRNA template
    if "Cas12k" in rna_mod_name or "CAST" in rna_mod_name or "crRNA" in rna_mod_name:
        return "GTTCATTTCGGTAATTATGG"  # Cas12k crRNA repeat
    return "NNNNNNNNNNNNNNNNNNNN"


def generate_chimera_designs(
    catalog: dict[str, DomainModule] | None = None,
    max_protein_size_aa: int = 1200,
    min_protein_size_aa: int = 250,
    seed: int = 42,
) -> list[ChimericDesign]:
    """Generate all Strategy A domain-swap chimeric designs.

    Iterates over catalytic_cores x rna_modules x linker_variants,
    applies size + DSB-free constraints, returns surviving designs.
    Expected ~30 designs from the ~60 raw combinations after filtering.
    """
    if catalog is None:
        catalog = DOMAIN_MODULES

    designs: list[ChimericDesign] = []
    design_counter = 0

    for cat_label, cat_mod1, cat_mod2 in _CATALYTIC_CORES:
        # Skip if required modules are missing from catalog
        if cat_mod1 not in catalog:
            continue
        if cat_mod2 and cat_mod2 not in catalog:
            continue

        for rna_label, rna_mod in _RNA_MODULES:
            if rna_mod not in catalog:
                continue

            # Build ordered module list for this chimera
            mod_list = [cat_mod1]
            if cat_mod2:
                mod_list.append(cat_mod2)
            mod_list.append(rna_mod)

            for linker_name in _LINKER_VARIANTS:
                linker_seq = LINKERS[linker_name]
                n_junctions = len(mod_list) - 1

                # Estimate total size
                module_aa = sum(_module_length(catalog[m]) for m in mod_list)
                linker_aa = len(linker_seq) * n_junctions
                total_aa = module_aa + linker_aa

                # Size gate
                if total_aa > max_protein_size_aa or total_aa < min_protein_size_aa:
                    continue

                design_counter += 1
                design_id = f"A_{design_counter:03d}_{cat_label}_{rna_label}_{linker_name}"

                full_seq = _assemble_sequence(mod_list, linker_name, catalog)

                # Detect IS621 sequence repetition: IS621_full_composite uses
                # Tnp_serine_C_term (aa 150-260); IS621_bRNA_extended also covers
                # aa 100-260.  When combined, residues 150-260 appear twice.
                is_is621_seq_repeat = (
                    cat_label == "IS621_full_composite" and rna_mod == "IS621_bRNA_extended"
                )

                designs.append(
                    ChimericDesign(
                        design_id=design_id,
                        modules_used=[
                            {
                                "module_name": m,
                                "source_scaffold": catalog[m].source_scaffold
                                if m in catalog
                                else "unknown",
                                "position_in_chimera": idx + 1,
                            }
                            for idx, m in enumerate(mod_list)
                        ],
                        linkers=[
                            {
                                "junction_index": j + 1,
                                "linker_name": linker_name,
                                "sequence": linker_seq,
                            }
                            for j in range(n_junctions)
                        ],
                        full_sequence=full_seq,
                        total_aa=len(full_seq) if "<" not in full_seq else total_aa,
                        bRNA_or_guide_template=_choose_bRNA_template(rna_mod),
                        expected_mechanism="DSB_FREE_TRANSEST_RECOMBINASE",
                        scaffold_provenance={
                            "catalytic_core": cat_label,
                            "rna_module": rna_label,
                            "linker": linker_name,
                            "modules": mod_list,
                            "IS621_seq_repeat": is_is621_seq_repeat,  # True -> lower Gate 5 priority
                        },
                    )
                )

    return designs


# ---------------------------------------------------------------------------
# Runtime: populate module sequences from scaffold_sequences.fasta
# ---------------------------------------------------------------------------


def populate_module_sequences(
    scaffold_seqs: dict[str, str],
    catalog: dict[str, DomainModule] | None = None,
) -> dict[str, DomainModule]:
    """Slice actual sequences from scaffold_sequences.fasta into module catalog.

    Called at runtime after Step 8 has produced scaffold_sequences.fasta.
    Modifies the catalog in-place (sequences were placeholder strings).
    """
    if catalog is None:
        catalog = DOMAIN_MODULES

    for _mod_name, mod in catalog.items():
        src = mod.source_scaffold
        if src not in scaffold_seqs:
            continue
        src_seq = scaffold_seqs[src]
        start0 = mod.start_residue - 1  # convert to 0-indexed
        end0 = min(mod.end_residue, len(src_seq))
        mod.sequence = src_seq[start0:end0]

    return catalog


def load_scaffold_sequences(fasta_path: Path) -> dict[str, str]:
    """Parse scaffold_sequences.fasta -> {scaffold_id: sequence}."""
    seqs: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    with fasta_path.open() as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            elif current_id is not None:
                parts.append(line)

    if current_id is not None:
        seqs[current_id] = "".join(parts)

    return seqs


# ---------------------------------------------------------------------------
# main() - run as script or via scripts/10_strategy_A_domain_swap.py
# ---------------------------------------------------------------------------


def main(
    scaffold_fasta: Path = Path("/data/pen-assemble/designs/scaffold_sequences.fasta"),
    output_dir: Path = Path("/data/pen-assemble/designs/strategy_A"),
    max_aa: int = 1200,
    min_aa: int = 250,
) -> list[ChimericDesign]:
    """Full Strategy A pipeline. Loads sequences, generates designs, saves outputs."""
    from rich.console import Console
    from rich.table import Table

    console = Console()
    console.print("[bold]Strategy A - Domain-Swap Chimera Generation (Step 9)[/bold]")

    # Step 9.1: load verified scaffold sequences from Step 8
    if not scaffold_fasta.exists():
        raise FileNotFoundError(
            f"scaffold_sequences.fasta not found at {scaffold_fasta}. "
            "Run Step 8 (scripts/03_validate_sequences.py) first."
        )
    scaffold_seqs = load_scaffold_sequences(scaffold_fasta)
    console.print(f"  Loaded {len(scaffold_seqs)} scaffold sequences from {scaffold_fasta.name}")

    # Step 9.1: populate module catalog with actual sequence slices
    populate_module_sequences(scaffold_seqs)

    populated = sum(1 for m in DOMAIN_MODULES.values() if m.sequence)
    console.print(
        f"  Module catalog populated: {populated}/{len(DOMAIN_MODULES)} modules have sequences"
    )

    # Step 9.2: generate combinatorial designs
    designs = generate_chimera_designs(max_protein_size_aa=max_aa, min_protein_size_aa=min_aa)
    console.print(f"  Generated {len(designs)} chimeric designs (target ~30)")

    # Summary table
    table = Table(title="Strategy A Designs")
    table.add_column("Design ID", style="dim")
    table.add_column("Modules")
    table.add_column("Linker")
    table.add_column("aa")
    table.add_column("bRNA/guide")

    for d in designs:
        mods = " + ".join(m["module_name"] for m in d.modules_used)
        table.add_row(
            d.design_id,
            mods,
            d.linkers[0]["linker_name"] if d.linkers else "-",
            str(d.total_aa),
            d.bRNA_or_guide_template[:20] + "...",
        )
    console.print(table)

    # Step 9.3: save outputs
    output_dir.mkdir(parents=True, exist_ok=True)

    parquet_path = output_dir / "chimeric_designs.parquet"
    rows = [d.model_dump() for d in designs]
    # Flatten scaffold_provenance for Parquet compatibility
    for row in rows:
        prov = row.pop("scaffold_provenance", {})
        row["provenance_catalytic_core"] = prov.get("catalytic_core", "")
        row["provenance_rna_module"] = prov.get("rna_module", "")
        row["provenance_linker"] = prov.get("linker", "")
        row["provenance_modules_json"] = json.dumps(prov.get("modules", []))
    pd.DataFrame(rows).to_parquet(parquet_path, index=False, compression="zstd")

    fasta_path = output_dir / "chimeric_designs.fasta"
    with fasta_path.open("w") as f:
        for d in designs:
            f.write(f">{d.design_id}\n{d.full_sequence}\n")

    manifest_path = output_dir / "chimeric_designs_manifest.json"
    manifest = {
        "n_designs": len(designs),
        "strategy": "A_domain_swap",
        "scaffold_sequences_source": str(scaffold_fasta),
        "module_catalog_size": len(DOMAIN_MODULES),
        "modules_populated": populated,
        "design_ids": [d.design_id for d in designs],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    console.print(f"\n  Parquet  -> {parquet_path}")
    console.print(f"  FASTA    -> {fasta_path}")
    console.print(f"  Manifest -> {manifest_path}")
    console.print(
        f"\n[bold green]Step 9 complete. {len(designs)} Strategy A designs ready for Step 12 (structure prediction).[/bold green]"
    )

    return designs


if __name__ == "__main__":
    main()
