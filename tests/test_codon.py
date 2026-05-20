"""
Tests for pen_assemble.codon
"""
import pytest
from pen_assemble.codon import (
    CODON_TABLE_HUMAN,
    RESTRICTION_SITES,
    codon_optimise,
    gc_content,
    check_restriction_sites,
    build_expression_orf,
)


class TestCodonTable:
    def test_all_20_amino_acids_present(self):
        aa = set("ACDEFGHIKLMNPQRSTVWY")
        assert aa.issubset(set(CODON_TABLE_HUMAN.keys()))

    def test_all_codons_length_3(self):
        for aa, codon in CODON_TABLE_HUMAN.items():
            if aa == "*":
                continue
            assert len(codon) == 3, f"{aa}: {codon}"

    def test_all_codons_valid_bases(self):
        valid = set("ACGT")
        for aa, codon in CODON_TABLE_HUMAN.items():
            assert set(codon).issubset(valid), f"Bad codon for {aa}: {codon}"

    def test_met_is_ATG(self):
        assert CODON_TABLE_HUMAN["M"] == "ATG"

    def test_trp_is_TGG(self):
        # Trp has only one codon
        assert CODON_TABLE_HUMAN["W"] == "TGG"


class TestCodonOptimise:
    def test_met_ala(self):
        assert codon_optimise("MA") == "ATG" + CODON_TABLE_HUMAN["A"]

    def test_length(self):
        seq = "MARKG"
        dna = codon_optimise(seq)
        assert len(dna) == len(seq) * 3

    def test_case_insensitive(self):
        assert codon_optimise("ma") == codon_optimise("MA")

    def test_stop_codon(self):
        dna = codon_optimise("M*")
        assert dna == "ATG" + CODON_TABLE_HUMAN["*"]

    def test_unknown_amino_acid(self):
        dna = codon_optimise("MXA")
        assert dna[3:6] == "NNN"

    def test_empty_string(self):
        assert codon_optimise("") == ""


class TestGCContent:
    def test_pure_gc(self):
        assert gc_content("GCGC") == pytest.approx(1.0)

    def test_pure_at(self):
        assert gc_content("ATAT") == pytest.approx(0.0)

    def test_mixed(self):
        # GCGCAT: G,C,G,C,A,T -> 4/6
        assert gc_content("GCGCAT") == pytest.approx(4 / 6)

    def test_empty(self):
        assert gc_content("") == 0.0

    def test_case_insensitive(self):
        assert gc_content("gcgcat") == pytest.approx(gc_content("GCGCAT"))

    def test_typical_human_range(self):
        # Human-preferred codon table should yield GC ~50-60%
        from pen_assemble.codon import codon_optimise
        seq = "MDRFFPVIRICKVGFTMEHELHYIGICTAKEKLDVD"
        dna = codon_optimise(seq)
        gc = gc_content(dna)
        assert 0.45 <= gc <= 0.70, f"GC out of expected range: {gc:.2%}"


class TestRestrictionSites:
    def test_ecori_detected(self):
        # EcoRI site: GAATTC
        assert "EcoRI" in check_restriction_sites("AAAGAATTCAAA")

    def test_no_sites(self):
        assert check_restriction_sites("ACGTACGTACGT") == []

    def test_multiple_sites(self):
        # EcoRI + BamHI
        dna = "GAATTCNNNNGGATCC"
        hits = check_restriction_sites(dna)
        assert "EcoRI" in hits
        assert "BamHI" in hits

    def test_returns_sorted(self):
        dna = "GGATCCGAATTC"
        hits = check_restriction_sites(dna)
        assert hits == sorted(hits)

    def test_case_insensitive(self):
        assert check_restriction_sites("gaattc") == check_restriction_sites("GAATTC")


class TestBuildExpressionORF:
    def test_no_kozak_no_stop(self):
        orf = build_expression_orf("MA", kozak=False, stop=False)
        assert orf == "ATGGCC"

    def test_kozak_prepended(self):
        orf = build_expression_orf("MA", kozak=True, stop=False)
        assert orf.startswith("GCCACC")

    def test_stop_appended(self):
        orf = build_expression_orf("MA", kozak=False, stop=True)
        assert orf.endswith("TGA")

    def test_full_orf(self):
        orf = build_expression_orf("MA", kozak=True, stop=True)
        assert orf.startswith("GCCACC")
        assert orf.endswith("TGA")
        # Kozak(6) + ATG(3) + GCC(3) + TGA(3) = 15
        assert len(orf) == 15

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            build_expression_orf("")

    def test_non_met_raises(self):
        with pytest.raises(ValueError, match="Met"):
            build_expression_orf("ARKG")

    def test_length_consistency(self):
        seq = "MDRFFPVIR"
        orf = build_expression_orf(seq, kozak=False, stop=False)
        assert len(orf) == len(seq) * 3
