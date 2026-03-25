import os
import pytest

from splace.read_genbank.genbank_handler import genbank_to_fasta
from splace.read_fasta.fasta_handler import fasta_extractor


# ---------------------------------------------------------------------------
# Gene filter loading logic
# ---------------------------------------------------------------------------


def test_load_genes_from_text_file(genes_list_file):
    """Parse genes_list.txt using the same logic as splace.py."""
    with open(genes_list_file, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    assert len(genes) == 83
    assert "rbcL" in genes
    assert "matK" in genes
    assert "psbA" in genes
    assert "accD" in genes


def test_parse_comma_separated_genes():
    """Parse a comma-separated gene string."""
    raw = "COI, CYTB, ND1"
    genes = [g.strip() for g in raw.split(",") if g.strip()]

    assert len(genes) == 3
    assert genes == ["COI", "CYTB", "ND1"]


def test_parse_comma_separated_no_spaces():
    """Parse comma-separated without spaces."""
    raw = "12S,16S,COI"
    genes = [g.strip() for g in raw.split(",") if g.strip()]

    assert len(genes) == 3
    assert genes == ["12S", "16S", "COI"]


# ---------------------------------------------------------------------------
# Cross-type filter scenarios
# ---------------------------------------------------------------------------


async def test_cp_filter_on_mt_genbank(sample_genbank_mt, tmp_output_dir, genes_list_file):
    """Using chloroplast gene names to filter mitochondrial GenBank should yield nothing."""
    with open(genes_list_file, "r") as f:
        cp_genes = [line.strip() for line in f if line.strip()]

    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=cp_genes,
    )
    assert result == []


async def test_gene_filter_from_file_with_fasta(sample_fasta_mt, tmp_output_dir):
    """Simulate loading a filter from file and using it with FASTA extraction."""
    genes_filter = ["ND1"]
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=genes_filter,
    )
    assert isinstance(result, list)
    if result:
        assert os.path.exists(os.path.join(tmp_output_dir, "ND1.fasta"))
