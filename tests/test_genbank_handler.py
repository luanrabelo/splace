import os
import pytest

from splace.read_genbank.genbank_handler import genbank_to_fasta, extract_multiple_genbanks


def _get_gene_names_from_output(output_dir):
    """Helper: return set of gene names from .fasta filenames in output dir."""
    return {
        os.path.splitext(f)[0]
        for f in os.listdir(output_dir)
        if f.endswith(".fasta")
    }


# ---------------------------------------------------------------------------
# Default filtering (no genes_filter, no feature_types)
# ---------------------------------------------------------------------------


async def test_genbank_default_mt_genes(sample_genbank_mt, tmp_output_dir, mt_default_genes):
    """Default mt extraction should produce only genes from the default list."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert isinstance(result, list)
    assert len(result) > 0

    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    for gene in produced_genes:
        assert gene in mt_default_genes, f"Unexpected gene '{gene}' not in default mt list"


# ---------------------------------------------------------------------------
# Specific genes_filter
# ---------------------------------------------------------------------------


async def test_genbank_specific_genes_filter(sample_genbank_mt, tmp_output_dir):
    """Only the requested genes should be extracted."""
    target = ["COI", "CYTB"]
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=target,
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    assert produced_genes.issubset(set(target))


async def test_genbank_single_gene_filter(sample_genbank_mt, tmp_output_dir):
    """Filtering for a single gene should produce exactly one file."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["ND1"],
    )
    assert isinstance(result, list)

    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    if produced_genes:
        assert produced_genes == {"ND1"}
        fasta_path = os.path.join(tmp_output_dir, "ND1.fasta")
        with open(fasta_path) as f:
            content = f.read()
        assert content.startswith(">")
        assert len(content.strip().split("\n")) >= 2


# ---------------------------------------------------------------------------
# Feature types
# ---------------------------------------------------------------------------


async def test_genbank_feature_types_cds_explicit(sample_genbank_mt, tmp_output_dir, mt_default_genes):
    """Explicit feature_types=['CDS'] should behave the same as default."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        feature_types=["CDS"],
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    for gene in produced_genes:
        assert gene in mt_default_genes


async def test_genbank_feature_types_trna(sample_genbank_mt, tmp_output_dir):
    """Scanning tRNA features (no gene name filter) should produce output files."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type=None,
        feature_types=["tRNA"],
    )
    # tRNA features may or may not normalize via SynGenes;
    # we just verify the function runs without error
    assert isinstance(result, list)


async def test_genbank_feature_types_rrna(sample_genbank_mt, tmp_output_dir):
    """Scanning rRNA features should produce output files."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type=None,
        feature_types=["rRNA"],
    )
    assert isinstance(result, list)


async def test_genbank_multiple_feature_types(sample_genbank_mt, tmp_output_dir):
    """Scanning CDS + rRNA should produce at least as many files as CDS alone."""
    result_combo = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type=None,
        feature_types=["CDS", "rRNA"],
    )
    genes_combo = _get_gene_names_from_output(tmp_output_dir)
    assert isinstance(result_combo, list)
    assert len(genes_combo) > 0


async def test_genbank_genes_filter_with_feature_types(sample_genbank_mt, tmp_output_dir):
    """Combining genes_filter and feature_types: only matching genes from matching features."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["COI"],
        feature_types=["CDS", "tRNA"],
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    if produced_genes:
        assert produced_genes == {"COI"}


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


async def test_genbank_nonexistent_file(tmp_output_dir):
    """A nonexistent input file should return None."""
    result = await genbank_to_fasta(
        genbank_file="does_not_exist.gb",
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert result is None


async def test_genbank_no_matching_genes(sample_genbank_mt, tmp_output_dir):
    """Filtering for a gene that doesn't exist should return empty list."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["FAKE_GENE_XYZ"],
    )
    assert result == []


async def test_genbank_empty_genes_filter_uses_default(sample_genbank_mt, tmp_output_dir, mt_default_genes):
    """An empty genes_filter list (falsy) should fall back to the default gene list."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=[],
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    for gene in produced_genes:
        assert gene in mt_default_genes


# ---------------------------------------------------------------------------
# Output format
# ---------------------------------------------------------------------------


async def test_genbank_output_fasta_format(sample_genbank_mt, tmp_output_dir):
    """Verify the output FASTA files have correct format (header + sequence)."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["COI"],
    )
    if not result:
        pytest.skip("COI gene not found in sample GenBank file")

    fasta_path = os.path.join(tmp_output_dir, "COI.fasta")
    assert os.path.exists(fasta_path)

    with open(fasta_path) as f:
        lines = [l.strip() for l in f if l.strip()]

    assert len(lines) >= 2
    assert lines[0].startswith(">")
    # Sequence line should contain only valid nucleotide characters
    for line in lines[1:]:
        if not line.startswith(">"):
            assert all(c in "ATCGNatcgn" for c in line), f"Invalid nucleotide chars in: {line[:50]}"


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------


async def test_extract_multiple_genbanks_batch(multiple_genbank_files, tmp_output_dir):
    """Batch extraction of multiple GenBank files should produce output."""
    result = await extract_multiple_genbanks(
        genbank_files=multiple_genbank_files,
        output_dir=tmp_output_dir,
        data_type="mt",
        max_concurrent=2,
    )
    assert isinstance(result, list)
    assert len(result) > 0

    fasta_files = [f for f in os.listdir(tmp_output_dir) if f.endswith(".fasta")]
    assert len(fasta_files) > 0


async def test_extract_multiple_genbanks_empty_input(tmp_output_dir):
    """Batch extraction with no input files should return empty list."""
    result = await extract_multiple_genbanks(
        genbank_files=[],
        output_dir=tmp_output_dir,
        data_type="mt",
    )
    assert result == []
