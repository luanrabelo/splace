import os
import pytest

from splace.read_fasta.fasta_handler import fasta_extractor, extract_multiple_fastas


def _get_gene_names_from_output(output_dir):
    """Helper: return set of gene names from .fasta filenames in output dir."""
    return {
        os.path.splitext(f)[0]
        for f in os.listdir(output_dir)
        if f.endswith(".fasta")
    }


# ---------------------------------------------------------------------------
# Default filtering
# ---------------------------------------------------------------------------


async def test_fasta_mt_default_genes(sample_fasta_mt, tmp_output_dir, mt_default_genes):
    """Default mt extraction from FASTA with [gene=...] headers."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert isinstance(result, list)
    assert len(result) > 0

    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    for gene in produced_genes:
        assert gene in mt_default_genes, f"Unexpected gene '{gene}'"


# ---------------------------------------------------------------------------
# Specific genes_filter
# ---------------------------------------------------------------------------


async def test_fasta_mt_specific_genes_filter(sample_fasta_mt, tmp_output_dir):
    """Filtering for specific genes should produce only those genes."""
    target = ["ND1", "ND2"]
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=target,
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    assert produced_genes.issubset(set(target))


async def test_fasta_mt_single_gene(sample_fasta_mt, tmp_output_dir):
    """Filtering for COI should produce at most one file."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["COI"],
    )
    assert isinstance(result, list)
    produced_genes = _get_gene_names_from_output(tmp_output_dir)
    if produced_genes:
        assert produced_genes == {"COI"}


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


async def test_fasta_nonexistent_file(tmp_output_dir):
    """A nonexistent input file should return empty list."""
    result = await fasta_extractor(
        fasta_file="nonexistent.fasta",
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert result == []


async def test_fasta_no_matching_genes(sample_fasta_mt, tmp_output_dir):
    """Filtering for a nonexistent gene should return empty list."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["FAKE_GENE"],
    )
    assert result == []


# ---------------------------------------------------------------------------
# Output format
# ---------------------------------------------------------------------------


async def test_fasta_output_content(sample_fasta_mt, tmp_output_dir):
    """Output FASTA should have valid header and sequence."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
    )
    if not result:
        pytest.skip("No genes extracted from sample FASTA")

    # Pick the first output file
    fasta_files = [f for f in os.listdir(tmp_output_dir) if f.endswith(".fasta")]
    assert len(fasta_files) > 0

    fasta_path = os.path.join(tmp_output_dir, fasta_files[0])
    with open(fasta_path) as f:
        lines = [l.strip() for l in f if l.strip()]

    assert len(lines) >= 2
    assert lines[0].startswith(">")


# ---------------------------------------------------------------------------
# Chloroplast FASTA (simple headers - known limitation)
# ---------------------------------------------------------------------------


async def test_fasta_cp_simple_headers(sample_fasta_cp, tmp_output_dir):
    """Chloroplast FASTAs with simple >geneName headers (no [gene=...] tags).

    The fasta_extractor may not parse these simple headers successfully.
    This test documents the current behavior.
    """
    result = await fasta_extractor(
        fasta_file=sample_fasta_cp,
        output_path=tmp_output_dir,
        data_type="cp",
    )
    # The result depends on whether the header parsing can match genes.
    # This test simply verifies no crash occurs.
    assert isinstance(result, list)


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------


async def test_extract_multiple_fastas_mt(multiple_fasta_mt, tmp_output_dir):
    """Batch extraction of multiple FASTA files should produce output."""
    result = await extract_multiple_fastas(
        fasta_files=multiple_fasta_mt,
        output_dir=tmp_output_dir,
        data_type="mt",
        max_concurrent=2,
    )
    assert isinstance(result, list)
    assert len(result) > 0


async def test_extract_multiple_fastas_empty_input(tmp_output_dir):
    """Batch extraction with no input files should return empty list."""
    result = await extract_multiple_fastas(
        fasta_files=[],
        output_dir=tmp_output_dir,
        data_type="mt",
    )
    assert result == []
