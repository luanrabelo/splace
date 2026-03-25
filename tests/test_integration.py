import os
import pytest

from splace.read_genbank.genbank_handler import genbank_to_fasta, extract_multiple_genbanks
from splace.read_fasta.fasta_handler import fasta_extractor, extract_multiple_fastas


def _list_fasta_files(directory):
    """Helper: list .fasta files in a directory."""
    if not os.path.isdir(directory):
        return []
    return [f for f in os.listdir(directory) if f.endswith(".fasta")]


def _gene_names(directory):
    """Helper: get gene names from .fasta filenames."""
    return {os.path.splitext(f)[0] for f in _list_fasta_files(directory)}


def _validate_fasta_format(filepath):
    """Helper: verify a file is valid FASTA."""
    with open(filepath, "r") as f:
        content = f.read().strip()
    assert content, f"File is empty: {filepath}"
    for block in content.split(">")[1:]:
        lines = block.strip().split("\n")
        assert len(lines) >= 2, f"Invalid FASTA block in {filepath}"
        assert lines[0].strip(), f"Empty header in {filepath}"
        seq = "".join(lines[1:]).strip()
        assert seq, f"Empty sequence in {filepath}"


# ---------------------------------------------------------------------------
# Scenario 1: mt GenBank default genes
# ---------------------------------------------------------------------------
async def test_integration_genbank_mt_default(sample_genbank_mt, tmp_output_dir, mt_default_genes):
    """Extract mt GenBank with no filter -> only default 13 genes."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert isinstance(result, (list, set))
    genes = _gene_names(tmp_output_dir)
    assert genes, "No gene files produced"
    for g in genes:
        assert g in mt_default_genes, f"Unexpected gene: {g}"
        _validate_fasta_format(os.path.join(tmp_output_dir, f"{g}.fasta"))


# ---------------------------------------------------------------------------
# Scenario 2: mt GenBank + specific genes
# ---------------------------------------------------------------------------
async def test_integration_genbank_specific_genes(sample_genbank_mt, tmp_output_dir):
    """Extract only COI and CYTB from mt GenBank."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["COI", "CYTB"],
    )
    genes = _gene_names(tmp_output_dir)
    assert genes == {"COI", "CYTB"} or genes.issubset({"COI", "CYTB"})
    for g in genes:
        _validate_fasta_format(os.path.join(tmp_output_dir, f"{g}.fasta"))


# ---------------------------------------------------------------------------
# Scenario 3: mt GenBank + feature types CDS + rRNA
# ---------------------------------------------------------------------------
async def test_integration_genbank_cds_rrna(sample_genbank_mt, tmp_output_dir):
    """Extract CDS and rRNA features from mt GenBank."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type=None,
        feature_types=["CDS", "rRNA"],
    )
    genes = _gene_names(tmp_output_dir)
    assert len(genes) > 0, "No gene files produced for CDS+rRNA"
    for g in genes:
        _validate_fasta_format(os.path.join(tmp_output_dir, f"{g}.fasta"))


# ---------------------------------------------------------------------------
# Scenario 4: mt GenBank + genes from text file
# ---------------------------------------------------------------------------
async def test_integration_genbank_genes_from_file(sample_genbank_mt, tmp_output_dir, genes_list_file):
    """Load cp genes from text file and use as filter on mt GenBank (should produce nothing)."""
    with open(genes_list_file, "r") as f:
        cp_genes = [line.strip() for line in f if line.strip()]
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=cp_genes,
    )
    genes = _gene_names(tmp_output_dir)
    assert len(genes) == 0, "cp gene filter should not match mt data"


# ---------------------------------------------------------------------------
# Scenario 5: mt GenBank + tRNA feature type only
# ---------------------------------------------------------------------------
async def test_integration_genbank_trna_only(sample_genbank_mt, tmp_output_dir):
    """Extract only tRNA features."""
    result = await genbank_to_fasta(
        genbank_file=sample_genbank_mt,
        output_path=tmp_output_dir,
        data_type=None,
        feature_types=["tRNA"],
    )
    genes = _gene_names(tmp_output_dir)
    # tRNA features exist in mitochondrial genomes
    assert isinstance(result, (list, set))


# ---------------------------------------------------------------------------
# Scenario 6: mt FASTA default genes
# ---------------------------------------------------------------------------
async def test_integration_fasta_mt_default(sample_fasta_mt, tmp_output_dir, mt_default_genes):
    """Extract mt FASTA with default gene list."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
    )
    assert isinstance(result, list)
    if result:
        genes = _gene_names(tmp_output_dir)
        for g in genes:
            assert g in mt_default_genes, f"Unexpected gene: {g}"
            _validate_fasta_format(os.path.join(tmp_output_dir, f"{g}.fasta"))


# ---------------------------------------------------------------------------
# Scenario 7: mt FASTA + specific genes
# ---------------------------------------------------------------------------
async def test_integration_fasta_specific_genes(sample_fasta_mt, tmp_output_dir):
    """Extract only ND1 from mt FASTA."""
    result = await fasta_extractor(
        fasta_file=sample_fasta_mt,
        output_path=tmp_output_dir,
        data_type="mt",
        genes_filter=["ND1"],
    )
    genes = _gene_names(tmp_output_dir)
    assert genes.issubset({"ND1"})
    if genes:
        _validate_fasta_format(os.path.join(tmp_output_dir, "ND1.fasta"))


# ---------------------------------------------------------------------------
# Scenario 8: Batch GenBank extraction (3 files)
# ---------------------------------------------------------------------------
async def test_integration_batch_genbank(multiple_genbank_files, tmp_output_dir, mt_default_genes):
    """Batch extract 3 GenBank files. Gene files should accumulate sequences."""
    result = await extract_multiple_genbanks(
        genbank_files=multiple_genbank_files,
        output_dir=tmp_output_dir,
        data_type="mt",
        max_concurrent=2,
    )
    assert isinstance(result, list)
    assert len(result) > 0, "Batch extraction produced no files"
    genes = _gene_names(tmp_output_dir)
    for g in genes:
        assert g in mt_default_genes
        filepath = os.path.join(tmp_output_dir, f"{g}.fasta")
        _validate_fasta_format(filepath)
        # Multiple species should be in the same gene file
        with open(filepath, "r") as f:
            headers = [line for line in f if line.startswith(">")]
        assert len(headers) >= 2, f"Expected multiple species in {g}.fasta, got {len(headers)}"


# ---------------------------------------------------------------------------
# Scenario 9: Batch FASTA extraction (2 files)
# ---------------------------------------------------------------------------
async def test_integration_batch_fasta(multiple_fasta_mt, tmp_output_dir, mt_default_genes):
    """Batch extract 2 mt FASTA files."""
    result = await extract_multiple_fastas(
        fasta_files=multiple_fasta_mt,
        output_dir=tmp_output_dir,
        data_type="mt",
        max_concurrent=2,
    )
    assert isinstance(result, list)
    if result:
        genes = _gene_names(tmp_output_dir)
        for g in genes:
            assert g in mt_default_genes
            _validate_fasta_format(os.path.join(tmp_output_dir, f"{g}.fasta"))


# ---------------------------------------------------------------------------
# Scenario 10: Gene coverage check across batch
# ---------------------------------------------------------------------------
async def test_integration_gene_coverage(multiple_genbank_files, tmp_output_dir):
    """Verify which genes are shared across all species in batch extraction."""
    result = await extract_multiple_genbanks(
        genbank_files=multiple_genbank_files,
        output_dir=tmp_output_dir,
        data_type="mt",
        max_concurrent=2,
    )
    genes = _gene_names(tmp_output_dir)
    num_species = len(multiple_genbank_files)

    complete_genes = []
    for gene in genes:
        filepath = os.path.join(tmp_output_dir, f"{gene}.fasta")
        with open(filepath, "r") as f:
            headers = [line for line in f if line.startswith(">")]
        if len(headers) == num_species:
            complete_genes.append(gene)

    # At least some genes should be present in all species
    assert len(complete_genes) > 0, "No genes are shared across all species"
