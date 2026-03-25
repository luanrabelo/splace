import os
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATASETS_DIR = os.path.join(PROJECT_ROOT, "datasets")


@pytest.fixture(autouse=True)
def clear_gene_cache():
    """Clear global gene name caches between tests to ensure isolation."""
    from splace.read_genbank.genbank_handler import GENE_NAME_CACHE as gb_cache
    from splace.read_fasta.fasta_handler import GENE_NAME_CACHE as fa_cache
    gb_cache.clear()
    fa_cache.clear()


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Provide a unique temporary output directory for each test."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return str(output_dir)


@pytest.fixture
def sample_genbank_mt():
    """Path to a single mitochondrial GenBank file."""
    path = os.path.join(
        DATASETS_DIR, "mitochondrial", "genbank", "Pteropodidae",
        "Balionycteris", "Balionycteris_maculata",
        "NC_061537.1_Balionycteris_maculata.gb"
    )
    assert os.path.exists(path), f"Sample GenBank file not found: {path}"
    return path


@pytest.fixture
def multiple_genbank_files():
    """List of 3 GenBank files from Testudinidae for batch tests."""
    base = os.path.join(DATASETS_DIR, "mitochondrial", "genbank", "Testudinidae")
    files = [
        os.path.join(base, "DQ080041.1_Stigmochelys_pardalis.gb"),
        os.path.join(base, "DQ080040.1_Manouria_emys.gb"),
        os.path.join(base, "DQ080049.1_Testudo_graeca.gb"),
    ]
    for f in files:
        assert os.path.exists(f), f"GenBank file not found: {f}"
    return files


@pytest.fixture
def sample_fasta_mt():
    """Path to a mitochondrial FASTA file with [gene=...] headers."""
    path = os.path.join(
        DATASETS_DIR, "mitochondrial", "genbank", "Testudinidae",
        "Chrysemys picta.fasta"
    )
    assert os.path.exists(path), f"Sample FASTA file not found: {path}"
    return path


@pytest.fixture
def multiple_fasta_mt():
    """List of mitochondrial FASTA files for batch tests."""
    base = os.path.join(DATASETS_DIR, "mitochondrial", "genbank", "Testudinidae")
    files = [
        os.path.join(base, "Chrysemys picta.fasta"),
        os.path.join(base, "Graptemys pseudogeographica.fasta"),
    ]
    for f in files:
        assert os.path.exists(f), f"FASTA file not found: {f}"
    return files


@pytest.fixture
def sample_fasta_cp():
    """Path to a chloroplast FASTA file."""
    path = os.path.join(
        DATASETS_DIR, "choroplast", "fastas",
        "MK479229_Acer_yangbiense.fasta"
    )
    assert os.path.exists(path), f"Sample FASTA file not found: {path}"
    return path


@pytest.fixture
def multiple_fasta_cp():
    """List of 3 chloroplast FASTA files for batch tests."""
    base = os.path.join(DATASETS_DIR, "choroplast", "fastas")
    files = [
        os.path.join(base, "MK479229_Acer_yangbiense.fasta"),
        os.path.join(base, "KJ566590_Panax_notoginseng.fasta"),
        os.path.join(base, "KM088018_Panax_quinquefolius.fasta"),
    ]
    for f in files:
        assert os.path.exists(f), f"FASTA file not found: {f}"
    return files


@pytest.fixture
def genes_list_file():
    """Path to the chloroplast genes list text file."""
    path = os.path.join(DATASETS_DIR, "genes_list.txt")
    assert os.path.exists(path), f"Genes list file not found: {path}"
    return path


@pytest.fixture
def mt_default_genes():
    """Default mitochondrial gene list (same as in genbank_handler.py)."""
    return [
        "COI", "COII", "COIII", "CYTB",
        "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
        "ATP6", "ATP8",
    ]
