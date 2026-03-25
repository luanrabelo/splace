import os
import pytest

from splace.phylogeny.run_iqtree import prepare_supermatrix


@pytest.fixture
def synthetic_gene_files(tmp_path):
    """
    Create synthetic FASTA files simulating 3 genes and 3 taxa.

    - gene_A.fasta: has TaxonA, TaxonB, TaxonC (complete)
    - gene_B.fasta: has TaxonA, TaxonB, TaxonC (complete)
    - gene_C.fasta: has TaxonA, TaxonB only (TaxonC missing)
    """
    genes_dir = tmp_path / "genes"
    genes_dir.mkdir()

    # gene_A: 10bp, all 3 taxa
    with open(genes_dir / "gene_A.fasta", "w") as f:
        f.write(">TaxonA\nATCGATCGAT\n")
        f.write(">TaxonB\nGCTAGCTAGC\n")
        f.write(">TaxonC\nTTTTTAAAAA\n")

    # gene_B: 8bp, all 3 taxa
    with open(genes_dir / "gene_B.fasta", "w") as f:
        f.write(">TaxonA\nAAAAAAAA\n")
        f.write(">TaxonB\nCCCCCCCC\n")
        f.write(">TaxonC\nGGGGGGGG\n")

    # gene_C: 6bp, only TaxonA and TaxonB (TaxonC is missing)
    with open(genes_dir / "gene_C.fasta", "w") as f:
        f.write(">TaxonA\nATATAT\n")
        f.write(">TaxonB\nGCGCGC\n")

    files = [
        str(genes_dir / "gene_A.fasta"),
        str(genes_dir / "gene_B.fasta"),
        str(genes_dir / "gene_C.fasta"),
    ]
    return files


def _read_nexus_matrix(nexus_path):
    """Parse NEXUS matrix into dict: {taxon: sequence}."""
    matrix = {}
    in_matrix = False
    with open(nexus_path, "r") as f:
        for line in f:
            line = line.strip()
            if line == "MATRIX":
                in_matrix = True
                continue
            if in_matrix:
                if line.startswith(";"):
                    break
                if line:
                    # Taxon name is in single quotes, followed by spaces and sequence
                    parts = line.split("'")
                    if len(parts) >= 3:
                        taxon = parts[1]
                        seq = parts[2].strip()
                        matrix[taxon] = seq
    return matrix


def _read_nexus_charsets(nexus_path):
    """Parse CHARSET lines from NEXUS."""
    charsets = {}
    with open(nexus_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("CHARSET"):
                # CHARSET gene_A = 1-10;
                parts = line.replace(";", "").split()
                name = parts[1]
                range_str = parts[3]
                charsets[name] = range_str
    return charsets


# ---------------------------------------------------------------------------
# Test: allow_missing=True -> all genes kept, missing filled with '?'
# ---------------------------------------------------------------------------
def test_allow_missing_true(synthetic_gene_files, tmp_path):
    """With allow_missing=True, gene_C should be included and TaxonC filled with '?'."""
    output_dir = str(tmp_path / "output_missing")
    os.makedirs(output_dir, exist_ok=True)

    result = prepare_supermatrix(
        trimmed_files=synthetic_gene_files,
        output_dir=output_dir,
        allow_missing=True,
    )

    assert result is not None
    assert os.path.exists(result)

    matrix = _read_nexus_matrix(result)
    charsets = _read_nexus_charsets(result)

    # All 3 taxa should be present
    assert len(matrix) == 3
    assert "TaxonA" in matrix
    assert "TaxonB" in matrix
    assert "TaxonC" in matrix

    # All 3 genes should be in charsets
    assert "gene_A" in charsets
    assert "gene_B" in charsets
    assert "gene_C" in charsets

    # Total length: 10 + 8 + 6 = 24
    for taxon, seq in matrix.items():
        assert len(seq) == 24, f"{taxon} has length {len(seq)}, expected 24"

    # TaxonC should have '?' for gene_C (last 6 chars)
    taxon_c_seq = matrix["TaxonC"]
    gene_c_region = taxon_c_seq[18:]  # positions 18-23 (0-indexed)
    assert gene_c_region == "??????", f"Expected '??????' for TaxonC gene_C, got '{gene_c_region}'"

    # TaxonA should have real data for gene_C
    taxon_a_seq = matrix["TaxonA"]
    gene_c_a = taxon_a_seq[18:]
    assert gene_c_a == "ATATAT", f"Expected 'ATATAT' for TaxonA gene_C, got '{gene_c_a}'"


# ---------------------------------------------------------------------------
# Test: allow_missing=False -> incomplete gene removed
# ---------------------------------------------------------------------------
def test_allow_missing_false(synthetic_gene_files, tmp_path):
    """With allow_missing=False, gene_C should be excluded (missing TaxonC)."""
    output_dir = str(tmp_path / "output_strict")
    os.makedirs(output_dir, exist_ok=True)

    result = prepare_supermatrix(
        trimmed_files=synthetic_gene_files,
        output_dir=output_dir,
        allow_missing=False,
    )

    assert result is not None
    assert os.path.exists(result)

    matrix = _read_nexus_matrix(result)
    charsets = _read_nexus_charsets(result)

    # All 3 taxa should still be present (they all have gene_A and gene_B)
    assert len(matrix) == 3

    # Only gene_A and gene_B should be in charsets (gene_C removed)
    assert "gene_A" in charsets
    assert "gene_B" in charsets
    assert "gene_C" not in charsets

    # Total length: 10 + 8 = 18 (no gene_C)
    for taxon, seq in matrix.items():
        assert len(seq) == 18, f"{taxon} has length {len(seq)}, expected 18"

    # No '?' characters should be present
    for taxon, seq in matrix.items():
        assert "?" not in seq, f"Found '?' in {taxon} sequence (strict mode)"


# ---------------------------------------------------------------------------
# Test: all genes complete -> no difference between modes
# ---------------------------------------------------------------------------
def test_all_genes_complete(tmp_path):
    """When all genes have all taxa, both modes produce identical results."""
    genes_dir = tmp_path / "complete_genes"
    genes_dir.mkdir()

    with open(genes_dir / "geneX.fasta", "w") as f:
        f.write(">Sp1\nAAAA\n>Sp2\nCCCC\n")
    with open(genes_dir / "geneY.fasta", "w") as f:
        f.write(">Sp1\nGGGG\n>Sp2\nTTTT\n")

    files = [str(genes_dir / "geneX.fasta"), str(genes_dir / "geneY.fasta")]

    out_allow = str(tmp_path / "out_allow")
    out_strict = str(tmp_path / "out_strict")
    os.makedirs(out_allow)
    os.makedirs(out_strict)

    r1 = prepare_supermatrix(files[:], out_allow, allow_missing=True)
    r2 = prepare_supermatrix(files[:], out_strict, allow_missing=False)

    m1 = _read_nexus_matrix(r1)
    m2 = _read_nexus_matrix(r2)

    assert m1 == m2

    c1 = _read_nexus_charsets(r1)
    c2 = _read_nexus_charsets(r2)

    assert c1 == c2


# ---------------------------------------------------------------------------
# Test: no complete genes in strict mode -> returns None
# ---------------------------------------------------------------------------
def test_no_complete_genes_returns_none(tmp_path):
    """If no gene has all taxa and allow_missing=False, return None."""
    genes_dir = tmp_path / "all_missing"
    genes_dir.mkdir()

    # gene1: only Sp1
    with open(genes_dir / "gene1.fasta", "w") as f:
        f.write(">Sp1\nAAAA\n")
    # gene2: only Sp2
    with open(genes_dir / "gene2.fasta", "w") as f:
        f.write(">Sp2\nCCCC\n")

    files = [str(genes_dir / "gene1.fasta"), str(genes_dir / "gene2.fasta")]
    out = str(tmp_path / "out_none")
    os.makedirs(out)

    result = prepare_supermatrix(files, out, allow_missing=False)
    assert result is None


# ---------------------------------------------------------------------------
# Test: default behavior (no allow_missing param) is strict
# ---------------------------------------------------------------------------
def test_default_is_strict(synthetic_gene_files, tmp_path):
    """Default (no allow_missing argument) should behave as strict mode."""
    output_dir = str(tmp_path / "output_default")
    os.makedirs(output_dir, exist_ok=True)

    result = prepare_supermatrix(
        trimmed_files=synthetic_gene_files,
        output_dir=output_dir,
    )

    assert result is not None
    charsets = _read_nexus_charsets(result)
    assert "gene_C" not in charsets, "Default mode should exclude incomplete genes"
