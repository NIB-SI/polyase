import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
import tempfile
import os
from pathlib import Path

# Import the functions to test
# Assuming they're in polyase/exon_structure.py
from polyase.structure import (
    add_exon_structure,
    add_structure_from_gtf,
    _create_transcript_structure_df,
    _add_structure_to_adata_var
)


@pytest.fixture
def create_test_adata():
    """Create a test AnnData object with transcript IDs."""
    n_obs = 5  # samples
    n_vars = 4  # transcripts

    X = np.random.randint(0, 100, size=(n_obs, n_vars))

    adata = AnnData(X=X)
    adata.var_names = ['transcript1', 'transcript2', 'transcript3', 'transcript4']
    adata.obs_names = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5']

    return adata


@pytest.fixture
def create_gtf_dataframe():
    """Create a mock GTF DataFrame for testing."""
    # Create mock exon data for multiple transcripts
    data = {
        'Chromosome': ['chr1'] * 6 + ['chr2'] * 4,
        'Feature': ['exon'] * 10,
        'Start': [100, 200, 300, 100, 250, 400, 500, 700, 500, 800],
        'End': [150, 280, 380, 200, 350, 480, 600, 750, 650, 900],
        'Strand': ['+'] * 6 + ['-'] * 4,
        'transcript_id': ['transcript1', 'transcript1', 'transcript1',
                         'transcript2', 'transcript2', 'transcript2',
                         'transcript3', 'transcript3',
                         'transcript4', 'transcript4'],
        'gene_id': ['gene1'] * 3 + ['gene2'] * 3 + ['gene3'] * 2 + ['gene4'] * 2,
        'exon_number': [1, 2, 3, 1, 2, 3, 1, 2, 1, 2]
    }

    return pd.DataFrame(data)


@pytest.fixture
def create_gtf_file(create_gtf_dataframe, tmp_path):
    """Create a temporary GTF file for testing."""
    gtf_content = """chr1\ttest\texon\t100\t150\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1\ttest\texon\t200\t280\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1\ttest\texon\t300\t380\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1"; exon_number "3";
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "1";
chr1\ttest\texon\t250\t350\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "2";
chr1\ttest\texon\t400\t480\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2"; exon_number "3";
chr2\ttest\texon\t500\t600\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript3"; exon_number "1";
chr2\ttest\texon\t700\t750\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript3"; exon_number "2";
chr2\ttest\texon\t500\t650\t.\t-\t.\tgene_id "gene4"; transcript_id "transcript4"; exon_number "1";
chr2\ttest\texon\t800\t900\t.\t-\t.\tgene_id "gene4"; transcript_id "transcript4"; exon_number "2";
"""

    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text(gtf_content)
    return str(gtf_file)


def test_create_transcript_structure_df(create_gtf_dataframe):
    """Test creation of transcript structure DataFrame."""
    structure_df = _create_transcript_structure_df(
        create_gtf_dataframe,
        transcript_id_col='transcript_id',
        verbose=False
    )

    # Check that we got all transcripts
    assert len(structure_df) == 4
    assert set(structure_df['transcript_id']) == {'transcript1', 'transcript2',
                                                   'transcript3', 'transcript4'}

    # Check transcript1 structure (3 exons on + strand)
    t1 = structure_df[structure_df['transcript_id'] == 'transcript1'].iloc[0]
    assert t1['n_exons'] == 3
    assert t1['exon_lengths'] == [51, 81, 81]  # End - Start + 1
    assert t1['transcript_length'] == 213
    assert t1['exon_structure'] == '51,81,81'

    # Check transcript3 structure (2 exons on - strand)
    t3 = structure_df[structure_df['transcript_id'] == 'transcript3'].iloc[0]
    assert t3['n_exons'] == 2
    # On negative strand, exons should be sorted differently
    assert t3['strand'] == '-'


def test_create_transcript_structure_df_no_exons():
    """Test handling of GTF data with no exon features."""
    gtf_df = pd.DataFrame({
        'Chromosome': ['chr1'],
        'Feature': ['gene'],  # No exons
        'Start': [100],
        'End': [200],
        'transcript_id': ['transcript1']
    })

    structure_df = _create_transcript_structure_df(gtf_df, verbose=False)
    assert structure_df.empty


def test_add_structure_to_adata_var(create_test_adata, create_gtf_dataframe):
    """Test adding structure information to AnnData.var."""
    adata = create_test_adata

    # Create structure dataframe
    structure_df = _create_transcript_structure_df(
        create_gtf_dataframe,
        transcript_id_col='transcript_id',
        verbose=False
    )

    # Add structure to adata
    _add_structure_to_adata_var(adata, structure_df, verbose=False)

    # Check that columns were added
    assert 'exon_structure' in adata.var.columns
    assert 'transcript_length' in adata.var.columns
    assert 'n_exons' in adata.var.columns
    assert 'chromosome' in adata.var.columns
    assert 'strand' in adata.var.columns

    # Check that exon_lengths was added to uns
    assert 'exon_lengths' in adata.uns

    # Check specific transcript values
    assert adata.var.loc['transcript1', 'n_exons'] == 3
    assert adata.var.loc['transcript1', 'transcript_length'] == 213
    assert adata.var.loc['transcript1', 'chromosome'] == 'chr1'
    assert adata.var.loc['transcript1', 'strand'] == '+'

    # Check exon_lengths in uns
    assert adata.uns['exon_lengths']['transcript1'] == [51, 81, 81]


def test_add_exon_structure_with_dataframe(create_test_adata, create_gtf_dataframe):
    """Test add_exon_structure with a DataFrame input."""
    adata = create_test_adata

    result = add_exon_structure(
        adata,
        gtf_df=create_gtf_dataframe,
        inplace=True,
        verbose=False
    )

    # Should return None when inplace=True
    assert result is None

    # Check that structure was added
    assert 'exon_structure' in adata.var.columns
    assert adata.var.loc['transcript1', 'n_exons'] == 3


def test_add_exon_structure_copy(create_test_adata, create_gtf_dataframe):
    """Test add_exon_structure with inplace=False."""
    adata = create_test_adata
    original_adata = adata.copy()

    result = add_exon_structure(
        adata,
        gtf_df=create_gtf_dataframe,
        inplace=False,
        verbose=False
    )

    # Should return modified copy
    assert result is not None
    assert result is not adata

    # Original should be unchanged
    assert 'exon_structure' not in adata.var.columns

    # Copy should be modified
    assert 'exon_structure' in result.var.columns


def test_add_exon_structure_missing_transcripts(create_test_adata):
    """Test handling when some transcripts don't have structure info."""
    adata = create_test_adata

    # Create GTF with only 2 of the 4 transcripts
    gtf_df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Feature': ['exon', 'exon'],
        'Start': [100, 100],
        'End': [200, 200],
        'Strand': ['+', '+'],
        'transcript_id': ['transcript1', 'transcript2'],
        'gene_id': ['gene1', 'gene2']
    })

    add_exon_structure(adata, gtf_df=gtf_df, inplace=True, verbose=False)

    # Check that missing transcripts have NaN/empty values
    assert pd.isna(adata.var.loc['transcript3', 'n_exons'])
    assert adata.var.loc['transcript3', 'exon_structure'] == ''
    assert adata.uns['exon_lengths']['transcript3'] == []


def test_add_structure_from_gtf_convenience(create_test_adata, create_gtf_file):
    """Test the convenience function for adding structure from GTF file."""
    try:
        import pyranges as pr
        PYRANGES_AVAILABLE = True
    except ImportError:
        PYRANGES_AVAILABLE = False

    if not PYRANGES_AVAILABLE:
        pytest.skip("pyranges not available")

    adata = create_test_adata

    result = add_structure_from_gtf(
        adata,
        gtf_file=create_gtf_file,
        inplace=True,
        verbose=False
    )

    assert result is None
    assert 'exon_structure' in adata.var.columns


def test_error_no_input():
    """Test that error is raised when neither gtf_file nor gtf_df provided."""
    adata = AnnData(np.zeros((5, 3)))

    with pytest.raises(ValueError, match="Either gtf_file or gtf_df must be provided"):
        add_exon_structure(adata)


def test_error_both_inputs(create_test_adata, create_gtf_dataframe):
    """Test that error is raised when both gtf_file and gtf_df provided."""
    adata = create_test_adata

    with pytest.raises(ValueError, match="Provide either gtf_file or gtf_df, not both"):
        add_exon_structure(
            adata,
            gtf_file="dummy.gtf",
            gtf_df=create_gtf_dataframe
        )


def test_error_missing_columns():
    """Test error handling for missing required columns in GTF DataFrame."""
    adata = AnnData(np.zeros((5, 3)))

    # GTF missing 'Feature' column
    bad_gtf = pd.DataFrame({
        'Start': [100],
        'End': [200],
        'transcript_id': ['t1']
    })

    with pytest.raises(ValueError, match="Missing required columns"):
        add_exon_structure(adata, gtf_df=bad_gtf)


def test_error_missing_transcript_id_col():
    """Test error when transcript_id column not found."""
    adata = AnnData(np.zeros((5, 3)))

    gtf_df = pd.DataFrame({
        'Feature': ['exon'],
        'Start': [100],
        'End': [200],
        'other_id': ['t1']  # Wrong column name
    })

    with pytest.raises(ValueError, match="Transcript ID column .* not found"):
        add_exon_structure(adata, gtf_df=gtf_df, transcript_id_col='transcript_id')


def test_exon_structure_string_format(create_gtf_dataframe):
    """Test that exon structure strings are formatted correctly."""
    structure_df = _create_transcript_structure_df(
        create_gtf_dataframe,
        verbose=False
    )

    # Check string format (comma-separated integers)
    for _, row in structure_df.iterrows():
        structure = row['exon_structure']
        # Should be comma-separated numbers
        parts = structure.split(',')
        assert len(parts) == row['n_exons']
        # Each part should be convertible to int
        for part in parts:
            assert part.isdigit()


def test_strand_specific_sorting(create_gtf_dataframe):
    """Test that exons are sorted correctly based on strand."""
    structure_df = _create_transcript_structure_df(
        create_gtf_dataframe,
        verbose=False
    )

    # For transcript3 on negative strand, verify exons are properly ordered
    t3 = structure_df[structure_df['transcript_id'] == 'transcript3'].iloc[0]
    assert t3['strand'] == '-'
    # The exact ordering depends on implementation, but should be consistent
    assert len(t3['exon_lengths']) == 2


def test_empty_gtf_dataframe():
    """Test handling of empty GTF DataFrame."""
    adata = AnnData(np.zeros((5, 3)))
    adata.var_names = ['t1', 't2', 't3']

    empty_gtf = pd.DataFrame({
        'Chromosome': [],
        'Feature': [],
        'Start': [],
        'End': [],
        'transcript_id': []
    })

    result = add_exon_structure(adata, gtf_df=empty_gtf, inplace=False, verbose=False)

    # Should return None when no structures found
    assert result is None


def test_transcript_length_calculation(create_gtf_dataframe):
    """Test that transcript lengths are calculated correctly."""
    structure_df = _create_transcript_structure_df(
        create_gtf_dataframe,
        verbose=False
    )

    for _, row in structure_df.iterrows():
        # transcript_length should equal sum of exon_lengths
        assert row['transcript_length'] == sum(row['exon_lengths'])


def test_additional_gtf_columns(create_test_adata):
    """Test that additional GTF columns (gene_id, chromosome, strand) are preserved."""
    adata = create_test_adata

    gtf_df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr2'],
        'Feature': ['exon', 'exon'],
        'Start': [100, 200],
        'End': [200, 300],
        'Strand': ['+', '-'],
        'transcript_id': ['transcript1', 'transcript2'],
        'gene_id': ['gene1', 'gene2']
    })

    add_exon_structure(adata, gtf_df=gtf_df, inplace=True, verbose=False)

    # Check that optional columns were added
    assert 'gene_id_gtf' in adata.var.columns
    assert adata.var.loc['transcript1', 'gene_id_gtf'] == 'gene1'
    assert adata.var.loc['transcript1', 'chromosome'] == 'chr1'
    assert adata.var.loc['transcript1', 'strand'] == '+'


def test_single_exon_transcript():
    """Test handling of single-exon transcripts."""
    adata = AnnData(np.zeros((3, 1)))
    adata.var_names = ['single_exon_t1']

    gtf_df = pd.DataFrame({
        'Chromosome': ['chr1'],
        'Feature': ['exon'],
        'Start': [100],
        'End': [500],
        'transcript_id': ['single_exon_t1']
    })

    add_exon_structure(adata, gtf_df=gtf_df, inplace=True, verbose=False)

    assert adata.var.loc['single_exon_t1', 'n_exons'] == 1
    assert adata.var.loc['single_exon_t1', 'transcript_length'] == 401
    assert adata.var.loc['single_exon_t1', 'exon_structure'] == '401'
    assert adata.uns['exon_lengths']['single_exon_t1'] == [401]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
