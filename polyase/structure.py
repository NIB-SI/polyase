"""
Module for adding exon structure information to AnnData objects.
"""

import pandas as pd
import numpy as np
import anndata as ad
from typing import List, Tuple, Optional, Union
try:
    import pyranges as pr
    PYRANGES_AVAILABLE = True
except ImportError:
    PYRANGES_AVAILABLE = False
    print("Warning: pyranges not available. GTF reading will be limited.")


def add_exon_structure(
    adata: ad.AnnData,
    gtf_file: Optional[str] = None,
    gtf_df: Optional[pd.DataFrame] = None,
    transcript_id_col: str = 'transcript_id',
    inplace: bool = True,
    verbose: bool = True
) -> Optional[ad.AnnData]:
    """
    Add exon structure information to AnnData.var from GTF/GFF data.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing transcript data
    gtf_file : str, optional
        Path to GTF/GFF file. Either gtf_file or gtf_df must be provided.
    gtf_df : pd.DataFrame, optional
        DataFrame with GTF/GFF data. Either gtf_file or gtf_df must be provided.
    transcript_id_col : str, default='transcript_id'
        Column name in GTF data containing transcript identifiers
    inplace : bool, default=True
        If True, modify the AnnData object in place. If False, return a copy.
    verbose : bool, default=True
        Whether to print progress information

    Returns
    -------
    AnnData or None
        If inplace=False, returns modified copy of AnnData object.
        If inplace=True, returns None and modifies the input object.

    Raises
    ------
    ValueError
        If neither gtf_file nor gtf_df is provided, or if required columns are missing
    """

    # Input validation
    if gtf_file is None and gtf_df is None:
        raise ValueError("Either gtf_file or gtf_df must be provided")

    if gtf_file is not None and gtf_df is not None:
        raise ValueError("Provide either gtf_file or gtf_df, not both")

    # Work on copy if not inplace
    if not inplace:
        adata = adata.copy()

    # Load GTF data if file path provided
    if gtf_file is not None:
        if not PYRANGES_AVAILABLE:
            raise ImportError("pyranges is required to read GTF files. Install with: pip install pyranges")

        if verbose:
            print(f"Loading GTF file: {gtf_file}")
        try:
            gtf_ranges = pr.read_gtf(gtf_file)
            gtf_df = gtf_ranges.df
        except Exception as e:
            raise ValueError(f"Error reading GTF file {gtf_file}: {str(e)}")

    # Validate GTF dataframe
    required_cols = ['Feature', 'Start', 'End']
    missing_cols = [col for col in required_cols if col not in gtf_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in GTF data: {missing_cols}")

    if transcript_id_col not in gtf_df.columns:
        raise ValueError(f"Transcript ID column '{transcript_id_col}' not found in GTF data")

    # Create structure dataframe
    if verbose:
        print("Processing exon structures...")

    structure_df = _create_transcript_structure_df(gtf_df, transcript_id_col, verbose)

    if structure_df.empty:
        print("Warning: No exon structures could be extracted")
        return None if inplace else adata

    # Map structure information to AnnData.var
    if verbose:
        print("Adding structure information to AnnData.var...")

    _add_structure_to_adata_var(adata, structure_df, verbose)

    if verbose:
        print(f"Successfully added exon structure information for {len(structure_df)} transcripts")

    return None if inplace else adata


def _create_transcript_structure_df(
    gtf_df: pd.DataFrame,
    transcript_id_col: str = 'transcript_id',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Create a DataFrame with transcript structures (exon lengths) from GTF/GFF data.

    Parameters
    ----------
    gtf_df : pd.DataFrame
        DataFrame with genomic coordinates containing required columns
    transcript_id_col : str, default='transcript_id'
        Column name containing transcript identifiers
    verbose : bool, default=True
        Whether to print progress information

    Returns
    -------
    pd.DataFrame
        DataFrame with transcript structure information
    """

    # Filter for exon features only
    exon_df = gtf_df[gtf_df['Feature'] == 'exon'].copy()

    if exon_df.empty:
        if verbose:
            print("Warning: No exon features found in the data")
        return pd.DataFrame()

    # Calculate exon lengths
    exon_df['exon_length'] = exon_df['End'] - exon_df['Start'] + 1

    # Group by transcript_id to get structure for each transcript
    results = []

    for transcript_id, group in exon_df.groupby(transcript_id_col):
        # Sort exons by position (accounting for strand if available)
        if 'Strand' in group.columns:
            if group['Strand'].iloc[0] == '-':
                # For negative strand, sort by decreasing start position
                sorted_group = group.sort_values('Start', ascending=False)
            else:
                # For positive strand, sort by increasing start position
                sorted_group = group.sort_values('Start', ascending=True)
        else:
            # Default to sorting by start position if strand info not available
            sorted_group = group.sort_values('Start', ascending=True)

        # If exon_number column exists and is not all NaN, use it for sorting
        if 'exon_number' in sorted_group.columns and not sorted_group['exon_number'].isna().all():
            sorted_group = group.sort_values('exon_number')

        # Extract exon lengths in order
        exon_lengths = sorted_group['exon_length'].tolist()

        # Create structure string
        structure_string = ','.join(map(str, exon_lengths))

        # Calculate total length
        total_length = sum(exon_lengths)

        # Get additional info if available
        gene_id = group['gene_id'].iloc[0] if 'gene_id' in group.columns else None
        chromosome = group['Chromosome'].iloc[0] if 'Chromosome' in group.columns else None
        strand = group['Strand'].iloc[0] if 'Strand' in group.columns else None

        results.append({
            'transcript_id': transcript_id,
            'exon_structure': structure_string,
            'exon_lengths': exon_lengths,
            'transcript_length': total_length,
            'n_exons': len(exon_lengths),
            'gene_id': gene_id,
            'chromosome': chromosome,
            'strand': strand
        })

    structure_df = pd.DataFrame(results)

    if verbose:
        print(f"Processed {len(structure_df)} transcripts")
        print("Exon count distribution:")
        print(structure_df['n_exons'].value_counts().sort_index().head(10))

    return structure_df


def _add_structure_to_adata_var(
    adata: ad.AnnData,
    structure_df: pd.DataFrame,
    verbose: bool = True
) -> None:
    """
    Add structure information to AnnData.var.

    Parameters
    ----------
    adata : AnnData
        AnnData object to modify
    structure_df : pd.DataFrame
        DataFrame containing structure information
    verbose : bool
        Whether to print progress information
    """

    # Set transcript_id as index for easier mapping
    structure_df = structure_df.set_index('transcript_id')

    # Initialize new columns with default values
    n_transcripts = adata.n_vars
    adata.var['exon_structure'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')
    adata.var['transcript_length'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='float64')
    adata.var['n_exons'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='Int64')

    # Add optional columns if they exist in structure_df
    if 'chromosome' in structure_df.columns:
        adata.var['chromosome'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')

    if 'strand' in structure_df.columns:
        adata.var['strand'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')

    if 'gene_id' in structure_df.columns:
        adata.var['gene_id_gtf'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')

    # Map structure information to AnnData transcripts
    matched_transcripts = 0

    for transcript_id in adata.var_names:
        if transcript_id in structure_df.index:
            row = structure_df.loc[transcript_id]

            adata.var.loc[transcript_id, 'exon_structure'] = row['exon_structure']
            adata.var.loc[transcript_id, 'transcript_length'] = row['transcript_length']
            adata.var.loc[transcript_id, 'n_exons'] = row['n_exons']

            # Add optional columns
            if 'chromosome' in structure_df.columns and not pd.isna(row['chromosome']):
                adata.var.loc[transcript_id, 'chromosome'] = row['chromosome']

            if 'strand' in structure_df.columns and not pd.isna(row['strand']):
                adata.var.loc[transcript_id, 'strand'] = row['strand']

            if 'gene_id' in structure_df.columns and not pd.isna(row['gene_id']):
                adata.var.loc[transcript_id, 'gene_id_gtf'] = row['gene_id']

            matched_transcripts += 1

    # Store exon_lengths in uns for more complex analysis
    exon_lengths_dict = {}
    for transcript_id in adata.var_names:
        if transcript_id in structure_df.index:
            exon_lengths_dict[transcript_id] = structure_df.loc[transcript_id, 'exon_lengths']
        else:
            exon_lengths_dict[transcript_id] = []

    # Store detailed exon lengths information in uns
    adata.uns['exon_lengths'] = exon_lengths_dict

    if verbose:
        print(f"Matched structure information for {matched_transcripts}/{len(adata.var_names)} transcripts")
        if matched_transcripts < len(adata.var_names):
            print(f"Warning: {len(adata.var_names) - matched_transcripts} transcripts had no structure information")


# Convenience function for common use case
def add_structure_from_gtf(
    adata: ad.AnnData,
    gtf_file: str,
    inplace: bool = True,
    verbose: bool = True
) -> Optional[ad.AnnData]:
    """
    Convenience function to add exon structure from GTF file.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing transcript data
    gtf_file : str
        Path to GTF/GFF file
    inplace : bool, default=True
        If True, modify the AnnData object in place
    verbose : bool, default=True
        Whether to print progress information

    Returns
    -------
    AnnData or None
        Modified AnnData object if inplace=False, otherwise None
    """
    return add_exon_structure(
        adata=adata,
        gtf_file=gtf_file,
        inplace=inplace,
        verbose=verbose
    )
