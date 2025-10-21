"""
Module for adding exon and intron structure information to AnnData objects.
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
    include_introns: bool = True,
    inplace: bool = True,
    verbose: bool = True
) -> Optional[ad.AnnData]:
    """
    Add exon and intron structure information to AnnData.var from GTF/GFF data.

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
    include_introns : bool, default=True
        Whether to calculate and include intron structure information
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

    structure_df = _create_transcript_structure_df(
        gtf_df,
        transcript_id_col,
        include_introns=include_introns,
        verbose=verbose
    )

    if structure_df.empty:
        print("Warning: No exon structures could be extracted")
        return None if inplace else adata

    # Map structure information to AnnData.var
    if verbose:
        print("Adding structure information to AnnData.var...")

    _add_structure_to_adata_var(adata, structure_df, include_introns=include_introns, verbose=verbose)

    if verbose:
        print(f"Successfully added exon structure information for {len(structure_df)} transcripts")
        if include_introns:
            print(f"  - Intron structures calculated for multi-exon transcripts")

    return None if inplace else adata


def _create_transcript_structure_df(
    gtf_df: pd.DataFrame,
    transcript_id_col: str = 'transcript_id',
    include_introns: bool = True,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Create a DataFrame with transcript structures (exon and intron lengths) from GTF/GFF data.

    Parameters
    ----------
    gtf_df : pd.DataFrame
        DataFrame with genomic coordinates containing required columns
    transcript_id_col : str, default='transcript_id'
        Column name containing transcript identifiers
    include_introns : bool, default=True
        Whether to calculate intron lengths
    verbose : bool, default=True
        Whether to print progress information

    Returns
    -------
    pd.DataFrame
        DataFrame with transcript structure information including introns
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

        # Calculate intron lengths if requested and there are multiple exons
        intron_lengths = []
        if include_introns and len(sorted_group) > 1:
            # Get the sorted starts and ends
            starts = sorted_group['Start'].tolist()
            ends = sorted_group['End'].tolist()

            # Calculate intron lengths (gap between consecutive exons)
            for i in range(len(ends) - 1):
                # Intron is the gap between end of exon i and start of exon i+1
                intron_length = starts[i + 1] - ends[i] - 1
                intron_lengths.append(intron_length)

        # Create structure strings
        exon_structure_string = ','.join(map(str, exon_lengths))
        intron_structure_string = ','.join(map(str, intron_lengths)) if intron_lengths else ''

        # Calculate total lengths
        total_exon_length = sum(exon_lengths)
        total_intron_length = sum(intron_lengths) if intron_lengths else 0

        # Get additional info if available
        gene_id = group['gene_id'].iloc[0] if 'gene_id' in group.columns else None
        chromosome = group['Chromosome'].iloc[0] if 'Chromosome' in group.columns else None
        strand = group['Strand'].iloc[0] if 'Strand' in group.columns else None

        result_dict = {
            'transcript_id': transcript_id,
            'exon_structure': exon_structure_string,
            'exon_lengths': exon_lengths,
            'transcript_length': total_exon_length,
            'n_exons': len(exon_lengths),
            'gene_id': gene_id,
            'chromosome': chromosome,
            'strand': strand
        }

        # Add intron information if calculated
        if include_introns:
            result_dict.update({
                'intron_structure': intron_structure_string,
                'intron_lengths': intron_lengths,
                'total_intron_length': total_intron_length,
                'n_introns': len(intron_lengths)
            })

        results.append(result_dict)

    structure_df = pd.DataFrame(results)

    if verbose:
        print(f"Processed {len(structure_df)} transcripts")
        print("Exon count distribution:")
        print(structure_df['n_exons'].value_counts().sort_index().head(10))

        if include_introns:
            multi_exon = structure_df[structure_df['n_exons'] > 1]
            if len(multi_exon) > 0:
                print(f"Calculated introns for {len(multi_exon)} multi-exon transcripts")
                print("Intron count distribution:")
                print(multi_exon['n_introns'].value_counts().sort_index().head(10))

    return structure_df


def _add_structure_to_adata_var(
    adata: ad.AnnData,
    structure_df: pd.DataFrame,
    include_introns: bool = True,
    verbose: bool = True
) -> None:
    """
    Add structure information (exons and introns) to AnnData.var.

    Parameters
    ----------
    adata : AnnData
        AnnData object to modify
    structure_df : pd.DataFrame
        DataFrame containing structure information
    include_introns : bool, default=True
        Whether to include intron structure information
    verbose : bool
        Whether to print progress information
    """

    # Set transcript_id as index for easier mapping
    structure_df = structure_df.set_index('transcript_id')

    # Initialize new columns with default values
    n_transcripts = adata.n_vars

    # Exon columns
    adata.var['exon_structure'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')
    adata.var['transcript_length'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='float64')
    adata.var['n_exons'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='Int64')

    # Intron columns
    if include_introns:
        adata.var['intron_structure'] = pd.Series([''] * n_transcripts, index=adata.var_names, dtype='object')
        adata.var['total_intron_length'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='float64')
        adata.var['n_introns'] = pd.Series([np.nan] * n_transcripts, index=adata.var_names, dtype='Int64')

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

            # Exon information
            adata.var.loc[transcript_id, 'exon_structure'] = row['exon_structure']
            adata.var.loc[transcript_id, 'transcript_length'] = row['transcript_length']
            adata.var.loc[transcript_id, 'n_exons'] = row['n_exons']

            # Intron information
            if include_introns:
                adata.var.loc[transcript_id, 'intron_structure'] = row['intron_structure']
                adata.var.loc[transcript_id, 'total_intron_length'] = row['total_intron_length']
                adata.var.loc[transcript_id, 'n_introns'] = row['n_introns']

            # Add optional columns
            if 'chromosome' in structure_df.columns and not pd.isna(row['chromosome']):
                adata.var.loc[transcript_id, 'chromosome'] = row['chromosome']

            if 'strand' in structure_df.columns and not pd.isna(row['strand']):
                adata.var.loc[transcript_id, 'strand'] = row['strand']

            if 'gene_id' in structure_df.columns and not pd.isna(row['gene_id']):
                adata.var.loc[transcript_id, 'gene_id_gtf'] = row['gene_id']

            matched_transcripts += 1

    # Store detailed information in uns for more complex analysis
    exon_lengths_dict = {}
    intron_lengths_dict = {}

    for transcript_id in adata.var_names:
        if transcript_id in structure_df.index:
            exon_lengths_dict[transcript_id] = structure_df.loc[transcript_id, 'exon_lengths']
            if include_introns:
                intron_lengths_dict[transcript_id] = structure_df.loc[transcript_id, 'intron_lengths']
        else:
            exon_lengths_dict[transcript_id] = []
            if include_introns:
                intron_lengths_dict[transcript_id] = []

    # Store in uns
    adata.uns['exon_lengths'] = exon_lengths_dict
    if include_introns:
        adata.uns['intron_lengths'] = intron_lengths_dict

    if verbose:
        print(f"Matched structure information for {matched_transcripts}/{len(adata.var_names)} transcripts")
        if matched_transcripts < len(adata.var_names):
            print(f"Warning: {len(adata.var_names) - matched_transcripts} transcripts had no structure information")


def calculate_structure_similarity(
    structure1: List[int],
    structure2: List[int],
    mode: str = 'exon'
) -> float:
    """
    Calculate similarity between two transcript structures.

    Parameters
    ----------
    structure1 : List[int]
        First structure as list of exon/intron lengths
    structure2 : List[int]
        Second structure as list of exon/intron lengths
    mode : str, default='exon'
        Type of structure ('exon' or 'intron')

    Returns
    -------
    float
        Similarity score between 0 and 1 (Jaccard index)
    """
    if not structure1 or not structure2:
        return 0.0

    set1 = set(structure1)
    set2 = set(structure2)

    if not set1 or not set2:
        return 0.0

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    return intersection / union if union > 0 else 0.0


def calculate_combined_structure_similarity(
    exon_structure1: List[int],
    exon_structure2: List[int],
    intron_structure1: List[int],
    intron_structure2: List[int],
    exon_weight: float = 0.6,
    intron_weight: float = 0.4
) -> float:
    """
    Calculate combined similarity using both exon and intron structures.

    Parameters
    ----------
    exon_structure1 : List[int]
        First transcript's exon structure
    exon_structure2 : List[int]
        Second transcript's exon structure
    intron_structure1 : List[int]
        First transcript's intron structure
    intron_structure2 : List[int]
        Second transcript's intron structure
    exon_weight : float, default=0.6
        Weight for exon similarity (must sum with intron_weight to 1.0)
    intron_weight : float, default=0.4
        Weight for intron similarity (must sum with exon_weight to 1.0)

    Returns
    -------
    float
        Combined similarity score between 0 and 1
    """
    if abs(exon_weight + intron_weight - 1.0) > 1e-6:
        raise ValueError("exon_weight and intron_weight must sum to 1.0")

    exon_sim = calculate_structure_similarity(exon_structure1, exon_structure2, mode='exon')

    # Only calculate intron similarity if both transcripts have introns
    if intron_structure1 and intron_structure2:
        intron_sim = calculate_structure_similarity(intron_structure1, intron_structure2, mode='intron')
        combined_sim = exon_weight * exon_sim + intron_weight * intron_sim
    else:
        # If one or both transcripts lack introns, use only exon similarity
        combined_sim = exon_sim

    return combined_sim


# Convenience function for common use case
def add_structure_from_gtf(
    adata: ad.AnnData,
    gtf_file: str,
    include_introns: bool = True,
    inplace: bool = True,
    verbose: bool = True
) -> Optional[ad.AnnData]:
    """
    Convenience function to add exon and intron structure from GTF file.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing transcript data
    gtf_file : str
        Path to GTF/GFF file
    include_introns : bool, default=True
        Whether to calculate and include intron structures
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
        include_introns=include_introns,
        inplace=inplace,
        verbose=verbose
    )
