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


def _add_structure_to_adata_var(
    adata: ad.AnnData,
    structure_df: pd.DataFrame,
    include_introns: bool = True,
    verbose: bool = True
) -> None:
    """
    Add structure information (exons and introns) to AnnData.var.
    Optimized version using vectorized operations.
    """
    # Set transcript_id as index
    structure_df = structure_df.set_index('transcript_id')
    
    # Reindex to match adata.var_names (handles missing transcripts automatically)
    structure_aligned = structure_df.reindex(adata.var_names)
    
    # Assign columns in bulk (vectorized operations)
    # Convert to string first to avoid categorical dtype issues
    adata.var['exon_structure'] = structure_aligned['exon_structure'].astype(str).replace('nan', '')
    adata.var['transcript_length'] = structure_aligned['transcript_length']
    adata.var['n_exons'] = structure_aligned['n_exons'].astype('Int64')
    
    if include_introns:
        adata.var['intron_structure'] = structure_aligned['intron_structure'].astype(str).replace('nan', '')
        adata.var['total_intron_length'] = structure_aligned['total_intron_length']
        adata.var['n_introns'] = structure_aligned['n_introns'].astype('Int64')
    
    # Optional columns - batch check and assign
    optional_cols = {
        'chromosome': 'chromosome',
        'strand': 'strand',
        'gene_id': 'gene_id_gtf'
    }
    
    for src_col, dest_col in optional_cols.items():
        if src_col in structure_df.columns:
            # Convert to string first to avoid categorical issues
            col_data = structure_aligned[src_col].astype(str)
            # Replace 'nan' string with empty string
            col_data = col_data.replace('nan', '')
            adata.var[dest_col] = col_data
    
    # Store in uns - more efficient conversion
    adata.uns['exon_lengths'] = structure_aligned['exon_lengths'].apply(
        lambda x: x if isinstance(x, list) else []
    ).to_dict()
    
    if include_introns:
        adata.uns['intron_lengths'] = structure_aligned['intron_lengths'].apply(
            lambda x: x if isinstance(x, list) else []
        ).to_dict()
    
    if verbose:
        matched = structure_aligned['exon_structure'].notna().sum()
        print(f"Matched structure information for {matched}/{len(adata.var_names)} transcripts")
        if matched < len(adata.var_names):
            print(f"Warning: {len(adata.var_names) - matched} transcripts had no structure information")


def _create_transcript_structure_df(
    gtf_df: pd.DataFrame,
    transcript_id_col: str = 'transcript_id',
    include_introns: bool = True,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Create a DataFrame with transcript structures (exon and intron lengths) from GTF/GFF data.
    Optimized version with vectorized operations where possible.
    """
    # Filter for exon features only
    exon_df = gtf_df[gtf_df['Feature'] == 'exon'].copy()

    if exon_df.empty:
        if verbose:
            print("Warning: No exon features found in the data")
        return pd.DataFrame()

    # Calculate exon lengths (vectorized)
    exon_df['exon_length'] = exon_df['End'] - exon_df['Start'] + 1

    # Determine sort column once per transcript using vectorized operations
    has_strand = 'Strand' in exon_df.columns
    has_exon_number = 'exon_number' in exon_df.columns and not exon_df['exon_number'].isna().all()
    
    # Pre-sort the dataframe for more efficient groupby
    if has_exon_number:
        # If exon_number exists, use it (fastest)
        exon_df = exon_df.sort_values([transcript_id_col, 'exon_number'])
        sort_by_exon_num = True
    else:
        # Otherwise, sort by position accounting for strand
        if has_strand:
            # Create a sort key that accounts for strand
            exon_df['_sort_key'] = exon_df.apply(
                lambda row: -row['Start'] if row['Strand'] == '-' else row['Start'], 
                axis=1
            )
            exon_df = exon_df.sort_values([transcript_id_col, '_sort_key'])
        else:
            exon_df = exon_df.sort_values([transcript_id_col, 'Start'])
        sort_by_exon_num = False

    # Use groupby with aggregation functions (much faster than iterating)
    agg_dict = {
        'exon_length': ['sum', 'count', list],
        'Start': list,
        'End': list,
    }
    
    # Add optional columns
    if 'gene_id' in exon_df.columns:
        agg_dict['gene_id'] = 'first'
    if 'Chromosome' in exon_df.columns:
        agg_dict['Chromosome'] = 'first'
    if 'Strand' in exon_df.columns:
        agg_dict['Strand'] = 'first'
    
    grouped = exon_df.groupby(transcript_id_col, sort=False).agg(agg_dict)
    
    # Flatten multi-level columns
    grouped.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col 
                       for col in grouped.columns]
    
    # Rename aggregated columns
    column_mapping = {
        'exon_length_sum': 'transcript_length',
        'exon_length_count': 'n_exons',
        'exon_length_list': 'exon_lengths',
        'Start_list': 'starts',
        'End_list': 'ends',
    }
    grouped = grouped.rename(columns=column_mapping)
    
    # Create exon structure strings (vectorized)
    grouped['exon_structure'] = grouped['exon_lengths'].apply(
        lambda x: ','.join(map(str, x))
    )
    
    # Calculate intron information if requested
    if include_introns:
        def calc_introns(row):
            """Calculate intron lengths for a transcript"""
            if row['n_exons'] <= 1:
                return [], '', 0, 0
            
            starts = row['starts']
            ends = row['ends']
            
            # Calculate intron lengths (gap between consecutive exons)
            intron_lengths = [starts[i+1] - ends[i] - 1 for i in range(len(ends) - 1)]
            intron_structure = ','.join(map(str, intron_lengths))
            total_intron_length = sum(intron_lengths)
            n_introns = len(intron_lengths)
            
            return intron_lengths, intron_structure, total_intron_length, n_introns
        
        # Apply intron calculation (still need apply here, but only once per transcript)
        intron_data = grouped.apply(calc_introns, axis=1)
        
        # Unpack results
        grouped['intron_lengths'] = intron_data.apply(lambda x: x[0])
        grouped['intron_structure'] = intron_data.apply(lambda x: x[1])
        grouped['total_intron_length'] = intron_data.apply(lambda x: x[2])
        grouped['n_introns'] = intron_data.apply(lambda x: x[3])
    
    # Drop temporary columns
    grouped = grouped.drop(columns=['starts', 'ends'], errors='ignore')
    
    # Reset index to make transcript_id a column
    structure_df = grouped.reset_index()
    structure_df = structure_df.rename(columns={transcript_id_col: 'transcript_id'})
    
    # Ensure correct column names for optional fields
    rename_map = {}
    if 'gene_id_first' in structure_df.columns:
        rename_map['gene_id_first'] = 'gene_id'
    if 'Chromosome_first' in structure_df.columns:
        rename_map['Chromosome_first'] = 'chromosome'
    if 'Strand_first' in structure_df.columns:
        rename_map['Strand_first'] = 'strand'
    
    if rename_map:
        structure_df = structure_df.rename(columns=rename_map)
    
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

