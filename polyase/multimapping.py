"""Utilities for calculating mulitmapping ratios per syntelog."""

import numpy as np
import pandas as pd

class MultimappingRatioCalculator:
    """
    Class for calculating multimapping ratios in AnnData objects.
    """

    def __init__(self, adata=None):
        """
        Initialize the calculator with an optional AnnData object.

        Parameters:
        -----------
        adata : AnnData, optional
            AnnData object containing transcript data
        """
        self.adata = adata

    def set_data(self, adata):
        """
        Set or update the AnnData object.

        Parameters:
        -----------
        adata : AnnData
            AnnData object containing transcript data
        """
        self.adata = adata

    def calculate_ratios(self, unique_layer='unique_counts', multi_layer='ambiguous_counts'):
        """
        Calculate multimapping ratios for each transcript grouped by Synt_id.

        Parameters:
        -----------
        unique_layer : str, optional (default: 'unique_counts')
            Layer containing unique counts to use for ratio calculations
        multi_layer : str, optional (default: 'ambiguous_counts')
            Layer containing multimapping counts to use for ratio calculations

        Returns:
        --------
        adata : AnnData
            Updated AnnData object with multimapping ratio layer added
        """
        if self.adata is None:
            raise ValueError("No AnnData object has been set")

        # Make sure Synt_id is in var
        if 'Synt_id' not in self.adata.var:
            raise ValueError("'Synt_id' not found in var")

        # Get counts from specified layer
        if multi_layer not in self.adata.layers:
            raise ValueError(f"Layer '{multi_layer}' not found")

        if unique_layer not in self.adata.layers:
            raise ValueError(f"Layer '{unique_layer}' not found")

        # Get unique Synt_ids
        unique_synt_ids = pd.unique(self.adata.var['Synt_id'])

        # Initialize ratio array with zeros
        #ratio_matrix = np.zeros_like(self.adata.layers[multi_layer], dtype=float)
        multi_ratio = np.zeros(self.adata.shape[1], dtype=float)
        # Calculate ratios for each Synt_id group
        for synt_id in unique_synt_ids:
            # Skip groups with Synt_id = 0 or None if needed
            if synt_id == 0 or synt_id is None:
                continue

            # Create mask for current Synt_id
            mask = self.adata.var['Synt_id'] == synt_id

            # Get counts for this group
            multi_group_counts = self.adata.layers[multi_layer][:,mask]
            unique_group_counts = self.adata.layers[unique_layer][:,mask]

            # Calculate total unique and ambiguous counts for this Synt_id
            multi_total_counts = np.sum(multi_group_counts)
            unique_total_counts = np.sum(unique_group_counts)

            # Calculate ratio if total_counts > 0
            if unique_total_counts > 0:
                multi_ratio[mask] = multi_total_counts / (unique_total_counts + multi_total_counts)

        layer_name = f'multimapping_ratio'

        # Add the ratios as a new layer in the AnnData object
        self.adata.var[layer_name] = multi_ratio

        return self.adata


    def get_ratios_for_synt_id(self, synt_id, multi_layer='multimapping_ratio'):
        """
        Get multimapping ratios for a specific Synt_id.

        Parameters:
        -----------
        synt_id : int or str
            The Synt_id to get ratios for
        multi_layer : str, optional
            Name of the layer containing the ratio data

        Returns:
        --------
        ratios : numpy array
            Array of mulitmapping values for the specified Synt_id
        """
        if multi_layer not in self.adata.layers:
            raise ValueError(f"Mulitmapping layer '{multi_layer}' not found. Calculate ratios first.")

        mask = self.adata.var['Synt_id'] == synt_id
        return self.adata.layers[multi_layer][:,mask]


def calculate_multi_ratios(adata, unique_layer='unique_counts', multi_layer='ambiguous_counts'):
    """
    Calculat emultimapping ratios for each transcript grouped by Synt_id.

    Parameters:
    -----------
    adata : AnnData
        AnnData object containing transcript data
    unique_layer : str, optional (default: 'unique_counts')
        Layer containing counts to use for ratio calculations
    multi_layer : str, optional (default: 'ambiguous_counts')
        Layer containing counts to use for ratio calculations

    Returns:
    --------
    adata : AnnData
        Updated AnnData object with 'multimapping_ratio' layer added
    """
    calculator = MultimappingRatioCalculator(adata)
    return calculator.calculate_ratios(unique_layer=unique_layer, multi_layer=multi_layer)

def calculate_per_allele_ratios(adata, unique_layer='unique_counts', multi_layer='ambiguous_counts',
                               gene_grouping_column='gene_id', inplace=True, count_scaling=True,
                               min_counts_threshold=10, scaling_method='weighted_average'):
    """
    Calculate multimapping ratios for each individual allele/transcript.
    For transcripts from the same gene, assigns ratios based on count-weighted calculations.

    Parameters:
    -----------
    adata : AnnData
        AnnData object containing transcript data
    unique_layer : str, optional (default: 'unique_counts')
        Layer containing unique counts to use for ratio calculations
    multi_layer : str, optional (default: 'ambiguous_counts')
        Layer containing multimapping counts to use for ratio calculations
    gene_grouping_column : str, optional (default: 'gene_id')
        Column in adata.var to group transcripts by (e.g., 'Synt_id', 'gene_id')
    inplace : bool, optional (default: True)
        Whether to modify the input AnnData object or return a copy
    count_scaling : bool, optional (default: True)
        Whether to scale multimapping ratios by count abundance
    min_counts_threshold : int, optional (default: 10)
        Minimum total counts threshold - transcripts below this get reduced weight
    scaling_method : str, optional (default: 'weighted_average')
        Method for combining ratios within genes:
        - 'weighted_average': Weight by total counts
        - 'max_weighted': Take max ratio but weight by counts
        - 'abundance_filtered': Only consider transcripts above threshold

    Returns:
    --------
    adata : AnnData
        Updated AnnData object with per-allele multimapping ratio added to var
    """
    import numpy as np
    import pandas as pd

    if adata is None:
        raise ValueError("No AnnData object provided")

    # Work on a copy if not inplace
    if not inplace:
        adata = adata.copy()

    # Get counts from specified layers
    if multi_layer not in adata.layers:
        raise ValueError(f"Layer '{multi_layer}' not found")
    if unique_layer not in adata.layers:
        raise ValueError(f"Layer '{unique_layer}' not found")

    # Check if gene grouping column exists
    if gene_grouping_column not in adata.var.columns:
        raise ValueError(f"Gene grouping column '{gene_grouping_column}' not found in adata.var")

    # Get counts for each transcript/allele
    unique_counts = adata.layers[unique_layer]
    multi_counts = adata.layers[multi_layer]

    # Sum across all samples for each transcript
    unique_totals = np.sum(unique_counts, axis=0)
    multi_totals = np.sum(multi_counts, axis=0)
    total_counts = unique_totals + multi_totals

    # Calculate per-transcript multimapping ratios
    per_transcript_ratios = np.zeros(len(total_counts), dtype=float)
    non_zero_mask = total_counts > 0
    per_transcript_ratios[non_zero_mask] = multi_totals[non_zero_mask] / total_counts[non_zero_mask]

    # Initialize per-allele ratios
    per_allele_ratios = np.zeros(len(per_transcript_ratios), dtype=float)
    gene_ids = adata.var[gene_grouping_column]
    unique_genes = gene_ids.unique()

    for gene_id in unique_genes:
        # Skip NaN or zero gene IDs
        if pd.isna(gene_id) or gene_id == 0:
            gene_mask = pd.isna(gene_ids) | (gene_ids == 0)
            per_allele_ratios[gene_mask] = per_transcript_ratios[gene_mask]
            continue

        # Get indices of transcripts belonging to this gene
        gene_mask = gene_ids == gene_id
        gene_transcript_indices = np.where(gene_mask)[0]

        # Get ratios and counts for transcripts in this gene
        gene_ratios = per_transcript_ratios[gene_transcript_indices]
        gene_counts = total_counts[gene_transcript_indices]

        if not count_scaling:
            # Original behavior: use maximum ratio
            max_ratio = np.max(gene_ratios)
            per_allele_ratios[gene_transcript_indices] = max_ratio
        else:
            # Apply count-based scaling
            if scaling_method == 'weighted_average':
                # Calculate weighted average ratio, giving more weight to highly expressed transcripts
                if np.sum(gene_counts) > 0:
                    weights = gene_counts / np.sum(gene_counts)
                    weighted_ratio = np.sum(gene_ratios * weights)
                else:
                    weighted_ratio = 0
                per_allele_ratios[gene_transcript_indices] = weighted_ratio

            elif scaling_method == 'max_weighted':
                # Take maximum ratio but reduce it if it comes from low-count transcripts
                max_idx = np.argmax(gene_ratios)
                max_ratio = gene_ratios[max_idx]
                max_counts = gene_counts[max_idx]

                # Apply scaling factor based on count abundance
                if max_counts < min_counts_threshold:
                    # Reduce the ratio for low-count transcripts
                    scaling_factor = max_counts / min_counts_threshold
                    scaled_ratio = max_ratio * scaling_factor
                else:
                    scaled_ratio = max_ratio

                per_allele_ratios[gene_transcript_indices] = scaled_ratio

            elif scaling_method == 'abundance_filtered':
                # Only consider transcripts above the count threshold
                high_count_mask = gene_counts >= min_counts_threshold

                if np.any(high_count_mask):
                    # Use maximum ratio among high-count transcripts
                    filtered_ratios = gene_ratios[high_count_mask]
                    max_ratio = np.max(filtered_ratios)
                else:
                    # If no transcripts meet threshold, use weighted average
                    if np.sum(gene_counts) > 0:
                        weights = gene_counts / np.sum(gene_counts)
                        max_ratio = np.sum(gene_ratios * weights)
                    else:
                        max_ratio = 0

                per_allele_ratios[gene_transcript_indices] = max_ratio

            else:
                raise ValueError(f"Unknown scaling_method: {scaling_method}")

    # Add to var with descriptive column name
    if count_scaling:
        column_name = f'multimapping_ratio_per_allele_{scaling_method}'
    else:
        column_name = 'multimapping_ratio_per_allele'

    adata.var[column_name] = per_allele_ratios

    # Also store some diagnostic information
    adata.var['transcript_total_counts'] = total_counts
    adata.var['transcript_multimapping_ratio'] = per_transcript_ratios

    return adata
