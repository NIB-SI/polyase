"""
AnnData Filtering Module
------------------------
This module provides functions for filtering AnnData objects based on group expression patterns.
It supports filtering by group expression levels with various normalization methods and thresholds.
"""

import numpy as np
import pandas as pd
import scipy.sparse
import anndata as ad
from typing import Dict, List, Optional, Set, Tuple, Union, Literal, Callable


def _get_group_mapping(
    adata: ad.AnnData,
    group_col: str,
    group_source: str
) -> pd.Series:
    """
    Extract group mapping from AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing group information
    group_col : str or tuple
        Column name containing group IDs. For obsm/varm, should be a tuple (key, column_index)
    group_source : str
        Location of the group column in AnnData ('obs', 'var', 'obsm', 'varm')

    Returns
    -------
    pandas.Series
        Series mapping indices to their group IDs

    Raises
    ------
    ValueError
        If the specified group column or source is not found in the AnnData object
    """
    if group_source == 'obs':
        if group_col not in adata.obs:
            raise ValueError(f"Group column '{group_col}' not found in adata.obs")
        return adata.obs[group_col]

    elif group_source == 'var':
        if group_col not in adata.var:
            raise ValueError(f"Group column '{group_col}' not found in adata.var")
        return adata.var[group_col]

    elif group_source == 'obsm':
        if isinstance(group_col, tuple) and len(group_col) == 2:
            obsm_key, col_idx = group_col
            if obsm_key not in adata.obsm:
                raise ValueError(f"obsm key '{obsm_key}' not found")
            return pd.Series(adata.obsm[obsm_key][:, col_idx], index=adata.obs_names)
        else:
            raise ValueError("For group_source='obsm', group_col must be a tuple (key, column_index)")

    elif group_source == 'varm':
        if isinstance(group_col, tuple) and len(group_col) == 2:
            varm_key, col_idx = group_col
            if varm_key not in adata.varm:
                raise ValueError(f"varm key '{varm_key}' not found")
            return pd.Series(adata.varm[varm_key][:, col_idx], index=adata.var_names)
        else:
            raise ValueError("For group_source='varm', group_col must be a tuple (key, column_index)")

    else:
        raise ValueError("group_source must be one of: 'obs', 'var', 'obsm', 'varm'")


def filter_low_expressed_genes(
    adata: ad.AnnData,
    min_expression: Union[float, Dict[str, float], Callable[[float], float]] = 1.0,
    library_size_dependent: bool = False,
    lib_size_normalization: Optional[str] = 'cpm',  # Options: 'cpm', 'rpkm', 'tpm', None
    layer: Optional[str] = None,
    group_col: Union[str, Tuple[str, int]] = 'Synt_id',
    group_source: str = 'var',
    mode: str = 'any',
    return_dropped: bool = False,
    copy: bool = True,
    filter_axis: Literal[0, 1] = 1,  # 1 = filter rows, 0 = filter columns
    verbose: bool = True
) -> Union[ad.AnnData, Tuple[ad.AnnData, List[str]]]:
    """
    Filter an AnnData object to remove groups with low expression across samples.

    Parameters
    ----------
    adata : AnnData
        AnnData object with group IDs and expression data
    min_expression : float, dict, or callable, default=1.0
        Minimum expression threshold for groups:
        - float: same threshold applied to all samples/features
        - dict: {name: threshold} for sample/feature-specific thresholds
        - callable: function that takes library size and returns threshold
          (e.g., lambda lib_size: lib_size * 1e-6 for 0.0001% of lib size)
    library_size_dependent : bool, default=False
        If True, scale thresholds by library size for each sample
    lib_size_normalization : str or None, default='cpm'
        How to normalize for library size:
        - 'cpm': Counts Per Million (divide by lib_size/1e6)
        - 'rpkm': Reads Per Kilobase Million (not implemented yet)
        - 'tpm': Transcripts Per Million (not implemented yet)
        - None: No normalization
    layer : str or None, default=None
        Layer to use for expression values. If None, use .X
    group_col : str or tuple, default='Synt_id'
        Column name containing group IDs. For obsm/varm, should be a tuple (key, column_index)
    group_source : str, default='var'
        Location of the group column in AnnData ('obs', 'var', 'obsm', 'varm')
    mode : str, default='any'
        'any': Keep groups that pass threshold in any sample/feature
        'all': Keep groups that pass threshold in all samples/features
        'mean': Keep groups that pass threshold on average across samples/features
    return_dropped : bool, default=False
        If True, also return list of dropped group IDs
    copy : bool, default=True
        If True, return a copy of the filtered AnnData object
        If False, filter the AnnData object in place
    filter_axis : int, default=1
        1: Filter rows (obs) based on group expression across columns (var)
        0: Filter columns (var) based on group expression across rows (obs)
    verbose : bool, default=True
        Whether to print additional information during filtering

    Returns
    -------
    AnnData or tuple
        Filtered AnnData object, and optionally a list of dropped group IDs

    Raises
    ------
    ValueError
        If parameters are invalid or required data is missing
    """
    # Check parameters
    if mode not in ['any', 'all', 'mean']:
        raise ValueError("mode must be one of 'any', 'all', or 'mean'")

    if lib_size_normalization not in ['cpm', None]:
        if lib_size_normalization in ['rpkm', 'tpm']:
            raise NotImplementedError(f"Normalization method '{lib_size_normalization}' not yet implemented")
        else:
            raise ValueError(f"Unsupported normalization method: {lib_size_normalization}, use 'cpm' or None")

    # Get group mapping
    group_mapping = _get_group_mapping(adata, group_col, group_source)

    # Get expression matrix from layer or .X
    if layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"Layer {layer} not found in AnnData object")
        expr_matrix = adata.layers[layer]
    else:
        expr_matrix = adata.X

    # Convert to dense if sparse
    if scipy.sparse.issparse(expr_matrix):
        expr_matrix = expr_matrix.toarray()

    # Calculate library sizes if needed
    if library_size_dependent or lib_size_normalization:
        if filter_axis == 1:  # When filtering rows
            lib_sizes = np.sum(expr_matrix, axis=1)
            adata.obs['lib_size'] = lib_sizes
        else:  # When filtering columns
            lib_sizes = np.sum(expr_matrix, axis=0)
            adata.var['lib_size'] = lib_sizes

    # Create a dataframe with expression data
    if filter_axis == 0:  # Filter rows based on group expression
        expr_df = pd.DataFrame(
            expr_matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )
        # Add group column
        expr_df['group'] = group_mapping

        # Determine which dimension to sum over for thresholds
        threshold_dim = list(adata.var_names)

        # Get the library sizes for normalization if needed
        if library_size_dependent:
            sample_lib_sizes = adata.var.get('lib_size', None)
            if sample_lib_sizes is None:
                sample_lib_sizes = {var: 1.0 for var in threshold_dim}
            else:
                sample_lib_sizes = sample_lib_sizes.to_dict()

    elif filter_axis == 1:  # Filter columns based on group expression
        expr_df = pd.DataFrame(
            expr_matrix.T,  # Transpose for column filtering
            index=adata.var_names,
            columns=adata.obs_names
        )
        # Add group column
        expr_df['group'] = group_mapping

        # Determine which dimension to sum over for thresholds
        threshold_dim = list(adata.obs_names)

        # Get the library sizes for normalization if needed
        if library_size_dependent:
            sample_lib_sizes = adata.obs.get('lib_size', None)
            if sample_lib_sizes is None:
                sample_lib_sizes = {obs: 1.0 for obs in threshold_dim}
            else:
                sample_lib_sizes = sample_lib_sizes.to_dict()

    else:
        raise ValueError("filter_axis must be 1 (filter rows) or 0 (filter columns)")

    # Apply normalization if specified
    if lib_size_normalization:
        norm_factor = {}

        if lib_size_normalization == 'cpm':
            for item in threshold_dim:
                if filter_axis == 1:
                    norm_factor[item] = adata.obs.loc[item, 'lib_size'] / 1e6
                else:
                    norm_factor[item] = adata.var.loc[item, 'lib_size'] / 1e6

        # Apply normalization to the data
        for item in threshold_dim:
            if item in norm_factor and norm_factor[item] > 0:
                expr_df[item] = expr_df[item] / norm_factor[item]

    # Group by group ID and sum expression
    grouped_expr = expr_df.groupby('group').sum()

    # Convert threshold to dictionary if it's a scalar
    thresholds = {}

    if callable(min_expression) and library_size_dependent:
        # Use the function to calculate thresholds based on library size
        for item in threshold_dim:
            if filter_axis == 1:
                lib_size = adata.obs.loc[item, 'lib_size']
            else:
                lib_size = adata.var.loc[item, 'lib_size']
            thresholds[item] = min_expression(lib_size)

    elif isinstance(min_expression, (int, float)):
        if library_size_dependent:
            # Scale the threshold by library size
            for item in threshold_dim:
                if filter_axis == 1:
                    scaling = adata.obs.loc[item, 'lib_size'] / 1e6  # CPM scaling
                else:
                    scaling = adata.var.loc[item, 'lib_size'] / 1e6  # CPM scaling
                thresholds[item] = min_expression * scaling
        else:
            # Same threshold for all
            thresholds = {item: min_expression for item in threshold_dim}
    else:
        # Dictionary of thresholds provided
        thresholds = min_expression

    if verbose:
        if library_size_dependent:
            print("Using library-size adjusted thresholds:")
            for k, v in list(thresholds.items())[:5]:
                print(f"  {k}: {v:.2f}")
            if len(thresholds) > 5:
                print(f"  ... and {len(thresholds)-5} more")
        else:
            print(f"Using uniform threshold: {min_expression}")

    # Determine which groups to keep based on mode
    if mode == 'any':
        # Keep groups that pass threshold in any sample/feature
        keep_groups = grouped_expr.apply(lambda row: any(row[item] >= thresholds[item]
                                            for item in threshold_dim), axis=1)
    elif mode == 'all':
        # Keep groups that pass threshold in all samples/features
        keep_groups = grouped_expr.apply(lambda row: all(row[item] >= thresholds[item]
                                            for item in threshold_dim), axis=1)
    elif mode == 'mean':
        # Keep groups that pass threshold on average
        avg_expr = grouped_expr.mean(axis=1)
        avg_threshold = np.mean(list(thresholds.values()))
        keep_groups = avg_expr >= avg_threshold

    # Get group IDs to keep
    keep_group_ids = keep_groups[keep_groups].index.tolist()

    # Store dropped IDs for reference
    dropped_group_ids = list(set(group_mapping.unique()) - set(keep_group_ids))

    # Filter the AnnData object
    if filter_axis == 0:  # Filter rows
        keep_indices = group_mapping.isin(keep_group_ids)
    else:  # Filter columns
        keep_indices = group_mapping.isin(keep_group_ids)

    if verbose:
        print(f"Filtered out {len(dropped_group_ids)} groups")
        print(f"Kept {keep_indices.sum()} / {len(keep_indices)} items")

    if copy:
        if filter_axis == 0:
            filtered_adata = adata[keep_indices].copy()
        else:
            filtered_adata = adata[:, keep_indices].copy()
    else:
        # Filter in place
        if filter_axis == 0:
            adata._inplace_subset_obs(keep_indices)
        else:
            adata._inplace_subset_var(keep_indices)
        filtered_adata = adata

    if return_dropped:
        return filtered_adata, dropped_group_ids

    return filtered_adata
