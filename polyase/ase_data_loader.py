import pandas as pd
import anndata as ad
import os
from pathlib import Path
import numpy as np
from scipy import sparse

def load_ase_data(
    var_obs_file,
    isoform_counts_dir,
    tx_to_gene_file,
    sample_info=None,
    counts_file=None,
    fillna=0,
    calculate_cpm=True,
    quant_dir=None
):
    """
    Load allele-specific expression data from long-read RNAseq at isoform level.

    Parameters
    -----------
    var_obs_file : str
        Path to the variant observations file
    isoform_counts_dir : str
        Directory containing the isoform counts files
    tx_to_gene_file : str
        Path to TSV file mapping transcript_id to gene_id
    sample_info : dict, optional
        Dictionary mapping sample IDs to their conditions (e.g., {'SRR14993892': 'leaf', 'SRR14993896': 'tuber'})
        If None, all TSV files in isoform_counts_dir will be used and conditions will be extracted from filenames
    counts_file : str, optional
        Path to additional counts file (salmon merged transcript counts). Optional.
    fillna : int or float, optional
        Value to fill NA values with
    calculate_cpm : bool, optional
        Whether to calculate CPM (Counts Per Million) from EM counts (default: True)
    quant_dir : str, optional
        Directory containing quantification files with EM counts (quant.sf files in {sample_id}/quant.sf structure)
        If None, uses the parent directory of isoform_counts_dir

    Returns
    --------
    anndata.AnnData
        AnnData object containing the processed isoform-level data with EM counts and CPM layers
    """
    # Load variant observations file
    var_obs = pd.read_csv(var_obs_file, delimiter="\t", index_col=0)

    # Load transcript to gene mapping
    tx_to_gene = pd.read_csv(tx_to_gene_file, delimiter="\t")
    # Ensure the mapping has the expected columns
    if 'transcript_id' not in tx_to_gene.columns or 'gene_id' not in tx_to_gene.columns:
        raise ValueError("tx_to_gene_file must contain 'transcript_id' and 'gene_id' columns")
    tx_to_gene_dict = tx_to_gene.set_index('transcript_id')['gene_id'].to_dict()

    # Load additional counts file if provided
    if counts_file:
        additional_counts = pd.read_csv(counts_file, delimiter="\t")
        # You can add code here to use the additional_counts data if needed

    # Find sample files and their conditions if sample_info not provided
    if sample_info is None:
        sample_info = {}
        counts_files = list(Path(isoform_counts_dir).glob("*.counts.tsv"))

        for file_path in counts_files:
            # Extract sample ID and condition from filename (assuming format like "SRR14993892_leaf.counts.tsv")
            filename = file_path.stem  # Gets filename without extension
            parts = filename.split('_')
            sample_id = parts[0]
            condition = parts[1] if len(parts) > 1 else "unknown"
            sample_info[sample_id] = condition

    sample_ids = list(sample_info.keys())

    # Set quant_dir if not provided (assume it's in the parent directory of isoform_counts_dir)
    if quant_dir is None:
        quant_dir = Path(isoform_counts_dir).parent

    # Load ambiguous, unique counts, and EM counts for each sample
    ambig_counts_dfs = []
    unique_counts_dfs = []
    em_counts_dfs = []

    for sample_id in sample_ids:
        condition = sample_info[sample_id]

        # Look for isoform-specific files with multiple naming patterns
        file_patterns = [
            f"{sample_id}_{condition}.isoform.counts.tsv",
            f"{sample_id}_{condition}.transcript.counts.tsv",
            f"{sample_id}_{condition}.tx.counts.tsv",
            f"{sample_id}_{condition}.counts.tsv"
        ]

        file_path = None
        for pattern in file_patterns:
            potential_path = os.path.join(isoform_counts_dir, pattern)
            if os.path.exists(potential_path):
                file_path = potential_path
                break

        # Check if the file exists, otherwise try without condition
        if not file_path:
            # ambig_info.tsv output from salmon/oarfish
            # transcript ids need to be added in the first columns!
            alternate_files = list(Path(isoform_counts_dir).glob(f"{sample_id}*.ambig_info.tsv"))
            if alternate_files:
                file_path = str(alternate_files[0])
            else:
                print(f"Warning: No file found for sample {sample_id}. The counts dir should contain files named like '{sample_id}_condition.counts.tsv' or '{sample_id}*.ambig_info.tsv'.")
                continue

        counts_df = pd.read_csv(file_path, delimiter="\t", index_col=0)

        # Try different possible column names
        ambig_col = None
        unique_col = None

        # Check for AmbigCount/UniqueCount first, then fallback to alternatives
        if 'AmbigCount' in counts_df.columns:
            ambig_col = 'AmbigCount'
        elif 'ambig_reads' in counts_df.columns:
            ambig_col = 'ambig_reads'

        if 'UniqueCount' in counts_df.columns:
            unique_col = 'UniqueCount'
        elif 'unique_reads' in counts_df.columns:
            unique_col = 'unique_reads'

        if ambig_col and unique_col:
            ambig_counts_dfs.append(counts_df[ambig_col])
            unique_counts_dfs.append(counts_df[unique_col])
        else:
            print(f"Warning: Expected columns not found in {file_path}")
            print(f"Available columns: {counts_df.columns.tolist()}")

        # Load EM counts from quant.sf file
        quant_file_path = os.path.join(quant_dir, sample_id, "quant.sf")
        if os.path.exists(quant_file_path):
            try:
                # Read the quant.sf file (tab-separated with header)
                # Note: quant.sf uses 'Name' as the first column, not index
                em_df = pd.read_csv(quant_file_path, delimiter="\t")

                # Check for expected columns
                if 'tname' not in em_df.columns:
                    raise ValueError(f"Expected 'tname' column not found in {quant_file_path}")
                if 'num_reads' not in em_df.columns:
                    raise ValueError(f"Expected 'num_reads' column not found in {quant_file_path}")

                # Set transcript name as index and extract num_reads
                em_df = em_df.set_index('tname')
                em_counts_dfs.append(em_df['num_reads'])

            except Exception as e:
                print(f"Error reading EM counts from {quant_file_path}: {str(e)}")
                # Add empty series if file can't be read
                em_counts_dfs.append(pd.Series(dtype=float, name='num_reads'))
        else:
            print(f"Warning: EM counts file not found at {quant_file_path}")
            # Add empty series if file doesn't exist
            em_counts_dfs.append(pd.Series(dtype=float, name='num_reads'))

    # Concatenate ambiguous counts across samples
    ambig_counts = pd.concat(ambig_counts_dfs, axis=1)
    ambig_counts = ambig_counts.loc[~ambig_counts.index.duplicated(keep='first')]

    # Concatenate unique counts across samples
    unique_counts = pd.concat(unique_counts_dfs, axis=1)
    unique_counts = unique_counts.loc[~unique_counts.index.duplicated(keep='first')]

    # Concatenate EM counts across samples
    em_counts = pd.concat(em_counts_dfs, axis=1)
    em_counts = em_counts.loc[~em_counts.index.duplicated(keep='first')]

    # Rename columns
    ambig_counts.columns = sample_ids
    unique_counts.columns = sample_ids
    em_counts.columns = sample_ids

    # Create isoform metadata
    # Use the union of all transcript IDs from all count matrices
    all_transcript_ids = set(unique_counts.index) | set(ambig_counts.index) | set(em_counts.index)
    isoform_var = pd.DataFrame(index=list(all_transcript_ids))
    isoform_var['transcript_id'] = isoform_var.index
    isoform_var['gene_id'] = isoform_var['transcript_id'].map(tx_to_gene_dict)
    isoform_var['feature_type'] = 'transcript'

    # Match var_obs data based on gene_id
    # First, get gene_ids that exist in both var_obs and isoform data
    var_obs_genes = set(var_obs.index)
    isoform_genes = set(isoform_var['gene_id'].dropna())
    matching_genes = var_obs_genes.intersection(isoform_genes)

    print(f"Found {len(matching_genes)} genes with matching annotations")

    if len(matching_genes) > 0:
        # Add var_obs columns to isoform_var, matching by gene_id
        for col in var_obs.columns:
            # Determine the appropriate dtype for this column
            col_dtype = var_obs[col].dtype
            if col_dtype == 'object' or pd.api.types.is_string_dtype(col_dtype):
                # For string/object columns, use None instead of np.nan
                isoform_var[col] = None
                isoform_var[col] = isoform_var[col].astype('object')
            else:
                # For numeric columns, use np.nan
                isoform_var[col] = np.nan

            # For each transcript, use the var_obs data from its parent gene
            for transcript_id in isoform_var.index:
                gene_id = isoform_var.loc[transcript_id, 'gene_id']
                if pd.notna(gene_id) and gene_id in var_obs.index:
                    isoform_var.loc[transcript_id, col] = var_obs.loc[gene_id, col]

    # Filter to only keep transcripts that have matching gene_ids in var_obs
    transcripts_with_gene_match = isoform_var['gene_id'].isin(var_obs.index)
    isoform_var = isoform_var[transcripts_with_gene_match]

    print(f"Keeping {len(isoform_var)} transcripts with gene matches")

    # Reindex counts to match filtered isoform_var (only transcripts with gene matches)
    ambig_counts = ambig_counts.reindex(isoform_var.index, fill_value=fillna)
    unique_counts = unique_counts.reindex(isoform_var.index, fill_value=fillna)
    em_counts = em_counts.reindex(isoform_var.index, fill_value=fillna)

    # Fill NA values
    ambig_counts.fillna(fillna, inplace=True)
    unique_counts.fillna(fillna, inplace=True)
    em_counts.fillna(fillna, inplace=True)

    # Create AnnData object (use EM counts as main X matrix if available, otherwise unique counts)
    if not em_counts.empty and em_counts.sum().sum() > 0:
        main_counts = em_counts
        print("Using EM counts as main expression matrix (X)")
    else:
        main_counts = unique_counts
        print("Using unique counts as main expression matrix (X) - no EM counts available")

    adata = ad.AnnData(
        X=main_counts[sample_ids].T,
        var=isoform_var,
        obs=pd.DataFrame(index=sample_ids),
    )

    # Add conditions to obs
    conditions = [sample_info[sample_id] for sample_id in sample_ids]
    adata.obs["condition"] = conditions

    # Add sample IDs to obs_names
    adata.obs_names = sample_ids

    # Add layers of unique, ambiguous, and EM counts
    adata.layers["unique_counts"] = unique_counts.T.to_numpy()
    adata.layers["ambiguous_counts"] = ambig_counts.T.to_numpy()
    adata.layers["em_counts"] = em_counts.T.to_numpy()

    # Calculate CPM if requested (using EM counts as the primary source)
    if calculate_cpm:
        # Calculate library sizes from EM counts (total EM counts per sample)
        em_lib_sizes = np.sum(adata.layers["em_counts"], axis=1)
        unique_lib_sizes = np.sum(adata.layers["unique_counts"], axis=1)
        ambig_lib_sizes = np.sum(adata.layers["ambiguous_counts"], axis=1)

        # Store library sizes in obs
        adata.obs["em_lib_size"] = em_lib_sizes
        adata.obs["unique_lib_size"] = unique_lib_sizes
        adata.obs["ambig_lib_size"] = ambig_lib_sizes
        adata.obs["total_lib_size"] = unique_lib_sizes + ambig_lib_sizes

        # Calculate CPM for EM counts (primary CPM calculation)
        em_cpm = np.zeros_like(adata.layers["em_counts"], dtype=float)
        for i in range(len(sample_ids)):
            if em_lib_sizes[i] > 0:
                em_cpm[i, :] = (adata.layers["em_counts"][i, :] / em_lib_sizes[i]) * 1e6

        # Calculate CPM for unique counts (based on EM library size for consistency)
        unique_cpm = np.zeros_like(adata.layers["unique_counts"], dtype=float)
        for i in range(len(sample_ids)):
            if em_lib_sizes[i] > 0:
                unique_cpm[i, :] = (adata.layers["unique_counts"][i, :] / em_lib_sizes[i]) * 1e6

        # Calculate CPM for ambiguous counts (based on EM library size for consistency)
        ambig_cpm = np.zeros_like(adata.layers["ambiguous_counts"], dtype=float)
        for i in range(len(sample_ids)):
            if em_lib_sizes[i] > 0:
                ambig_cpm[i, :] = (adata.layers["ambiguous_counts"][i, :] / em_lib_sizes[i]) * 1e6

        # Add CPM layers
        adata.layers["em_cpm"] = em_cpm
        adata.layers["unique_cpm"] = unique_cpm
        adata.layers["ambiguous_cpm"] = ambig_cpm

        print(f"Added CPM layers (normalized by EM count library sizes):")
        print(f"  - em_cpm: EM counts normalized to CPM")
        print(f"  - unique_cpm: Unique counts normalized to CPM (using EM lib size)")
        print(f"  - ambiguous_cpm: Ambiguous counts normalized to CPM (using EM lib size)")

        # Print library size statistics
        print(f"\nLibrary size statistics:")
        print(f"  - EM counts mean lib size: {em_lib_sizes.mean():.0f}")
        print(f"  - Unique counts mean lib size: {unique_lib_sizes.mean():.0f}")
        print(f"  - Ambiguous counts mean lib size: {ambig_lib_sizes.mean():.0f}")

    return adata


def aggregate_transcripts_to_genes(adata_tx):
    """
    Aggregate transcript-level AnnData to gene-level AnnData.
    Optimized version using vectorized operations and sparse matrices.

    Parameters
    ----------
    adata_tx : AnnData
        Transcript-level AnnData object

    Returns
    -------
    AnnData
        Gene-level AnnData object
    """
    # Get unique genes and create mapping
    gene_ids = adata_tx.var['gene_id'].dropna()
    unique_genes = gene_ids.unique()
    n_genes = len(unique_genes)
    n_obs = adata_tx.n_obs

    print(f"Aggregating {adata_tx.n_vars} transcripts to {n_genes} genes")

    # Create mapping from gene_id to index
    gene_to_idx = {gene: idx for idx, gene in enumerate(unique_genes)}

    # Create transcript to gene mapping matrix (sparse for efficiency)
    # This matrix has shape (n_transcripts, n_genes)
    tx_indices = []
    gene_indices = []

    for tx_idx, gene_id in enumerate(gene_ids):
        if pd.notna(gene_id) and gene_id in gene_to_idx:
            tx_indices.append(tx_idx)
            gene_indices.append(gene_to_idx[gene_id])

    # Create sparse mapping matrix
    mapping_matrix = sparse.csr_matrix(
        (np.ones(len(tx_indices)), (tx_indices, gene_indices)),
        shape=(len(gene_ids), n_genes)
    )

    # Vectorized aggregation
    # Convert to sparse if not already
    X_sparse = sparse.csr_matrix(adata_tx.X) if not sparse.issparse(adata_tx.X) else adata_tx.X
    unique_counts_sparse = sparse.csr_matrix(adata_tx.layers['unique_counts']) if not sparse.issparse(adata_tx.layers['unique_counts']) else adata_tx.layers['unique_counts']
    ambiguous_counts_sparse = sparse.csr_matrix(adata_tx.layers['ambiguous_counts']) if not sparse.issparse(adata_tx.layers['ambiguous_counts']) else adata_tx.layers['ambiguous_counts']

    # Aggregate using matrix multiplication
    X_gene = X_sparse @ mapping_matrix
    unique_counts_gene = unique_counts_sparse @ mapping_matrix

    # For ambiguous counts (mean), we need to divide by the number of transcripts per gene
    ambiguous_sum = ambiguous_counts_sparse @ mapping_matrix
    transcripts_per_gene = mapping_matrix.sum(axis=0).A1  # Convert to 1D array
    transcripts_per_gene[transcripts_per_gene == 0] = 1  # Avoid division by zero
    ambiguous_counts_gene = ambiguous_sum / transcripts_per_gene[np.newaxis, :]

    # Handle EM counts if present
    if 'em_counts' in adata_tx.layers:
        em_counts_sparse = sparse.csr_matrix(adata_tx.layers['em_counts']) if not sparse.issparse(adata_tx.layers['em_counts']) else adata_tx.layers['em_counts']
        em_counts_gene = em_counts_sparse @ mapping_matrix
        em_counts_gene = em_counts_gene.toarray() if sparse.issparse(em_counts_gene) else em_counts_gene

    # Convert back to dense arrays if needed
    X_gene = X_gene.toarray() if sparse.issparse(X_gene) else X_gene
    unique_counts_gene = unique_counts_gene.toarray() if sparse.issparse(unique_counts_gene) else unique_counts_gene
    ambiguous_counts_gene = ambiguous_counts_gene.toarray() if sparse.issparse(ambiguous_counts_gene) else ambiguous_counts_gene

    # Create gene-level var DataFrame
    gene_var = pd.DataFrame(index=unique_genes)
    gene_var['gene_id'] = unique_genes
    gene_var['feature_type'] = 'gene'

    # Aggregate metadata efficiently using groupby
    gene_metadata_cols = [
        'transcript_id',
        'Synt_id', 'synteny_category', 'syntenic_genes', 'haplotype',
        'CDS_length_category', 'CDS_percent_difference',
        'CDS_haplotype_with_longest_annotation'
    ]

    # Filter to only transcripts with valid gene_ids
    valid_tx_mask = gene_ids.notna()
    tx_var_valid = adata_tx.var[valid_tx_mask].copy()

    for col in gene_metadata_cols:
        if col in adata_tx.var.columns:
            # For other metadata, take the first non-null value for each gene
            gene_metadata = tx_var_valid.groupby('gene_id')[col].first()
            gene_var[col] = gene_metadata.reindex(unique_genes)

    # Create gene-level AnnData
    adata_gene = ad.AnnData(
        X=X_gene,
        obs=adata_tx.obs.copy(),
        var=gene_var
    )

    # Add layers
    adata_gene.layers['unique_counts'] = unique_counts_gene
    adata_gene.layers['ambiguous_counts'] = ambiguous_counts_gene

    # Add EM counts layer if present in original data
    if 'em_counts' in adata_tx.layers:
        adata_gene.layers['em_counts'] = em_counts_gene

    # Aggregate CPM layers if they exist
    cpm_layers = ['em_cpm', 'unique_cpm', 'ambiguous_cpm']
    for layer_name in cpm_layers:
        if layer_name in adata_tx.layers:
            layer_sparse = sparse.csr_matrix(adata_tx.layers[layer_name]) if not sparse.issparse(adata_tx.layers[layer_name]) else adata_tx.layers[layer_name]
            layer_gene = layer_sparse @ mapping_matrix
            adata_gene.layers[layer_name] = layer_gene.toarray() if sparse.issparse(layer_gene) else layer_gene

    # Add summary statistics
    n_transcripts_per_gene = gene_ids.value_counts()
    adata_gene.var['n_transcripts'] = n_transcripts_per_gene.reindex(
        adata_gene.var_names, fill_value=0
    )

    print(f"Created gene-level AnnData: {adata_gene.n_obs} × {adata_gene.n_vars}")
    print(f"Average transcripts per gene: {adata_gene.var['n_transcripts'].mean():.2f}")

    return adata_gene
