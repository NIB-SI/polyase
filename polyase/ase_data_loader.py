import pandas as pd
import anndata as ad
import os
from pathlib import Path
import numpy as np
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings

def _load_sample_counts(sample_id, condition, isoform_counts_dir, quant_dir, fillna=0):
    """
    Load counts for a single sample - designed for parallel execution.
    """
    result = {
        'sample_id': sample_id,
        'ambig_counts': None,
        'unique_counts': None,
        'em_counts': None
    }

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
        alternate_files = list(Path(isoform_counts_dir).glob(f"{sample_id}*.ambig_info.tsv"))
        if alternate_files:
            file_path = str(alternate_files[0])
        else:
            print(f"Warning: No file found for sample {sample_id}")
            return result

    try:
        counts_df = pd.read_csv(file_path, delimiter="\t", index_col=0)

        # Try different possible column names
        ambig_col = None
        unique_col = None

        if 'AmbigCount' in counts_df.columns:
            ambig_col = 'AmbigCount'
        elif 'ambig_reads' in counts_df.columns:
            ambig_col = 'ambig_reads'

        if 'UniqueCount' in counts_df.columns:
            unique_col = 'UniqueCount'
        elif 'unique_reads' in counts_df.columns:
            unique_col = 'unique_reads'

        if ambig_col and unique_col:
            result['ambig_counts'] = counts_df[ambig_col]
            result['unique_counts'] = counts_df[unique_col]
        else:
            print(f"Warning: Expected columns not found in {file_path}")
            return result
    except Exception as e:
        print(f"Error reading counts from {file_path}: {str(e)}")
        return result

    # Load EM counts from quant.sf file
    quant_file_path = os.path.join(quant_dir, sample_id, "quant.sf")
    if os.path.exists(quant_file_path):
        try:
            em_df = pd.read_csv(quant_file_path, delimiter="\t")

            if 'tname' in em_df.columns and 'num_reads' in em_df.columns:
                print("Oarfish quant.sf format detected")
                em_df = em_df.set_index('tname')
                result['em_counts'] = em_df['num_reads']
            elif 'Name' in em_df.columns and 'NumReads' in em_df.columns:
                print("Salmon quant.sf format detected")
                em_df = em_df.set_index('Name')
                result['em_counts'] = em_df['NumReads']
            else:
                print(f"Warning: Expected columns not found in {quant_file_path}")
        except Exception as e:
            print(f"Error reading EM counts from {quant_file_path}: {str(e)}")
    else:
        print(f"Warning: EM counts file not found at {quant_file_path}")

    return result


def load_ase_data(
    var_obs_file,
    isoform_counts_dir,
    tx_to_gene_file,
    sample_info=None,
    counts_file=None,
    fillna=0,
    calculate_cpm=True,
    quant_dir=None,
    n_jobs=4  # Number of parallel jobs for loading samples
):
    """
    Load allele-specific expression data from long-read RNAseq at isoform level.
    Optimized version with parallel loading and vectorized operations.

    Parameters
    -----------
    var_obs_file : str
        Path to the variant observations file
    isoform_counts_dir : str
        Directory containing the isoform counts files
    tx_to_gene_file : str
        Path to TSV file mapping transcript_id to gene_id
    sample_info : dict, optional
        Dictionary mapping sample IDs to their conditions
    counts_file : str, optional
        Path to additional counts file (salmon merged transcript counts). Optional.
    fillna : int or float, optional
        Value to fill NA values with
    calculate_cpm : bool, optional
        Whether to calculate CPM (Counts Per Million) from EM counts (default: True)
    quant_dir : str, optional
        Directory containing quantification files with EM counts
    n_jobs : int, optional
        Number of parallel jobs for loading samples (default: 4)

    Returns
    --------
    anndata.AnnData
        AnnData object containing the processed isoform-level data with EM counts and CPM layers.
        Includes all transcripts from expression matrix and tx2gene mapping, with NaN values
        for var_obs data when genes are not found in var_obs_file.
    """
    print("Loading metadata files...")

    # Load variant observations file
    var_obs = pd.read_csv(var_obs_file, delimiter="\t", index_col=0)

    # Load transcript to gene mapping
    tx_to_gene = pd.read_csv(tx_to_gene_file, delimiter="\t")
    if 'transcript_id' not in tx_to_gene.columns or 'gene_id' not in tx_to_gene.columns:
        raise ValueError("tx_to_gene_file must contain 'transcript_id' and 'gene_id' columns")
    tx_to_gene_dict = tx_to_gene.set_index('transcript_id')['gene_id'].to_dict()

    # Find sample files and their conditions if sample_info not provided
    if sample_info is None:
        sample_info = {}
        counts_files = list(Path(isoform_counts_dir).glob("*.counts.tsv"))
        for file_path in counts_files:
            filename = file_path.stem
            parts = filename.split('_')
            sample_id = parts[0]
            condition = parts[1] if len(parts) > 1 else "unknown"
            sample_info[sample_id] = condition

    sample_ids = list(sample_info.keys())

    # Set quant_dir if not provided
    if quant_dir is None:
        quant_dir = Path(isoform_counts_dir).parent

    print(f"Loading counts for {len(sample_ids)} samples in parallel...")

    # Load samples in parallel
    sample_results = {}
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        future_to_sample = {
            executor.submit(_load_sample_counts, sample_id, sample_info[sample_id],
                          isoform_counts_dir, quant_dir, fillna): sample_id
            for sample_id in sample_ids
        }

        for future in as_completed(future_to_sample):
            result = future.result()
            if result['unique_counts'] is not None:
                sample_results[result['sample_id']] = result

    # Filter sample_ids to only those successfully loaded
    successful_samples = list(sample_results.keys())
    print(f"Successfully loaded {len(successful_samples)} out of {len(sample_ids)} samples")

    if not successful_samples:
        raise ValueError("No samples were successfully loaded")

    print("Concatenating count matrices...")

    # Prepare data for concatenation
    ambig_data = {}
    unique_data = {}
    em_data = {}

    for sample_id in successful_samples:
        result = sample_results[sample_id]
        if result['ambig_counts'] is not None:
            ambig_data[sample_id] = result['ambig_counts']
        if result['unique_counts'] is not None:
            unique_data[sample_id] = result['unique_counts']
        if result['em_counts'] is not None:
            em_data[sample_id] = result['em_counts']

    # Concatenate count matrices efficiently
    ambig_counts = pd.DataFrame(ambig_data) if ambig_data else pd.DataFrame()
    unique_counts = pd.DataFrame(unique_data) if unique_data else pd.DataFrame()
    em_counts = pd.DataFrame(em_data) if em_data else pd.DataFrame()

    # Get all transcript IDs from expression matrices
    all_transcript_ids = set()
    for df in [ambig_counts, unique_counts, em_counts]:
        if not df.empty:
            all_transcript_ids.update(df.index)

    print("Creating isoform metadata...")

    # Create isoform metadata for ALL transcripts found in expression data
    isoform_var = pd.DataFrame(index=list(all_transcript_ids))
    isoform_var['transcript_id'] = isoform_var.index
    isoform_var['gene_id'] = isoform_var['transcript_id'].map(tx_to_gene_dict)
    isoform_var['feature_type'] = 'transcript'

    # Add var_obs data for ALL transcripts, using NaN for missing genes
    print(f"Adding var_obs annotations for {len(all_transcript_ids)} transcripts...")

    var_obs_genes = set(var_obs.index)
    isoform_genes = set(isoform_var['gene_id'].dropna())
    matching_genes = var_obs_genes.intersection(isoform_genes)
    missing_genes = isoform_genes - var_obs_genes

    print(f"Found {len(matching_genes)} genes with var_obs annotations")
    print(f"Found {len(missing_genes)} genes without var_obs annotations (will be filled with NaN)\n e.g {list(missing_genes)[:3]}\
           {list(missing_genes)[-3:]}.... Could be that they are non-coding or novel..")

    # Create a mapping for efficient lookup
    gene_to_var_obs = var_obs.to_dict('index')

    # Add all var_obs columns to isoform_var, filling with NaN for missing genes
    for col in var_obs.columns:
        col_dtype = var_obs[col].dtype
        if col_dtype == 'object' or pd.api.types.is_string_dtype(col_dtype):
            default_value = None
        else:
            default_value = np.nan

        # Vectorized mapping using gene_id - keeps ALL transcripts
        isoform_var[col] = isoform_var['gene_id'].map(
            lambda gene_id: gene_to_var_obs.get(gene_id, {}).get(col, default_value)
            if pd.notna(gene_id) else default_value
        )

    # Keep ALL transcripts (no filtering based on var_obs matches)
    print(f"Keeping all {len(isoform_var)} transcripts from expression data")

    # Report statistics
    transcripts_with_gene_annotation = isoform_var['gene_id'].notna().sum()
    transcripts_with_var_obs = isoform_var['gene_id'].isin(var_obs.index).sum()

    print(f"Statistics:")
    print(f"  - Total transcripts: {len(isoform_var)}")
    print(f"  - Transcripts with gene_id annotation: {transcripts_with_gene_annotation}")
    print(f"  - Transcripts with var_obs data: {transcripts_with_var_obs}")
    print(f"  - Transcripts with NaN var_obs data: {len(isoform_var) - transcripts_with_var_obs}")

    # Reindex count matrices efficiently
    print("Aligning count matrices...")
    transcript_index = isoform_var.index

    ambig_counts = ambig_counts.reindex(transcript_index, fill_value=fillna)
    unique_counts = unique_counts.reindex(transcript_index, fill_value=fillna)
    em_counts = em_counts.reindex(transcript_index, fill_value=fillna)

    # Ensure all samples are present in all matrices
    for sample_id in successful_samples:
        if sample_id not in ambig_counts.columns:
            ambig_counts[sample_id] = fillna
        if sample_id not in unique_counts.columns:
            unique_counts[sample_id] = fillna
        if sample_id not in em_counts.columns:
            em_counts[sample_id] = fillna

    # Fill remaining NA values
    ambig_counts.fillna(fillna, inplace=True)
    unique_counts.fillna(fillna, inplace=True)
    em_counts.fillna(fillna, inplace=True)

    print("Creating AnnData object...")

    main_counts = unique_counts[successful_samples]
    print("Using unique counts as main expression matrix (X)")

    # Create AnnData object
    adata = ad.AnnData(
        X=main_counts.T.values,  # Use .values for faster conversion
        var=isoform_var,
        obs=pd.DataFrame(index=successful_samples),
    )

    # Add conditions to obs
    conditions = [sample_info[sample_id] for sample_id in successful_samples]
    adata.obs["condition"] = conditions
    adata.obs_names = successful_samples

    # Add layers efficiently using numpy arrays
    adata.layers["unique_counts"] = unique_counts[successful_samples].T.values
    adata.layers["ambiguous_counts"] = ambig_counts[successful_samples].T.values
    adata.layers["em_counts"] = em_counts[successful_samples].T.values

    # Calculate CPM if requested
    if calculate_cpm:
        print("Calculating CPM...")

        # Vectorized library size calculations
        em_lib_sizes = np.sum(adata.layers["em_counts"], axis=1)
        unique_lib_sizes = np.sum(adata.layers["unique_counts"], axis=1)
        ambig_lib_sizes = np.sum(adata.layers["ambiguous_counts"], axis=1)

        # Store library sizes
        adata.obs["em_lib_size"] = em_lib_sizes
        adata.obs["unique_lib_size"] = unique_lib_sizes
        adata.obs["ambig_lib_size"] = ambig_lib_sizes
        adata.obs["total_lib_size"] = unique_lib_sizes + ambig_lib_sizes

        # Vectorized CPM calculations
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='divide by zero')

            # Calculate CPM using broadcasting
            em_lib_sizes_expanded = em_lib_sizes[:, np.newaxis]
            em_lib_sizes_expanded[em_lib_sizes_expanded == 0] = 1  # Avoid division by zero

            adata.layers["em_cpm"] = (adata.layers["em_counts"] / em_lib_sizes_expanded) * 1e6
            adata.layers["unique_cpm"] = (adata.layers["unique_counts"] / em_lib_sizes_expanded) * 1e6
            adata.layers["ambiguous_cpm"] = (adata.layers["ambiguous_counts"] / em_lib_sizes_expanded) * 1e6

        print(f"Added CPM layers (normalized by EM count library sizes)")
        print(f"Library size statistics:")
        print(f"  - EM counts mean lib size: {em_lib_sizes.mean():.0f}")
        print(f"  - Unique counts mean lib size: {unique_lib_sizes.mean():.0f}")
        print(f"  - Ambiguous counts mean lib size: {ambig_lib_sizes.mean():.0f}")

    print("Data loading completed!")
    return adata

def aggregate_transcripts_to_genes(adata_tx):
    """
    Aggregate transcript-level AnnData to gene-level AnnData.
    Optimized version using vectorized operations and sparse matrices.
    """
    # Get unique genes and create mapping
    gene_ids = adata_tx.var['gene_id'].dropna()
    unique_genes = gene_ids.unique()
    n_genes = len(unique_genes)

    print(f"Aggregating {adata_tx.n_vars} transcripts to {n_genes} genes")

    # Create mapping from gene_id to index
    gene_to_idx = {gene: idx for idx, gene in enumerate(unique_genes)}

    # Create transcript to gene mapping matrix (sparse for efficiency)
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

    # Vectorized aggregation for all layers
    layers_to_aggregate = ['unique_counts', 'ambiguous_counts', 'em_counts',
                          'em_cpm', 'unique_cpm', 'ambiguous_cpm']

    aggregated_layers = {}

    # Aggregate X matrix
    X_sparse = sparse.csr_matrix(adata_tx.X) if not sparse.issparse(adata_tx.X) else adata_tx.X
    X_gene = (X_sparse @ mapping_matrix).toarray()

    # Aggregate all layers
    for layer_name in layers_to_aggregate:
        if layer_name in adata_tx.layers:
            layer_sparse = sparse.csr_matrix(adata_tx.layers[layer_name]) if not sparse.issparse(adata_tx.layers[layer_name]) else adata_tx.layers[layer_name]

            if 'ambiguous' in layer_name and 'cpm' not in layer_name:
                # For ambiguous counts, take mean
                layer_sum = layer_sparse @ mapping_matrix
                transcripts_per_gene = mapping_matrix.sum(axis=0).A1
                transcripts_per_gene[transcripts_per_gene == 0] = 1
                aggregated_layers[layer_name] = (layer_sum / transcripts_per_gene[np.newaxis, :]).toarray()
            else:
                # For others, sum
                aggregated_layers[layer_name] = (layer_sparse @ mapping_matrix).toarray()

    # Create gene-level var DataFrame efficiently
    gene_var = pd.DataFrame(index=unique_genes)
    gene_var['gene_id'] = unique_genes
    gene_var['feature_type'] = 'gene'

    # Aggregate metadata using groupby (more efficient than loops)
    gene_metadata_cols = [
        'transcript_id', 'Synt_id', 'synteny_category', 'syntenic_genes',
        'haplotype', 'CDS_length_category', 'CDS_percent_difference',
        'CDS_haplotype_with_longest_annotation', 'functional_annotation'
    ]

    valid_tx_mask = gene_ids.notna()
    tx_var_valid = adata_tx.var[valid_tx_mask].copy()

    for col in gene_metadata_cols:
        if col in adata_tx.var.columns:
            gene_metadata = tx_var_valid.groupby('gene_id')[col].first()
            gene_var[col] = gene_metadata.reindex(unique_genes)

    # Create gene-level AnnData
    adata_gene = ad.AnnData(
        X=X_gene,
        obs=adata_tx.obs.copy(),
        var=gene_var
    )

    # Add aggregated layers
    for layer_name, layer_data in aggregated_layers.items():
        adata_gene.layers[layer_name] = layer_data

    # Add summary statistics
    n_transcripts_per_gene = gene_ids.value_counts()
    adata_gene.var['n_transcripts'] = n_transcripts_per_gene.reindex(
        adata_gene.var_names, fill_value=0
    )

    print(f"Created gene-level AnnData: {adata_gene.n_obs} × {adata_gene.n_vars}")
    print(f"Average transcripts per gene: {adata_gene.var['n_transcripts'].mean():.2f}")

    return adata_gene
