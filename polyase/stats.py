def test_allelic_ratios_within_conditions(adata, layer="unique_counts", test_condition="control", inplace=True):
    """
    Test if alleles of a gene have unequal expression and store results in AnnData object.

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    test_condition : str, optional
        Variable column name containing condition for testing within (default: "control")
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)

    Returns
    --------
    AnnData or None
        If inplace=False, returns modified copy of AnnData; otherwise returns None
        Results are stored in:
        - adata.uns['allelic_ratio_test']: Complete test results as DataFrame
        - adata.var['allelic_ratio_pval']: P-values for each allele
        - adata.var['allelic_ratio_FDR']: FDR-corrected p-values for each allele
    pd.DataFrame
        Results of statistical tests for each syntelog
    """
    import pandas as pd
    import numpy as np
    import re
    from statsmodels.stats.multitest import multipletests
    from isotools._transcriptome_stats import betabinom_lr_test
    from anndata import AnnData

    # Validate inputs
    if not isinstance(adata, AnnData):
        raise ValueError("Input adata must be an AnnData object")

    # Check if layer exists
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    # Work on a copy if not inplace
    if not inplace:
        adata = adata.copy()

    # Get counts and metadata
    counts = adata.layers[layer].copy()  # Create a copy to avoid modifying original

    # Ensure allelic ratio layer exists
    if "allelic_ratio_unique_counts" not in adata.layers:
        raise ValueError("Layer 'allelic_ratio_unique_counts' not found in AnnData object")
    allelic_ratio_counts = adata.layers["allelic_ratio_unique_counts"].copy()

    # Get CPM data if available
    cpm_layer_name = layer.replace('counts', 'cpm')  # e.g., unique_counts -> unique_cpm
    cpm_counts = None
    if cpm_layer_name in adata.layers:
        cpm_counts = adata.layers[cpm_layer_name].copy()
        print(f"Using CPM data from layer: {cpm_layer_name}")
    else:
        print(f"CPM layer '{cpm_layer_name}' not found, CPM values will not be included")

    # Check for syntelog IDs
    if "Synt_id" not in adata.var:
        raise ValueError("'Synt_id' not found in adata.var")
    synt_ids = adata.var["Synt_id"]

    # Check for transcript IDs
    if not adata.var_names.any():
        raise ValueError("'transcript_id' not found in adata.var_names")
    gene_ids = adata.var_names
    transcript_ids = adata.var['transcript_id']

    # Check conditions
    if test_condition not in adata.obs['condition'].unique() and test_condition != "all":
        raise ValueError(f"Condition '{test_condition}' not found in adata.obs['condition']")

    unique_synt_ids = np.unique(synt_ids)

    # Prepare results dataframe
    results = []

    # Create empty arrays for storing p-values in adata.var
    pvals = np.full(adata.n_vars, np.nan)
    fdr_pvals = np.full(adata.n_vars, np.nan)
    ratio_diff = np.full(adata.n_vars, np.nan)

    # Create empty arrays for mean ratios per condition
    mean_ratio_cond1 = np.full(adata.n_vars, np.nan)
    mean_ratio_cond2 = np.full(adata.n_vars, np.nan)

    # Track progress
    total_syntelogs = len(unique_synt_ids)
    processed = 0

    # Process each syntelog
    for synt_id in unique_synt_ids:
        processed += 1
        if processed % 100 == 0:
            print(f"Processing syntelog {processed}/{total_syntelogs}")

        # Find alleles (observations) belonging to this syntelog
        allele_indices = np.where(synt_ids == synt_id)[0]

        # Skip if fewer than 2 alleles found (need at least 2 for ratio testing)
        if len(allele_indices) < 2:
            continue

        for allele_idx, allele_pos in enumerate(allele_indices):
            allele_counts = []
            condition_total = []
            allelic_ratios = {}

            # Get samples for this condition
            if test_condition == "all":
                condition_indices = np.arange(counts.shape[0])
            else:
                # Get samples for this condition
                condition_indices = np.where(adata.obs['condition'] == test_condition)[0]

            # Extract counts for these alleles and samples
            condition_counts = counts[np.ix_(condition_indices, allele_indices)]

            # Sum across samples to get total counts per allele for this condition
            total_counts = np.sum(condition_counts, axis=1)

            # Get allelic ratios for this condition
            condition_ratios = allelic_ratio_counts[np.ix_(condition_indices, allele_indices)]

            # Get CPM values for this condition if available
            condition_cpm = None
            if cpm_counts is not None:
                condition_cpm = cpm_counts[np.ix_(condition_indices, allele_indices)]

            # Append arrays for total counts
            condition_total.append(total_counts)

            # Append array for this specific allele's counts
            allele_counts.append(condition_counts[:,allele_idx])

            # Store ratios for this test condition
            allelic_ratios = condition_ratios[:,allele_idx]

            # generate balanced allele counts based on condition total counts
            # balanced counts need to be integers for the test
            balanced_counts = [np.round(x * 1/len(allele_indices)) for x in condition_total]
            allele_counts.append(balanced_counts[0])
            # add the total counts again for the balanced counts
            condition_total.append(total_counts)

            # Run the beta-binomial likelihood ratio test
            try:
                test_result = betabinom_lr_test(allele_counts, condition_total)
                p_value, ratio_stats = test_result[0], test_result[1]
                # Calculate absolute difference in mean ratios between conditions
                ratio_difference = abs(ratio_stats[0] - ratio_stats[2])
            except Exception as e:
                print(f"Error testing syntelog {synt_id}, allele {allele_idx}: {str(e)}")
                continue

            # Get gene ID and parse allele info
            gene_id = gene_ids[allele_pos]
            # Get transcript ID
            transcript_id = transcript_ids[allele_pos]

            haplotype = adata.var['haplotype'].iloc[allele_indices[allele_idx]]
            # Extract allele number from haplotype

            try:
                allele_match = re.search(r'hap(\d+)', haplotype)  # Capture the number
                if allele_match:
                    allele_num = allele_match.group(1)  # Get the captured number directly
                else:
                    allele_num = f"{allele_idx+1}"  # Fallback if regex fails
                    print(f"No match found, using fallback: {allele_num}")
            except Exception as e:
                print(f"Error: {e}")
                allele_num = f"{allele_idx+1}"  # Fallback if any error occurs

            # Store p-value in the arrays we created
            pvals[allele_pos] = p_value
            ratio_diff[allele_pos] = ratio_difference
            mean_ratio_cond1[allele_pos] = ratio_stats[0]
            mean_ratio_cond2[allele_pos] = ratio_stats[2]

            # Prepare result dictionary
            result_dict = {
                'Synt_id': synt_id,
                'allele': allele_num,
                'transcript_id': transcript_id,
                'p_value': p_value,
                'ratio_difference': ratio_difference,
                'n_alleles': len(allele_indices),
                f'ratios_{test_condition}_mean': ratio_stats[0],
                f'ratios_rep_{test_condition}': allelic_ratios
            }

            # Add CPM values if available
            if condition_cpm is not None:
                # Calculate mean CPM for this allele in this condition
                allele_cpm_values = condition_cpm[:, allele_idx]
                mean_cpm = np.mean(allele_cpm_values)

                result_dict[f'cpm_{test_condition}_mean'] = mean_cpm
                result_dict[f'cpm_rep_{test_condition}'] = allele_cpm_values

            results.append(result_dict)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Multiple testing correction if we have results
    if len(results_df) > 0:
        # PROBLEM: p_value is nan sometimes, replace with 1 for now
        results_df['p_value'] = results_df['p_value'].fillna(1)
        results_df['FDR'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df = results_df.sort_values('p_value')

        # Map FDR values back to the individual alleles
        # Group by transcript_id and take the first FDR value (they should be the same for all replicates)
        fdr_map = results_df.groupby('transcript_id')['FDR'].first().to_dict()

        # Update the FDR array
        for i, transcript_id in enumerate(transcript_ids):
            if transcript_id in fdr_map:
                fdr_pvals[i] = fdr_map[transcript_id]

    # Store results in the AnnData object
    adata.uns['allelic_ratio_test'] = results_df
    adata.var['allelic_ratio_pval'] = pvals
    adata.var['allelic_ratio_FDR'] = fdr_pvals
    adata.var['allelic_ratio_difference'] = ratio_diff
    adata.var[f'allelic_ratio_mean_{test_condition}'] = mean_ratio_cond1

    # Group by Synt_id and take minimum FDR value and max ratio difference
    if len(results_df) > 0:
        grouped_results = results_df.groupby('Synt_id').agg({
            'FDR': 'min',
            'ratio_difference': 'max'
        })
        # Print summary
        significant_results = grouped_results[(grouped_results['FDR'] < 0.005) & (grouped_results['ratio_difference'] > 0.1)]
        print(f"Found {len(significant_results)} from {len(grouped_results)} syntelogs with at least one significantly different allele (FDR < 0.005 and ratio difference > 0.1)")

    # Return AnnData object if not inplace
    if not inplace:
        return adata
    else:
        return results_df


def test_allelic_ratios_between_conditions(adata, layer="unique_counts", group_key="condition", inplace=True):
    """
    Test if allelic ratios change between conditions and store results in AnnData object.

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    group_key : str, optional
        Variable column name containing condition information (default: "condition")
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)

    Returns
    --------
    AnnData or None
        If inplace=False, returns modified copy of AnnData; otherwise returns None
        Results are stored in:
        - adata.uns['allelic_ratio_test']: Complete test results as DataFrame
        - adata.var['allelic_ratio_pval']: P-values for each allele
        - adata.var['allelic_ratio_FDR']: FDR-corrected p-values for each allele
    pd.DataFrame
        Results of statistical tests for each syntelog
    """
    import pandas as pd
    import numpy as np
    import re
    from statsmodels.stats.multitest import multipletests
    from isotools._transcriptome_stats import betabinom_lr_test
    from anndata import AnnData

    # Validate inputs
    if not isinstance(adata, AnnData):
        raise ValueError("Input adata must be an AnnData object")

    # Check if layer exists
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    # Check if group_key exists in obs
    if group_key not in adata.obs:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")

    # Work on a copy if not inplace
    if not inplace:
        adata = adata.copy()

    # Get counts and metadata
    counts = adata.layers[layer].copy()  # Create a copy to avoid modifying original

    # Ensure allelic ratio layer exists
    if "allelic_ratio_unique_counts" not in adata.layers:
        raise ValueError("Layer 'allelic_ratio_unique_counts' not found in AnnData object")
    allelic_ratio_counts = adata.layers["allelic_ratio_unique_counts"].copy()

    # Get CPM data if available
    cpm_layer_name = layer.replace('counts', 'cpm')  # e.g., unique_counts -> unique_cpm
    cpm_counts = None
    if cpm_layer_name in adata.layers:
        cpm_counts = adata.layers[cpm_layer_name].copy()
        print(f"Using CPM data from layer: {cpm_layer_name}")
    else:
        print(f"CPM layer '{cpm_layer_name}' not found, CPM values will not be included")

    # Check for syntelog IDs
    if "Synt_id" not in adata.var:
        raise ValueError("'Synt_id' not found in adata.var")
    synt_ids = adata.var["Synt_id"]

    # Check for transcript IDs
    if not adata.var_names.any():
        raise ValueError("'transcript_id' not found in adata.var_names")
    gene_ids = adata.var_names
    transcript_ids = adata.var['transcript_id']

    # Check conditions
    if group_key not in adata.obs:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")
    conditions = adata.obs[group_key].values

    # Get unique conditions and syntelog IDs
    unique_conditions = np.unique(conditions)
    if len(unique_conditions) != 2:
        raise ValueError(f"Need exactly 2 conditions, found {len(unique_conditions)}: {unique_conditions}")

    unique_synt_ids = np.unique(synt_ids)

    # Prepare results dataframe
    results = []

    # Create empty arrays for storing p-values in adata.var
    pvals = np.full(adata.n_vars, np.nan)
    fdr_pvals = np.full(adata.n_vars, np.nan)
    ratio_diff = np.full(adata.n_vars, np.nan)

    # Create empty arrays for mean ratios per condition
    mean_ratio_cond1 = np.full(adata.n_vars, np.nan)
    mean_ratio_cond2 = np.full(adata.n_vars, np.nan)

    # Track progress
    total_syntelogs = len(unique_synt_ids)
    processed = 0

    # Process each syntelog
    for synt_id in unique_synt_ids:
        processed += 1
        if processed % 100 == 0:
            print(f"Processing syntelog {processed}/{total_syntelogs}")

        # Find alleles (observations) belonging to this syntelog
        allele_indices = np.where(synt_ids == synt_id)[0]

        # Skip if fewer than 2 alleles found (need at least 2 for ratio testing)
        if len(allele_indices) < 2:
            continue

        for allele_idx, allele_pos in enumerate(allele_indices):
            allele_counts = []
            condition_total = []
            allelic_ratios = {}
            cpm_values = {}

            for condition_idx, condition in enumerate(unique_conditions):
                # Get samples for this condition
                condition_indices = np.where(conditions == condition)[0]

                # Extract counts for these alleles and samples
                condition_counts = counts[np.ix_(condition_indices, allele_indices)]

                # Sum across samples to get total counts per allele for this condition
                total_counts = np.sum(condition_counts, axis=1)

                # Get allelic ratios for this condition
                condition_ratios = allelic_ratio_counts[np.ix_(condition_indices, allele_indices)]

                # Get CPM values for this condition if available
                if cpm_counts is not None:
                    condition_cpm = cpm_counts[np.ix_(condition_indices, allele_indices)]
                    cpm_values[condition] = condition_cpm[:, allele_idx]

                # Append arrays for total counts
                condition_total.append(total_counts)

                # Append array for this specific allele's counts
                allele_counts.append(condition_counts[:,allele_idx])

                # Store ratios for this condition
                allelic_ratios[condition] = condition_ratios[:,allele_idx]

            # Run the beta-binomial likelihood ratio test
            try:
                test_result = betabinom_lr_test(allele_counts, condition_total)
                p_value, ratio_stats = test_result[0], test_result[1]

                # Calculate absolute difference in mean ratios between conditions
                ratio_difference = abs(ratio_stats[0] - ratio_stats[2])
            except Exception as e:
                print(f"Error testing syntelog {synt_id}, allele {allele_idx}: {str(e)}")
                continue

            # Get transcript ID and parse allele info
            gene_id = gene_ids[allele_pos]
            transcript_id = transcript_ids[allele_pos]

            haplotype = adata.var['haplotype'].iloc[allele_indices[allele_idx]]
            # Extract allele number from haplotype
            try:
                allele_match = re.search(r'hap(\d+)', haplotype)  # Capture the number
                if allele_match:
                    allele_num = allele_match.group(1)  # Get the captured number directly
                else:
                    allele_num = f"{allele_idx+1}"  # Fallback if regex fails
                    print(f"No match found, using fallback: {allele_num}")
            except Exception as e:
                print(f"Error: {e}")
                allele_num = f"{allele_idx+1}"  # Fallback if any error occurs

            # Store p-value in the arrays we created
            pvals[allele_pos] = p_value
            ratio_diff[allele_pos] = ratio_difference
            mean_ratio_cond1[allele_pos] = ratio_stats[0]
            mean_ratio_cond2[allele_pos] = ratio_stats[2]

            # Prepare result dictionary
            result_dict = {
                'Synt_id': synt_id,
                'gene_id': gene_id,
                'transcript_id': transcript_id,
                'allele': allele_num,
                'p_value': p_value,
                'ratio_difference': ratio_difference,
                'n_alleles': len(allele_indices),
                f'ratios_{unique_conditions[0]}_mean': ratio_stats[0],
                f'ratios_rep_{unique_conditions[0]}': allelic_ratios[unique_conditions[0]],
                f'ratios_{unique_conditions[1]}_mean': ratio_stats[2],
                f'ratios_rep_{unique_conditions[1]}': allelic_ratios[unique_conditions[1]]
            }

            # Add CPM values if available
            if cpm_counts is not None:
                for condition in unique_conditions:
                    if condition in cpm_values:
                        mean_cpm = np.mean(cpm_values[condition])
                        result_dict[f'cpm_{condition}_mean'] = mean_cpm
                        result_dict[f'cpm_rep_{condition}'] = cpm_values[condition]

            results.append(result_dict)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Multiple testing correction if we have results
    if len(results_df) > 0:
        # PROBLEM: p_value is nan sometimes, replace with 1 for now
        results_df['p_value'] = results_df['p_value'].fillna(1)
        results_df['FDR'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df = results_df.sort_values('p_value')

        # Map FDR values back to the individual alleles
        # Group by transcript_id and take the first FDR value (they should be the same for all replicates)
        fdr_map = results_df.groupby('transcript_id')['FDR'].first().to_dict()

        # Update the FDR array
        for i, transcript_id in enumerate(transcript_ids):
            if transcript_id in fdr_map:
                fdr_pvals[i] = fdr_map[transcript_id]

    # Store results in the AnnData object
    adata.uns['allelic_ratio_test'] = results_df
    adata.var['allelic_ratio_pval'] = pvals
    adata.var['allelic_ratio_FDR'] = fdr_pvals
    adata.var['allelic_ratio_difference'] = ratio_diff
    adata.var[f'allelic_ratio_mean_{unique_conditions[0]}'] = mean_ratio_cond1
    adata.var[f'allelic_ratio_mean_{unique_conditions[1]}'] = mean_ratio_cond2

    # Group by Synt_id and take minimum FDR value
    if len(results_df) > 0:
        grouped_results = results_df.groupby('Synt_id').min("FDR")
        # Print summary
        significant_results = grouped_results[(grouped_results['FDR'] < 0.05)]
        print(f"Found {len(significant_results)} from {len(grouped_results)} syntelogs with at least one significantly different allelic ratio (FDR < 0.05)")

    # Return AnnData object if not inplace
    if not inplace:
        return adata
    else:
        return results_df

def get_top_differential_syntelogs(results_df, n=5, sort_by='p_value', fdr_threshold=0.05, ratio_threshold=0.1):
    """
    Get the top n syntelogs with differential allelic ratios.

    Parameters
    -----------
    results_df : pd.DataFrame
        Results dataframe from test_allelic_ratios function
    n : int, optional
        Number of top syntelogs to return (default: 5)
    sort_by : str, optional
        Column to sort results by ('p_value', 'FDR', or 'ratio_difference') (default: 'p_value')
    fdr_threshold : float, optional
        Maximum FDR to consider a result significant (default: 0.05)

    Returns
    --------
    pd.DataFrame
        Filtered dataframe containing only the top n syntelogs
    """
    if len(results_df) == 0:
        print("No results to filter")
        return results_df

    # Validate sort_by parameter
    if sort_by not in ['p_value', 'FDR', 'ratio_difference']:
        print(f"Invalid sort_by parameter '{sort_by}'. Using 'p_value' instead.")
        sort_by = 'p_value'

    if sort_by == 'ratio_difference':
        sort_bool = False
    else:
        sort_bool = True

    # Apply FDR filter if column exists
    if 'FDR' in results_df.columns:
        sig_results = results_df[(results_df['FDR'] <= fdr_threshold) & (results_df['ratio_difference'] >= ratio_threshold)]

        if len(sig_results) == 0:
            print(f"No results with FDR <= {fdr_threshold} and ratio_difference >= {ratio_threshold}. Using all results.")
            sig_results = results_df
    else:
        sig_results = results_df

    # Get top n syntelogs
    top_syntelogs = sig_results.sort_values(sort_by, ascending=sort_bool).drop_duplicates('Synt_id').head(n)['Synt_id'].unique()

    # Return filtered dataframe
    return results_df[results_df['Synt_id'].isin(top_syntelogs)]


def test_isoform_DIU_between_conditions(adata, layer="unique_counts", group_key="condition", gene_id_key="gene_id", inplace=True):
    """
    Test if isoform usage ratios change between conditions and store results in AnnData object.

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    group_key : str, optional
        Variable column name containing condition information (default: "condition")
    gene_id_key : str, optional
        Variable column name containing gene ID information (default: "gene_id")
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)

    Returns
    --------
    AnnData or None
        If inplace=False, returns modified copy of AnnData; otherwise returns None
        Results are stored in:
        - adata.uns['isoform_usage_test']: Complete test results as DataFrame
        - adata.var['isoform_usage_pval']: P-values for each isoform
        - adata.var['isoform_usage_FDR']: FDR-corrected p-values for each isoform
    pd.DataFrame
        Results of statistical tests for each gene
    """
    import pandas as pd
    import numpy as np
    import re
    from statsmodels.stats.multitest import multipletests
    from isotools._transcriptome_stats import betabinom_lr_test
    from anndata import AnnData

    # Validate inputs
    if not isinstance(adata, AnnData):
        raise ValueError("Input adata must be an AnnData object")

    # Check if layer exists
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    # Check if group_key exists in obs
    if group_key not in adata.obs:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")

    # Check if gene_id_key exists in var
    if gene_id_key not in adata.var:
        raise ValueError(f"Gene ID key '{gene_id_key}' not found in adata.var")

    # Work on a copy if not inplace
    if not inplace:
        adata = adata.copy()

    no_counts_isoform = 0


    # Get counts and metadata
    counts = adata.layers[layer].copy()  # Create a copy to avoid modifying original
    gene_ids = adata.var[gene_id_key]
    transcript_ids = adata.var_names
    conditions = adata.obs[group_key].values

    # Get unique conditions and gene IDs
    unique_conditions = np.unique(conditions)
    if len(unique_conditions) != 2:
        raise ValueError(f"Need exactly 2 conditions, found {len(unique_conditions)}: {unique_conditions}")

    unique_gene_ids = np.unique(gene_ids)

    # Calculate isoform ratios for each gene
    print("Calculating isoform ratios...")
    isoform_ratios = np.zeros_like(counts, dtype=float)


    for gene_id in unique_gene_ids:
        # Find isoforms (variables) belonging to this gene
        isoform_indices = np.where(gene_ids == gene_id)[0]

        if len(isoform_indices) < 2:
            continue  # Skip genes with only one isoform

        # Calculate total gene expression for each sample
        gene_totals = np.sum(counts[:, isoform_indices], axis=1, keepdims=True)

        # Avoid division by zero
        gene_totals[gene_totals == 0] = 1

        # Calculate isoform ratios
        isoform_ratios[:, isoform_indices] = counts[:, isoform_indices] / gene_totals

    # Store isoform ratios in a new layer
    adata.layers['isoform_ratios'] = isoform_ratios

    # Prepare results dataframe
    results = []

    # Create empty arrays for storing p-values in adata.var
    pvals = np.full(adata.n_vars, np.nan)
    fdr_pvals = np.full(adata.n_vars, np.nan)
    ratio_diff = np.full(adata.n_vars, np.nan)

    # Create empty arrays for mean ratios per condition
    mean_ratio_cond1 = np.full(adata.n_vars, np.nan)
    mean_ratio_cond2 = np.full(adata.n_vars, np.nan)

    # Track progress
    total_genes = len(unique_gene_ids)
    processed = 0

    # Process each gene
    for gene_id in unique_gene_ids:
        processed += 1
        if processed % 100 == 0:
            print(f"Processing gene {processed}/{total_genes}")

        # Find isoforms (variables) belonging to this gene
        isoform_indices = np.where(gene_ids == gene_id)[0]

        # Skip if fewer than 2 isoforms found (need at least 2 for ratio testing)
        if len(isoform_indices) < 2:
            continue

        # Test each isoform within this gene
        for isoform_idx, isoform_pos in enumerate(isoform_indices):
            isoform_counts = []
            gene_total_counts = []
            isoform_ratios_per_condition = {}

            for condition_idx, condition in enumerate(unique_conditions):
                # Get samples for this condition
                condition_indices = np.where(conditions == condition)[0]

                # Extract counts for all isoforms of this gene in this condition
                condition_gene_counts = counts[np.ix_(condition_indices, isoform_indices)]

                # Get total gene counts per sample (sum across all isoforms)
                condition_gene_totals = np.sum(condition_gene_counts, axis=1)

                # Get this specific isoform's counts
                condition_isoform_counts = counts[np.ix_(condition_indices, [isoform_pos])].flatten()

                # Store data for beta-binomial test
                isoform_counts.append(condition_isoform_counts)
                gene_total_counts.append(condition_gene_totals)

                # Calculate isoform ratios for this condition
                condition_ratios = np.divide(condition_isoform_counts, condition_gene_totals,
                                           out=np.zeros_like(condition_isoform_counts, dtype=float),
                                           where=condition_gene_totals!=0)
                isoform_ratios_per_condition[condition] = condition_ratios

            # Run the beta-binomial likelihood ratio test
            # if isoform counts or gene total counts are 0, skip this isoform
            isoform_counts = np.array(isoform_counts)
            gene_total_counts = np.array(gene_total_counts)
            if np.any(isoform_counts == 0) or np.any(gene_total_counts == 0):
                no_counts_isoform = no_counts_isoform + 1
                # print(f"Skipping gene {gene_id}, isoform {isoform_idx + 1} due to zero counts")
                continue
            try:
                test_result = betabinom_lr_test(isoform_counts, gene_total_counts)
                p_value, ratio_stats = test_result[0], test_result[1]

                # Calculate absolute difference in mean ratios between conditions
                ratio_difference = abs(ratio_stats[0] - ratio_stats[2])
            except Exception as e:
                print(f"Error testing gene {gene_id}, isoform {isoform_idx}: {str(e)}")
                continue

            # Get transcript ID
            transcript_id = transcript_ids[isoform_pos]

            # Store p-value in the arrays we created
            pvals[isoform_pos] = p_value
            ratio_diff[isoform_pos] = ratio_difference
            mean_ratio_cond1[isoform_pos] = ratio_stats[0]
            mean_ratio_cond2[isoform_pos] = ratio_stats[2]

            # Store results
            results.append({
                'gene_id': gene_id,
                'isoform_number': isoform_idx + 1,
                'transcript_id': transcript_id,
                'p_value': p_value,
                'ratio_difference': ratio_difference,
                'n_isoforms': len(isoform_indices),
                f'ratios_{unique_conditions[0]}_mean': ratio_stats[0],
                f'ratios_rep_{unique_conditions[0]}': isoform_ratios_per_condition[unique_conditions[0]],
                f'ratios_{unique_conditions[1]}_mean': ratio_stats[2],
                f'ratios_rep_{unique_conditions[1]}': isoform_ratios_per_condition[unique_conditions[1]]
            })

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Multiple testing correction if we have results
    if len(results_df) > 0:
        # Handle NaN p-values
        results_df['p_value'] = results_df['p_value'].fillna(1)
        results_df['FDR'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df = results_df.sort_values('p_value')

        # Map FDR values back to the individual isoforms
        fdr_map = results_df.groupby('transcript_id')['FDR'].first().to_dict()

        # Update the FDR array
        for i, transcript_id in enumerate(transcript_ids):
            if transcript_id in fdr_map:
                fdr_pvals[i] = fdr_map[transcript_id]

    # Store results in the AnnData object
    adata.uns['isoform_usage_test'] = results_df
    adata.var['isoform_usage_pval'] = pvals
    adata.var['isoform_usage_FDR'] = fdr_pvals
    adata.var['isoform_usage_difference'] = ratio_diff
    adata.var[f'isoform_usage_mean_{unique_conditions[0]}'] = mean_ratio_cond1
    adata.var[f'isoform_usage_mean_{unique_conditions[1]}'] = mean_ratio_cond2

    # Group by gene_id and take minimum FDR value
    grouped_results = results_df.groupby('gene_id').agg({
        'FDR': 'min',
        'p_value': 'min',
        'n_isoforms': 'first'
    }).reset_index()

    # Print summary
    significant_results = grouped_results[grouped_results['FDR'] < 0.05]
    print(f"Found {len(significant_results)} from {len(grouped_results)} genes with at least one significantly different isoform usage (FDR < 0.05)")
    print(f"Skipped {no_counts_isoform} isoforms due to zero counts")

    # Return AnnData object if not inplace
    if not inplace:
        return adata
    else:
        return results_df
def test_isoform1_DIU_between_alleles(adata, layer="unique_counts", test_condition="control", inplace=True):
    """
    Test if alleles have different isoform usage and store results in AnnData object.

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    test_condition : str, optional
        Variable column name containing condition for testing within (default: "control")
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)

    Returns
    --------
    pd.DataFrame
        Results of statistical tests for each syntelog
    """
    import pandas as pd
    import numpy as np
    import re
    from statsmodels.stats.multitest import multipletests
    from isotools._transcriptome_stats import betabinom_lr_test
    from anndata import AnnData

    # Validate inputs
    if not isinstance(adata, AnnData):
        raise ValueError("Input adata must be an AnnData object")

    # Check if layer exists
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    # Work on a copy if not inplace
    if not inplace:
        adata = adata.copy()

    # Get counts and metadata
    counts = adata.layers[layer].copy()

    # Check for required columns
    if "Synt_id" not in adata.var:
        raise ValueError("'Synt_id' not found in adata.var")
    synt_ids = adata.var["Synt_id"]

    if "haplotype" not in adata.var:
        raise ValueError("'haplotype' not found in adata.var")
    haplotypes = adata.var["haplotype"]

    # Check for transcript IDs
    if not adata.var_names.any():
        raise ValueError("'transcript_id' not found in adata.var_names")
    transcript_ids = adata.var_names

    # Check conditions
    if test_condition not in adata.obs['condition'].unique() and test_condition != "all":
        raise ValueError(f"Condition '{test_condition}' not found in adata.obs['condition']")

    unique_synt_ids = np.unique(synt_ids)

    # Remove NaN and 0 values from unique_synt_ids
    unique_synt_ids = unique_synt_ids[~pd.isna(unique_synt_ids)]
    unique_synt_ids = unique_synt_ids[unique_synt_ids != 0]

    # Prepare results dataframe
    results = []

    # Track progress
    total_syntelogs = len(unique_synt_ids)
    processed = 0

    # Process each syntelog
    for synt_id in unique_synt_ids:
        processed += 1
        if processed % 100 == 0:
            print(f"Processing syntelog {processed}/{total_syntelogs}")

        # Find all transcripts belonging to this syntelog
        synt_mask = synt_ids == synt_id
        synt_indices = np.where(synt_mask)[0]

        # Skip if no transcripts found
        if len(synt_indices) == 0:
            continue

        # Get unique haplotypes for this syntelog
        synt_haplotypes = haplotypes.iloc[synt_indices]
        unique_haplotypes = synt_haplotypes.unique()

        # Skip if fewer than 2 haplotypes and not more than one isoform
        if len(unique_haplotypes) < 2 and len(synt_indices) < (len(unique_haplotypes) *2):
            print(f"Skipping syntelog {synt_id} with less than 2 haplotypes or less than 2 isoforms")
            continue


        # Get sample indices for the test condition
        if test_condition == "all":
            condition_indices = np.arange(counts.shape[0])
        else:
            condition_indices = np.where(adata.obs['condition'] == test_condition)[0]

        # Find the most expressed isoform across all haplotypes for this syntelog
        synt_counts = counts[np.ix_(condition_indices, synt_indices)]
        total_counts_per_transcript = np.sum(synt_counts, axis=0)

        # Skip if all counts are zero
        if np.sum(total_counts_per_transcript) == 0:
            continue

        max_count_local_idx = np.argmax(total_counts_per_transcript)
        max_count_global_idx = synt_indices[max_count_local_idx]
        max_count_transcript_id = transcript_ids[max_count_global_idx]

        # Extract the isoform number from the most expressed transcript
        try:
            isoform_match = re.search(r'\.(\d+)\.', max_count_transcript_id)
            if isoform_match:
                target_isoform = isoform_match.group(1)
            elif isoform_match is None:
                # Try to extract the transcript name without isoform number
                target_isoform = max_count_transcript_id.split('.')[0]

                if target_isoform is None or target_isoform == "":
                    # If we couldn't extract a valid isoform number, skip this syntelog
                    print(f"Could not extract isoform number from {max_count_transcript_id}, skipping syntelog {synt_id}")
                    continue
        except Exception as e:
            print(f"Error extracting isoform from {max_count_transcript_id}: {e}")
            continue

        # Find the same isoform in each haplotype
        haplotype_isoform_data = {}

        for hap in unique_haplotypes:
            # Get indices for this haplotype within the syntelog
            hap_mask = synt_haplotypes == hap
            hap_indices_local = np.where(hap_mask)[0]
            hap_indices_global = synt_indices[hap_indices_local]

            # Find the target isoform in this haplotype
            target_isoform_idx = None
            for idx in hap_indices_global:
                transcript_id = transcript_ids[idx]

                try:
                    # First try to match the isoform number pattern
                    isoform_match = re.search(r'\.(\d+)\.', transcript_id)
                    if isoform_match and isoform_match.group(1) == target_isoform:
                        target_isoform_idx = idx
                        break
                    else:
                        # Fallback: try different patterns or matching strategies
                        # This depends on your specific transcript ID format
                        transcript_parts = transcript_id.split('.')
                        if len(transcript_parts) >= 2 and transcript_parts[1] == target_isoform:
                            target_isoform_idx = idx
                            break
                except:
                    continue

            if target_isoform_idx is not None:
                # Get counts for this specific isoform
                isoform_counts = counts[np.ix_(condition_indices, [target_isoform_idx])][:, 0]

                # Get total counts for all isoforms of this haplotype in this syntelog
                hap_total_counts = np.sum(counts[np.ix_(condition_indices, hap_indices_global)], axis=1)

                haplotype_isoform_data[hap] = {
                    'isoform_counts': isoform_counts,
                    'total_counts': hap_total_counts,
                    'transcript_id': transcript_ids[target_isoform_idx]
                }

        # Skip if we don't have the target isoform in at least 2 haplotypes
        if len(haplotype_isoform_data) < 2:
            print(f"Target isoform {target_isoform} not found in enough haplotypes for syntelog {synt_id}")
            continue

        # Calculate average ratios for each haplotype
        haplotype_ratios = {}
        for hap, data in haplotype_isoform_data.items():
            # Avoid division by zero
            valid_samples = data['total_counts'] > 0
            if np.sum(valid_samples) > 0:
                ratios = np.zeros_like(data['total_counts'], dtype=float)
                ratios[valid_samples] = data['isoform_counts'][valid_samples] / data['total_counts'][valid_samples]
                haplotype_ratios[hap] = np.mean(ratios)
            else:
                haplotype_ratios[hap] = 0.0

        # Find haplotypes with max and min ratios
        if len(haplotype_ratios) < 2:
            continue

        sorted_haps = sorted(haplotype_ratios.keys(), key=lambda x: haplotype_ratios[x])
        min_hap = sorted_haps[0]
        max_hap = sorted_haps[-1]

        # Skip if ratios are the same (no difference to test)
        if haplotype_ratios[min_hap] == haplotype_ratios[max_hap]:
            continue

        # Prepare data for statistical test
        min_hap_data = haplotype_isoform_data[min_hap]
        max_hap_data = haplotype_isoform_data[max_hap]

        allele_counts = [min_hap_data['isoform_counts'], max_hap_data['isoform_counts']]
        condition_total = [min_hap_data['total_counts'], max_hap_data['total_counts']]

        # Skip if any total counts are zero
        if np.any([np.sum(ct) == 0 for ct in condition_total]):
            print(f"Skipping syntelog {synt_id} with zero total counts in some haplotypes.")
            continue

        # Run the beta-binomial likelihood ratio test
        try:
            test_result = betabinom_lr_test(allele_counts, condition_total)
            p_value, ratio_stats = test_result[0], test_result[1]

            # Calculate absolute difference in mean ratios between haplotypes
            ratio_difference = abs(ratio_stats[0] - ratio_stats[2]) if len(ratio_stats) >= 3 else abs(haplotype_ratios[min_hap] - haplotype_ratios[max_hap])

        except Exception as e:
            print(f"Error testing syntelog {synt_id}: {str(e)}")
            continue

        # Prepare result dictionary
        result_dict = {
            'Synt_id': synt_id,
            'target_isoform': target_isoform,
            'min_ratio_haplotype': min_hap,
            'max_ratio_haplotype': max_hap,
            'min_ratio_transcript_id': min_hap_data['transcript_id'],
            'max_ratio_transcript_id': max_hap_data['transcript_id'],
            'p_value': p_value,
            'ratio_difference': ratio_difference,
            'n_haplotypes': len(haplotype_isoform_data),
            f'ratio_{min_hap}_mean': haplotype_ratios[min_hap],
            f'ratio_{max_hap}_mean': haplotype_ratios[max_hap]
        }

        # Add mean ratios for ALL haplotypes (including those not tested)
        for hap, ratio in haplotype_ratios.items():
            result_dict[f'ratio_{hap}_mean'] = ratio

        # Add transcript IDs for all haplotypes that have the target isoform
        for hap, data in haplotype_isoform_data.items():
            result_dict[f'transcript_id_{hap}'] = data['transcript_id']

        # Add list of all haplotypes for this syntelog
        result_dict['all_haplotypes'] = list(haplotype_ratios.keys())

        # Store results
        results.append(result_dict)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Multiple testing correction if we have results
    if len(results_df) > 0:
        # Handle NaN p-values
        results_df['p_value'] = results_df['p_value'].fillna(1)
        results_df['FDR'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df = results_df.sort_values('p_value')

        # Print summary
        significant_results = results_df[(results_df['FDR'] < 0.05) & (results_df['ratio_difference'] > 0.2)]
        print(f"Found {len(significant_results)} from {len(results_df)} syntelogs with significantly different isoform usage between alleles (FDR < 0.05 and ratio difference > 0.2)")

        # Store results in AnnData object if inplace
        if inplace:
            adata.uns['isoform_diu_test'] = results_df
    else:
        print("No results found")

    return results_df
