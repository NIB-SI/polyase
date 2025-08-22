"""
stats.py
====================================
The stats module of polyase project
"""
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

    # Select the right layer for allelic ratios for plotting
    if layer == "unique_counts":
        allelic_ratio_counts = adata.layers["allelic_ratio_unique_counts"].copy()
    elif layer == "em_counts":
        allelic_ratio_counts = adata.layers["allelic_ratio_em_counts"].copy()
    else:
        raise ValueError("Layer must be either 'allelic_ratio_unique_counts' or 'allelic_ratio_em_counts'")

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

    # Get allelic ratio counts
    if layer == "unique_counts":
        allelic_ratio_counts = adata.layers["allelic_ratio_unique_counts"].copy()
    elif layer == "em_counts":
        allelic_ratio_counts = adata.layers["allelic_ratio_em_counts"].copy()
    else:
        raise ValueError("Layer must be either 'allelic_ratio_unique_counts' or 'allelic_ratio_em_counts'")

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

    if "functional_annotation" in adata.var:
        functional_annotations = adata.var["functional_annotation"]
    else:
        functional_annotations = None
        print("No functional annotations found in adata.var, skipping functional annotation processing.")

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
            transcript_id = transcript_ids.iloc[allele_pos]
            functional_annotation = functional_annotations.iloc[allele_pos]

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
                'functional_annotation': functional_annotation,
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

    if "functional_annotation" in adata.var:
        functional_annotations = adata.var["functional_annotation"]
    else:
        functional_annotations = None
        print("No functional annotations found in adata.var, skipping functional annotation processing.")

    if "gene_id" in adata.var:
        gene_ids = adata.var['gene_id']

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
            functional_annotation = functional_annotations.iloc[allele_pos]

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
                'functional_annotation': functional_annotation,
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
    pd.DataFrame
        Plotting results table with one row per replicate, condition, isoform ratio, and transcript
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

    # Calculate library sizes (total counts per sample) for CPM calculation
    library_sizes = np.sum(counts, axis=1)

    # Calculate CPM (Counts Per Million)
    cpm = np.zeros_like(counts, dtype=float)
    for i, lib_size in enumerate(library_sizes):
        if lib_size > 0:
            cpm[i, :] = (counts[i, :] / lib_size) * 1e6
        else:
            cpm[i, :] = 0

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
    plotting_results = []  # New list for plotting table

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

            # Create plotting results table - one row per replicate, condition, ratio, transcript
            for condition in unique_conditions:
                condition_indices = np.where(conditions == condition)[0]
                sample_names = adata.obs_names[condition_indices]
                ratios = isoform_ratios_per_condition[condition]

                # Get CPM values for this isoform in this condition
                isoform_cpm_values = cpm[condition_indices, isoform_pos]

                # Get raw counts for this isoform in this condition
                isoform_count_values = counts[condition_indices, isoform_pos]

                for rep_idx, (sample_name, ratio_value, cpm_value, count_value) in enumerate(zip(sample_names, ratios, isoform_cpm_values, isoform_count_values)):
                    plotting_results.append({
                        'gene_id': gene_id,
                        'transcript_id': transcript_id,
                        'isoform_number': isoform_idx + 1,
                        'condition': condition,
                        'replicate': rep_idx + 1,
                        'sample_name': sample_name,
                        'isoform_ratio': ratio_value,
                        'unique_counts': count_value,  # Raw counts
                        'unique_counts_cpm': cpm_value,  # CPM values
                        'p_value': p_value,
                        'ratio_difference': ratio_difference,
                        'n_isoforms': len(isoform_indices)
                    })

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    plotting_df = pd.DataFrame(plotting_results)

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

        # Add FDR to plotting dataframe
        plotting_df['FDR'] = plotting_df['transcript_id'].map(fdr_map).fillna(1.0)

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
    print(f"Created plotting table with {len(plotting_df)} rows (one per replicate, condition, isoform ratio, and transcript)")

    # Return AnnData object if not inplace
    if not inplace:
        return adata, results_df, plotting_df
    else:
        return results_df, plotting_df


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


def test_isoform_DIU_between_alleles_with_major_minor_plotting(
    adata, layer="unique_counts", test_condition="control",
    structure_similarity_threshold=0.6, min_similarity_for_matching=0.4,
    inplace=True, verbose=False, return_plotting_data=True
):
    """
    Test if alleles have different isoform usage and return plotting data with both
    major (most expressed) and minor (second most expressed from same haplotype) isoforms.

    This function performs structure-based DIU analysis and returns data in long format
    suitable for plotting, including both the primary isoform and secondary isoform
    from the same haplotype as the reference.

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data with exon structure information
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    test_condition : str, optional
        Variable column name containing condition for testing within (default: "control")
    structure_similarity_threshold : float, optional
        Minimum similarity score (0-1) for considering structures as similar (default: 0.6)
    min_similarity_for_matching : float, optional
        Minimum similarity required for a transcript to be considered a match (default: 0.4)
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)
    verbose : bool, optional
        Whether to print detailed progress information (default: True)
    return_plotting_data : bool, optional
        Whether to return the plotting-ready long format data (default: True)

    Returns
    --------
    tuple or pd.DataFrame
        If return_plotting_data=True: returns (results_df, plotting_df)
        If return_plotting_data=False: returns results_df only

        plotting_df includes rows for both major and minor isoforms:
            - isoform_rank: "major" or "minor"
            - isoform_id: reference ID for major, actual transcript ID for minor
            - <additional_info>
    """
    import pandas as pd
    import numpy as np
    from statsmodels.stats.multitest import multipletests
    from anndata import AnnData

    # First run the main analysis to get statistical results
    results_df = test_isoform_DIU_between_alleles_by_structure(
        adata, layer, test_condition, structure_similarity_threshold,
        min_similarity_for_matching, inplace, verbose
    )

    if not return_plotting_data:
        return results_df

    if len(results_df) == 0:
        print("No results to format for plotting")
        return results_df, pd.DataFrame()

    # Now create the enhanced plotting format with major and minor isoforms
    if verbose:
        print("Creating plotting data with major and minor isoforms from same haplotype...")

    plotting_data = []

    # Get necessary data
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    counts = adata.layers[layer]
    synt_ids = adata.var["Synt_id"]
    functional_annotation = adata.var["functional_annotation"]
    haplotypes = adata.var["haplotype"]
    transcript_ids = adata.var_names
    exon_lengths_dict = adata.uns['exon_lengths']

    # Get gene_id if available
    gene_ids = adata.var.get("gene_id", pd.Series([None] * len(adata.var), index=adata.var_names))

    # Get condition indices
    if test_condition == "all":
        condition_indices = np.arange(counts.shape[0])
        condition_values = adata.obs['condition'].values
    else:
        condition_indices = np.where(adata.obs['condition'] == test_condition)[0]
        condition_values = [test_condition] * len(condition_indices)

    # Get sample names for replicates
    sample_names = adata.obs_names[condition_indices]

    # Process each syntelog from the results
    for _, result in results_df.iterrows():
        synt_id = result['Synt_id']
        reference_transcript_id = result['reference_transcript']

        # Get all transcripts for this syntelog
        synt_mask = synt_ids == synt_id
        synt_indices = np.where(synt_mask)[0]

        if len(synt_indices) == 0:
            continue

        synt_haplotypes = haplotypes.iloc[synt_indices]
        unique_haplotypes = synt_haplotypes.dropna().unique()

        # Get the reference structure
        reference_structure = exon_lengths_dict.get(reference_transcript_id, [])

        # Find which haplotype the reference transcript belongs to
        reference_haplotype = None
        for idx in synt_indices:
            if transcript_ids[idx] == reference_transcript_id:
                reference_haplotype = haplotypes.iloc[idx]
                break

        if reference_haplotype is None:
            if verbose:
                print(f"Warning: Could not find haplotype for reference transcript {reference_transcript_id}")
            continue

        # Find major and minor isoforms from the reference haplotype
        ref_hap_mask = synt_haplotypes == reference_haplotype
        ref_hap_indices_local = np.where(ref_hap_mask)[0]
        ref_hap_indices_global = synt_indices[ref_hap_indices_local]

        # Calculate expression for all transcripts in the reference haplotype
        ref_hap_expressions = []
        for idx in ref_hap_indices_global:
            transcript_id = transcript_ids[idx]
            total_expr = np.sum(counts[np.ix_(condition_indices, [idx])])
            structure = exon_lengths_dict.get(transcript_id, [])
            similarity = _calculate_structure_similarity(reference_structure, structure)

            ref_hap_expressions.append({
                'transcript_id': transcript_id,
                'transcript_idx': idx,
                'total_expression': total_expr,
                'structure': structure,
                'similarity_score': similarity
            })

        # Sort by expression (descending)
        ref_hap_expressions.sort(key=lambda x: x['total_expression'], reverse=True)

        # Get major (reference) and minor (second in same haplotype) isoforms
        major_isoform = None
        minor_isoform = None

        # Find the reference transcript in the list
        for i, transcript_data in enumerate(ref_hap_expressions):
            if transcript_data['transcript_id'] == reference_transcript_id:
                major_isoform = transcript_data
                # Get the second most expressed transcript from same haplotype
                if len(ref_hap_expressions) > 1:
                    # Find next transcript that isn't the reference
                    for j, other_transcript in enumerate(ref_hap_expressions):
                        if other_transcript['transcript_id'] != reference_transcript_id:
                            minor_isoform = other_transcript
                            break
                break

        if major_isoform is None:
            continue

        # For each haplotype, create plotting data
        for hap in unique_haplotypes:
            # Get transcripts for this haplotype
            hap_mask = synt_haplotypes == hap
            hap_indices_local = np.where(hap_mask)[0]
            hap_indices_global = synt_indices[hap_indices_local]

            if len(hap_indices_global) == 0:
                continue

            # Find the matching transcript for this haplotype (for major isoform)
            # This is the transcript selected by the DIU analysis
            major_transcript_col = f'transcript_id_{hap}'
            if major_transcript_col not in result or pd.isna(result[major_transcript_col]):
                continue

            selected_transcript_id = result[major_transcript_col]

            # Find the transcript index
            try:
                major_transcript_idx = transcript_ids.get_loc(selected_transcript_id)
            except KeyError:
                continue

            # Get the gene_id for this transcript
            gene_id = gene_ids.get(selected_transcript_id, None)

            # Get major isoform data
            major_structure = exon_lengths_dict.get(selected_transcript_id, [])
            major_similarity = _calculate_structure_similarity(reference_structure, major_structure)

            # For minor isoform: if this is the reference haplotype, use the minor from reference haplotype
            # If this is a different haplotype, find the second most expressed transcript in this haplotype
            if hap == reference_haplotype and minor_isoform is not None:
                # Use the pre-identified minor isoform from reference haplotype
                minor_transcript_idx = minor_isoform['transcript_idx']
                minor_transcript_id = minor_isoform['transcript_id']
                minor_structure = minor_isoform['structure']
                minor_similarity = minor_isoform['similarity_score']
            else:
                # Find second most expressed in this haplotype
                hap_expressions = []
                for idx in hap_indices_global:
                    transcript_id = transcript_ids[idx]
                    if transcript_id == selected_transcript_id:  # Skip the major isoform
                        continue
                    total_expr = np.sum(counts[np.ix_(condition_indices, [idx])])
                    structure = exon_lengths_dict.get(transcript_id, [])
                    similarity = _calculate_structure_similarity(reference_structure, structure)

                    hap_expressions.append({
                        'transcript_id': transcript_id,
                        'transcript_idx': idx,
                        'total_expression': total_expr,
                        'structure': structure,
                        'similarity_score': similarity
                    })

                if hap_expressions:
                    # Sort by expression and take the most expressed (excluding major)
                    hap_expressions.sort(key=lambda x: x['total_expression'], reverse=True)
                    minor_data = hap_expressions[0]
                    minor_transcript_idx = minor_data['transcript_idx']
                    minor_transcript_id = minor_data['transcript_id']
                    minor_structure = minor_data['structure']
                    minor_similarity = minor_data['similarity_score']
                else:
                    # No minor isoform available for this haplotype
                    minor_transcript_idx = None
                    minor_transcript_id = None
                    minor_structure = []
                    minor_similarity = 0.0

            # Create plotting data for both major and minor isoforms
            isoforms_to_process = [
                ('major', {
                    'transcript_idx': major_transcript_idx,
                    'transcript_id': selected_transcript_id,
                    'structure': major_structure,
                    'similarity_score': major_similarity
                }, reference_transcript_id)  # Use reference ID for major across ALL haplotypes
            ]

            if minor_transcript_idx is not None:
                # For minor isoform: use the minor transcript ID from reference haplotype for all haplotypes
                minor_reference_id = minor_isoform['transcript_id'] if minor_isoform else minor_transcript_id
                isoforms_to_process.append(
                    ('minor', {
                        'transcript_idx': minor_transcript_idx,
                        'transcript_id': minor_transcript_id,
                        'structure': minor_structure,
                        'similarity_score': minor_similarity
                    }, minor_reference_id)  # Use reference minor ID for all haplotypes
                )

            # Process each isoform (major and minor)
            for isoform_rank, isoform_data, isoform_id in isoforms_to_process:
                transcript_idx = isoform_data['transcript_idx']
                transcript_id = isoform_data['transcript_id']
                structure = isoform_data['structure']
                similarity_score = isoform_data['similarity_score']

                # For each replicate/sample
                for sample_idx, (condition_idx, sample_name, condition) in enumerate(
                    zip(condition_indices, sample_names, condition_values)
                ):

                    # Get isoform counts
                    isoform_counts = counts[condition_idx, transcript_idx]

                    # Get total counts for this haplotype
                    hap_total_counts = np.sum(counts[condition_idx, hap_indices_global])

                    # Calculate ratio
                    isoform_ratio = isoform_counts / hap_total_counts if hap_total_counts > 0 else 0.0

                    # Add row to plotting data
                    plotting_data.append({
                        'Synt_id': synt_id,
                        'gene_id': gene_id,  # Actual gene_id from adata.var['gene_id']
                        'haplotype': hap,
                        'sample': sample_name,
                        'isoform_rank': isoform_rank,  # "major" or "minor"
                        'isoform_id': isoform_id,  # Reference ID for major, actual ID for minor
                        'transcript_id': transcript_id,  # Actual transcript used
                        'isoform_counts': int(isoform_counts),
                        'total_counts': int(hap_total_counts),
                        'isoform_ratio': float(isoform_ratio),
                        'similarity_score': float(similarity_score),
                        'structure': ','.join(map(str, structure)),
                        'reference_structure': ','.join(map(str, reference_structure)),
                        'reference_haplotype': reference_haplotype,  # New column
                        'is_reference_haplotype': hap == reference_haplotype,  # New column
                        'condition': condition,
                        'p_value': result['p_value'],
                        'FDR': result['FDR'],
                        'ratio_difference': result['ratio_difference'],
                        'matching_quality': result['matching_quality'],
                        'significance': 'significant' if (result['FDR'] < 0.05 and result['ratio_difference'] > 0.2) else 'not_significant',
                        'n_haplotypes': result['n_haplotypes'],
                        'mean_similarity_score': result['mean_similarity_score']
                    })

    plotting_df = pd.DataFrame(plotting_data)

    if verbose:
        print(f"Created plotting data with {len(plotting_df)} rows")
        print(f"Covering {plotting_df['Synt_id'].nunique()} syntelogs")
        print(f"Including both major and minor isoforms from same reference haplotype")

        # Show breakdown by isoform rank
        rank_counts = plotting_df['isoform_rank'].value_counts()
        print(f"\nIsoform rank distribution:")
        for rank, count in rank_counts.items():
            print(f"  {rank}: {count} measurements")

        # Show reference haplotype info
        ref_hap_counts = plotting_df['is_reference_haplotype'].value_counts()
        print(f"\nReference haplotype distribution:")
        print(f"  Reference haplotype: {ref_hap_counts.get(True, 0)} measurements")
        print(f"  Other haplotypes: {ref_hap_counts.get(False, 0)} measurements")

        # Show example data
        if len(plotting_df) > 0:
            print(f"\nExample rows:")
            example_cols = ['Synt_id', 'haplotype', 'isoform_rank', 'isoform_ratio', 'is_reference_haplotype']
            print(plotting_df[example_cols].head(6))

    return results_df, plotting_df

def test_isoform_DIU_between_alleles_by_structure(
    adata, layer="unique_counts", test_condition="control",
    structure_similarity_threshold=0.6, min_similarity_for_matching=0.8,
    inplace=True, verbose=True
):
    """
    Test if alleles have different isoform usage using expression-aware adaptive structure-based matching.

    This function handles cases where haplotypes have different numbers of isoforms by:
    1. Finding the most expressed isoform structure across all haplotypes
    2. For each haplotype, finding the most similar transcript structure
    3. When multiple transcripts have similar structures, choosing the most expressed one
    4. Comparing usage of the best matching transcripts across haplotypes

    Parameters
    -----------
    adata : AnnData
        AnnData object containing expression data with exon structure information
    layer : str, optional
        Layer containing count data (default: "unique_counts")
    test_condition : str, optional
        Variable column name containing condition for testing within (default: "control")
    structure_similarity_threshold : float, optional
        Minimum similarity score (0-1) for considering structures as similar (default: 0.6)
    min_similarity_for_matching : float, optional
        Minimum similarity required for a transcript to be considered a match (default: 0.4)
    inplace : bool, optional
        Whether to modify the input AnnData object or return a copy (default: True)
    verbose : bool, optional
        Whether to print detailed progress information (default: True)

    Returns
    --------
    pd.DataFrame
        Results of statistical tests for each syntelog with additional columns:
        - expression_based_choices: number of haplotypes where expression was used to choose
        - total_candidates_found: total number of candidate transcripts found
        - expression_levels: dictionary of expression levels for selected transcripts
    """
    import pandas as pd
    import numpy as np
    from statsmodels.stats.multitest import multipletests
    from isotools._transcriptome_stats import betabinom_lr_test
    from anndata import AnnData

    # Validate inputs
    if not isinstance(adata, AnnData):
        raise ValueError("Input adata must be an AnnData object")

    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in AnnData object")

    required_structure_cols = ['exon_structure', 'n_exons']
    missing_cols = [col for col in required_structure_cols if col not in adata.var.columns]
    if missing_cols:
        raise ValueError(f"Missing required structure columns: {missing_cols}. "
                        "Please run add_exon_structure() first.")

    if 'exon_lengths' not in adata.uns:
        raise ValueError("'exon_lengths' not found in adata.uns. Please run add_exon_structure() first.")

    if not inplace:
        adata = adata.copy()

    counts = adata.layers[layer].copy()

    if "Synt_id" not in adata.var:
        raise ValueError("'Synt_id' not found in adata.var")
    synt_ids = adata.var["Synt_id"]

    if "haplotype" not in adata.var:
        raise ValueError("'haplotype' not found in adata.var")
    haplotypes = adata.var["haplotype"]

    if not adata.var_names.any():
        raise ValueError("'transcript_id' not found in adata.var_names")
    transcript_ids = adata.var_names

    if test_condition not in adata.obs['condition'].unique() and test_condition != "all":
        raise ValueError(f"Condition '{test_condition}' not found in adata.obs['condition']")

    exon_lengths_dict = adata.uns['exon_lengths']

    unique_synt_ids = np.unique(synt_ids)
    unique_synt_ids = unique_synt_ids[~pd.isna(unique_synt_ids)]
    unique_synt_ids = unique_synt_ids[unique_synt_ids != 0]

    results = []
    total_syntelogs = len(unique_synt_ids)
    processed = 0
    successful_matches = 0
    failed_matches = 0
    expression_based_choices = 0

    if verbose:
        print(f"Processing {total_syntelogs} syntelogs for expression-aware adaptive structure-based DIU analysis...")

    # Process each syntelog
    for synt_id in unique_synt_ids:
        processed += 1
        if verbose and processed % 100 == 0:
            print(f"Processing syntelog {processed}/{total_syntelogs}")

        synt_mask = synt_ids == synt_id
        synt_indices = np.where(synt_mask)[0]

        if len(synt_indices) == 0:
            continue

        synt_haplotypes = haplotypes.iloc[synt_indices]
        unique_haplotypes = synt_haplotypes.dropna().unique()

        if len(unique_haplotypes) < 2:
            continue

        if test_condition == "all":
            condition_indices = np.arange(counts.shape[0])
        else:
            condition_indices = np.where(adata.obs['condition'] == test_condition)[0]

        # Use improved reference structure selection
        reference_structure, reference_transcript = _find_reference_structure(
            synt_indices, transcript_ids, exon_lengths_dict, counts, condition_indices
        )

        if reference_structure is None:
            continue

        # Use expression-aware matching
        haplotype_matches = _find_best_matches_per_haplotype(
            synt_indices, unique_haplotypes, synt_haplotypes, transcript_ids,
            exon_lengths_dict, reference_structure, min_similarity_for_matching,
            counts, condition_indices, verbose and processed <= 5
        )

        if len(haplotype_matches) < 2:
            failed_matches += 1
            continue

        successful_matches += 1

        # Count cases where expression was used to break ties
        expression_choices = sum(1 for match in haplotype_matches.values()
                               if match.get('n_candidates_best_similarity', 1) > 1)
        expression_based_choices += expression_choices

        # Perform statistical test with enhanced results
        try:
            test_results = _perform_statistical_test(
                haplotype_matches, synt_id, reference_structure, reference_transcript
            )
            if test_results:
                results.append(test_results)
        except Exception as e:
            if verbose and processed <= 5:
                print(f"  Error testing syntelog {synt_id}: {str(e)}")
            continue

    if verbose:
        print(f"\nCompleted processing:")
        print(f"  Total syntelogs: {total_syntelogs}")
        print(f"  Successful matches: {successful_matches}")
        print(f"  Failed matches: {failed_matches}")
        print(f"  Expression-based choices: {expression_based_choices}")
        print(f"  Results generated: {len(results)}")

    # Convert results to DataFrame and apply multiple testing correction
    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        results_df['p_value'] = results_df['p_value'].fillna(1)
        results_df['FDR'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df = results_df.sort_values('p_value')

        significant_results = results_df[
            (results_df['FDR'] < 0.05) &
            (results_df['ratio_difference'] > 0.2)
        ]


        print(f"\nFound {len(significant_results)} from {len(results_df)} syntelogs with "
              f"significantly different isoform usage between alleles (FDR < 0.05, ratio diff > 0.2)")

        if inplace:
            adata.uns['adaptive_structure_diu_test_v2'] = results_df
    else:
        if verbose:
            print("No results found")

    return results_df


def _calculate_structure_similarity(structure1, structure2):
    """
    Enhanced structure similarity calculation that handles different numbers of exons better.
    """
    import numpy as np

    if len(structure1) == 0 and len(structure2) == 0:
        return 1.0

    if len(structure1) == 0 or len(structure2) == 0:
        return 0.0

    # Number of exons similarity (penalize differences less harshly)
    max_exons = max(len(structure1), len(structure2))
    min_exons = min(len(structure1), len(structure2))
    exon_count_sim = min_exons / max_exons

    # Total length similarity
    total1, total2 = sum(structure1), sum(structure2)
    if total1 == 0 and total2 == 0:
        length_sim = 1.0
    elif total1 == 0 or total2 == 0:
        length_sim = 0.0
    else:
        length_sim = min(total1, total2) / max(total1, total2)

    # Pattern similarity (compare overlapping exons)
    pattern_sim = 0.0
    if min_exons > 0:
        overlapping_diffs = []
        for i in range(min_exons):
            e1, e2 = structure1[i], structure2[i]
            if max(e1, e2) > 0:
                diff = abs(e1 - e2) / max(e1, e2)
                overlapping_diffs.append(1 - diff)

        if overlapping_diffs:
            pattern_sim = np.mean(overlapping_diffs)

    # Weighted combination (adjusted for adaptive matching)
    similarity = (0.3 * exon_count_sim + 0.3 * length_sim + 0.4 * pattern_sim)

    return max(0.0, min(1.0, similarity))

def _find_reference_structure(synt_indices, transcript_ids, exon_lengths_dict,
                                               counts, condition_indices):
    """
    Find the reference structure more intelligently by considering both expression and structure complexity.
    """
    import numpy as np

    transcript_candidates = []

    for idx in synt_indices:
        transcript_id = transcript_ids[idx]
        structure = exon_lengths_dict.get(transcript_id, [])

        if not structure:
            continue

        # Calculate total expression for this transcript
        transcript_counts = counts[np.ix_(condition_indices, [idx])]
        total_expression = np.sum(transcript_counts)

        # Calculate structure complexity metrics
        n_exons = len(structure)
        total_length = sum(structure)
        length_variance = np.var(structure) if len(structure) > 1 else 0

        transcript_candidates.append({
            'transcript_id': transcript_id,
            'transcript_idx': idx,
            'structure': structure,
            'total_expression': total_expression,
            'n_exons': n_exons,
            'total_length': total_length,
            'length_variance': length_variance
        })

    if not transcript_candidates:
        return None, None

    # Sort by expression (descending)
    transcript_candidates.sort(key=lambda x: x['total_expression'], reverse=True)

    # Take the top 20% by expression or at least top 3 transcripts
    top_n = max(3, len(transcript_candidates) // 5)
    top_expressed = transcript_candidates[:top_n]

    # Among highly expressed transcripts, prefer those with moderate complexity
    def complexity_score(candidate):
        exon_score = 1.0
        if candidate['n_exons'] < 2:
            exon_score = 0.5  # Single exon is less representative
        elif candidate['n_exons'] > 8:
            exon_score = 0.8  # Very complex might be outlier

        length_score = 1.0
        if candidate['total_length'] < 1000:
            length_score = 0.7
        elif candidate['total_length'] > 30000:
            length_score = 0.8

        return exon_score * length_score

    # Score top expressed transcripts by complexity
    for candidate in top_expressed:
        candidate['complexity_score'] = complexity_score(candidate)

    # Choose the one with best complexity score among top expressed
    best_candidate = max(top_expressed, key=lambda x: x['complexity_score'])

    return best_candidate['structure'], best_candidate['transcript_id']


def _find_best_matches_per_haplotype(
    synt_indices, unique_haplotypes, synt_haplotypes, transcript_ids,
    exon_lengths_dict, reference_structure, min_similarity,
    counts, condition_indices, verbose=False
):
    """
    Find the best matching transcript for each haplotype, considering both structure similarity
    and expression levels. When multiple transcripts have similar structures, choose the most expressed one.
    """
    import numpy as np

    haplotype_matches = {}

    if verbose:
        print(f"    Reference structure: {','.join(map(str, reference_structure))}")

    for hap in unique_haplotypes:
        # Get transcript indices for this haplotype
        hap_mask = synt_haplotypes == hap
        hap_indices_local = np.where(hap_mask)[0]
        hap_indices_global = synt_indices[hap_indices_local]

        # Find all transcripts that meet minimum similarity threshold
        candidate_matches = []

        for idx in hap_indices_global:
            transcript_id = transcript_ids[idx]
            structure = exon_lengths_dict.get(transcript_id, [])

            if not structure:
                continue

            similarity = _calculate_structure_similarity(reference_structure, structure)

            if similarity >= min_similarity:
                # Calculate expression for this transcript
                transcript_counts = counts[np.ix_(condition_indices, [idx])]
                total_expression = np.sum(transcript_counts)

                candidate_matches.append({
                    'transcript_id': transcript_id,
                    'transcript_idx': idx,
                    'similarity_score': similarity,
                    'structure': structure,
                    'total_expression': total_expression
                })

        if not candidate_matches:
            if verbose:
                print(f"    {hap}: No suitable matches found")
            continue

        # Group candidates by similarity score (rounded to avoid floating point issues)
        similarity_groups = {}
        for candidate in candidate_matches:
            sim_rounded = round(candidate['similarity_score'], 3)
            if sim_rounded not in similarity_groups:
                similarity_groups[sim_rounded] = []
            similarity_groups[sim_rounded].append(candidate)

        # Find the group with highest similarity
        best_similarity = max(similarity_groups.keys())
        best_group = similarity_groups[best_similarity]

        # if the next best similarity is within 0.1, consider it as well
        next_best_similarity = max(
            [sim for sim in similarity_groups.keys() if sim < best_similarity],
            default=None
        )
        # append
        if next_best_similarity is not None and best_similarity - next_best_similarity < 0.05:
            best_group.extend(similarity_groups[next_best_similarity])

        # Within the best similarity group, choose the most expressed transcript
        best_match = max(best_group, key=lambda x: x['total_expression'])

        # Within the best similarity group, choose the most expressed transcript
        best_match = max(best_group, key=lambda x: x['total_expression'])

        # Get counts for this transcript
        isoform_counts = counts[np.ix_(condition_indices, [best_match['transcript_idx']])][:, 0]

        # Get total counts for all transcripts of this haplotype in this syntelog
        hap_total_counts = np.sum(counts[np.ix_(condition_indices, hap_indices_global)], axis=1)

        haplotype_matches[hap] = {
            'transcript_id': best_match['transcript_id'],
            'transcript_idx': best_match['transcript_idx'],
            'similarity_score': best_match['similarity_score'],
            'structure': best_match['structure'],
            'total_expression': best_match['total_expression'],
            'isoform_counts': isoform_counts,
            'total_counts': hap_total_counts,
            'n_candidates': len(candidate_matches),
            'n_candidates_best_similarity': len(best_group)
        }

        if verbose:
            structure_str = ','.join(map(str, best_match['structure']))
            expr_info = f", expr: {best_match['total_expression']:.0f}"
            candidates_info = f", {len(candidate_matches)} candidates"
            if len(best_group) > 1:
                candidates_info += f" ({len(best_group)} with best similarity)"

            print(f"    {hap}: {best_match['transcript_id']} "
                  f"(similarity: {best_match['similarity_score']:.3f}{expr_info}{candidates_info})")
            print(f"         Structure: {structure_str}")

    return haplotype_matches


def _perform_statistical_test(haplotype_matches, synt_id, reference_structure, reference_transcript):
    """
    Statistical test that includes information about expression-based choices.
    """
    import numpy as np

    if len(haplotype_matches) < 2:
        return None

    # Calculate ratios for each haplotype
    haplotype_ratios = {}
    expression_levels = {}

    for hap, data in haplotype_matches.items():
        valid_samples = data['total_counts'] > 0
        if np.sum(valid_samples) > 0:
            ratios = np.zeros_like(data['total_counts'], dtype=float)
            ratios[valid_samples] = (data['isoform_counts'][valid_samples] /
                                   data['total_counts'][valid_samples])
            haplotype_ratios[hap] = np.mean(ratios)
        else:
            haplotype_ratios[hap] = 0.0

        expression_levels[hap] = data['total_expression']

    if len(haplotype_ratios) < 2:
        return None

    sorted_haps = sorted(haplotype_ratios.keys(), key=lambda x: haplotype_ratios[x])
    min_hap = sorted_haps[0]
    max_hap = sorted_haps[-1]

    if haplotype_ratios[min_hap] == haplotype_ratios[max_hap]:
        return None

    # Prepare data for statistical test
    min_hap_data = haplotype_matches[min_hap]
    max_hap_data = haplotype_matches[max_hap]

    allele_counts = [min_hap_data['isoform_counts'], max_hap_data['isoform_counts']]
    condition_total = [min_hap_data['total_counts'], max_hap_data['total_counts']]

    if np.any([np.sum(ct) == 0 for ct in condition_total]):
        return None

    # Run statistical test
    try:
        from isotools._transcriptome_stats import betabinom_lr_test
        test_result = betabinom_lr_test(allele_counts, condition_total)
        p_value, ratio_stats = test_result[0], test_result[1]

        ratio_difference = abs(ratio_stats[0] - ratio_stats[2]) if len(ratio_stats) >= 3 else \
                          abs(haplotype_ratios[min_hap] - haplotype_ratios[max_hap])
    except Exception:
        return None

    # Assess matching quality
    similarities = [data['similarity_score'] for data in haplotype_matches.values()]
    mean_similarity = np.mean(similarities)
    min_similarity = np.min(similarities)

    if mean_similarity >= 0.8:
        matching_quality = "excellent"
    elif mean_similarity >= 0.6:
        matching_quality = "good"
    elif mean_similarity >= 0.4:
        matching_quality = "fair"
    else:
        matching_quality = "poor"

    # Calculate additional metrics
    total_candidates = sum(data.get('n_candidates', 1) for data in haplotype_matches.values())
    expression_based_choices = sum(1 for data in haplotype_matches.values()
                                 if data.get('n_candidates_best_similarity', 1) > 1)

    # Build enhanced result dictionary
    result = {
        'Synt_id': synt_id,
        'reference_structure': ','.join(map(str, reference_structure)),
        'reference_transcript': reference_transcript,
        'min_ratio_haplotype': min_hap,
        'max_ratio_haplotype': max_hap,
        'min_ratio_transcript_id': min_hap_data['transcript_id'],
        'max_ratio_transcript_id': max_hap_data['transcript_id'],
        'p_value': p_value,
        'ratio_difference': ratio_difference,
        'n_haplotypes': len(haplotype_matches),
        'matching_quality': matching_quality,
        'mean_similarity_score': mean_similarity,
        'min_similarity_score': min_similarity,
        'total_candidates_found': total_candidates,
        'expression_based_choices': expression_based_choices,
        f'ratio_{min_hap}_mean': haplotype_ratios[min_hap],
        f'ratio_{max_hap}_mean': haplotype_ratios[max_hap]
    }

    # Add detailed information for each haplotype
    similarity_scores = {}
    candidate_counts = {}

    for hap, data in haplotype_matches.items():
        result[f'ratio_{hap}_mean'] = haplotype_ratios[hap]
        result[f'transcript_id_{hap}'] = data['transcript_id']
        result[f'similarity_{hap}'] = data['similarity_score']
        result[f'expression_{hap}'] = data['total_expression']
        result[f'n_candidates_{hap}'] = data.get('n_candidates', 1)

        similarity_scores[hap] = data['similarity_score']
        candidate_counts[hap] = data.get('n_candidates', 1)

    result['similarity_scores'] = similarity_scores
    result['candidate_counts'] = candidate_counts
    result['expression_levels'] = expression_levels
    result['all_haplotypes'] = list(haplotype_matches.keys())

    return result
