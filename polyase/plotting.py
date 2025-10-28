import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Union, Tuple

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Union, Tuple

def plot_allelic_ratios(
    adata,
    synteny_category: str,
    sample: Union[str, List[str]] = "all",
    multimapping_threshold: float = 0.5,
    ratio_type: str = "both",
    bins: int = 30,
    figsize: Tuple[int, int] = (12, 6),
    kde: bool = True,
    save_path: Optional[str] = None
    ):
    """
    Plot allelic ratios for transcripts in a specific synteny category.
    Filters out Synt_ids where all allelic ratios are 0 for a given sample.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing transcript data
    synteny_category : str
        Synteny category to filter for
    sample : str or List[str], default="all"
        Sample(s) to plot. If list, will plot each sample separately
    multimapping_threshold : float, default=0.5
        Threshold for high multimapping ratio
    ratio_type : str, default="both"
        Type of ratio to plot: "unique", "em", or "both"
    bins : int, default=30
        Number of bins for the histogram
    figsize : tuple, default=(12, 6)
        Figure size (width, height) in inches
    kde : bool, default=True
        Whether to show KDE curve on histogram
    save_path : str, optional
        Path to save the plot. If None, plot is shown but not saved

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the plot(s)
    """
    # Validate parameters
    valid_ratio_types = ["unique", "em", "both"]
    if ratio_type not in valid_ratio_types:
        raise ValueError(f"ratio_type must be one of {valid_ratio_types}")

    # Make sample always a list for consistent processing
    if isinstance(sample, str):
        samples = [sample]
    else:
        samples = sample

    # Ensure all samples exist in var
    for s in samples:
        if s not in adata.var and s != "all":
            raise ValueError(f"Sample '{s}' not found in adata.var")

    # Filter the data for the specific synteny category
    filtered_data = adata[:,adata.var['synteny_category'] == synteny_category].copy()

    if len(filtered_data) == 0:
        print(f"No data found for synteny category: {synteny_category}")
        return None

    # Add a tag for high multimapping ratio
    filtered_data.var['ambiguous_counts'] = np.where(
    filtered_data.var['multimapping_ratio'] > multimapping_threshold,
    'high',
    'low')

    # Create figure with appropriate number of subplots
    if ratio_type == "both":
        n_ratio_plots = 2
    else:
        n_ratio_plots = 1

    n_sample_plots = len(samples)
    total_plots = n_ratio_plots * n_sample_plots

    # Calculate grid dimensions
    if total_plots <= 2:
        n_rows, n_cols = 1, total_plots
    else:
        n_cols = min(3, total_plots)
        n_rows = (total_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    plot_idx = 0

    # Create plots for each sample and ratio type
    for sample_idx, current_sample in enumerate(samples):
        # Extract data for current sample
        if current_sample != "all":
            sample_indices = np.where(filtered_data.obs_names == current_sample)[0]
        else:
            sample_indices = np.arange(filtered_data.shape[0])

        # Determine which ratio types to plot
        ratio_types_to_plot = ["unique", "em"] if ratio_type == "both" else [ratio_type]

        for rt in ratio_types_to_plot:
            if rt == "unique":
                layer_name = "allelic_ratio_unique_counts"
                title_suffix = "Unique Counts"
            else:
                layer_name = "allelic_ratio_em_counts"
                title_suffix = "EM Counts"

            # Check if layer exists
            if layer_name not in filtered_data.layers:
                print(f"Warning: Layer '{layer_name}' not found, skipping {title_suffix}")
                continue

            # Extract allelic ratios for the current sample and layer
            allelic_ratios = filtered_data.layers[layer_name][sample_indices]

            # Filter out Synt_ids where all ratios are 0 for this sample
            valid_synt_mask = []
            synt_ids = filtered_data.var['Synt_id'].values
            unique_synt_ids = np.unique(synt_ids)

            for synt_id in unique_synt_ids:
                # Get indices for this Synt_id
                synt_mask = synt_ids == synt_id
                synt_indices = np.where(synt_mask)[0]

                # Get ratios for this Synt_id and sample
                synt_ratios = allelic_ratios[:, synt_indices]

                # Check if all ratios for this Synt_id are 0 (excluding NaN)
                synt_ratios_clean = synt_ratios[~np.isnan(synt_ratios)]
                if len(synt_ratios_clean) > 0 and not np.all(synt_ratios_clean == 0):
                    # This Synt_id has non-zero ratios, keep it
                    valid_synt_mask.extend([True] * np.sum(synt_mask))
                else:
                    # This Synt_id has all zero ratios, exclude it
                    valid_synt_mask.extend([False] * np.sum(synt_mask))
                    print(f"Excluding Synt_id {synt_id} for sample {current_sample} ({rt}): all ratios are 0")

            valid_synt_mask = np.array(valid_synt_mask)

            # Apply the mask to filter data
            if not np.any(valid_synt_mask):
                print(f"No valid Synt_ids for sample {current_sample}, ratio type {rt}")
                continue

            filtered_ratios = allelic_ratios[:, valid_synt_mask]
            filtered_ambiguous = filtered_data.var['ambiguous_counts'].values[valid_synt_mask]

            # Create DataFrame for plotting
            plot_data = pd.DataFrame({
                'allelic_ratio': filtered_ratios.flatten(order='F'),
                'ambiguous_counts': np.repeat(filtered_ambiguous, len(sample_indices))
            })

            # Drop NaN values
            plot_data = plot_data.dropna()

            if len(plot_data) == 0:
                print(f"No valid data for sample {current_sample}, ratio type {rt}")
                continue

            # Plot histogram
            sns.histplot(
                data=plot_data,
                x='allelic_ratio',
                hue='ambiguous_counts',
                kde=kde,
                bins=bins,
                palette={'low':'green', 'high':'grey'},
                ax=axes[plot_idx]
            )

            axes[plot_idx].set_title(f"{title_suffix} Allelic Ratio ({current_sample})")
            axes[plot_idx].set_xlabel('Allelic Ratio')
            axes[plot_idx].set_ylabel('Count')
            axes[plot_idx].grid(True, linestyle='--', alpha=0.7)
            axes[plot_idx].set_xticks([0, 0.25, 0.5, 0.75, 1])

            plot_idx += 1

    # Remove any empty subplots
    for i in range(plot_idx, len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()

    # Save the plot if a path is provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig

def plot_top_differential_syntelogs(results_df, n=5, figsize=(16, 12), palette=None, jitter=0.2, alpha=0.7,
                                   ylim=None, sort_by='p_value', output_file=None, sig_threshold=0.05,
                                   difference_threshold=0.05, sig_color='red', plot_type='ratios'):
    """
    Plot the top n syntelogs with differential allelic ratios or CPM values in a grid layout (6 plots per row).
    Syntelogs with significant differences will have their titles highlighted in red.

    Parameters
    -----------
    results_df : pd.DataFrame
        Results dataframe from test_allelic_ratios function
    n : int, optional
        Number of top syntelogs to plot (default: 5)
    figsize : tuple, optional
        Figure size as (width, height) in inches (default: (16, 12))
    palette : dict or None, optional
        Color palette for conditions (default: None, uses seaborn defaults)
    jitter : float, optional
        Amount of jitter for strip plot (default: 0.2)
    alpha : float, optional
        Transparency of points (default: 0.7)
    ylim : tuple, optional
        Y-axis limits (default: None, auto-determined based on plot_type)
    sort_by : str, optional
        Column to sort results by ('p_value', 'FDR', or 'ratio_difference') (default: 'p_value')
    output_file : str, optional
        Path to save the figure (default: None, displays figure but doesn't save)
    sig_threshold : float, optional
        Significance threshold for p-value or FDR (default: 0.05)
    difference_threshold : float, optional
        Ratio difference threshold for significance (default: 0.05)
    sig_color : str, optional
        Color for titles of syntelogs with significant differences (default: 'red')
    plot_type : str, optional
        What to plot: 'ratios' for allelic ratios or 'cpm' for CPM values (default: 'ratios')

    Returns
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import math

    if len(results_df) == 0:
        print("No results to plot")
        return None

    # Validate plot_type parameter
    if plot_type not in ['ratios', 'cpm']:
        print(f"Invalid plot_type '{plot_type}'. Using 'ratios' instead.")
        plot_type = 'ratios'

    # Validate sort_by parameter
    if sort_by not in ['p_value', 'FDR', 'ratio_difference']:
        print(f"Invalid sort_by parameter '{sort_by}'. Using 'p_value' instead.")
        sort_by = 'p_value'

    # Ensure FDR column exists
    if 'FDR' not in results_df.columns and sort_by == 'FDR':
        print("FDR column not found. Using p_value for sorting.")
        sort_by = 'p_value'

    # Ensure ratio_difference column exists
    if 'ratio_difference' not in results_df.columns and sort_by == 'ratio_difference':
        print("ratio_difference column not found. Using p_value for sorting.")
        sort_by = 'p_value'

    ascending_bool = False if sort_by == 'ratio_difference' else True

    # Determine what columns to use based on plot_type
    if plot_type == 'ratios':
        value_prefix = 'ratios'
        y_label = 'Expression Ratio'
        default_ylim = (0, 1)
    else:  # plot_type == 'cpm'
        value_prefix = 'cpm'
        y_label = 'CPM'
        default_ylim = None  # Auto-scale for CPM

    # Set y-limits
    if ylim is None:
        ylim = default_ylim

    # Get the condition names from appropriate columns
    condition_columns = [col for col in results_df.columns if col.startswith(f'{value_prefix}_rep_')]
    if not condition_columns:
        print(f"No {value_prefix} columns found in dataframe")
        # Try to fall back to ratios if CPM was requested but not available
        if plot_type == 'cmp':
            print("Falling back to ratios...")
            condition_columns = [col for col in results_df.columns if col.startswith('ratios_rep_')]
            if condition_columns:
                plot_type = 'ratios'
                value_prefix = 'ratios'
                y_label = 'Expression Ratio'
                ylim = (0, 1) if ylim is None else ylim
            else:
                print("No ratio columns found either")
                return None
        else:
            return None

    conditions = [col.replace(f'{value_prefix}_rep_', '') for col in condition_columns]

    # Get top n syntelogs
    top_syntelogs = results_df.sort_values(sort_by, ascending=ascending_bool).drop_duplicates('Synt_id').head(n)['Synt_id'].unique()
    top_results = results_df[results_df['Synt_id'].isin(top_syntelogs)]

    # Calculate grid dimensions - 6 plots per row
    cols = 6
    rows = math.ceil(len(top_syntelogs) / cols)

    # Create the figure with grid layout
    fig, axes = plt.subplots(rows, cols, figsize=figsize)

    # Convert axes to flattened array for easier indexing
    if rows == 1 and cols == 1:
        axes = np.array([axes])
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    # Plot each syntelog
    for i, synt_id in enumerate(top_syntelogs):
        # Get data for this syntelog
        synt_data = top_results[top_results['Synt_id'] == synt_id].copy()
        synt_data = synt_data.sort_values('allele')

        # Get stats for this syntelog
        p_value = synt_data['p_value'].min()
        fdr = synt_data['FDR'].min() if 'FDR' in synt_data.columns else np.nan
        ratio_difference = (synt_data.loc[synt_data['FDR'] < sig_threshold, 'ratio_difference'].max()
                          if ('ratio_difference' in synt_data.columns and 'FDR' in synt_data.columns)
                          else np.nan)

        # Handle exploding columns with different lengths properly
        explode_cols = [col for col in synt_data.columns if col.startswith(f'{value_prefix}_rep_')]
        exploded_rows = []

        for idx, row in synt_data.iterrows():
            # Get the base row data (non-list columns)
            base_row = {col: row[col] for col in synt_data.columns if col not in explode_cols}

            # Find the maximum length among all explode columns for this row
            max_length = 0
            explode_data = {}
            for col in explode_cols:
                if isinstance(row[col], (list, np.ndarray)):
                    explode_data[col] = row[col]
                    max_length = max(max_length, len(row[col]))
                else:
                    # Handle single values by converting to list
                    explode_data[col] = [row[col]]
                    max_length = max(max_length, 1)

            # Create one row for each replicate
            for rep_idx in range(max_length):
                new_row = base_row.copy()
                for col in explode_cols:
                    if rep_idx < len(explode_data[col]):
                        new_row[col] = explode_data[col][rep_idx]
                    else:
                        # Fill missing values with NaN if one condition has fewer replicates
                        new_row[col] = np.nan
                exploded_rows.append(new_row)

        # Convert to DataFrame
        synt_data_exploded = pd.DataFrame(exploded_rows)

        # Reshape data for seaborn
        id_vars = ['Synt_id', 'allele', 'transcript_id']
        if 'FDR' in synt_data.columns:
            id_vars.append('FDR')
        if 'gene_id' in synt_data.columns:
            id_vars.append('gene_id')
        if 'functional_annotation' in synt_data.columns:
            id_vars.append('functional_annotation')

        synt_data_melted = pd.melt(
            synt_data_exploded,
            id_vars=id_vars,
            value_vars=condition_columns,
            var_name='condition',
            value_name='value'
        )

        # Clean up condition names and remove NaN values
        synt_data_melted['condition'] = synt_data_melted['condition'].str.replace(f'{value_prefix}_rep_', '')
        synt_data_melted = synt_data_melted.dropna(subset=['value'])

        # Create the stripplot
        ax = axes[i]
        sns.stripplot(
            x='allele',
            y='value',
            hue='condition',
            data=synt_data_melted,
            jitter=jitter,
            alpha=alpha,
            palette=palette,
            ax=ax
        )

        # Add mean values as horizontal lines for each allele and condition
        for allele in synt_data['allele'].unique():
            allele_pos = list(synt_data['allele'].unique()).index(allele)

            for j, cond in enumerate(conditions):
                mean_col = f'{value_prefix}_{cond}_mean'
                if mean_col in synt_data.columns:
                    mean_val = synt_data[synt_data['allele'] == allele][mean_col].iloc[0]
                    # Get color from the plot
                    if ax.get_legend() and len(ax.get_legend().get_lines()) > j:
                        line_color = ax.get_legend().get_lines()[j].get_color()
                    else:
                        # Fallback colors if legend not available
                        colors = plt.cm.Set1(np.linspace(0, 1, len(conditions)))
                        line_color = colors[j]

                    ax.hlines(
                        y=mean_val,
                        xmin=allele_pos-0.2,
                        xmax=allele_pos+0.2,
                        colors=line_color,
                        linewidth=2
                    )

        # Determine title components
        fdr_text = f", FDR={fdr:.2e}" if not np.isnan(fdr) else ""

        # Determine if this syntelog has a significant difference
        is_significant = False
        if 'FDR' in synt_data.columns and not np.isnan(fdr):
            is_significant = (fdr <= sig_threshold) & (ratio_difference > difference_threshold)
        else:
            is_significant = p_value <= sig_threshold

        title_color = sig_color if is_significant else 'black'

        # Get functional annotation
        function_annotation_text = "NA"
        if 'functional_annotation' in synt_data_melted.columns:
            function_annotation = synt_data_melted['functional_annotation'].iloc[0]
            if function_annotation is not None:
                function_annotation_text = str(function_annotation)[:20]  # Shorten long descriptions

        transcript_id = synt_data_melted['transcript_id'].iloc[0]

        # Set title with stats and color based on significance
        ax.set_title(f"{transcript_id}\n{function_annotation_text}{fdr_text}", color=title_color)
        ax.set_xlabel('Allele')
        ax.set_ylabel(y_label)

        # Set y-limits
        if ylim is not None:
            ax.set_ylim(ylim)

        # Adjust legend
        if ax.get_legend():
            ax.legend(title='Condition', loc='best')

    # Hide unused subplots if any
    for j in range(len(top_syntelogs), rows * cols):
        axes[j].set_visible(False)

    plt.tight_layout()

    # Save figure if requested
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')

    return fig


def plot_top_differential_isoforms(results_df, n=5, figsize=(16, 12), palette=None, jitter=0.2, alpha=0.7, ylim=(0, 1), sort_by='p_value', output_file=None, sig_threshold=0.05, difference_threshold=0.05, sig_color='red'):
    """
    Plot the top n genes with differential isoform usage in a grid layout (3 plots per row).
    Genes with significant differences will have their titles highlighted in red.

    Parameters
    -----------
    results_df : pd.DataFrame
        Results dataframe from test_allelic_ratios function
    n : int, optional
        Number of top genes to plot (default: 5)
    figsize : tuple, optional
        Figure size as (width, height) in inches (default: (16, 12))
    palette : dict or None, optional
        Color palette for conditions (default: None, uses seaborn defaults)
    jitter : float, optional
        Amount of jitter for strip plot (default: 0.2)
    alpha : float, optional
        Transparency of points (default: 0.7)
    ylim : tuple, optional
        Y-axis limits (default: (0, 1))
    sort_by : str, optional
        Column to sort results by ('p_value', 'FDR', or 'ratio_difference') (default: 'p_value')
    output_file : str, optional
        Path to save the figure (default: None, displays figure but doesn't save)
    sig_threshold : float, optional
        Significance threshold for p-value or FDR (default: 0.05)
    sig_color : str, optional
        Color for titles of genes with significant differences (default: 'red')

    Returns
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import math

    if len(results_df) == 0:
        print("No results to plot")
        return None

    # Validate sort_by parameter
    if sort_by not in ['p_value', 'FDR', 'ratio_difference']:
        print(f"Invalid sort_by parameter '{sort_by}'. Using 'p_value' instead.")
        sort_by = 'p_value'

    # Ensure FDR column exists
    if 'FDR' not in results_df.columns and sort_by == 'FDR':
        print("FDR column not found. Using p_value for sorting.")
        sort_by = 'p_value'

    # Ensure ratio_difference column exists
    if 'ratio_difference' not in results_df.columns and sort_by == 'ratio_difference':
        print("ratio_difference column not found. Using p_value for sorting.")
        sort_by = 'p_value'

    if sort_by == 'ratio_difference':
        ascending_bool = False
    else:
        ascending_bool = True

    # Get the condition names
    condition_columns = [col for col in results_df.columns if col.startswith('ratios_rep_')]
    if not condition_columns:
        print("No ratio columns found in dataframe")
        return None

    conditions = [col.replace('ratios_rep_', '') for col in condition_columns]

    # Get top n genes with lowest sort_by values
    top_genes = results_df.sort_values(sort_by, ascending=ascending_bool).drop_duplicates('gene_id').head(n)['gene_id'].unique()

    # Filter results to include only these syntelogs
    top_results = results_df[results_df['gene_id'].isin(top_genes)]

    # Calculate grid dimensions - 3 plots per row
    cols = 6
    rows = math.ceil(len(top_genes) / cols)

    # Create the figure with grid layout
    fig, axes = plt.subplots(rows, cols, figsize=figsize)

    # Convert axes to flattened array for easier indexing
    if rows == 1 and cols == 1:
        axes = np.array([axes])
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    # Plot each gene
    for i, gene_id in enumerate(top_genes):
        # Get data for this gene
        gene_data = top_results[top_results['gene_id'] == gene_id].copy()

        # Sort by isoform number for better visualization
        gene_data = gene_data.sort_values('isoform_number')

        # Get stats for this syntelog (take first row since they're the same for all replicates)
        p_value = gene_data['p_value'].min()
        fdr = gene_data['FDR'].min() if 'FDR' in gene_data.columns else np.nan
        ratio_difference = gene_data.loc[gene_data['FDR'] < sig_threshold, 'ratio_difference'].max() if ('ratio_difference' in gene_data.columns and 'FDR' in gene_data.columns) else np.nan
        n_isoforms = gene_data['n_isoforms'].iloc[0]

        # Explode the replicate ratio columns
        explode_cols = [col for col in gene_data.columns if col.startswith('ratios_rep_')]
        gene_data_exploded = gene_data.explode(explode_cols)

        # Reshape data for seaborn
        gene_data_melted = pd.melt(
            gene_data_exploded,
            id_vars=['gene_id', 'isoform_number', 'transcript_id', 'FDR'] if 'FDR' in gene_data.columns else ['gene_id', 'isoform_number', 'transcript_id'],
            value_vars=condition_columns,
            var_name='condition',
            value_name='ratio'
        )

        # Clean up condition names
        gene_data_melted['condition'] = gene_data_melted['condition'].str.replace('ratios_rep_', '')

        # Create the stripplot
        ax = axes[i]
        sns.stripplot(
            x='isoform_number',
            y='ratio',
            hue='condition',
            data=gene_data_melted,
            jitter=jitter,
            alpha=alpha,
            palette=palette,
            ax=ax
        )


        # Add mean values as horizontal lines for each isoform and condition
        i = 0
        for isoform in gene_data['isoform_number'].unique():
            isoform_pos = list(gene_data['isoform_number'].unique()).index(isoform)
            #xax.text(x=isoform_pos-0.1 , y=0.9, s=pvalue_asterisks[i])
            i = i+1
            for j, cond in enumerate(conditions):
                mean_col = f'ratios_{cond}_mean'
                if mean_col in gene_data.columns:
                    mean_val = gene_data[gene_data['isoform_number'] == isoform][mean_col].iloc[0]
                    isoform_pos = list(gene_data['isoform_number'].unique()).index(isoform)
                    ax.hlines(
                        y=mean_val,
                        xmin=isoform_pos-0.2,
                        xmax=isoform_pos+0.2,
                        colors=ax.get_legend().get_lines()[j].get_color(),
                        linewidth=2
                    )

        # Set title and labels
        fdr_text = f", FDR = {fdr:.2e}" if not np.isnan(fdr) else ""
        p_value_text = f", p = {p_value:.2e}"

        # Determine if this gebe has a significant difference
        is_significant = False
        if 'FDR' in gene_data.columns and not np.isnan(fdr):
            is_significant = (fdr <= sig_threshold) & (ratio_difference > difference_threshold)
        else:
            is_significant = p_value <= sig_threshold

        # Set title color based on significance
        title_color = sig_color if is_significant else 'black'

        # Add title with optional stats and color based on significance
        ax.set_title(f"{gene_id}{fdr_text}", color=title_color)
        ax.set_xlabel('isoform_number')
        ax.set_ylabel('Expression Ratio')

        # Set y-limits
        ax.set_ylim(ylim)

        # Adjust legend
        ax.legend(title='Condition', loc='best')

    # Hide unused subplots if any
    for j in range(len(top_genes), rows * cols):
        axes[j].set_visible(False)

    plt.tight_layout()

    # Save figure if requested
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')

    #return fig


def convert_pvalue_to_asterisks(pvalue):
    """Convert p-values to significance asterisks notation."""
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    else:
        return ""


def plot_isoform_diu_results(
    results_df,
    adata=None,
    n=6,
    figsize=(18, 12),
    palette=None,
    jitter=0.2,
    alpha=0.7,
    ylim=(0, 1),
    sort_by='p_value',
    output_file=None,
    sig_threshold=0.05,
    difference_threshold=0.1,
    sig_color='red',
    show_raw_counts=False,
    layer='unique_counts',
    test_condition='control'
):
    """
    Plot the top n syntelogs with differential isoform usage between alleles/haplotypes.

    Parameters
    -----------
    results_df : pd.DataFrame
        Results dataframe from test_isoform1_DIU_between_alleles function
    adata : AnnData, optional
        Original AnnData object containing expression data (needed if show_raw_counts=True)
    n : int, optional
        Number of top syntelogs to plot (default: 6)
    figsize : tuple, optional
        Figure size as (width, height) in inches (default: (18, 12))
    palette : dict or None, optional
        Color palette for haplotypes (default: None, uses seaborn defaults)
    jitter : float, optional
        Amount of jitter for strip plot (default: 0.2)
    alpha : float, optional
        Transparency of points (default: 0.7)
    ylim : tuple, optional
        Y-axis limits for isoform ratios (default: (0, 1))
    sort_by : str, optional
        Column to sort results by ('p_value', 'FDR', or 'ratio_difference') (default: 'p_value')
    output_file : str, optional
        Path to save the figure (default: None, displays figure but doesn't save)
    sig_threshold : float, optional
        Significance threshold for p-value or FDR (default: 0.05)
    difference_threshold : float, optional
        Minimum ratio difference threshold (default: 0.1)
    sig_color : str, optional
        Color for titles of syntelogs with significant differences (default: 'red')
    show_raw_counts : bool, optional
        Whether to show raw counts in addition to ratios (default: False)
    layer : str, optional
        Layer to use for raw counts if show_raw_counts=True (default: 'unique_counts')
    test_condition : str, optional
        Condition to filter samples for if show_raw_counts=True (default: 'control')

    Returns
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import math
    import re

    if len(results_df) == 0:
        print("No results to plot")
        return None

    # Validate sort_by parameter
    if sort_by not in ['p_value', 'FDR', 'ratio_difference']:
        print(f"Invalid sort_by parameter '{sort_by}'. Using 'p_value' instead.")
        sort_by = 'p_value'

    # Ensure FDR column exists
    if 'FDR' not in results_df.columns and sort_by == 'FDR':
        print("FDR column not found. Using p_value for sorting.")
        sort_by = 'p_value'

    # Determine sorting order
    ascending_bool = True if sort_by != 'ratio_difference' else False

    # Get top n syntelogs
    top_syntelogs = results_df.sort_values(sort_by, ascending=ascending_bool).head(n)['Synt_id'].unique()

    # Filter results to include only these syntelogs
    top_results = results_df[results_df['Synt_id'].isin(top_syntelogs)]

    # Calculate grid dimensions - 3 plots per row
    cols = 3
    rows = math.ceil(len(top_syntelogs) / cols)

    # Adjust figure size based on whether we're showing raw counts
    if show_raw_counts and adata is not None:
        figsize = (figsize[0], figsize[1] * 1.5)  # Make taller for subplots

    # Create the figure with grid layout
    fig, axes = plt.subplots(rows, cols, figsize=figsize)

    # Convert axes to flattened array for easier indexing
    if rows == 1 and cols == 1:
        axes = np.array([axes])
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    # Set default palette if none provided
    if palette is None:
        palette = sns.color_palette("Set2", n_colors=10)

    # Plot each syntelog
    for i, synt_id in enumerate(top_syntelogs):
        # Get data for this syntelog
        synt_data = top_results[top_results['Synt_id'] == synt_id].iloc[0]

        # Extract information
        target_isoform = synt_data['target_isoform']
        min_hap = synt_data['min_ratio_haplotype']
        max_hap = synt_data['max_ratio_haplotype']
        min_transcript = synt_data['min_ratio_transcript_id']
        max_transcript = synt_data['max_ratio_transcript_id']
        p_value = synt_data['p_value']
        fdr = synt_data['FDR'] if 'FDR' in synt_data else np.nan
        ratio_difference = synt_data['ratio_difference']

        # Get mean ratios
        min_ratio = synt_data[f'ratio_{min_hap}_mean']
        max_ratio = synt_data[f'ratio_{max_hap}_mean']

        ax = axes[i]

        # If we have access to the original data, plot individual sample points
        if show_raw_counts and adata is not None:
            try:
                # Get sample indices for the condition
                if test_condition == "all":
                    condition_indices = np.arange(adata.shape[0])
                else:
                    condition_indices = np.where(adata.obs['condition'] == test_condition)[0]

                # Find transcript indices
                min_transcript_idx = np.where(adata.var_names == min_transcript)[0]
                max_transcript_idx = np.where(adata.var_names == max_transcript)[0]

                if len(min_transcript_idx) > 0 and len(max_transcript_idx) > 0:
                    min_transcript_idx = min_transcript_idx[0]
                    max_transcript_idx = max_transcript_idx[0]

                    # Get counts for both transcripts
                    counts = adata.layers[layer]

                    # Get all transcripts for this syntelog to calculate total counts per haplotype
                    synt_mask = adata.var['Synt_id'] == synt_id
                    synt_indices = np.where(synt_mask)[0]

                    # Calculate ratios for each sample
                    plot_data = []

                    for hap, transcript_idx in [(min_hap, min_transcript_idx), (max_hap, max_transcript_idx)]:
                        # Get haplotype-specific transcript indices
                        hap_mask = (adata.var['Synt_id'] == synt_id) & (adata.var['haplotype'] == hap)
                        hap_indices = np.where(hap_mask)[0]

                        for sample_idx in condition_indices:
                            # Get counts for this specific transcript
                            isoform_count = counts[sample_idx, transcript_idx]

                            # Get total counts for this haplotype
                            hap_total = np.sum(counts[sample_idx, hap_indices])

                            # Calculate ratio
                            if hap_total > 0:
                                ratio = isoform_count / hap_total
                            else:
                                ratio = 0

                            plot_data.append({
                                'haplotype': hap,
                                'ratio': ratio,
                                'sample': adata.obs_names[sample_idx]
                            })

                    # Convert to DataFrame and plot
                    plot_df = pd.DataFrame(plot_data)

                    if len(plot_df) > 0:
                        sns.stripplot(
                            data=plot_df,
                            x='haplotype',
                            y='ratio',
                            jitter=jitter,
                            alpha=alpha,
                            palette=palette,
                            ax=ax,
                            size=6
                        )

                        # Add mean lines
                        hap_order = [min_hap, max_hap]
                        means = [min_ratio, max_ratio]

                        for j, (hap, mean_val) in enumerate(zip(hap_order, means)):
                            ax.hlines(
                                y=mean_val,
                                xmin=j-0.3,
                                xmax=j+0.3,
                                colors='red',
                                linewidth=3,
                                alpha=0.8
                            )

            except Exception as e:
                print(f"Could not plot individual points for syntelog {synt_id}: {e}")
                # Fall back to bar plot
                show_raw_counts = False

        # If we can't show individual points, show a bar plot of mean ratios
        if not show_raw_counts or adata is None:
            # Get all haplotypes and their ratios for this syntelog
            all_haplotypes = synt_data['all_haplotypes'] if 'all_haplotypes' in synt_data else [min_hap, max_hap]

            # Extract ratios for all haplotypes
            haplotype_ratios_dict = {}
            for hap in all_haplotypes:
                ratio_col = f'ratio_{hap}_mean'
                if ratio_col in synt_data:
                    haplotype_ratios_dict[hap] = synt_data[ratio_col]

            # If we don't have all haplotype data, fall back to min/max
            if not haplotype_ratios_dict:
                haplotype_ratios_dict = {min_hap: min_ratio, max_hap: max_ratio}

            # Sort haplotypes by ratio for better visualization
            sorted_haps = sorted(haplotype_ratios_dict.keys(), key=lambda x: haplotype_ratios_dict[x])
            haplotypes = sorted_haps
            ratios = [haplotype_ratios_dict[hap] for hap in sorted_haps]

            # Create color palette - highlight min and max haplotypes
            colors = []
            for hap in haplotypes:
                if hap == min_hap:
                    colors.append('lightcoral')  # Min ratio in red
                elif hap == max_hap:
                    colors.append('lightgreen')  # Max ratio in green
                else:
                    colors.append('lightgray')   # Other haplotypes in gray

            bars = ax.bar(haplotypes, ratios, color=colors, alpha=0.7, edgecolor='black')

            # Add value labels on bars
            for bar, ratio in zip(bars, ratios):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                       f'{ratio:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=8)

        # Determine if this syntelog is significant
        is_significant = False
        if 'FDR' in synt_data and not pd.isna(fdr):
            is_significant = (fdr <= sig_threshold) and (ratio_difference >= difference_threshold)
        else:
            is_significant = (p_value <= sig_threshold) and (ratio_difference >= difference_threshold)

        # Format title with statistics
        fdr_text = f", FDR = {fdr:.2e}" if not pd.isna(fdr) else ""
        title_color = sig_color if is_significant else 'black'

        ax.set_title(
            f"Syntelog {synt_id}, Isoform {target_isoform}\np = {p_value:.2e}{fdr_text}\nΔratio = {ratio_difference:.3f}",
            color=title_color,
            fontsize=10
        )

        ax.set_xlabel('Haplotype')
        ax.set_ylabel(f'Isoform {target_isoform} Usage Ratio')
        ax.set_ylim(ylim)
        ax.grid(True, alpha=0.3)

        # Rotate x-axis labels if they're long
        if max(len(str(min_hap)), len(str(max_hap))) > 8:
            ax.tick_params(axis='x', rotation=45)

    # Hide unused subplots
    for j in range(len(top_syntelogs), rows * cols):
        if j < len(axes):
            axes[j].set_visible(False)

    # Add overall title
    fig.suptitle(
        f'Top {len(top_syntelogs)} Syntelogs with Differential Isoform Usage Between Alleles\n'
        f'(sorted by {sort_by})',
        fontsize=16,
        y=0.98
    )

    plt.tight_layout()

    # Save figure if requested
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_file}")

    #return fig




def plot_differential_isoform_usage(
    results_df,
    annotation_df,
    fdr_threshold=0.05,
    ratio_difference_threshold=0.1,
    width=1300,
    height=350,
    template="simple_white"
):
    """
    Visualize differential isoform usage results using transcript structure plots.

    Filtering is performed at the gene level: if any isoform of a gene meets the
    significance criteria, all isoforms of that gene will be plotted.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame containing differential isoform usage results with columns:
        - gene_id: Gene identifier
        - transcript_id: Transcript identifier
        - sample_name: Sample name
        - {layer}_cpm: CPM normalized counts for the specified layer
        - isoform_ratio: Isoform usage ratio
        - condition: Experimental condition
        - FDR: False discovery rate (optional for filtering)
        - ratio_difference: Difference in ratios (optional for filtering)
    annotation_df : polars.DataFrame
        Polars DataFrame containing GTF annotation with:
        - gene_id: Gene identifier
        - functional_annotation: Gene annotation text (optional)
        Plus standard GTF columns for RNApysoforms with columns:
        - gene_id: Gene identifier
        - functional_annotation: Gene annotation text (optional)
        Plus standard GTF columns for RNApysoforms
    fdr_threshold : float, default=0.05
        FDR threshold for identifying significant genes (applied at gene level)
    ratio_difference_threshold : float, default=0.1
        Minimum ratio difference threshold for identifying significant genes (applied at gene level)
    width : int, default=1300
        Plot width in pixels
    height : int, default=350
        Plot height in pixels
    template : str, default="simple_white"
        Plotly template for styling

    Returns
    -------
    list
        List of generated plotly figures

    Raises
    ------
    ImportError
        If required packages (RNApysoforms) are not available
    ValueError
        If required columns are missing from input DataFrames
    """
    try:
        import RNApysoforms as RNApy
        import pandas as pd
        import polars as pl
    except ImportError as e:
        raise ImportError(f"Required package not available: {e}. "
                         "Please install RNApysoforms for transcript visualization.")

    # Detect the layer being used by looking for columns ending with '_cpm'
    layer = None
    layer_cpm_col = None

    for col in results_df.columns:
        if col.endswith('_cpm'):
            layer_cpm_col = col
            layer = col.replace('_cpm', '')
            break

    if layer is None:
        # Fallback to looking for 'unique_counts_cpm' for backward compatibility
        if 'unique_counts_cpm' in results_df.columns:
            layer = 'unique_counts'
            layer_cpm_col = 'unique_counts_cpm'
        else:
            raise ValueError("No CPM column found. Expected column ending with '_cpm'")

    print(f"Detected layer: {layer} (using column: {layer_cpm_col})")

    # Validate required columns
    required_cols = ['gene_id', 'transcript_id', 'sample_name', layer_cpm_col,
                     'isoform_ratio', 'condition']
    missing_cols = [col for col in required_cols if col not in results_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in results_df: {missing_cols}")

    # Filter results if thresholds provided - identify significant genes
    significant_genes = set()

    # Find genes with at least one significant isoform
    if 'FDR' in results_df.columns:
        significant_genes.update(
            results_df[results_df['FDR'] <= fdr_threshold]['gene_id'].unique()
        )

    if 'ratio_difference' in results_df.columns:
        significant_genes.update(
            results_df[results_df['ratio_difference'] >= ratio_difference_threshold]['gene_id'].unique()
        )

    # If no filtering columns are present, include all genes
    if 'FDR' not in results_df.columns and 'ratio_difference' not in results_df.columns:
        significant_genes = set(results_df['gene_id'].unique())

    # Filter to only include data for significant genes (but keep all isoforms of those genes)
    filtered_df = results_df[results_df['gene_id'].isin(significant_genes)].copy()

    print(f"Found {len(significant_genes)} genes with significant isoforms. Plotting all isoforms for these genes.")

    figures = []

    # Create dynamic subplot titles based on detected layer
    layer_display_name = layer.replace('_', ' ').title()
    subplot_titles = [
        "Transcript Structure",
        f"CPM ({layer_display_name})",
        "Isoform Usage Ratio"
    ]

    # Process each gene
    for gene in filtered_df['gene_id'].unique():
        print(f"Processing gene: {gene}")

        # Get annotation text (optional)
        annotation_txt = ""
        try:
            # Check if functional_annotation column exists
            if "functional_annotation" in annotation_df.columns:
                filtered_annotation = annotation_df.filter(pl.col("gene_id") == gene)
                if len(filtered_annotation) > 0:
                    annotation_txt = filtered_annotation.head(1)["functional_annotation"][0].split('/')
                    annotation_txt = list(dict.fromkeys(annotation_txt))
                    annotation_txt = ' '.join(annotation_txt)
                else:
                    print(f"Warning: No annotation found for gene {gene}")
            else:
                print(f"Info: No functional_annotation column found, using gene ID only")
        except (IndexError, KeyError, TypeError) as e:
            print(f"Warning: Could not extract functional annotation for {gene}: {e}")
            annotation_txt = ""

        # Subset counts matrix for this gene
        counts_matrix_gene = filtered_df[filtered_df['gene_id'] == gene].copy()

        # Skip if only one transcript
        if len(counts_matrix_gene['transcript_id'].unique()) < 2:
            print(f"Skipping {gene} as it has only one transcript")
            continue

        # Prepare data for RNApysoforms
        counts_matrix_gene['sample_id'] = counts_matrix_gene['sample_name']
        counts_matrix_gene['counts'] = counts_matrix_gene[layer_cpm_col].astype(float)
        counts_matrix_gene = pl.from_pandas(counts_matrix_gene)

        try:
            # Filter annotation and counts matrix
            sod1_annotation, sod1_counts_matrix = RNApy.gene_filtering(
                annotation=annotation_df,
                expression_matrix=counts_matrix_gene,
                target_gene=gene
            )

            # Create traces
            traces = RNApy.make_traces(
                annotation=sod1_annotation,
                expression_matrix=sod1_counts_matrix,
                y='transcript_id',
                annotation_hue="type",
                hover_start="start",
                hover_end="end",
                expression_columns=["counts", "isoform_ratio"],
                expression_hue="condition"
            )

            # Create figure with dynamic subplot titles
            fig = RNApy.make_plot(
                traces=traces,
                subplot_titles=subplot_titles,
                width=width,
                height=height,
                boxgap=0.1,
                boxgroupgap=0.5,
                horizontal_spacing=0.08,
                template=template,
                column_widths=[0.4, 0.3, 0.3]
            )

            # Update layout with title
            title_text = f"{gene}"
            if annotation_txt:
                title_text += f" {annotation_txt}"

            fig.update_layout(title=dict(
                text=title_text,
                x=0.5,
                xanchor="center"
            ))

            # Show the plot
            fig.show()
            figures.append(fig)

        except Exception as e:
            print(f"Error processing gene {gene}: {str(e)}")
            continue

    print(f"Generated {len(figures)} plots for differential isoform usage")
    return figures

def plot_allele_specific_isoform_structure(
    results_df,
    annotation_df,
    ratio_difference_threshold=0.2,
    width=1300,
    height=350,
    template="simple_white"
):
    """
    Visualize allele-specific isoform structure differences.

    Filtering is performed at the syntelog level: if any isoform of a syntelog meets the
    significance criteria, all isoforms of that syntelog will be plotted.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame containing allele-specific results with columns:
        - Synt_id: Syntelog identifier
        - isoform_id: Isoform/transcript identifier
        - gene_id: Gene identifier
        - haplotype: Haplotype identifier
        - reference_haplotype: Reference haplotype identifier
        - sample: Sample identifier
        - isoform_counts: Isoform counts (will be detected dynamically)
        - isoform_ratio: Isoform usage ratio
        - ratio_difference: Difference in ratios between haplotypes
    annotation_df : polars.DataFrame
        Polars DataFrame containing GTF annotation
    ratio_difference_threshold : float, default=0.2
        Minimum ratio difference threshold for identifying significant syntelogs (applied at syntelog level)
    width : int, default=1300
        Plot width in pixels
    height : int, default=350
        Plot height in pixels
    template : str, default="simple_white"
        Plotly template for styling

    Returns
    -------
    list
        List of generated plotly figures

    Raises
    ------
    ImportError
        If required packages are not available
    ValueError
        If required columns are missing
    """
    try:
        import RNApysoforms as RNApy
        import pandas as pd
        import polars as pl
    except ImportError as e:
        raise ImportError(f"Required package not available: {e}")

    # Detect the layer being used by looking for columns with count data
    layer = None
    layer_counts_col = None

    # Look for various count column patterns
    count_patterns = ['_cpm', '_counts', '_raw']
    for pattern in count_patterns:
        for col in results_df.columns:
            if col.endswith(pattern) and 'isoform' in col:
                layer_counts_col = col
                if pattern == '_cpm':
                    layer = col.replace('_cpm', '')
                elif pattern == '_counts':
                    layer = col.replace('_counts', '')
                elif pattern == '_raw':
                    layer = col.replace('_raw', '')
                break
        if layer is not None:
            break

    # Fallback to looking for 'isoform_counts' for backward compatibility
    if layer is None:
        if 'isoform_counts' in results_df.columns:
            layer = 'unique'  # Default assumption
            layer_counts_col = 'isoform_counts'
        else:
            raise ValueError("No count column found. Expected column with pattern '_cpm', '_counts', or '_raw', or 'isoform_counts'")

    print(f"Detected layer: {layer} (using column: {layer_counts_col})")

    # Validate required columns
    required_cols = ['Synt_id', 'isoform_id', 'gene_id', 'haplotype',
                     'reference_haplotype', 'sample', layer_counts_col, 'isoform_ratio']
    missing_cols = [col for col in required_cols if col not in results_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in results_df: {missing_cols}")

    # Filter by ratio difference threshold - identify significant syntelogs
    significant_syntelogs = set()

    # Find syntelogs with at least one significant isoform
    significant_syntelogs.update(
        results_df[results_df['ratio_difference'] >= ratio_difference_threshold]['Synt_id'].unique()
    )

    # Filter to only include data for significant syntelogs (but keep all isoforms of those syntelogs)
    filtered_df = results_df[results_df['Synt_id'].isin(significant_syntelogs)].copy()

    print(f"Found {len(significant_syntelogs)} syntelogs with significant allelic differences. Plotting all isoforms for these syntelogs.")

    figures = []

    # Create dynamic subplot titles based on detected layer
    layer_display_name = layer.replace('_', ' ').title()

    # Determine the type of data based on column pattern
    if layer_counts_col.endswith('_cpm'):
        data_type = f"CPM ({layer_display_name})"
    elif layer_counts_col.endswith('_raw'):
        data_type = f"Raw Counts ({layer_display_name})"
    elif layer_counts_col.endswith('_counts'):
        data_type = f"Counts ({layer_display_name})"
    else:
        data_type = f"Counts ({layer_display_name})"  # Default

    subplot_titles = [
        "Transcript Structure",
        data_type,
        "Isoform Usage Ratio"
    ]

    # Process each syntelog
    for synt_id in filtered_df['Synt_id'].unique():
        print(f"Processing Synt_id: {synt_id}")

        # Get gene ID from reference haplotype
        try:
            gene = (filtered_df[
                (filtered_df['Synt_id'] == synt_id) &
                (filtered_df['haplotype'] == filtered_df['reference_haplotype'])
            ]['gene_id'].unique()[0])
            print(f"Gene: {gene}")
        except IndexError:
            print(f"Warning: No reference haplotype found for Synt_id {synt_id}")
            continue

        # Get annotation text (optional)
        annotation_txt = ""
        try:
            # Check if functional_annotation column exists
            if "functional_annotation" in annotation_df.columns:
                filtered_annotation = annotation_df.filter(pl.col("gene_id") == gene)
                if len(filtered_annotation) > 0:
                    annotation_txt = filtered_annotation.head(1)["functional_annotation"][0].split('/')
                    annotation_txt = list(dict.fromkeys(annotation_txt))
                    annotation_txt = ' '.join(annotation_txt)
                else:
                    print(f"Warning: No annotation found for gene {gene}")
            else:
                print(f"Info: No functional_annotation column found, using gene ID only")
        except (IndexError, KeyError, TypeError) as e:
            print(f"Warning: Could not extract functional annotation for {gene}: {e}")
            annotation_txt = ""

        # Subset for this syntelog
        counts_matrix_gene = filtered_df[filtered_df['Synt_id'] == synt_id].copy()
        counts_matrix_gene['gene_id'] = gene
        counts_matrix_gene['transcript_id'] = counts_matrix_gene['isoform_id']

        # Skip if only one transcript
        if len(counts_matrix_gene['transcript_id'].unique()) < 2:
            print(f"Skipping {gene} (Synt_id: {synt_id}) as it has only one transcript")
            continue

        # Prepare data for RNApysoforms
        counts_matrix_gene['sample_id'] = counts_matrix_gene['sample']
        counts_matrix_gene['counts'] = counts_matrix_gene[layer_counts_col].astype(float)
        counts_matrix_gene = pl.from_pandas(counts_matrix_gene)
        # for rows where isoform_ratio is 0.0, set to None
        counts_matrix_gene = counts_matrix_gene.with_columns([
            pl.when(pl.col("isoform_ratio") == 0.0)
            .then(None)
            .otherwise(pl.col("isoform_ratio"))
            .alias("isoform_ratio")
        ])
        try:
            # Filter annotation and counts matrix
            sod1_annotation, sod1_counts_matrix = RNApy.gene_filtering(
                annotation=annotation_df,
                expression_matrix=counts_matrix_gene,
                target_gene=gene
            )

            #sod1_annotation = RNApy.shorten_gaps(sod1_annotation)
            # This is not working

            # Create traces
            traces = RNApy.make_traces(
                annotation=sod1_annotation,
                expression_matrix=sod1_counts_matrix,
                y='transcript_id',
                annotation_hue="type",
                hover_start="start",
                hover_end="end",
                expression_columns=["counts", "isoform_ratio"],
                expression_hue="haplotype",
                marker_size=2
            )

            

            # Create figure with dynamic subplot titles
            fig = RNApy.make_plot(
                traces=traces,
                subplot_titles=subplot_titles,
                width=width,
                height=height,
                boxgap=0.1,
                boxgroupgap=0.5,
                horizontal_spacing=0.08,
                template=template,
                column_widths=[0.4, 0.3, 0.3]
            )

            # Update layout with title
            title_text = f"{gene}"
            if annotation_txt:
                title_text += f" {annotation_txt}"
            title_text += f" (Synt_id: {synt_id})"

            fig.update_layout(title=dict(
                text=title_text,
                x=0.5,
                xanchor="center"
            ))

            # Show the plot
            fig.show()
            figures.append(fig)

        except Exception as e:
            print(f"Error processing Synt_id {synt_id} (gene {gene}): {str(e)}")
            continue

    print(f"Generated {len(figures)} plots for allele-specific isoform structure")
    return figures
