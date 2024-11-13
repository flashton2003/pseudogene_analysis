import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
import numpy as np

# [Previous mappings and helper functions remain unchanged...]
STRAIN_MAPPING = {
    'GCF_000020705.1': 'SL476',
    'GCF_000020745.1': 'CVM19633',
    'GCF_000020885.1': 'SL483',
    'GCF_000009505.1': 'P125109',
    'GCF_000018705.1': 'SPB7',
    'GCF_000195995.1': 'CT18',
    'GCF_000007545.1': 'Ty2',
    'GCF_000011885.1': 'ATCC 9150',
    'GCF_000020925.1': 'CT_02021853',
    'GCF_000009525.1': '287/91',
    'GCF_000008105.1': 'SC-B67',
    'GCF_000018385.1': 'RKS4594',
    'GCF_000026565.1': 'AKU_12601'
}

ei_gi_lookup = {
    'GCF_000007545.1': 'EI',
    'GCF_000008105.1': 'EI',
    'GCF_000009505.1': 'GI',
    'GCF_000009525.1': 'EI',
    'GCF_000011885.1': 'EI',
    'GCF_000018385.1': 'EI',
    'GCF_000018705.1': 'GI',
    'GCF_000020705.1': 'GI',
    'GCF_000020745.1': 'GI',
    'GCF_000020885.1': 'GI',
    'GCF_000020925.1': 'EI',
    'GCF_000195995.1': 'EI',
    'GCF_000026565.1': 'EI'
}

color_map = {
    'dbs': 'orange',
    'bakta_db': 'red',
    'ncbi': 'green',
    'salmonella': 'purple'
}

marker_map = {
    'EI': 'o',  # circle
    'GI': '^'   # triangle
}

def extract_accession(filename):
    """Extract GCF accession from filename."""
    match = re.search(r'(GCF_\d+\.\d+)', filename)
    return match.group(1) if match else None

def extract_tool(filename):
    """Extract tool name from filename."""
    if 'dbs' in filename.lower():
        return 'dbs'
    elif 'bakta_db' in filename.lower():
        return 'bakta_db'
    elif 'ncbi' in filename.lower():
        return 'ncbi'
    elif 'salmonella' in filename.lower():
        return 'salmonella'
    return 'unknown'

def analyze_files(file_pattern):
    results = []
    
    for file_path in glob.glob(file_pattern):
        try:
            accession = extract_accession(file_path)
            tool = extract_tool(file_path)
            
            if not accession:
                print(f"Could not extract accession from {file_path}")
                continue
                
            df = pd.read_excel(file_path)
            strain_col = STRAIN_MAPPING.get(accession)
            
            if not strain_col:
                print(f"No strain mapping found for accession {accession}")
                continue
            
            true_pseudogenes = df[df[strain_col].astype(str).str.startswith('2')].shape[0]
            pseudogene_calls = df['is_pseudogene'].sum()
            pathogen_type = ei_gi_lookup.get(accession)
            
            results.append({
                'file': file_path,
                'accession': accession,
                'tool': tool,
                'true_pseudogenes': true_pseudogenes,
                'pseudogene_calls': pseudogene_calls,
                'pathogen_type': pathogen_type
            })
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
    
    return pd.DataFrame(results)

def plot_correlation_by_tool(results_df, tool):
    """Create correlation plot for a specific tool."""
    tool_data = results_df[results_df['tool'] == tool]
    
    if tool_data.empty:
        return None
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot points for each pathogen type
    for pathogen_type in marker_map.keys():
        mask = tool_data['pathogen_type'] == pathogen_type
        data = tool_data[mask]
        
        if not data.empty:
            ax.scatter(
                data['true_pseudogenes'],
                data['pseudogene_calls'],
                c=color_map[tool],
                marker=marker_map[pathogen_type],
                label=f'{pathogen_type}',
                alpha=0.7,
                s=100
            )
    
    # Calculate and add correlation line
    slope, intercept = np.polyfit(tool_data['true_pseudogenes'], tool_data['pseudogene_calls'], 1)
    x_values = tool_data['true_pseudogenes']
    y_values = slope * x_values + intercept
    ax.plot(
        x_values,
        y_values,
        "k--",
        alpha=0.5,
        label='Trend line'
    )
    
    # Calculate correlation coefficient
    corr = tool_data['true_pseudogenes'].corr(tool_data['pseudogene_calls'])
    
    # Add correlation coefficient and equation to plot
    equation = f'y = {slope:.2f}x + {intercept:.2f}'
    stats_text = f'Correlation: {corr:.2f}\n{equation}'
    ax.text(
        0.05, 0.95,
        stats_text,
        transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.8),
        verticalalignment='top'
    )
    
    # Set axes to start at 0
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    # Calculate the maximum values for x and y axes with some padding
    max_x = tool_data['true_pseudogenes'].max() * 1.1
    max_y = tool_data['pseudogene_calls'].max() * 1.1
    
    # Set the upper limits
    ax.set_xlim(right=max_x)
    ax.set_ylim(top=max_y)
    
    ax.set_xlabel('True Pseudogenes')
    ax.set_ylabel('Pseudogene Calls')
    ax.set_title(f'Correlation for {tool.upper()}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

# Main execution
if __name__ == "__main__":
    # Analyze files
    results_df = analyze_files("2024.11.11/*xlsx")
    
    if not results_df.empty:
        # Create separate plots for each tool
        for tool in color_map.keys():
            plot = plot_correlation_by_tool(results_df, tool)
            if plot:
                plot.savefig(f'pseudogene_correlation_{tool}.png', bbox_inches='tight', dpi=300)
                plt.close(plot)
        
        # Print summary statistics
        print("\nSummary Statistics:")
        print(f"Total files processed: {len(results_df)}")
        print("\nBy Tool:")
        print(results_df.groupby('tool').size())
        print("\nBy Pathogen Type:")
        print(results_df.groupby('pathogen_type').size())
        
        # Print correlation coefficients and equations for each tool
        print("\nStatistics by tool:")
        for tool in color_map.keys():
            tool_data = results_df[results_df['tool'] == tool]
            if not tool_data.empty:
                slope, intercept = np.polyfit(tool_data['true_pseudogenes'], tool_data['pseudogene_calls'], 1)
                corr = tool_data['true_pseudogenes'].corr(tool_data['pseudogene_calls'])
                print(f"{tool}:")
                print(f"  Correlation: {corr:.3f}")
                print(f"  Equation: y = {slope:.2f}x + {intercept:.2f}")
    else:
        print("No files were processed successfully.")