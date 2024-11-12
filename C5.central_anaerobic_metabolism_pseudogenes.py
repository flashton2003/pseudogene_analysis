import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

# Existing mapping dictionaries remain the same
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
    'GCF_000007545.1': 'EI', 'GCF_000008105.1': 'EI', 
    'GCF_000009505.1': 'GI', 'GCF_000009525.1': 'EI', 
    'GCF_000011885.1': 'EI', 'GCF_000018385.1': 'EI', 
    'GCF_000018705.1': 'GI', 'GCF_000020705.1': 'GI', 
    'GCF_000020745.1': 'GI', 'GCF_000020885.1': 'GI', 
    'GCF_000020925.1': 'EI', 'GCF_000195995.1': 'EI', 
    'GCF_000026565.1': 'EI'
}

color_map = {
    'dbs': 'orange',
    'bakta_db': 'red',
    'ncbi': 'green',
    'salmonella': 'purple'
}

marker_map = {
    'EI': 'o',
    'GI': '^'
}

def get_tool_type(filename):
    if 'dbs' in filename.lower():
        return 'dbs'
    elif 'bakta_db' in filename.lower():
        return 'bakta_db'
    elif 'ncbi' in filename.lower():
        return 'ncbi'
    elif 'salmonella' in filename.lower():
        return 'salmonella'
    return None

def get_accession(filename):
    match = re.search(r'(GCF_\d{9}\.\d)', filename)
    return match.group(1) if match else None

def analyze_pseudogenes():
    # Read central anaerobic genes reference
    central_genes = pd.read_excel("2024.11.05b/mbo001141769st7.central_anaerobic_genes.xlsx")
    
    results = []
    
    for filepath in glob.glob("2024.11.11/*xlsx"):
        try:
            accession = get_accession(filepath)
            if not accession:
                continue
                
            tool_type = get_tool_type(filepath)
            pathogen_type = ei_gi_lookup.get(accession)
            strain = STRAIN_MAPPING.get(accession)
            
            df = pd.read_excel(filepath)
            merged = pd.merge(
                df,
                central_genes,
                on="Reference locus tag(s)",
                how="inner"
            )
            
            true_pseudogenes = merged[merged[strain].str.startswith('2', na=False)].shape[0]
            called_pseudogenes = merged[merged['is_pseudogene'] == 1].shape[0]
            
            results.append({
                'accession': accession,
                'tool': tool_type,
                'pathogen': pathogen_type,
                'true_pseudogenes': true_pseudogenes,
                'called_pseudogenes': called_pseudogenes
            })
            
        except Exception as e:
            print(f"Error processing {filepath}: {str(e)}")
    
    return pd.DataFrame(results)

def plot_tool_correlation(results_df, tool_name, output_dir='plots'):
    plt.figure(figsize=(10, 8))
    
    # Filter data for specific tool
    tool_data = results_df[results_df['tool'] == tool_name]
    
    # Plot points for each pathogen type
    for pathogen in marker_map.keys():
        mask = tool_data['pathogen'] == pathogen
        data = tool_data[mask]
        
        if not data.empty:
            plt.scatter(
                data['true_pseudogenes'],
                data['called_pseudogenes'],
                c=color_map[tool_name],
                marker=marker_map[pathogen],
                label=f"{pathogen}",
                s=100,
                alpha=0.7
            )
    
    # Add diagonal line
    max_val = max(
        tool_data['true_pseudogenes'].max(),
        tool_data['called_pseudogenes'].max()
    )
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    # Calculate correlation coefficient for this tool
    correlation = tool_data['true_pseudogenes'].corr(tool_data['called_pseudogenes'])
    
    # Customize plot
    plt.xlabel('True Pseudogenes (2* in strain column)')
    plt.ylabel('Called Pseudogenes (is_pseudogene=1)')
    plt.title(f'Correlation for {tool_name.upper()}\nPseudogenes in Central Anaerobic Metabolism\nr = {correlation:.3f}')
    plt.legend(title='Pathogen Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Ensure equal aspect ratio
    plt.axis('equal')
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save plot
    output_path = os.path.join(output_dir, f'pseudogenes_correlation_{tool_name}.png')
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    return correlation

def main():
    # Create results directory
    output_dir = 'correlation_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # Run analysis
    results_df = analyze_pseudogenes()
    
    # Create individual plots for each tool
    correlations = {}
    for tool in color_map.keys():
        correlation = plot_tool_correlation(results_df, tool, output_dir)
        correlations[tool] = correlation
        
    # Print correlation coefficients
    print("\nCorrelation coefficients by tool:")
    for tool, corr in correlations.items():
        print(f"{tool.upper()}: {corr:.3f}")

    results_df.to_excel('2024.11.11/central_anaerobic_metabolism_pseudogene_analysis_results.xlsx', index=False)
    
    return results_df

if __name__ == "__main__":
    main()