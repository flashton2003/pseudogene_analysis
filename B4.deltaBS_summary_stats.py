import pandas as pd
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Define strain mapping
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

def get_reference_name(accession):
    """Convert accession to reference name using the mapping."""
    if accession not in STRAIN_MAPPING:
        raise ValueError(f"Unknown accession: {accession}. Valid accessions are: {', '.join(STRAIN_MAPPING.keys())}")
    return STRAIN_MAPPING[accession]

def analyze_genes(data_file, accession, dbs_threshold):
    # Convert accession to reference name
    reference_name = get_reference_name(accession)
    
    # Read the Excel file
    df = pd.read_excel(data_file, sheet_name='Deduplicated_Data')
    
    # Filter out genes that are absent from the reference genome
    reference_col = df[reference_name].astype(str)
    mask_present = ~reference_col.str.startswith('3')
    present_genes = df[mask_present].copy()
    reference_col = reference_col[mask_present]
    
    # Calculate basic counts that don't depend on threshold
    total_genes = len(present_genes)
    genes_with_blast = present_genes['qseqid_fwd'].notna().sum()
    genes_with_dbs = present_genes['gene_1'].notna().sum()
    
    # Identify true HDCs (those starting with '2' in reference column)
    true_hdc_mask = reference_col.str.startswith('2')
    true_hdcs = present_genes[true_hdc_mask]
    total_true_hdcs = len(true_hdcs)
    
    # Calculate True HDCs with blast hits and dbs matches
    true_hdcs_with_blast = true_hdcs['qseqid_fwd'].notna().sum()
    true_hdcs_with_dbs = true_hdcs['gene_1'].notna().sum()
    
    # Identify predicted HDCs using the provided threshold
    predicted_hdcs_mask = (present_genes['delta-bitscore'] > dbs_threshold) & \
                         (present_genes['loss_of_function'] == 1)
    total_predicted_hdcs = predicted_hdcs_mask.sum()
    
    # Calculate true positives (correctly identified HDCs)
    true_positives = (predicted_hdcs_mask & reference_col.str.startswith('2')).sum()
    
    # Calculate sensitivity and PPV
    sensitivity = true_positives / total_true_hdcs if total_true_hdcs > 0 else 0
    ppv = true_positives / total_predicted_hdcs if total_predicted_hdcs > 0 else 0
    
    # Split the summary into two dictionaries
    summary_counts = {
        'Accession': accession,
        'Reference_Name': reference_name,
        'Total_Genes': total_genes,
        'Genes_with_Blast': genes_with_blast,
        'Genes_with_DBS': genes_with_dbs,
        'True_HDCs': total_true_hdcs,
        'True_HDCs_with_Blast': true_hdcs_with_blast,
        'True_HDCs_with_DBS': true_hdcs_with_dbs
    }
    
    summary_performance = {
        'DBS_Threshold': dbs_threshold,
        'True_HDCs': total_true_hdcs,
        'Predicted_HDCs': total_predicted_hdcs,
        'True_Positives': true_positives,
        'Sensitivity': sensitivity,
        'PPV': ppv
    }
    
    return summary_counts, summary_performance

def plot_sensitivity_ppv(summaries, output_dir):
    # Set the Seaborn style and context
    sns.set_style("whitegrid")
    sns.set_context("talk")
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'Sensitivity': [s['Sensitivity'] for s in summaries],
        'PPV': [s['PPV'] for s in summaries],
        'Threshold': [s['DBS_Threshold'] for s in summaries]
    })
    
    # Create figure with specified size
    plt.figure(figsize=(10, 10))
    
    # Create scatter plot with connected lines using Seaborn
    sns.scatterplot(data=plot_data, x='Sensitivity', y='PPV', s=100)
    plt.plot(plot_data['Sensitivity'], plot_data['PPV'], 'b-', alpha=0.6)
    
    # Add threshold labels to each point
    for idx, row in plot_data.iterrows():
        plt.annotate(f'DBS>{row["Threshold"]:.1f}', 
                    (row['Sensitivity'], row['PPV']),
                    xytext=(10, 5),
                    textcoords='offset points',
                    fontsize=10)
    
    # Add labels and title with Seaborn's style
    plt.xlabel('Sensitivity', fontsize=12)
    plt.ylabel('Positive Predictive Value', fontsize=12)
    plt.title('Sensitivity vs PPV for Different Delta-bitscore Thresholds', 
             fontsize=14, pad=20)
    
    # Set axis limits with some padding
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    
    # Add diagonal line for reference
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot with high DPI
    plot_path = os.path.join(output_dir, 'sensitivity_vs_ppv.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <excel_file> <accession>")
        print("\nValid accessions:")
        for acc, ref in STRAIN_MAPPING.items():
            print(f"  {acc} -> {ref}")
        sys.exit(1)
        
    excel_file = sys.argv[1]
    accession = sys.argv[2]
    
    try:
        # Get reference name from accession
        reference_name = get_reference_name(accession)
        
        # Create output directory based on accession and reference name
        output_dir = f"deltabs_analysis_{accession}_{reference_name}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Read data to calculate percentiles
        df = pd.read_excel(excel_file, sheet_name='Deduplicated_Data')
        percentile_97_5 = np.percentile(df['delta-bitscore'].dropna(), 97.5)
        percentile_99 = np.percentile(df['delta-bitscore'].dropna(), 99)
        
        # Define thresholds including both percentiles
        thresholds = [0, 2, 4, 6, 8, 10, percentile_97_5, percentile_99]
        
        # Calculate summaries for each threshold
        summaries_counts = []
        summaries_performance = []
        for threshold in thresholds:
            counts, performance = analyze_genes(excel_file, accession, threshold)
            summaries_counts.append(counts)
            summaries_performance.append(performance)
        
        # Create plot
        plot_sensitivity_ppv(summaries_performance, output_dir)
        
        # Prepare output files
        counts_file = os.path.join(output_dir, f"deltabs_gene_counts_{accession}.tsv")
        performance_file = os.path.join(output_dir, f"deltabs_performance_{accession}.tsv")
        
        # Write counts table
        with open(counts_file, 'w') as f:
            count_headers = ['Accession', 'Reference_Name', 'Total_Genes', 'Genes_with_Blast', 
                           'Genes_with_DBS', 'True_HDCs', 'True_HDCs_with_Blast', 
                           'True_HDCs_with_DBS']
            print("\nGene Counts Summary:")
            print('\t'.join(count_headers))
            f.write('\t'.join(count_headers) + '\n')
            
            # Only write the first row since counts don't change with threshold
            row = [str(summaries_counts[0][h]) for h in count_headers]
            print('\t'.join(row))
            f.write('\t'.join(row) + '\n')
        
        # Write performance table
        with open(performance_file, 'w') as f:
            perf_headers = ['DBS_Threshold', 'True_HDCs', 'Predicted_HDCs', 
                          'True_Positives', 'Sensitivity', 'PPV']
            print("\nPerformance at Different Thresholds:")
            print('\t'.join(perf_headers))
            f.write('\t'.join(perf_headers) + '\n')
            
            for summary in summaries_performance:
                row = [
                    f"{summary['DBS_Threshold']:.3f}",
                    str(summary['True_HDCs']),
                    str(summary['Predicted_HDCs']),
                    str(summary['True_Positives']),
                    f"{summary['Sensitivity']:.3f}",
                    f"{summary['PPV']:.3f}"
                ]
                print('\t'.join(row))
                f.write('\t'.join(row) + '\n')
                
        print(f"\nResults saved to directory: {output_dir}")
        print(f"Counts file: {counts_file}")
        print(f"Performance file: {performance_file}")
        print(f"Plot file: {os.path.join(output_dir, 'sensitivity_vs_ppv.png')}")
            
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
