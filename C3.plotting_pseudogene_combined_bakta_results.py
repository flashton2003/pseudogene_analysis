import pandas as pd
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns

def determine_analysis_type(filename):
    if 'dbs' in filename.lower():
        return 'delta_bit_score'
    elif 'bakta_db' in filename.lower():
        return 'bakta_db'
    elif 'ncbi' in filename.lower():
        return 'ncbi_db'
    elif 'salmonella' in filename.lower():
        return 'salmonella_db'
    print(filename)
    return 'unknown'

def extract_gcf(filename):
    # Regular expression to match GCF patterns
    match = re.search(r'(GCF_\d{9}\.\d)', filename)
    if match:
        return match.group(1)
    return None

def process_files(directory, ei_gi_lookup):
    # List to store all processed data
    all_data = []
    
    # Get all TSV files in the directory
    files = glob.glob(f"{directory}/*.tsv")
    
    for file in files:
        try:
            # Read TSV file
            df = pd.read_csv(file, sep='\t')
            
            # Process each row in the file
            for _, row in df.iterrows():
                analysis_type = determine_analysis_type(row['input_file'])
                gcf = extract_gcf(row['input_file'])
                
                # Get salmonella type from lookup
                salmonella_type = ei_gi_lookup.get(gcf, 'unknown')
                if salmonella_type == 'unknown':
                    print(f"Unknown salmonella type for GCF {gcf}")
                
                all_data.append({
                    'sensitivity': float(row['sensitivity'].strip('%')),
                    'ppv': float(row['ppv'].strip('%')),
                    'analysis_type': analysis_type,
                    'salmonella_type': salmonella_type,
                    'strain': row['strain']
                })
                
        except Exception as e:
            print(f"Error processing file {file}: {str(e)}")
    
    return pd.DataFrame(all_data)

def create_visualization(df):
    # Set the style
    # plt.style.use('seaborn')
    
    # Create color mapping for analysis types
    color_map = {
        'delta_bit_score': 'orange',
        'bakta_db': 'red',
        'ncbi_db': 'green',
        'salmonella_db': 'purple'
    }
    
    # Create marker mapping for salmonella types
    marker_map = {
        'EI': 'o',  # circle
        'GI': '^'   # triangle
    }
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Plot points for each analysis type and salmonella type
    for analysis in df['analysis_type'].unique():
        for salmonella in df['salmonella_type'].unique():
            mask = (df['analysis_type'] == analysis) & (df['salmonella_type'] == salmonella)
            data = df[mask]
            
            if not data.empty:
                plt.scatter(
                    data['sensitivity'],
                    data['ppv'],
                    c=color_map.get(analysis, 'gray'),
                    marker=marker_map.get(salmonella, 's'),
                    label=f'{analysis} - {salmonella}',
                    s=100,
                    alpha=0.7
                )
    
    # Customize the plot
    plt.xlabel('Sensitivity (%)')
    plt.ylabel('PPV (%)')
    plt.title('Sensitivity vs PPV by Analysis Type and Salmonella Type')
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    
    # Add a tight layout to prevent label cutoff
    plt.tight_layout()
    
    return plt

# Main execution
def main():
    # Your provided dictionary
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
    
    # Process the files
    df = process_files('2024.11.11', ei_gi_lookup)
    
    # Create and save the visualization
    plt = create_visualization(df)
    plt.savefig('sensitivity_ppv_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()