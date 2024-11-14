import pandas as pd
import os
from pathlib import Path

def calculate_metrics(df, method_col, truth_col):
    """Calculate PPV, sensitivity, and counts for a given method column"""
    # Convert truth values starting with '2' to 1, others to 0
    true_positives = ((df[truth_col].astype(str).str.startswith('2')) & (df[method_col] == 1)).sum()
    false_positives = ((~df[truth_col].astype(str).str.startswith('2')) & (df[method_col] == 1)).sum()
    false_negatives = ((df[truth_col].astype(str).str.startswith('2')) & (df[method_col] == 0)).sum()
    
    # Calculate total positives
    total_positives = (df[method_col] == 1).sum()
    
    # Calculate PPV and sensitivity
    ppv = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    
    return ppv, sensitivity, total_positives, true_positives

def count_central_metabolism_pseudogenes(df, method_col):
    """Count pseudogenes that are also central metabolism genes"""
    return ((df[method_col] == 1) & (df['central_anaerobic_metabolism'] == 1)).sum()

def analyze_file(excel_path, strain_mapping):
    """Analyze a single Excel file and return its metrics"""
    # Read the Excel file
    df = pd.read_excel(excel_path, sheet_name='Deduplicated_Data')
    
    # Get the GCF accession from the filename
    gcf_acc = excel_path.stem.split('.calls_vs')[0]
    strain = strain_mapping.get(gcf_acc)
    
    if not strain:
        print(f"Warning: No strain mapping found for {gcf_acc}")
        return None
    
    # The methods to analyze
    methods = [
        'bakta_pseudogene',
        'pseudofinder_baktadb_pseudogene',
        'pseudofinder_salmonella_pseudogene',
        'pseudofinder_ncbi_pseudogene',
        'dbs_pseudogene'
    ]
    
    results = {
        'strain': strain,
        'gcf_acc': gcf_acc
    }
    
    for method in methods:
        ppv, sens, total_pos, true_pos = calculate_metrics(df, method, strain)
        cam_count = count_central_metabolism_pseudogenes(df, method)
        
        results.update({
            f'{method}_ppv': ppv,
            f'{method}_sensitivity': sens,
            f'{method}_total_positives': total_pos,
            f'{method}_true_positives': true_pos,
            f'{method}_cam_count': cam_count
        })
    
    return results

def analyze_all_files(todo_list, strain_mapping):
    """Analyze all Excel files in the directory"""
    results = []
    
    # Process each Excel file
    for file_path in todo_list:
        file = Path(file_path)
        result = analyze_file(file, strain_mapping)
        if result:
            results.append(result)
    
    # Create a DataFrame with all results
    results_df = pd.DataFrame(results)
    
    return results_df

# Example usage:
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

list_of_excels = ['2024.11.14/GCF_000007545.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000008105.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000009505.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000009525.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000011885.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000018385.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000018705.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000020745.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000020885.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000020925.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000026565.1.calls_vs_nuccio.xlsx', '2024.11.14/GCF_000195995.1.calls_vs_nuccio.xlsx']

results = analyze_all_files(list_of_excels, STRAIN_MAPPING)

# add a column to results based on the GCF accession lookup in ei_gi_lookup
results['salm_type'] = results['gcf_acc'].map(ei_gi_lookup)
# move salm_type to third column in dataframe

# To move column 'D' to the third position (index 2):
cols = list(results.columns)
cols.remove('salm_type')             # Remove the column name you want to move
cols.insert(2, 'salm_type')         # Insert it at position 2 (third position)
results = results[cols]               # Reorder the DataFrame

results.to_csv('2024.11.14/2024.11.14.pseudogene_validation_results.csv', index=False)