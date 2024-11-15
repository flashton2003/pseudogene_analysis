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

def count_group_pseudogenes(df, method_col, group_id):
    """Count pseudogenes in a specific functional group"""
    return ((df[method_col] == 1) & (df['Cross-reference'].str.contains(group_id, na=False))).sum()

def analyze_file(excel_path, strain_mapping, coord_matching):
    """Analyze a single Excel file and return its metrics"""
    # Read the Excel file
    if coord_matching is True:
        df = pd.read_excel(excel_path)
    elif coord_matching is False:
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
    
    # Define functional groups
    functional_groups = {
        'fimbrae': 'GroupID:G01',
        'T3SS-1_effector': 'GroupID:G02',
        'T3SS-2_effector': 'GroupID:G03',
        'motility_chemotaxis': 'GroupID:G05'
    }
    
    positives_in_truth = df[strain].astype(str).str.startswith('2').sum()
    positives_in_cam_truth = df[(df[strain].astype(str).str.startswith('2')) & (df['central_anaerobic_metabolism'] == 1)].shape[0]

    results = {
        'strain': strain,
        'gcf_acc': gcf_acc,
        'total_positives_in_truth': positives_in_truth,
        'total_positives_in_cam_truth': positives_in_cam_truth
    }
    
    # Count truth positives for each functional group
    for group_name, group_id in functional_groups.items():
        truth_count = df[(df[strain].astype(str).str.startswith('2')) & 
                        (df['Cross-reference'].str.contains(group_id, na=False))].shape[0]
        results[f'{group_name}_truth'] = truth_count

    for method in methods:
        ppv, sens, total_pos, true_pos = calculate_metrics(df, method, strain)
        cam_count = count_central_metabolism_pseudogenes(df, method)
        
        # Add basic metrics
        results.update({
            f'{method}_ppv': ppv,
            f'{method}_sensitivity': sens,
            f'{method}_total_positives': total_pos,
            f'{method}_true_positives': true_pos,
            f'{method}_cam_count': cam_count
        })
        
        # Add functional group counts
        for group_name, group_id in functional_groups.items():
            group_count = count_group_pseudogenes(df, method, group_id)
            results[f'{method}_{group_name}_count'] = group_count
    
    return results

def analyze_all_files(todo_list, strain_mapping, coord_matching = False):
    """Analyze all Excel files in the directory"""
    results = []
    
    # Process each Excel file
    for file_path in todo_list:
        file = Path(file_path)
        result = analyze_file(file, strain_mapping, coord_matching)
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

results = analyze_all_files(list_of_excels, STRAIN_MAPPING, coord_matching = True)

# Add salm_type column based on GCF accession lookup
results['salm_type'] = results['gcf_acc'].map(ei_gi_lookup)

# Move salm_type to third column in dataframe
cols = list(results.columns)
cols.remove('salm_type')
cols.insert(2, 'salm_type')
results = results[cols]

# Save results
results.to_csv('2024.11.14/2024.11.14.pseudogene_validation_results.diamond_matching.csv', index=False)