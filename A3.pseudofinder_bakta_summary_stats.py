import pandas as pd
import argparse

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

def evaluate_predictions(data, input_filename, strain_col):
    # Extract strain column and is_pseudogene column
    is_pseudogene_col = 'pseudogene'
    cam_col = 'central_anaerobic_metabolism'
    
    # Initialize counters
    true_positives = 0
    false_positives = 0
    false_negatives = 0
    total_analyzed = 0
    pseudo_and_cam_count = 0  # New counter for genes that are both pseudogenes and CAM
    
    # Iterate through rows
    for idx, row in data.iterrows():
        strain_value = str(row[strain_col])
        is_pseudogene = int(row[is_pseudogene_col])
        is_cam = int(row[cam_col]) if pd.notna(row[cam_col]) else 0  # Handle potential NaN values
        
        # Count genes that are both pseudogenes and CAM
        if is_pseudogene == 1 and is_cam == 1:
            pseudo_and_cam_count += 1
        
        # Skip if strain value is "Absent" or empty
        if strain_value == "3|Absent" or pd.isna(strain_value):
            continue
            
        total_analyzed += 1
        
        # Check if it's a true pseudogene in strain (starts with "2|")
        is_true_pseudogene = strain_value.startswith("2|")
        
        # Count true positives, false positives, and false negatives
        if is_pseudogene and is_true_pseudogene:
            true_positives += 1
        elif is_pseudogene and not is_true_pseudogene:
            false_positives += 1
        elif not is_pseudogene and is_true_pseudogene:
            false_negatives += 1
    
    # Calculate sensitivity and PPV
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    ppv = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    
    # Create results dictionary
    results = {
        'input_file': input_filename,
        'strain': strain_col,
        'total_analyzed': total_analyzed,
        'true_positives': true_positives,
        'false_positives': false_positives,
        'false_negatives': false_negatives,
        'sensitivity': sensitivity,
        'ppv': ppv,
        'pseudo_and_cam_count': pseudo_and_cam_count  # Add new metric to results
    }
    
    return results

def create_tsv_output(results):
    # Create header and data rows
    header = ['input_file', 'strain', 'total_analyzed', 'true_positives', 'false_positives', 
              'false_negatives', 'sensitivity', 'ppv', 'pseudo_and_cam_count']  # Added new column
    
    data_row = [
        results['input_file'],
        results['strain'],
        str(results['total_analyzed']),
        str(results['true_positives']),
        str(results['false_positives']),
        str(results['false_negatives']),
        f"{results['sensitivity']:.2%}",
        f"{results['ppv']:.2%}",
        str(results['pseudo_and_cam_count'])  # Added new metric
    ]
    
    return '\t'.join(header) + '\n' + '\t'.join(data_row) + '\n'

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Analyze pseudogene predictions.')
    parser.add_argument('--input_file', help='Input Excel file path')
    parser.add_argument('--output_file', help='Output TSV file path')
    parser.add_argument('--accession', help='Accession number to compare against')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Get strain name from accession
    if args.accession not in STRAIN_MAPPING:
        raise ValueError(f"Unknown accession: {args.accession}")
    strain = STRAIN_MAPPING[args.accession]
    
    # Read the input file
    try:
        data = pd.read_excel(args.input_file)
    except Exception as e:
        raise Exception(f"Error reading input file: {str(e)}")
    
    # Check if required columns exist
    if strain not in data.columns:
        raise ValueError(f"Strain column '{strain}' not found in input file")
    if 'central_anaerobic_metabolism' not in data.columns:
        raise ValueError("Column 'central_anaerobic_metabolism' not found in input file")
    
    # Run the analysis
    results = evaluate_predictions(data, args.input_file, strain)
    
    # Create TSV output
    tsv_output = create_tsv_output(results)
    
    # Print to screen
    print(tsv_output)
    
    # Write to file
    try:
        with open(args.output_file, 'w') as f:
            f.write(tsv_output)
        print(f"\nResults have been written to {args.output_file}")
    except Exception as e:
        raise Exception(f"Error writing output file: {str(e)}")

if __name__ == "__main__":
    main()