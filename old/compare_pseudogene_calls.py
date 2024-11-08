import pandas as pd
import argparse
from typing import List, Tuple, Dict
import os
import re

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

SEQNAME_MAPPING = {
    'GCF_000007545.1': {'contig_1': 'NC_004631.1'},
    'GCF_000008105.1': {'contig_1': 'NC_006905.1', 'contig_2': 'NC_006855.1', 'contig_3': 'NC_006856.1'},
    'GCF_000009505.1': {'contig_1': 'NC_011294.1'},
    'GCF_000009525.1': {'contig_1': 'NC_011274.1'},
    'GCF_000011885.1': {'contig_1': 'NC_006511.1'},
    'GCF_000018385.1': {'contig_1': 'NC_012125.1', 'contig_2': 'NC_012124.1'},
    'GCF_000018705.1': {'contig_1': 'NC_010102.1'},
    'GCF_000020705.1': {'contig_1': 'NC_011083.1', 'contig_2': 'NC_011082.1', 'contig_3': 'NC_011081.1'},
    'GCF_000020745.1': {'contig_1': 'NC_011094.1', 'contig_2': 'NC_011092.1', 'contig_3': 'NC_011093.1'},
    'GCF_000020885.1': {'contig_1': 'NC_011149.1', 'contig_2': 'NC_011148.1'},
    'GCF_000020925.1': {'contig_1': 'NC_011205.1', 'contig_2': 'NC_011204.1'},
    'GCF_000026565.1': {'contig_1': 'NC_011147.1'},
    'GCF_000195995.1': {'contig_1': 'NC_003198.1', 'contig_2': 'NC_003384.1', 'contig_3': 'NC_003385.1'}
}

def sanitize_coordinate(value: str) -> float:
    """
    Clean coordinate values by removing non-numeric characters (except minus sign and decimal point)
    and converting to a float.
    """
    if pd.isna(value):
        return 0.0
    # Allow only digits, minus sign, and decimal point
    cleaned = re.sub(r'[^\d.-]', '', str(value))
    try:
        return float(cleaned)
    except ValueError:
        return 0.0  # Return 0.0 if conversion fails


def sanitize_seqname(value: str) -> str:
    """
    Clean seqname by removing unwanted characters and standardizing format
    """
    if pd.isna(value):
        return ''
    # Remove '>' and any leading/trailing whitespace
    cleaned = str(value).strip().lstrip('>')
    return cleaned

def convert_seqname(seqname: str, genome_accession: str) -> str:
    """
    Convert seqname based on the genome accession using SEQNAME_MAPPING
    """
    cleaned_seqname = sanitize_seqname(seqname)
    if genome_accession in SEQNAME_MAPPING and cleaned_seqname in SEQNAME_MAPPING[genome_accession]:
        return SEQNAME_MAPPING[genome_accession][cleaned_seqname]
    return cleaned_seqname

def infer_file_type(file_path: str) -> str:
    """
    Infer the file type from the file extension
    """
    _, extension = os.path.splitext(file_path.lower())
    if extension in ['.xls', '.xlsx']:
        return 'excel'
    elif extension == '.csv':
        return 'csv'
    elif extension in ['.gff', '.gff3']:
        return 'gff'
    else:
        raise ValueError(f"Unsupported file type: {extension}")

def is_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Helper function to check if two regions overlap"""
    return start1 <= end2 and end1 >= start2

def extract_orf_percentage(attribute: str) -> float:
    """Extract the ORF percentage from the attribute string."""
    match = re.search(r'ORF is (\d+\.?\d*)%', attribute)
    if match:
        return float(match.group(1))
    return None

def filter_pseudofinder_calls(df: pd.DataFrame, max_orf_percentage: float = 100.0) -> pd.DataFrame:
    """Filter pseudofinder calls based on ORF percentage."""
    if 'attribute' not in df.columns:
        return df
    
    def should_keep_row(row):
        if 'ORF is' not in str(row['attribute']):
            return True
        
        orf_percentage = extract_orf_percentage(row['attribute'])
        if orf_percentage is None:
            return True
        
        return orf_percentage <= max_orf_percentage
    
    return df[df.apply(should_keep_row, axis=1)]

def process_truth_data(file_path: str, genome_accession: str) -> pd.DataFrame:
    """
    Process truth data from Excel file based on strain mapping
    """
    # Get strain name from genome accession
    strain_name = STRAIN_MAPPING.get(genome_accession)
    if not strain_name:
        raise ValueError(f"No strain mapping found for accession: {genome_accession}")
    
    # Read Excel file
    truth_df = pd.read_excel(file_path)
    
    # Filter rows where the strain column starts with '2'
    strain_col = truth_df[strain_name].astype(str)
    truth_df = truth_df[strain_col.str.startswith('2')]
    
    # Process coordinates column
    coord_splits = truth_df[strain_name].str.split('|')
    
    # Create new dataframe with processed coordinates
    processed_df = pd.DataFrame({
        'seqname': coord_splits.str[2].apply(sanitize_seqname),
        'start': coord_splits.str[3].apply(sanitize_coordinate),
        'stop': coord_splits.str[4].apply(sanitize_coordinate)
    })
    
    # Convert seqnames
    processed_df['seqname'] = processed_df['seqname'].apply(lambda x: convert_seqname(x, genome_accession))
    
    return processed_df

def read_and_filter_data(file_path: str, genome_accession: str, is_pseudofinder: bool = False, max_orf_percentage: float = None) -> pd.DataFrame:
    """
    Read and filter data from input files
    """
    file_type = infer_file_type(file_path)
    
    if file_type == 'excel':
        df = pd.read_excel(file_path)
    elif file_type == 'csv':
        df = pd.read_csv(file_path)
    elif file_type == 'gff':
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")
    
    # Sanitize input data
    if 'seqname' in df.columns:
        df['seqname'] = df['seqname'].apply(sanitize_seqname)
        df['seqname'] = df['seqname'].apply(lambda x: convert_seqname(x, genome_accession))
    
    # print(df['start'])
    if 'start' in df.columns:
        df['start'] = df['start'].apply(sanitize_coordinate)
    # print(df['start'])
    if 'end' in df.columns:
        df['end'] = df['end'].apply(sanitize_coordinate)
    
    if is_pseudofinder and max_orf_percentage is not None:
        df = filter_pseudofinder_calls(df, max_orf_percentage)
    
    return df

def calculate_sensitivity_ppv(truth_df: pd.DataFrame, call_df: pd.DataFrame) -> Tuple[float, float]:
    """
    Calculate sensitivity and PPV based on truth and call data
    """
    true_positives = 0
    false_negatives = 0
    false_positives = 0
    
    for _, truth_row in truth_df.iterrows():
        truth_start = truth_row['start']
        truth_end = truth_row['stop']
        truth_seqname = truth_row['seqname']
        
        matched = any(
            is_overlap(truth_start, truth_end, call_row['start'], call_row['end'])
            and truth_seqname == call_row['seqname']
            for _, call_row in call_df.iterrows()
        )
        
        if matched:
            true_positives += 1
        else:
            false_negatives += 1
    
    for _, call_row in call_df.iterrows():
        call_start = call_row['start']
        call_end = call_row['end']
        call_seqname = call_row['seqname']
        
        if not any(
            is_overlap(truth_row['start'], truth_row['stop'], call_start, call_end)
            and call_seqname == truth_row['seqname']
            for _, truth_row in truth_df.iterrows()
        ):
            false_positives += 1
    
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    ppv = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    
    return sensitivity, ppv

def process_datasets(truth_file: str, genome_accession: str, call_files: List[Dict[str, str]], max_orf_percentage: float = None) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Process all datasets and calculate metrics
    """
    truth_df = process_truth_data(truth_file, genome_accession)
    
    results = []
    all_calls = {}
    
    for call_file in call_files:
        is_pseudofinder = 'pseudofinder' in call_file['name'].lower()
        call_df = read_and_filter_data(call_file['path'], genome_accession, is_pseudofinder, max_orf_percentage)
        
        if call_file['name'].lower() == 'bakta_pseudo':
            call_df = call_df[call_df['attribute'].notna() & call_df['attribute'].str.contains('pseudo=True')]
            # print(call_df)
        truth_df.to_csv('~/Desktop/truth.csv')
        sensitivity, ppv = calculate_sensitivity_ppv(truth_df, call_df)
        
        # Add filter status to dataset name if it's pseudofinder
        dataset_name = call_file['name']
        if is_pseudofinder and max_orf_percentage is not None:
            dataset_name = f"{dataset_name} (ORF≤{max_orf_percentage}%)"
        
        results.append({
            "Dataset": dataset_name,
            "Sensitivity": sensitivity,
            "PPV": ppv,
            "Total_Calls": len(call_df)
        })
        
        calls_key = f"{call_file['name']}_filtered_{max_orf_percentage}" if max_orf_percentage is not None else call_file['name']
        all_calls[calls_key] = call_df
    
    return pd.DataFrame(results), all_calls

def main():
    parser = argparse.ArgumentParser(description="Compare pseudogene calls against a truth set.")
    parser.add_argument("--truth", required=True, help="Path to the truth Excel file")
    parser.add_argument("--genome_accession", required=True, help="Genome accession number (e.g., GCF_000020705.1)")
    parser.add_argument("--calls", required=True, nargs='+', help="Paths to call files")
    parser.add_argument("--call_names", required=True, nargs='+', help="Names for each call dataset")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    parser.add_argument("--max_orf_percentage", type=float, default=100.0, 
                      help="Maximum ORF percentage to include for pseudofinder calls (default: 100.0)")
    
    args = parser.parse_args()
    
    if len(args.calls) != len(args.call_names):
        raise ValueError("The number of call files must match the number of call names")
    
    call_files = [{'path': path, 'name': name} 
                  for path, name in zip(args.calls, args.call_names)]
    
    print("Processing files...")
    print(call_files)

    # Run analysis without filtering
    print("\nRunning analysis without ORF filtering...")
    unfiltered_results, unfiltered_calls = process_datasets(args.truth, args.genome_accession, call_files)
    
    # Run analysis with filtering
    print(f"\nRunning analysis with ORF filtering (max {args.max_orf_percentage}%)...")
    filtered_results, filtered_calls = process_datasets(args.truth, args.genome_accession, call_files, args.max_orf_percentage)
    
    # Combine results
    combined_results = pd.concat([unfiltered_results, filtered_results])
    
    print("\nCombined Results:")
    print(combined_results)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # add a column for the genome accession
    combined_results['genome_accession'] = args.genome_accession

    # Save combined results
    combined_results.to_csv(os.path.join(args.output_dir, f'{args.genome_accession}_summary_results.csv'), index=False)
    
    # the below isn't needed, as is same as the input.
    # # Save unfiltered call datasets
    # for name, df in unfiltered_calls.items():
    #     output_file = os.path.join(args.output_dir, f'{args.genome_accession}_{name}_unfiltered_calls.csv')
    #     df.to_csv(output_file, index=False)
    #     print(f"Saved unfiltered calls for {name} to {output_file}")
    
    # # Save filtered call datasets
    # for name, df in filtered_calls.items():
    #     output_file = os.path.join(args.output_dir, f'{args.genome_accession}_{name}_filtered_calls.csv')
    #     df.to_csv(output_file, index=False)
    #     print(f"Saved filtered calls for {name} to {output_file}")

if __name__ == "__main__":
    main()