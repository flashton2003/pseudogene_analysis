import pandas as pd
import argparse
from typing import List, Tuple, Dict
import os
import re

def infer_file_type(file_path: str) -> str:
    _, extension = os.path.splitext(file_path.lower())
    if extension in ['.xls', '.xlsx']:
        return 'excel'
    elif extension == '.csv':
        return 'csv'
    elif extension in ['.gff', '.gff3']:
        return 'gff'
    else:
        raise ValueError(f"Unsupported file type: {extension}")

def is_overlap(truth_start: int, truth_end: int, call_start: int, call_end: int) -> bool:
    truth_length = truth_end - truth_start
    overlap_start = max(truth_start, call_start)
    overlap_end = min(truth_end, call_end)

    if overlap_start < overlap_end:
        overlap_length = overlap_end - overlap_start
        if overlap_length >= 0.1 * truth_length:
            return True
    return False

def calculate_sensitivity_ppv(truth_df: pd.DataFrame, call_df: pd.DataFrame) -> Tuple[float, float]:
    true_positives = 0
    false_negatives = 0
    false_positives = 0
    
    for _, truth_row in truth_df.iterrows():
        truth_start = int(truth_row['start'])
        truth_end = int(truth_row['stop'])
        
        matched = any(is_overlap(truth_start, truth_end, int(call_row['start']), int(call_row['end']))
                      for _, call_row in call_df.iterrows())
        
        if matched:
            true_positives += 1
        else:
            false_negatives += 1
    
    for _, call_row in call_df.iterrows():
        call_start = int(call_row['start'])
        call_end = int(call_row['end'])
        
        if not any(is_overlap(int(truth_row['start']), int(truth_row['stop']), call_start, call_end)
                   for _, truth_row in truth_df.iterrows()):
            false_positives += 1
    
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    ppv = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    
    return sensitivity, ppv

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

def read_and_filter_data(file_path: str, is_pseudofinder: bool = False, max_orf_percentage: float = None) -> pd.DataFrame:
    file_type = infer_file_type(file_path)
    if file_type == 'excel':
        df = pd.read_excel(file_path)
    elif file_type == 'csv':
        df = pd.read_csv(file_path)
    elif file_type == 'gff':
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
        df = df[df['seqname'] == 'contig_1']
    else:
        raise ValueError(f"Unsupported file type: {file_type}")
    
    if is_pseudofinder and max_orf_percentage is not None:
        df = filter_pseudofinder_calls(df, max_orf_percentage)
    
    return df

def process_datasets(truth_file: str, call_files: List[Dict[str, str]], max_orf_percentage: float = None) -> pd.DataFrame:
    truth_df = read_and_filter_data(truth_file)
    
    results = []
    all_calls = {}
    
    for call_file in call_files:
        is_pseudofinder = 'pseudofinder' in call_file['name'].lower()
        call_df = read_and_filter_data(call_file['path'], is_pseudofinder, max_orf_percentage)
        
        if call_file['name'].lower() == 'bakta_pseudo':
            call_df = call_df[call_df['attribute'].notna() & call_df['attribute'].str.contains('pseudo=True')]
        
        call_df['in_truth'] = call_df.apply(lambda row: any(is_overlap(int(truth_row['start']), int(truth_row['stop']), 
                                                                       int(row['start']), int(row['end'])) 
                                                            for _, truth_row in truth_df.iterrows()), axis=1)
        
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
        
        # Store the calls with a name that indicates if they were filtered
        calls_key = f"{call_file['name']}_filtered_{max_orf_percentage}" if max_orf_percentage is not None else call_file['name']
        all_calls[calls_key] = call_df
    
    return pd.DataFrame(results), all_calls

def main():
    parser = argparse.ArgumentParser(description="Compare pseudogene calls against a truth set.")
    parser.add_argument("--truth", required=True, help="Path to the truth file")
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
    unfiltered_results, unfiltered_calls = process_datasets(args.truth, call_files)
    
    # Run analysis with filtering
    print(f"\nRunning analysis with ORF filtering (max {args.max_orf_percentage}%)...")
    filtered_results, filtered_calls = process_datasets(args.truth, call_files, args.max_orf_percentage)
    
    # Combine results
    combined_results = pd.concat([unfiltered_results, filtered_results])
    
    print("\nCombined Results:")
    print(combined_results)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save combined results
    combined_results.to_csv(os.path.join(args.output_dir, 'summary_results.csv'), index=False)
    
    # Save unfiltered call datasets
    for name, df in unfiltered_calls.items():
        output_file = os.path.join(args.output_dir, f'{name}_unfiltered_calls.csv')
        df.to_csv(output_file, index=False)
        print(f"Saved unfiltered calls for {name} to {output_file}")
    
    # Save filtered call datasets
    for name, df in filtered_calls.items():
        output_file = os.path.join(args.output_dir, f'{name}_filtered_calls.csv')
        df.to_csv(output_file, index=False)
        print(f"Saved filtered calls for {name} to {output_file}")

if __name__ == "__main__":
    main()