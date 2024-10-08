import pandas as pd
import argparse
from typing import List, Tuple, Dict
import os

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

def read_and_filter_data(file_path: str) -> pd.DataFrame:
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
    return df

def process_datasets(truth_file: str, call_files: List[Dict[str, str]]) -> pd.DataFrame:
    truth_df = read_and_filter_data(truth_file)
    
    results = []
    for call_file in call_files:
        call_df = read_and_filter_data(call_file['path'])
        
        if call_file['name'].lower() == 'bakta_pseudo':
            call_df = call_df[call_df['attribute'].notna() & call_df['attribute'].str.contains('pseudo=True')]
        
        sensitivity, ppv = calculate_sensitivity_ppv(truth_df, call_df)
        results.append({
            "Dataset": call_file['name'],
            "Sensitivity": sensitivity,
            "PPV": ppv
        })
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Compare pseudogene calls against a truth set.")
    parser.add_argument("--truth", required=True, help="Path to the truth file")
    parser.add_argument("--calls", required=True, nargs='+', help="Paths to call files")
    parser.add_argument("--call_names", required=True, nargs='+', help="Names for each call dataset")
    
    args = parser.parse_args()
    
    if len(args.calls) != len(args.call_names):
        raise ValueError("The number of call files must match the number of call names")
    
    call_files = [{'path': path, 'name': name} 
                  for path, name in zip(args.calls, args.call_names)]
    
    results_df = process_datasets(args.truth, call_files)
    print(results_df)

if __name__ == "__main__":
    main()