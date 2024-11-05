import pandas as pd
import argparse
import subprocess
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Perform DIAMOND reciprocal best hits analysis.')
    parser.add_argument('--query_fasta', required=True,
                      help='Path to query protein FASTA file')
    parser.add_argument('--subject_fasta', required=True,
                      help='Path to subject protein FASTA file')
    parser.add_argument('--output', required=True,
                      help='Path for output TSV file of reciprocal best hits')
    parser.add_argument('--threads', type=int, default=4,
                      help='Number of threads for DIAMOND to use')
    parser.add_argument('--tmp_dir', default='diamond_tmp',
                      help='Directory for temporary files')
    return parser.parse_args()

def check_diamond_installation():
    """Check if DIAMOND is installed and accessible."""
    try:
        subprocess.run(['diamond', 'version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.exit("Error: DIAMOND is not installed or not in PATH. Please install DIAMOND first.")

def create_diamond_db(fasta_file, output_db, threads):
    """Create a DIAMOND database from a FASTA file."""
    cmd = [
        'diamond', 'makedb',
        '--in', fasta_file,
        '--db', output_db,
        '--threads', str(threads)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error creating DIAMOND database: {e.stderr}")

def run_diamond_search(query_fasta, db_path, output_file, threads):
    """Run DIAMOND search with sensitive mode and specified output format."""
    cmd = [
        'diamond', 'blastp',
        '--query', query_fasta,
        '--db', db_path,
        '--sensitive',
        '--out', output_file,
        '--outfmt', '6', 
        'qseqid', 'qlen', 'sseqid', 'slen', 'pident', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore', 'gaps',
        '--threads', str(threads)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error running DIAMOND search: {e.stderr}")

def process_diamond_results(forward_file, reverse_file, output_file):
    """Process DIAMOND results to find reciprocal best hits."""
    # Load the files into DataFrames with the correct column names
    column_names = [
        'qseqid', 'qlen', 'sseqid', 'slen', 'pident', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore', 'gaps'
    ]
    
    df_forward = pd.read_csv(forward_file, sep='\t', names=column_names)
    df_reverse = pd.read_csv(reverse_file, sep='\t', names=column_names)

    # Extract protein ID from sseqid by splitting on '|' and taking the second element
    df_forward['protein_id'] = df_forward['sseqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)
    df_reverse['protein_id'] = df_reverse['qseqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)

    # Calculate query coverage and subject coverage
    df_forward['query_cov'] = (df_forward['length'] - df_forward['gaps']) / df_forward['qlen'] * 100
    df_forward['subject_cov'] = (df_forward['length'] - df_forward['gaps']) / df_forward['slen'] * 100

    df_reverse['query_cov'] = (df_reverse['length'] - df_reverse['gaps']) / df_reverse['qlen'] * 100
    df_reverse['subject_cov'] = (df_reverse['length'] - df_reverse['gaps']) / df_reverse['slen'] * 100

    # Filter based on query coverage, subject coverage, and evalue
    filtered_forward = df_forward[
        (df_forward['query_cov'] > 70) &
        (df_forward['subject_cov'] > 70) &
        (df_forward['evalue'] < 1e-15)
    ]

    filtered_reverse = df_reverse[
        (df_reverse['query_cov'] > 70) &
        (df_reverse['subject_cov'] > 70) &
        (df_reverse['evalue'] < 1e-15)
    ]

    # Identify the best hits for each query
    best_hits_forward = filtered_forward.loc[filtered_forward.groupby('qseqid')['bitscore'].idxmax()]
    best_hits_reverse = filtered_reverse.loc[filtered_reverse.groupby('qseqid')['bitscore'].idxmax()]

    # Merge to find reciprocal best hits
    reciprocal_hits = best_hits_forward.merge(
        best_hits_reverse,
        left_on=['qseqid', 'protein_id'],
        right_on=['sseqid', 'protein_id'],
        suffixes=('_fwd', '_rev')
    )

    # Save results
    reciprocal_hits.to_csv(output_file, sep='\t', index=False)
    return reciprocal_hits

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Check if DIAMOND is installed
    check_diamond_installation()
    
    # Create temporary directory if it doesn't exist
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    # Define paths for temporary files
    query_db = os.path.join(args.tmp_dir, "query_db")
    subject_db = os.path.join(args.tmp_dir, "subject_db")
    forward_results = os.path.join(args.tmp_dir, "forward_results.tsv")
    reverse_results = os.path.join(args.tmp_dir, "reverse_results.tsv")
    
    # Create DIAMOND databases
    print("Creating DIAMOND databases...")
    create_diamond_db(args.query_fasta, query_db, args.threads)
    create_diamond_db(args.subject_fasta, subject_db, args.threads)
    
    # Run DIAMOND searches
    print("Running forward DIAMOND search...")
    run_diamond_search(args.query_fasta, subject_db, forward_results, args.threads)
    
    print("Running reverse DIAMOND search...")
    run_diamond_search(args.subject_fasta, query_db, reverse_results, args.threads)
    
    # Process results
    print("Processing reciprocal best hits...")
    reciprocal_hits = process_diamond_results(forward_results, reverse_results, args.output)
    
    print(f"Found {len(reciprocal_hits)} reciprocal best hits")
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()