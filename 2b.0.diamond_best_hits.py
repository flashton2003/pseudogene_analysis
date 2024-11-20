import pandas as pd
import argparse
import subprocess
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Perform one-way DIAMOND search analysis.')
    parser.add_argument('--query_fasta', required=True,
                      help='Path to query protein FASTA file')
    parser.add_argument('--subject_fasta', required=True,
                      help='Path to subject protein FASTA file')
    parser.add_argument('--output', required=True,
                      help='Path for output TSV file of best hits')
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

def process_diamond_results(results_file, output_file):
    """Process DIAMOND results to find best hits."""
    # Load the file into DataFrame with the correct column names
    column_names = [
        'qseqid', 'qlen', 'sseqid', 'slen', 'pident', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore', 'gaps'
    ]
    
    df = pd.read_csv(results_file, sep='\t', names=column_names)

    # Extract protein ID from sseqid by splitting on '|' and taking the second element
    df['protein_id'] = df['sseqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)

    # Calculate query coverage and subject coverage
    df['query_cov'] = (df['length'] - df['gaps']) / df['qlen'] * 100
    df['subject_cov'] = (df['length'] - df['gaps']) / df['slen'] * 100

    # Filter based on query coverage, and evalue
    filtered_results = df[
        (df['query_cov'] > 70) &
        (df['evalue'] < 1e-10)
    ]

    # Identify the best hits for each query
    best_hits = filtered_results.loc[filtered_results.groupby('qseqid')['bitscore'].idxmax()]

    # Save results
    best_hits.to_csv(output_file, sep='\t', index=False)
    return best_hits

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Check if DIAMOND is installed
    check_diamond_installation()
    
    # Create temporary directory if it doesn't exist
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    # Define paths for temporary files
    subject_db = os.path.join(args.tmp_dir, "subject_db")
    search_results = os.path.join(args.tmp_dir, "search_results.tsv")
    
    # Create DIAMOND database
    print("Creating DIAMOND database...")
    create_diamond_db(args.subject_fasta, subject_db, args.threads)
    
    # Run DIAMOND search
    print("Running DIAMOND search...")
    run_diamond_search(args.query_fasta, subject_db, search_results, args.threads)
    
    # Process results
    print("Processing best hits...")
    best_hits = process_diamond_results(search_results, args.output)
    
    print(f"Found {len(best_hits)} best hits")
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()