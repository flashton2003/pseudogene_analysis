import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process DIAMOND reciprocal best hits.')
    parser.add_argument('--forward', required=True,
                      help='Path to forward comparison file (query vs subject)')
    parser.add_argument('--reverse', required=True,
                      help='Path to reverse comparison file (subject vs query)')
    parser.add_argument('--output', required=True,
                      help='Path for output TSV file of reciprocal best hits')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Load the files into DataFrames
    df_forward = pd.read_csv(args.forward, sep='\t')
    df_reverse = pd.read_csv(args.reverse, sep='\t')

    # Extract protein ID from sseqid by splitting on '|' and taking the second element ([1])
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

    # Identify the best hits (highest bitscore) for each query using idxmax()
    best_hits_forward = filtered_forward.loc[filtered_forward.groupby('qseqid')['bitscore'].idxmax()]
    best_hits_reverse = filtered_reverse.loc[filtered_reverse.groupby('qseqid')['bitscore'].idxmax()]

    print(best_hits_forward)
    print(best_hits_reverse)

    # Merge to find reciprocal best hits
    reciprocal_hits = best_hits_forward.merge(
        best_hits_reverse,
        left_on=['qseqid', 'protein_id'],
        right_on=['sseqid', 'protein_id'],
        suffixes=('_fwd', '_rev')
    )

    # Display the reciprocal best hits
    print("Reciprocal Best Hits:")
    print(reciprocal_hits)
    
    # Save results
    reciprocal_hits.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()