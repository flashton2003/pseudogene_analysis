import pandas as pd
import re
from pathlib import Path
import argparse
import numpy as np

def extract_uniprot_id(cross_ref):
    """Extract UniProtKB ID from cross-reference string using regex"""
    if pd.isnull(cross_ref):
        return None
    match = re.search(r'UniProtKB:([A-Z0-9]+)', str(cross_ref))
    return match.group(1) if match else None

def process_nuccio(file_path):
    """Read and process the nuccio file"""
    df = pd.read_excel(file_path)
    df['UniProtKB_ID'] = df['Cross-reference'].apply(extract_uniprot_id)
    return df

def process_bakta(file_path):
    """Process the bakta GFF file to get pseudo genes"""
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                    names=['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'attribute'],
                    low_memory=False)
    
    df['locus_tag'] = df['attribute'].apply(
        lambda x: re.search(r'locus_tag=([^;]+)', str(x)).group(1) 
        if re.search(r'locus_tag=([^;]+)', str(x)) else None
    )
    
    return df[df['attribute'].str.contains('pseudo=True', na=False)]['locus_tag'].tolist()

def process_pseudofinder(file_path):
    """Process the pseudofinder calls to get matching locus tags"""
    if not file_path:
        return []
        
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                    names=['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'attribute'],
                    low_memory=False)
    
    def extract_locus_tag(attribute):
        if pd.isna(attribute):
            return None
        old_locus_match = re.search(r'old_locus_tag=([^;]+)', str(attribute))
        if not old_locus_match:
            return None
        
        locus_tags = old_locus_match.group(1).split(',')
        pattern = r'[A-Z]{6}_\d{5}'
        matching_tags = [tag for tag in locus_tags if re.match(pattern, tag.strip())]
        return matching_tags[0] if matching_tags else None
    
    df['matching_locus_tag'] = df['attribute'].apply(extract_locus_tag)
    return df['matching_locus_tag'].dropna().tolist()

def process_diamond(file_path):
    """Process reciprocal diamond file"""
    df = pd.read_csv(file_path, sep='\t')
    return df[['qseqid', 'protein_id']]

def process_dbs(file_path):
    """Process DBS results file"""
    return pd.read_csv(file_path, sep='\t', skiprows=1)

def process_anaerobic(file_path):
    """Process anaerobic genes file"""
    return pd.read_excel(file_path)['Reference locus tag(s)'].tolist()

def calculate_dbs_threshold(dbs_values):
    """Calculate the 97.5th percentile threshold for delta-bit-scores"""
    return np.percentile(dbs_values.dropna(), 97.5)

def is_dbs_pseudogene(row, threshold):
    """
    Determine if a gene is a pseudogene based on DBS criteria:
    1. Delta-bit-score > threshold OR
    2. loss_of_function = 1
    """
    if pd.isnull(row['delta-bitscore']):
        return 0
    return 1 if (row['delta-bitscore'] > threshold) else 0

def create_consensus_row(group):
    """
    Create a consensus row from a group of rows based on pseudogene columns.
    If any row has a 1 in a pseudogene column, the consensus row gets a 1.
    """
    # Get the first row as the base for non-pseudogene columns
    consensus = group.iloc[0].copy()
    
    # List of pseudogene columns to check
    pseudo_columns = [
        'bakta_pseudogene',
        'pseudofinder_baktadb_pseudogene',
        'pseudofinder_salmonella_pseudogene',
        'pseudofinder_ncbi_pseudogene',
        'dbs_pseudogene'
    ]
    
    # For each pseudogene column, set to 1 if any row in the group has a 1
    for col in pseudo_columns:
        if col in group.columns:  # Check if column exists
            consensus[col] = 1 if (group[col] == 1).any() else 0
    
    return pd.Series(consensus)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pseudogene data from multiple sources')
    parser.add_argument('--nuccio', required=True, help='Path to Nuccio Excel file')
    parser.add_argument('--bakta', help='Path to Bakta GFF3 file')
    parser.add_argument('--pseudofinder-baktadb', help='Path to Pseudofinder GFF file (BaktaDB)')
    parser.add_argument('--pseudofinder-salmonella', help='Path to Pseudofinder GFF file (Salmonella)')
    parser.add_argument('--pseudofinder-ncbi', help='Path to Pseudofinder GFF file (NCBI)')
    parser.add_argument('--diamond', required=True, help='Path to Diamond TSV file')
    parser.add_argument('--dbs', help='Path to DBS TSV file')
    parser.add_argument('--anaerobic', required=True, help='Path to anaerobic genes Excel file')
    parser.add_argument('--output', required=True, help='Path for output Excel file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Process common files
    nuccio_df = process_nuccio(args.nuccio)
    diamond_df = process_diamond(args.diamond)
    anaerobic_genes = process_anaerobic(args.anaerobic)
    
    # Initial merge of nuccio and diamond data
    merged_df = pd.merge(nuccio_df, diamond_df, 
                        left_on='UniProtKB_ID', 
                        right_on='protein_id', 
                        how='left')
    
    # Add anaerobic metabolism column
    merged_df['central_anaerobic_metabolism'] = merged_df['Reference locus tag(s)'].isin(anaerobic_genes).astype(int)
    
    # Process Bakta if provided
    if args.bakta:
        bakta_pseudos = process_bakta(args.bakta)
        merged_df['bakta_pseudogene'] = merged_df['qseqid'].isin(bakta_pseudos).astype(int)
    
    # Process each Pseudofinder result
    pseudofinder_results = {
        'baktadb': process_pseudofinder(args.pseudofinder_baktadb),
        'salmonella': process_pseudofinder(args.pseudofinder_salmonella),
        'ncbi': process_pseudofinder(args.pseudofinder_ncbi)
    }
    
    # Add columns for each Pseudofinder result
    for source, pseudos in pseudofinder_results.items():
        merged_df[f'pseudofinder_{source}_pseudogene'] = merged_df['qseqid'].isin(pseudos).astype(int)
    
    # Process DBS if provided
    if args.dbs:
        dbs_results = process_dbs(args.dbs)
        
        # Merge DBS results
        merged_df = pd.merge(merged_df, dbs_results[['gene_2', 'delta-bitscore', 'loss_of_function']],
                           left_on='qseqid',
                           right_on='gene_2',
                           how='left')
        
        # Calculate DBS threshold and add DBS pseudogene column
        dbs_threshold = calculate_dbs_threshold(merged_df['delta-bitscore'])
        merged_df['dbs_pseudogene'] = merged_df.apply(
            lambda row: is_dbs_pseudogene(row, dbs_threshold), axis=1
        )
        
        # Clean up temporary DBS columns
        merged_df.drop(['gene_2', 'delta-bitscore', 'loss_of_function'], axis=1, inplace=True)
    
    # Create consensus-based deduplicated dataframe
    dedup_df = (merged_df.groupby('Index', group_keys=True)
                .apply(create_consensus_row, include_groups=False)
                .reset_index())
        
    # Save both complete and deduplicated data
    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        merged_df.to_excel(writer, sheet_name='Complete_Data', index=False)
        dedup_df.to_excel(writer, sheet_name='Deduplicated_Data', index=False)
        
    print(f"Processing complete. Output saved to: {args.output}")
    if args.dbs:
        print(f"DBS threshold (97.5th percentile): {dbs_threshold:.2f}")
    print("Sheets created: Complete_Data and Deduplicated_Data")

if __name__ == "__main__":
    main()