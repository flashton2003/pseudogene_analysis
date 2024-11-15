import pandas as pd
import re
import argparse

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

def process_anaerobic(file_path):
    """Process anaerobic genes file"""
    return pd.read_excel(file_path)['Reference locus tag(s)'].tolist()

def create_consensus_row(group):
    """Create a consensus row from a group of rows"""
    consensus = group.iloc[0].copy()
    
    if 'pseudofinder_baktadb_pseudogene' in group.columns:
        consensus['pseudofinder_baktadb_pseudogene'] = 1 if (group['pseudofinder_baktadb_pseudogene'] == 1).any() else 0
    
    return consensus

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pseudogene data from Pseudofinder BaktaDB')
    parser.add_argument('--nuccio', required=True, help='Path to Nuccio Excel file')
    parser.add_argument('--pseudofinder-baktadb', required=True, help='Path to Pseudofinder GFF file (BaktaDB)')
    parser.add_argument('--diamond', required=True, help='Path to Diamond TSV file')
    parser.add_argument('--anaerobic', required=True, help='Path to anaerobic genes Excel file')
    parser.add_argument('--output', required=True, help='Path for output Excel file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Process files
    nuccio_df = process_nuccio(args.nuccio)
    diamond_df = process_diamond(args.diamond)
    pseudofinder_baktadb = process_pseudofinder(args.pseudofinder_baktadb)
    anaerobic_genes = process_anaerobic(args.anaerobic)
    
    # Merge nuccio and diamond data
    merged_df = pd.merge(nuccio_df, diamond_df, 
                        left_on='UniProtKB_ID', 
                        right_on='protein_id', 
                        how='left')
    
    # Add Pseudofinder BaktaDB results
    merged_df['pseudofinder_baktadb_pseudogene'] = merged_df['qseqid'].isin(pseudofinder_baktadb).astype(int)
    
    # Add anaerobic metabolism column
    merged_df['central_anaerobic_metabolism'] = merged_df['Reference locus tag(s)'].isin(anaerobic_genes).astype(int)
    
    # Create consensus-based deduplicated dataframe
    dedup_df = merged_df.groupby('Index').apply(create_consensus_row)
    dedup_df = dedup_df.reset_index(drop=True)
    
    # Save both complete and deduplicated data
    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        merged_df.to_excel(writer, sheet_name='Complete_Data', index=False)
        dedup_df.to_excel(writer, sheet_name='Deduplicated_Data', index=False)
    
    print(f"Processing complete. Output saved to: {args.output}")

if __name__ == "__main__":
    main()