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
    return df[['qseqid_fwd', 'protein_id']]

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
    # print(threshold, row['delta-bitscore'])
    if pd.isnull(row['delta-bitscore']):
        return 0
    return 1 if (row['delta-bitscore'] > threshold) else 0

def select_row_for_index(group):
    """
    Select the appropriate row for each Index group based on DBS logic:
    1. If no loss_of_function=1, take highest delta-bitscore row
    2. If highest positive delta-bitscore has loss_of_function=1, take that row
    3. If lowest delta-bitscore has loss_of_function=1, take that row
    4. Otherwise, take highest delta-bitscore row
    """
    sorted_group = group.sort_values('delta-bitscore', ascending=False)
    has_loss_of_function = (sorted_group['loss_of_function'] == 1).any()
    
    if not has_loss_of_function:
        return sorted_group.iloc[0]
    
    positive_dbs_with_loss = sorted_group[
        (sorted_group['delta-bitscore'] > 0) & 
        (sorted_group['loss_of_function'] == 1)
    ]
    
    if not positive_dbs_with_loss.empty:
        return positive_dbs_with_loss.iloc[0]
    elif sorted_group.iloc[-1]['loss_of_function'] == 1:
        return sorted_group.iloc[-1]
    else:
        return sorted_group.iloc[0]

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
        merged_df['bakta_pseudogene'] = merged_df['qseqid_fwd'].isin(bakta_pseudos).astype(int)
    
    # Process each Pseudofinder result
    pseudofinder_results = {
        'baktadb': process_pseudofinder(args.pseudofinder_baktadb),
        'salmonella': process_pseudofinder(args.pseudofinder_salmonella),
        'ncbi': process_pseudofinder(args.pseudofinder_ncbi)
    }
    
    # Add columns for each Pseudofinder result
    for source, pseudos in pseudofinder_results.items():
        merged_df[f'pseudofinder_{source}_pseudogene'] = merged_df['qseqid_fwd'].isin(pseudos).astype(int)
    
    
    
    # Process DBS if provided
    if args.dbs:
        dbs_results = process_dbs(args.dbs)
        
        # do this the same as the above, filter the DBS and then return a list for
        # checking if it is in?
        merged_df = pd.merge(merged_df, dbs_results[['gene_2', 'delta-bitscore', 'loss_of_function']],
                           left_on='qseqid_fwd',
                           right_on='gene_2',
                           how='left')
        # merged_df = pd.merge(merged_df, dbs_results,
        #                    left_on='qseqid_fwd',
        #                    right_on='gene_2',
        #                    how='left')
        
        # Create deduplicated version using DBS logic
        dedup_rows = []
        for idx, group in merged_df.groupby('Index'):
            selected_row = select_row_for_index(group)
            dedup_rows.append(selected_row)
        dedup_df = pd.DataFrame(dedup_rows)
        
        # Calculate DBS threshold and add DBS pseudogene column
        dbs_threshold = calculate_dbs_threshold(dedup_df['delta-bitscore'])        
        dedup_df['dbs_pseudogene'] = dedup_df.apply(
            lambda row: is_dbs_pseudogene(row, dbs_threshold), axis=1
        )



        dedup_df.to_csv('dedup.csv')
        # Add DBS pseudogene column with default of 0
        dedup_df['dbs_pseudogene'] = dedup_df['dbs_pseudogene'].fillna(0).astype(int)


        merged_df.drop(['gene_2', 'delta-bitscore', 'loss_of_function'], axis=1, inplace=True)
        dedup_df.drop(['gene_2', 'delta-bitscore', 'loss_of_function'], axis=1, inplace=True)
        
        # Save both complete and deduplicated data
        with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
            merged_df.to_excel(writer, sheet_name='Complete_Data', index=False)
            dedup_df.to_excel(writer, sheet_name='Deduplicated_Data', index=False)
            
        print(f"Processing complete. Output saved to: {args.output}")
        print(f"DBS threshold (97.5th percentile): {dbs_threshold:.2f}")
        print("Sheets created: Complete_Data and Deduplicated_Data")
    else:
        # Save only the complete data if no DBS file provided
        merged_df.to_excel(args.output, index=False, engine='openpyxl')
        print(f"Processing complete. Output saved to: {args.output}")

if __name__ == "__main__":
    main()