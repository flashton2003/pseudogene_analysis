import pandas as pd
import re
from pathlib import Path
import argparse
import shutil


# Read and process the nuccio file
def process_nuccio(file_path):
    df = pd.read_excel(file_path)
    # Extract UniProtKB ID using regex
    df['UniProtKB_ID'] = df['Cross-reference'].apply(
        lambda x: re.search(r'UniProtKB:([A-Z0-9]+)', str(x)).group(1) 
        if re.search(r'UniProtKB:([A-Z0-9]+)', str(x)) else None
    )
    return df

# Process the bakta GFF file
def process_bakta(file_path):
    # Read GFF file, skipping comments
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                    names=['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'attribute'],
                          low_memory=False)
    
    # Extract locus tag
    df['locus_tag'] = df['attribute'].apply(
        lambda x: re.search(r'locus_tag=([^;]+)', str(x)).group(1) 
        if re.search(r'locus_tag=([^;]+)', str(x)) else None
    )
    
    # Filter for pseudo=True
    pseudo_locus_tags = df[df['attribute'].str.contains('pseudo=True', na=False)]['locus_tag'].tolist()
    return pseudo_locus_tags

# Process the pseudofinder calls
def process_pseudofinder(file_path):
    # Read GFF file, skipping comments
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                    names=['seqname', 'source', 'feature', 'start', 'end', 
                          'score', 'strand', 'frame', 'attribute'],
                          low_memory=False)
    
    # Extract locus tags matching the pattern: 6 uppercase letters, underscore, 5 digits
    def extract_locus_tag(attribute):
        if pd.isna(attribute):
            return None
        # Find the old_locus_tag section
        old_locus_match = re.search(r'old_locus_tag=([^;]+)', str(attribute))
        if not old_locus_match:
            return None
        
        # Split multiple locus tags and find ones matching the pattern
        locus_tags = old_locus_match.group(1).split(',')
        pattern = r'[A-Z]{6}_\d{5}'
        matching_tags = [tag for tag in locus_tags if re.match(pattern, tag.strip())]
        
        # Return the first matching tag found, or None if none found
        return matching_tags[0] if matching_tags else None
    # df['matching_locus_tag', 'attribute'].to_csv('pseudofinder.csv')
    
    df['matching_locus_tag'] = df['attribute'].apply(extract_locus_tag)
    df.to_csv('pseudofinder.csv')
    matching_tags = df['matching_locus_tag'].dropna().tolist()
    # print(matching_tags)
    # print(f"Found {len(matching_tags)} matching locus tags")
    return matching_tags

# Process reciprocal diamond file
def process_diamond(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df[['qseqid_fwd', 'protein_id']]

# Process anaerobic genes file
def process_anaerobic(file_path):
    return pd.read_excel(file_path)['Reference locus tag(s)'].tolist()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pseudogene data from multiple sources')
    parser.add_argument('--nuccio', required=True, help='Path to Nuccio Excel file')
    parser.add_argument('--bakta', required=True, help='Path to Bakta GFF3 file')
    parser.add_argument('--pseudofinder', required=True, help='Path to Pseudofinder GFF file')
    parser.add_argument('--diamond', required=True, help='Path to Diamond TSV file')
    parser.add_argument('--anaerobic', required=True, help='Path to anaerobic genes Excel file')
    parser.add_argument('--output', default='processed_results.xlsx', 
                      help='Path for output Excel file (default: processed_results.xlsx)')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Read all files
    nuccio_df = process_nuccio(args.nuccio)
    bakta_pseudos = process_bakta(args.bakta)
    pseudofinder_pseudos = process_pseudofinder(args.pseudofinder)
    diamond_df = process_diamond(args.diamond)
    anaerobic_genes = process_anaerobic(args.anaerobic)
    
    # First join: nuccio and diamond
    merged_df = pd.merge(nuccio_df, diamond_df, 
                        left_on='UniProtKB_ID', 
                        right_on='protein_id', 
                        how='left')
    
    # Add bakta pseudogene column
    merged_df['bakta_pseudogene'] = merged_df['qseqid_fwd'].isin(bakta_pseudos).astype(int)
    
    # Add pseudofinder pseudogene column
    merged_df['pseudofinder_pseudogene'] = merged_df['qseqid_fwd'].isin(pseudofinder_pseudos).astype(int)
    
    # Add anaerobic metabolism column
    merged_df['central_anaerobic_metabolism'] = merged_df['Reference locus tag(s)'].isin(anaerobic_genes).astype(int)
    
    # Create combined pseudogene column
    merged_df['pseudogene'] = ((merged_df['bakta_pseudogene'] == 1) | 
                              (merged_df['pseudofinder_pseudogene'] == 1)).astype(int)
    
    # Save the results
    merged_df.to_excel(args.output, index=False, engine='openpyxl')

if __name__ == "__main__":
    main()