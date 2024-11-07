import pandas as pd
import argparse

def extract_uniprot_id(cross_ref):
    if pd.isnull(cross_ref):
        return None
    parts = cross_ref.split('|')
    for part in parts:
        if part.startswith('UniProtKB:'):
            return part.replace('UniProtKB:', '')
    return None

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process and merge Nuccio analysis with DBS results')
    parser.add_argument('--nuccio', required=True, help='Path to Nuccio analysis Excel file')
    parser.add_argument('--lookup', required=True, help='Path to reciprocal best hits lookup TSV file')
    parser.add_argument('--dbs', required=True, help='Path to DBS TSV file')
    parser.add_argument('--anaerobic', required=True, help='Path to central anaerobic genes Excel file')
    parser.add_argument('--output', required=True, help='Path for output Excel file')
    
    args = parser.parse_args()

    # Load the Excel files
    nuccio_analysis = pd.read_excel(args.nuccio)
    nuccio_isangi_lookup = pd.read_csv(args.lookup, sep = '\t')
    dbs_results = pd.read_csv(args.dbs, sep = '\t', skiprows=1)
    anaerobic_genes = pd.read_excel(args.anaerobic)
    
    # Create UniProtKB_ID column
    nuccio_analysis['UniProtKB_ID'] = nuccio_analysis['Cross-reference'].apply(extract_uniprot_id)
    
    # select only the columns we need from nuccio_isangi_lookup
    nuccio_isangi_lookup = nuccio_isangi_lookup[['qseqid_fwd', 'protein_id']]

    # Join nuccio_analysis with nuccio_isangi_lookup using qseqid_fwd instead of protein_id
    merged_df = pd.merge(nuccio_analysis, nuccio_isangi_lookup, 
                        left_on='UniProtKB_ID', right_on='protein_id', 
                        how='left')
    
    # Join with dbs_results
    #merged_df.to_csv('merged_df.csv', sep =',')
    #dbs_results.to_csv('dbs_results.csv')
    final_df = pd.merge(merged_df, dbs_results, 
                       left_on='qseqid_fwd', right_on='gene_2', 
                       how='left')
    
    # Create the central anaerobic metabolism column
    # Convert anaerobic genes' reference locus tags to a set for faster lookup
    anaerobic_locus_tags = set(anaerobic_genes['Reference locus tag(s)'].dropna())
    
    # Create new column marking central anaerobic metabolism genes
    final_df['central anaerobic metabolism gene?'] = final_df['Reference locus tag(s)'].isin(anaerobic_locus_tags)
    
    # Save the final dataframe
    final_df.to_excel(args.output, index=False)
    
    print(f"Processing complete. Output saved to: {args.output}")
    print("\nFirst few rows of the final dataframe:")
    print(final_df.head())

if __name__ == "__main__":
    main()
