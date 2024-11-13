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

def select_row_for_index(group):
    """
    Select the appropriate row for each Index group based on the specified logic:
    1. If no loss_of_function=1, take highest delta-bitscore row
    2. If highest positive delta-bitscore has loss_of_function=1, take that row
    3. If lowest delta-bitscore has loss_of_function=1, take that row
    4. Otherwise, take highest delta-bitscore row
    """
    # Sort by delta-bitscore in descending order
    sorted_group = group.sort_values('delta-bitscore', ascending=False)
    
    # Check if any rows have loss_of_function = 1
    has_loss_of_function = (sorted_group['loss_of_function'] == 1).any()
    
    if not has_loss_of_function:
        # If no loss of function, take highest delta-bitscore row
        return sorted_group.iloc[0]
    else:
        # Get rows with positive delta-bitscore and loss_of_function=1
        positive_dbs_with_loss = sorted_group[
            (sorted_group['delta-bitscore'] > 0) & 
            (sorted_group['loss_of_function'] == 1)
        ]
        
        if not positive_dbs_with_loss.empty:
            # Take the highest delta-bitscore row among positive values with loss_of_function=1
            return positive_dbs_with_loss.iloc[0]
        elif sorted_group.iloc[-1]['loss_of_function'] == 1:
            # If no positive delta-bitscore with loss_of_function=1,
            # check if lowest delta-bitscore has loss_of_function=1
            return sorted_group.iloc[-1]
        else:
            # If neither condition is met, take highest delta-bitscore row
            return sorted_group.iloc[0]
        
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
    nuccio_isangi_lookup = pd.read_csv(args.lookup, sep='\t')
    dbs_results = pd.read_csv(args.dbs, sep='\t', skiprows=1)
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
    final_df = pd.merge(merged_df, dbs_results, 
                       left_on='qseqid_fwd', right_on='gene_2', 
                       how='left')
    
    # Create the central anaerobic metabolism column
    # Convert anaerobic genes' reference locus tags to a set for faster lookup
    anaerobic_locus_tags = set(anaerobic_genes['Reference locus tag(s)'].dropna())
    
    # Create new column marking central anaerobic metabolism genes
    final_df['central anaerobic metabolism gene?'] = final_df['Reference locus tag(s)'].isin(anaerobic_locus_tags)
    
    # Create deduplicated version using the new logic
    dedup_rows = []
    for idx, group in final_df.groupby('Index'):
        selected_row = select_row_for_index(group)
        dedup_rows.append(selected_row)
    
    dedup_df = pd.DataFrame(dedup_rows)
    
    # Save both dataframes to different sheets in the same Excel file
    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        final_df.to_excel(writer, sheet_name='Complete_Data', index=False)
        dedup_df.to_excel(writer, sheet_name='Deduplicated_Data', index=False)
    
    print(f"Processing complete. Output saved to: {args.output}")
    print("Sheet names:")
    print("- Complete_Data: contains all rows")
    print("- Deduplicated_Data: contains one row per Index based on selection criteria")
    print("\nSelection criteria for deduplicated data:")
    print("1. If no loss_of_function=1 rows, take highest delta-bitscore")
    print("2. If highest positive delta-bitscore has loss_of_function=1, take that")
    print("3. If lowest delta-bitscore has loss_of_function=1, take that")
    print("4. Otherwise, take highest delta-bitscore")
    print("\nFirst few rows of the deduplicated data:")
    print(dedup_df.head())

if __name__ == "__main__":
    main()
