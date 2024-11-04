import pandas as pd

# Load the Excel files
nuccio_analysis = pd.read_excel('/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.04/mbo001141769st1.adding_isangi.xlsx')
nuccio_isangi_lookup = pd.read_excel('/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.04/reciprocal_best_hits.lookup.xlsx')
dbs_results = pd.read_excel('/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.01/typmu_vs_isangi/results.dbs.xlsx')

# Function to extract UniProtKB ID from the 'Cross-reference' column
def extract_uniprot_id(cross_ref):
    if pd.isnull(cross_ref):
        return None
    parts = cross_ref.split('|')
    for part in parts:
        if part.startswith('UniProtKB:'):
            return part.replace('UniProtKB:', '')
    return None

# Apply the function to create a new column 'UniProtKB_ID'
nuccio_analysis['UniProtKB_ID'] = nuccio_analysis['Cross-reference'].apply(extract_uniprot_id)

# Join nuccio_analysis with nuccio_isangi_lookup based on 'UniProtKB_ID' and 'protein_id'
merged_df = pd.merge(nuccio_analysis, nuccio_isangi_lookup, left_on='UniProtKB_ID', right_on='protein_id', how='left')

# Join the result with dbs_results based on 'gene_2' and 'qseqid_isang'
final_df = pd.merge(merged_df, dbs_results, left_on='qseqid_isang', right_on='gene_2', how='left')

# Optionally, save the final dataframe to a new Excel file
final_df.to_excel('/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.04/2024.11.04.combined_nuccio_and_isangi_dbs.xlsx', index=False)

# Display the first few rows of the final dataframe
print(final_df.head())
