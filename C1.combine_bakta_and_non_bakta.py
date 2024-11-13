import pandas as pd
import sys
import numpy as np
import re
import os

def create_pseudogene_column(df):
    """Create pseudogene column based on delta-bitscore threshold"""
    # Calculate 97.5th percentile of delta-bitscore
    threshold = np.percentile(df['delta-bitscore'].dropna(), 97.5)
    # print(threshold)
    # Create pseudogene column where 1 indicates above threshold
    df['is_pseudogene'] = (df['delta-bitscore'] > threshold).astype(int)
    return df

def process_files(bakta_file, nonbakta_file):
    """Process and combine pseudogene information from both files"""
    # Read Excel files
    bakta_df = pd.read_excel(bakta_file)
    
    # Check if nonbakta file needs pseudogene column creation
    if nonbakta_file.endswith('dbs.xlsx'):
        nonbakta_df = pd.read_excel(nonbakta_file, sheet_name='Deduplicated_Data')
        nonbakta_df = create_pseudogene_column(nonbakta_df)
    else:
        nonbakta_df = pd.read_excel(nonbakta_file)
    # print(nonbakta_df['delta-bitscore'], nonbakta_df[['is_pseudogene']])
    # nonbakta_df.to_csv('nonbakta_df.csv')
    
    # Ensure both dataframes have is_pseudogene column
    if 'is_pseudogene' not in bakta_df.columns or 'is_pseudogene' not in nonbakta_df.columns:
        raise ValueError("Missing is_pseudogene column in one or both files")
    
    # Combine pseudogene information
    # select only 'Index' and 'is_pseudogene' column from nonbakta_df
    combined_df = pd.merge(bakta_df, nonbakta_df[['Index', 'is_pseudogene']], 
                          left_on='Index', right_on='Index', 
                          suffixes=('_bakta', '_nonbakta'))
    
    # Create combined pseudogene column (1 if either source has 1)
    combined_df['is_pseudogene'] = ((combined_df['is_pseudogene_bakta'] == 1) | 
                                           (combined_df['is_pseudogene_nonbakta'] == 1)).astype(int)
    
    return combined_df


def extract_gcf(filename):
        match = re.search(r'(GCF_\d{9}\.\d)', filename)
        return match.group(1) if match else None    

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <bakta_file.xlsx> <nonbakta_file.xlsx> <output_root_dir>")
        sys.exit(1)
    
    bakta_file = sys.argv[1]
    nonbakta_file = sys.argv[2]
    output_root_dir = sys.argv[3]

    # i want to extract the accession, which looks like this "GCF_000007545.1" from the bakta_file 
    # which looks like this 2024.11.08/GCF_000007545.1.bakta.gff3.vs_nuccio.xlsx
    accession = extract_gcf(bakta_file)
    # nonbakta filename looks like this 2024.11.08b/2024.11.08.combined_nuccio_and_GCF_000007545_dbs.xlsx
    # i want something like this 2024.11.08.combined_nuccio_and_GCF_000007545_dbs
    nonbakta_filename_clean = os.path.splitext(os.path.basename(nonbakta_file))[0]

    
    # try:
    result_df = process_files(bakta_file, nonbakta_file)
    
    # Save to new Excel file
    output_filename = f'{output_root_dir}/{accession}_bakta_{nonbakta_filename_clean}_combined.xlsx'
    result_df.to_excel(output_filename, index=True)
    print(f"Combined results saved to {output_filename}")
        
    # except Exception as e:
    #     print(f"Error processing files: {str(e)}")
    #     sys.exit(1)

if __name__ == "__main__":
    main()