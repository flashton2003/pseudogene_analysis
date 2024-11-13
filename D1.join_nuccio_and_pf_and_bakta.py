import pandas as pd
import re
from pathlib import Path

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
                          'score', 'strand', 'frame', 'attribute'])
    
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
    
    # Extract BIANFB locus tags
    def extract_bianfb(attribute):
        if pd.isna(attribute):
            return None
        # Find the old_locus_tag section
        old_locus_match = re.search(r'old_locus_tag=([^;]+)', str(attribute))
        if not old_locus_match:
            return None
        
        # Split multiple locus tags and find BIANFB ones
        locus_tags = old_locus_match.group(1).split(',')
        bianfb_tags = [tag for tag in locus_tags if 'BIANFB_' in tag]
        
        # Return the first BIANFB tag found, or None if none found
        return bianfb_tags[0] if bianfb_tags else None

    df['bianfb_locus_tag'] = df['attribute'].apply(extract_bianfb)

    df.to_csv('pseudofinder_pseudos.csv')

    return df['bianfb_locus_tag'].dropna().tolist()

# Process reciprocal diamond file
def process_diamond(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df[['qseqid_fwd', 'protein_id']]

# Process anaerobic genes file
def process_anaerobic(file_path):
    return pd.read_excel(file_path)['Reference locus tag(s)'].tolist()

def main():
    # Read all files
    nuccio_df = process_nuccio("/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.05b/mbo001141769st1.adding_isangi.xlsx")
    bakta_pseudos = process_bakta("/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.12/CNS1F3.gff3")
    pseudofinder_pseudos = process_pseudofinder("/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.12/CNS1F3_pseudofinder_pseudos.gff")
    diamond_df = process_diamond("/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.12/CNS1F3_vs_nuccio.reciprocal_diamond.tsv")
    anaerobic_genes = process_anaerobic("/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.05b/mbo001141769st7.central_anaerobic_genes.xlsx")
    
    # write pseudofinder_pseudos to file
    # with open("pseudofinder_pseudos.txt", "w") as f:
    #     for item in pseudofinder_pseudos:
    #         f.write("%s\n" % item)

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
    
    # Count anaerobic metabolism genes that are pseudogenes
    pseudo_anaerobic = merged_df[
        (merged_df['central_anaerobic_metabolism'] == 1) & 
        (merged_df['pseudogene'] == 1)
    ].shape[0]
    
    print(f"Number of central anaerobic metabolism genes that are pseudogenes: {pseudo_anaerobic}")
    
    # Optionally save the results
    merged_df.to_excel("processed_results.xlsx", index=False)

if __name__ == "__main__":
    main()