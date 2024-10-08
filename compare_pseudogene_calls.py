import pandas as pd

# Function to check overlap by 50% rule
def is_overlap(truth_start, truth_end, call_start, call_end):
    # Get the length of the truth pseudogene
    truth_length = truth_end - truth_start
    # Determine the overlap region
    # This line ensures that the starting point of the overlap is the later of the two starts 
    # (either the truth pseudogene's start or the call pseudogene's start).
    # This is because an overlap can only begin after both genes have started.
    overlap_start = max(truth_start, call_start)
    # This line ensures that the end point of the overlap is the earlier of the two ends 
    # (either the truth pseudogene's end or the call pseudogene's end). 
    # An overlap can only extend up to the point where both genes have not yet ended.
    overlap_end = min(truth_end, call_end)

    if overlap_start < overlap_end:  # valid overlap exists
        overlap_length = overlap_end - overlap_start
        # print('overlap_length:', overlap_length)
        # Check if the overlap is at least 50% of the truth gene length
        if overlap_length >= 0.1 * truth_length:
            return True
    return False

# Function to calculate sensitivity and specificity
def calculate_sensitivity_ppv(truth_df, call_df):
    true_positives = 0
    false_negatives = 0
    false_positives = 0
    
    # Loop over truth and call data
    for _, truth_row in truth_df.iterrows():
        truth_start = int(truth_row['start'])
        truth_end = int(truth_row['stop'])
        
        # Check if there is any matching pseudogene in the call data
        matched = False
        for _, call_row in call_df.iterrows():
            call_start = int(call_row['start'])
            call_end = int(call_row['end'])
            
            if is_overlap(truth_start, truth_end, call_start, call_end):
                matched = True
                # print('Matched')
                # print(truth_start, truth_end, call_start, call_end)
                break
        
        if matched:
            true_positives += 1
        else:
            false_negatives += 1
    
    # For ppv, check calls that do not match the truth
    for _, call_row in call_df.iterrows():
        call_start = int(call_row['start'])
        call_end = int(call_row['end'])
        
        matched = False
        for _, truth_row in truth_df.iterrows():
            truth_start = int(truth_row['start'])
            truth_end = int(truth_row['stop'])
            
            if is_overlap(truth_start, truth_end, call_start, call_end):
                matched = True
                break
        
        if not matched:
            false_positives += 1
    
    # Sensitivity = True Positives / (True Positives + False Negatives)
    sensitivity = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    # do positive predictive value
    ppv = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    
    return sensitivity, ppv


def read_and_filter_data(file_path, file_type):
    if file_type == 'excel':
        df = pd.read_excel(file_path)
    elif file_type == 'csv':
        df = pd.read_csv(file_path)
    elif file_type == 'gff':
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
        df = df[df['seqname'] == 'contig_1']
    else:
        raise ValueError(f"Unsupported file type: {file_type}")
    return df

# Load the provided files
# truth_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/pbio.3000059.s002.pseudogenes.csv"
truth_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/nuccio_pseudo.NC_003197.xlsx"
pseudofinder_baktadb_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/GCF_000006945.2_bakta_db_pseudos.gff"
pseudofinder_singlegenomedb_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/GCF_000006945.2_ncbi_pseudos.gff"
pseudofinder_salmonellapangenomedb_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/GCF_000006945.2_salmonella_pseudos.gff"
bakta_annotation = '/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.10.07/GCF_000006945.2.gff3'


truth_df = read_and_filter_data(truth_file, 'excel')
pseudofinder_baktadb_df = read_and_filter_data(pseudofinder_baktadb_file, 'gff')
pseudofinder_singlegenomedb_df = read_and_filter_data(pseudofinder_singlegenomedb_file, 'gff')
pseudofinder_salmonellapangenomedb_df = read_and_filter_data(pseudofinder_salmonellapangenomedb_file, 'gff')
bakta_annotation = read_and_filter_data(bakta_annotation, 'gff')

bakta_pseudogene_df = bakta_annotation[bakta_annotation['attribute'].notna() & bakta_annotation['attribute'].str.contains('pseudo=True')]

# Calculating sensitivity and specificity for each pseudogene file
pseudofinder_baktadb_sensitivity, pseudofinder_baktadb_ppv = calculate_sensitivity_ppv(truth_df, pseudofinder_baktadb_df)
pseudofinder_singlegenomedb_sensitivity, pseudofinder_singlegenomedb_ppv = calculate_sensitivity_ppv(truth_df, pseudofinder_singlegenomedb_df)
pseudofinder_salmonellapangenomedb_sensitivity, pseudofinder_salmonellapangenomedb_ppv = calculate_sensitivity_ppv(truth_df, pseudofinder_salmonellapangenomedb_df)
bakta_pseudogene_sensitivity, bakta_pseudogene_ppv = calculate_sensitivity_ppv(truth_df, bakta_pseudogene_df)

# Prepare the results for display
results = {
    "Dataset": ["pseudofinder_baktadb", "pseudofinder_singlegenomedb", "pseudofinder_salmonellapangenomedb", "Bakta_pseudogene"],
    "Sensitivity": [pseudofinder_baktadb_sensitivity, pseudofinder_singlegenomedb_sensitivity, pseudofinder_salmonellapangenomedb_sensitivity, bakta_pseudogene_sensitivity],
    "PPV": [pseudofinder_baktadb_ppv, pseudofinder_singlegenomedb_ppv, pseudofinder_salmonellapangenomedb_ppv, bakta_pseudogene_ppv]
}

# Convert the results to a DataFrame for better presentation
results_df = pd.DataFrame(results)
print(results_df)
