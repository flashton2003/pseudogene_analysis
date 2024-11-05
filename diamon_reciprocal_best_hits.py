import pandas as pd

# Define file names
isang_vs_nucciocogs = '2024.11.04/CNS1F3_vs_2024.11.04.nuccio_baumler_uniprotkb.diamond.tsv'
nucciocogs_vs_isangi = '2024.11.04/2024.11.04.nuccio_baumler_uniprotkb_vs_CNS1F3.diamond.tsv'

# Load the files into DataFrames
df_isangi_vs_nucciocogs = pd.read_csv(isang_vs_nucciocogs, sep='\t')
df_nucciocogs_vs_isangi = pd.read_csv(nucciocogs_vs_isangi, sep='\t')

# Extract protein ID from sseqid by splitting on '|' and taking the second element ([1])
df_isangi_vs_nucciocogs['protein_id'] = df_isangi_vs_nucciocogs['sseqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)
df_nucciocogs_vs_isangi['protein_id'] = df_nucciocogs_vs_isangi['qseqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)

# Calculate query coverage and subject coverage
df_isangi_vs_nucciocogs['query_cov'] = (df_isangi_vs_nucciocogs['length'] - df_isangi_vs_nucciocogs['gaps']) / df_isangi_vs_nucciocogs['qlen'] * 100
df_isangi_vs_nucciocogs['subject_cov'] = (df_isangi_vs_nucciocogs['length'] - df_isangi_vs_nucciocogs['gaps']) / df_isangi_vs_nucciocogs['slen'] * 100

df_nucciocogs_vs_isangi['query_cov'] =(df_nucciocogs_vs_isangi['length'] - df_nucciocogs_vs_isangi['gaps']) / df_nucciocogs_vs_isangi['qlen'] * 100
df_nucciocogs_vs_isangi['subject_cov'] = (df_nucciocogs_vs_isangi['length'] - df_nucciocogs_vs_isangi['gaps']) / df_nucciocogs_vs_isangi['slen'] * 100

# Filter based on query coverage, subject coverage, and evalue
filtered_isangi_vs_nucciocogs = df_isangi_vs_nucciocogs[
    (df_isangi_vs_nucciocogs['query_cov'] > 70) &
    (df_isangi_vs_nucciocogs['subject_cov'] > 70) &
    (df_isangi_vs_nucciocogs['evalue'] < 1e-15)
]

filtered_nucciocogs_vs_isangi = df_nucciocogs_vs_isangi[
    (df_nucciocogs_vs_isangi['query_cov'] > 70) &
    (df_nucciocogs_vs_isangi['subject_cov'] > 70) &
    (df_nucciocogs_vs_isangi['evalue'] < 1e-15)
]

filtered_isangi_vs_nucciocogs.to_csv('~/Desktop/filtered_isangi_vs_nucciocogs.csv', index=False)

# Identify the best hits (highest bitscore) for each query in each filtered DataFrame
best_hits_isangi_vs_nucciocogs = filtered_isangi_vs_nucciocogs.groupby('qseqid').apply(lambda x: x.loc[x['bitscore'].idxmax()])
best_hits_nucciocogs_vs_isangi = filtered_nucciocogs_vs_isangi.groupby('qseqid').apply(lambda x: x.loc[x['bitscore'].idxmax()])

best_hits_isangi_vs_nucciocogs.to_csv('~/Desktop/best_hits_isangi_vs_nucciocogs.csv', index=False)

print(best_hits_isangi_vs_nucciocogs)

print(best_hits_nucciocogs_vs_isangi)

# Reset index to flatten the grouped DataFrames
best_hits_isangi_vs_nucciocogs = best_hits_isangi_vs_nucciocogs.reset_index(drop=True)
best_hits_nucciocogs_vs_isangi = best_hits_nucciocogs_vs_isangi.reset_index(drop=True)

# Merge to find reciprocal best hits
# We look for matches where qseqid in file1 matches sseqid in file2 and vice versa
reciprocal_hits = best_hits_isangi_vs_nucciocogs.merge(
    best_hits_nucciocogs_vs_isangi,
    left_on=['qseqid', 'protein_id'],
    right_on=['sseqid', 'protein_id'],
    suffixes=('_isang', '_nuccio')
)

# Display the reciprocal best hits
print("Reciprocal Best Hits:")
# print(reciprocal_hits[['qseqid', 'sseqid_isang', 'sseqid_nuccio', 'protein_id', 'bitscore_isang', 'bitscore_nuccio', 'pident_isang', 'query_cov_isang', 'subject_cov_isang', 'pident_nuccio', 'query_cov_nuccio', 'subject_cov_nuccio']])
print(reciprocal_hits)
reciprocal_hits.to_csv('2024.11.04/reciprocal_best_hits.tsv', sep='\t', index=False)