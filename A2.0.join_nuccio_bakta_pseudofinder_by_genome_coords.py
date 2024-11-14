import pandas as pd
import argparse
from collections import defaultdict
import re
import os

# Mapping dictionaries remain unchanged
STRAIN_MAPPING = {
    'GCF_000020705.1': 'SL476',
    'GCF_000020745.1': 'CVM19633',
    'GCF_000020885.1': 'SL483',
    'GCF_000009505.1': 'P125109',
    'GCF_000018705.1': 'SPB7',
    'GCF_000195995.1': 'CT18',
    'GCF_000007545.1': 'Ty2',
    'GCF_000011885.1': 'ATCC 9150',
    'GCF_000020925.1': 'CT_02021853',
    'GCF_000009525.1': '287/91',
    'GCF_000008105.1': 'SC-B67',
    'GCF_000018385.1': 'RKS4594',
    'GCF_000026565.1': 'AKU_12601'
}

SEQNAME_MAPPING = {
    'GCF_000007545.1': {'contig_1': 'NC_004631.1'},
    'GCF_000008105.1': {'contig_1': 'NC_006905.1', 'contig_2': 'NC_006855.1', 'contig_3': 'NC_006856.1'},
    'GCF_000009505.1': {'contig_1': 'NC_011294.1'},
    'GCF_000009525.1': {'contig_1': 'NC_011274.1'},
    'GCF_000011885.1': {'contig_1': 'NC_006511.1'},
    'GCF_000018385.1': {'contig_1': 'NC_012125.1', 'contig_2': 'NC_012124.1'},
    'GCF_000018705.1': {'contig_1': 'NC_010102.1'},
    'GCF_000020705.1': {'contig_1': 'NC_011083.1', 'contig_2': 'NC_011082.1', 'contig_3': 'NC_011081.1'},
    'GCF_000020745.1': {'contig_1': 'NC_011094.1', 'contig_2': 'NC_011092.1', 'contig_3': 'NC_011093.1'},
    'GCF_000020885.1': {'contig_1': 'NC_011149.1', 'contig_2': 'NC_011148.1'},
    'GCF_000020925.1': {'contig_1': 'NC_011205.1', 'contig_2': 'NC_011204.1'},
    'GCF_000026565.1': {'contig_1': 'NC_011147.1'},
    'GCF_000195995.1': {'contig_1': 'NC_003198.1', 'contig_2': 'NC_003384.1', 'contig_3': 'NC_003385.1'}
}

def sanitize_coordinate(value: str) -> float:
    """
    Clean coordinate values by removing non-numeric characters (except minus sign and decimal point)
    and converting to a float.
    """
    if pd.isna(value):
        return None
    # Allow only digits, minus sign, and decimal point
    cleaned = re.sub(r'[^\d.-]', '', str(value))
    return float(cleaned)

def sanitize_seqname(value: str) -> str:
    """
    Clean seqname by removing unwanted characters and standardizing format
    """
    if pd.isna(value):
        return ''
    # Remove '>' and any leading/trailing whitespace
    cleaned = str(value).strip().lstrip('>')
    return cleaned

def parse_gff(gff_file, accession):
    """Parse GFF file and convert seqnames using mapping"""
    calls = []
    seqname_map = SEQNAME_MAPPING[accession]
    is_bakta = gff_file.lower().endswith('bakta.gff3')
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            # For Bakta files, only process records marked as pseudogenes
            if is_bakta and not "pseudo=True" in fields[8]:
                continue
                
            # Sanitize the seqname before mapping
            orig_seqname = sanitize_seqname(fields[0])
            if orig_seqname in seqname_map:
                seqname = seqname_map[orig_seqname]
                # Sanitize coordinates
                start = sanitize_coordinate(fields[3])
                end = sanitize_coordinate(fields[4])
                
                # Only add valid coordinates
                if start > 0 and end > 0:
                    calls.append({
                        'seqname': seqname,
                        'start': int(start),  # Convert to int after validation
                        'end': int(end),      # Convert to int after validation
                        'strand': fields[6],
                        'attributes': fields[8]
                    })
    
    return calls

def parse_strain_column(strain_data):
    """Parse the strain column format to extract coordinates"""
    parts = strain_data.split('|')
    if len(parts) >= 5:
        # Sanitize coordinates from strain data
        return {
            'seqname': sanitize_seqname(parts[2]),
            'start': int(sanitize_coordinate(parts[3])),
            'end': int(sanitize_coordinate(parts[4]))
        }
    return None

def check_overlap(call, truth_coords):
    """Check if two regions overlap"""
    return (call['seqname'] == truth_coords['seqname'] and 
            call['start'] <= truth_coords['end'] and 
            call['end'] >= truth_coords['start'])

def main():
    parser = argparse.ArgumentParser(description='Match GFF calls to truth dataset')
    parser.add_argument('--call', help='Input GFF file')
    parser.add_argument('--truth', help='Truth dataset file')
    parser.add_argument('--accession', help='Genome accession (e.g. GCF_000007545.1)')
    parser.add_argument('--output', help='Output file name')
    
    args = parser.parse_args()
    
    # Get strain name from accession
    if args.accession not in STRAIN_MAPPING:
        raise ValueError(f"Unknown accession: {args.accession}")
    strain = STRAIN_MAPPING[args.accession]
    
    # Parse input files
    calls = parse_gff(args.call, args.accession)
    
    # Read truth dataset
    truth = pd.read_excel(args.truth)
    
    # Add is_pseudogene column if it doesn't exist
    if 'is_pseudogene' not in truth.columns:
        truth['is_pseudogene'] = 0
    
    # Process each row in the truth dataset
    for idx, row in truth.iterrows():
        strain_data = row[strain]
        if pd.notna(strain_data) and strain_data != '3|Absent':
            coords = parse_strain_column(strain_data)
            if coords:
                # Check for overlap with any call
                for call in calls:
                    if check_overlap(call, coords):
                        if row['Reference locus tag(s)'] == 'STM0212':
                            print(call)
                            print(strain_data)
                        truth.at[idx, 'is_pseudogene'] = 1
                        break
    
    # Write modified truth dataset
    truth.to_excel(args.output, index=False)

if __name__ == "__main__":
    main()