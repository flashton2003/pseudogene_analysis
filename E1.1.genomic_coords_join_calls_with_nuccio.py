import pandas as pd
import re
import argparse
import pprint

# Mapping dictionary for strain lookups
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

def convert_seqname(seqname, gcf):
    """Convert contig names to NC accessions based on GCF"""
    if gcf in SEQNAME_MAPPING:
        return SEQNAME_MAPPING[gcf].get(seqname, seqname)
    return seqname

def sanitize_coordinate(value: str) -> float:
    """
    Clean coordinate values by removing non-numeric characters (except minus sign and decimal point)
    and converting to a float.
    """
    if pd.isna(value):
        return None
    # Allow only digits, minus sign, and decimal point
    cleaned = re.sub(r'[^\d.-]', '', str(value))
    try:
        return float(cleaned)
    except ValueError:
        return None

def sanitize_seqname(value: str) -> str:
    """
    Clean seqname by removing unwanted characters and standardizing format
    """
    if pd.isna(value):
        return ''
    # Remove '>' and any leading/trailing whitespace
    cleaned = str(value).strip().lstrip('>')
    return cleaned

class GenomicRegion:
    def __init__(self, seqname, start, end, gcf=None):
        self.seqname = convert_seqname(sanitize_seqname(seqname), gcf) if gcf else sanitize_seqname(seqname)
        start_coord = sanitize_coordinate(start)
        end_coord = sanitize_coordinate(end)
        
        if start_coord is None or end_coord is None:
            raise ValueError(f"Invalid coordinates: start={start}, end={end}")
            
        # Ensure start is less than end
        self.start = int(min(start_coord, end_coord))
        self.end = int(max(start_coord, end_coord))

def check_overlap(region1, region2):
    """Check if two genomic regions overlap"""
    if region1.seqname != region2.seqname:
        return False
    return (region1.start <= region2.end) and (region2.start <= region1.end)

def parse_strain_column(strain_data):
    """Parse the strain column to extract genomic coordinates"""
    try:
        parts = strain_data.split('|')
        return GenomicRegion(parts[2], parts[3], parts[4])
    except (IndexError, ValueError):
        return None

def process_nuccio(file_path):
    """Read and process the nuccio file to get truth data"""
    df = pd.read_excel(file_path)
    return df

def extract_gene_id(attribute_string):
    """Extract gene ID from GFF attribute string"""
    match = re.search(r'ID=([^;]+)', attribute_string)
    if match:
        return match.group(1)
    return None

def extract_gene_id(attribute_string):
    """Extract gene ID from GFF attribute string"""
    match = re.search(r'ID=([^;]+)', attribute_string)
    if match:
        return match.group(1)
    return None

def process_bakta_gff(file_path, gcf=None):
    """
    Process a Bakta GFF file and return:
    1. List of pseudogene regions 
    2. Dictionary mapping gene IDs to GenomicRegion objects for all CDS/gene features
    
    Args:
        file_path: Path to Bakta GFF file
        gcf: GCF accession for coordinate conversion
        
    Returns:
        tuple: (list of pseudogene GenomicRegions, dict of gene_id:GenomicRegion)
    """
    if not file_path:
        return [], {}
        
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                     names=['seqname', 'source', 'feature', 'start', 'end', 
                           'score', 'strand', 'frame', 'attribute'],
                     low_memory=False)
    
    # Create coordinate mapping dictionary for all CDS/gene features
    coord_dict = {}
    pseudo_regions = []
    
    for _, row in df.iterrows():
        # Process all CDS/gene features for coordinate dictionary
        if row['feature'] == 'CDS' or row['feature'] == 'gene':
            gene_id = extract_gene_id(row['attribute'])
            if gene_id:
                try:
                    coord_dict[gene_id] = GenomicRegion(row['seqname'], row['start'], row['end'], gcf)
                except ValueError:
                    print(f"Warning: Invalid coordinates for gene {gene_id}")
                    continue
        
        # Track pseudogene regions
        if 'pseudo=True' in str(row['attribute']):
            try:
                pseudo_regions.append(GenomicRegion(row['seqname'], row['start'], row['end'], gcf))
            except ValueError as e:
                print(f"Warning: Skipping invalid pseudogene region: {e}")
                continue
    
    return pseudo_regions, coord_dict

def process_pseudofinder_gff(file_path, gcf=None):
    """
    Process a Pseudofinder GFF file and return list of pseudogene regions.
    Only returns genomic regions marked as pseudogenes, no gene ID mapping needed.
    
    Args:
        file_path: Path to Pseudofinder GFF file
        gcf: GCF accession for coordinate conversion
        
    Returns:
        list: List of GenomicRegion objects for pseudogenes
    """
    if not file_path:
        return []
        
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None,
                     names=['seqname', 'source', 'feature', 'start', 'end', 
                           'score', 'strand', 'frame', 'attribute'],
                     low_memory=False)
    
    pseudo_regions = []
    
    for _, row in df.iterrows():
        try:
            pseudo_regions.append(GenomicRegion(row['seqname'], row['start'], row['end'], gcf))
        except ValueError as e:
            print(f"Warning: Skipping invalid Pseudofinder region: {e}")
            continue
    
    return pseudo_regions

def process_dbs(file_path, coord_dict):
    """Process DBS results file and extract coordinates using gene ID lookup"""
    if not file_path or not coord_dict:
        return []
    
    regions = []
    df = pd.read_csv(file_path, sep='\t', skiprows=1)
    
    # Calculate 97.5th percentile threshold of delta-bitscore
    threshold = df['delta-bitscore'].quantile(0.975)
    
    for _, row in df.iterrows():
        if pd.notna(row['gene_2']) and row['delta-bitscore'] > threshold:
            gene_id = row['gene_2'].split('|')[0] if '|' in row['gene_2'] else row['gene_2']
            
            if gene_id in coord_dict:
                try:
                    regions.append(coord_dict[gene_id])
                except (ValueError, AttributeError) as e:
                    print(f"Warning: Skipping invalid DBS region for gene {gene_id}: {e}")
                    continue
            else:
                print(f"Warning: Gene ID {gene_id} not found in coordinate dictionary")
    
    return regions

def calculate_overlaps(truth_df, strain_column, call_sets):
    """Calculate overlaps between truth data and call sets"""
    # Initialize pseudogene columns with zeros
    for call_set_name in call_sets.keys():
        truth_df[f'{call_set_name}_pseudogene'] = 0
    
    for idx, row in truth_df.iterrows():
        strain_data = row[strain_column]
        if pd.notna(strain_data) and strain_data != '3|Absent':
            truth_coords = parse_strain_column(strain_data)
            if truth_coords:
                # Check each call set
                for call_set_name, calls in call_sets.items():
                    for call in calls:
                        if check_overlap(call, truth_coords):
                            truth_df.at[idx, f'{call_set_name}_pseudogene'] = 1
                            break
    
    return truth_df

def create_consensus_row(group):
    """Create a consensus row from a group of rows based on pseudogene columns"""
    consensus = group.iloc[0].copy()
    
    pseudo_columns = [col for col in group.columns if col.endswith('_pseudogene')]
    
    for col in pseudo_columns:
        if col in group.columns:
            consensus[col] = 1 if (group[col] == 1).any() else 0
    
    return pd.Series(consensus)

def process_anaerobic(file_path):
    """Process anaerobic genes file"""
    return pd.read_excel(file_path)['Reference locus tag(s)'].tolist()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pseudogene data using coordinate overlap')
    parser.add_argument('--nuccio', required=True, help='Path to Nuccio Excel file')
    parser.add_argument('--gcf', required=True, help='GCF accession for strain lookup')
    parser.add_argument('--bakta', help='Path to Bakta GFF3 file')
    parser.add_argument('--pseudofinder-baktadb', help='Path to Pseudofinder GFF file (BaktaDB)')
    parser.add_argument('--pseudofinder-salmonella', help='Path to Pseudofinder GFF file (Salmonella)')
    parser.add_argument('--pseudofinder-ncbi', help='Path to Pseudofinder GFF file (NCBI)')
    parser.add_argument('--dbs', help='Path to DBS TSV file')
    parser.add_argument('--anaerobic', required=True, help='Path to anaerobic genes Excel file')
    parser.add_argument('--output', required=True, help='Path for output Excel file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Get strain name from GCF accession
    if args.gcf not in STRAIN_MAPPING:
        raise ValueError(f"Unknown GCF accession: {args.gcf}")
    strain = STRAIN_MAPPING[args.gcf]
    
    # Process truth data and anaerobic
    truth_df = process_nuccio(args.nuccio)
    anaerobic_genes = process_anaerobic(args.anaerobic)

    # Add anaerobic metabolism column
    truth_df['central_anaerobic_metabolism'] = truth_df['Reference locus tag(s)'].isin(anaerobic_genes).astype(int)
    

    # Process Bakta GFF file to get both pseudogenes and coordinate dictionary
    bakta_regions, coord_dict = process_bakta_gff(args.bakta, gcf=args.gcf)

    # Process all call sets and store their coordinates
    call_sets = {
        'bakta': bakta_regions,
        'pseudofinder_baktadb': process_pseudofinder_gff(args.pseudofinder_baktadb, gcf=args.gcf),
        'pseudofinder_salmonella': process_pseudofinder_gff(args.pseudofinder_salmonella, gcf=args.gcf),
        'pseudofinder_ncbi': process_pseudofinder_gff(args.pseudofinder_ncbi, gcf=args.gcf),
        'dbs': process_dbs(args.dbs, coord_dict)
    }
    
    # Remove empty call sets
    call_sets = {k: v for k, v in call_sets.items() if v}
    
    # Calculate overlaps and update truth dataframe
    results_df = calculate_overlaps(truth_df, strain, call_sets)
    
    # Write output to Excel file
    results_df.to_excel(args.output, index=False)
    
    print(f"Processing complete. Output saved to: {args.output}")
    print(f"Processed strain: {strain}")
    print("Sheets created: Complete_Data and Deduplicated_Data")

if __name__ == "__main__":
    main()