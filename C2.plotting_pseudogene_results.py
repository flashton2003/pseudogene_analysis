import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

# Function to convert percentage strings to floats
def convert_percentage(x):
    if isinstance(x, str) and '%' in x:
        return float(x.strip('%'))
    return x

# Create EI/GI lookup dictionary
ei_gi_lookup = {
    'GCF_000007545.1': 'EI',
    'GCF_000008105.1': 'EI',
    'GCF_000009505.1': 'GI',
    'GCF_000009525.1': 'EI',
    'GCF_000011885.1': 'EI',
    'GCF_000018385.1': 'EI',
    'GCF_000018705.1': 'GI',
    'GCF_000020705.1': 'GI',
    'GCF_000020745.1': 'GI',
    'GCF_000020885.1': 'GI',
    'GCF_000020925.1': 'EI',
    'GCF_000195995.1': 'EI'
}

# Function to extract GCF accession from filename
def extract_gcf(filename):
    match = re.search(r'(GCF_\d{9}\.\d)', filename)
    return match.group(1) if match else None

# Read all files
base_dir = "2024.11.08"
deltaBS_pattern = "2024.11.08b/*/*performance*"
files = glob.glob(f"{base_dir}/*.tsv")
deltaBS_files = glob.glob(deltaBS_pattern)

# Initialize dictionary to store dataframes
data_dict = {
    'bakta': {'EI': [], 'GI': []},
    'bakta_db': {'EI': [], 'GI': []},
    'ncbi': {'EI': [], 'GI': []},
    'salmonella': {'EI': [], 'GI': []},
    'deltaBS': {'EI': [], 'GI': []}
}

# Process regular files
for file in files:
    # Read the TSV file
    df = pd.read_csv(file, sep='\t')
    
    # Convert percentage columns to float
    df['sensitivity'] = df['sensitivity'].apply(convert_percentage)
    df['ppv'] = df['ppv'].apply(convert_percentage)
    
    # Extract GCF accession and get EI/GI classification
    gcf = extract_gcf(file)
    classification = ei_gi_lookup.get(gcf, 'Unknown')
    
    if classification != 'Unknown':
        # Determine which category the file belongs to
        if 'bakta.gff3' in file:
            data_dict['bakta'][classification].append(df)
        elif 'bakta_db_pseudos' in file:
            data_dict['bakta_db'][classification].append(df)
        elif 'ncbi_pseudos' in file:
            data_dict['ncbi'][classification].append(df)
        elif 'salmonella_pseudos' in file:
            data_dict['salmonella'][classification].append(df)

# Process deltaBS files
for file in deltaBS_files:
    # Read the TSV file
    df = pd.read_csv(file, sep='\t')
    
    # Extract GCF accession and get EI/GI classification
    gcf = extract_gcf(file)
    classification = ei_gi_lookup.get(gcf, 'Unknown')
    
    if classification != 'Unknown':
        # Get the second to last row
        row = df.iloc[-2]
        
        # Create a new dataframe with just this row's sensitivity and PPV
        deltaBS_df = pd.DataFrame({
            'input_file': [file],
            'sensitivity': [row['Sensitivity'] * 100],  # Convert to percentage
            'ppv': [row['PPV'] * 100]  # Convert to percentage
        })
        
        data_dict['deltaBS'][classification].append(deltaBS_df)

# Combine dataframes within each category and classification
for category in data_dict:
    for classification in ['EI', 'GI']:
        if data_dict[category][classification]:
            data_dict[category][classification] = pd.concat(
                data_dict[category][classification], 
                ignore_index=True
            )

# Set the font sizes
plt.rc('font', size=14)          # controls default text size
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the x tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the y tick labels
plt.rc('legend', fontsize=12)    # fontsize of the legend
plt.rc('figure', titlesize=18)   # fontsize of the figure title

# Create scatter plot
plt.figure(figsize=(12, 8))
colors = ['blue', 'red', 'green', 'purple', 'orange']
markers = {'EI': 'o', 'GI': 's'}  # Circle for EI, Square for GI

for idx, (category, color) in enumerate(zip(data_dict.keys(), colors)):
    for classification in ['EI', 'GI']:
        if isinstance(data_dict[category][classification], pd.DataFrame) and \
           not data_dict[category][classification].empty:
            plt.scatter(
                data_dict[category][classification]['sensitivity'],
                data_dict[category][classification]['ppv'],
                label=f'{category} ({classification})',
                alpha=0.6,
                c=color,
                marker=markers[classification],
                s=100
            )

plt.xlabel('Sensitivity (%)', fontsize=14, labelpad=10)
plt.ylabel('PPV (%)', fontsize=14, labelpad=10)
plt.title('Sensitivity vs PPV by Category and Classification', fontsize=16, pad=20)

# Increase legend text size and adjust position
plt.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')

# Add grid with increased line width
plt.grid(True, linestyle='--', alpha=0.7, linewidth=0.8)

# Add a bit of padding to the axes
plt.margins(0.1)

# Adjust layout to prevent legend cutoff
plt.tight_layout()

# Save the plot with increased figure size and DPI
plt.savefig('pseudo_analysis_scatter.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSummary Statistics:")
print("-" * 50)
for category in data_dict:
    for classification in ['EI', 'GI']:
        if isinstance(data_dict[category][classification], pd.DataFrame) and \
           not data_dict[category][classification].empty:
            print(f"\n{category.upper()} - {classification}:")
            stats = data_dict[category][classification][['sensitivity', 'ppv']].describe()
            print(stats)