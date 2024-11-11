import pandas as pd
import matplotlib.pyplot as plt
import glob
import pprint

# Function to convert percentage strings to floats
def convert_percentage(x):
    if isinstance(x, str) and '%' in x:
        return float(x.strip('%'))
    return x

# Read all files
base_dir = "2024.11.08"
deltaBS_pattern = "2024.11.08b/*/*performance*"
files = glob.glob(f"{base_dir}/*.tsv")
deltaBS_files = glob.glob(deltaBS_pattern)

# Initialize dictionary to store dataframes
data_dict = {
    'bakta': [],
    'bakta_db': [],
    'ncbi': [],
    'salmonella': [],
    'deltaBS': []  # Added deltaBS category
}

# Process regular files
for file in files:
    # Read the TSV file
    df = pd.read_csv(file, sep='\t')
    
    # Convert percentage columns to float
    df['sensitivity'] = df['sensitivity'].apply(convert_percentage)
    df['ppv'] = df['ppv'].apply(convert_percentage)
    
    # Determine which category the file belongs to
    if 'bakta.gff3' in file:
        data_dict['bakta'].append(df)
    elif 'bakta_db_pseudos' in file:
        data_dict['bakta_db'].append(df)
    elif 'ncbi_pseudos' in file:
        data_dict['ncbi'].append(df)
    elif 'salmonella_pseudos' in file:
        data_dict['salmonella'].append(df)

# Process deltaBS files
for file in deltaBS_files:
    # Read the TSV file
    df = pd.read_csv(file, sep='\t')
    
    # Get the second to last row
    row = df.iloc[-2]
    
    # Create a new dataframe with just this row's sensitivity and PPV
    deltaBS_df = pd.DataFrame({
        'input_file': [file],
        'sensitivity': [row['Sensitivity'] * 100],  # Convert to percentage
        'ppv': [row['PPV'] * 100]  # Convert to percentage
    })
    
    data_dict['deltaBS'].append(deltaBS_df)



# Combine dataframes within each category
for category in data_dict:
    if data_dict[category]:
        data_dict[category] = pd.concat(data_dict[category], ignore_index=True)

for category in data_dict:
    print(category)
    print(data_dict[category])

# Create scatter plot
plt.figure(figsize=(12, 8))
colors = ['blue', 'red', 'green', 'purple', 'orange']  # Added color for deltaBS
markers = ['o', 's', '^', 'D', '*']  # Added marker for deltaBS

for (category, color, marker) in zip(data_dict.keys(), colors, markers):
    if not data_dict[category].empty:
        plt.scatter(
            data_dict[category]['sensitivity'],
            data_dict[category]['ppv'],
            label=category,
            alpha=0.6,
            c=color,
            marker=marker,
            s=100
        )

plt.xlabel('Sensitivity (%)')
plt.ylabel('PPV (%)')
plt.title('Sensitivity vs PPV by Category')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

# Add a bit of padding to the axes
plt.margins(0.1)

# Save the plot
plt.savefig('pseudo_analysis_scatter.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSummary Statistics:")
print("-" * 50)
for category in data_dict:
    if not data_dict[category].empty:
        print(f"\n{category.upper()}:")
        stats = data_dict[category][['sensitivity', 'ppv']].describe()
        print(stats)