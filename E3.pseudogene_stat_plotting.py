import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Read data
df = pd.read_csv('2024.11.14/2024.11.14.pseudogene_validation_results.diamond_matching.csv')
# df = pd.read_csv('2024.11.14b/2024.11.14.pseudogene_validation_results.coordinate_matching.csv')

# Define markers
salm_markers = {'EI': 'o', 'GI': 's'}

# Create and save PPV vs Sensitivity plot
plt.figure(figsize=(10, 6))
tool_colors = {
    'bakta': '#1f77b4',
    'pseudofinder_baktadb': '#ff7f0e',
    'pseudofinder_salmonella': '#2ca02c',
    'pseudofinder_ncbi': '#d62728',
    'dbs': '#9467bd'
}

for tool in ['bakta', 'pseudofinder_baktadb', 'pseudofinder_salmonella', 'pseudofinder_ncbi', 'dbs']:
    ppv_col = f'{tool}_pseudogene_ppv'
    sens_col = f'{tool}_pseudogene_sensitivity'
    
    for salm_type in df['salm_type'].unique():
        mask = df['salm_type'] == salm_type
        plt.scatter(
            df[mask][sens_col], 
            df[mask][ppv_col],
            c=[tool_colors[tool]],
            marker=salm_markers[salm_type],
            label=f'{tool} ({salm_type})',
            s=100
        )

plt.xlabel('Sensitivity')
plt.ylabel('PPV')
plt.title('PPV vs Sensitivity by Tool and Strain Type')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('ppv_sensitivity_plot.png', bbox_inches='tight', dpi=300)
plt.close()

# Create and save Truth vs Total Positives correlation plot
plt.figure(figsize=(10, 6))
x1 = df['total_positives_in_truth']
y1 = df['pseudofinder_baktadb_pseudogene_total_positives']
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1, y1)
r_squared1 = r_value1 ** 2

for salm_type in df['salm_type'].unique():
    mask = df['salm_type'] == salm_type
    plt.scatter(
        df[mask]['total_positives_in_truth'],
        df[mask]['pseudofinder_baktadb_pseudogene_total_positives'],
        marker=salm_markers[salm_type],
        label=f'{salm_type}',
        s=100
    )

plt.plot(x1, slope1 * x1 + intercept1, color='red', label=f'y = {slope1:.2f}x + {intercept1:.2f}')
plt.xlabel('Total Positives in Truth')
plt.ylabel('Pseudofinder Baktadb Total Positives')
plt.title(f'Correlation Analysis (r² = {r_squared1:.3f}, p = {p_value1:.3e})')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.xlim(0, max(x1) * 1.1)
plt.ylim(0, max(y1) * 1.1)
plt.tight_layout()
plt.savefig('truth_vs_total_positives.png', bbox_inches='tight', dpi=300)
plt.close()

# Create and save CAM Truth vs CAM Count correlation plot
plt.figure(figsize=(10, 6))
x2 = df['total_positives_in_cam_truth']
y2 = df['pseudofinder_baktadb_pseudogene_cam_count']
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2, y2)
r_squared2 = r_value2 ** 2

for salm_type in df['salm_type'].unique():
    mask = df['salm_type'] == salm_type
    plt.scatter(
        df[mask]['total_positives_in_cam_truth'],
        df[mask]['pseudofinder_baktadb_pseudogene_cam_count'],
        marker=salm_markers[salm_type],
        label=f'{salm_type}',
        s=100
    )

plt.plot(x2, slope2 * x2 + intercept2, color='red', label=f'y = {slope2:.2f}x + {intercept2:.2f}')
plt.xlabel('Total Positives in CAM Truth')
plt.ylabel('Pseudofinder Baktadb CAM Count')
plt.title(f'Correlation Analysis (r² = {r_squared2:.3f}, p = {p_value2:.3e})')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.xlim(0, max(x2) * 1.1)
plt.ylim(0, max(y2) * 1.1)
plt.tight_layout()
plt.savefig('cam_truth_vs_cam_count.png', bbox_inches='tight', dpi=300)
plt.close()