import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from datetime import datetime

# Path configuration
db_path = 'data/ATLAS2025.db'
base_output_dir = 'postprocessing/results2025/motif_distribution'

# Create output directories
large_dir = os.path.join(base_output_dir, 'large')
medium_dir = os.path.join(base_output_dir, 'medium')
small_dir = os.path.join(base_output_dir, 'small')
data_dir = os.path.join(base_output_dir, 'data')

for dir_path in [large_dir, medium_dir, small_dir, data_dir]:
    os.makedirs(dir_path, exist_ok=True)

# Plot styling configuration - Three versions
# Large version (original)
large_base_font_size = 28 * 1.5
large_bar_label_font_size = 18
large_x_tick_font_size = 30

# Medium version (75% of large)
medium_base_font_size = int(large_base_font_size * 0.75)
medium_bar_label_font_size = int(large_bar_label_font_size * 0.75)
medium_x_tick_font_size = int(large_x_tick_font_size * 0.75)

# Small version
small_base_font_size = 16
small_bar_label_font_size = 12
small_x_tick_font_size = 14

colors = ['#264653', '#287171','#2a9d8f', '#e9c46a', '#f4a261', '#e76f51']
figure_size = (12, 9)

# Connect to database and extract data from both tables
conn = sqlite3.connect(db_path)

# Extract regular motif data
data_query = "SELECT motif_type, nt_number FROM data"
df = pd.read_sql_query(data_query, conn)

# Extract PK data
pk_query = "SELECT motif_type FROM PK"
pk_df = pd.read_sql_query(pk_query, conn)

conn.close()

# Calculate motif size from nucleotide numbers
def calculate_motif_size(nt_numbers):
    return len(nt_numbers.split(','))

df['motif_size'] = df['nt_number'].apply(calculate_motif_size)

# Count motifs by type and add PK count
motif_counts = df['motif_type'].value_counts()
motif_counts['PK'] = len(pk_df)  # Add total PK count

# Calculate size distributions for regular motifs
distribution = df.groupby('motif_type')['motif_size'].value_counts().unstack().fillna(0)

# === NEW: Generate CSV outputs ===
# 1. Save motif type counts
motif_counts_df = pd.DataFrame({
    'motif_type': motif_counts.index,
    'count': motif_counts.values
})
motif_counts_df.to_csv(os.path.join(data_dir, 'motif_type_counts.csv'), index=False)

# 2. Save PK distribution
pk_types = pk_df['motif_type'].str.replace('pseudoknot_', '', regex=False)
pk_counts = pk_types.value_counts()
pk_counts_df = pd.DataFrame({
    'pk_type': pk_counts.index,
    'count': pk_counts.values
})
pk_counts_df.to_csv(os.path.join(data_dir, 'pk_distribution.csv'), index=False)

# 3. Save size distributions for each motif type
for motif_type in distribution.index:
    size_dist = distribution.loc[motif_type]
    size_dist = size_dist[size_dist > 0]  # Only include non-zero counts
    size_dist_df = pd.DataFrame({
        'motif_size': size_dist.index,
        'count': size_dist.values
    })
    size_dist_df.to_csv(os.path.join(data_dir, f'{motif_type}_size_distribution.csv'), index=False)

# === NEW: Generate statistical report ===
report_lines = []
report_lines.append(f"RNA Motif Library Statistical Report")
report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append("=" * 50)
report_lines.append("")

# Total counts
total_regular_motifs = len(df)
total_pk_motifs = len(pk_df)
total_motifs = total_regular_motifs + total_pk_motifs

report_lines.append(f"TOTAL MOTIFS: {total_motifs}")
report_lines.append(f"  Regular motifs: {total_regular_motifs}")
report_lines.append(f"  Pseudoknot motifs: {total_pk_motifs}")
report_lines.append("")

# Regular motif type distribution
report_lines.append("REGULAR MOTIF TYPE DISTRIBUTION:")
for motif_type, count in motif_counts.drop('PK').items():
    percentage = (count / total_regular_motifs) * 100
    report_lines.append(f"  {motif_type}: {count} ({percentage:.1f}%)")
report_lines.append("")

# PK distribution
report_lines.append("PSEUDOKNOT TYPE DISTRIBUTION:")
for pk_type, count in pk_counts.items():
    percentage = (count / total_pk_motifs) * 100
    report_lines.append(f"  {pk_type}: {count} ({percentage:.1f}%)")
report_lines.append("")

# Size statistics for each motif type
report_lines.append("SIZE STATISTICS BY MOTIF TYPE:")
for motif_type in ['internal', 'hairpin', 'bulge', '3_junction', '4_junction', '5_junction']:
    if motif_type in df['motif_type'].values:
        sizes = df[df['motif_type'] == motif_type]['motif_size']
        report_lines.append(f"  {motif_type}:")
        report_lines.append(f"    Count: {len(sizes)}")
        report_lines.append(f"    Size range: {sizes.min()} - {sizes.max()}")
        report_lines.append(f"    Mean size: {sizes.mean():.1f}")
        report_lines.append(f"    Median size: {sizes.median():.1f}")
        report_lines.append("")

# Save report
with open(os.path.join(data_dir, 'statistical_report.txt'), 'w') as f:
    f.write('\n'.join(report_lines))

# Function to generate plots with specified font sizes
def generate_plots(version="large"):
    # Set font sizes and output directory based on version
    if version == "large":
        current_base_font_size = large_base_font_size
        current_bar_label_font_size = large_bar_label_font_size
        current_x_tick_font_size = large_x_tick_font_size
        output_dir = large_dir
    elif version == "medium":
        current_base_font_size = medium_base_font_size
        current_bar_label_font_size = medium_bar_label_font_size
        current_x_tick_font_size = medium_x_tick_font_size
        output_dir = medium_dir
    else:  # small
        current_base_font_size = small_base_font_size
        current_bar_label_font_size = small_bar_label_font_size
        current_x_tick_font_size = small_x_tick_font_size
        output_dir = small_dir

    # Generate motif type distribution plot (including PK)
    plt.figure(figsize=figure_size, dpi=300)

    # Sort data in descending order
    sorted_data = sorted(zip(motif_counts.values, motif_counts.index), reverse=True)
    sorted_counts, sorted_motif_types = zip(*sorted_data)

    bars = plt.bar(sorted_motif_types, sorted_counts, color=colors)
    plt.xlabel("Motif Type", fontsize=current_base_font_size)
    plt.ylabel("Number of Motifs", fontsize=current_base_font_size)
    plt.xticks(fontsize=current_x_tick_font_size, rotation=45, ha='right')
    plt.yticks(fontsize=current_base_font_size - 2)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add data labels on bars
    max_count = max(sorted_counts)
    plt.ylim(0, max_count * 1.4)
    for bar, count in zip(bars, sorted_counts):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                 ha='center', va='bottom', fontsize=current_bar_label_font_size)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'motif_type_vs_number.png'), dpi=300)
    plt.close()

    # Generate PK distribution plot
    plt.figure(figsize=figure_size, dpi=300)

    # Sort PK data in descending order
    sorted_pk_data = sorted(zip(pk_counts.values, pk_counts.index), reverse=True)
    sorted_pk_counts, sorted_pk_types = zip(*sorted_pk_data)

    bars = plt.bar(sorted_pk_types, sorted_pk_counts, color=colors)
    plt.xlabel("Pseudoknot Type", fontsize=current_base_font_size)
    plt.ylabel("Number of Motifs", fontsize=current_base_font_size)
    plt.xticks(fontsize=current_x_tick_font_size, rotation=45, ha='right')
    plt.yticks(fontsize=current_base_font_size - 2)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add data labels on bars
    max_pk_count = max(sorted_pk_counts)
    plt.ylim(0, max_pk_count * 1.4)
    for bar, count in zip(bars, sorted_pk_counts):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                 ha='center', va='bottom', fontsize=current_bar_label_font_size)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pk_distribution.png'), dpi=300)
    plt.close()

    # Generate size distribution plots for each motif type (regular motifs only)
    motif_types = ['internal', 'hairpin', 'bulge', '3_junction', '4_junction', '5_junction']
    x_ranges = [(3, 25), (3, 22), (3, 10), (3, 50), (3, 50), (3, 50)]
    color_map = {'internal': colors[2], 'hairpin': colors[4], 'bulge': colors[5], 
                 '3_junction': colors[3], '4_junction': colors[1], '5_junction': colors[0]}

    for i, motif_type in enumerate(motif_types):
        if motif_type not in distribution.index:
            continue
            
        plt.figure(figsize=figure_size, dpi=300)
        
        # Filter data within specified range
        distribution_subset = distribution.loc[motif_type]
        distribution_subset = distribution_subset[distribution_subset.index >= x_ranges[i][0]]
        distribution_subset = distribution_subset[distribution_subset.index <= x_ranges[i][1]]
        
        # Adjust motif size values based on type
        if motif_type == 'hairpin':
            adjusted_motif_size = [size - 2 for size in distribution_subset.index]
        elif motif_type == '3_junction':
            adjusted_motif_size = [size - 6 for size in distribution_subset.index]
        else:
            adjusted_motif_size = [size - 4 for size in distribution_subset.index]
        
        plt.bar(adjusted_motif_size, distribution_subset.values, color=color_map[motif_type])
        
        # Format plot labels and axes
        display_name = motif_type.replace('_junction', '-way junction')
        plt.xlabel(f"{display_name.capitalize()} Size", fontsize=current_base_font_size)
        plt.ylabel("Number of Motifs", fontsize=current_base_font_size)
        
        # Set x-ticks to multiples of 5
        max_x = max(adjusted_motif_size) if adjusted_motif_size else 0
        xticks = np.arange(0, max_x + 5, 5)
        xticks = xticks[xticks <= max_x]
        if 0 not in xticks:
            xticks = np.append([0], xticks)
        
        plt.xticks(xticks, fontsize=current_base_font_size - 2)
        plt.yticks(fontsize=current_base_font_size - 2)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Set axis limits
        plt.xlim(left=-0.5, right=max_x + 0.5 if max_x > 0 else 0.5)
        plt.ylim(bottom=0)
        
        plt.savefig(os.path.join(output_dir, f'motif_size_distribution_{motif_type}.png'), dpi=300)
        plt.close()

# Generate all three versions
generate_plots("large")
generate_plots("medium") 
generate_plots("small")

print("All plots and statistical data generated successfully!")
print(f"Output base directory: {base_output_dir}")
print("Generated files:")
print("  - PNG plots in large/, medium/, small/ folders")
print("  - CSV data files and statistical_report.txt in data/ folder")
