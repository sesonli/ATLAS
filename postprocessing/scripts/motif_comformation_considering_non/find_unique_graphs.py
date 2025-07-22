"""
Find unique graphs incorporating non-WC interactions, clustering the same type of motifs based on non-WC interactions.
This script can analyze motifs with either a specific size constraint or no size constraint.
"""

import sqlite3
import networkx as nx
from networkx.algorithms import isomorphism
import json
import matplotlib.pyplot as plt
import os
from datetime import datetime
import sys
from io import StringIO
import numpy as np

# Connect to the SQLite database
conn = sqlite3.connect('data/ATLAS2025.db')
cursor = conn.cursor()

def plot_graph_distribution_pie(frequencies, output_dir, motif_type, size_constraint):
    """
    Plot a pie chart showing the percentage distribution of graphs based on their frequencies.
    
    Parameters:
    frequencies (list): List of frequencies for each graph, ordered from most to least frequent.
    output_dir (str): Directory to save the output PNG file.
    motif_type (str): Type of motif for naming the output file.
    size_constraint (str): Size constraint information for the filename.
    """
    # Calculate the total frequency
    total = sum(frequencies)
    
    # Calculate percentages
    percentages = [(freq / total) * 100 for freq in frequencies]
    
    # Prepare data for the pie chart
    labels = []
    sizes = []
    
    # Process the top 4 graphs individually
    for i in range(min(4, len(percentages))):
        labels.append(f"Graph {i+1}")
        sizes.append(percentages[i])
    
    # Group graphs 5-8 as requested
    if len(percentages) > 4:
        g5_8_sum = sum(percentages[4:8]) if len(percentages) >= 8 else sum(percentages[4:])
        if g5_8_sum > 0:
            labels.append("Graphs 5-8")
            sizes.append(g5_8_sum)
    
    # Add "Other" category for all remaining graphs (9+)
    if len(percentages) > 8:
        other_sum = sum(percentages[8:])
        if other_sum > 0:
            labels.append("Other")
            sizes.append(other_sum)
    
    # Create the pie chart
    plt.figure(figsize=(12, 9))
    explode = [0.05] * len(sizes)  # Explode all slices slightly
    # Use consistent color theme with Fig2
    fig2_colors = ['#264653', '#287171','#2a9d8f', '#e9c46a', '#f4a261', '#e76f51']
    colors = [fig2_colors[i % len(fig2_colors)] for i in range(len(sizes))]
    
    # Use a doubled font size and only show percentages for slices ≥15%
    wedges, texts, autotexts = plt.pie(
        sizes, 
        explode=explode, 
        labels=labels, 
        colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p >= 15 else '',  # Only show percentage if ≥15%
        shadow=False,  # Removed shadow as requested
        startangle=90,
        textprops={'fontsize': 36}  # Doubled font size from 18 to 36
    )
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_weight('bold')
    
    # Add title and adjust layout
    plt.axis('equal')  # Equal aspect ratio ensures the pie chart is circular
    plt.title(f"Graph Distribution for {motif_type}", fontsize=48, pad=20)  # Doubled font size from 24 to 48
    
    # Save figure with descriptive filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"pie_chart_{motif_type}_{size_constraint}_{timestamp}.png"
    filepath = os.path.join(output_dir, filename)
    
    plt.tight_layout()
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Pie chart saved to: {filepath}")

# Function to calculate the number of nucleotides based on comma-separated string
def calculate_motif_size(nt_numbers):
    """Count the number of nucleotides, assuming they are separated by commas"""
    return len(nt_numbers.split(','))

# Function to parse the edges from the non_graph_edges field
def parse_edges(edge_list):
    """Parse edge list into a format suitable for NetworkX graph"""
    edges = []
    for edge in edge_list:
        try:
            node1, node2, attributes = edge
            edges.append((node1, node2, {'attribute': attributes['attribute']}))
        except ValueError:
            print(f"Error parsing edge: {edge}")
    return edges

# Function to plot a graph and display edge attributes
def plot_graph_with_attributes(graph, ax, graph_id, frequency):
    pos = nx.spring_layout(graph)  # Use spring layout for better visual separation of nodes
    #pos = nx.planar_layout(graph)  # Use spring layout for better visual separation of nodes
    nx.draw(graph, pos, ax=ax, with_labels=True, node_color="lightblue", node_size=500, font_size=10)

    # Draw edge labels (attributes)
    edge_labels = nx.get_edge_attributes(graph, 'attribute')
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_color='red', ax=ax)

    ax.set_title(f"Graph {graph_id} (Freq: {frequency})")
    
# Custom output stream that writes to both console and file
class TeeOutput(object):
    def __init__(self, file_stream):
        self.file_stream = file_stream
        self.console = sys.stdout
        
    def write(self, message):
        self.console.write(message)
        self.file_stream.write(message)
        
    def flush(self):
        self.console.flush()
        self.file_stream.flush()
    
# Main function to find unique graphs and plot them
def find_unique_graphs_and_plot(motif_type, nt_number_size=None, top_n=12, show_plot=False):
    """
    Find unique graphs based on motif_type and nt_number size, and plot the top N most frequent graphs.
    
    Parameters:
    motif_type (str): Type of motif to analyze (e.g., 'internal', 'hairpin', etc.)
    nt_number_size (int, optional): If specified, only include motifs of this exact size.
                                   If None or 0, include motifs of all sizes.
    top_n (int): Number of top frequent graphs to plot
    show_plot (bool): Whether to display the plot interactively
    """
    
    # Create output directory if it doesn't exist
    output_dir = 'postprocessing/results2025/unique_graphs'
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a descriptive filename
    size_str = f"size_{nt_number_size}" if nt_number_size else "all_sizes"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"unique_graphs_{motif_type}_{size_str}_{timestamp}.txt")
    
    # Open file for writing results
    with open(output_file, 'w') as f:
        # Redirect stdout to both console and file
        original_stdout = sys.stdout
        sys.stdout = TeeOutput(f)
        
        print(f"Unique Graphs Analysis for {motif_type}")
        print(f"Motif Size Constraint: {nt_number_size if nt_number_size else 'All Sizes'}")
        print(f"Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*50 + "\n")
        
        try:
            # Query to get all rows with the specified motif_type
            query = f"""
            SELECT nt_number, non_graph_edges FROM data
            WHERE motif_type = ?
            """
            cursor.execute(query, (motif_type,))
            rows = cursor.fetchall()
            print(f"Retrieved {len(rows)} records of type '{motif_type}' from database")

            # A list to store unique graphs and their frequencies
            unique_graphs = []
            graph_frequencies = []  # To track frequency of each graph
            graph_num = 0
            skipped = 0
            
            # Iterate through each row and filter based on the size of nt_number
            for row in rows:
                nt_number = row[0]  # Retrieve the nt_number field
                non_graph_edges = row[1]  # Retrieve the non_graph_edges field
                
                # Calculate the size of the motif
                motif_size = calculate_motif_size(nt_number)
                
                # Filter by size if specified
                if nt_number_size is not None and nt_number_size > 0:
                    if motif_size != nt_number_size:
                        skipped += 1
                        continue
                else:
                    # Always require at least 2 nucleotides (motif_size > 1)
                    if motif_size <= 1:
                        skipped += 1
                        continue
                
                graph_num += 1
                # Ensure non_graph_edges is deserialized if it's stored as a string
                if isinstance(non_graph_edges, str):
                    try:
                        non_graph_edges = json.loads(non_graph_edges)
                    except json.JSONDecodeError as e:
                        print(f"Error decoding non_graph_edges: {e}")
                        continue

                # Construct the graph using the parsed edges
                graph = nx.Graph()
                edges = parse_edges(non_graph_edges)
                graph.add_edges_from(edges)

                # Check if this graph is isomorphic to any graph in unique_graphs
                is_unique = True
                for idx, unique_graph in enumerate(unique_graphs):
                    GM = isomorphism.GraphMatcher(graph, unique_graph, 
                                                node_match=lambda n1, n2: True,  # Ignore node labels
                                                edge_match=lambda e1, e2: e1['attribute'] == e2['attribute'])  # Match edge attributes
                    if GM.is_isomorphic():
                        is_unique = False
                        graph_frequencies[idx] += 1  # Increment the frequency of the corresponding graph
                        break

                # If the graph is unique, add it to the list of unique graphs
                if is_unique:
                    unique_graphs.append(graph)
                    graph_frequencies.append(1)  # Initialize the frequency of this new graph as 1

            # Sort the graphs by frequency in descending order
            sorted_graphs_and_freqs = sorted(zip(unique_graphs, graph_frequencies), key=lambda x: x[1], reverse=True)

            # Select the top N most frequent graphs
            top_graphs_and_freqs = sorted_graphs_and_freqs[:top_n]
            
            # Calculate total graphs
            total_graphs = graph_num
            top_graphs_count = sum(freq for _, freq in top_graphs_and_freqs)
            other_graphs_count = total_graphs - top_graphs_count
            
            # Print the analysis results
            print(f"\nAnalysis Results:")
            print(f"Total graphs processed: {graph_num}")
            print(f"Total graphs skipped: {skipped}")
            print(f"Number of unique graphs found: {len(unique_graphs)}")
            print(f"Number of top graphs plotted (Top {top_n}): {top_graphs_count}")
            print(f"Number of other graphs not plotted: {other_graphs_count}")
            
            print("\nDetailed Frequency Information:")
            for i, (graph, frequency) in enumerate(top_graphs_and_freqs):
                print(f"Graph {i+1} has been found {frequency} time(s) ({frequency/total_graphs*100:.2f}%)")
                
                # Print edge information for each graph
                print(f"  Edges in Graph {i+1}:")
                for edge in graph.edges(data=True):
                    print(f"    {edge[0]} -- {edge[1]} : {edge[2]['attribute']}")
                print()

            print(f"\nResults saved to: {output_file}")

            # Plot the top N most frequent graphs
            rows_needed = (min(len(top_graphs_and_freqs), top_n) + 3) // 4  # Calculate rows needed
            fig, axes = plt.subplots(rows_needed, 4, figsize=(15, 5 * rows_needed))
            
            # Make axes 2D if it's 1D
            if rows_needed == 1:
                axes = [axes]
                
            # Flatten if we have a 2D array of axes
            if rows_needed > 1:
                axes = axes.flatten()
                
            for i, (graph, frequency) in enumerate(top_graphs_and_freqs):
                if i < top_n:
                    ax_idx = i
                    plot_graph_with_attributes(graph, axes[ax_idx], i+1, frequency)
            
            # Hide empty subplots
            for j in range(len(top_graphs_and_freqs), rows_needed * 4):
                if j < len(axes):
                    axes[j].axis('off')
            
            # Save the plot
            plt.tight_layout()
            plot_file = os.path.join(output_dir, f"unique_graphs_{motif_type}_{size_str}_{timestamp}.png")
            plt.savefig(plot_file, dpi=300)
            print(f"Plot saved to: {plot_file}")
            if show_plot:
                plt.show()
            else:
                plt.close()
                
            # Create and save the pie chart
            # Extract frequencies from sorted_graphs_and_freqs for the pie chart
            all_frequencies = [freq for _, freq in sorted_graphs_and_freqs]
            print("\nGenerating pie chart of graph distribution...")
            plot_graph_distribution_pie(all_frequencies, output_dir, motif_type, size_str)
            
        finally:
            # Restore stdout
            sys.stdout = original_stdout

# Run the analyses sequentially with show_plot=False for all but the last one
find_unique_graphs_and_plot(motif_type='hairpin', nt_number_size=6, top_n=8, show_plot=False)
find_unique_graphs_and_plot(motif_type='hairpin', nt_number_size=None, top_n=8, show_plot=False)
find_unique_graphs_and_plot(motif_type='internal', nt_number_size=None, top_n=8, show_plot=False)

# Close the database connection
conn.close()