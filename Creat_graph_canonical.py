# -*- coding: utf-8 -*-
"""
Changed from creatgraphv4. All the noncanonical bps will be ignored.
Created on Fri Mar  8 11:40:06 2024
For an edge, attribute=[a1,a2,a3] ai=0,1 1 means have, 0 means have not. 
I make the attribute to integral. The attribute of vertex are 0.
The attribute of e is 1 for covalent, 10 for noncanonical bp, 100 for canonical bp. If not only one, sum them. 
a1 means bp_cononical a2 means bp_noncononical a3 means covalent.
The nodes have no attribute.
@author: JingyiLi
"""
import networkx as nx
import os
import re
import pickle

# Keep node labels as-is without encoding
def encode_label(label):
    return label

# Convert a list of integers into a single concatenated integer
def list_to_number(lst):
    return int(''.join(str(num) for num in lst))

# Read and create a one-hot mapping from a text file
def read_data(data_doc_path):
    char_to_vector = {}
    try:
        with open(data_doc_path, 'r') as file:
            lines = file.readlines()
            for index, char in enumerate(lines):
                clean_char = char.strip()
                if clean_char:
                    char_to_vector[clean_char] = [int(i == index) for i in range(len(lines))]
    except FileNotFoundError:
        print(f"Error: The file {data_doc_path} was not found.")
    return char_to_vector

# Convert a string into a vector representation based on the mapping provided
def string2vector(line, char_to_vector, current_section):
    vector_sum = [0] * len(char_to_vector)
    parts = line.split(':')
    if len(parts) > 1:
        part1 = parts[1].strip()
        elements = part1.split(' ')
        if current_section == "base_pairs":
            filtered_elements = [elem for elem in elements if '/' in elem]
        else:
            filtered_elements = elements
        for elem in filtered_elements:
            vector = char_to_vector.get(elem, None)
            if vector:
                vector_sum = [sum(x) for x in zip(vector_sum, vector)]
    return vector_sum

# Extract node information from a line
def write_nt(line, current_section, filename):
    parts = line.strip().split()
    if len(parts) >= 4:
        vertex = parts[0]
        return vertex
    return None

# Extract edge information and attributes from a line
def write_interaction(line, current_section, filename, vertices):
    if ':' in line:
        encoded_results = []
        try:
            edge_info, attribute_info = line.strip().split(' : ')
            if len(edge_info.split('-')) == 2:
                parts = edge_info.split('-')
            else:
                parts = re.findall(r"[A-Za-z]-?\d+|'\d+'-?\d+", edge_info)
            if len(parts) != 2:
                raise ValueError("Error extracting vertices.")
            vertex1, vertex2 = parts
            if vertex1 in vertices and vertex2 in vertices:
                encoded_vector = string2vector(line, interaction_mapping, current_section)
                temp = sum(encoded_vector[21:])
                temp = min(temp, 1)
                encoded_vector = encoded_vector[:1]
                encoded_vector.extend([0, 0])
                encoded_results.append((vertex1, vertex2, encoded_vector))
        except ValueError as e:
            print(f"File '{filename}' error: '{line.strip()}' {e}")
    return encoded_results

# Parse the input file to extract vertices and edges
def get_graph_info(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
        return set(), []
    vertices = set()
    edges = []
    current_section = None
    for line in lines:
        line = line.strip()
        if line.startswith("Residue conformations"):
            current_section = "residues"
        elif line.startswith("Adjacent stackings"):
            current_section = "adjacent_stackings" 
        elif line.startswith("Base-pairs"):
            current_section = "base_pairs"
        elif current_section == "residues" and line:
            vertex = write_nt(line, current_section, filename)
            if vertex:
                vertices.add(vertex)
        elif current_section in ["base_pairs"] and '-' in line:
            edge_data = write_interaction(line, current_section, filename, vertices)
            edges.append(edge_data)
    return vertices, edges

# Create a graph object from vertices and edges
def create_graph(vertices, edges):
    G = nx.Graph()
    for vertex in vertices:
        G.add_node(vertex)
    for edge_list in edges:
        for edge in edge_list:
            vertex1, vertex2, edge_attribute = edge
            if edge_attribute != [0, 0, 0]:
                G.add_edge(vertex1, vertex2, attribute=edge_attribute)
    return G

# Generate next sequential nucleotide name
def nextnode(node_name):
    match = re.match(r"^(.*?)(\d+)$", node_name)
    if match:
        prefix = match.group(1)
        number = match.group(2)
        new_number = int(number) + 1  
        return f"{prefix}{new_number}"
    else:
        return None

# Add covalent edges between nodes based on naming convention
def add_covalent(G, filename):
    for node in G.nodes():
        next_node = nextnode(node)
        if next_node in G.nodes():
            if G.has_edge(node, next_node):
                if G[node][next_node]['attribute']:
                    G[node][next_node]['attribute'][-1] = 1
            else:
                G.add_edge(node, next_node, attribute=[0, 0, 1])
    return G

# Save the graph to a custom TXT format
def save_graph_to_txt(G, out_file, graph_number, filename):
    out_file.write(f't # {graph_number}\n')
    for node in G.nodes():
        out_file.write(f'v {encode_label(node)} 0\n')
    for u, v, data in G.edges(data=True):
        out_file.write(f'e {encode_label(u)} {encode_label(v)} {list_to_number(data["attribute"])}\n')

# Paths and setup for file processing
script_dir = os.path.dirname(os.path.abspath(__file__))
nt_list_path = os.path.join(script_dir,'data', 'unique_residue_conformations.txt')
interaction_list_path = os.path.join(script_dir,'data','unique_interactions_reduced.txt')
nt_mapping = read_data(nt_list_path)
interaction_mapping = read_data(interaction_list_path)

path = os.path.join(script_dir,'data','rna_annotation')
output_path = os.path.join(script_dir,'data','canonical_output_new.txt')
# Main processing loop to build graphs
graphs = {}
with open(output_path, 'w') as out_file:
    graph_number = 0
    for filename in os.listdir(path):
        if filename.endswith("_annotation.txt"):
            filepath = os.path.join(path, filename)
            vertices, edges = get_graph_info(filepath)
            G = create_graph(vertices, edges)
            G = add_covalent(G, filename)
            save_graph_to_txt(G, out_file, graph_number, filename)
            graphs[filename] = G
            graph_number += 1

print("Graphs have been saved to TXT successfully.")

# Save the graph data using pickle for later use
file_path_to_save = os.path.join(script_dir,'data','saved_graphs_canonical_new.pickle')
with open(file_path_to_save, 'wb') as file:
    pickle.dump(graphs, file)

print("Graphs have been saved successfully.")
