# -*- coding: utf-8 -*-
"""
The attribute of vertices is limited to 4 dim (only A T G C), and the edge character is limted into 3 dim. 
Created on Tue Mar 26 15:06:53 2024

@author: JingyiLi

Created on Fri Mar  8 11:40:06 2024
For an edge, attribute=[a1,a2,a3] ai=0,1 1 means have, 0 means have not. 
I make the attribute to integral. The attribute of vertex are 0.
The attribute of e is 1 for covalent, 10 for noncanonical bp, 100 for canonical bp. If not only one, sum them. 
a1 means bp_cononical a2 means bp_noncononical a3 means covalent.
The nodes have no attribute.
@author: JingyiLi
"""
import networkx as nx
import matplotlib.pyplot as plt
import os
import re
import pickle

"Possible A1, A-1,'1'1,'1'-1"
def encode_label(label):
    cleaned_label = label.replace("'", "")
    has_dash = '-' in cleaned_label
    import re
    parts = re.split('(-|\d+)', cleaned_label)
    parts = [p for p in parts if p]  
    prefix = '2' if has_dash else '1'
    encoded_parts = []
    for part in parts:
        if part.isalpha(): 
            ascii_value = ord(part)
            encoded_part = f'9{ascii_value:03}'
        elif part.isdigit():  
            encoded_part = part
        else: 
            continue
        encoded_parts.append(encoded_part)
    final_encoding = prefix + ''.join(encoded_parts)
    return final_encoding
def list_to_number(lst):
    number_str = ''.join(str(num) for num in lst)
    return int(number_str)


def read_data(data_doc_path): ########generate one-hot encoding
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
    except Exception as e:
        print(f"An error occurred while reading {data_doc_path}: {e}")
    return char_to_vector

# Input is the list of nt or interactions. return to a attribute of nt or interactions in vectior form
# The sequence doesn't matter, we represent each element as a vector. If the object has the element, 
# the correspnonding vacant is 1, if not, 0. 
def string2vector(line, char_to_vector, current_section):
    vector_sum = [0] * len(char_to_vector) # label vector
    parts = line.split(':') # divided the string by :
    if len(parts) > 1:
        part1 = parts[1].strip() # remove leading and trailing whitespaces
        elements = part1.split(' ')
        if current_section == "base_pairs":
            filtered_elements = [elem for elem in elements if '/' in elem]
        else:
            filtered_elements = elements
        for elem in filtered_elements: 
            vector = char_to_vector.get(elem, None)
            if vector:
                #Sum for all meet vectors, because it is treated no order.
                vector_sum = [sum(x) for x in zip(vector_sum, vector)] 
    return vector_sum

# Process a line to get the residue in vector form, including node and attribute.
def write_nt(line,current_section, filename):
    parts = line.strip().split()
    if len(parts) >= 4:  # the lengthe is to get the residue. keep :"A309 : A C3p_endo anti" skip:"F178 : ARG"
        vertex = parts[0]  # The name for nodes
        return vertex 
    else:
        return None
# write all interactions, I just don't want to change the name
def write_interaction(line, current_section, filename,vertices):
    if ':' in line:
        # Example of line: '0'79-'0'97 : G-G Ww/O2' pairing 
        encoded_results = []
        try:
            # edge info: '0'79-'0'97
            # attribute: G-G Ww/O2' pairing 
            edge_info, attribute_info = line.strip().split(' : ') 
            # Nodes can be A1 or A-1 or '0'97 or '0'-97
            if len(edge_info.split('-'))== 2: # example: A1-A2
                parts = edge_info.split('-')  # get two vertices
            else : 
                parts = re.findall(r"[A-Za-z]-?\d+|'\d+'-?\d+", edge_info)
            if len(parts) != 2:
                raise ValueError("error in getting 2 vertices。")
            vertex1, vertex2 = parts
            if vertex1 in vertices and vertex2 in vertices:
                encoded_vector = string2vector(line, interaction_mapping, current_section)
                # This part should be changed for different attribute form
                # a = [0,1,0] a[0] means whether bp or not a[1] means whether non-bp or not.
                # temp=1 if non-bp interactions 
                temp = sum(encoded_vector[21:])
                temp = min(temp,1)
                encoded_vector = encoded_vector[:1]
                encoded_vector.append(temp)
                # Need to append more if we want 2 more elements, the third is for covalent interaction
                encoded_vector.append(0)
                #print(encoded_vector)
                encoded_results.append((vertex1, vertex2, encoded_vector))  
        except ValueError as e:
            print(f"file '{filename}' error '{line.strip()}'{e}")
    return encoded_results

def get_graph_info(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
        return set(), []
    vertices = set()  # dictionary to store the nodes.
    edges = []
    current_section = None  # To record the current section.
    for line in lines:
        line = line.strip()
        #print(line)
        if line.startswith("Residue conformations"):
            current_section = "residues"
        elif line.startswith("Adjacent stackings"):
            current_section = "adjacent_stackings"
        elif line.startswith("Non-Adjacent stackings"):
            current_section = "non_adjacent_stackings"
        elif line.startswith("Base-pairs"):
            current_section = "base_pairs"
        elif current_section == "residues" and line:
            vertex = write_nt(line,current_section, filename)
            if vertex:
                vertices.add(vertex)
        #elif current_section in ["adjacent_stackings", "non_adjacent_stackings", "base_pairs"] and '-' in line:
        elif current_section in ["base_pairs"] and '-' in line:
            edge_data = write_interaction(line,current_section, filename,vertices)
            #print("edgedata:",edge_data)
            edges.append(edge_data)            
    return vertices, edges

def create_graph(vertices, edges):
    G = nx.Graph()
    for vertex in vertices:
        G.add_node(vertex)
    for edge_list in edges: # edges is a list containing a edge_data
        for edge in edge_list:  # in this code, only one edge in edge_list, check line 95.
            vertex1, vertex2, edge_attribute = edge
            G.add_edge(vertex1, vertex2, attribute = edge_attribute)
    return G

def nextnode(node_name):
    match = re.match(r"^(.*?)(\d+)$", node_name)
    if match:
        prefix = match.group(1)
        number = match.group(2)
        new_number = int(number) + 1  
        return f"{prefix}{new_number}"
    else:
        return None  
    
def add_covalent(G,filename):
    for node in G.nodes():
        next_node = nextnode(node)
        if next_node in G.nodes():  
            if G.has_edge(node, next_node):
                if G[node][next_node]['attribute']:
                    G[node][next_node]['attribute'][-1] = 1
            else:
                G.add_edge(node, next_node, attribute=[0,0,1])
    return G

def visualize_graph(G):
    plt.figure(figsize=(10, 8))  
    pos = nx.spring_layout(G)  
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, edge_color='k', linewidths=1, font_size=15, alpha=0.7)

    plt.title("Graph of Residues and Interactions", fontsize=40)
    plt.show()
    
def save_graph_to_txt(G, out_file, graph_number, filename):
    out_file.write(f't # {graph_number}\n')
    for node, data in G.nodes(data=True):
        attribute = data.get("list_to_number(attribute)", "DEFAULT") 
        out_file.write(f'v {encode_label(node)} {0}\n')
    for u, v, data in G.edges(data=True):
        attribute = data.get("attribute", "DEFAULT")  
        out_file.write(f'e {encode_label(u)} {encode_label(v)} {list_to_number(attribute)}\n')
        
# Genetate the mappings
script_dir = os.path.dirname(os.path.abspath(__file__))
nt_list_path = os.path.join(script_dir, 'data', 'unique_residue_conformations.txt')
nt_mapping = read_data(nt_list_path)
interaction_list_path = os.path.join(script_dir, 'data', 'unique_interactions.txt')
interaction_mapping = read_data(interaction_list_path)

#This is the path to read MC output date to build graphs
path = os.path.join(script_dir, 'data', 'rna_annotation')
# path = 'D:/test'
output_path = os.path.join(script_dir, 'data', 'created_graphs_noncanonical.txt')

graphs = {}
batch_size = 100
batch_count = 0

with open(output_path, 'w') as out_file:
    graph_number = 0
    for filename in os.listdir(path):
        if filename.endswith(".txt"):
            filepath = os.path.join(path, filename)
            vertices, edges = get_graph_info(filepath)
            G = create_graph(vertices, edges)
            G = add_covalent(G, filename)
            save_graph_to_txt(G, out_file, graph_number, filename)

            graphs[filename] = G
            graph_number += 1
            
            # Batch processing logic
            if graph_number % batch_size == 0:
                batch_file = os.path.join(script_dir, 'data', f'batch_{batch_count:04d}_graphs.pickle')
                with open(batch_file, 'wb') as f:
                    pickle.dump(graphs, f)
                print(f"Saved batch {batch_count} with {len(graphs)} graphs")
                graphs.clear()  # Clear memory
                batch_count += 1
            
            # Print progress every 100 files
            if graph_number % 100 == 0:
                print(f"Progress: {graph_number} files processed")
# Save final batch if there are remaining graphs
if graphs:
    batch_file = os.path.join(script_dir, 'data', f'batch_{batch_count:04d}_graphs.pickle')
    with open(batch_file, 'wb') as f:
        pickle.dump(graphs, f)
    print(f"Saved final batch {batch_count} with {len(graphs)} graphs")

print("Graphs have been saved to TXT successfully.")

# Function to merge all batch files
def merge_all_batches():
    print("Merging all batch files...")
    all_graphs = {}
    batch_num = 0
    
    while True:
        batch_file = os.path.join(script_dir, 'data', f'batch_{batch_num:04d}_graphs.pickle')
        if not os.path.exists(batch_file):
            break
        
        with open(batch_file, 'rb') as f:
            batch_graphs = pickle.load(f)
            all_graphs.update(batch_graphs)
        
        print(f"Loaded batch {batch_num} with {len(batch_graphs)} graphs")
        batch_num += 1
    
    # Save merged file
    final_file = os.path.join(script_dir, 'data', 'saved_graphs_noncanonical20250627.pickle')
    with open(final_file, 'wb') as f:
        pickle.dump(all_graphs, f)
    
    print(f"Merged {len(all_graphs)} graphs into final file")
    return all_graphs

# Merge all batches and save final file
merged_graphs = merge_all_batches()
print("All graphs have been saved successfully.")

