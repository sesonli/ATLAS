# -*- coding: utf-8 -*-
"""
Created from existing find_motif.py and find_hairpin.py scripts
Combined to find all motif types: hairpins, bulges, internal loops, and junctions
This script reuses the existing tested functions and adds functionality to process multiple motif types
@author: JingyiLi
"""

from networkx.algorithms import isomorphism
import time
import pickle
import os
import sqlite3
import re
import networkx as nx

# ========== UTILITY FUNCTIONS FROM EXISTING SCRIPTS ==========

def get_chain_type(node_name):
    if re.match(r'^[A-Z]\d+$', node_name):
        return 'A'
    elif re.match(r'^[A-Z]-\d+$', node_name):
        return 'A-'
    elif re.match(r"^'\d+'-\d+$", node_name):
        return "0'-"
    elif re.match(r"^'\d+'\d+$", node_name):
        return "0'"
    else:
        return None

def natural_key(node_name):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', node_name)]

def is_paired_node(G, node):
    for neighbor in G.neighbors(node):
        if G[node][neighbor]['attribute'] != [0, 0, 1]:
            return True
    return False

def is_unpaired_node(G, node):
    for neighbor in G.neighbors(node):
        if G[node][neighbor]['attribute'] != [0, 0, 1]:
            return False
    return True

def get_chains(lst):
    chains_dict = {}
    current_chain = [lst[0]]

    def parse_item(item):
        match = re.match(r"([A-Za-z]*)'?-?(\d+)", item)
        if match:
            letter_part = match.group(1)
            number_part = int(match.group(2))
            return letter_part, number_part
        return None, None

    chain_index = 1

    for i in range(1, len(lst)):
        prev_letter, prev_number = parse_item(lst[i-1])
        curr_letter, curr_number = parse_item(lst[i])

        if prev_letter == curr_letter and curr_number == prev_number + 1:
            current_chain.append(lst[i])
        else:
            chains_dict[f"Chain {chain_index}"] = current_chain
            current_chain = [lst[i]]
            chain_index += 1

    chains_dict[f"Chain {chain_index}"] = current_chain
    return chains_dict

# ========== GRAPH COMPRESSION FUNCTIONS FROM find_motif.py ==========

def compress_graph(G, graph_name, problematic_graphs):
    compressed_G = nx.Graph()
    lst = []
    try:
        for node in G.nodes():            
            lst.append(node)
        lst = sorted(lst, key=natural_key)
        chains = get_chains(lst)

        for u, v, data in G.edges(data=True):
            if data['attribute'] != [0, 0, 1]:
                compressed_G.add_edge(u, v, attribute=data['attribute'])
                
        for chain_type, chain in chains.items():
            chain = sorted(chain, key=natural_key)
            
            i = 0
            while i < len(chain) - 1:
                node1 = chain[i]

                if is_paired_node(G, node1) and is_unpaired_node(G, chain[i + 1]):
                    j = i + 1
                    while j < len(chain) - 1 and is_unpaired_node(G, chain[j]):
                        j += 1
                    if is_paired_node(G, chain[j]):
                        nodes_to_compress = chain[i + 1:j]
                        compressed_G.add_edge(node1, chain[j], weight=nodes_to_compress, attribute=[0, 0, 2])
                        i = j
                    else:
                        i += 1
                else:
                    if not G.has_edge(node1, chain[i + 1]):
                        print(f"No edge between {node1} and {chain[i + 1]}, skipping.")
                    else:
                        compressed_G.add_edge(node1, chain[i + 1], attribute=G[node1][chain[i + 1]]['attribute'])
                    i += 1

            if len(chain) > 1:
                node1 = chain[-2]
                node2 = chain[-1]
                if G.has_edge(node1, node2):
                    compressed_G.add_edge(node1, node2, attribute=G[node1][node2]['attribute'])
                else:
                    print(f"No edge between last two nodes {node1} and {node2}, skipping.")

    except (KeyError, IndexError) as e:
        print(f"Skipping graph {graph_name} due to error: {e}")
        problematic_graphs.append(graph_name)
        return None    

    return compressed_G

def depress_graph(compressed_G):
    decompressed_G = nx.Graph()

    for edge in compressed_G.edges(data=True):
        node1, node2, data = edge

        if natural_key(node1) > natural_key(node2):
            node1, node2 = node2, node1
        
        if 'weight' in data:
            compressed_nodes = data['weight']            
            sorted_compressed_nodes = sorted(compressed_nodes, key=natural_key)            
            recover = [node1] + sorted_compressed_nodes + [node2]           
            for i in range(len(recover) - 1):
                decompressed_G.add_edge(recover[i], recover[i+1], attribute=[0, 0, 1])
        else:
            decompressed_G.add_edge(node1, node2, attribute=data['attribute'])
    return decompressed_G

# ========== PDB READING FUNCTIONS FROM EXISTING SCRIPTS ==========

def read_pdb_and_find_nt(pdb_filepath, node_id):
    node_id = node_id.replace("'", "")
    content = []
    with open(pdb_filepath, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                current_chain = line[21:22].strip()
                current_nt = line[22:27].strip()
                node_id_in_file = current_chain + current_nt
                if node_id_in_file == node_id:
                    content.append(line.strip())
    return content

def extract_node_info(node_id, pdb_filename):
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    pdbrna_dir = os.path.join(data_dir, 'rna_chain_corrected_2')
    pdb_filename_base = pdb_filename[:4]  
    correct_pdb_filename = pdb_filename_base + '.pdb'
    pdb_filepath = os.path.join(pdbrna_dir, correct_pdb_filename)
    if not os.path.exists(pdb_filepath):
        return f"File {correct_pdb_filename} not found in rna_chain_corrected_2 folder."
    content = read_pdb_and_find_nt(pdb_filepath, node_id)
    return "\n".join(content) if content else "No matching nodes found"

# ========== DATABASE FUNCTIONS FROM EXISTING SCRIPTS ==========

def save_to_database(motif_type, graph_id, reduced_nt_numbers, combined_nt_numbers, combined_filecontent, conn):
    cursor = conn.cursor()
    pdbid = graph_id[:4]
    cursor.execute('''
    INSERT INTO files (motif_type, pdbid, paired_nt_number, nt_number, filecontent)
    VALUES (?, ?, ?, ?,?)
    ''', (motif_type, pdbid, reduced_nt_numbers, combined_nt_numbers, combined_filecontent))
    conn.commit()

# ========== ISOMORPHISM FUNCTIONS FROM EXISTING SCRIPTS ==========

def is_subgraph_isomorphic(graph, target_graph):
    GM = isomorphism.GraphMatcher(graph, target_graph, 
                                  node_match=lambda n1,n2: True,
                                  edge_match=lambda e1, e2: e1['attribute'] == e2['attribute'])
    
    matches = []
    matched_subgraphs = []

    for mapping in GM.subgraph_isomorphisms_iter():
        sorted_mapping = tuple(sorted(mapping.items()))

        if sorted_mapping not in matches:
            matches.append(sorted_mapping)
            subgraph_nodes = mapping.keys()
            matched_subgraph = graph.subgraph(subgraph_nodes).copy()
            
            def make_hashable(attr):
                return {k: tuple(v) if isinstance(v, list) else v for k, v in attr.items()}
            
            sorted_edges = tuple(sorted((u, v) for u, v, attr in matched_subgraph.edges(data=True)))
            matched_subgraphs.append(sorted_edges)
    
    if matches:
        unique_subgraphs = [graph.edge_subgraph(edges) for edges in matched_subgraphs]
        return matches, unique_subgraphs
    return None, None

def sort_key(node):
    node = node.replace('-', '')
    match = re.match(r"([A-Za-z]+)(\d+)", node)
    if not match:
        match = re.match(r"'(\d+)'(\d+)", node)
    if match:
        letter = match.group(1)
        number = int(match.group(2))
        return (letter, number)
    return ('', 0)

# ========== MODIFIED SEARCH FUNCTIONS FOR ALL MOTIFS ==========

def find_hairpin_motifs(graphs, target_graphs, conn):
    """Search for hairpin motifs (from find_hairpin.py logic)"""
    print("Searching for hairpin motifs...")
    matching_counts = {}
    target_match_counts = {tg_id: 0 for tg_id in target_graphs.keys()}
    total_graphs = len(graphs)
    graph_processed = 0
    matches_info = []
    problematic_graphs = []
    
    for graph_id, graph in graphs.items():   
        graph_processed += 1
        count = 0
        start_time = time.time()        
        for target_graph_id, target_graph in target_graphs.items():
            mappings, matched_subgraphs = is_subgraph_isomorphic(graph, target_graph)
            if matched_subgraphs:
                for mapping, matched_graph in zip(mappings, matched_subgraphs):
                    combined_filecontent = []
                    reduced_nt_numbers = []
                    combined_nt_numbers = []
                    mapping = sorted(mapping, key=lambda x: sort_key(x[1]))
                    for a_node, t_node in mapping:
                        reduced_nt_numbers.append(a_node)
                    for nodes in matched_graph:
                        nt_number = nodes
                        filecontent = extract_node_info(nt_number, graph_id)
                        combined_nt_numbers.append(nt_number)
                        combined_filecontent.append(filecontent)
                    combined_nt_numbers = sorted(combined_nt_numbers, key=sort_key)
                    reduced_nt_numbers = [combined_nt_numbers[0], combined_nt_numbers[-1]]
                    combined_nt_numbers_str = ','.join(combined_nt_numbers)
                    reduced_nt_numbers_str = ','.join(reduced_nt_numbers)
                    combined_filecontent_str = '\n'.join(combined_filecontent)
                    
                    save_to_database('hairpin', graph_id, reduced_nt_numbers_str, combined_nt_numbers_str, combined_filecontent_str, conn)
                    count += 1
                    target_match_counts[target_graph_id] += 1
                    
        matching_counts[graph_id] = count
        elapsed_time = time.time() - start_time
        if graph_processed % 100 == 0:
            print(f"Finished searching graph {graph_id} ({graph_processed}/{total_graphs}). Time spent: {elapsed_time:.2f} seconds.")
   
    print(f"Total graphs processed: {graph_processed}/{total_graphs}")
    return matching_counts, matches_info, target_match_counts, problematic_graphs

def find_other_motifs(graphs, target_graphs, conn, motif_type_name):
    """Search for other motifs (from find_motif.py logic)"""
    print(f"Searching for {motif_type_name} motifs...")
    matching_counts = {}
    target_match_counts = {tg_id: 0 for tg_id in target_graphs.keys()}
    total_graphs = len(graphs)
    graph_processed = 0
    matches_info = []
    problematic_graphs = []
    
    for graph_id, graph in graphs.items():   
        compressed_G = compress_graph(graph, graph_id, problematic_graphs)
        if compressed_G is None:
            continue
        graph_processed += 1
        count = 0
        start_time = time.time()        
        for target_graph_id, target_graph in target_graphs.items():
            mappings, matched_subgraphs = is_subgraph_isomorphic(compressed_G, target_graph)
            if matched_subgraphs:
                for mapping, matched_graph in zip(mappings, matched_subgraphs):
                    combined_filecontent = []
                    reduced_nt_numbers = []
                    combined_nt_numbers = []
                    mapping = sorted(mapping, key=lambda x: sort_key(x[1]))
                    for a_node, t_node in mapping:
                        reduced_nt_numbers.append(a_node)
                    matched_graph_depressed = depress_graph(matched_graph)                    
                    for nodes in matched_graph_depressed:
                        nt_number = nodes
                        filecontent = extract_node_info(nt_number, graph_id)
                        combined_nt_numbers.append(nt_number)
                        combined_filecontent.append(filecontent)                    
                    combined_nt_numbers_str = ','.join(combined_nt_numbers)
                    reduced_nt_numbers_str = ','.join(reduced_nt_numbers)
                    combined_filecontent_str = '\n'.join(combined_filecontent)
                    
                    save_to_database(target_graph_id, graph_id, reduced_nt_numbers_str, combined_nt_numbers_str, combined_filecontent_str, conn)
                    count += 1
                    target_match_counts[target_graph_id] += 1
                    
        matching_counts[graph_id] = count
        elapsed_time = time.time() - start_time
        if graph_processed % 100 == 0:
            print(f"Finished searching graph {graph_id} ({graph_processed}/{total_graphs}). Time spent: {elapsed_time:.2f} seconds.")
   
    print(f"Total graphs processed: {graph_processed}/{total_graphs}")
    return matching_counts, matches_info, target_match_counts, problematic_graphs

# ========== MAIN EXECUTION ==========

if __name__ == "__main__":
    # Setup directories
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, 'data')
    os.makedirs(data_dir, exist_ok=True)

    # Create database
    db_path = os.path.join(data_dir, 'all_motifs.db')
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS files (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            motif_type TEXT,
            pdbid TEXT,
            paired_nt_number TEXT, 
            nt_number TEXT,
            filecontent TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        ''')

        # Load graphs
        file_path = os.path.join(data_dir, "saved_graphs_canonical_new.pickle")
        if not os.path.exists(file_path):
            print(f"Error: {file_path} not found!")
            exit(1)
            
        with open(file_path, 'rb') as f:
            graphs = pickle.load(f)
        print(f"Loaded {len(graphs)} graphs")

        # Define motif files to process
        motif_files = {
            'hairpin': 'hairpin.pickle',
            'bulge_internal_junction': 'bulge_internal_junction.pickle'
        }

        all_results = {}
        total_start_time = time.time()

        # Process each motif type
        for motif_name, pickle_file in motif_files.items():
            target_file_path = os.path.join(data_dir, pickle_file)
            if not os.path.exists(target_file_path):
                print(f"Warning: {target_file_path} not found, skipping {motif_name}")
                continue
                
            print(f"\n{'='*50}")
            print(f"Processing {motif_name} motifs...")
            print(f"{'='*50}")
            
            with open(target_file_path, 'rb') as f:
                target_graphs = pickle.load(f)
            print(f"Loaded {len(target_graphs)} target graphs for {motif_name}")

            motif_start_time = time.time()
            
            if motif_name == 'hairpin':
                matching_counts, matches_info, target_match_counts, problematic_graphs = find_hairpin_motifs(graphs, target_graphs, conn)
            else:
                matching_counts, matches_info, target_match_counts, problematic_graphs = find_other_motifs(graphs, target_graphs, conn, motif_name)
            
            motif_end_time = time.time()
            
            all_results[motif_name] = {
                'matching_counts': matching_counts,
                'target_match_counts': target_match_counts,
                'problematic_graphs': problematic_graphs,
                'time_taken': motif_end_time - motif_start_time
            }
            
            print(f"Completed {motif_name} in {motif_end_time - motif_start_time:.2f} seconds")
            print(f"Target graph matches: {target_match_counts}")

        total_end_time = time.time()
        
        # Write summary results
        summary_file = os.path.join(data_dir, "all_motifs_summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"All Motifs Search Summary\n")
            f.write(f"Total time: {total_end_time - total_start_time:.2f} seconds\n")
            f.write(f"Total graphs processed: {len(graphs)}\n\n")
            
            for motif_name, results in all_results.items():
                f.write(f"{motif_name.upper()} RESULTS:\n")
                f.write(f"Time taken: {results['time_taken']:.2f} seconds\n")
                f.write(f"Target graph matches: {results['target_match_counts']}\n")
                f.write(f"Problematic graphs: {len(results['problematic_graphs'])}\n")
                f.write(f"\n")

        print(f"\nAll motif searches completed in {total_end_time - total_start_time:.2f} seconds")
        print(f"Results summary written to {summary_file}") 