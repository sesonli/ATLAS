# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:54:08 2024

@author: JingyiLi
"""
import re
import pickle
import sqlite3

def parse_node_name(item):
    match = re.match( r"([A-Za-z]|'?\d+'?)'?-?(\d+)", item)
    if match:
        letter_part = match.group(1)
        number_part = int(match.group(2))
        return letter_part, number_part

def graph_to_sequence_and_pairs(G):
    chains = {}
    for node in G.nodes:
        chain_id = node[0]
        if chain_id not in chains:
            chains[chain_id] = []
        chains[chain_id].append(node)

    sequences_and_pairs = []
    for chain_id, chain_nodes in chains.items():
        chain_nodes.sort(key=lambda x: parse_node_name(x))
        sequence = []
        pairs = []
        node_to_index = {}
        index = 1

        for node in chain_nodes:
            sequence.append(node)
            node_to_index[node] = index
            pairs.append([index, 0])
            index += 1

        for u, v, attr in G.edges(data="attribute"):
            if u in chain_nodes and v in chain_nodes and attr == [1, 0, 0]:
                u_index = node_to_index[u]
                v_index = node_to_index[v]
                pairs[u_index - 1][1] = v_index
                pairs[v_index - 1][1] = u_index

        sequences_and_pairs.append((sequence, pairs))
    return sequences_and_pairs

def determine_level(ss, a, b, pair_symbols):
    if ss[a] != '.' or ss[b] != '.':
        return -1
    level = 0
    while True:
        score = 0
        for i in range(a + 1, b):
            if ss[i] == pair_symbols[level][0]:
                score += 1
            elif ss[i] == pair_symbols[level][1]:
                score -= 1
        if score == 0:
            return level
        level += 1
        if level >= len(pair_symbols):
            return -1

def convert_to_dot_bracket(sequence, pairs):
    n = len(sequence)
    ss = ['.' for _ in range(n)]
    pair_symbols = [('(', ')'), ('[', ']')]
    id2index = {pair[0]: i for i, pair in enumerate(pairs)}
    for p in pairs:
        if p[1] != 0:
            i = id2index[p[0]]
            j = id2index[p[1]]
            level = determine_level(ss, i, j, pair_symbols)
            if level != -1:
                ss[i] = pair_symbols[level][0]
                ss[j] = pair_symbols[level][1]

    dot_bracket_str = ''.join(ss)
    return sequence, dot_bracket_str

def compress(lst):
    if not lst:
        return []
    compressed = [lst[0]]
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1]:
            compressed.append(lst[i])
    return compressed

def type_kp(lst):
    lst_c = compress(lst)
    if lst_c == [1]:
        return "H"
    elif lst_c == [1, -1]:
        return "HHH"
    elif lst_c == [1, 0]:
        return "Hlout"
    elif lst_c == [1, 0, 1]:
        return "Hlin"
    elif lst_c == [0, 1]:
        return "LL"
    else:
        return "LR"

def find_all_subsequences(dot_bracket_str):
    def is_balanced(subsequence):
        stack = []
        for char in subsequence:
            if char == '(':
                stack.append('(')
            elif char == ')':
                if not stack:
                    return False
                stack.pop()
        return not stack  # 如果栈为空，说明括号成对
    used_indices = set()
    results = []
    types = []
    judge_1 = 1
    judge_2 = 1
    i = 0
    while i < len(dot_bracket_str):
        
        if dot_bracket_str[i] == '[' and i not in used_indices:
            start_index = i
            end_index = None
            stack = ['[']
            is_continuous = True

            for k in range(i + 1, len(dot_bracket_str)):
                if k in used_indices:
                    continue
                if dot_bracket_str[k] == '[':
                    if dot_bracket_str[k - 1] != '[':
                        is_continuous = False
                        i = k - 1
                        break
                    stack.append('[')
                elif dot_bracket_str[k] == ']':
                    if stack and stack[-1] == '[':
                        stack.pop()
                        if not stack:
                            end_index = k
                            break

            if not is_continuous or end_index is None:
                used_indices.add(i)
                i += 1
                continue

            subsequence = dot_bracket_str[start_index:end_index + 1]
            
            stackrb = []
            stacklb = []
            start_index_extend = start_index
            end_index_extend = end_index
            ii = 0
            type_k = []
            while ii < len(subsequence):
                if subsequence[ii] == '(':
                    stacklb.append('(')
                    type_k.append(-1)
                elif subsequence[ii] == ')':
                    if stacklb:
                        stacklb.pop()
                        type_k.remove(-1)
                        type_k.append(0)
                    else:
                        stackrb.append(')')
                        type_k.append(1)
                ii += 1

            # Checking and extending on right
            if stacklb:
                judge_1 = 0
                for jj in range(end_index + 1, len(dot_bracket_str)):
                    if dot_bracket_str[jj] == '(':
                        stacklb.append('(')
                    elif dot_bracket_str[jj] == ')':
                        stacklb.pop()
                        if not stacklb:
                            end_index_extend = jj
                            judge_1 = 1
                            break
                    elif dot_bracket_str[jj] == '[' or dot_bracket_str[jj] == ']':
                        judge_1 = 0
                        break

            # Checking and extending on left
            if stackrb:
                judge_2 = 0
                for kk in range(start_index - 1, -1, -1):
                    if dot_bracket_str[kk] == ')':
                        stackrb.append(')')
                    elif dot_bracket_str[kk] == '(':
                        stackrb.pop()
                        if not stackrb:
                            start_index_extend = kk
                            judge_2 = 1
                            break
                    elif dot_bracket_str[kk] == '[' or dot_bracket_str[kk] == ']':
                        judge_2 = 0
                        break
                    else:
                        continue
            extended_subsequence = []
            if judge_1 == 1 or judge_2 == 1:
                extended_subsequence = dot_bracket_str[start_index_extend:end_index_extend + 1]

            if extended_subsequence and is_balanced(extended_subsequence):
                results.append(extended_subsequence)
                type_k1 = type_kp(type_k)
                types.append(type_k1)
            used_indices.update(range(start_index_extend, end_index_extend + 1))
        i += 1
    return results, types

def find_kps(G, key):
    sequences_and_pairs = graph_to_sequence_and_pairs(G)
    sequences_dot_brackets = []
    kps_nodes = {}

    for sequence, pairs in sequences_and_pairs:
        seq_list, dot_bracket_str = convert_to_dot_bracket(sequence, pairs)
        sequences_dot_brackets.append((seq_list, dot_bracket_str))

    kp_index = 1
    for seq_list, dot_bracket_str in sequences_dot_brackets:
        kp_subsequences, types = find_all_subsequences(dot_bracket_str)
        
        for kp, kp_type in zip(kp_subsequences, types):
            if kp:
                start_index = dot_bracket_str.index(kp)
                end_index = start_index + len(kp) - 1
                kp_nodes = seq_list[start_index:end_index + 1]
                kps_nodes[f"G{kp_index}"] = {
                    "nodes": kp_nodes,
                    "kp_type": f"pseudoknot_{kp_type}",
                    "dot_bracket": kp
                }
                kp_index += 1
    return kps_nodes

# Function to create and populate the database
def create_and_populate_db(data):
    conn = sqlite3.connect("rna_kps_11.db")
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS data (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        motif_type TEXT,
                        pdbid TEXT,
                        paired_nt_number TEXT,
                        nt_number TEXT,
                        dot_bracket TEXT)''')

    for entry in data:
        for pdbid, kps in entry.items():
            for kp_key, kp_data in kps.items():
                nodes = kp_data["nodes"]
                paired_nt_number = f"{nodes[0]},{nodes[-1]}" if nodes else ""
                nt_number = ','.join(nodes)
                motif_type = kp_data["kp_type"]
                dot_bracket = kp_data["dot_bracket"]
                pdbid_r = pdbid[3:7]
                if "." not in nt_number:
                    cursor.execute('''INSERT INTO data (motif_type, pdbid, paired_nt_number, nt_number, dot_bracket)
                                  VALUES (?, ?, ?, ?, ?)''', (motif_type, pdbid_r, paired_nt_number, nt_number, dot_bracket))

    conn.commit()
    conn.close()

with open("D:/RNA_design_py/data/saved_graphs_canonical20241010.pickle", "rb") as f:
    graph_data = pickle.load(f)
import os
pdb_folder = r"D:\RNA_design_py\data\PDBRNA"
DNA_complex = []
for filename in os.listdir(pdb_folder):
    if filename.endswith(".ent.pdb"):  # Check if the file has the correct extension
        filepath = os.path.join(pdb_folder, filename)
        with open(filepath, 'r') as file:
            content = file.read()
        if any(nucleotide in content for nucleotide in [" DA ", " DU ", " DG ", " DC "]):
            # Generate the key and append it to the DNA_complex list
            key = filename.replace("ent.pdb", "ent_output.txt")
            DNA_complex.append(key)
# print("DNA Complex Keys:", DNA_complex)
print("Number of DNA Complexes:", len(DNA_complex))
kps_graphs = []
count = 1
for key, graph in graph_data.items():
    if count % 100 == 0:
        print(f"finish {count} RNA searching")
    count += 1
    if count >= 4000:
        break
    if key not in DNA_complex:# and key =="pdb1ffk.ent_output.txt": 
        kps_nodes = find_kps(graph, key)
        if kps_nodes:
            kps_graphs.append({key: kps_nodes})

create_and_populate_db(kps_graphs)

print("Database creation and population complete.")
