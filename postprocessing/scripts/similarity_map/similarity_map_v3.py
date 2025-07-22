import sqlite3
import json
import itertools
from collections import defaultdict
import numpy as np
from scipy.optimize import linear_sum_assignment
import networkx as nx
import time
from grakel import Graph
from grakel.kernels import WeisfeilerLehman, EdgeHistogram, VertexHistogram

def compute_motif_similarity_4(motif1, motif2):
    """
    Compute S_mn using Weisfeiler-Lehman Graph Kernel.
    """
    G1 = build_motif_graph(motif1)
    G2 = build_motif_graph(motif2)
    
    # Convert to Grakel format
    g1 = Graph(nx.adjacency_matrix(G1).todense().tolist())
    g2 = Graph(nx.adjacency_matrix(G2).todense().tolist())

    # Compute WL kernel similarity
    wl_kernel = WeisfeilerLehman(n_iter=3, base_kernel=VertexHistogram)
    similarity_matrix = wl_kernel.fit_transform([g1, g2])
    
    return similarity_matrix[0, 1]  # Similarity score in [0,1]


# -------------------------
# Database Loading
# -------------------------
def load_data(db_path):
    """Load RNA motif data from SQLite database"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM data")
    rows = cursor.fetchall()
    columns = [desc[0] for desc in cursor.description]
    data = [dict(zip(columns, row)) for row in rows]
    conn.close()
    return data

# -------------------------
# Group Motifs by RNA
# -------------------------
def group_motifs_by_RNA(data):
    """Group all motifs (all types) by RNA (using pdbid as key)"""
    RNA_motifs = defaultdict(list)
    for motif in data:
        RNA_motifs[motif['pdbid'].lower()].append(motif)
    return RNA_motifs

# -------------------------
# Motif Weight Calculation
# -------------------------
def compute_motif_weight(motif):
    """Compute weight as: w(m) = (# nucleotides) * (1 + # noncanonical bp)"""
    if not motif.get('nt_number'):
        return 0
    nucleotides = motif['nt_number'].split(',')
    L = len(nucleotides)
    NC = 0
    try:
        edges = json.loads(motif['non_graph_edges'])
    except Exception:
        edges = []
    for edge in edges:
        if isinstance(edge[2], dict) and edge[2].get("attribute") == [0, 1, 0]:
            NC += 1
    return L * (1 + NC)

# -------------------------
# Build Motif Graph
# -------------------------
def build_motif_graph(motif):
    """Build a graph from the motif's non_graph_edges (all bond types)"""
    G = nx.Graph()
    if not motif.get('nt_number'):
        return G
    nucleotides = [nt.strip() for nt in motif['nt_number'].split(',')]
    G.add_nodes_from(nucleotides)
    try:
        edges = json.loads(motif['non_graph_edges'])
        for edge in edges:
            node1 = str(edge[0]).strip()
            node2 = str(edge[1]).strip()
            # Store edge attributes
            attr = edge[2].get("attribute", []) if isinstance(edge[2], dict) else []
            if node1 in G.nodes and node2 in G.nodes:
                G.add_edge(node1, node2, attr=attr)
    except Exception:
        pass
    return G

# -------------------------
# Motif Similarity Functions
# -------------------------
def compute_motif_similarity_1(motif1, motif2):
    """
    Compute S_mn = 1 - (D_mn / D_max) using graph_edit_distance.
    D_max = (N1 + N2) + (E1 + E2) with unit costs.
    """
    G1 = build_motif_graph(motif1)
    G2 = build_motif_graph(motif2)
    N1, E1 = G1.number_of_nodes(), G1.number_of_edges()
    N2, E2 = G2.number_of_nodes(), G2.number_of_edges()
    D_max = (N1 + N2) + (E1 + E2)
    try:
        D_mn = nx.graph_edit_distance(G1, G2, timeout=1)
    except Exception:
        D_mn = D_max / 2.0
    if D_mn is None:
        D_mn = D_max / 2.0
    return 1 - (D_mn / D_max) if D_max > 0 else 0.0 

def compute_motif_similarity_2(motif1, motif2):
    """Heuristic: Return 1 if non_graph_edges strings match exactly, else 0."""
    return 1.0 if motif1['non_graph_edges'] == motif2['non_graph_edges'] else 0.0

def compute_motif_similarity_3(motif1, motif2):
    """
    Lightweight heuristic comparing basic features:
    S = 1 - (|L1 - L2| + |E1 - E2|) / (L1 + L2 + E1 + E2)
    """
    L1 = len(motif1['nt_number'].split(',')) if motif1.get('nt_number') else 0
    L2 = len(motif2['nt_number'].split(',')) if motif2.get('nt_number') else 0
    try:
        E1 = len(json.loads(motif1['non_graph_edges']))
    except:
        E1 = 0
    try:
        E2 = len(json.loads(motif2['non_graph_edges']))
    except:
        E2 = 0
    numerator = abs(L1 - L2) + abs(E1 - E2)
    denominator = (L1 + L2 + E1 + E2)
    return 1 - (numerator / denominator) if denominator > 0 else 0.0

# -------------------------
# New Hybrid Graph Similarity Function
# -------------------------
def compute_hybrid_graph_similarity(motif1, motif2):
    """
    Compute graph similarity using a hybrid approach that considers:
    1. Structural similarity (based on graph properties)
    2. Edge attribute similarity
    """
    G1 = build_motif_graph(motif1)
    G2 = build_motif_graph(motif2)
    
    # If either graph is empty, return appropriate similarity
    if G1.number_of_nodes() == 0 and G2.number_of_nodes() == 0:
        return 1.0
    if G1.number_of_nodes() == 0 or G2.number_of_nodes() == 0:
        return 0.0
    
    # Check for trivial cases first
    if e1 == 0 and e2 == 0:
        return 1.0  # Both have no edges, so they are identical (trivial case)
    elif e1 == 0 or e2 == 0:
        return 0.0  # One has edges, one doesn't - they're fundamentally different
    
    # For RNA motifs, there are two key aspects to consider:
    # 1. The structural topology (connections between nucleotides)
    # 2. The base-pairing types (edge attributes)
    
    # Both graphs have same number of nodes, check for structural equivalence
    if nx.is_isomorphic(G1, G2):
        # Perfect structural match
        structure_sim = 1.0
    else:
        # No perfect structural match, compute approximate similarity using Jaccard similarity
        # of node-pair sets (which nucleotides are connected)
        edge_set1 = {frozenset((str(e[0]), str(e[1]))) for e in G1.edges()}
        edge_set2 = {frozenset((str(e[0]), str(e[1]))) for e in G2.edges()}
        
        # Jaccard similarity: |intersection| / |union|
        intersection = len(edge_set1.intersection(edge_set2))
        union = len(edge_set1.union(edge_set2))
        structure_sim = intersection / union if union > 0 else 0.0
    
    # For edge attributes, we are specifically interested in RNA bond types
    # Count edge attributes in each graph
    attr_count1 = {}
    attr_count2 = {}
    
    for _, _, attr in G1.edges(data=True):
        attr_tuple = tuple(attr.get('attr', []))
        attr_count1[attr_tuple] = attr_count1.get(attr_tuple, 0) + 1
    
    for _, _, attr in G2.edges(data=True):
        attr_tuple = tuple(attr.get('attr', []))
        attr_count2[attr_tuple] = attr_count2.get(attr_tuple, 0) + 1
    
    # Compute Cosine similarity between attribute count vectors
    # This measure is invariant to the absolute counts and focuses on the proportions
    all_attrs = set(attr_count1.keys()).union(attr_count2.keys())
    
    # Build vectors from the attribute counts
    vec1 = [attr_count1.get(attr, 0) for attr in all_attrs]
    vec2 = [attr_count2.get(attr, 0) for attr in all_attrs]
    
    # Compute cosine similarity: dot(vec1, vec2) / (||vec1|| * ||vec2||)
    dot_product = sum(a * b for a, b in zip(vec1, vec2))
    norm1 = sum(a * a for a in vec1) ** 0.5
    norm2 = sum(b * b for b in vec2) ** 0.5
    attr_sim = dot_product / (norm1 * norm2) if norm1 > 0 and norm2 > 0 else 0.0
    
    # Final similarity is a geometric mean of structure and attribute similarity
    # This represents that both aspects are equally important and multiplicative
    # If either is 0, the result is 0 (can't have equivalent RNA motifs with completely
    # different structures or completely different bond types)
    return (structure_sim * attr_sim) ** 0.5

# -------------------------
# Graph Similarity Function using Adjacency Matrix
# -------------------------
def compute_spectral_similarity(motif1, motif2):
    """
    Compute graph similarity using eigendecomposition of adjacency matrices.
    This is a well-established method that captures structural similarity.
    
    The similarity is based on comparing the set of eigenvalues of the adjacency matrices,
    which characterize the graph structure.
    """
    G1 = build_motif_graph(motif1)
    G2 = build_motif_graph(motif2)
    
    # Handle trivial cases
    if G1.number_of_nodes() == 0 and G2.number_of_nodes() == 0:
        return 1.0
    if G1.number_of_nodes() == 0 or G2.number_of_nodes() == 0:
        return 0.0
    
    # Get adjacency matrices
    A1 = nx.adjacency_matrix(G1).todense()
    A2 = nx.adjacency_matrix(G2).todense()
    
    # Check for exact graph isomorphism first (fastest path for identical graphs)
    # This handles the case where the graphs are structurally identical
    if A1.shape == A2.shape and np.all(A1 == A2):
        return 1.0
    
    # Now compute eigenvalues
    try:
        eig1 = np.linalg.eigvals(A1)
        eig2 = np.linalg.eigvals(A2)
        
        # Sort the eigenvalues
        eig1 = sorted(abs(eig1))
        eig2 = sorted(abs(eig2))
        
        # Pad the smaller set with zeros to make lengths equal
        if len(eig1) < len(eig2):
            eig1 = eig1 + [0] * (len(eig2) - len(eig1))
        elif len(eig2) < len(eig1):
            eig2 = eig2 + [0] * (len(eig1) - len(eig2))
        
        # Calculate similarity based on eigenvalue difference
        # This is a standard approach in spectral graph theory
        diff_sum = sum((e1 - e2)**2 for e1, e2 in zip(eig1, eig2))
        max_sum = sum(max(e1**2, e2**2) for e1, e2 in zip(eig1, eig2))
        
        if max_sum == 0:
            return 1.0
            
        similarity = 1 - (diff_sum / max_sum)**0.5
        return max(0, min(1, similarity))  # Ensure in range [0,1]
        
    except Exception as e:
        # Fallback to a simpler size-based similarity metric
        n1, e1 = G1.number_of_nodes(), G1.number_of_edges()
        n2, e2 = G2.number_of_nodes(), G2.number_of_edges()
        
        if n1 + n2 + e1 + e2 == 0:
            return 1.0
            
        return 1 - (abs(n1-n2) + abs(e1-e2)) / (n1 + n2 + e1 + e2)

def enhanced_wl_kernel_similarity(motif1, motif2):
    """
    Use GraKeL's graph kernels to compare graphs with edge attributes.
    This is an established algorithm that properly handles both structure and edge attributes.
    """
    G1 = build_motif_graph(motif1)
    G2 = build_motif_graph(motif2)
    
    # Handle trivial cases
    if G1.number_of_nodes() == 0 and G2.number_of_nodes() == 0:
        return 1.0
    if G1.number_of_nodes() == 0 or G2.number_of_nodes() == 0:
        return 0.0
    
    try:
        wl_kernel = WeisfeilerLehman(n_iter=3, normalize=True)
        edge_kernel = EdgeHistogram(normalize=True)
        
        # Properly convert to GraKeL format
        node_labels1 = {node: G1.degree(node) for node in G1.nodes()}
        node_labels2 = {node: G2.degree(node) for node in G2.nodes()}
        
        edge_labels1 = {}
        for u, v, attr in G1.edges(data=True):
            edge_key = (u, v) if u < v else (v, u)
            attr_tuple = tuple(attr.get('attr', []))
            edge_labels1[edge_key] = hash(str(attr_tuple)) % 1000000
        
        edge_labels2 = {}
        for u, v, attr in G2.edges(data=True):
            edge_key = (u, v) if u < v else (v, u)
            attr_tuple = tuple(attr.get('attr', []))
            edge_labels2[edge_key] = hash(str(attr_tuple)) % 1000000
        
        # Create GraKeL graph objects
        g1 = Graph(G1.edges(), node_labels=node_labels1, edge_labels=edge_labels1)
        g2 = Graph(G2.edges(), node_labels=node_labels2, edge_labels=edge_labels2)
        
        # Calculate similarity using properly formatted graph objects
        edge_similarity = edge_kernel.fit_transform([g1, g2])[0, 1]
        
        structure_similarity = wl_kernel.fit_transform([g1, g2])[0, 1]
        
        # Combine both metrics (geometric mean ensures both are important)
        similarity = (edge_similarity * structure_similarity) ** 0.5
        
        return similarity
    
    except Exception as e:
        print(f"Warning: Graph kernel error - {str(e)}. Falling back to basic similarity.")
        # Fall back to a simpler similarity measure
        if nx.is_isomorphic(G1, G2):
            return 1.0
        else:
            return compute_motif_similarity_3(motif1, motif2)  # Use the lightweight heuristic as fallback

# Choose the similarity function here:
similarity_function = enhanced_wl_kernel_similarity

# -------------------------
# Optimal Matching using Hungarian Algorithm
# -------------------------
def optimal_matching(motifs_A, motifs_B):
    n = len(motifs_A)
    m = len(motifs_B)
    sim_matrix = np.zeros((n, m))
    weight_matrix = np.zeros((n, m))
    
    for i in range(n):
        for j in range(m):
            sim = similarity_function(motifs_A[i], motifs_B[j])
            sim_matrix[i, j] = sim
            w_i = compute_motif_weight(motifs_A[i])
            w_j = compute_motif_weight(motifs_B[j])
            weight_matrix[i, j] = (w_i + w_j) / 2.0
            
    cost_matrix = 1 - sim_matrix
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    total_weighted_similarity = 0.0
    total_matched_weight = 0.0
    for i, j in zip(row_ind, col_ind):
        total_weighted_similarity += weight_matrix[i, j] * sim_matrix[i, j]
        total_matched_weight += weight_matrix[i, j]
    return total_weighted_similarity, total_matched_weight

# -------------------------
# Overall RNA Similarity Aggregation
# -------------------------
def compute_RNA_similarity(RNA1_motifs, RNA2_motifs):
    if not RNA1_motifs or not RNA2_motifs:
        return 0.0
    total_weight_RNA1 = sum(compute_motif_weight(m) for m in RNA1_motifs)
    total_weight_RNA2 = sum(compute_motif_weight(m) for m in RNA2_motifs)
    matched_sim, matched_weight = optimal_matching(RNA1_motifs, RNA2_motifs)
    S_ij = (2 * matched_sim) / (total_weight_RNA1 + total_weight_RNA2)
    return S_ij

# -------------------------
# Compute Similarity Map (Serial Version)
# -------------------------
def compute_similarity_map(RNA_motifs):
    RNA_ids = list(RNA_motifs.keys())
    similarity_map = {}
    pairs = [(i, j) for i in range(len(RNA_ids)) for j in range(i+1, len(RNA_ids))]
    total_pairs = len(pairs)
    processed = 0
    start_time = time.time()
    
    for pair in pairs:
        i, j = pair
        rna1 = RNA_ids[i]
        rna2 = RNA_ids[j]
        sim = compute_RNA_similarity(RNA_motifs[rna1], RNA_motifs[rna2])
        if sim > 0:  # Only save non-zero similarities
            similarity_map[(rna1, rna2)] = sim
        
        processed += 1
        if processed % 100 == 0:
            elapsed = time.time() - start_time
            print(f"Processed {processed}/{total_pairs} pairs... Elapsed time: {elapsed:.2f}s, Current similarity: {sim:.2f}")
    
    return similarity_map

# -------------------------
# Save Similarity Map to File
# -------------------------
def save_similarity_map(similarity_map, filename="similarity_results_structure_attr.txt"):
    with open(filename, "w") as f:
        for (rna1, rna2), sim in similarity_map.items():
            f.write(f"Similarity between {rna1} and {rna2}: {sim:.2f}\n")
    print(f"Full similarity map with non-zero scores saved to {filename}")

# -------------------------
# Main Pipeline
# -------------------------
def main():
    print("Starting RNA similarity analysis (serial version)...")
    db_path = "../../../data/ATLAS2025_TEST.db"
    data = load_data(db_path)
    RNA_motifs = group_motifs_by_RNA(data)
    print(f"Total RNAs: {len(RNA_motifs)}")
    
    print("Calculating full RNA similarity map (serial version)...")
    similarity_map = compute_similarity_map(RNA_motifs)
    save_similarity_map(similarity_map, filename="../../../postprocessing/results2025/similarity_map/similarity_results_ATLAS2025.txt")

if __name__ == "__main__":
    main()
