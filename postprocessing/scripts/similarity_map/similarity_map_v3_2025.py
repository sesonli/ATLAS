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
        # 创建基于边属性的核函数
        base_kernel = VertexHistogram(normalize=True)
        wl_kernel = WeisfeilerLehman(n_iter=3, base_kernel=base_kernel, normalize=True)
        
        # 为节点创建标签（使用度数）
        node_labels1 = {node: str(G1.degree(node)) for node in G1.nodes()}
        node_labels2 = {node: str(G2.degree(node)) for node in G2.nodes()}
        
        # 为边创建标签（将RNA键类型转换为字符串标签）
        edge_labels1 = {}
        for u, v, attr in G1.edges(data=True):
            edge_key = (u, v) if u <= v else (v, u)  # 确保一致的键顺序
            attr_list = attr.get('attr', [])
            
            # 将RNA键类型数组转换为字符串标签
            if len(attr_list) >= 3:
                # [0,0,1] = covalent, [1,0,0] = watson_crick, [0,1,0] = non_watson_crick
                if attr_list == [0, 0, 1]:
                    edge_label = "covalent"
                elif attr_list == [1, 0, 0]:
                    edge_label = "watson_crick"
                elif attr_list == [0, 1, 0]:
                    edge_label = "non_watson_crick"
                else:
                    edge_label = f"custom_{attr_list[0]}_{attr_list[1]}_{attr_list[2]}"
            else:
                edge_label = "unknown"
            
            edge_labels1[edge_key] = edge_label
        
        edge_labels2 = {}
        for u, v, attr in G2.edges(data=True):
            edge_key = (u, v) if u <= v else (v, u)  # 确保一致的键顺序
            attr_list = attr.get('attr', [])
            
            # 将RNA键类型数组转换为字符串标签
            if len(attr_list) >= 3:
                # [0,0,1] = covalent, [1,0,0] = watson_crick, [0,1,0] = non_watson_crick
                if attr_list == [0, 0, 1]:
                    edge_label = "covalent"
                elif attr_list == [1, 0, 0]:
                    edge_label = "watson_crick"
                elif attr_list == [0, 1, 0]:
                    edge_label = "non_watson_crick"
                else:
                    edge_label = f"custom_{attr_list[0]}_{attr_list[1]}_{attr_list[2]}"
            else:
                edge_label = "unknown"
            
            edge_labels2[edge_key] = edge_label
        
        # 创建GraKeL图对象，确保提供正确的格式
        try:
            # 转换为边列表格式
            edge_list1 = [(u, v) for u, v in G1.edges()]
            edge_list2 = [(u, v) for u, v in G2.edges()]
            
            # 创建图对象，包含节点和边标签
            g1 = Graph(edge_list1, node_labels=node_labels1, edge_labels=edge_labels1)
            g2 = Graph(edge_list2, node_labels=node_labels2, edge_labels=edge_labels2)
            
            # 计算相似性
            similarity_matrix = wl_kernel.fit_transform([g1, g2])
            similarity = similarity_matrix[0, 1]
            
            return max(0.0, min(1.0, similarity))  # 确保在[0,1]范围内
            
        except Exception as kernel_error:
            print(f"Warning: GraKeL kernel calculation failed - {str(kernel_error)}. Using edge histogram fallback.")
            
            # 使用边直方图作为后备方案
            try:
                edge_kernel = EdgeHistogram(normalize=True)
                edge_similarity = edge_kernel.fit_transform([g1, g2])[0, 1]
                return max(0.0, min(1.0, edge_similarity))
            except:
                # 如果所有图核都失败，使用基本相似性
                print(f"Warning: All graph kernels failed. Using basic structural similarity.")
                return compute_motif_similarity_3(motif1, motif2)
    
    except Exception as e:
        print(f"Warning: Graph kernel error - {str(e)}. Falling back to basic similarity.")
        # 回退到简单的相似性度量
        return compute_motif_similarity_3(motif1, motif2)  # 使用轻量级启发式作为回退

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
# 并行分块处理测试版本
# -------------------------
import multiprocessing as mp
import os
import shutil
import math
from functools import partial

def create_test_database(source_db, test_db, sample_rnas=100):
    """从原数据库抽取部分数据创建测试数据库"""
    print(f"创建测试数据库：从{sample_rnas}个RNA中抽取数据...")
    
    # 连接原数据库
    source_conn = sqlite3.connect(source_db)
    source_cursor = source_conn.cursor()
    
    # 获取指定数量的RNA列表
    source_cursor.execute("SELECT DISTINCT pdbid FROM data ORDER BY RANDOM() LIMIT ?", (sample_rnas,))
    test_rnas = [row[0] for row in source_cursor.fetchall()]
    
    print(f"选择的RNA数量: {len(test_rnas)}")
    
    # 获取原表的列信息
    source_cursor.execute("PRAGMA table_info(data)")
    columns_info = source_cursor.fetchall()
    print(f"原表结构：{len(columns_info)}列")
    for col in columns_info:
        print(f"  列{col[0]}: {col[1]} ({col[2]})")
    
    # 获取这些RNA的所有motif数据，明确指定列名
    column_names = [col[1] for col in columns_info]  # 获取列名
    columns_str = ", ".join(column_names)
    
    placeholders = ','.join(['?' for _ in test_rnas])
    query = f"SELECT {columns_str} FROM data WHERE pdbid IN ({placeholders})"
    print(f"查询SQL: {query}")
    
    source_cursor.execute(query, test_rnas)
    all_rows = source_cursor.fetchall()
    
    print(f"提取的motif数量: {len(all_rows)}")
    if all_rows:
        print(f"第一行数据长度: {len(all_rows[0])}")
        print(f"第一行数据类型: {[type(x).__name__ for x in all_rows[0]]}")
    
    # 创建测试数据库
    test_conn = sqlite3.connect(test_db)
    test_cursor = test_conn.cursor()
    
    # 手动创建表结构
    create_sql = "CREATE TABLE data ("
    column_definitions = []
    for col_info in columns_info:
        col_name = col_info[1]  # 列名
        col_type = col_info[2]  # 数据类型
        column_definitions.append(f"{col_name} {col_type}")
    create_sql += ", ".join(column_definitions) + ")"
    
    print(f"创建表SQL: {create_sql}")
    test_cursor.execute(create_sql)
    
    # 插入数据，使用正确的列数
    placeholders_insert = ','.join(['?' for _ in column_names])
    insert_sql = f"INSERT INTO data ({columns_str}) VALUES ({placeholders_insert})"
    
    print(f"插入SQL: {insert_sql}")
    print(f"期望列数: {len(column_names)}, 实际数据列数: {len(all_rows[0]) if all_rows else 0}")
    
    # 验证数据完整性
    if all_rows and len(all_rows[0]) != len(column_names):
        print(f"❌ 数据列数不匹配！期望{len(column_names)}列，实际{len(all_rows[0])}列")
        # 尝试重新查询，使用*选择所有列
        source_cursor.execute(f"SELECT * FROM data WHERE pdbid IN ({placeholders}) LIMIT 1", test_rnas)
        sample_row = source_cursor.fetchone()
        print(f"使用SELECT *的样本行长度: {len(sample_row) if sample_row else 0}")
        
        # 重新获取完整数据
        source_cursor.execute(f"SELECT * FROM data WHERE pdbid IN ({placeholders})", test_rnas)
        all_rows = source_cursor.fetchall()
        print(f"重新获取的数据行长度: {len(all_rows[0]) if all_rows else 0}")
    
    # 插入数据
    test_cursor.executemany(insert_sql, all_rows)
    
    test_conn.commit()
    source_conn.close()
    test_conn.close()
    
    print(f"测试数据库创建完成：{test_db}")
    print(f"包含 {len(test_rnas)} 个RNA，{len(all_rows)} 个motif")
    return test_rnas

def load_lightweight_index(db_path):
    """轻量级加载：只读取RNA-motif索引信息"""
    print("构建轻量级索引...")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # 读取pdbid，保持原始大小写
    cursor.execute("SELECT DISTINCT pdbid FROM data")
    rna_list = [row[0] for row in cursor.fetchall()]  # 移除.lower()
    
    conn.close()
    
    print(f"索引构建完成：{len(rna_list)} 个RNA")
    return rna_list

def create_rna_blocks(rna_list, block_size):
    """将RNA列表分块"""
    blocks = [rna_list[i:i+block_size] for i in range(0, len(rna_list), block_size)]
    print(f"分块完成：{len(rna_list)} 个RNA分成 {len(blocks)} 块，每块约 {block_size} 个")
    return blocks

def load_block_motifs(db_path, target_rna_list):
    """按需加载：只加载指定RNA的motif数据"""
    if not target_rna_list:
        return []
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # 使用IN查询只获取需要的RNA的motif
    placeholders = ','.join(['?' for _ in target_rna_list])
    query = f"SELECT * FROM data WHERE pdbid IN ({placeholders})"
    cursor.execute(query, target_rna_list)
    
    rows = cursor.fetchall()
    columns = [desc[0] for desc in cursor.description]
    data = [dict(zip(columns, row)) for row in rows]
    
    print(f"加载了 {len(data)} 个motif，来自 {len(target_rna_list)} 个RNA")
    
    conn.close()
    return data

def compute_block_pairs_similarity(block1_rnas, block2_rnas, all_motif_data):
    """计算两个块之间的RNA相似性"""
    # 使用现有函数进行分组
    RNA_motifs = group_motifs_by_RNA(all_motif_data)
    
    results = {}
    total_pairs = 0
    non_zero_pairs = 0
    
    for rna1 in block1_rnas:
        for rna2 in block2_rnas:
            # 避免重复计算（只在块内时需要）
            if block1_rnas == block2_rnas and rna1 >= rna2:
                continue
            
            total_pairs += 1
            
            # 统一使用小写进行匹配（与group_motifs_by_RNA保持一致）
            rna1_key = rna1.lower()
            rna2_key = rna2.lower()
            
            # 检查RNA是否存在于数据中
            if rna1_key not in RNA_motifs or rna2_key not in RNA_motifs:
                print(f"  跳过：{rna1}({rna1_key}) 或 {rna2}({rna2_key}) 不在数据中")
                continue
                
            # 直接调用现有函数，零修改
            sim = compute_RNA_similarity(RNA_motifs[rna1_key], RNA_motifs[rna2_key])
            
            if sim > 0:
                results[(rna1, rna2)] = sim
                non_zero_pairs += 1
                print(f"  相似性：{rna1} vs {rna2} = {sim:.3f}")
            else:
                print(f"  零相似性：{rna1} vs {rna2} (motif数: {len(RNA_motifs[rna1_key])}, {len(RNA_motifs[rna2_key])})")
    
    print(f"  块对统计：总计{total_pairs}对，非零{non_zero_pairs}对，保存{len(results)}对")
    return results

def distribute_block_pairs(total_blocks, num_processes):
    """将块对分配给多个进程"""
    # 生成所有块对组合
    all_pairs = []
    for i in range(total_blocks):
        for j in range(i, total_blocks):
            all_pairs.append((i, j))
    
    # 平均分配给进程
    pairs_per_process = len(all_pairs) // num_processes
    process_tasks = []
    
    for p in range(num_processes):
        start_idx = p * pairs_per_process
        if p == num_processes - 1:  # 最后一个进程处理剩余的
            end_idx = len(all_pairs)
        else:
            end_idx = (p + 1) * pairs_per_process
        
        process_tasks.append(all_pairs[start_idx:end_idx])
    
    return process_tasks

def process_worker_with_file_output(assigned_pairs, db_path, rna_blocks, output_file, worker_id):
    """Worker进程：计算并直接写文件"""
    print(f"Worker {worker_id}: 开始处理 {len(assigned_pairs)} 个块对...")
    
    processed_pairs = 0
    total_similarities = 0
    
    with open(output_file, "w") as f:
        for pair_idx, (block_i, block_j) in enumerate(assigned_pairs):
            print(f"Worker {worker_id}: 处理块对 ({block_i},{block_j}) - {pair_idx+1}/{len(assigned_pairs)}")
            
            # 确定需要加载的RNA
            target_rnas = list(set(rna_blocks[block_i] + rna_blocks[block_j]))
            
            # 加载当前块对的数据
            block_data = load_block_motifs(db_path, target_rnas)
            
            if not block_data:
                print(f"Worker {worker_id}: 块对 ({block_i},{block_j}) 无数据，跳过")
                continue
            
            # 计算相似性
            results = compute_block_pairs_similarity(rna_blocks[block_i], rna_blocks[block_j], block_data)
            
            # 直接写入文件
            for (rna1, rna2), sim in results.items():
                f.write(f"Similarity between {rna1} and {rna2}: {sim:.2f}\n")
            
            total_similarities += len(results)
            processed_pairs += 1
            
            # 实时flush，确保数据写入
            f.flush()
            
            # 释放内存
            del block_data
    
    print(f"Worker {worker_id}: 完成！处理了 {processed_pairs} 个块对，计算了 {total_similarities} 个相似性")
    return total_similarities

def merge_result_files(temp_files, final_output):
    """合并所有worker的结果文件"""
    print("合并结果文件...")
    total_lines = 0
    
    with open(final_output, "w") as outf:
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                with open(temp_file, "r") as inf:
                    lines = inf.readlines()
                    outf.writelines(lines)
                    total_lines += len(lines)
                os.remove(temp_file)  # 删除临时文件
    
    print(f"合并完成：{len(temp_files)} 个文件，总计 {total_lines} 行结果")
    return total_lines

def parallel_blockwise_similarity_map(db_path, block_size=10, num_processes=2, output_file="test_results.txt"):
    """并行分块计算相似性图谱（策略1：分别写文件）"""
    print("🚀 启动并行分块相似性计算...")
    start_time = time.time()
    
    # 第1步：轻量级索引
    rna_list = load_lightweight_index(db_path)
    rna_blocks = create_rna_blocks(rna_list, block_size)
    
    total_blocks = len(rna_blocks)
    total_pairs = total_blocks * (total_blocks + 1) // 2
    print(f"配置：{len(rna_list)} 个RNA，{total_blocks} 块，{num_processes} 进程")
    print(f"预计计算 {total_pairs} 个块对")
    
    # 第2步：任务分配
    process_tasks = distribute_block_pairs(total_blocks, num_processes)
    
    for i, task in enumerate(process_tasks):
        print(f"进程 {i}: 分配到 {len(task)} 个块对")
    
    # 第3步：创建临时目录
    temp_dir = "temp_results"
    os.makedirs(temp_dir, exist_ok=True)
    
    # 第4步：并行计算
    print("启动并行计算...")
    
    try:
        with mp.Pool(processes=num_processes) as pool:
            async_results = []
            temp_files = []
            
            for i, assigned_pairs in enumerate(process_tasks):
                if not assigned_pairs:  # 跳过空任务
                    continue
                    
                # 每个worker有自己的输出文件
                worker_output = f"{temp_dir}/worker_{i}_results.txt"
                temp_files.append(worker_output)
                
                result = pool.apply_async(
                    process_worker_with_file_output,
                    (assigned_pairs, db_path, rna_blocks, worker_output, i)
                )
                async_results.append(result)
            
            # 等待所有worker完成
            total_similarities = 0
            for async_result in async_results:
                worker_similarities = async_result.get()  # 等待完成
                total_similarities += worker_similarities
        
        # 第5步：合并结果文件
        final_lines = merge_result_files(temp_files, output_file)
        
        # 清理临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        
        # 统计结果
        elapsed = time.time() - start_time
        print(f"✅ 并行计算完成！")
        print(f"⏱️  总耗时：{elapsed:.1f}秒")
        print(f"📈 计算了 {total_similarities} 个相似性对")
        print(f"📁 结果保存到：{output_file}")
        
        return total_similarities
        
    except Exception as e:
        print(f"❌ 并行计算出错：{e}")
        # 清理临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        raise

def quick_test_main():
    """快速测试主程序：100个RNA，2个worker，几分钟完成"""
    print("🧪 快速测试程序启动...")
    print("配置：100个RNA，10块，2个worker，预期2-5分钟完成")
    
    start_time = time.time()
    
    # 配置参数 - 修复：使用完整的ATLAS数据库
    source_db = "../../../data/ATLAS.db"  # 使用完整的正常数据库
    test_db = "../../../data/ATLAS2025_QUICK_TEST_FIXED.db"  # 创建新的测试数据库
    output_file = "../../../postprocessing/results2025/similarity_map/quick_test_results.txt"
    
    try:
        # 第1步：创建小测试数据库（如果不存在）
        print("📊 第1步：检查测试数据库...")
        if os.path.exists(test_db):
            print(f"测试数据库已存在：{test_db}")
            # 读取现有数据库的RNA数量
            conn = sqlite3.connect(test_db)
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(DISTINCT pdbid) FROM data")
            rna_count = cursor.fetchone()[0]
            conn.close()
            print(f"现有数据库包含 {rna_count} 个RNA")
        else:
            test_rnas = create_test_database(source_db, test_db, sample_rnas=100)
        
        # 第2步：并行计算
        print("⚡ 第2步：启动并行计算...")
        total_similarities = parallel_blockwise_similarity_map(
            db_path=test_db,
            block_size=10,      # 每块10个RNA（共约10块）
            num_processes=2,    # 2个worker
            output_file=output_file
        )
        
        # 第3步：验证结果
        print("🔍 第3步：验证结果...")
        if os.path.exists(output_file):
            with open(output_file, "r") as f:
                lines = f.readlines()
            print(f"结果文件行数：{len(lines)}")
            
            # 显示前几行作为样本
            print("结果样本：")
            for i, line in enumerate(lines[:5]):
                print(f"  {line.strip()}")
            if len(lines) > 5:
                print(f"  ... (还有 {len(lines)-5} 行)")
        
        # 统计总结
        elapsed = time.time() - start_time
        print(f"\n✅ 快速测试完成！")
        print(f"⏱️  总耗时：{elapsed:.1f}秒")
        print(f"📈 计算相似性对：{total_similarities}")
        print(f"📁 结果文件：{output_file}")
        
        return total_similarities
        
    except Exception as e:
        print(f"❌ 测试失败：{e}")
        import traceback
        traceback.print_exc()
        return None

# -------------------------
# Main Pipeline
# -------------------------
def main():
    print("Starting RNA similarity analysis (serial version)...")
    db_path = "../../../data/ATLAS.db"  # 修复：使用完整的ATLAS数据库
    data = load_data(db_path)
    RNA_motifs = group_motifs_by_RNA(data)
    print(f"Total RNAs: {len(RNA_motifs)}")
    
    print("Calculating full RNA similarity map (serial version)...")
    similarity_map = compute_similarity_map(RNA_motifs)
    save_similarity_map(similarity_map, filename="../../../postprocessing/results2025/similarity_map/similarity_results_ATLAS2025.txt")

if __name__ == "__main__":
    # 运行快速测试
    quick_test_main()
