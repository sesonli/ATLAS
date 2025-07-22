#!/usr/bin/env python3
"""
Parallel version of find_all_motifs.py for test.
Uses 8 processes to split motif search across chunks of graphs.
Generates test DBs under 'parallel_test_db' directory.
"""

import os
import sqlite3
import pickle
import time
from multiprocessing import Pool

import find_all_motifs


def init_db(db_path):
    conn = sqlite3.connect(db_path)
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
    conn.commit()
    return conn


def process_chunk(args):
    idx, graph_items = args
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parallel_dir = os.path.join(current_dir, 'data', 'parallel_test_db')
    os.makedirs(parallel_dir, exist_ok=True)
    db_path = os.path.join(parallel_dir, f'all_motifs_test_part_{idx}.db')
    conn = init_db(db_path)

    graphs = dict(graph_items)

    data_dir = os.path.join(current_dir, 'data')
    motif_files = {
        'hairpin': 'hairpin.pickle',
        'bulge_internal_junction': 'bulge_internal_junction.pickle'
    }

    for motif_name, pickle_filename in motif_files.items():
        target_path = os.path.join(data_dir, pickle_filename)
        if not os.path.exists(target_path):
            print(f"[Process {idx}] Warning: {target_path} not found, skipping {motif_name}")
            continue
        with open(target_path, 'rb') as f:
            target_graphs = pickle.load(f)
        print(f"[Process {idx}] Loaded {len(target_graphs)} target graphs for {motif_name}")
        start_time = time.time()

        if motif_name == 'hairpin':
            find_all_motifs.find_hairpin_motifs(graphs, target_graphs, conn)
        else:
            find_all_motifs.find_other_motifs(graphs, target_graphs, conn, motif_name)

        elapsed = time.time() - start_time
        print(f"[Process {idx}] Completed {motif_name} in {elapsed:.2f}s")

    conn.close()
    print(f"[Process {idx}] Done")


def chunkify(seq, n):
    """Split sequence into n chunks (as equal as possible)."""
    k, m = divmod(len(seq), n)
    return [seq[i*k + min(i, m):(i+1)*k + min(i+1, m)] for i in range(n)]


if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, 'data')
    graph_pickle = os.path.join(data_dir, 'saved_graphs_canonical_new.pickle')
    if not os.path.exists(graph_pickle):
        print(f"Error: Graph pickle {graph_pickle} not found!")
        exit(1)
    with open(graph_pickle, 'rb') as f:
        graphs = pickle.load(f)
    graph_items = list(graphs.items())
    n_procs = 8
    chunks = chunkify(graph_items, n_procs)
    tasks = [(i, chunks[i]) for i in range(len(chunks))]

    print(f"Starting parallel motif search with {n_procs} processes")
    start_total = time.time()
    with Pool(processes=n_procs) as pool:
        pool.map(process_chunk, tasks)
    end_total = time.time()
    print(f"Parallel test completed in {end_total - start_total:.2f}s")
    print("Test DB files are under 'data/parallel_test_db' directory") 