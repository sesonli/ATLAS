# -*- coding: utf-8 -*-
"""
RNA Motif Database Post-Processing Script
Created on 2025-01-26

This script:
1. Merges motif and pseudoknot databases
2. Adds graph connectivity information
3. Reorders atomic coordinates
4. Removes duplicates
5. Creates final ATLAS database

@author: JingyiLi
"""
import sqlite3
import pickle
import networkx as nx
import json
import os
import re
from datetime import datetime

# Global batch counter for test mode
GLOBAL_BATCH_COUNTER = 0
TEST_MODE = False

def log_progress(message):
    """Log progress with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_message = f"[{timestamp}] {message}"
    print(log_message)

def increment_batch_counter():
    """Increment global batch counter and check if should stop in test mode"""
    global GLOBAL_BATCH_COUNTER
    GLOBAL_BATCH_COUNTER += 1
    if TEST_MODE and GLOBAL_BATCH_COUNTER > 1:
        log_progress(f"TEST MODE: Stopping after processing {GLOBAL_BATCH_COUNTER - 1} batch(es)")
        return True
    return False

def reset_batch_counter():
    """Reset global batch counter"""
    global GLOBAL_BATCH_COUNTER
    GLOBAL_BATCH_COUNTER = 0

def load_pickle_file(pickle_file_path):
    """Load pickle file containing graphs"""
    log_progress(f"Loading pickle file: {pickle_file_path}")
    with open(pickle_file_path, 'rb') as f:
        graph_data = pickle.load(f)
    log_progress(f"Loaded {len(graph_data)} graphs from pickle file")
    return graph_data

def get_graph_from_noncanonical_db(db_path, pdbid):
    """Get a specific graph from noncanonical database by PDB ID"""
    possible_keys = [
        f"{pdbid}_annotation.txt",
        f"pdb{pdbid}.ent_output.txt"
    ]
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    for key in possible_keys:
        cursor.execute("SELECT graph_data FROM graphs WHERE key = ?", (key,))
        row = cursor.fetchone()
        if row:
            try:
                graph = pickle.loads(row[0])
                conn.close()
                if isinstance(graph, nx.Graph):
                    nodes = list(graph.nodes)
                    edges = list(graph.edges(data=True))
                    return {'nodes': nodes, 'edges': edges}
                else:
                    log_progress(f"Warning: Unexpected graph type for {pdbid}: {type(graph)}")
                    return None
            except Exception as e:
                log_progress(f"Warning: Failed to load graph for {pdbid}: {e}")
                conn.close()
                return None
    
    conn.close()
    return None

def get_batch_graphs_from_noncanonical_db(db_path, pdbids):
    """Get multiple graphs from noncanonical database by batch of PDB IDs"""
    graphs = {}
    
    # Create all possible keys for the batch
    key_to_pdbid = {}
    for pdbid in pdbids:
        possible_keys = [
            f"{pdbid}_annotation.txt",
            f"pdb{pdbid}.ent_output.txt"
        ]
        for key in possible_keys:
            key_to_pdbid[key] = pdbid
    
    if not key_to_pdbid:
        return graphs
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Build IN clause for batch query
    placeholders = ','.join(['?'] * len(key_to_pdbid))
    cursor.execute(f"SELECT key, graph_data FROM graphs WHERE key IN ({placeholders})", 
                   list(key_to_pdbid.keys()))
    
    for key, blob_data in cursor.fetchall():
        try:
            graph = pickle.loads(blob_data)
            pdbid = key_to_pdbid[key]
            if isinstance(graph, nx.Graph):
                nodes = list(graph.nodes)
                edges = list(graph.edges(data=True))
                graphs[pdbid] = {'nodes': nodes, 'edges': edges}
            else:
                log_progress(f"Warning: Unexpected graph type for {pdbid}: {type(graph)}")
        except Exception as e:
            pdbid = key_to_pdbid[key]
            log_progress(f"Warning: Failed to load graph for {pdbid}: {e}")
    
    conn.close()
    return graphs

def sort_atomic_coordinates(filecontent):
    """Sort atomic coordinates by chain ID and residue number"""
    if not filecontent or filecontent.strip() == "":
        return filecontent
    
    lines = filecontent.strip().split("\n")
    blocks = {}

    # Group lines by (chain_id, residue_number)
    for line in lines:
        if line.strip() == "":
            continue
        try:
            # Extract key (column 21 and trimmed column 22-26 as integer)
            if len(line) >= 26:
                chain_id = line[21]
                residue_num = int(line[22:26].strip())
                key = (chain_id, residue_num)
            else:
                # If line is too short, use the line itself as key
                key = (line, 0)
        except (ValueError, IndexError):
            # If parsing fails, use the line itself as key
            key = (line, 0)

        if key not in blocks:
            blocks[key] = []
        blocks[key].append(line)

    # Sort the blocks based on the extracted key
    sorted_blocks = sorted(blocks.items(), key=lambda x: (x[0][0], x[0][1]))

    # Reassemble the sorted blocks into a single string
    result = []
    for key, block in sorted_blocks:
        result.extend(block)

    return "\n".join(result)

def extract_graph_data_for_pdbid(graph_data, pdbid):
    """Extract graph information for a specific PDB ID"""
    # Try different key formats
    possible_keys = [
        f"{pdbid}_annotation.txt",
        f"pdb{pdbid}.ent_output.txt"
    ]
    
    for key in possible_keys:
        if key in graph_data:
            graph = graph_data[key]
            if isinstance(graph, nx.Graph):
                nodes = list(graph.nodes)
                edges = list(graph.edges(data=True))
                return {'nodes': nodes, 'edges': edges}
            else:
                log_progress(f"Warning: Unexpected graph type for {pdbid}: {type(graph)}")
                return None
    
    log_progress(f"Warning: No graph found for {pdbid}")
    return None

def filter_edges(motif_nodes, all_edges):
    """Filter edges to include only those connecting nodes in the motif"""
    motif_node_set = set(motif_nodes)
    filtered_edges = []
    
    for edge in all_edges:
        u, v, attr = edge
        if u in motif_node_set and v in motif_node_set:
            # Convert edge to serializable format
            filtered_edges.append([u, v, attr])
    
    return filtered_edges

def add_graph_info_to_table_batch(db_path, table_name, canonical_graphs, noncanonical_db_path, batch_size=1000):
    """Add graph connectivity information to database table using batch processing"""
    log_progress(f"Adding graph information to table: {table_name} (batch size: {batch_size})")
    if TEST_MODE:
        log_progress("TEST MODE: Will stop after 1 batch")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Add columns if they don't exist
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN non_graph_edges TEXT")
    except sqlite3.OperationalError:
        pass  # Column already exists
    
    try:
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN cano_graph_edges TEXT")
    except sqlite3.OperationalError:
        pass  # Column already exists
    
    # Get total count
    cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
    total_rows = cursor.fetchone()[0]
    log_progress(f"Processing {total_rows} rows in batches of {batch_size}")
    
    processed = 0
    for offset in range(0, total_rows, batch_size):
        # Check batch counter for test mode
        if increment_batch_counter():
            break
        # Get batch of rows
        cursor.execute(f"SELECT id, pdbid, nt_number FROM {table_name} LIMIT ? OFFSET ?", 
                      (batch_size, offset))
        batch_rows = cursor.fetchall()
        
        # Extract unique pdbids from this batch
        batch_pdbids = list(set([row[1] for row in batch_rows]))
        
        # Batch load noncanonical graphs for this batch only
        noncanonical_graphs_batch = get_batch_graphs_from_noncanonical_db(noncanonical_db_path, batch_pdbids)
        
        # Process each row in the batch
        for row_id, pdbid, nt_number in batch_rows:
            processed += 1
            if processed % 100 == 0:
                log_progress(f"Processed {processed}/{total_rows} rows for graph info")
            
            motif_nodes = nt_number.split(',')
            
            # Add noncanonical graph edges (from batch-loaded data)
            if pdbid in noncanonical_graphs_batch:
                noncanonical_graph = noncanonical_graphs_batch[pdbid]
                noncanonical_edges = filter_edges(motif_nodes, noncanonical_graph['edges'])
                cursor.execute(f"UPDATE {table_name} SET non_graph_edges = ? WHERE id = ?", 
                              (json.dumps(noncanonical_edges), row_id))
            
            # Add canonical graph edges (from pre-loaded data)
            canonical_graph = extract_graph_data_for_pdbid(canonical_graphs, pdbid)
            if canonical_graph:
                canonical_edges = filter_edges(motif_nodes, canonical_graph['edges'])
                cursor.execute(f"UPDATE {table_name} SET cano_graph_edges = ? WHERE id = ?", 
                              (json.dumps(canonical_edges), row_id))
        
        # Commit after each batch
        conn.commit()
        log_progress(f"Completed batch {offset//batch_size + 1}/{(total_rows + batch_size - 1)//batch_size}")
    
    conn.close()
    log_progress(f"Graph information added to {processed} rows in {table_name}")

def reorder_nucleotides_in_table_batch(db_path, table_name, batch_size=1000):
    """Reorder atomic coordinates in filecontent column using batch processing"""
    log_progress(f"Reordering nucleotides in table: {table_name} (batch size: {batch_size})")
    if TEST_MODE:
        log_progress("TEST MODE: Will stop after 1 batch")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Get total count
    cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
    total_rows = cursor.fetchone()[0]
    log_progress(f"Processing {total_rows} rows in batches of {batch_size}")
    
    processed = 0
    for offset in range(0, total_rows, batch_size):
        # Check batch counter for test mode
        if increment_batch_counter():
            break
        # Get batch of rows
        cursor.execute(f"SELECT id, filecontent FROM {table_name} LIMIT ? OFFSET ?", 
                      (batch_size, offset))
        batch_rows = cursor.fetchall()
        
        # Process each row in the batch
        for row_id, filecontent in batch_rows:
            processed += 1
            if processed % 100 == 0:
                log_progress(f"Reordered {processed}/{total_rows} entries")
            
            if filecontent:
                sorted_filecontent = sort_atomic_coordinates(filecontent)
                cursor.execute(f"UPDATE {table_name} SET filecontent = ? WHERE id = ?", 
                              (sorted_filecontent, row_id))
        
        # Commit after each batch
        conn.commit()
        log_progress(f"Completed batch {offset//batch_size + 1}/{(total_rows + batch_size - 1)//batch_size}")
    
    conn.close()
    log_progress(f"Nucleotides reordered in {processed} rows")

def remove_duplicates_from_table(db_path, table_name):
    """Remove duplicate entries based on motif_type, pdbid, nt_number"""
    log_progress(f"Removing duplicates from table: {table_name}")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Count duplicates before removal
    cursor.execute(f"""
        SELECT COUNT(*) FROM {table_name}
        WHERE rowid NOT IN (
            SELECT MIN(rowid) 
            FROM {table_name} 
            GROUP BY motif_type, pdbid, nt_number
        )
    """)
    duplicate_count = cursor.fetchone()[0]
    
    if duplicate_count > 0:
        log_progress(f"Found {duplicate_count} duplicate entries to remove")
        
        # Remove duplicates, keeping the first occurrence
        cursor.execute(f"""
            DELETE FROM {table_name} 
            WHERE rowid NOT IN (
                SELECT MIN(rowid) 
                FROM {table_name} 
                GROUP BY motif_type, pdbid, nt_number
            )
        """)
        
        conn.commit()
        log_progress(f"Removed {duplicate_count} duplicate entries")
    else:
        log_progress("No duplicates found")
    
    conn.close()

def create_final_database(motif_db_path, pk_db_path, output_db_path):
    """Create final database by merging motif and PK databases"""
    log_progress("Creating final merged database")
    if TEST_MODE:
        log_progress("TEST MODE: Creating database with limited data")
    
    # Connect to all databases
    motif_conn = sqlite3.connect(motif_db_path)
    pk_conn = sqlite3.connect(pk_db_path)
    output_conn = sqlite3.connect(output_db_path)
    
    motif_cursor = motif_conn.cursor()
    pk_cursor = pk_conn.cursor()
    output_cursor = output_conn.cursor()
    
    # Drop existing tables
    output_cursor.execute('DROP TABLE IF EXISTS data')
    output_cursor.execute('DROP TABLE IF EXISTS PK')
    
    # Create data table (exactly matching ATLAS.db structure)
    output_cursor.execute('''
        CREATE TABLE data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            motif_type TEXT,
            pdbid TEXT,
            paired_nt_number TEXT,
            nt_number TEXT,
            filecontent TEXT,
            non_graph_edges TEXT,
            cano_graph_edges TEXT
        )
    ''')
    
    # Create PK table (exactly matching ATLAS.db structure - note the special format)
    output_cursor.execute('''
        CREATE TABLE PK (
            id,
            motif_type,
            pdbid,
            paired_nt_number,
            nt_number,
            file_content,
            non_graph_edges,
            cano_graph_edges,
            dot_bracket
        )
    ''')
    
    # Copy motif data (exclude created_at column)
    try:
        motif_cursor.execute('''
            SELECT motif_type, pdbid, paired_nt_number, nt_number, filecontent, 
                   non_graph_edges, cano_graph_edges 
            FROM files
        ''')
        motif_rows = motif_cursor.fetchall()
        
        output_cursor.executemany('''
            INSERT INTO data (motif_type, pdbid, paired_nt_number, nt_number, filecontent, 
                            non_graph_edges, cano_graph_edges)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', motif_rows)
        
        log_progress(f"Copied {len(motif_rows)} motif entries to data table")
    except sqlite3.OperationalError as e:
        log_progress(f"Warning: Could not copy all motif columns: {e}")
        # Try with basic columns only
        motif_cursor.execute('SELECT motif_type, pdbid, paired_nt_number, nt_number, filecontent FROM files')
        motif_rows = motif_cursor.fetchall()
        output_cursor.executemany('''
            INSERT INTO data (motif_type, pdbid, paired_nt_number, nt_number, filecontent)
            VALUES (?, ?, ?, ?, ?)
        ''', motif_rows)
        log_progress(f"Copied {len(motif_rows)} motif entries with basic columns")
    
    # Copy PK data (map filecontent to file_content)
    try:
        pk_cursor.execute('''
            SELECT motif_type, pdbid, paired_nt_number, nt_number, filecontent, 
                   non_graph_edges, cano_graph_edges, dot_bracket 
            FROM PK
        ''')
        pk_rows = pk_cursor.fetchall()
        
        output_cursor.executemany('''
            INSERT INTO PK (motif_type, pdbid, paired_nt_number, nt_number, file_content, 
                          non_graph_edges, cano_graph_edges, dot_bracket)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ''', pk_rows)
        
        log_progress(f"Copied {len(pk_rows)} pseudoknot entries to PK table")
    except sqlite3.OperationalError as e:
        log_progress(f"Warning: Could not copy all PK columns: {e}")
        # Try with basic columns only
        pk_cursor.execute('SELECT motif_type, pdbid, paired_nt_number, nt_number, filecontent, dot_bracket FROM PK')
        pk_rows = pk_cursor.fetchall()
        output_cursor.executemany('''
            INSERT INTO PK (motif_type, pdbid, paired_nt_number, nt_number, file_content, dot_bracket)
            VALUES (?, ?, ?, ?, ?, ?)
        ''', pk_rows)
        log_progress(f"Copied {len(pk_rows)} pseudoknot entries with basic columns")
    
    # Commit and close connections
    output_conn.commit()
    motif_conn.close()
    pk_conn.close()
    output_conn.close()
    
    log_progress("Final database created successfully")

def main(test_mode=False):
    """Main post-processing workflow"""
    global TEST_MODE
    TEST_MODE = test_mode
    
    log_progress("Starting RNA motif database post-processing")
    if TEST_MODE:
        log_progress("RUNNING IN TEST MODE - Will process only 1 batch per operation")
        log_progress("TEST MODE: Original databases will NOT be modified")
    
    # Reset batch counter
    reset_batch_counter()
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, 'data')
    
    # Fixed input and output paths
    motif_db_path = os.path.join(data_dir, 'all_motifs_merged.db')
    pk_db_path = os.path.join(data_dir, 'ATLAS20250703_pk.db')
    
    # Use different output path for test mode to avoid overwriting
    if TEST_MODE:
        output_db_path = os.path.join(data_dir, 'ATLAS2025_TEST.db')
    else:
        output_db_path = os.path.join(data_dir, 'ATLAS2025.db')
    
    log_progress(f"Input motif database: {motif_db_path}")
    log_progress(f"Input PK database: {pk_db_path}")
    log_progress(f"Output database: {output_db_path}")
    
    # Check if input databases exist
    if not os.path.exists(motif_db_path):
        log_progress(f"Error: Motif database not found: {motif_db_path}")
        return
    
    if not os.path.exists(pk_db_path):
        log_progress(f"Error: PK database not found: {pk_db_path}")
        return
    
    # Setup graph data sources
    canonical_pickle = os.path.join(data_dir, "saved_graphs_canonical_new.pickle")
    noncanonical_db = os.path.join(data_dir, "noncanonical_graphs.db")
    
    if not os.path.exists(canonical_pickle):
        log_progress(f"Error: Canonical graphs not found: {canonical_pickle}")
        return
    
    if not os.path.exists(noncanonical_db):
        log_progress(f"Error: Noncanonical graphs database not found: {noncanonical_db}")
        return
    
    # Load canonical graphs (small enough to keep in memory)
    canonical_graphs = load_pickle_file(canonical_pickle)
    
    # In test mode, work with copies to avoid modifying original databases
    if TEST_MODE:
        import shutil
        test_motif_db = os.path.join(data_dir, 'test_motifs.db')
        test_pk_db = os.path.join(data_dir, 'test_pk.db')
        
        log_progress("TEST MODE: Creating temporary copies of databases...")
        shutil.copy2(motif_db_path, test_motif_db)
        shutil.copy2(pk_db_path, test_pk_db)
        
        motif_db_work = test_motif_db
        pk_db_work = test_pk_db
    else:
        motif_db_work = motif_db_path
        pk_db_work = pk_db_path
    
    # Process motif database with batch processing (noncanonical graphs loaded on-demand)
    log_progress("Processing motif database...")
    reset_batch_counter()
    add_graph_info_to_table_batch(motif_db_work, "files", canonical_graphs, noncanonical_db, batch_size=1000)
    reset_batch_counter()
    reorder_nucleotides_in_table_batch(motif_db_work, "files", batch_size=1000)
    if not TEST_MODE:
        remove_duplicates_from_table(motif_db_work, "files")
    
    # Process PK database with batch processing (noncanonical graphs loaded on-demand)
    log_progress("Processing PK database...")
    reset_batch_counter()
    add_graph_info_to_table_batch(pk_db_work, "PK", canonical_graphs, noncanonical_db, batch_size=1000)
    reset_batch_counter()
    reorder_nucleotides_in_table_batch(pk_db_work, "PK", batch_size=1000)
    if not TEST_MODE:
        remove_duplicates_from_table(pk_db_work, "PK")
    
    # Create final merged database
    create_final_database(motif_db_work, pk_db_work, output_db_path)
    
    # Final deduplication across merged data
    if not TEST_MODE:
        log_progress("Performing final deduplication...")
        remove_duplicates_from_table(output_db_path, "data")
        remove_duplicates_from_table(output_db_path, "PK")
    
    # Final statistics
    conn = sqlite3.connect(output_db_path)
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM data")
    motif_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM PK")
    pk_count = cursor.fetchone()[0]
    
    conn.close()
    
    log_progress("=" * 50)
    if TEST_MODE:
        log_progress("TEST MODE COMPLETE!")
        log_progress(f"Test database created: {output_db_path}")
        log_progress(f"Processed 1 batch from each operation")
    else:
        log_progress("POST-PROCESSING COMPLETE!")
        log_progress(f"Final database: {output_db_path}")
    log_progress(f"Regular motifs: {motif_count}")
    log_progress(f"Pseudoknots: {pk_count}")
    log_progress(f"Total entries: {motif_count + pk_count}")
    log_progress("=" * 50)
    
    # Clean up test files
    if TEST_MODE:
        try:
            os.remove(test_motif_db)
            os.remove(test_pk_db)
            log_progress("Cleaned up temporary test files")
        except:
            pass

if __name__ == "__main__":
    import sys
    test_mode = len(sys.argv) > 1 and sys.argv[1] == "--test"
    main(test_mode=test_mode) 