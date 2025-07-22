#!/usr/bin/env python3
"""
Convert pickle batch files to SQLite database
This is used for noncanonical graphs
Uses streaming approach to avoid memory issues
"""

import sqlite3
import pickle
import glob
import os
import time
from pathlib import Path

def get_memory_usage():
    """Get current memory usage in MB (simplified version)"""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        return 0.0  # Return 0 if psutil not available

def convert_pickles_to_db(pickle_pattern, output_db):
    """
    Convert pickle batch files to SQLite database
    
    Args:
        pickle_pattern: Pattern to match pickle files (e.g., 'data/noncanonical_graphs/batch_*.pickle')
        output_db: Output SQLite database path
    """
    
    # Get all pickle files
    pickle_files = glob.glob(pickle_pattern)
    pickle_files.sort()  # Process in order
    
    if not pickle_files:
        print(f"No pickle files found matching pattern: {pickle_pattern}")
        return
    
    print(f"Found {len(pickle_files)} pickle files to process")
    print(f"Output database: {output_db}")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_db), exist_ok=True)
    
    # Remove existing database if it exists
    if os.path.exists(output_db):
        print(f"Removing existing database: {output_db}")
        os.remove(output_db)
    
    # Connect to database
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()
    
    # Create table
    cursor.execute('''
        CREATE TABLE graphs (
            key TEXT PRIMARY KEY,
            graph_data BLOB,
            batch_file TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Process each pickle file
    total_keys = 0
    start_time = time.time()
    
    for i, pickle_file in enumerate(pickle_files):
        file_start_time = time.time()
        batch_name = os.path.basename(pickle_file)
        
        print(f"\nProcessing [{i+1}/{len(pickle_files)}]: {batch_name}")
        print(f"Memory usage: {get_memory_usage():.1f} MB")
        
        try:
            # Load pickle file
            with open(pickle_file, 'rb') as f:
                data = pickle.load(f)
            
            if not isinstance(data, dict):
                print(f"Warning: {pickle_file} does not contain a dictionary, skipping")
                continue
            
            # Insert data into database
            batch_keys = 0
            for key, graph in data.items():
                # Serialize the graph object
                graph_blob = pickle.dumps(graph, protocol=pickle.HIGHEST_PROTOCOL)
                
                # Insert into database
                cursor.execute(
                    'INSERT OR REPLACE INTO graphs (key, graph_data, batch_file) VALUES (?, ?, ?)',
                    (key, graph_blob, batch_name)
                )
                batch_keys += 1
                total_keys += 1
                
                # Commit every 1000 records to free memory
                if batch_keys % 1000 == 0:
                    conn.commit()
            
            # Final commit for this file
            conn.commit()
            
            # Clear the data from memory
            del data
            
            file_time = time.time() - file_start_time
            print(f"  Processed {batch_keys} keys in {file_time:.1f}s")
            
        except Exception as e:
            print(f"Error processing {pickle_file}: {e}")
            continue
    
    # Create index for faster lookups
    print("\nCreating index...")
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_key ON graphs(key)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_batch_file ON graphs(batch_file)')
    conn.commit()
    
    # Get final statistics
    cursor.execute('SELECT COUNT(*) FROM graphs')
    final_count = cursor.fetchone()[0]
    
    cursor.execute('SELECT COUNT(DISTINCT batch_file) FROM graphs')
    batch_count = cursor.fetchone()[0]
    
    conn.close()
    
    total_time = time.time() - start_time
    db_size = os.path.getsize(output_db) / 1024 / 1024  # MB
    
    print(f"\n=== Conversion Complete ===")
    print(f"Total keys processed: {total_keys}")
    print(f"Keys in database: {final_count}")
    print(f"Batch files processed: {batch_count}")
    print(f"Total time: {total_time:.1f}s")
    print(f"Database size: {db_size:.1f} MB")
    print(f"Final memory usage: {get_memory_usage():.1f} MB")

def test_database(db_path):
    """Test the created database"""
    print(f"\n=== Testing Database: {db_path} ===")
    
    if not os.path.exists(db_path):
        print("Database file does not exist!")
        return False
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check table structure
        cursor.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='graphs'")
        table_info = cursor.fetchone()
        print(f"Table structure: {table_info[0] if table_info else 'Table not found'}")
        
        # Count total records
        cursor.execute('SELECT COUNT(*) FROM graphs')
        total_count = cursor.fetchone()[0]
        print(f"Total records: {total_count}")
        
        # Count distinct batch files
        cursor.execute('SELECT COUNT(DISTINCT batch_file) FROM graphs')
        batch_count = cursor.fetchone()[0]
        print(f"Distinct batch files: {batch_count}")
        
        # Sample a few records
        cursor.execute('SELECT key, batch_file FROM graphs LIMIT 5')
        samples = cursor.fetchall()
        print(f"Sample records:")
        for key, batch_file in samples:
            print(f"  {key} (from {batch_file})")
        
        # Test loading a graph object
        cursor.execute('SELECT key, graph_data FROM graphs LIMIT 1')
        result = cursor.fetchone()
        if result:
            key, graph_blob = result
            try:
                graph = pickle.loads(graph_blob)
                print(f"Successfully loaded graph for key '{key}': {type(graph)}")
                # If it's a networkx graph, show some info
                if hasattr(graph, 'number_of_nodes'):
                    print(f"  Nodes: {graph.number_of_nodes()}, Edges: {graph.number_of_edges()}")
            except Exception as e:
                print(f"Error loading graph data: {e}")
        
        conn.close()
        print("Database test completed successfully!")
        return True
        
    except Exception as e:
        print(f"Database test failed: {e}")
        return False

if __name__ == "__main__":
    # Configuration
    pickle_pattern = "data/noncanonical_graphs/batch_*.pickle"
    output_db = "data/noncanonical_graphs.db"
    
    print("=== Pickle to Database Converter ===")
    print(f"Input pattern: {pickle_pattern}")
    print(f"Output database: {output_db}")
    print(f"Initial memory usage: {get_memory_usage():.1f} MB")
    
    # Convert pickles to database
    convert_pickles_to_db(pickle_pattern, output_db)
    
    # Test the database
    test_database(output_db)
    
    print("\nConversion and testing complete!") 