#!/usr/bin/env python3
"""
SQLite Database Merger
Merge SQLite database shards in motif_database directory
"""

import sqlite3
import os
import time
from pathlib import Path
from datetime import datetime

class MergeManager:
    def __init__(self):
        self.start_time = time.time()
        
    def log(self, message):
        """Print timestamped log message"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")
    
    def get_memory_usage(self):
        """Get current memory usage"""
        import tracemalloc
        try:
            tracemalloc.start()
            current, peak = tracemalloc.get_traced_memory()
            return current / 1024 / 1024
        except:
            return 0.0
    
    def merge_sqlite_databases(self):
        """Merge SQLite database files"""
        self.log("Starting SQLite database merge...")
        
        # Locate source files
        db_dir = Path("data/motif_database")
        source_files = sorted(db_dir.glob("all_motifs_test_part_*.db"))
        target_file = "all_motifs_merged.db"
        
        if not source_files:
            self.log("Error: No source database files found")
            return False
        
        self.log(f"Found {len(source_files)} source database files")
        
        # Remove existing target file
        if os.path.exists(target_file):
            os.remove(target_file)
            self.log(f"Removed existing target file: {target_file}")
        
        try:
            # 创建目标数据库
            target_conn = sqlite3.connect(target_file)
            target_cursor = target_conn.cursor()
            
            # 创建表结构 (基于检查结果)
            target_cursor.execute('''
                CREATE TABLE files (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    motif_type TEXT,
                    pdbid TEXT,
                    paired_nt_number TEXT,
                    nt_number TEXT,
                    filecontent TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            self.log("Target database table structure created")
            
            # Count total records for progress display
            total_records = 0
            for source_file in source_files:
                conn = sqlite3.connect(source_file)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM files")
                count = cursor.fetchone()[0]
                total_records += count
                conn.close()
            
            self.log(f"Total records: {total_records:,}")
            
            # Start merging
            processed_records = 0
            current_id = 1
            
            for i, source_file in enumerate(source_files):
                self.log(f"Processing file {i+1}/{len(source_files)}: {source_file.name}")
                
                source_conn = sqlite3.connect(source_file)
                source_cursor = source_conn.cursor()
                
                # Get current file record count
                source_cursor.execute("SELECT COUNT(*) FROM files")
                file_record_count = source_cursor.fetchone()[0]
                
                # Read and insert data in batches
                batch_size = 10000
                offset = 0
                
                while offset < file_record_count:
                    # Read a batch of data (excluding id field)
                    source_cursor.execute('''
                        SELECT motif_type, pdbid, paired_nt_number, nt_number, filecontent, created_at
                        FROM files 
                        LIMIT ? OFFSET ?
                    ''', (batch_size, offset))
                    
                    rows = source_cursor.fetchall()
                    if not rows:
                        break
                    
                    # Prepare insert data (reassign IDs)
                    insert_data = []
                    for row in rows:
                        insert_data.append((current_id,) + row)
                        current_id += 1
                    
                    # Batch insert
                    target_cursor.executemany('''
                        INSERT INTO files (id, motif_type, pdbid, paired_nt_number, nt_number, filecontent, created_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', insert_data)
                    
                    target_conn.commit()
                    
                    processed_records += len(rows)
                    offset += batch_size
                    
                    # Show progress
                    progress = (processed_records / total_records) * 100
                    self.log(f"  Progress: {processed_records:,}/{total_records:,} ({progress:.1f}%)")
                
                source_conn.close()
                self.log(f"  Completed file: {source_file.name} ({file_record_count:,} records)")
            
            # Create indexes to optimize query performance
            self.log("Creating indexes...")
            target_cursor.execute("CREATE INDEX idx_motif_type ON files(motif_type)")
            target_cursor.execute("CREATE INDEX idx_pdbid ON files(pdbid)")
            target_conn.commit()
            
            target_conn.close()
            
            # Verify results
            final_conn = sqlite3.connect(target_file)
            final_cursor = final_conn.cursor()
            final_cursor.execute("SELECT COUNT(*) FROM files")
            final_count = final_cursor.fetchone()[0]
            final_conn.close()
            
            self.log(f"Merge completed! Final record count: {final_count:,}")
            
            if final_count == total_records:
                self.log("✓ Data validation passed")
                return True
            else:
                self.log("✗ Data validation failed")
                return False
                
        except Exception as e:
            self.log(f"Error merging SQLite databases: {e}")
            return False
    

    
    def validate_results(self):
        """Validate merge results"""
        self.log("Validating merge results...")
        
        # Validate SQLite file
        sqlite_file = "all_motifs_merged.db"
        if os.path.exists(sqlite_file):
            size_mb = os.path.getsize(sqlite_file) / (1024**2)
            self.log(f"✓ SQLite merge file: {sqlite_file} ({size_mb:.1f} MB)")
        else:
            self.log(f"✗ SQLite merge file not found: {sqlite_file}")
        

    
    def run(self):
        """Run complete merge workflow"""
        self.log("="*60)
        self.log("Starting data merge program")
        self.log("="*60)
        
        # Show system info
        self.log("Starting merge operations (System RAM: 64GB)")
        
        success_count = 0
        
        # Merge SQLite databases
        if self.merge_sqlite_databases():
            success_count += 1
        
        # Validate results
        self.validate_results()
        
        # Summary
        elapsed_time = time.time() - self.start_time
        self.log("="*60)
        self.log(f"Merge program completed! Success: {success_count}/1")
        self.log(f"Total time: {elapsed_time/60:.1f} minutes")
        self.log("="*60)

def main():
    """Main program entry point"""
    merger = MergeManager()
    merger.run()

if __name__ == "__main__":
    main() 