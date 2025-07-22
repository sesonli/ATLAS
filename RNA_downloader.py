#!/usr/bin/env python3
"""
RNA Processing Pipeline - Space-Efficient Combined Script
=========================================================
Combines 4 steps in one script to save disk space:
1. Download PDB files
2. Extract RNA chains
3. Apply corrections (atom number overflow)
4. Fix chain IDs (two-letter to single-letter)

Processes one file at a time to minimize disk usage.
Only keeps final corrected files.
"""

import requests
import os
import time
from pathlib import Path
import string
from datetime import datetime
import tempfile
import argparse
import json

class RNAPipeline:
    def __init__(self, output_dir="rna_final", log_dir="logs", delay=0.1, timeout=30):
        self.output_dir = Path(output_dir)
        self.log_dir = Path(log_dir)
        self.delay = delay
        self.timeout = timeout
        
        # Create directories
        self.output_dir.mkdir(exist_ok=True)
        self.log_dir.mkdir(exist_ok=True)
        
        # Statistics
        self.stats = {
            'total': 0,
            'downloaded': 0,
            'rna_found': 0,
            'corrected': 0,
            'failed': 0,
            'already_exists': 0
        }
        
        # Setup logging
        self.log_file = self.log_dir / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    def log(self, message, level="INFO"):
        """Log message to both console and file"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] [{level}] {message}"
        print(log_entry)
        
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(log_entry + '\n')
    
    def get_rna_structure_list(self):
        """Get list of all RNA-containing structures from PDB API"""
        self.log("Fetching RNA structure list from PDB API...")
        
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": "RNA"
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": 10000
                }
            }
        }
        
        try:
            response = requests.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=query,
                timeout=self.timeout
            )
            response.raise_for_status()
            
            data = response.json()
            structure_ids = [item['identifier'] for item in data['result_set']]
            
            self.log(f"Found {len(structure_ids)} RNA structures from PDB API")
            return structure_ids
            
        except Exception as e:
            self.log(f"Error fetching RNA structure list: {e}", "ERROR")
            return []
    
    def download_pdb(self, pdb_id):
        """Download PDB content to memory"""
        try:
            # Try PDB format first
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(pdb_url, timeout=self.timeout)
            
            if response.status_code == 200:
                return response.text, 'pdb'
            
            # Try CIF format if PDB failed
            cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            response = requests.get(cif_url, timeout=self.timeout)
            
            if response.status_code == 200:
                return response.text, 'cif'
            
            self.log(f"Failed to download {pdb_id}: HTTP {response.status_code}", "ERROR")
            return None, None
            
        except Exception as e:
            self.log(f"Error downloading {pdb_id}: {e}", "ERROR")
            return None, None
    
    def find_rna_chains_pdb(self, pdb_content):
        """Find chains that contain A, U, G, or C residues in PDB format"""
        rna_chains = set()
        rna_residues = {'A', 'U', 'G', 'C'}
        
        lines = pdb_content.strip().split('\n')
        in_model_1 = True
        model_found = False
        
        for line in lines:
            if line.startswith('MODEL'):
                model_found = True
                model_num = line.split()[1] if len(line.split()) > 1 else '1'
                in_model_1 = (model_num == '1')
                continue
            
            if line.startswith('ENDMDL'):
                if model_found and in_model_1:
                    break
                in_model_1 = False
                continue
            
            if line.startswith(('ATOM', 'HETATM')) and in_model_1:
                parts = line.split()
                if len(parts) >= 5:
                    residue_name = parts[3]
                    chain_and_residue = parts[4]
                    chain_id = chain_and_residue[0]
                    
                    if residue_name in rna_residues:
                        rna_chains.add(chain_id)
        
        return rna_chains
    
    def find_rna_chains_cif(self, cif_content):
        """Find chains that contain A, U, G, or C residues in CIF format"""
        rna_chains = set()
        rna_residues = {'A', 'U', 'G', 'C'}
        
        lines = cif_content.strip().split('\n')
        in_atom_section = False
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('_atom_site.pdbx_PDB_model_num'):
                in_atom_section = True
                continue
            
            if in_atom_section and (not line or line.startswith('_') or line.startswith('#')):
                break
            
            if in_atom_section and line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 7:
                    residue_name = parts[5]
                    chain_id = parts[6]
                    
                    if residue_name in rna_residues:
                        rna_chains.add(chain_id)
        
        return rna_chains
    
    def extract_rna_chains_pdb(self, pdb_content, rna_chains):
        """Extract RNA chains from PDB content"""
        lines = pdb_content.strip().split('\n')
        output_lines = []
        
        in_model_1 = True
        model_found = False
        
        for line in lines:
            if line.startswith('MODEL'):
                model_found = True
                model_num = line.split()[1] if len(line.split()) > 1 else '1'
                in_model_1 = (model_num == '1')
                continue
            
            if line.startswith('ENDMDL'):
                if model_found and in_model_1:
                    break
                in_model_1 = False
                continue
            
            if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK')):
                output_lines.append(line)
            
            elif line.startswith(('ATOM', 'HETATM')) and in_model_1:
                parts = line.split()
                if len(parts) >= 5:
                    chain_and_residue = parts[4]
                    chain_id = chain_and_residue[0]
                    if chain_id in rna_chains:
                        output_lines.append(line)
            
            elif line.startswith(('TER', 'CONECT')) and in_model_1:
                output_lines.append(line)
            
            elif line.startswith('END'):
                output_lines.append(line)
        
        return '\n'.join(output_lines)
    
    def extract_rna_chains_cif(self, cif_content, rna_chains):
        """Extract RNA chains from CIF content and convert to PDB format"""
        lines = cif_content.strip().split('\n')
        output_lines = ["HEADER    RNA STRUCTURE"]
        
        in_atom_section = False
        atom_count = 1
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('_atom_site.pdbx_PDB_model_num'):
                in_atom_section = True
                continue
            
            if in_atom_section and (not line or line.startswith('_') or line.startswith('#')):
                break
            
            if in_atom_section and line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 15:
                    chain_id = parts[6]
                    
                    if chain_id in rna_chains:
                        atom_name = parts[3]
                        residue_name = parts[5]
                        residue_num = parts[8]
                        x, y, z = parts[10], parts[11], parts[12]
                        occupancy, b_factor = parts[13], parts[14]
                        element = parts[2]
                        
                        pdb_line = f"ATOM  {atom_count:5d}  {atom_name:<4s} {residue_name:>3s} {chain_id}{residue_num:>4s}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}{float(occupancy):6.2f}{float(b_factor):6.2f}           {element:>2s}"
                        output_lines.append(pdb_line)
                        atom_count += 1
        
        output_lines.append("END")
        return '\n'.join(output_lines)
    
    def detect_and_map_two_char_chains(self, pdb_content):
        """Detect two-character chain IDs and create mapping to single characters"""
        lines = pdb_content.strip().split('\n')
        two_char_chains = set()
        existing_single_chains = set()
        
        # Process lines to fix atom overflow first, then detect chains
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Apply atom overflow fix
                temp_line = line
                if len(temp_line) > 13:
                    try:
                        atom_num_field = temp_line[6:12].strip()
                        if atom_num_field.isdigit() and int(atom_num_field) >= 100000:
                            temp_line = temp_line[:13] + temp_line[14:]
                    except (ValueError, IndexError):
                        pass
                
                if len(temp_line) > 21:
                    # Check for single-character chain ID
                    if len(temp_line) > 22 and temp_line[22] == ' ':
                        single_chain = temp_line[21]
                        if single_chain != ' ' and single_chain.isalpha():
                            existing_single_chains.add(single_chain)
                    
                    # Check for two-character chain ID
                    elif len(temp_line) > 22 and temp_line[22] != ' ':
                        two_char_chain = temp_line[21:23]
                        if (temp_line[21] != ' ' and 
                            temp_line[21].isalpha() and temp_line[22].isalpha()):
                            two_char_chains.add(two_char_chain)
        
        # Create mapping
        available_names = list(string.ascii_uppercase)
        for existing_chain in existing_single_chains:
            if existing_chain in available_names:
                available_names.remove(existing_chain)
        
        chain_mapping = {}
        for two_char_chain in sorted(two_char_chains):
            if available_names:
                single_char = available_names.pop(0)
                chain_mapping[two_char_chain] = single_char
        
        return chain_mapping
    
    def correct_pdb_line(self, line, chain_mapping):
        """Apply both atom overflow and chain ID corrections"""
        if not line.startswith(('ATOM', 'HETATM')):
            return line
        
        if len(line) <= 22:
            return line
        
        # Fix atom number overflow (≥100000)
        if len(line) > 13:
            try:
                atom_num_field = line[6:12].strip()
                if atom_num_field.isdigit() and int(atom_num_field) >= 100000:
                    line = line[:13] + line[14:]
            except (ValueError, IndexError):
                pass
        
        # Fix two-letter chain IDs
        if len(line) > 22 and line[22] != ' ' and len(line) > 21:
            two_char_chain = line[21:23]
            if (line[21] != ' ' and 
                line[21].isalpha() and line[22].isalpha() and 
                two_char_chain in chain_mapping):
                return line[:21] + chain_mapping[two_char_chain] + ' ' + line[23:]
        
        return line
    
    def apply_corrections(self, pdb_content):
        """Apply all corrections to PDB content"""
        # Detect chain mappings
        chain_mapping = self.detect_and_map_two_char_chains(pdb_content)
        
        # Apply corrections line by line
        lines = pdb_content.strip().split('\n')
        corrected_lines = []
        corrections_made = 0
        
        for line in lines:
            original_line = line
            corrected_line = self.correct_pdb_line(line, chain_mapping)
            
            if original_line != corrected_line:
                corrections_made += 1
            
            corrected_lines.append(corrected_line)
        
        return '\n'.join(corrected_lines), corrections_made, len(chain_mapping)
    
    def process_single_pdb(self, pdb_id):
        """Process a single PDB through the entire pipeline"""
        start_time = time.time()
        
        # Check if already processed
        output_file = self.output_dir / f"{pdb_id}.pdb"
        if output_file.exists():
            self.stats['already_exists'] += 1
            return True, "Already exists"
        
        try:
            # Step 1: Download
            self.log(f"Processing {pdb_id}: Downloading...")
            pdb_content, file_format = self.download_pdb(pdb_id)
            
            if pdb_content is None:
                self.stats['failed'] += 1
                return False, "Download failed"
            
            self.stats['downloaded'] += 1
            
            # Step 2: Find RNA chains
            if file_format == 'cif':
                rna_chains = self.find_rna_chains_cif(pdb_content)
            else:
                rna_chains = self.find_rna_chains_pdb(pdb_content)
            
            if not rna_chains:
                self.stats['failed'] += 1
                return False, "No RNA chains found"
            
            self.log(f"Processing {pdb_id}: Found RNA chains: {sorted(rna_chains)}")
            
            # Step 3: Extract RNA chains
            if file_format == 'cif':
                rna_content = self.extract_rna_chains_cif(pdb_content, rna_chains)
            else:
                rna_content = self.extract_rna_chains_pdb(pdb_content, rna_chains)
            
            self.stats['rna_found'] += 1
            
            # Step 4: Apply corrections
            final_content, corrections_made, chain_mappings = self.apply_corrections(rna_content)
            
            # Save final result
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(final_content)
            
            self.stats['corrected'] += 1
            
            elapsed = time.time() - start_time
            self.log(f"Processing {pdb_id}: ✓ Complete ({corrections_made} corrections, {chain_mappings} chain mappings, {elapsed:.1f}s)")
            
            return True, f"Success ({corrections_made} corrections)"
            
        except Exception as e:
            self.stats['failed'] += 1
            self.log(f"Processing {pdb_id}: ✗ Error: {e}", "ERROR")
            return False, f"Error: {e}"
    
    def run_pipeline(self, pdb_list=None, max_files=None):
        """Run the complete pipeline"""
        self.log("🧬 RNA Processing Pipeline Started")
        self.log("=" * 50)
        
        start_time = time.time()
        
        # Get PDB list
        if pdb_list is None:
            pdb_list = self.get_rna_structure_list()
        
        if not pdb_list:
            self.log("No PDB structures to process", "ERROR")
            return
        
        if max_files:
            pdb_list = pdb_list[:max_files]
            self.log(f"Limited to first {max_files} files for testing")
        
        self.stats['total'] = len(pdb_list)
        self.log(f"Processing {len(pdb_list)} PDB structures...")
        
        # Process each PDB
        for i, pdb_id in enumerate(pdb_list, 1):
            success, message = self.process_single_pdb(pdb_id)
            
            # Progress update every 100 files
            if i % 100 == 0:
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                eta = (self.stats['total'] - i) / rate if rate > 0 else 0
                
                self.log(f"Progress: {i}/{self.stats['total']} | "
                        f"Success: {self.stats['corrected']} | "
                        f"Failed: {self.stats['failed']} | "
                        f"Rate: {rate:.1f}/s | "
                        f"ETA: {eta/60:.1f} min")
            
            # Rate limiting
            time.sleep(self.delay)
        
        # Final summary
        total_time = time.time() - start_time
        self.log("\n🎉 PIPELINE COMPLETE!")
        self.log("=" * 50)
        self.log(f"Total processed: {self.stats['total']}")
        self.log(f"Downloaded: {self.stats['downloaded']}")
        self.log(f"RNA found: {self.stats['rna_found']}")
        self.log(f"Successfully corrected: {self.stats['corrected']}")
        self.log(f"Already existed: {self.stats['already_exists']}")
        self.log(f"Failed: {self.stats['failed']}")
        self.log(f"Total time: {total_time/60:.1f} minutes")
        self.log(f"Average rate: {self.stats['total']/total_time:.1f} files/second")
        self.log(f"Output directory: {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(description="RNA Processing Pipeline - Space-Efficient")
    parser.add_argument('--output-dir', default='rna_final', help='Output directory for final files')
    parser.add_argument('--log-dir', default='logs', help='Log directory')
    parser.add_argument('--pdb-list', help='File containing PDB IDs to process (one per line)')
    parser.add_argument('--max-files', type=int, help='Maximum number of files to process (for testing)')
    parser.add_argument('--delay', type=float, default=0.1, help='Delay between downloads (seconds)')
    parser.add_argument('--timeout', type=int, default=30, help='Request timeout (seconds)')
    
    args = parser.parse_args()
    
    # Create pipeline
    pipeline = RNAPipeline(
        output_dir=args.output_dir,
        log_dir=args.log_dir,
        delay=args.delay,
        timeout=args.timeout
    )
    
    # Load PDB list if provided
    pdb_list = None
    if args.pdb_list:
        try:
            with open(args.pdb_list, 'r') as f:
                pdb_list = [line.strip().upper() for line in f if line.strip()]
            pipeline.log(f"Loaded {len(pdb_list)} PDB IDs from {args.pdb_list}")
        except Exception as e:
            pipeline.log(f"Error loading PDB list: {e}", "ERROR")
            return
    
    # Run pipeline
    pipeline.run_pipeline(pdb_list=pdb_list, max_files=args.max_files)

if __name__ == "__main__":
    main() 