# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:40:07 2024
@author: JingyiLi
"""

from flask import Flask, render_template, request, send_file
import sqlite3
import pandas as pd
import os
import zipfile
from io import BytesIO
from flask import send_from_directory
# import logging

# # Configure logging
# logging.basicConfig(
#     filename='app.log',        # Log file name
#     level=logging.DEBUG,       # Log level (DEBUG shows all details)
#     format='%(asctime)s [%(levelname)s] %(message)s'  # Log format
# )

# logger = logging.getLogger(__name__)

app = Flask(__name__)
app = Flask(__name__, static_folder='static', static_url_path='/static')

# Connect to SQLite database
def get_db_connection():
    base_dir = os.path.dirname(__file__)
    db_path = os.path.join(base_dir, 'ATLAS.db')
    conn = sqlite3.connect(db_path)
    return conn

# Search RNA motifs based on motif_type and nt_number or subtype
def search_motif(motif_type, nt_number=None):
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        if motif_type in ["hairpin", "internal", "bulge"]:
            # Query for hairpin, internal, and bulge from "data" table
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type,))
            rows = cursor.fetchall()
            # Filter rows by nt_number (length of nt_number list)
            filtered_rows = [row for row in rows if len(row[3].split(',')) == int(nt_number)]

        elif 'junction' in motif_type:
            # Query for multiway junctions
            motif_type_db = motif_type.replace("-way", "").replace(" ", "_")
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db,))
            filtered_rows = cursor.fetchall()

        elif 'pseudoknot' in motif_type:
            # Query for pseudoknot subtypes from "PK" table
            motif_type_db = motif_type
            query = "SELECT id, motif_type, pdbid, nt_number, file_content FROM PK WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db.strip(),))
            filtered_rows = cursor.fetchall()

        else:
            # If motif_type is unknown
            filtered_rows = []

        conn.close()

        # If results are found, save them and return
        if filtered_rows:
            # Save results to CSV
            df = pd.DataFrame(filtered_rows, columns=['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Filecontent'])
            base_dir = os.path.dirname(os.path.abspath(__file__))
            csv_path = os.path.join(base_dir, 'search_results.csv')
            df[['ID', 'Motif Type', 'PDB ID', 'NT Number']].to_csv(csv_path, index=False)

            # Save PDB files
            pdb_dir = os.path.join(base_dir, 'pdb_files')
            if os.path.exists(pdb_dir):
                for file in os.listdir(pdb_dir):
                    os.remove(os.path.join(pdb_dir, file))
            else:
                os.makedirs(pdb_dir)

            for row in filtered_rows:
                motif_id = row[0]
                pdb_content = row[4]  # Filecontent column
                pdb_filename = os.path.join(pdb_dir, f'{motif_id}.pdb')
                with open(pdb_filename, 'w') as pdb_file:
                    pdb_file.write(pdb_content)

            # Create a ZIP archive
            zip_buffer = BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for file in os.listdir(pdb_dir):
                    file_path = os.path.join(pdb_dir, file)
                    zip_file.write(file_path, arcname=file)
            zip_buffer.seek(0)

            return filtered_rows, csv_path, zip_buffer

        # No results found
        return [], None, None

    except sqlite3.Error as e:
        print(f"Error during search: {e}")
        return [], None, None

# Route to the home page
@app.route('/')
def home():
    return render_template('index.html')

# Route to handle motif search requests
@app.route('/search', methods=['POST'])
def search():
    motif_type = request.form.get('motif_type')
    nt_number = request.form.get('nt_number')
    if motif_type == 'multiway-junction':
        # Here nt_number will be a string like '3-way junction'
        motif_type = nt_number
        nt_number_count = None
    elif motif_type == 'pseudoknot':
        motif_type = f'pseudoknot_{nt_number}'
        nt_number_count = None
    else:
        # For other types, it's a number, so convert to integer
        try:
            nt_number_count = int(nt_number)
        except ValueError:
            nt_number_count = None

    search_results, csv_path, pdb_zip_buffer = search_motif(motif_type, nt_number_count)
    
    if search_results:
        return render_template('results.html', search_results=search_results, csv_path=csv_path)
    else:
        return render_template('no_results.html')
    
# Route to download CSV file
@app.route('/download_csv', methods=['GET'])
def download_csv():
    csv_path = request.args.get('csv_path')

    if not csv_path or not os.path.exists(csv_path):
        return "CSV file not found.", 404

    try:
        return send_file(csv_path, mimetype='text/csv', as_attachment=True, download_name='search_results.csv')
    except Exception as e:
        return f"An error occurred during CSV download: {str(e)}", 500

# Route to download ZIP of PDB file
@app.route('/download_zip', methods=['GET'])
def download_zip():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_dir = os.path.join(base_dir, 'pdb_files')

    if not os.path.exists(pdb_dir) or not os.listdir(pdb_dir):
        return "PDB files not found.", 404

    # Create buffer
    zip_buffer = BytesIO()

    # Zip the files
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for root, dirs, files in os.walk(pdb_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, start=pdb_dir)
                zip_file.write(file_path, arcname=arcname)

    zip_buffer.seek(0)

    return send_file(
        zip_buffer,
        mimetype='application/zip',
        as_attachment=True,
        download_name='pdb_files.zip'
    )

@app.route('/download-database', methods=['GET'])
def download_database():
    # Path to the 'data' folder
    root_folder = '.'
    file_name = "ATLAS.db"  # Replace with your database file name

    # Use send_from_directory to serve the file
    return send_from_directory(directory=root_folder, path=file_name, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
