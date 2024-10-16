from flask import Flask, request, jsonify, render_template, url_for, redirect
import requests
from bs4 import BeautifulSoup
import xmltodict
from Bio import SeqIO 
import pytest
import xml.etree.ElementTree as ET
from io import StringIO
from unittest.mock import patch
import subprocess
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, UserMixin, login_user, logout_user, login_required, current_user
from flask_bcrypt import Bcrypt
import os
from datetime import datetime



app = Flask(__name__)

# BOLD API base URL
BOLD_API_URL = "http://v3.boldsystems.org/index.php/Ids_xml"
# Define the BLAST API endpoint
blast_endpoint = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
api_key = "4650e1877eb97708a2e62f181bc9749c5508"
# Define the BLAST API endpoint
blast_endpoint = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.config['SECRET_KEY'] = 'supersecretkey'

db = SQLAlchemy(app)
bcrypt = Bcrypt(app)

login_manager = LoginManager(app)
login_manager.login_view = 'login'

UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# In-memory storage for uploaded files
files = []

# User model
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(150), unique=True, nullable=False)
    password = db.Column(db.String(150), nullable=False)
    is_admin = db.Column(db.Boolean, default=False)

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

@app.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        hashed_password = bcrypt.generate_password_hash(password).decode('utf-8')

        new_user = User(username=username, password=hashed_password, is_admin=True)
        db.session.add(new_user)
        db.session.commit()

        return redirect(url_for('login'))
    return render_template('register.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        user = User.query.filter_by(username=username).first()

        if user and bcrypt.check_password_hash(user.password, password):
            login_user(user)
            return redirect(url_for('admin'))

        return "Invalid credentials"
    
    return render_template('login.html')

@app.route('/admin')
def admin():
    return render_template('admin_panel.html')

@app.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))

# Define a function to make a BLAST search request
def blast_search(query, database, program):
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": query
    }
    response = requests.post(blast_endpoint, params=params)
    return response.text

# Define a route for the BLAST search
@app.route("/blast", methods=["POST"])
def blast():
    query = request.form["queries"]
    database = request.form["database"]
    program = request.form["program"]
    result = blast_search(query, database, program)
    soup = BeautifulSoup(result, 'html.parser')
    results = soup.find('div', {'id': 'content'})
    return render_template('results.html', results=results)

# Define a route for the main page
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

# Define a route for Species Identification 
@app.route("/identifyspecies")
def identifyspecies():
    return render_template("identifyspecies.html")

@app.route('/blastsearch')
def blastsearch():
    return render_template('blast.html') 

@app.route('/restriction')
def restriction():
    return render_template('restriction.html')

# Define a route for the main page
@app.route("/about")
def about():
    return render_template("about.html")

# Route for uploading a file
@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']

    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file:
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], file.filename))
        upload_date = datetime.now().strftime('%Y-%m-%d')
        files.append({'filename': file.filename, 'date': upload_date})
        return redirect(url_for('admin_panel'))

# Route to delete a file
@app.route('/delete/<filename>', methods=['POST'])
def delete_file(filename):
    global files
    files = [f for f in files if f['filename'] != filename]

    # Delete the file from the server if needed
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if os.path.exists(file_path):
        os.remove(file_path)

    return redirect(url_for('admin_panel'))

# Function to query BOLD ID Engine API
def query_bold_id_engine(sequence):
    api_url = f"http://v3.boldsystems.org/index.php/Ids_xml?db=COX1_SPECIES_PUBLIC&sequence={sequence}"
    response = requests.get(api_url)
    
    # Check if the request was successful
    if response.status_code == 200:
        return parse_bold_response(response.text)
    else:
        return None

# Function to parse XML response from BOLD API
def parse_bold_response(xml_data):
    import xml.etree.ElementTree as ET

    results = []
    root = ET.fromstring(xml_data)
    
    for match in root.findall(".//match"):
        result = {
            'species': match.find('identification').text if match.find('identification') is not None else "Unknown",
            'similarity': match.find('similarity').text if match.find('similarity') is not None else "N/A",
            'specimen': {
                'publicpage': match.find('specimen_url').text if match.find('specimen_url') is not None else "N/A"
            }
        }
        results.append(result)
    
    return {'identifications': {'match': results}}

# Define the route that processes the DNA sequence and queries the BOLD API
@app.route('/identify', methods=['POST'])
def identify_sequence():
    dna_sequence = request.form.get('sequence')  # Get DNA sequence from the form
    print("DNA Sequence received:", dna_sequence)  # Debugging log
    
    if not dna_sequence:
        return jsonify({"error": "No DNA sequence provided"}), 400
    
    # Query the BOLD API using the DNA sequence
    result = query_bold_id_engine(dna_sequence)
    
    if result:
        print("BOLD API Result:", result)  # Log the API result for debugging
        
        # Extract relevant data from the response
        matches = result.get('identifications', {}).get('match', [])
        extracted_data = [
            {
                "species": match['species'],
                "similarity": match['similarity'],
                "specimen_page": match['specimen']['publicpage']
            }
            for match in matches
        ]
        
        return jsonify(extracted_data)
    else:
        print("Failed to retrieve data from BOLD API")  # Log failure for debugging
        return jsonify({"error": "Failed to retrieve data from BOLD API"}), 500


@app.route("/genedna")
def genedna():
    return render_template("genedna.html")

@app.route("/proteindata")
def proteindata():
    return render_template("protein.html")

@app.route("/submission")
def submission():
    return render_template("submission.html")

@app.route("/dnafulldata1")
def dnafulldata1():
    return render_template("dnafulldata1.html")

@app.route("/dnafulldata2")
def dnafulldata2():
    return render_template("dnafulldata2.html")

@app.route("/dnafulldata3")
def dnafulldata3():
    return render_template("dnafulldata3.html")

@app.route("/dnafulldata4")
def dnafulldata4():
    return render_template("dnafulldata4.html")

@app.route("/dnafulldata5")
def dnafulldata5():
    return render_template("dnafulldata5.html")

@app.route("/rnadata")
def rnadata():
    return render_template("rna.html")

if __name__ == "__main__":
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)