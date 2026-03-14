from flask import Flask, render_template, send_from_directory
import os

app = Flask(__name__)

# Route for the main report page
@app.route('/')
def index():
    # Read the report content to display it in the app
    report_path = 'REPORT.md'
    report_content = ""
    if os.path.exists(report_path):
        with open(report_path, 'r') as f:
            report_content = f.read()
    
    # List of figures to display
    figures = [
        {'title': '1. Technical Validation', 'path': 'images/01_qc_library_size.png', 'desc': 'Library size consistency across cohorts.'},
        {'title': '2. Global Ordinance', 'path': 'images/02_pcoa_discovery.png', 'desc': 'Bray-Curtis separation of biomes.'},
        {'title': '3. Mechanistic Drivers', 'path': 'images/03_volcano_mechanisms.png', 'desc': 'Top biomarkers of microbial colonization.'},
        {'title': '4. Functional Interpretation', 'path': 'images/04_pathway_interpretation.png', 'desc': 'Metabolic pathway enrichment scores.'},
        {'title': '5. Group Boundaries', 'path': 'images/05_pcoa_boundaries.png', 'desc': 'PCoA with emphasized group boundaries for better cluster interpretation.'}
    ]
    
    return render_template('index.html', report=report_content, figures=figures)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
