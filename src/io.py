import pandas as pd
import requests
import os
import gzip
import shutil
import io
import time
from Bio import Entrez

# Configure NCBI Entrez email
Entrez.email = "example@example.com"  # Please set your own email if needed

def fetch_geo_counts(accession, tissue='Colon', out_dir='data/raw'):
    """
    Downloads tissue-specific count matrices from GSE278548.
    """
    os.makedirs(out_dir, exist_ok=True)
    print(f"Fetching {tissue} counts for {accession}...")
    
    file_name = f"{accession}_counts_{tissue}.tsv.gz"
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file&file={file_name}"
    
    target_path = os.path.join(out_dir, file_name)
    extracted_path = target_path.replace('.gz', '')
    
    if not os.path.exists(extracted_path):
        print(f"  Downloading {file_name}...")
        try:
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(target_path, 'wb') as f:
                    f.write(response.content)
                
                print(f"  Extracting...")
                with gzip.open(target_path, 'rb') as f_in:
                    with open(extracted_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(target_path)
            else:
                print(f"  Error: Could not download file (Status: {response.status_code})")
                return None
        except Exception as e:
            print(f"  Error downloading {accession}: {e}")
            return None
    
    print(f"  Loading count matrix...")
    df = pd.read_csv(extracted_path, sep='\t', index_col=0)
    return df

def generate_metadata_mock(sample_ids):
    """
    Generates experimental design metadata based on GSE278548 sample structure.
    """
    print("Generating experimental design metadata...")
    meta = pd.DataFrame(index=sample_ids)
    
    # Study GSE278548 has 24 samples: Conv vs GF, Young vs Old
    # Pattern: 6 replicates per group
    microbiome = []
    age = []
    
    for i in range(len(sample_ids)):
        if i < 12:
            microbiome.append('Conv')
        else:
            microbiome.append('GF')
            
        if (i % 12) < 6:
            age.append('Young')
        else:
            age.append('Old')
            
    meta['Microbiome'] = microbiome
    meta['Age'] = age
    meta['Condition'] = meta['Microbiome'] + "_" + meta['Age']
    return meta

def fetch_geo_metadata_biopython(accession):
    """
    Uses Bio.Entrez to fetch real metadata for a GEO series.
    """
    print(f"Fetching Biopython-powered metadata for {accession}...")
    try:
        # Search for the GSE in gds database
        handle = Entrez.esearch(db="gds", term=f"{accession}[ACCN]", retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        ids = record["IdList"]
        if not ids:
            return None
            
        # Get details for these IDs
        handle = Entrez.esummary(db="gds", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        
        metadata_list = []
        for summary in summaries:
            # GSE records are often types like 'gse' or 'gsm'
            if summary['EntryType'] == 'GSM':
                # Parse sample characteristics
                chars = summary.get('ExtRelation', [])
                metadata_list.append({
                    'GSM': summary['Accession'],
                    'Title': summary['Title'],
                    'Summary': summary['Summary']
                })
        
        return pd.DataFrame(metadata_list)
    except Exception as e:
        print(f"Biopython metadata fetch error: {e}")
        return None

def get_mouse_gene_map():
    """
    Fetches a simple mapping from Ensembl ID to Gene Symbol.
    """
    print("Fetching gene mapping (Ensembl -> Symbol)...")
    url = "https://raw.githubusercontent.com/dpryan79/Answers/master/phylogeny/ensembl_to_symbol.mouse.txt"
    try:
        mapping = pd.read_csv(url, sep='\t', header=None, names=['Ensembl', 'Symbol'])
        return dict(zip(mapping['Ensembl'], mapping['Symbol']))
    except:
        return {}
