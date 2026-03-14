import pandas as pd
from pydeseq2.preprocessing import deseq2_norm
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import numpy as np
import requests

def get_mouse_gene_map():
    print("Fetching gene mapping (Ensembl -> Symbol)...")
    url = "https://raw.githubusercontent.com/dpryan79/Answers/master/phylogeny/ensembl_to_symbol.mouse.txt"
    try:
        mapping = pd.read_csv(url, sep='\t', header=None, names=['Ensembl', 'Symbol'])
        return dict(zip(mapping['Ensembl'], mapping['Symbol']))
    except:
        return {}

def run_deseq2(adata, design_factor='Microbiome'):
    print(f"Running DESeq2 by {design_factor}...")
    # Use integer counts
    counts_df = pd.DataFrame(adata.X.astype(int), index=adata.obs_names, columns=adata.var_names)
    
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=adata.obs,
        design_factors=design_factor
    )
    dds.deseq2()
    
    # Contrast: [factor, treated, control]
    stat_res = DeseqStats(dds, contrast=[design_factor, 'Conv', 'GF'])
    stat_res.summary()
    
    # In pydeseq2 0.5.x, results are accessed via stat_res.results_df
    # If that fails, try accessing the attribute after summary()
    try:
        res_df = stat_res.results_df
    except AttributeError:
        # Fallback for different versions
        print("Attribute results_df not found, attempting fallback...")
        res_df = stat_res.summary() 
    
    gene_map = get_mouse_gene_map()
    res_df['Symbol'] = res_df.index.map(lambda x: gene_map.get(x, x))
    
    return res_df

def run_gsea(de_results, gene_set='KEGG_2019_Mouse'):
    print(f"Running GSEA using {gene_set}...")
    if 'Symbol' not in de_results.columns:
        print("No Symbol column found for GSEA.")
        return pd.DataFrame()
        
    ranking = de_results[['Symbol', 'log2FoldChange']].dropna()
    ranking['Symbol'] = ranking['Symbol'].str.upper()
    ranking = ranking.sort_values(by='log2FoldChange', ascending=False)
    ranking = ranking.groupby('Symbol').first()
    
    try:
        pre_res = gp.prerank(rnk=ranking, 
                             gene_sets=gene_set,
                             threads=4,
                             min_size=5,
                             max_size=500,
                             permutation_num=100, 
                             outdir=None, 
                             seed=42)
        return pre_res.res2d
    except Exception as e:
        print(f"GSEA failed: {e}")
        return pd.DataFrame()
