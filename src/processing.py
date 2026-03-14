import pandas as pd
import numpy as np
import scanpy as sc
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform

def process_counts(counts_df, metadata_df, config):
    """
    Normalization, Filtering and HVG selection using Scanpy
    """
    print("Normalizing and Filtering...")
    # Create AnnData
    # counts_df is Gene x Sample, so we transpose to Sample x Gene
    adata = sc.AnnData(X=counts_df.values.T, obs=metadata_df, var=pd.DataFrame(index=counts_df.index))
    
    # Filter genes
    sc.pp.filter_genes(adata, min_counts=config['filtering_threshold'])
    
    # Normalize
    if config['normalization_method'] == 'CPM':
        sc.pp.normalize_total(adata, target_sum=1e6)
    
    sc.pp.log1p(adata)
    
    # HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=config['n_hvgs'])
    
    return adata

def run_ordination(adata, metric='braycurtis'):
    """
    Runs PCoA using scikit-bio
    """
    print(f"Running PCoA using {metric} metric...")
    # Extract data for HVGs
    hvg_data = adata[:, adata.var.highly_variable].X
    
    # Distance matrix (pdist handles 'braycurtis' as a string)
    dist_matrix = pdist(hvg_data, metric=metric)
    dist_matrix = squareform(dist_matrix)
    
    # PCoA
    pcoa_results = pcoa(dist_matrix)
    
    # Add coordinates back to adata.obsm
    adata.obsm['X_pcoa'] = pcoa_results.samples[['PC1', 'PC2']].values
    adata.uns['pcoa_variance'] = pcoa_results.proportion_explained[:2].values
    
    return adata
