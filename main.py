import yaml
import argparse
import os
from src.io import fetch_geo_counts, generate_metadata_mock
from src.processing import process_counts, run_ordination
from src.stats import run_deseq2, run_gsea
from src.viz import (plot_qc_library_size, plot_improved_pcoa, 
                     plot_volcano, plot_gsea_pathways, plot_pcoa_boundaries)

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="RNA-seq Post-Quant Analysis Pipeline")
    parser.add_argument('--config', default='config/study_config.yaml', help='Path to YAML config')
    args = parser.parse_args()

    # Load parameters
    config = load_config(args.config)
    print(f"Starting analysis for study: {config['study_id']}")

    # 1. Fetch Data
    counts = fetch_geo_counts(config['study_id'])
    if counts is None:
        return
    
    # 2. Generate Metadata
    meta = generate_metadata_mock(counts.columns.tolist())

    # 3. Pre-processing & QC
    adata = process_counts(counts, meta, config)
    plot_qc_library_size(adata) # STORY PLOT 1
    
    # 4. Ordination (PCoA)
    adata = run_ordination(adata)
    plot_improved_pcoa(adata)   # STORY PLOT 2
    plot_pcoa_boundaries(adata) # STORY PLOT 5
    
    # 5. Differential Expression
    de_results = run_deseq2(adata)
    plot_volcano(de_results)    # STORY PLOT 3
    
    # 6. Pathway Enrichment
    gsea_res = run_gsea(de_results)
    plot_gsea_pathways(gsea_res) # STORY PLOT 4

    print("Pipeline finished successfully. Figures generated in results/figures/")

if __name__ == "__main__":
    main()
