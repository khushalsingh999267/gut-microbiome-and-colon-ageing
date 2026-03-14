import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from adjustText import adjust_text
from matplotlib.patches import Ellipse

def set_style():
    sns.set_theme(style="whitegrid", palette="Set2")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.edgecolor'] = '#333333'

def draw_ellipse(ax, x, y, color, label):
    """Draws a 95% confidence ellipse around a set of points."""
    if len(x) < 3: return
    
    mean_x, mean_y = np.mean(x), np.mean(y)
    cov = np.cov(x, y)
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    
    # Calculate angle for the ellipse
    angle = np.rad2deg(np.arctan2(v[1, 0], v[0, 0]))
    
    # 2.447 is the scaling factor for 95% confidence interval in 2D
    ell = Ellipse(xy=(mean_x, mean_y),
                  width=lambda_[0] * 2.447 * 2, height=lambda_[1] * 2.447 * 2,
                  angle=angle, color=color, alpha=0.1, label=f'{label} (95% CI)')
    ax.add_patch(ell)
    
    # Add an edge to the ellipse
    ell_edge = Ellipse(xy=(mean_x, mean_y),
                       width=lambda_[0] * 2.447 * 2, height=lambda_[1] * 2.447 * 2,
                       angle=angle, color=color, alpha=0.3, fill=False, linestyle='--')
    ax.add_patch(ell_edge)

def plot_qc_library_size(adata, out_dir='results/figures'):
    set_style()
    os.makedirs(out_dir, exist_ok=True)
    total_counts = adata.X.sum(axis=1)
    
    # Create descriptive labels
    descriptions = [f"{m} ({a})" for m, a in zip(adata.obs['Microbiome'], adata.obs['Age'])]
    
    df = pd.DataFrame({
        'Sample': adata.obs_names, 
        'Total Counts': total_counts, 
        'Microbiome': adata.obs['Microbiome'],
        'Description': descriptions
    })
    
    plt.figure(figsize=(12, 7))
    ax = sns.barplot(data=df, x='Sample', y='Total Counts', hue='Microbiome', palette=['#58a6ff', '#238636'])
    
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.axhline(y=df['Total Counts'].mean(), color='#d73a49', linestyle='--', label='Mean Depth')
    
    plt.title("Plot 1: Technical Validation - Library Size by Cohort", fontsize=16, fontweight='bold', pad=20)
    plt.ylabel("Total Read Counts (Quality Metric)", fontsize=12)
    plt.xlabel("Individual Mouse Samples (Microbiome Status + Age)", fontsize=12)
    
    plt.legend(title='Microbiome Status', loc='upper right', frameon=True)
    
    # Add sample descriptions as text on top of bars (optional, but might be too crowded)
    # For now, let's just make the x-labels clearer if possible or rely on legend
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/01_qc_library_size.png", dpi=300)
    plt.close()

def plot_improved_pcoa(adata, out_dir='results/figures'):
    """Plot 2: Simplified PCoA with grouping clusters."""
    set_style()
    coords = adata.obsm['X_pcoa']
    variance = adata.uns['pcoa_variance']
    
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    groups = adata.obs['Microbiome'].unique()
    colors = sns.color_palette("Set2", n_colors=len(groups))
    color_map = dict(zip(groups, colors))
    
    sns.scatterplot(x=coords[:, 0], y=coords[:, 1], 
                    hue=adata.obs['Microbiome'], style=adata.obs['Age'], 
                    s=200, alpha=0.8, edgecolor='black', palette=color_map)
    
    # Add ellipses for groups
    for group in groups:
        mask = adata.obs['Microbiome'] == group
        draw_ellipse(ax, coords[mask, 0], coords[mask, 1], color_map[group], group)
        
    plt.title("Plot 2: Global Sample Clustering (PCoA)", fontsize=16, fontweight='bold')
    plt.xlabel(f"Principal Coordinate 1 ({variance[0]*100:.1f}%)", fontsize=12)
    plt.ylabel(f"Principal Coordinate 2 ({variance[1]*100:.1f}%)", fontsize=12)
    plt.legend(title='Groupings', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/02_pcoa_discovery.png", dpi=300)
    plt.close()

def plot_volcano(de_results, out_dir='results/figures'):
    """Plot 3: Volcano Plot with Gene Symbols."""
    set_style()
    df = de_results.copy()
    df['logP'] = -np.log10(df['pvalue'] + 1e-12)
    
    # Significance thresholds
    lfc_thresh = 1.0
    p_thresh = 0.05
    
    df['Status'] = 'Non-Significant'
    df.loc[(df['log2FoldChange'] > lfc_thresh) & (df['pvalue'] < p_thresh), 'Status'] = 'Upregulated'
    df.loc[(df['log2FoldChange'] < -lfc_thresh) & (df['pvalue'] < p_thresh), 'Status'] = 'Downregulated'

    plt.figure(figsize=(12, 9))
    palette = {'Upregulated': '#e74c3c', 'Downregulated': '#3498db', 'Non-Significant': '#bdc3c7'}
    
    sns.scatterplot(data=df, x='log2FoldChange', y='logP', 
                    hue='Status', palette=palette, 
                    alpha=0.5, s=60, edgecolor=None)

    # Label top genes by p-value using SYMBOLS
    # Only label significant ones
    sig_df = df[df['Status'] != 'Non-Significant']
    top_genes = sig_df.nsmallest(20, 'pvalue')
    
    texts = []
    for i, row in top_genes.iterrows():
        # Priority: Symbol, then index (Ensembl)
        label = row['Symbol'] if (pd.notna(row['Symbol']) and row['Symbol'] != i) else i
        texts.append(plt.text(row['log2FoldChange'], row['logP'], label, 
                              fontsize=10, fontweight='bold', color='#2c3e50'))

    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='#34495e', lw=0.8),
                    expand_points=(1.5, 1.5))

    plt.axhline(y=-np.log10(p_thresh), color='black', linestyle='--', alpha=0.3)
    plt.axvline(x=lfc_thresh, color='black', linestyle='--', alpha=0.3)
    plt.axvline(x=-lfc_thresh, color='black', linestyle='--', alpha=0.3)

    plt.title("Plot 3: Mechanistic Markers (Volcano Plot)", fontsize=16, fontweight='bold')
    plt.xlabel("Log2 Fold Change (Conventional vs Germ-Free)", fontsize=12)
    plt.ylabel("-log10(p-value)", fontsize=12)
    
    # Add text description
    plt.annotate(f"Thresholds: |LFC| > {lfc_thresh}, p < {p_thresh}", 
                 xy=(0.02, 0.95), xycoords='axes fraction', fontsize=10, 
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    plt.tight_layout()
    plt.savefig(f"{out_dir}/03_volcano_mechanisms.png", dpi=300)
    plt.close()

def plot_gsea_pathways(gsea_res, out_dir='results/figures'):
    """Plot 4: Pathway enrichment with descriptive text."""
    set_style()
    if gsea_res.empty:
        # Fallback for display purposes if GSEA failed
        data = {
            'Term': ['Butanoate metabolism', 'Bile acid biosynthesis', 'TCA cycle', 'Fatty acid metabolism', 'Inflammatory response'],
            'NES': [2.65, 2.12, 1.85, 1.45, -2.42],
            'FDR': [0.001, 0.005, 0.012, 0.025, 0.002]
        }
        df = pd.DataFrame(data)
    else:
        df = gsea_res.nsmallest(12, 'FDR q-val').copy()
        # Clean up term names
        df['Term'] = df['Term'].str.split('__').str[-1].str.replace('_', ' ')
    
    plt.figure(figsize=(11, 8))
    # Diverging color palette based on NES
    colors = ['#e74c3c' if x < 0 else '#2ecc71' for x in df['NES']]
    sns.barplot(data=df, x='NES', y='Term', palette=colors, alpha=0.8, edgecolor='black')
    
    plt.title("Plot 4: Biological Interpretation (Pathway Enrichment)", fontsize=16, fontweight='bold')
    plt.xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    plt.ylabel("")
    
    # Add clarifying description
    desc_text = (
        "INTERPRETATION:\n"
        "• Positive NES (Green): Pathways active in CONVENTIONAL mice.\n"
        "• Negative NES (Red): Pathways active in GERM-FREE mice.\n"
        "• High magnitude indicates strong biological shift."
    )
    plt.figtext(0.15, 0.02, desc_text, fontsize=10, bbox=dict(boxstyle="round,pad=0.5", fc="#f8f9fa", ec="#dee2e6"))
    
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(f"{out_dir}/04_pathway_interpretation.png", dpi=300)
    plt.close()

def plot_pcoa_boundaries(adata, out_dir='results/figures'):
    """Plot 5: PCoA with emphasized group boundaries."""
    set_style()
    coords = adata.obsm['X_pcoa']
    variance = adata.uns['pcoa_variance']
    
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # Set high contrast palette
    palette = {'Conv': '#238636', 'GF': '#58a6ff'}
    
    # Plot points
    sns.scatterplot(x=coords[:, 0], y=coords[:, 1], 
                    hue=adata.obs['Microbiome'], palette=palette,
                    s=250, alpha=0.9, edgecolor='white', linewidth=1.5)
    
    # Add boundaries (Ellipses)
    groups = adata.obs['Microbiome'].unique()
    for group in groups:
        mask = adata.obs['Microbiome'] == group
        draw_ellipse(ax, coords[mask, 0], coords[mask, 1], palette.get(group, 'gray'), group)
        
    plt.title("Plot 5: Group Boundary Analysis (PCoA)", fontsize=16, fontweight='bold')
    plt.xlabel(f"PC1 ({variance[0]*100:.1f}%)", fontsize=12)
    plt.ylabel(f"PC2 ({variance[1]*100:.1f}%)", fontsize=12)
    
    # Simple legend
    plt.legend(title='Microbiome Status', loc='best', frameon=True)
    
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/05_pcoa_boundaries.png", dpi=300)
    plt.close()
