import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from src.io import fetch_geo_counts, generate_metadata_mock, get_mouse_gene_map
from src.processing import process_counts
import os

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="RNA-Seq Explorer | Microbiome-Ageing Axis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- MINIMALIST DARK THEME CSS ---
st.markdown("""
    <style>
    .stApp { background-color: #0e1117; color: #e0e0e0; }
    section[data-testid="stSidebar"] { background-color: #161b22; border-right: 1px solid #30363d; }
    h1, h2, h3 { color: #ffffff; font-weight: 700; }
    .metric-card {
        background-color: #161b22;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #30363d;
        text-align: center;
    }
    .plot-container {
        border: 1px solid #30363d;
        border-radius: 10px;
        padding: 10px;
        margin-bottom: 20px;
        background-color: #0d1117;
    }
    </style>
    """, unsafe_allow_html=True)

# --- SIDEBAR & CONFIG ---
with st.sidebar:
    st.image("https://img.icons8.com/ios-filled/100/ffffff/dna-helix.png", width=60)
    st.title("BLOCK-O-MICS")
    st.markdown("---")
    
    geo_id = st.text_input("GEO Accession", value="GSE278548")
    tissue = st.selectbox("Tissue Target", ["Colon", "Ileum", "Lung"], index=0)
    
    st.markdown("### Interactive Controls")
    norm_method = st.radio("Normalization", ["CPM", "Log1p(CPM)", "Raw"])
    
    st.markdown("---")
    if st.button("🔄 Refresh Data"):
        st.cache_data.clear()
        st.rerun()

# --- DATA LOADING ---
@st.cache_data
def load_data(accession, tissue_type):
    try:
        counts = fetch_geo_counts(accession, tissue=tissue_type)
        if counts is not None:
            meta = generate_metadata_mock(counts.columns.tolist())
            gene_map = get_mouse_gene_map()
            return counts, meta, gene_map
    except Exception as e:
        st.error(f"Error loading data: {e}")
    return None, None, None

# --- MAIN PAGE ---
st.title("🧬 RNA-Seq Post-Quant Explorer")
st.markdown(f"**Current View:** {geo_id} | {tissue} Analysis Pipeline")

counts, meta, gene_map = load_data(geo_id, tissue)

if counts is not None:
    # --- METRICS ---
    m1, m2, m3, m4 = st.columns(4)
    with m1: st.markdown(f'<div class="metric-card"><b>Samples</b><br><h2>{len(counts.columns)}</h2></div>', unsafe_allow_html=True)
    with m2: st.markdown(f'<div class="metric-card"><b>Total Genes</b><br><h2>{len(counts):,}</h2></div>', unsafe_allow_html=True)
    with m3: st.markdown(f'<div class="metric-card"><b>Condition</b><br><h2>{tissue}</h2></div>', unsafe_allow_html=True)
    with m4: st.markdown(f'<div class="metric-card"><b>Status</b><br><h2>Ready</h2></div>', unsafe_allow_html=True)

    st.markdown("---")

    # --- TABS ---
    tab_report, tab_interactive = st.tabs(["📄 Publication Report", "🔍 Interactive Data"])

    with tab_report:
        st.header("Analytical Story & Findings")
        st.info("The following plots are generated using the high-end Matplotlib pipeline for maximum scientific accuracy.")
        
        # Display the 5 plots we generated
        plot_files = [
            ("1. Technical Validation (Library Size)", "results/figures/01_qc_library_size.png", "Uniform sequencing depth across all cohorts ensures unbiased comparison."),
            ("2. Global Discovery (PCoA)", "results/figures/02_pcoa_discovery.png", "Principal Coordinate Analysis reveals the microbiome as the dominant driver of transcriptomic variance."),
            ("3. Mechanistic Drivers (Volcano)", "results/figures/03_volcano_mechanisms.png", "Identification of key biomarkers (e.g., Reg3g, Reg3b) differentially expressed in colonized colons."),
            ("4. Pathway Enrichment", "results/figures/04_pathway_interpretation.png", "Functional analysis highlights Butanoate Metabolism as a critical microbial-driven process."),
            ("5. Group Boundary Analysis", "results/figures/05_pcoa_boundaries.png", "Emphasized clustering showing distinct biological states between Conventional and Germ-Free mice.")
        ]
        
        for title, path, desc in plot_files:
            with st.container():
                st.subheader(title)
                if os.path.exists(path):
                    st.image(path, width="stretch")
                else:
                    st.warning(f"Plot not found at {path}. Please run `python3 main.py` to generate it.")
                st.markdown(f"**Interpretation:** {desc}")
                st.markdown("<br>", unsafe_allow_html=True)

    with tab_interactive:
        st.header("Real-time Data Exploration")
        
        # Apply normalization based on parameter
        processed_counts = counts.copy()
        if norm_method == "CPM":
            processed_counts = (counts / counts.sum()) * 1e6
        elif norm_method == "Log1p(CPM)":
            processed_counts = np.log1p((counts / counts.sum()) * 1e6)
        
        col_left, col_right = st.columns([1, 2])
        
        with col_left:
            st.markdown("### Metadata")
            st.dataframe(meta, width="stretch")
            
        with col_right:
            st.markdown(f"### Normalization: {norm_method}")
            # Show interactive bar chart for library sizes
            lib_sizes = counts.sum().reset_index()
            lib_sizes.columns = ['Sample', 'Total Counts']
            lib_sizes = lib_sizes.merge(meta, left_on='Sample', right_index=True)
            
            fig = px.bar(lib_sizes, x='Sample', y='Total Counts', color='Microbiome',
                         template="plotly_dark", color_discrete_sequence=['#58a6ff', '#238636'])
            fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', margin=dict(t=10, b=10))
            st.plotly_chart(fig, width="stretch")

        st.markdown("---")
        st.markdown("### Top Variable Genes (Live Heatmap)")
        # Show top 25 most variable genes
        top_var_genes = processed_counts.var(axis=1).sort_values(ascending=False).head(25).index
        heatmap_df = processed_counts.loc[top_var_genes]
        # Map to symbols
        heatmap_df.index = [gene_map.get(g, g) for g in heatmap_df.index]
        
        fig_heat = px.imshow(heatmap_df, 
                            color_continuous_scale='Viridis',
                            template="plotly_dark",
                            aspect="auto")
        fig_heat.update_layout(margin=dict(t=10, b=10))
        st.plotly_chart(fig_heat, width="stretch")

else:
    st.error("Failed to initialize data. Please verify your internet connection and GEO ID.")

st.markdown("---")
