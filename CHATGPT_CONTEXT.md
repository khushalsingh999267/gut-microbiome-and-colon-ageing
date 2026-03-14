# 🧬 BLOCK-O-MICS: Project Reference for ChatGPT

This document provides a comprehensive overview of the **BLOCK-O-MICS** bioinformatics project. Use this as context for ChatGPT to help it understand the experiment, technical stack, and biological findings.

---

## 1. Project Overview
**Project Name:** BLOCK-O-MICS (The Microbiome-Ageing Axis in the Colon)  
**Study ID:** [GSE278548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278548) (Murine Colon Transcriptomics)  
**Core Question:** Does the gut microbiome (Conventional vs. Germ-Free) have a greater impact on colon health and gene expression than the natural ageing process (Young vs. Old)?

---

## 2. Experimental Design
The study uses a 2x2 factorial design with 24 samples:
- **Microbiome Status:**
  - **Conventional (Conv):** Mice with a normal gut microbiome.
  - **Germ-Free (GF):** Mice raised in sterile environments without any microbes.
- **Age Groups:**
  - **Young:** Early adulthood.
  - **Old:** Natural senescence.

---

## 3. Technical Architecture (The Pipeline)
The pipeline follows a modular "Unix Philosophy" and is orchestrated by **Snakemake**.

### Key Components:
- **Data Ingestion (`src/io.py`):** 
  - Fetches raw count matrices from NCBI GEO using `Biopython` and `Requests`.
  - Simulates metadata based on the experimental design.
- **Processing (`src/processing.py`):**
  - **Normalization:** Counts Per Million (CPM) + Log1p transformation.
  - **Ordination:** Principal Coordinate Analysis (PCoA) using Bray-Curtis dissimilarity (via `scikit-bio`).
- **Statistical Modeling (`src/stats.py`):**
  - **Differential Expression:** Uses **PyDESeq2** (industry standard) to calculate Fold Change and P-values.
  - **ID Mapping:** Translates Ensembl Gene IDs to human-readable Gene Symbols.
- **Pathway Analysis:**
  - **GSEA (Gene Set Enrichment Analysis):** Identifies biological processes (KEGG Pathways) significantly affected by the experimental variables.
- **Visualization (`src/viz.py`):**
  - Generates high-resolution (300 DPI) plots including Library Size QC, PCoA, Volcano Plots, and GSEA Pathway maps.
- **Deployment:**
  - **Streamlit (`streamlit_app.py`):** An interactive dashboard for real-time data exploration and visualization.

---

## 4. Key Biological Findings

### A. Microbiome > Ageing
PCoA analysis (Global Discovery) shows that **Microbiome Status** (Conventional vs. Germ-Free) is the primary driver of transcriptomic variance, separating samples much more distinctly than Age.

### B. Biomarkers of Colonization
Differential expression analysis identified key genes "switched on" by bacteria:
- **Antimicrobial Peptides:** *Reg3b*, *Reg3g* (essential for managing gut-microbe boundaries).
- **Immune Modulators:** Genes involved in the response to bacterial stimuli.

### C. Metabolic Starvation in GF Mice
GSEA revealed that **Butanoate (Butyrate) Metabolism** is the most enriched pathway in Conventional mice.
- **Insight:** Butyrate is a Short-Chain Fatty Acid (SCFA) produced by gut bacteria and is the primary energy source for colonocytes.
- **Conclusion:** Without a microbiome, the colon enters a state of **metabolic starvation**, which likely accelerates age-related decline.

---

## 5. Potential Discussion Points for ChatGPT
1. **Bioinformatics Strategy:** Ask about optimizing the PyDESeq2 parameters or alternative normalization methods like TMM or RLE.
2. **Biological Interpretation:** Discuss the role of *Reg3g* in gut homeostasis or how butyrate deficiency impacts the "Ageing Clock" of the colon.
3. **Pipeline Enhancement:** Explore adding Single-Cell RNA-seq (scRNA-seq) integration or Multi-omics (Metagenomics + Transcriptomics) correlation.
4. **Code Review:** Ask ChatGPT to analyze the modular structure of the `src/` directory or the efficiency of the `Snakemake` workflow.

---
*Created by Gemini CLI to assist in advanced bioinformatics research.*
