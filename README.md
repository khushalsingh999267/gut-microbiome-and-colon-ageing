# Gut Microbiome and Colon Ageing: A Transcriptomic Study (GSE278548)

This repository contains a professional-grade bioinformatics pipeline designed to investigate the relative impact of microbial colonization and biological ageing on the colon's transcriptomic landscape.

Using RNA-sequencing (RNA-seq) analysis on the **GSE278548** dataset, we provide a detailed map of how the presence of bacteria (Conventional) versus their absence (Germ-Free) reshapes the biological processes of the colon across different age groups.

---

## Project Overview

The primary objective of this study is to identify the **transcriptomic signatures** of microbial colonization. Key research questions include:
1.  **Microbial Impact:** How does the presence of a healthy microbiome alter gene expression compared to sterile (Germ-Free) environments?
2.  **Ageing vs. Microbiome:** Does biological ageing or microbial status have a more profound effect on the colon's molecular profile?
3.  **Metabolic Regulation:** Which specific metabolic pathways are regulated by the gut microbiome?

---

## Pipeline & Methodology

The analysis is implemented as a modular and reproducible pipeline:

*   **Data Acquisition:** Automated retrieval of raw count matrices from NCBI GEO.
*   **Normalization:** Conversion of raw counts to **Counts Per Million (CPM)** for cross-sample comparability.
*   **Differential Expression Analysis:** Utilizing **PyDESeq2** to calculate Fold Change and statistical significance (P-values).
*   **Pathway Enrichment:** Conducting **GSEA (Gene Set Enrichment Analysis)** to identify biological processes, such as Butanoate metabolism and inflammatory responses.
*   **Workflow Orchestration:** Managed by **Snakemake** to ensure reproducibility and scalability.

---

## Results & Biological Insights

### 1. Technical Quality Control
Library size analysis confirms consistent sequencing depth across all experimental groups (Young, Old, Conventional, and Germ-Free), ensuring that biological comparisons are not biased by technical variation.
![Library Size](results/figures/01_qc_library_size.png)

### 2. Global Transcriptomic Variance (PCoA)
Principal Coordinate Analysis (PCoA) reveals that microbial status (Conventional vs. Germ-Free) is the dominant driver of transcriptomic variation, accounting for significantly more variance than biological age.
![PCoA Discovery](results/figures/02_pcoa_discovery.png)

### 3. Differential Gene Expression
Volcano plots highlight key biomarkers of microbial colonization, including antimicrobial peptides like *Reg3b* and *Reg3g*, which are upregulated in the presence of bacteria.
![Volcano Plot](results/figures/03_volcano_mechanisms.png)

### 4. Functional Pathway Enrichment
GSEA identifies **Butanoate (Butyrate) Metabolism** as a highly enriched pathway. As butyrate is the primary energy source for colonocytes, its absence in Germ-Free mice suggests a state of metabolic deprivation.
![Pathway Enrichment](results/figures/04_pathway_interpretation.png)

### 5. Group Cluster Boundaries
95% confidence ellipses demonstrate distinct biological states between Conventional and Germ-Free mice, with no overlap, confirming the fundamental reprogramming of the colon by the microbiome.
![Group Boundaries](results/figures/05_pcoa_boundaries.png)

---

## Conclusion

This study demonstrates that the **gut microbiome is a primary determinant of colon health**, often overshadowing the effects of chronological ageing. The findings emphasize the importance of microbial-derived metabolites, particularly butyrate, in maintaining colonic homeostasis and suggest that targeting microbial pathways may be a viable strategy for mitigating age-related gut decline.

---

## Getting Started

### Prerequisites
Install the required Python packages:
```bash
pip install biopython streamlit plotly pandas requests pyyaml snakemake scikit-bio pydeseq2 adjustText
```

### Usage
1.  **Execute the Analysis Pipeline:**
    ```bash
    python3 main.py
    ```
2.  **Launch the Interactive Dashboard:**
    ```bash
    streamlit run streamlit_app.py
    ```

---
*Project maintained by Khushal Singh.*
