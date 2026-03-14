rule all:
    input:
        "results/figures/01_qc_library_size.png",
        "results/figures/02_pcoa_discovery.png",
        "results/figures/03_volcano_mechanisms.png",
        "results/figures/04_pathway_interpretation.png",
        "results/figures/05_pcoa_boundaries.png"

rule run_analysis:
    output:
        "results/figures/01_qc_library_size.png",
        "results/figures/02_pcoa_discovery.png",
        "results/figures/03_volcano_mechanisms.png",
        "results/figures/04_pathway_interpretation.png",
        "results/figures/05_pcoa_boundaries.png"
    script:
        "main.py"
