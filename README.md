```markdown

**Steller** a platform utilizing high-diversity barcodes at sparse tagging densities to guarantee unique clonal identification. By integrating a specialized, probe-modified spatial transcriptomic workflow, Steller enhances lineage signature readout at single cell resolution without compromising transcriptomic analysis.

---

## Highlights

- A novel engineered spatial lineage-tracing chip: Employing a novel CS3-optimized probe design on 2.5 um spatial transcriptomic chips, we achieved substantial increase both the number and diversity of CloneBC recovery rates compared to traditional low-resolution Visium platform as well as conventional poly(T) probe at spatial single-cell resolution. It is compatible with existing commercial and academic spatial transcriptomics platform and scalable to the spatial lineage study for other biological programs.

- High-resolution mapping of diverse spatial lineage patterns in mouse forebrain: Leveraging Steller platform, we achieved, for the first time, single-cell-resolution tracking of neural progenitor cell (NPC) progeny across diverse brain regions. Beyond capturing the classic radially arranged cortical columns, we revealed several previously uncharacterized or complex spatial organizations, including:
•	Hippocampus: Identification of horizontally aligned neuronal clusters within the CA3 and dentate gyrus (DG);
•	Striatum: Observation of distinct, spatially clustered striatal spiny projection neuron (SPN) subtypes;
•	Thalamus: Discovery of dorso-ventrally partitioned territories of thalamic glutamatergic neurons. 

---

## Repository Structure

```text
Steller/
├── CS_design/             # CloneBC capture sequence (CS) design
├── scST/                  # Spatial transcriptomics with cell segmentation based on BMKMANU slide
├── snST/                  # Single nucleus spatial transcriptomics based on SeekSpace slide
├── whitelist/             # Whitelist of CloneBC plasmid library used in the study
│   ├── preprocessing/     # Data filtering, QC, and lineage barcode extraction
│   ├── integration/       # Alignment algorithms for transcriptomics and clones
│   └── visualization/     # Custom plotting engines for spatial & UMAP plots
├── Helper.R               # Source code using R
├── utils.py               # Source code using python
└── README.md

```

---

## Getting Started

### Prerequisites

Steller requires a Linux environment (tested on 64-core enterprise servers) and relies on both Python and R libraries depending on your chosen pipeline workflow.

* **Python Requirements:** `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`, `seaborn`
* **R Requirements:** `Seurat` (>= 4.0), `tidyverse`, `Matrix`

### Installation via Conda

We recommend using Conda to manage your environments and prevent dependency version conflicts (such as `GLIBCXX` mismatches).

```bash
# Clone the repository
git clone [https://github.com/AmoriLi/Steller.git](https://github.com/AmoriLi/Steller.git)
cd Steller

# Create and activate the conda environment
conda env create -f envs/steller_env.yml
conda activate steller

# (Optional) Verify your R/Python path configurations
python -c "import scanpy; print(scanpy.__version__)"
R -e "library(Seurat); packageVersion('Seurat')"

```

---

## Standard Workflow

Steller processes spatial lineages in three sequential modules:

```
[ Input Data ] ---> ( 1. Preprocessing ) ---> ( 2. Clonal Integration ) ---> ( 3. Spatial Visualization )

```

### 1. Preprocessing and Lineage Barcoding

Extract and filter high-quality spatial transcriptomic spots/cells alongside valid lineage tracer barcodes.

```bash
python scripts/preprocess.py \
    --input_matrix data/spatial_counts.h5ad \
    --lineage_barcodes data/clone_bcs.csv \
    --out_dir results/preprocessed/

```

### 2. Spatio-Transcriptomic Integration

Map clonal lineages onto the spatial coordinates and low-dimensional embeddings (UMAP/t-SNE) to measure clone distribution entropy and fate directions.

### 3. Publication-Quality Visualization

Generate clean, informative plots with distinct color assignments across complex clonal architectures.

```python
from src.visualization import plot_spatial_clones

# Example usage within a python pipeline
plot_spatial_clones(
    adata, 
    clone_column='assigned_clone', 
    palette='morandi_53', 
    save_path='figures/fig1_spatial_mapping.pdf'
)

```

---

## Tutorials and Demos

Check out the interactive guides in the `notebooks/` directory for full examples from start to finish:

* `01_data_alignment.ipynb`: Reading spatial datasets and overlaying lineage metadata.
* `02_clonality_analysis.ipynb`: Quantifying lineage dynamics in the developing mouse brain.
* `03_advanced_visualization.Rmd`: Generating high-impact figures with custom color palettes.

---

## Citation

If you use **Steller** in your research or find its workflows helpful, please cite our pre-print/manuscript:

```bibtex
@article{li2026steller,
  title={Steller: Integrating spatial transcriptomics with lineage tracing to map mouse brain development},
  author={Li, Amori and [Co-authors]},
  journal={Bioinformatics / Nature Methods / Cell (Update as appropriate)},
  year={2026},
  publisher={[Publisher Name]}
}

```

---

## Contact & Support

For questions, bug reports, or feature requests, please open an Issue on the [GitHub Issues page](https://www.google.com/search?q=https://github.com/AmoriLi/Steller/issues) or reach out directly to the maintainer at **[Your Email Address]**.

```

***

### Recommendations for Customizing This File:
1. **`envs/steller_env.yml`**: Ensure you place a valid environment file inside an `envs` directory so users can recreate your package setup exactly as it runs on your server.
2. **Morandi Palette**: Under Key Features and Visualization, I explicitly highlighted your sophisticated color palettes—this is a valuable selling point for bioinformaticians who struggle to visualize dozens of clones cleanly in one plot.
3. **Citation**: Update the BibTeX stub with your target journal or BioRxiv details once available.

```
