**Steller** a platform utilizing high-diversity barcodes at sparse tagging densities to guarantee unique clonal identification. By integrating a specialized, probe-modified spatial transcriptomic workflow, Steller enhances lineage signature readout at single cell resolution without compromising transcriptomic analysis.

---

## 🌟 Highlights

- **A novel engineered spatial lineage-tracing chip:**<br> Employing a novel CS3-optimized probe design on 2.5 um spatial transcriptomic chips, we achieved substantial increase both the number and diversity of CloneBC recovery rates compared to traditional low-resolution Visium platform as well as conventional poly(T) probe at spatial single-cell resolution. It is compatible with existing commercial and academic spatial transcriptomics platform and scalable to the spatial lineage study for other biological programs.

- **High-resolution mapping of diverse spatial lineage patterns in mouse forebrain:**<br> Leveraging Steller platform, we achieved, for the first time, single-cell-resolution tracking of neural progenitor cell (NPC) progeny across diverse brain regions. Beyond capturing the classic radially arranged cortical columns, we revealed several previously uncharacterized or complex spatial organizations, including:
<br>•	Hippocampus: Identification of horizontally aligned neuronal clusters within the CA3 and dentate gyrus (DG);
<br>•	Striatum: Observation of distinct, spatially clustered striatal spiny projection neuron (SPN) subtypes;
<br>•	Thalamus: Discovery of dorso-ventrally partitioned territories of thalamic glutamatergic neurons. 

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
[ 1. CS probe selection ] --->
( 2. CloneBC whitelist generation ) --->
( 3. Brain region definition ) --->
( 4. CloneBC extraction & CloneCalling ) --->
( 5. Cell type annotation ) --->
( 6. Spatial lineage fate analysis & visualization )

```
### 1. CS probe selection

#### 1a. CS probe selection
Firstly, we designed three types of CS probes for CloneBC targeted capture. Based on bulk amplicon dataset, we can identify which probe has the highest capture efficiency for CloneBC transcripts from transfected 4T1 total RNA.

```bash
python CS_design/Read1_CS123.py \
    CS_design/CS123/fastq #the directory path of test CS123 fastq data
```
You will see cs1, cs2, cs3 and other diresctory generated under the "CS_design/CS123/split_fastq" directory. There are readnames enriched by each probes (other means undetermined reads).

Then split original PE fastq file into subfiles of each probes.
```bash
sbatch -a 0-3 \ #specify slurm array range for running multiple tasks
    CS_design/2_splitfastq.sh \
    CS_design/CS123 #the directory path of test CS123
```
New PE fastq.gz are generated in cs1/cs2/cs3/other directories.

Finally, we need to determine the gene identity of each read sequence by running STAR mapping based on manual customized genome reference with mCherry-CloneBC information, which can be created via 'STAR_generate_index_v2.sh'.

```bash
sbatch STAR_generate_index_v2.sh \
    $(path of downloaded mouse_gencode_vM23) \
    $(new output directory of customized genome)

sbatch -a 0-3 \
    CS_design/3_mapping.sh \
    CS_design/CS123 #the directory path of test CS123
    $(path of customized genome reference) #star/mouse/P045_genocode_vM23/P045_mcherry_BC
```
Each probe's directory will contain a 'STAR' result, we can quantify and compare the read type and gene type captured by different probes.


#### 1b. CS3 probe quality control
Next, we modified CS3 probe on conventional polyT spatial beads, resulting about polyT:CS3 = 1:1. Based on bulk amplicon dataset, we can further compare whether CS3 probes show better performance than polyT on CloneBC transcripts enrichment from transfected 4T1 total RNA.

In order to identify the read type accurately, run fastp first for removing fuzzy reads which may be incorrectly assigned into polyT captured type due to very closed and similar 5'-end pattern of read1 sequence between polyT and CS3
<img width="704" height="121" alt="Screenshot 2026-06-03 at 10 43 18" src="https://github.com/user-attachments/assets/5a518153-8ed3-4222-b4e7-a01d85e35263" />

```bash
sbatch CS_design/0_fastp_v2.sh \
    CS_design/CS3_PT #the directory path of test CS3_polyT fastq data
```
'CS_design/CS3_PT/fastp' directory was added, which contains read sequence with base quanlity > Q25 and length > 140.

Followed by read type and gene type detection
```bash
python CS_design/1_extract_readname_v3.py \
    CS_design/CS3_PT

sbatch -a 0-2 \
    CS_design/2_splitfastq.sh \
    CS_design/CS3_PT


sbatch -a 0-2 \
    CS_design/3_mapping.sh \
    CS_design/CS3_PT \
    $(path of customized genome reference) 
```

We can confirm the probe type of top read sequences by running '2_2_read_summary_top50.sh', to check the details base composition adjacent polyT.
```bash
sbatch -a 0-2 CS_design/2_2_read_summary_top50.sh \
    CS_design/CS3_PT
```
We can see the differences between read end sequence captured by CS3 and conventional polyT probes, where CS3-captured-read started with "TCA".
<img width="461" height="153" alt="Screenshot 2026-06-03 at 15 40 04" src="https://github.com/user-attachments/assets/c0ee4875-c9e0-4205-b117-5bcfeaba4d93" />
<img width="457" height="156" alt="Screenshot 2026-06-03 at 15 39 51" src="https://github.com/user-attachments/assets/c27f9910-3822-4469-8b1a-df135188439c" />


### Preprocessing and Lineage Barcoding

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

If you use **Steller** in your research or find its workflows helpful, please cite our pre-print/manuscript.

---

## Contact & Support

For questions, bug reports, or feature requests, please open an Issue on the [GitHub Issues page](https://www.google.com/search?q=https://github.com/AmoriLi/Steller/issues) or reach out directly to the maintainer at **[li_zhuxia@gibh.ac.cn]**.

```

***

### Recommendations for Customizing This File:
1. **`envs/steller_env.yml`**: Ensure you place a valid environment file inside an `envs` directory so users can recreate your package setup exactly as it runs on your server.
2. **Morandi Palette**: Under Key Features and Visualization, I explicitly highlighted your sophisticated color palettes—this is a valuable selling point for bioinformaticians who struggle to visualize dozens of clones cleanly in one plot.
3. **Citation**: Update the BibTeX stub with your target journal or BioRxiv details once available.

```
