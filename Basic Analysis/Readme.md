# Spatial Transcriptomics Analysis with Scanpy

**Author:** Giovanni Palla  
**Notebook:** `basic_analysis.ipynb`  
**Framework:** [Scanpy](https://scanpy.readthedocs.io/) — Single-Cell Analysis in Python

---

## Overview

This notebook demonstrates a complete, end-to-end workflow for analyzing and visualizing spatial transcriptomics data using Scanpy. It is aimed at researchers and bioinformaticians who work with spatially-resolved gene expression data and want to integrate tissue morphology information with transcriptomic clustering. The notebook covers two major spatial transcriptomics platforms: 10x Genomics Visium (the primary focus) and MERFISH (as a secondary example), and walks through every major step from raw data loading to marker gene visualization.

The core idea behind spatial transcriptomics is that traditional single-cell RNA sequencing loses the positional information of where a cell exists within a tissue. Spatial transcriptomics platforms like Visium capture gene expression while preserving the spatial location of each measurement spot on the tissue section, enabling researchers to directly relate gene expression patterns to tissue structure, morphology, and cell organization.

---

## Dataset

### Primary Dataset: 10x Genomics Visium — Human Lymph Node

The main dataset used throughout this notebook is a publicly available Visium spatial transcriptomics dataset of a human lymph node, provided by 10x Genomics. It can be accessed at the [10x Genomics website](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node).

This dataset is loaded directly using Scanpy's built-in downloader function:

```python
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
```

The returned object is an `AnnData` structure that contains:
- Raw gene expression count matrix (spots x genes)
- High-resolution H&E (Hematoxylin and Eosin) histology image
- Spatial coordinates for each capture spot on the tissue section
- Metadata per spot including total counts and gene counts

### Secondary Dataset: MERFISH — U2-OS Cells (Xia et al., 2019)

The notebook also includes a MERFISH example using data from [Xia et al. 2019](https://www.pnas.org/content/116/39/19490.abstract), which measured gene expression in cultured U2-OS cells. This example demonstrates how to work with FISH-based spatial data by manually assigning spatial coordinates to an AnnData object. The MERFISH section is more compact and serves as a reference for users with non-Visium data formats.

---

## Environment and Dependencies

The notebook installs and imports the following key libraries:

```python
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scanpy as sc
```

Additional installations required:
- `scanpy` — core analysis library for single-cell and spatial data
- `python-igraph` — required for Leiden graph-based clustering
- `openpyxl` — required for reading the MERFISH coordinate Excel file

Scanpy's verbosity is set to level 3 (full output), and figure parameters are set to a white background at 8x8 inches for clean, publication-ready visualizations.

---

## Workflow

The notebook is organized into the following sections:

### 1. Reading the Data

The Visium dataset is loaded into an `AnnData` object. After loading, gene names are made unique to avoid conflicts, and mitochondrial genes are flagged using the `MT-` prefix convention:

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```

QC metrics are computed for each spot, including:
- `total_counts`: total number of RNA molecules detected per spot
- `n_genes_by_counts`: number of genes detected per spot
- `pct_counts_mt`: percentage of counts originating from mitochondrial genes

The AnnData object printed at this stage reveals the full dimensionality of the dataset — number of spots (observations) and genes (variables) — as well as all the embedded metadata layers including images and spatial coordinates.

---

### 2. QC and Preprocessing

#### Quality Control Histograms

A panel of four histograms is generated to inspect the distribution of key QC metrics:

- `total_counts` (all spots): the full distribution of RNA molecule counts per spot
- `total_counts` filtered to below 10,000: a zoomed view to inspect the lower end of the distribution where low-quality spots tend to cluster
- `n_genes_by_counts` (all spots): distribution of the number of distinct genes detected per spot
- `n_genes_by_counts` filtered to below 4,000: a zoomed version to inspect the lower tail

These histograms are critical for choosing appropriate filtering thresholds. Spots with very few counts or genes are likely empty or poor quality, while spots with extremely high counts can indicate doublets or technical artifacts.

#### Filtering

Based on the QC distributions, the following filters are applied:

```python
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=10)
```

- Spots with fewer than 5,000 total counts are removed (too low quality)
- Spots with more than 35,000 total counts are removed (likely artifacts)
- Spots where more than 20% of counts come from mitochondrial genes are removed (dying or damaged cells)
- Genes detected in fewer than 10 spots are removed

The notebook prints the number of spots that remain after mitochondrial filtering, giving a concrete count of how many spots passed quality control.

#### Normalization and Feature Selection

Counts are library-size normalized using `normalize_total`, followed by a log1p transformation to stabilize variance:

```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

The top 2,000 highly variable genes are identified using the Seurat method. These genes carry the most biological signal and are used as input for dimensionality reduction and clustering in subsequent steps.

---

### 3. Manifold Embedding and Clustering

#### Dimensionality Reduction

PCA is run first to reduce the dimensionality from thousands of genes to a manageable number of principal components:

```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

A k-nearest neighbor graph is then constructed in PCA space, which forms the basis for both clustering and UMAP embedding.

#### Clustering

Leiden community detection is applied to the neighbor graph to identify transcriptionally similar groups of spots:

```python
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
```

This produces cluster labels stored in `adata.obs["clusters"]`.

#### UMAP Visualization

A UMAP plot is generated showing three features side by side:
- `total_counts` per spot
- `n_genes_by_counts` per spot
- `clusters` — Leiden cluster identity

This triple UMAP is an important quality check. If total counts or gene counts align suspiciously well with cluster boundaries, it can indicate that technical variation rather than true biology is driving the clustering. In a well-processed dataset, the biological cluster structure should not strongly follow technical covariates.

---

### 4. Visualization in Spatial Coordinates

This section is the core of what makes spatial transcriptomics uniquely powerful. The `sc.pl.spatial` function overlays spot-level information directly on top of the H&E tissue image.

#### Total Counts and Gene Counts in Space

```python
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
```

This produces a two-panel figure where each panel shows the H&E image with spots colored by either total RNA content or number of detected genes. This allows the researcher to see whether RNA-rich regions correlate with morphologically distinct tissue structures.

#### Cluster Visualization in Space

```python
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
```

This is one of the most important outputs of the notebook. Each spot on the tissue is colored according to its Leiden cluster assignment. The critical observation made in the notebook is that spots belonging to the same cluster in gene expression space tend to occupy coherent, contiguous regions in tissue space — a hallmark of true biological organization. Specifically, spots in cluster 5 are noted to be frequently surrounded by spots from cluster 0, suggesting spatial co-organization of distinct cell populations.

#### Zoomed Cluster View with Transparency

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)
```

This zoomed-in visualization focuses on clusters 5 and 9 within a specific crop of the tissue image. By setting alpha to 0.5, the spot colors become semi-transparent, allowing the underlying H&E histology to be visible through the colored overlay. This enables a direct visual link between gene expression cluster identity and tissue morphology.

The `sc.pl.spatial` function supports additional parameters discussed in the notebook:
- `img_key`: which stored image to use (e.g., `"hires"` or `"lowres"`)
- `crop_coord`: pixel coordinates to crop the image (left, right, top, bottom)
- `alpha_img`: transparency of the background image itself
- `bw`: flag to convert image to grayscale
- `size`: scales spot sizes (behaves differently than in other Scanpy plotting functions — here it is a scaling factor, not an absolute size)

---

### 5. Cluster Marker Genes

#### Differential Expression

Marker genes are computed for each cluster using a t-test:

```python
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
```

This ranks genes by their specificity to each cluster, identifying genes that are significantly upregulated in one cluster relative to all others.

#### Heatmap of Top Marker Genes

```python
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")
```

A heatmap is generated showing the top 10 marker genes for cluster 9 across all clusters. This provides a compact but dense view of how the selected marker genes are expressed across the entire dataset, making it possible to assess the specificity of the identified markers.

#### Spatial Expression of CR2

A key finding highlighted in the notebook is that the gene CR2, identified as a top marker for cluster 9, also shows a spatially coherent expression pattern:

```python
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])
```

This side-by-side plot confirms that the transcriptional cluster identity is recapitulated in physical space — the region where cluster 9 spots are concentrated also shows high CR2 expression. CR2 (Complement Receptor 2) is a well-known marker of mature B cells, which is biologically expected in lymph node tissue.

#### Spatial Expression of COL1A2 and SYPL1

```python
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)
```

Two additional genes — COL1A2 (Collagen Type I Alpha 2 Chain) and SYPL1 (Synaptophysin-Like Protein 1) — are visualized spatially with 70% transparency. COL1A2 is associated with connective tissue and fibroblasts, while SYPL1 is involved in vesicle transport. Plotting them side by side allows a spatial comparison of their expression territories across the lymph node section.

---

### 6. MERFISH Example

The final section of the notebook demonstrates how to work with FISH-based spatial transcriptomics data, using the MERFISH dataset from Xia et al. 2019.

#### Data Loading

Unlike Visium, MERFISH data does not come bundled with a built-in Scanpy loader. Coordinates and counts are loaded separately and combined manually:

```python
coordinates = pd.read_excel("./data/pnas.1912459116.sd15.xlsx", index_col=0)
counts = sc.read_csv("./data/pnas.1912459116.sd12.csv").transpose()
adata_merfish = counts[coordinates.index, :].copy()
adata_merfish.obsm["spatial"] = coordinates.to_numpy()
```

This pattern — assigning a NumPy array of coordinates to `adata.obsm["spatial"]` — is the standard way to store spatial information for any custom or non-Visium dataset in the Scanpy/Squidpy ecosystem.

#### Preprocessing and Clustering

Standard preprocessing is applied with per-cell normalization to 1 million counts:

```python
sc.pp.normalize_per_cell(adata_merfish, counts_per_cell_after=1e6)
sc.pp.log1p(adata_merfish)
sc.pp.pca(adata_merfish, n_comps=15)
sc.pp.neighbors(adata_merfish)
sc.tl.umap(adata_merfish)
sc.tl.leiden(adata_merfish, key_added="clusters", resolution=0.5, n_iterations=2, flavor="igraph", directed=False)
```

Fifteen PCA components are used (fewer than the Visium dataset, reflecting the smaller gene panel typical of FISH-based methods).

#### MERFISH Cluster Visualization

Two plots are generated:
1. UMAP colored by Leiden cluster
2. Spatial embedding colored by Leiden cluster using `sc.pl.embedding` with `basis="spatial"`

An important biological note made in the notebook: since this MERFISH experiment measured U2-OS cells in culture (not in a tissue), the clusters represent cell cycle states rather than distinct cell types. No coherent spatial structure is expected or observed — cells at different stages of the cell cycle are randomly distributed across the culture dish. This makes it a useful contrast to the Visium lymph node dataset, where strong spatial organization is present.

---

## Key Results and Biological Insights

The notebook produces several important findings:

- After QC filtering, a clean set of spots from the human lymph node tissue is retained, with total counts between 5,000 and 35,000 and mitochondrial content below 20%.
- Leiden clustering identifies distinct transcriptional populations, and these populations show strong spatial organization within the lymph node tissue section.
- Cluster 5 spots are frequently surrounded by cluster 0 spots, suggesting a spatial niche or communication relationship between these two populations.
- CR2 expression recapitulates the spatial distribution of cluster 9, and because CR2 is a known B cell marker, this confirms that cluster 9 corresponds to a B cell-rich region of the lymph node — which is biologically expected in lymphoid tissue.
- COL1A2 and SYPL1 mark spatially distinct regions of the tissue, consistent with connective tissue and vesicle-trafficking-related processes respectively.
- The MERFISH example demonstrates that the same Scanpy workflow can be applied to non-Visium platforms with minor adaptations, specifically the manual assignment of spatial coordinates.

---

## File Structure

```
.
├── basic_analysis.ipynb        # Main notebook
└── data/
    ├── pnas.1912459116.sd15.xlsx   # MERFISH spatial coordinates
    └── pnas.1912459116.sd12.csv    # MERFISH gene expression counts
```

The Visium dataset is downloaded automatically by Scanpy and cached locally when `sc.datasets.visium_sge` is first called.

---

## How to Run

1. Install the required packages:
   ```bash
   pip install scanpy python-igraph openpyxl seaborn
   ```

2. Download the MERFISH data files from the Xia et al. 2019 supplementary materials and place them in the `./data/` directory with the filenames shown above.

3. Open and run `basic_analysis.ipynb` sequentially from top to bottom. The Visium dataset will be downloaded automatically on first run.

4. For up-to-date spatial analysis workflows and extended functionality including spatial statistics, neighborhood enrichment analysis, and ligand-receptor interaction analysis, refer to the [Squidpy tutorials](https://squidpy.readthedocs.io/en/stable/tutorials.html).

---

## References

- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Mol Syst Biol.* 2019. [Link](https://www.embopress.org/doi/full/10.15252/msb.20188746)
- Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression (SCTransform). *Genome Biol.* 2019. [Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
- Townes FW et al. Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model (GLM-PCA). *Genome Biol.* 2019. [Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1861-6)
- Xia C et al. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. *PNAS.* 2019. [Link](https://www.pnas.org/content/116/39/19490.abstract)

---

## Notes

- This notebook was originally authored as part of the Scanpy documentation tutorials.
- For production use with your own Visium data, use `squidpy.read.visium` instead of `sc.datasets.visium_sge`.
- The notebook is a foundational tutorial. More advanced analyses (spatial autocorrelation, neighborhood enrichment, co-occurrence analysis) are covered in the Squidpy documentation.
