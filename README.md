# Spatial Transcriptomics Analysis

A collection of end-to-end spatial transcriptomics analysis notebooks covering multiple platforms, tools, and biological systems — from basic Scanpy workflows to advanced Squidpy spatial statistics and subcellular-resolution Xenium analysis.

---

## Repository Structure

```
Spatial-Transcriptomics/
│
├── Basic Analysis/
│   ├── basic_analysis.ipynb          # Scanpy-based workflow (Visium + MERFISH)
│   ├── Readme.md
│   └── Images/
│
├── Squidpy Visium hne/
│   ├── tutorial_visium_hne.ipynb     # Squidpy H&E mouse brain analysis
│   ├── Readme.md
│   └── Images/
│
├── Squidpy Visium Fluorescence/
│   ├── tutorial_visium_fluo.ipynb    # Squidpy fluorescence image feature analysis
│   ├── Readme.md
│   └── Images/
│
└── Xenium Analysis/
    ├── tutorial_xenium_colab_v3.ipynb  # 10x Xenium human lung cancer analysis
    ├── Readme.md
    └── Images/
```

---

## Modules

### 1. Basic Analysis — Scanpy Visium & MERFISH
**Notebook:** `Basic Analysis/basic_analysis.ipynb`  
**Platform:** 10x Genomics Visium + MERFISH  
**Tissue:** Human Lymph Node + U2-OS Cells (in culture)

A complete beginner-to-intermediate spatial transcriptomics workflow using Scanpy as the primary framework. This notebook covers two platforms in a single end-to-end tutorial:

**Visium (Human Lymph Node):**
- Loading the `V1_Human_Lymph_Node` dataset via `sc.datasets.visium_sge`
- QC with per-spot histograms for total counts, gene counts, and mitochondrial fraction
- Filtering by count range (5,000–35,000), mitochondrial content (<20%), and gene prevalence (≥10 spots)
- Library-size normalization, log1p transformation, and selection of 2,000 highly variable genes
- PCA → KNN graph → UMAP → Leiden clustering
- Spatial overlays of total counts, gene counts, cluster labels, and individual gene expression directly on the H&E tissue image
- Differential expression using t-test; top marker gene identification (e.g., CR2 for B cells)
- Multi-gene spatial visualization (COL1A2, SYPL1) with transparency control

**MERFISH (U2-OS Cells):**
- Manual spatial coordinate assignment to an AnnData object — the standard pattern for non-Visium data
- Normalized clustering and UMAP on a 15-PC reduction (smaller panel, FISH-based)
- Demonstrates that spatial structure is absent in cell culture data (clusters reflect cell cycle, not anatomy)

**Key outputs:** `qc_histograms.png`, `umap_clusters.png`, `spatial_clusters_full.png`, `spatial_clusters_zoomed.png`, `spatial_CR2.png`, `spatial_COL1A2_SYPL1.png`, `marker_heatmap.png`, `merfish_umap.png`

---

### 2. Squidpy Visium H&E — Advanced Spatial Statistics
**Notebook:** `Squidpy Visium hne/tutorial_visium_hne.ipynb`  
**Platform:** 10x Genomics Visium  
**Tissue:** Mouse Brain (coronal section, H&E stained)

A comprehensive Squidpy-based analysis that goes well beyond standard clustering. This notebook uses a pre-annotated mouse brain Visium dataset (with cluster labels from the Allen Brain Atlas and Linnarson Lab Mouse Brain Atlas) to demonstrate five spatial analysis methods:

| Analysis | Method | Output |
|---|---|---|
| Image-based clustering | Summary features → PCA → Leiden | Comparison of gene-space vs. image-space clusters |
| Spatial neighbors graph | Delaunay / radius-based adjacency | Foundation for all graph statistics |
| Neighborhood enrichment | Permutation-based co-localization z-scores | Pairwise cluster spatial interaction heatmap |
| Co-occurrence analysis | Conditional probability vs. distance | Interaction distance curves per cluster pair |
| Ligand-receptor interactions | CellPhoneDB reimplementation (Omnipath DB) | Significant interaction dot plot |
| Spatially variable genes | Moran's I autocorrelation | Top SVGs: Olfm1, Plp1, Itpka |

**Key biological findings:**
- Hippocampus, Pyramidal Layer, and Pyramidal Layer DG clusters form a tightly co-localized spatial unit (confirmed by both neighborhood enrichment and co-occurrence)
- Image morphology and gene expression overlap on large-scale anatomy (white matter vs. gray matter) but diverge on fine cortical layering — gene expression is required to resolve laminar structure
- Plp1 (oligodendrocyte myelin marker) and Olfm1 (neural marker) are the most spatially autocorrelated genes

**Key outputs:** `spatial_clusters.png`, `image_vs_gene_clusters.png`, `spatially_variable_genes.png`, `co_occurrence_hippocampus.png`

---

### 3. Squidpy Visium Fluorescence — Image Feature Extraction
**Notebook:** `Squidpy Visium Fluorescence/tutorial_visium_fluo.ipynb`  
**Platform:** 10x Genomics Visium (fluorescence imaging)  
**Tissue:** Mouse Brain (coronal section, multiplex fluorescence: DAPI, NeuN, GFAP)

This notebook focuses on the image analysis capabilities of Squidpy, using a three-channel fluorescence image (DAPI, anti-NeuN, anti-GFAP) alongside gene expression data. The central question is whether image-derived features alone can recapitulate or even improve on transcriptomic cluster annotations.

**Pipeline:**
1. Gaussian smoothing of the fluorescence image for noise reduction
2. Watershed nucleus segmentation on the DAPI channel; results stored as a label image
3. Segmentation-based features per spot: nucleus count (cell density proxy), mean NEUN intensity, mean GFAP intensity
4. Multi-scale image feature extraction:
   - Full-resolution spot circle: summary + texture + histogram features
   - Full-resolution square crop: summary + histogram features (includes local context)
   - 25%-downsampled: summary + histogram features (captures broader tissue patterns)
5. Concatenation of all feature sets into a single per-spot feature matrix
6. Leiden clustering independently from summary, histogram, and texture features

**Key biological findings:**
- The hippocampal pyramidal layer is detectable as a high-cell-density region from segmentation alone — a resolution not available from gene-expression clustering
- NEUN intensity is highest in cortical gene-expression clusters; GFAP intensity is highest in the Fiber tract and lateral ventricle clusters — both consistent with known biology
- All three image feature types subdivide the hippocampus and cortex into finer anatomical units than gene-expression Leiden clustering
- Image analysis provides spatial resolution approaching the pixel level, complementing spot-level transcriptomic resolution

**Key outputs:** `01_spatial_cluster_annotation.png`, `02_fluorescence_channels.png`, `03_nucleus_segmentation_crop.png`, `04_segmentation_features_spatial.png`, `05_feature_cluster_comparison.png`

---

### 4. Xenium Analysis — Subcellular-Resolution Human Lung Cancer
**Notebook:** `Xenium Analysis/tutorial_xenium_colab_v3.ipynb`  
**Platform:** 10x Genomics Xenium (in situ, single-molecule resolution)  
**Tissue:** Human Lung Cancer (2-FOV dataset, 161,000 cells, 480 genes)  
**Runtime:** Google Colab, fully automated (no manual data handling required)

The most technically advanced notebook in this repository. Xenium resolves individual RNA molecules at their exact two-dimensional position within fixed tissue, enabling analysis at true single-cell and subcellular resolution. The full pipeline covers data ingestion through spatial statistics:

**Data ingestion:**
- Selective download of seven required Xenium output files (265 MB), skipping large morphology TIFFs
- Parsing with `spatialdata-io.xenium` reader into a SpatialData object containing cell boundaries, nucleus boundaries, a 40M-point transcript cloud, and the 161,000 × 480 count matrix
- Zarr-format on-disk storage for efficient repeated access

**QC:**
- Negative control probe rate: ~0.005% (far below 1% threshold, confirming high signal-to-noise)
- Distribution inspection of total transcripts per cell, unique genes per cell, cell area, and nucleus-to-cell area ratio
- Filtering: ≥10 total transcripts per cell; ≥5 cells per gene

**Single-cell analysis:**
- Normalization → log1p → PCA → KNN graph → UMAP → Leiden clustering (~20 populations)
- UMAP confirms transcriptionally distinct tumor epithelial, stromal, immune, endothelial, and other lung-resident populations

**Spatial statistics:**

| Analysis | Key Finding |
|---|---|
| Spatial scatter | Tumor nests, stroma, and immune infiltrates are spatially coherent |
| Centrality scores | Stromal/vascular clusters are most spatially central; tumor clusters have highest clustering coefficient |
| Co-occurrence | Interaction distance peaks at 20–50 µm (direct contact) and 100–200 µm (compartment-level) |
| Neighborhood enrichment | Full pairwise interaction map; reveals tumor-immune boundary contacts |
| Moran's I | Top SVGs: AREG (0.696), MET (0.683), EPCAM (0.633), IGKC, IGHG1, IDO1, SPARC, APOE |

**Key biological findings:**
- AREG and MET (lung adenocarcinoma-associated) are the two most spatially autocorrelated genes, concentrated in tumor epithelial nests
- Immunoglobulin genes (IGKC, IGHG1) mark focal immune aggregates within the tumor microenvironment
- IDO1 (immunosuppression) and SPARC (stroma) define distinct spatial zones of the tumor microenvironment

**Key outputs:** `qc_histograms.png`, `umap_clusters.png`, `spatial_leiden.png`, `centrality_scores.png`, `co_occurrence.png`, `nhood_enrichment.png`, `spatial_areg_met.png`

---

## Technologies and Tools

| Tool | Role |
|---|---|
| [Scanpy](https://scanpy.readthedocs.io/) | Core single-cell preprocessing, PCA, UMAP, Leiden clustering, differential expression |
| [Squidpy](https://squidpy.readthedocs.io/) | Spatial transcriptomics analysis: neighborhood enrichment, co-occurrence, Moran's I, image features, LR interactions |
| [SpatialData](https://spatialdata.scverse.org/) | Universal multi-modal spatial omics data container (Zarr-backed) |
| [spatialdata-io](https://spatialdata-io.scverse.org/) | Platform-specific readers (Xenium, Visium, etc.) |
| [AnnData](https://anndata.readthedocs.io/) | Core data structure for observations × features matrices |
| Leiden / python-igraph | Graph-based community detection for clustering |
| CellPhoneDB (via Squidpy) | Permutation-based ligand-receptor interaction testing |
| [Omnipath](https://omnipathdb.org/) | Curated ligand-receptor interaction database |
| Moran's I | Spatial autocorrelation statistic for spatially variable gene detection |
| Watershed | Classical image segmentation for nucleus detection |
| Pandas / NumPy | Data manipulation and numerical computation |
| Matplotlib / Seaborn | Visualization |

---

## Getting Started

### Basic Analysis & Squidpy Notebooks (local)

```bash
# Install core dependencies
pip install scanpy squidpy python-igraph leidenalg openpyxl seaborn

# Launch the notebook of your choice
jupyter notebook "Basic Analysis/basic_analysis.ipynb"
jupyter notebook "Squidpy Visium hne/tutorial_visium_hne.ipynb"
jupyter notebook "Squidpy Visium Fluorescence/tutorial_visium_fluo.ipynb"
```

For reproducible Squidpy environments, use the official conda environment file from the [squidpy_notebooks repository](https://github.com/scverse/squidpy_notebooks/blob/main/environment.yml).

### Xenium Analysis (Google Colab — recommended)

1. Open [colab.research.google.com](https://colab.research.google.com)
2. Upload `Xenium Analysis/tutorial_xenium_colab_v3.ipynb`
3. Select **Runtime → Run all**
4. The notebook will install all packages, download the dataset automatically, and run the full pipeline
5. Estimated runtime: ~30–45 minutes on a standard Colab CPU instance

---

## Datasets

| Module | Dataset | Source |
|---|---|---|
| Basic Analysis | 10x Visium Human Lymph Node (`V1_Human_Lymph_Node`) | [10x Genomics](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node) |
| Basic Analysis | MERFISH U2-OS cells (Xia et al., 2019) | [PNAS 2019](https://www.pnas.org/content/116/39/19490.abstract) |
| Squidpy H&E | 10x Visium Mouse Brain (coronal, H&E) | Via `sq.datasets` / [10x Genomics](https://support.10xgenomics.com/spatial-gene-expression/datasets) |
| Squidpy Fluorescence | 10x Visium Mouse Brain (coronal, fluorescence) | Via `sq.datasets` |
| Xenium Analysis | Xenium Human Lung Cancer 2-FOV (XOA v2.0.0) | [10x Genomics](https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Lung_2fov/Xenium_V1_human_Lung_2fov_outs.zip) |

All datasets are publicly available and downloaded automatically by the notebooks on first run. No manual data preparation is required.

---

## References

- Wolf F A, Angerer P, Theis F J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology.* 2018.
- Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nature Methods.* 2022.
- Marconato L, Palla G, Yamauchi K A, et al. SpatialData: an open and universal data framework for spatial omics. *Nature Methods.* 2024.
- Efremova M, et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. *Nature Protocols.* 2020.
- Turei D, et al. Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. *Molecular Systems Biology.* 2021.
- Xia C, et al. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. *PNAS.* 2019.
- Luecken M D, Theis F J. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology.* 2019.
- 10x Genomics Visium dataset portal: https://support.10xgenomics.com/spatial-gene-expression/datasets
- 10x Genomics Xenium platform: https://www.10xgenomics.com/platforms/xenium
- Allen Mouse Brain Atlas: https://mouse.brain-map.org/
- Linnarson Lab Mouse Brain Atlas: http://mousebrain.org/
