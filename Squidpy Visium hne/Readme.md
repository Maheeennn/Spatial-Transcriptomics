# Spatial Transcriptomics Analysis of Mouse Brain Visium H&E Data

## Overview

This notebook presents a complete end-to-end spatial transcriptomics analysis of a 10x Genomics Visium slide derived from a coronal section of the mouse brain. The analysis is built on top of the **Squidpy** library, a Python framework specifically designed for the analysis of spatial molecular data. The workflow integrates both gene expression data and high-resolution tissue histology images to extract biologically meaningful insights about tissue organization, cell-cell communication, and spatially variable gene expression patterns.

The dataset is sourced from the [10x Genomics dataset portal](https://support.10xgenomics.com/spatial-gene-expression/datasets) and is provided in pre-processed form via Squidpy's built-in dataset loader. Cluster annotations were established using multiple reference resources including the Allen Brain Atlas, the Linnarson Lab Mouse Brain gene expression atlas, and recent published literature.

---

## Table of Contents

- [Background and Motivation](#background-and-motivation)
- [Dependencies and Environment Setup](#dependencies-and-environment-setup)
- [Dataset Description](#dataset-description)
- [Notebook Structure and Cell-by-Cell Analysis](#notebook-structure-and-cell-by-cell-analysis)
  - [1. Package Import and Data Loading](#1-package-import-and-data-loading)
  - [2. Spatial Visualization of Gene Expression Clusters](#2-spatial-visualization-of-gene-expression-clusters)
  - [3. Image Feature Extraction](#3-image-feature-extraction)
  - [4. Image-Based Clustering](#4-image-based-clustering)
  - [5. Spatial Neighbors Graph Construction](#5-spatial-neighbors-graph-construction)
  - [6. Neighborhood Enrichment Analysis](#6-neighborhood-enrichment-analysis)
  - [7. Co-occurrence Analysis](#7-co-occurrence-analysis)
  - [8. Ligand-Receptor Interaction Analysis](#8-ligand-receptor-interaction-analysis)
  - [9. Spatially Variable Genes with Moran's I](#9-spatially-variable-genes-with-morans-i)
- [Key Results and Biological Interpretation](#key-results-and-biological-interpretation)
- [Technologies Used](#technologies-used)
- [References](#references)

---

## Background and Motivation

Spatial transcriptomics is a cutting-edge genomics technology that allows researchers to measure gene expression while preserving the spatial organization of cells within a tissue. The 10x Genomics Visium platform captures mRNA from tissue sections placed on a spatially barcoded array, enabling the mapping of thousands of genes simultaneously across the tissue.

This notebook goes beyond simple gene expression visualization. It combines:

- **Histological image analysis** using convolutional-style image feature extraction
- **Spatial graph statistics** including neighborhood enrichment and co-occurrence scoring
- **Cell-cell communication inference** using a reimplementation of CellPhoneDB
- **Spatial autocorrelation analysis** using Moran's I to detect spatially variable genes

The mouse brain is an ideal tissue for this type of analysis because of its well-defined anatomical regions, each with distinct cell type compositions and gene expression signatures.

---

## Dependencies and Environment Setup

The notebook is designed to be run in a conda environment. The environment can be set up using:

The environment file is available in the [squidpy_notebooks GitHub repository](https://github.com/scverse/squidpy_notebooks/blob/main/environment.yml).

The following core packages are used throughout the notebook:

- `numpy` and `pandas` for numerical and tabular data manipulation
- `anndata` for the AnnData data structure, which is the standard for single-cell and spatial genomics data
- `scanpy` for single-cell analysis utilities including PCA, neighbor computation, and clustering
- `squidpy` for spatial transcriptomics-specific analysis functions and plotting
- `leidenalg` and `python-igraph` for Leiden community detection used in clustering
- `matplotlib` for plot rendering (invoked inline via `%matplotlib inline`)

---

## Dataset Description

The dataset consists of two main objects:

**1. AnnData object (`adata`)**

This is the primary data container, storing:
- A matrix of gene expression counts across all Visium capture spots
- Spatial coordinates of each spot on the tissue slide stored in `adata.obsm['spatial']`
- Pre-computed cluster annotations stored in `adata.obs['cluster']`
- Pre-computed highly variable gene flags in `adata.var`
- Various embedding results added during the analysis in `adata.obsm`

The cluster annotations represent distinct anatomical brain regions including the Hippocampus, Pyramidal layer, Pyramidal layer dentate gyrus, Fiber tract, and cortical layers, among others.

**2. ImageContainer object (`img`)**

This is a Squidpy-specific data structure that wraps the high-resolution H&E stained tissue image. The H&E (Hematoxylin and Eosin) staining highlights cell nuclei in blue-purple and cytoplasm and extracellular matrix in pink. This image provides morphological context that is analyzed in parallel with the gene expression data.

---

## Notebook Structure and Cell-by-Cell Analysis

### 1. Package Import and Data Loading

The first executable cell sets up the plotting backend with `%matplotlib inline`, then installs `squidpy` via pip. Subsequently, all required packages are imported and the pre-processed dataset is downloaded and loaded using:

**Output and what it means:**

This cell prints the Scanpy header, which displays version information for all key single-cell analysis packages (`scanpy`, `anndata`, `umap`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `statsmodels`). This output serves as an environment verification step and is critical for reproducibility. The squidpy version is also printed. This diagnostic output ensures that the analysis was performed with a known, consistent software stack.

---

### 2. Spatial Visualization of Gene Expression Clusters

This single line calls Squidpy's spatial scatter plot function, which overlays the cluster annotations directly onto the H&E tissue image. Each spot is colored by its assigned cluster label, and the underlying histology image provides anatomical context.

**Output and what it means:**

The output is a scatter plot rendered on top of the H&E stained tissue image of the mouse brain coronal section. Each dot represents one Visium capture spot (approximately 55 micrometers in diameter), and each color corresponds to a distinct annotated brain region. The result visually confirms that the gene expression-based cluster annotations align with known anatomical structures of the mouse brain. For example, the hippocampus, cortex, and white matter tracts are clearly delineated by different clusters. This is the foundational visualization of the entire analysis, establishing that the pre-processing and annotation were performed correctly.

---

### 3. Image Feature Extraction

Summary features include channel-wise statistics (mean, standard deviation, percentiles) of the pixel intensities within each spot's crop from the tissue image. By computing these at two scales (`scale=1.0` and `scale=2.0`), the analysis captures both local fine-grained morphology and broader tissue context surrounding each spot.

The resulting feature matrices from both scales are then concatenated into a single `obs x features` matrix and stored in `adata.obsm["features"]`. Duplicate column names are resolved using `ad.utils.make_index_unique`.

**Output and what it means:**

This cell does not produce a visual output directly. Instead, it populates `adata.obsm["features_summary_scale1.0"]`, `adata.obsm["features_summary_scale2.0"]`, and the combined `adata.obsm["features"]`. The successful execution of this cell (without errors) means that for every Visium spot, a numerical descriptor of the tissue morphology has been computed. This is conceptually important: it means the notebook is now working with two parallel representations of each spot — one from gene expression and one from image morphology.

---

### 4. Image-Based Clustering

**What the code does:**

A helper function `cluster_features` is defined and applied. This function takes the image feature matrix, scales it, computes PCA with up to 10 components, builds a nearest-neighbor graph, and then applies Leiden community detection:

**Output and what it means:**

The output is a side-by-side spatial scatter plot showing two cluster annotations on the tissue:
- **Left panel:** `features_cluster` — clusters derived purely from tissue image morphology
- **Right panel:** `cluster` — clusters derived from gene expression

This comparison is one of the most intellectually significant outputs in the notebook. It allows for a direct visual assessment of how much information is shared between the histological appearance of the tissue and its molecular identity.

The key biological finding visible in this output: some regions show strong concordance between image-based and gene-based clusters, particularly the Fiber tract and Hippocampus regions. This makes sense because these regions have highly distinctive morphological appearances that reflect their unique gene expression profiles. In the cortex, however, the two clustering schemes diverge. Gene expression clusters reveal the layered laminar structure of the cortex, while image features cluster the cortex more by broad cortical region rather than layer. This highlights both the power and the limitation of image-based analysis — it captures large-scale tissue organization well but misses molecular subtype differences that are not morphologically distinct.

---

### 5. Spatial Neighbors Graph Construction

**What the code does:**

```python
sq.gr.spatial_neighbors(adata)
```

This function constructs a spatial connectivity graph based on the physical positions of spots on the tissue. Spots that are spatially adjacent (neighbors on the Visium array) are connected by an edge. The resulting adjacency matrix is stored in `adata.obsp['spatial_connectivities']` and the distances are stored in `adata.obsp['spatial_distances']`.

**Output and what it means:**

This cell produces no visual output. It is a prerequisite step for all subsequent graph-based analyses (neighborhood enrichment, co-occurrence). The spatial neighbor graph formalizes the concept of "proximity" in the tissue by defining which spots are physically close to one another. Unlike transcriptomic neighbor graphs (which connect spots with similar gene expression), this graph connects spots that are physically adjacent, regardless of their molecular similarity.

---

### 6. Neighborhood Enrichment Analysis

**What the code does:**

```python
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

The `nhood_enrichment` function computes a permutation-based enrichment score for every pair of clusters. The score reflects how often spots from two different clusters appear as spatial neighbors compared to what would be expected by chance (based on random permutations of cluster labels). A high positive score means two clusters are spatially co-located more than expected; a negative score means they are more spatially segregated than expected.

**Output and what it means:**

The output is a heatmap (matrix plot) where rows and columns represent cluster labels, and the color of each cell encodes the enrichment score for that cluster pair. The color scale typically runs from blue (depletion/segregation) through white (no enrichment) to red (enrichment/co-localization).

The key biological finding from this plot is the strong neighborhood enrichment between the Hippocampus cluster and both the Pyramidal layer and Pyramidal layer dentate gyrus clusters. This is anatomically expected — in the mouse brain, the pyramidal cell layers of the hippocampus are directly embedded within the broader hippocampal region. The quantitative enrichment score confirms what is visible to the eye: these clusters are spatially interleaved in a non-random manner. Conversely, some cluster pairs show depletion, indicating that they tend to occupy spatially segregated domains of the tissue.

---

### 7. Co-occurrence Analysis

Unlike neighborhood enrichment (which operates on the discrete connectivity graph), co-occurrence analysis works directly on continuous spatial coordinates. For a given source cluster, it computes the conditional probability of observing spots from each other cluster within increasing radii. The score is the ratio of the conditional probability to the marginal probability:

`score = p(target cluster | within radius of source cluster) / p(target cluster)`

A score greater than 1 indicates the target cluster appears more frequently near the source than expected by chance.

**Output and what it means:**

The output is a line plot with distance (in pixels of the source image) on the x-axis and the co-occurrence score on the y-axis. A separate line is drawn for each cluster, conditioned on proximity to the Hippocampus cluster.

The key finding is that the Pyramidal layer cluster shows the highest co-occurrence score at short distances from the Hippocampus, and this score decreases as the radius increases. This confirms the tight spatial association between these anatomically related regions and demonstrates that the relationship is strongest at very short distances — consistent with the Pyramidal layer being physically embedded in or adjacent to the Hippocampus. This analysis is complementary to neighborhood enrichment but provides an additional dimension: the distance at which co-occurrence is strongest.

---

### 8. Ligand-Receptor Interaction Analysis

This analysis uses Squidpy's reimplementation of CellPhoneDB, combined with the Omnipath ligand-receptor interaction database. For each pair of clusters, it tests whether annotated ligand-receptor pairs are co-expressed at levels higher than expected by chance (assessed via permutation testing). The results are filtered to show only interactions with high mean expression (`means_range=(3, np.inf)`) and a very stringent significance threshold (adjusted p-value < 0.0001, `alpha=1e-4`).

**Output and what it means:**

The output is a dot plot where:
- The x-axis represents ligand-receptor interaction pairs
- The y-axis represents the cluster pair (source: Hippocampus, target: Pyramidal layer or Pyramidal layer dentate gyrus)
- The size of each dot encodes the mean expression level of the interaction pair
- The color of each dot encodes the statistical significance (smaller p-value = more intense color)

This is arguably the most hypothesis-generating output in the notebook. It identifies specific molecular signaling interactions that may be driving communication between hippocampal and pyramidal layer cells. These are candidate interactions that could be validated experimentally or investigated further with single-cell deconvolution methods to understand which specific cell types within each spot are responsible for these interactions.

---

### 9. Spatially Variable Genes with Moran's I

Moran's I is a classical spatial autocorrelation statistic. For each gene, it measures whether spots with high expression of that gene tend to cluster spatially (positive autocorrelation, Moran's I close to +1), are randomly distributed (I close to 0), or are dispersed in a checkerboard-like pattern (negative I). The analysis is run on a subset of 1000 highly variable genes for computational efficiency. Results (Moran's I statistic, p-value, adjusted p-value) are stored in `adata.uns['moranI']`.

The top 10 genes are then displayed:

**Output and what it means:**

**Table output (`adata.uns["moranI"].head(10)`):**

This is a pandas DataFrame showing the top 10 genes sorted by Moran's I statistic in descending order. Each row contains the gene name, the Moran's I score, the raw p-value, and the Benjamini-Hochberg adjusted p-value. Genes with Moran's I close to 1 and very small adjusted p-values are genes whose expression patterns are highly spatially structured — they are not randomly distributed across the tissue, but instead concentrated in specific anatomical regions. These are the most biologically interesting candidates for spatially restricted gene expression programs.

**Spatial scatter plot output:**

The final visualization shows the spatial expression patterns of three top-ranked spatially variable genes (Olfm1, Plp1, Itpka) alongside the cluster annotation:

- **Olfm1** (Olfactomedin 1): A secreted glycoprotein involved in neural development. Its spatially restricted expression pattern, visible in this plot, associates it with specific brain subregions, most likely related to pyramidal or hippocampal populations.
- **Plp1** (Proteolipid Protein 1): A major myelin constituent. Its spatial pattern is expected to align tightly with white matter tracts and fiber tract clusters, consistent with its known role in oligodendrocytes.
- **Itpka** (Inositol 1,4,5-trisphosphate 3-kinase A): A calcium signaling enzyme. Its spatial pattern in the context of the mouse brain has been associated with pyramidal cell layers, consistent with high neuronal activity in those regions.

By overlaying these expression patterns with the cluster annotation (fourth panel), it becomes immediately clear how the top spatially variable genes define and reinforce the cluster boundaries identified by the transcriptomic analysis.

---

## Key Results and Biological Interpretation

Taken together, the analyses in this notebook converge on several consistent biological findings:

**Hippocampal organization is a dominant spatial signal.** Across multiple independent analyses — neighborhood enrichment, co-occurrence, and ligand-receptor interactions — the Hippocampus, Pyramidal layer, and Pyramidal layer dentate gyrus clusters emerge as a tightly associated spatial unit. This reflects the known anatomical organization of the hippocampal formation in the mouse brain.

**Image morphology and gene expression are partially redundant.** The comparison of gene-derived and image-derived clustering shows that large-scale tissue organization (white matter vs gray matter, major anatomical boundaries) is captured by histology, but finer molecular distinctions (e.g., cortical layers) require gene expression information.

**Spatially variable genes reflect known brain cell type biology.** The top Moran's I genes include markers of myelinating oligodendrocytes (Plp1) and neurons (Olfm1, Itpka), and their spatial patterns are anatomically coherent. This validates the Moran's I approach as a principled method for gene discovery in spatial data.

**Ligand-receptor analysis points to specific signaling hypotheses.** The CellPhoneDB-style analysis identifies statistically significant interaction pairs between hippocampal and pyramidal cell populations, providing a ranked list of candidate intercellular communication mechanisms for follow-up investigation.

---

## Technologies Used

| Tool | Purpose |
|---|---|
| Squidpy | Spatial transcriptomics analysis and visualization |
| Scanpy | Single-cell preprocessing, PCA, UMAP, clustering |
| AnnData | Core data structure for storing all analysis results |
| Leiden algorithm | Graph-based community detection for clustering |
| CellPhoneDB (via Squidpy) | Ligand-receptor interaction permutation testing |
| Omnipath | Ligand-receptor interaction database |
| Moran's I | Spatial autocorrelation for spatially variable gene detection |
| Pandas / NumPy | Data manipulation and numerical computation |
| Matplotlib | Figure rendering |

---

## References

- Palla, G. et al. Squidpy: a scalable framework for spatial omics analysis. *Nature Methods* (2022).
- Efremova, M. et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. *Nature Protocols* (2020).
- Turei, D. et al. Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. *Molecular Systems Biology* (2021). [Omnipath]
- Wolf, F.A. et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* (2018).
- 10x Genomics Visium Spatial Gene Expression dataset: Mouse Brain Coronal Section. [https://support.10xgenomics.com/spatial-gene-expression/datasets](https://support.10xgenomics.com/spatial-gene-expression/datasets)
- Allen Brain Atlas: [https://mouse.brain-map.org](https://mouse.brain-map.org)
- Linnarson Lab Mouse Brain Atlas: [http://mousebrain.org](http://mousebrain.org)
