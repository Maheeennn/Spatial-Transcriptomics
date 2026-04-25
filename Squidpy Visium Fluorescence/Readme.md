# Squidpy Visium Fluorescence Analysis Tutorial

This notebook demonstrates how to perform image-based spatial analysis of 10x Visium fluorescence data using the Squidpy library. The analysis focuses on extracting quantitative image features from a high-resolution fluorescence tissue image, clustering spots based on those features, and comparing the resulting clusters with gene-expression-derived cluster annotations. The goal is to show that image features can provide complementary or even higher-resolution biological information than gene expression alone.

---

## Table of Contents

- [Background](#background)
- [Dataset](#dataset)
- [Environment Setup](#environment-setup)
- [Analysis Pipeline](#analysis-pipeline)
  - [1. Package Imports and Data Loading](#1-package-imports-and-data-loading)
  - [2. Spatial Cluster Visualization](#2-spatial-cluster-visualization)
  - [3. Fluorescence Channel Visualization](#3-fluorescence-channel-visualization)
  - [4. Image Preprocessing and Nucleus Segmentation](#4-image-preprocessing-and-nucleus-segmentation)
  - [5. Segmentation Feature Extraction](#5-segmentation-feature-extraction)
  - [6. Multi-Scale Image Feature Extraction](#6-multi-scale-image-feature-extraction)
  - [7. Feature-Based Leiden Clustering](#7-feature-based-leiden-clustering)
- [Results and Biological Interpretation](#results-and-biological-interpretation)
- [Key Data Outputs](#key-data-outputs)
- [Dependencies](#dependencies)
- [References](#references)

---

## Background

10x Visium is a spatial transcriptomics platform that captures gene expression measurements at spatially defined spots across a tissue section, paired with a high-resolution tissue image. Each Visium spot covers a roughly 55-micron diameter area on the tissue, and the accompanying image can be used as an additional source of information beyond gene counts.

Squidpy is a Python library built on top of Scanpy and AnnData that provides tools for spatial omics analysis, including graph-based neighborhood analysis, spatial statistics, and image feature extraction. This tutorial specifically focuses on Squidpy's image analysis module (`squidpy.im`), which enables users to compute quantitative features from the tissue image and integrate them with the transcriptomic data stored in an `AnnData` object.

The central hypothesis of this tutorial is that image-derived features — such as cell density, fluorescence intensity per channel, texture, and pixel intensity distributions — can recapitulate and even refine the biological structure captured by gene-expression clustering.

---

## Dataset

The dataset used in this tutorial is a coronal section of the mouse brain, originally published by 10x Genomics and available through their [dataset portal](https://support.10xgenomics.com/spatial-gene-expression/datasets). For this tutorial, Squidpy provides a pre-processed, pre-cropped version of the dataset that can be loaded directly via `sq.datasets`.

The pre-processing pipeline mirrors the one described in the official [Scanpy Visium tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html), and cluster annotations were derived using multiple resources:

- The [Allen Brain Atlas](https://mouse.brain-map.org/)
- The [Linnarson lab Mouse Brain gene expression atlas](http://mousebrain.org/)
- Supporting published literature

The dataset is provided as two objects:

- `adata`: An `AnnData` object containing normalized gene expression counts, spatial coordinates for each Visium spot, dimensionality reductions, and pre-annotated cluster labels in `adata.obs["cluster"]`.
- `img`: A `squidpy.im.ImageContainer` object wrapping the high-resolution fluorescence tissue image.

The tissue image contains three fluorescence channels:

| Channel Index | Marker | Biological Target |
|---------------|--------|------------------|
| 0 | DAPI | DNA — marks all cell nuclei |
| 1 | anti-NEUN | NeuN protein — marks mature neurons |
| 2 | anti-GFAP | GFAP protein — marks glial cells (astrocytes) |

The crop provided covers a sub-region of the full brain section, including parts of the hippocampus, cortex, fiber tracts, and lateral ventricles.

---

## Environment Setup

To reproduce this analysis, create a conda environment using the `environment.yml` file from the [squidpy_notebooks repository](https://github.com/scverse/squidpy_notebooks/blob/main/environment.yml):

---

## Analysis Pipeline

### 1. Package Imports and Data Loading

This cell imports all required libraries and loads the dataset. `sc.logging.print_header()` prints a summary of the environment including the versions of Scanpy, AnnData, and their core dependencies, which is useful for reproducibility. The dataset is downloaded automatically from Squidpy's hosted datasets on first run and cached locally.

After this cell, `adata` contains the full transcriptomic and spatial metadata for the Visium spots, and `img` contains the fluorescence tissue image, both ready for downstream analysis.

---

### 2. Spatial Cluster Visualization

This produces a spatial scatter plot where each Visium spot is positioned at its physical location on the tissue and colored according to its pre-annotated gene-expression cluster label. The underlying tissue image is rendered in the background.

The plot establishes the baseline biological organization of the tissue section as determined by transcriptomics. Visible anatomical regions in this crop include cortical layers (Cortex_1 through Cortex_5), hippocampal regions, fiber tracts, and lateral ventricles. This serves as the reference annotation against which all image-derived cluster results will later be compared.

Because the dataset is a crop rather than the full section, the plot shows a smaller region than a typical Visium experiment, but the anatomical structure is still clearly identifiable.

---

### 3. Fluorescence Channel Visualization

This renders the tissue image one channel at a time, producing three grayscale images displayed side by side. Each image represents the raw fluorescence signal from one marker:

- **Channel 0 (DAPI):** Shows nuclei across the entire tissue. Pixel brightness corresponds to DNA content, so regions with more cells appear brighter. This channel is also used as the basis for nucleus segmentation in the following step.
- **Channel 1 (anti-NEUN):** Shows neuron distribution. Bright regions correspond to areas enriched in mature neurons, such as the cortex and hippocampus.
- **Channel 2 (anti-GFAP):** Shows glial cell distribution. Elevated signal appears in fiber tracts and around the lateral ventricles, consistent with the known biology of astrocytes and other glial populations in those regions.

This visualization is important because it grounds the subsequent quantitative analysis in observable biology. The patterns visible here will reappear in the segmentation and feature extraction results.

---

### 4. Image Preprocessing and Nucleus Segmentation

This section performs two sequential image processing steps:

**Smoothing (`sq.im.process`):** A Gaussian smoothing filter is applied to the raw image to reduce high-frequency noise. The smoothed result is stored as a new layer `img["image_smooth"]` inside the `ImageContainer`, leaving the original image layer intact. Smoothing is a standard preprocessing step before segmentation because it prevents the algorithm from treating noise spikes as separate objects.

**Watershed segmentation (`sq.im.segment`):** The Watershed algorithm is applied to the DAPI channel of the smoothed image (`channel=0`). Watershed is a classical image segmentation method that treats pixel intensity as a topographic surface and floods it from local intensity minima, creating watershed lines between distinct objects. In this context, it identifies and separates individual cell nuclei. The `chunks=1000` parameter controls how the image is tiled during processing to manage memory usage.

The result is stored as `img["segmented_watershed"]`, a label image of the same spatial dimensions as the original tissue image, where background pixels are labeled 0 and each detected nucleus is assigned a unique positive integer.

**Visualization:** A 500x500 pixel crop of the tissue (starting at pixel coordinates 2000, 2000) is shown as a two-panel figure. The left panel shows the raw DAPI channel for that crop — bright elliptical blobs are individual nuclei. The right panel shows the corresponding segmentation label image — each uniquely colored region is one detected nucleus.

This visualization confirms that the Watershed algorithm successfully delineated individual nuclei. The quality of segmentation here directly determines the accuracy of the cell-count and intensity features computed in the next step.

---

### 5. Segmentation Feature Extraction

`sq.im.calculate_image_features` iterates over each Visium spot, crops the image to that spot's region, and computes the specified features from the crop. When `features="segmentation"` is used, it derives the following per-spot measurements from the segmented label image:

- **`segmentation_label`:** The number of detected nuclei (segmented objects) within the spot boundary — effectively an estimate of cell count per spot.
- **`segmentation_ch-X_mean_intensity_mean`:** For each fluorescence channel X, the average pixel intensity within segmented nuclei masks. This provides a channel-specific signal gated to cells, not background.

All features are stored in `adata.obsm["features_segmentation"]` as a DataFrame with one row per spot. The `sq.pl.extract` helper function temporarily copies these features into `adata.obs` so they can be used as color variables in the spatial scatter plot.

**The four-panel output shows:**

- **Top-left — `segmentation_label` (cell count):** Spots are colored by how many nuclei were detected within them. The hippocampal pyramidal layer is visibly denser than surrounding regions. This fine-grained spatial heterogeneity in cell density is not captured by the gene-space clusters, which group the entire hippocampus into a single cluster.

- **Top-right — `cluster` (gene-space reference):** The same spatial plot colored by transcriptomic cluster, shown here as a reference for direct comparison. The hippocampus appears as a single color, illustrating the resolution difference.

- **Bottom-left — `segmentation_ch-1_mean_intensity_mean` (NEUN intensity):** Spots labeled as Cortex_1 and Cortex_3 in the gene cluster annotation show elevated average NEUN intensity, indicating higher neuronal content in those cortical regions.

- **Bottom-right — `segmentation_ch-2_mean_intensity_mean` (GFAP intensity):** Spots in the Fiber_tracts and lateral_ventricles gene clusters show elevated average GFAP intensity, consistent with the known enrichment of astrocytes and myelin-associated glial cells in those regions.

This step demonstrates that image-derived features can independently recapitulate known cell-type spatial organization and provide quantitative support for gene-expression cluster annotations.

---

### 6. Multi-Scale Image Feature Extraction

This section computes a richer and more diverse set of image features using three different parameter configurations, then combines them into a single feature matrix per spot.

**Feature configurations:**

- **`features_orig`:** Extracts summary, texture, and histogram features at full resolution, restricted to the circular area covered by the Visium spot (`mask_circle=True`). This is the most faithful representation of what is directly measured under each spot.

- **`features_context`:** Extracts summary and histogram features at full resolution from the full square crop around each spot (no circle mask). This includes some pixels from tissue outside the strict spot boundary, providing a small amount of spatial context from the immediate neighborhood.

- **`features_lowres`:** Extracts summary and histogram features after downsampling the image to 25% of its original resolution (`scale=0.25`). Operating at lower resolution means the crop effectively covers a larger physical area of tissue, capturing broader spatial patterns.

**Feature types explained:**

- **Summary features:** Descriptive statistics of pixel intensities within the crop, computed per channel — mean, standard deviation, and multiple percentiles. These summarize the overall brightness and spread of the fluorescence signal.

- **Histogram features:** The distribution of pixel intensities is binned into a fixed number of bins per channel. The resulting bin counts form a feature vector that describes the shape of the intensity distribution, capturing whether the signal is uniform, bimodal, skewed, and so on.

- **Texture features:** Haralick texture features computed from gray-level co-occurrence matrices (GLCM). These capture the spatial relationships between pixel intensities — properties like contrast, correlation, and energy that describe tissue microstructure patterns beyond simple intensity levels.

All three feature sets are concatenated along the column axis into `adata.obsm["features"]`, yielding a combined feature matrix that describes each spot across multiple spatial scales and feature types. Duplicate column names that arise from computing summary and histogram features across multiple configurations are made unique using `ad.utils.make_index_unique`.

---

### 7. Feature-Based Leiden Clustering

The `cluster_features` helper function takes the combined feature matrix and applies a standard single-cell clustering pipeline to it:

1. **Feature filtering:** Only features whose column names match the `like` argument are retained (e.g., only summary-type features).
2. **Scaling:** Feature values are z-score normalized across spots so that features with larger absolute ranges do not dominate the PCA.
3. **PCA:** Principal component analysis reduces the feature matrix to at most 10 components, capturing the main axes of variation across spots.
4. **Neighbor graph construction:** A k-nearest-neighbor graph is built in PCA space.
5. **Leiden clustering:** Community detection is applied to the neighbor graph to produce cluster assignments.

This function is called three times — once each for summary, histogram, and texture features — producing three independent cluster annotations stored in `adata.obs`.

**The final four-panel comparison plot shows:**

- **Summary feature clusters:** Spots grouped based on average pixel intensity statistics. These clusters tend to align with broad tissue compartments and are sensitive to overall fluorescence levels per region.

- **Histogram feature clusters:** Spots grouped based on the shape of their pixel intensity distributions. These clusters are sensitive to the mix of bright and dark regions within a spot, which can reflect differences in cell density and staining heterogeneity.

- **Texture feature clusters:** Spots grouped based on tissue microstructure patterns. These clusters often reflect differences in tissue organization — for example, the dense, regularly arranged nuclei of the hippocampal pyramidal layer versus the more dispersed arrangement in cortical layers.

- **Gene-space clusters (reference):** The original transcriptomic annotation shown for comparison.

All three image-derived cluster sets are spatially coherent — nearby spots tend to be assigned to the same cluster. More importantly, all three subdivide the hippocampus into multiple clusters corresponding to its distinct structural layers, while the gene-expression annotation assigns the whole hippocampus to a single cluster. Similarly, the cortex is subdivided with finer granularity by image features than by gene expression.

---

## Results and Biological Interpretation

The analysis produces four categories of results:

**Cell density mapping:** The segmentation feature `segmentation_label` reveals that the hippocampal pyramidal layer is significantly more cell-dense than surrounding tissue, a structural feature that is biologically well-established but not captured in gene-expression clusters. The image-derived cell count per spot provides a quantitative proxy for local cell density without requiring single-cell dissociation.

**Cell-type-specific fluorescence:** Average channel intensities within segmented nuclei serve as spatial maps of cell type enrichment. Elevated NEUN signal in cortical clusters confirms their neuronal character, while elevated GFAP signal in fiber tracts and ventricle-adjacent regions confirms the presence of glia, consistent with white matter and ependymal cell biology.

**Image features resolve finer anatomical structure:** Across all three image feature types (summary, histogram, texture), the resulting Leiden clusters subdivide the hippocampus and cortex into more regions than gene-expression clustering. This demonstrates that morphological and photometric information in the fluorescence image carries biologically meaningful signal that is distinct from or complementary to transcriptomic signal.

**Multi-scale features are informative:** The combination of features at different spatial scales (spot-level, neighborhood-level, low-resolution) and different feature types produces cluster assignments that each reflect different aspects of tissue organization, none of which are strictly redundant with gene-space clusters.

The overall conclusion is that image feature analysis is a valuable complement to spatial transcriptomics. It is computationally inexpensive relative to additional sequencing assays, it provides spatial resolution approaching the image pixel level rather than the spot level, and it captures continuous morphological variation that gene expression may only discretize coarsely into clusters.

---

## Key Data Outputs

| Object | Location | Description |
|--------|----------|-------------|
| Smoothed image | `img["image_smooth"]` | Gaussian-smoothed version of the raw fluorescence image |
| Segmentation label image | `img["segmented_watershed"]` | Per-pixel nucleus labels from Watershed segmentation |
| Segmentation features | `adata.obsm["features_segmentation"]` | Cell count and per-channel mean intensities per spot |
| Full-res masked features | `adata.obsm["features_orig"]` | Summary, texture, and histogram features under spot circle at full resolution |
| Full-res context features | `adata.obsm["features_context"]` | Summary and histogram features from square crop at full resolution |
| Low-res features | `adata.obsm["features_lowres"]` | Summary and histogram features at 0.25x resolution |
| Combined feature matrix | `adata.obsm["features"]` | Concatenation of all three feature sets above |
| Summary clusters | `adata.obs["features_summary_cluster"]` | Leiden clusters derived from summary image features |
| Histogram clusters | `adata.obs["features_histogram_cluster"]` | Leiden clusters derived from histogram image features |
| Texture clusters | `adata.obs["features_texture_cluster"]` | Leiden clusters derived from texture image features |

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `scanpy` | Single-cell analysis framework; PCA, neighbor graph, Leiden clustering, plotting |
| `squidpy` | Spatial omics analysis; image feature extraction, segmentation, spatial plotting |
| `anndata` | Core data structure for storing expression matrices, observations, and embeddings |
| `leidenalg` | Backend implementation of the Leiden community detection algorithm |
| `pandas` | Feature matrix manipulation and concatenation |
| `matplotlib` | Figure rendering |

---

## References

- Palla, G., Spitzer, H., Klein, M., et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- 10x Genomics Visium dataset portal: https://support.10xgenomics.com/spatial-gene-expression/datasets
- Allen Mouse Brain Atlas: https://mouse.brain-map.org/
- Mouse Brain gene expression atlas (Linnarson lab): http://mousebrain.org/
