# Spatial Transcriptomics Analysis of Human Lung Cancer Using 10x Genomics Xenium

## Overview

This repository contains a Jupyter notebook that performs a complete spatial transcriptomics analysis pipeline on a real human lung cancer tissue sample profiled with the 10x Genomics Xenium in situ platform. The notebook is designed to run in Google Colab without any manual file handling. It automatically downloads the dataset, parses it using the official Xenium reader, and walks through the full analysis from raw data ingestion to spatial statistics.

Xenium is a subcellular-resolution in situ gene expression technology. Unlike bulk or standard single-cell RNA sequencing, Xenium physically images individual RNA molecules directly inside fixed tissue sections. This preserves the exact two-dimensional position of each detected transcript within the tissue, making it possible to ask not just what genes are expressed in a cell, but where those cells are located within the tissue and how they relate spatially to their neighbors. The result is a dataset that combines the gene expression depth of single-cell sequencing with the spatial context of histology.

The analysis in this notebook follows the official squidpy Xenium tutorial and uses the scverse ecosystem throughout, specifically the `spatialdata`, `spatialdata-io`, `spatialdata-plot`, `scanpy`, and `squidpy` libraries.

---

## Dataset

The dataset is the 10x Genomics Xenium Human Lung Cancer 2-FOV dataset (XOA version 2.0.0). It is publicly available at:

```
https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Lung_2fov/Xenium_V1_human_Lung_2fov_outs.zip
```

This is a compact 265 MB version of the full dataset that contains two fields of view from a formalin-fixed paraffin-embedded human lung cancer section. It has the same file structure and output format as the full 7.5 GB version. The full version includes large morphology TIFF images which are not required for the transcript-level analysis performed here. The notebook automatically downloads and extracts only the seven files needed, skipping the TIFF images entirely.

The seven files used are:

- `experiment.xenium` contains the run metadata including the software version and the pixel-to-micron conversion factor, which is required by the reader to register coordinates.
- `cells.parquet` contains per-cell summary statistics including centroid coordinates in micrometers, cell area, nucleus area, nucleus count, and counts of negative control probes and codewords.
- `transcripts.parquet` contains one row per detected RNA molecule with its exact x, y, and z coordinates, gene name, quality value, and cell assignment.
- `cell_boundaries.parquet` contains polygon vertex coordinates defining the outer boundary of each segmented cell.
- `nucleus_boundaries.parquet` contains polygon vertex coordinates defining the boundary of each segmented nucleus.
- `cell_feature_matrix.h5` contains the gene-by-cell sparse count matrix in the standard 10x HDF5 format.
- `cells.zarr.zip` contains a Zarr store with z-level assignments and nucleus counts per cell.

The full dataset contains 161,000 cells and 480 genes across two fields of view from a single lung cancer tissue section. The tissue contains a mixture of tumor epithelial cells, stromal fibroblasts, immune infiltrates, endothelial cells, and other lung-resident populations.

---

## Requirements

The notebook runs in Google Colab. All dependencies are installed automatically in the first cell. The packages required are `spatialdata`, `spatialdata-io`, `spatialdata-plot`, `squidpy`, `scanpy`, `seaborn`, and `matplotlib`. No local installation is needed. No data needs to be manually downloaded or uploaded. The notebook handles everything.

---

## How to Run

1. Go to colab.research.google.com
2. Select File then Upload notebook and upload the `.ipynb` file from this repository
3. Select Runtime then Run all
4. The notebook will install packages, download the dataset, and run the full pipeline automatically
5. Total runtime is approximately 30 to 45 minutes on a standard Colab CPU instance, with the co-occurrence computation being the longest step

---

## Notebook Structure and Detailed Explanation

### Package Installation

The first code cell installs all required packages using pip. This is necessary in Colab because none of the spatial omics libraries are pre-installed in the default environment. Installation takes approximately two to three minutes.

---

### Imports

The imports follow the official squidpy Xenium tutorial exactly. `spatialdata` provides the data container. `spatialdata_io.xenium` provides the file reader. `scanpy` provides all single-cell preprocessing and clustering routines. `squidpy` provides all spatial analysis routines. `seaborn` and `matplotlib` are used for all plotting.

---

### Data Download and Extraction

The notebook downloads the 265 MB zip archive directly from the 10x Genomics public server using `wget`. It then opens the zip file using Python's built-in `zipfile` module and extracts only the seven required files, skipping the large TIFF images. The download and extraction are both guarded by existence checks so that re-running the notebook after a kernel restart does not re-download anything.

---

### Data Loading and Zarr Conversion

The `xenium()` reader from `spatialdata-io` parses all seven files simultaneously and assembles them into a single `SpatialData` object. It reads the pixel size from `experiment.xenium` and converts all coordinates from pixels to micrometers. It constructs polygon shapes from the boundary parquet files, registers all spatial elements to a shared coordinate system called `global`, and populates `obsm["spatial"]` from the cell centroid coordinates automatically.

The SpatialData object is written to disk as a Zarr store. Zarr is a chunked, compressed array format optimized for large multi-dimensional data. Writing once and reading from the Zarr store subsequently is more efficient than re-parsing the raw files on every run. All downstream analysis reads exclusively from the Zarr store.

**Output of `sdata` after reading from Zarr:** The printed summary shows all spatial elements in the object. It confirms the presence of the `cell_boundaries` shapes layer with 161,000 polygons, the `nucleus_boundaries` shapes layer, the `transcripts` points layer with over 40 million individual RNA molecule detections as a three-dimensional point cloud, the `cell_circles` shapes layer, and the `table` AnnData layer with the 161,000-by-480 count matrix. All elements are registered to the `global` coordinate system.

---

### AnnData Object Inspection

The count matrix and all cell-level metadata are stored in the `table` slot of the SpatialData object as an `AnnData` object. This is the primary object used for all single-cell analyses.

**Output of `adata`:** The AnnData summary confirms 161,000 cells and 480 genes. The `obs` slot contains per-cell metadata columns including `cell_id`, `transcript_counts`, `control_probe_counts`, `control_codeword_counts`, `unassigned_codeword_counts`, `deprecated_codeword_counts`, `total_counts`, `cell_area`, `nucleus_area`, `region`, `z_level`, `nucleus_count`, and `cell_labels`. The presence of `spatial` in `obsm` confirms that coordinates were set by the reader.

**Output of `adata.obs`:** A 161,000-row dataframe where each row is one segmented cell. The most important columns for quality control are `total_counts`, `cell_area`, and `nucleus_area`. The control probe and codeword count columns are assay quality indicators. In the official output, the first few rows show cells with total counts ranging from roughly 55 to 188 transcripts and cell areas ranging from roughly 45 to 200 square micrometers, illustrating the considerable heterogeneity in cell size and transcriptional activity across the tissue.

**Output of `adata.obsm["spatial"]`:** A 161,000-by-2 numpy array of floating-point centroid coordinates in micrometers. The x coordinates span roughly 100 to 7,200 micrometers corresponding to the two fields of view. This array is used by squidpy for all spatial distance calculations and graph construction.

---

### Quality Control

**Output of the control probe percentage cell:** Two percentages are printed. The official result shows a negative DNA probe count of approximately 0.0051 percent and a negative decoding count of approximately 0.0025 percent. Both are far below 1 percent, confirming a high signal-to-noise assay. Control probes are synthetic sequences that should not bind any real transcript, so their detection represents pure instrument background. Values this low indicate that the overwhelming majority of detected transcripts are genuine gene expression signal.

**Output of the QC histogram figure:** A four-panel figure that is one of the most important outputs in the notebook because it informs filtering.

The total transcripts per cell panel shows a right-skewed distribution. Most cells have between 50 and 300 total transcript detections. The tail on the right represents highly transcriptionally active cells. A small spike near zero represents empty segmentation regions that will be removed by filtering.

The unique transcripts per cell panel shows the number of distinct genes detected per cell. Comparing the shape of this distribution to the total counts distribution reveals whether the data is driven by a broad or narrow set of genes.

The area of segmented cells panel shows the physical size distribution of all cells as determined by the Xenium segmentation algorithm in square micrometers. The distribution is roughly log-normal. Very small values at the left correspond to segmentation fragments or debris. Very large values may indicate merged cells or large cell types such as macrophages.

The nucleus ratio panel shows the fraction of each cell's total area that is nucleus. A ratio close to 1 indicates a cell where almost no cytoplasm was detected outside the nucleus boundary, which can occur in tightly packed tissue. A low ratio indicates a cell with extensive cytoplasm, characteristic of fibroblasts and certain epithelial cells.

---

### Filtering

Cells with fewer than 10 total transcript counts are removed using `scanpy.pp.filter_cells`. Genes detected in fewer than 5 cells are removed using `scanpy.pp.filter_genes`. These thresholds are intentionally conservative and are chosen by visual inspection of the histogram distributions. Alternative filtering criteria that can be applied include minimum cell area, DAPI signal intensity, or minimum unique transcript count.

---

### Normalization, Dimensionality Reduction, and Clustering

The raw counts are saved to `adata.layers["counts"]` before normalization to preserve the original integers for any analysis that requires them.

`scanpy.pp.normalize_total` scales every cell to the same total count, correcting for differences in transcript capture efficiency between cells.

`scanpy.pp.log1p` applies a log-plus-one transformation, compressing the dynamic range of expression values and making the data more suitable for linear methods.

`scanpy.pp.pca` reduces the 480-gene expression space to principal components that capture the major axes of transcriptional variation.

`scanpy.pp.neighbors` builds a k-nearest neighbor graph in PCA space, connecting each cell to its most transcriptionally similar neighbors.

`scanpy.tl.umap` computes a two-dimensional embedding of the neighbor graph that preserves local neighborhood structure.

`scanpy.tl.leiden` partitions the neighbor graph into communities. Each cell receives a cluster label in `adata.obs["leiden"]` representing its transcriptional population, which serves as a cell-type proxy throughout all spatial analyses.

---

### UMAP Visualization

**Output of the UMAP plot:** A three-panel figure where each point is a cell positioned by its UMAP coordinates.

The first panel, colored by total transcript count, checks whether library size is a major driver of the transcriptional structure. In a well-processed dataset this should be relatively uniform across the UMAP.

The second panel, colored by unique gene count, shows whether clusters differ in transcriptional complexity.

The third panel, colored by Leiden cluster, is the key result. In the official output from the full dataset this reveals approximately 20 distinct transcriptional populations forming separated islands, reflecting the diverse cell types present in lung cancer tissue. The number, size, and separation of islands in UMAP space reflects how transcriptionally distinct each cell population is from the others.

---

### Spatial Scatter Visualization

**Output of the spatial scatter plot:** Cells plotted at their true tissue x-y coordinates, colored by Leiden cluster. This is the central visualization of the notebook. It shows whether the transcriptionally defined populations have spatial organization.

In the official output, the spatial scatter clearly shows spatially coherent regions corresponding to different tissue compartments. Tumor epithelial clusters form compact regions consistent with tumor nests visible histologically. Stromal clusters fill the interstitial spaces. Immune clusters appear as focal infiltrates or distributed along boundaries. The alignment between transcriptional cluster identity and spatial tissue architecture confirms that Leiden clustering is capturing biologically real cell types rather than technical artifacts.

---

### Spatial Neighborhood Graph Construction

`squidpy.gr.spatial_neighbors` with `coord_type="generic"` and `delaunay=True` constructs a Delaunay triangulation of the cell centroid coordinates. Delaunay triangulation connects each cell to its natural geometric neighbors without any arbitrary distance cutoff. The resulting sparse connectivity and distance matrices stored in `adata.obsp` form the mathematical foundation for all three spatial statistics that follow.

---

### Centrality Scores

**Output of the centrality scores plot:** A three-panel bar chart with one bar per Leiden cluster.

The closeness centrality panel shows which clusters occupy central positions in the tissue. A high score means cells of that cluster are on average a short path length away from cells of all other types in the spatial graph. In lung tissue, stromal and vascular populations typically score highest because they are distributed throughout the tissue rather than concentrated in one region.

The degree centrality panel shows the fraction of non-cluster cells that directly border cells of each cluster. Clusters at the interface between major tissue compartments, such as immune cells at tumor boundaries or peritumoral fibroblasts, tend to score highest here.

The clustering coefficient panel shows how tightly cells of each cluster aggregate with each other. A score near 1 means a cell's spatial neighbors are almost exclusively other cells from the same cluster, indicating a compact, homogeneous aggregate. Tumor epithelial clusters almost always show the highest clustering coefficients because they form solid nests.

---

### Co-occurrence Probability

A 50 percent random subsample is taken before computing co-occurrence to reduce the computational cost of this analysis. The subsample is stored as a separate table in the SpatialData object, demonstrating the multi-table capability of the SpatialData format where different analytical subsets can coexist.

**Output of the co-occurrence plot for cluster 12:** A multi-line plot where the x-axis is distance in micrometers and the y-axis is the co-occurrence ratio for each cluster relative to cluster 12. The ratio at any given distance answers the question of how much more or less likely you are to find a cell from that cluster near cluster 12 compared to what would be expected if cells were randomly distributed.

Values above 1 indicate spatial enrichment at that distance. Peaks in the curves indicate characteristic interaction distances. A peak at 20 to 50 micrometers reflects direct cellular contact. A peak at 100 to 200 micrometers reflects compartment-level co-localization. Values consistently below 1 across all distances indicate spatial exclusion between that cluster and cluster 12.

**Output of the spatial scatter of the subsample:** The subsampled cells at tissue coordinates colored by Leiden cluster, providing the spatial context in which to interpret the co-occurrence curves.

---

### Neighborhood Enrichment

**Output of the neighborhood enrichment figure:** A side-by-side layout with a heatmap on the left and a spatial scatter on the right.

The heatmap shows every Leiden cluster on both axes. The color of each cell encodes a z-score derived by comparing observed spatial neighbor contacts to the distribution from 1,000 permutations of cluster labels. A warm color (positive z-score) means those two clusters are spatial neighbors more often than chance. A cool color (negative z-score) means they avoid each other spatially.

This heatmap is the most information-dense output in the notebook. It encodes the complete pairwise spatial interaction landscape of all cell types simultaneously. Strong positive off-diagonal values reveal specific cell-type pairs that share physical boundaries in the tissue, which in a tumor context reflects direct tumor-immune contacts, epithelial-stromal interfaces, and vascular niches. Strong negative values reflect segregated populations that occupy distinct compartments.

The spatial scatter on the right provides the visual confirmation of the patterns encoded in the heatmap.

---

### Moran's I Spatial Autocorrelation

**Output of the Moran's I table:** A dataframe of the top 10 most spatially autocorrelated genes. The columns include the Moran's I statistic, p-values under normality and permutation assumptions, and Benjamini-Hochberg corrected p-values.

In the official output from the full dataset, the top 10 genes are AREG at 0.696, MET at 0.683, ANXA1 at 0.667, EPCAM at 0.633, DMBT1 at 0.588, IGKC at 0.549, IGHG1 at 0.518, IDO1 at 0.496, SPARC at 0.446, and APOE at 0.437. All have p-values of zero under both testing frameworks, indicating these spatial patterns are far beyond chance.

The biological interpretation of these results is significant. AREG and MET are receptor-ligand signaling molecules associated with lung adenocarcinoma. EPCAM is a canonical epithelial cell marker concentrated in tumor nests. IGKC and IGHG1 are immunoglobulin genes expressed by B cells and plasma cells concentrated in immune aggregates. IDO1 is an immunosuppressive enzyme expressed in the tumor microenvironment. SPARC is a stromal extracellular matrix protein. APOE is expressed by macrophages. The fact that these genes span multiple cell types and compartments, and that all are highly spatially autocorrelated, confirms that the spatial structure of the tissue is biologically meaningful and faithfully captured by the Xenium assay.

**Output of the AREG and MET spatial scatter:** The normalized expression level of AREG and MET plotted at each cell's tissue coordinates using a continuous color scale. Both genes show sharply localized high expression in specific spatial regions corresponding to the tumor epithelial compartment visible in the spatial Leiden scatter. The spatial concentration of these oncologically relevant genes directly reflects the spatial organization of the tumor in the tissue section.

---

### Morphology Image Overlay

This cell uses `spatialdata-plot` to render AREG and MET expression overlaid on the raw morphology focus image from the Xenium instrument. The cell circles are colored by expression level and composited on the grayscale morphology image. This requires the morphology TIFF images which are not included in the 265 MB download used in this notebook. The cell should be skipped when running with the 2-FOV dataset. The full dataset produces figures showing expression concentrated in glandular tumor structures consistent with lung adenocarcinoma morphology.

---

### Interactive Visualization with napari-spatialdata

`napari-spatialdata` provides an interactive GUI for exploring SpatialData objects. Calling `Interactive(sdata)` opens a multi-layer viewer where spatial elements can be toggled, panned, zoomed, and inspected at any scale. Gene expression can be overlaid on morphology images. This requires a desktop display server and cannot be run in Colab. The official tutorial includes a screenshot showing AREG expression visualized across all cells in the napari viewer.

---

## Key Results Summary

The control probe percentages confirm a high-quality Xenium run with background signal below 0.01 percent of all detected transcripts.

UMAP embedding and Leiden clustering reveal approximately 20 distinct transcriptional populations reflecting the known cellular heterogeneity of human lung cancer tissue.

The spatial scatter visualization confirms that the transcriptionally defined populations have strong spatial coherence, with distinct compartments corresponding to tumor nests, stroma, and immune infiltrates.

Centrality analysis quantifies which populations are spatially central, boundary-dwelling, or tightly aggregated, providing a summary of the architectural role of each cell type.

Co-occurrence analysis quantifies the spatial interaction distances between cell populations, revealing which types are in direct contact and which co-localize at the compartment level.

Neighborhood enrichment provides a complete pairwise spatial interaction map showing which cell types preferentially neighbor each other and which are spatially excluded from each other.

Moran's I identifies AREG and MET as the two most spatially autocorrelated genes with values above 0.68, with EPCAM, DMBT1, and immunoglobulin genes also in the top 10, all confirmed by permutation testing to be far beyond chance.

---

## Repository Structure

```
.
├── README.md
├── tutorial_xenium_colab_v3.ipynb
└── Images/
    ├── qc_histograms.png
    ├── umap_clusters.png
    ├── spatial_leiden.png
    ├── centrality_scores.png
    ├── co_occurrence.png
    ├── spatial_subsample.png
    ├── nhood_enrichment.png
    └── spatial_areg_met.png
```

---

## References

- Wolf F A, Angerer P, Theis F J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018.
- Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. Nature Methods. 2022.
- Marconato L, Palla G, Yamauchi K A, et al. SpatialData: an open and universal data framework for spatial omics. Nature Methods. 2024.
- 10x Genomics Xenium platform: https://www.10xgenomics.com/platforms/xenium
- Official squidpy Xenium tutorial: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html
