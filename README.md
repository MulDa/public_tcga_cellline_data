Tumor-Type Specific MTAP Methylation Data
=========================================

This repository accompanies the manuscript *“Tumor-Type Specific Methylation Patterns of MTAP in human samples and cell lines”* and provides processed TCGA/TARGET and cell line data used for the DNA methylation and expression analyses.

The goal is to make the key processed data publicly available in a simple structure so that all figures and correlation analyses can be reproduced without re-running the full GDC/GEO pipelines.

## Repository structure

- `tcga/`
  - One subfolder per TCGA or TARGET project, for example:
    - `tcga/TCGA-COAD/`
      - `TCGA-COADnormalized_expression_data_picked_genes.RData`
      - `TCGA-COAD_MTAP.RData`
    - `tcga/TCGA-BRCA/`
      - `TCGA-BRCAnormalized_expression_data_picked_genes.RData`
      - `TCGA-BRCA_MTAP.RData`

  Each of these subfolders contains:

  - Variance-normalized RNA-seq expression matrix for a selected gene set (including MTAP), as used in the paper.
  - Processed methylation beta values at MTAP-associated CpG probes (e.g. `cg00230302`, `cg13492671`, `cg14548963`, `cg25162921`) for matched samples in that project.

- `cell_lines/`
  - `cell_lines.csv`  
    HM450K β-values for the four MTAP-associated CpG sites across 63 cancer and non-malignant cell lines.
  - `cell_line_info.csv`  
    Annotated table linking the methylation measurements to ENCODE/GEO metadata, including cell line name, tissue type, karyotype and other descriptive fields.

## TCGA data: download and processing

### Source

TCGA data were obtained from the Genomic Data Commons (GDC) using the `TCGAbiolinks` R package. For each project, raw data were downloaded into a project-specific folder under:

```text
GDCdata/TCGA-PROJECT/
```

Only samples with both Illumina HumanMethylation450K (HM450K) array data and RNA-seq gene expression quantification were retained for the integrated analyses.

### Expression data (`TCGA-PROJECTnormalized_expression_data_picked_genes.RData`)


- RNA-seq count data were queried with `TCGAbiolinks::GDCquery` and downloaded per project.
- Counts were normalized in R using **DESeq2** (variance-stabilizing size factor normalization).
- Ensembl IDs were mapped to HGNC symbols via **EnsDb.Hsapiens.v79**.
- A selected set of genes (including `MTAP` and extracellular matrix genes listed in the manuscript) was retained.
- For each project, the resulting normalized expression matrix was saved as:

```r
save(normalized_counts,
     file = paste0(projectname, "normalized_expression_data_picked_genes.RData"))
```

In all downstream analyses and figures, `MTAP_log2` was defined as `log2(MTAP + 1)` from this normalized matrix.

### Methylation data (`TCGA-PROJECT_MTAP.RData`)

- HM450K IDAT files were downloaded via `TCGAbiolinks` for each TCGA project.
- Raw intensities were imported using **minfi** (`read.metharray.exp`), and poor-quality probes were filtered using detection p-values (`detectionP < 0.01` in at least 80% of samples).
- Beta values were calculated as:

```text
β = M / (M + U + 100)
```

- Probe-type bias (Type I vs Type II probes) was corrected using **BMIQ** from the `wateRmelon` package.
- Probes on sex chromosomes, probes overlapping known SNPs, and known cross-reactive probes (Pidsley et al.) were removed.
- MTAP-associated CpG probes (`cg00230302`, `cg13492671`, `cg14548963`, `cg25162921`) were selected using manifest annotation.
- The final per-sample β-values for these probes were saved per project as:

```r
save(cpg_means, file = paste0(projectname, "_MTAP.RData"))
```

These `.RData` objects were then merged with the normalized MTAP expression to create the `correlation_table.RData` used for all correlation and plotting scripts.

### Sample filtering and matching

- Only `Primary Tumor` and `Solid Tissue Normal` samples were retained for the case–control comparisons in the boxplots and correlation plots.
- Sample IDs were harmonized using matching barcodes across expression and methylation tables.
- Matched samples were combined into a single per-project table containing:
  - `sample.submitter_id`
  - `project`
  - `sample_type` (`Primary Tumor` vs `Solid Tissue Normal`)
  - MTAP-associated CpG β-values
  - `MTAP_log2` expression

These combined tables were used to generate the figures showing β-value distributions and their correlation with MTAP expression, with clear case/control labelling and color coding (`Primary Tumor` in red, `Solid Tissue Normal` in blue).

## Cell line data (GEO GSE40699)

### Source

Cell line DNA methylation data were obtained from GEO series **GSE40699**, which provides HM450K profiles for 63 human cancer and non-malignant cell lines (ENCODE).

### Processing pipeline

The scripts `cell_lines.R` and `cell_line_info.R` describe the full workflow:

- Raw IDAT files were imported with **minfi** (`read.metharray.exp`).
- Detection p-value filtering (`detectionP < 0.01`) and BMIQ normalization (via `wateRmelon`) were applied as for TCGA.
- Sex-chromosome, SNP-affected and cross-reactive probes were removed using the same annotation lists as for TCGA.
- The four MTAP-associated CpG probes (`cg00230302`, `cg13492671`, `cg14548963`, `cg25162921`) were extracted and β-values compiled across all 63 cell lines, saved as:
  - `cell_lines.csv` – simple matrix of CpG β-values per cell line.
- GEO series metadata (`GSE40699_series_matrix.txt`) were parsed to obtain cell-line level annotation (cell name, lineage, karyotype, treatment, etc.), which were merged with `cell_lines.csv` in `cell_line_info.R` to produce:
  - `cell_line_info.csv` – annotated cell line metadata plus MTAP CpG β-values.

### Note on cell line plots

Each cell line is represented by a single HM450K sample. Because there is no within-line replication, detailed correlation plots analogous to the TCGA/TARGET figures (with Pearson r and regression lines) are not statistically meaningful for cell lines and were therefore not generated.

## Using the data

### Loading TCGA project data in R

```r
setwd("path/to/public_tcga_cellline_data")

# Example: load TCGA-COAD methylation and expression
load("tcga/TCGA-COAD/TCGA-COAD_MTAP.RData")                      # MTAP β-values
load("tcga/TCGA-COAD/TCGA-COADnormalized_expression_data_picked_genes.RData")  # expression
```

You can then merge β-values and expression on sample IDs (e.g. `sample.submitter_id`) to reconstruct the `correlation_table` used for plotting.

### Loading cell line data

```r
cell_lines <- read.csv("cell_lines/cell_lines.csv")
cell_line_info <- read.csv("cell_lines/cell_line_info.csv")
```

These tables can be used to reproduce the cell-line methylation distributions (e.g. β-value boxplots for the four CpG sites) reported in the manuscript.

## Reproducibility notes

- All statistical analyses were performed in **R 4.4.0**.
- Key packages used include: `TCGAbiolinks`, `DESeq2`, `wateRmelon`, `RPMM`, `EnsDb.Hsapiens.v79`, `minfi`, and `ggplot2`.
- The exact scripts used to generate the `.RData` objects are the ones already present in the original analysis directory (`process_tcga_data.R`, `process_tcga_data_idat.R`, `process_tcga_correlation.R`, `cell_lines.R`, `cell_line_info.R`).

For questions about the data or code, please contact the corresponding authors listed in the manuscript.

