# ARG_HPC_snRNAseq
This project contains R code for analysis of single-nucleus RNA-sequencing data from mouse and human hippocampal neurons. 

# Description of analyses

## 00_check_read_lengths
- Script to check read 1 files for discrepant read lengths. Should be (and are) exactly 28bp = 16 [BC] + 12 [UMI]. 

## 01_align
- Scripts to align sequencing data from all samples using 10x Genomics’ cellranger count version 3.0.2. [here] 

## 02_processing
- Script to build initial SingleCellExperiment (sce) object, filter empty droplets, perform quality control, remove doublets, and build new HPC-only sce object.

## 03_clustering_annotation
- Script to cluster HPC nuclei, annotate cluster, and build figure 1 plots.

## 04_broad_de_analysis
- Script to perform initial, sample-level (non-cell type specific) pseudobulk differential expression analysis, perform GO over-representation analysis (ORA), and build figure 2 plots.

## 05_cell_type_specific_de_analysis
- Script to perform cell type specific pseudobulk differential expression analysis and build figure 3 plots.

## 06_gsea
- Script to run gene set enrichment analysis (GSEA) on cell type-specific DE analysis results for GO terms within biological process (BP), cellular component (CC), and molecular function (MF) domains.
- Script to process GSEA results, perform semantic similarity clustering of GO terms, and build figure 4 plots.
- Script to make heatmaps for figure 4 plots.

## 07_nmf
- Script to perform non-negative matrix factorization (NMF) using CoGAPS.
- Script to build figure 5 plots.

# Processed data
- File containing SingleCellExperiment object with all nuclei (pre-thalamic neuron removal) [here] https://github.com/Erik-D-Nelson/ARG_HPC_snRNAseq/blob/main/processed_data/sce_total.rda.xz
- File containing SingleCellExperiment object with HPC nuclei only [here] https://github.com/Erik-D-Nelson/ARG_HPC_snRNAseq/blob/main/processed_data/sce_subset.rda.xz
