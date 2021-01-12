## When running compare_cellranger_alevin_results.R, I noticed that the same
## samples had very different percent mitochondrial expression. This script
## is investigating that. The conclusion is that, because there were fewer
## mitochondrial genes in alevin than cellranger, the reads to the missing mito
## genes were instead mapping to chromosomal mimic pseudogenes and thus not
## getting counted for percent mitochondrial reads. Ultimately a problem of
## different reference transcriptomes. Should only need to be a one-off script.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(data.table)
  library(Seurat)
  library(scater)
  library(dplyr)
})

source("config.R")

#Get alevin and cellranger runs of a single sample, standardize barcode format
alevin <- readRDS("sce_objects/pre_filtering/H117_pTum1_alevin.rds")
alevin$barcode <- gsub(".*_", "", colnames(alevin))
cellranger <- readRDS("sce_objects/pre_filtering/H117_pTum1_cellranger.rds")
cellranger$barcode <- gsub("-[0-9]", "", colnames(cellranger))

#Pick a barcode with different but fairly high mitochondrial expression in both
cell_barcode <- "AAAGCAACAAACCCAT"
colData(alevin)[alevin$barcode == "AAAGCAACAAACCCAT", ]
colData(cellranger)[cellranger$barcode == "AAAGCAACAAACCCAT", ]

#Look at all non-zero expressed genes in these cells
which(alevin$barcode == "AAAGCAACAAACCCAT")
alevin_cell <- assay(alevin[, 954])
alevin_cell <- as.data.frame(alevin_cell)
setnames(alevin_cell, "expression")
alevin_cell <- subset(alevin_cell, alevin_cell$expression > 0)
sum(alevin_cell$expression)

#Look at which mito genes are most expressed in cellranger
which(cellranger$barcode == "AAAGCAACAAACCCAT")
cellranger_cell <- assay(cellranger[, 2])
cellranger_cell <- as.data.frame(cellranger_cell)
setnames(cellranger_cell, "expression")
cellranger_cell <- subset(cellranger_cell, cellranger_cell$expression > 0)
sum(cellranger_cell$expression)

#Try with another cell, one that gets treated differently by the miQC runs
which(alevin$barcode == "CTCACACCAGACGCCT")
alevin_cell <- assay(alevin[, 559])
alevin_cell <- as.data.frame(alevin_cell)
setnames(alevin_cell, "expression")
alevin_cell <- subset(alevin_cell, alevin_cell$expression > 0)
sum(alevin_cell$expression)

which(cellranger$barcode == "CTCACACCAGACGCCT")
cellranger_cell <- assay(cellranger[, 651])
cellranger_cell <- as.data.frame(cellranger_cell)
setnames(cellranger_cell, "expression")
cellranger_cell <- subset(cellranger_cell, cellranger_cell$expression > 0)
sum(cellranger_cell$expression)
