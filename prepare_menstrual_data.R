## Casey stumbled across this public dataset with a lot of dead cells and we
## decided to run miQC on it. Since it was downloaded from a different way
## than the public datasets in the scRNAseq package, I made a new script.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tximport)
  library(scater)
  library(scran)
  library(ggplot2)
  library(flexmix)
  library(scRNAseq)
  library(biomaRt)
})

source("config.R")

counts <- readRDS("GSE111976_ct_endo_10x.rds")
sce <- SingleCellExperiment(list(counts = counts))
dim(sce)
rm(counts); gc()

# Get mito genes
mt_genes <- grepl("^MT-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sum(mt_genes)

# Run scater
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

sce <- addPerFeatureQC(sce, BPPARAM = BiocParallel::MulticoreParam())
rowData(sce)$gene <- rownames(rowData(sce))

# Remove worst cells
filt <- sce$sum > 500 & sce$detected > 100
length(which(!filt))
sce <- sce[, filt]

# Write to file
filepath <- paste(path_to_mito_filtering, "sce_objects",
                  "pre_filtering/Wang.rds", sep = "/")
saveRDS(sce, file = filepath)
