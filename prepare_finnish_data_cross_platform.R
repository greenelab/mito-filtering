## In order to compare the results from cellranger and alevin, a colleague at
## University of Helsinki resent their tumor data as quantified by the two
## platforms. This is converting the objects sent from Seurat to sce objects
## and processing them with scater for use in miQC.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(data.table)
  library(Seurat)
  library(scater)
  library(dplyr)
})

source("config.R")

# Open alevin results
alevin <- readRDS(paste(path_to_finnish_data,
                        "updated_finnish_data/alevin_counts_v2.Rds", sep = "/"))
alevin <- as.SingleCellExperiment(alevin)

# Convert ensembl ids to gene names to easily pull out mitochondrial genes
gene_map <- fread(paste(path_to_hgsc_data,
                        "index", "gene2symbol.tsv", sep = "/"), header = F)
setnames(gene_map, c("ensembl_id", "gene_symbol"))
gene_map$ensembl_id <- gsub("\\..*", "", gene_map$ensembl_id)
gene_map <- unique(gene_map)
rowData(alevin)$ensembl_id <- gsub("\\..*", "", rownames(alevin))
setdiff(rowData(alevin)$ensembl_id, gene_map$ensembl_id)
# TODO: there are 12 extra genes in their alevin run than in our gene map,
# which is weird since they used the index we gave them. Follow up on this.
alevin <- subset(alevin, rowData(alevin)$ensembl_id %in% gene_map$ensembl_id)
rowData(alevin) <- left_join(as.data.frame(rowData(alevin)), gene_map)
rownames(alevin) <- rowData(alevin)$gene_symbol
rowData(alevin)$gene_symbol <- NULL
rowData(alevin)$ensembl_id <- NULL

# Get mitochondrial genes
mt_genes <- grepl("^MT-",  rownames(alevin))
sum(mt_genes)
feature_ctrls <- list(mito = rownames(alevin)[mt_genes])

# Run scater QC calculations
alevin <- addPerCellQC(alevin, subsets = feature_ctrls,
                       BPPARAM = BiocParallel::MulticoreParam())
alevin <- addPerFeatureQC(alevin,
                          BPPARAM = BiocParallel::MulticoreParam())

# Filter out obviously bad cells
filt <- alevin$sum > 500 & alevin$detected > 100
length(which(!filt))
alevin <- alevin[, filt]

# Split data by sample of origin, save to file
table(alevin$orig.ident)
samples <- unique(alevin$orig.ident)

for (i in 1:4) {
  donor <- samples[i]
  x1 <- alevin[, alevin$orig.ident == donor]
  filename <- paste(donor, "_alevin.rds", sep = "")
  saveRDS(x1, file = paste(path_to_mito_filtering,
                           "/sce_objects/pre_filtering/",
                           filename, sep = ""))
}

# Load cellranger object
load(paste(path_to_finnish_data,
              "updated_finnish_data/cellranger_counts.RData", sep = "/"))
cellranger <- as.SingleCellExperiment(cellranger_counts)

# Get mitochondrial genes
mt_genes <- grepl("^MT-",  rownames(cellranger))
sum(mt_genes)
feature_ctrls <- list(mito = rownames(cellranger)[mt_genes])

# Run scater QC calculations
cellranger <- addPerCellQC(cellranger, subsets = feature_ctrls,
                           BPPARAM = BiocParallel::MulticoreParam())
cellranger <- addPerFeatureQC(cellranger,
                              BPPARAM = BiocParallel::MulticoreParam())

# Remove obviously bad cells
filt <- cellranger$sum > 500 & cellranger$detected > 100
length(which(!filt))
cellranger <- cellranger[, filt]

# Divide data by sample of origin, write to file
cellranger$sample <- gsub("[ATCG]+-", "", colnames(cellranger))
sample_map <- c("H094_pOme", "H095_pOvaR", "H117_pTum1", "H122_pTum1")

for (i in 1:4) {
  donor <- sample_map[i]
  x1 <- cellranger[, cellranger$sample == i]
  filename <- paste(donor, "_cellranger.rds", sep = "")
  saveRDS(x1, file = paste(path_to_mito_filtering,
                           "/sce_objects/pre_filtering/",
                           filename, sep = ""))
}
