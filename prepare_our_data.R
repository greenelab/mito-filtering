suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tximport)
  library(scater)
  library(scran)
})

source("config.R")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="16030X2")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load data into sce object
filename <- paste(path_to_hgsc_data, "tumors", opt$file, "alevin_output/alevin/quants_mat.gz",sep="/")
txia <- suppressMessages(tximport(
  paste(path_to_hgsc_data, "tumors", opt$file,
        "alevin_output/alevin/quants_mat.gz",sep="/"),
  type = "alevin"))
sce <- SingleCellExperiment(assays = list(counts = txia$counts))

# convert ensembl gene ids to gene symbols
gene_map <- fread(paste(path_to_hgsc_data,"index","gene2symbol.tsv",sep="/"), header = F)
setnames(gene_map, c("ensembl_id", "gene_symbol"))
rowData(sce)$ensembl_id <- rownames(sce)
rowData(sce) <- left_join(as.data.frame(rowData(sce)), gene_map)
rownames(sce) <- rowData(sce)$gene_symbol
rowData(sce)$gene_symbol <- NULL
rowData(sce)$ensembl_id <- NULL

# get list of mitochondrial genes
mt_genes <- grepl("^MT-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sum(mt_genes)

# calculate per-cell and per-feature metrics using scater
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

sce <- addPerFeatureQC(sce, BPPARAM = BiocParallel::MulticoreParam())
#rowData(sce)$Gene <- rownames(rowData(sce))
#sce

# remove unambiguously failed cells
filt <- sce$sum > 100 & sce$detected > 50
length(which(!filt))
sce<-sce[,filt]
sce

#save object
filename <- paste(opt$file, ".rds", sep="")
filepath <- paste(path_to_mito_filtering, "sce_objects",
                  "pre_filtering", filename, sep="/")
saveRDS(sce, file=filepath)
