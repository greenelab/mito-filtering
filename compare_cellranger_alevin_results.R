## Building off prepare_finnish_data_cross_platform.R. Different quantifying
## platforms like alevin and cellranger will return different results for
## number of quantified genes and percent mitochondrial reads. Here, we show
## that for (most) samples, using miQC as opposed to a flat cutoff will not
## only preserve more cells overall, but there will also be more concordance
## between which cells are kept across tools This script makes confusion
## matrices and calculates cross-platform consistency.

## Parameters: --file is the name of the tumor you want to look at, --model is
## the kind of mixture model used for miQC. At this point, it should only be
## the default linear, but we might try other models in the future.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(data.table)
  library(Seurat)
  library(scater)
  library(dplyr)
  library(caret)
})

source("config.R")

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "H094_pOme"),
  make_option(c("-m", "--model"), type = "character", default = "linear")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
rm(opt_parser, option_list)

# Get SingleCellExperiment objects of tumor samples with miQC filtering results
cellranger <- readRDS(paste(path_to_mito_filtering,
                            "/sce_objects/post_filtering/", opt$file,
                            "_cellranger_", opt$model, ".rds", sep = ""))
alevin <- readRDS(paste(path_to_mito_filtering,
                        "/sce_objects/post_filtering/", opt$file,
                        "_alevin_", opt$model, ".rds", sep = ""))

cellranger$barcode <- gsub("-[0-9]", "", colnames(cellranger))
cellranger_metrics <- as.data.frame(subset(colData(cellranger),
                                           select = c("barcode", "keep")))
setnames(cellranger_metrics, c("barcode", "cellranger_status"))
alevin$barcode <- gsub(".*_", "", colnames(alevin))
alevin_metrics <- as.data.frame(subset(colData(alevin),
                                       select = c("barcode", "keep")))
setnames(alevin_metrics, c("barcode", "alevin_status"))

# Make confusion matrix, where true positive = cell is not filtered out by
# miQC in either the cellranger or the alevin quantifications
twoway <- inner_join(alevin_metrics, cellranger_metrics)
twoway$cellranger_status <- as.factor(twoway$cellranger_status)
twoway$alevin_status <- as.factor(twoway$alevin_status)
miQC <- confusionMatrix(data = twoway$cellranger_status,
                        reference = twoway$alevin_status)


# Make basic cutoffs
cellranger$cutoff <- cellranger$subsets_mito_percent > 10
alevin$cutoff <- alevin$subsets_mito_percent > 10

cellranger_cutoff <- as.data.frame(subset(colData(cellranger),
                                          select = c("barcode", "cutoff")))
setnames(cellranger_cutoff, c("barcode", "cellranger_cutoff"))
alevin_cutoff <- as.data.frame(subset(colData(alevin),
                                      select = c("barcode", "cutoff")))
setnames(alevin_cutoff, c("barcode", "alevin_cutoff"))

# Make confusion matrix for cutoffs, where true positive = cell is below
# cutoff value in both cellranger and alevin quantifications
confusion <- inner_join(alevin_cutoff, cellranger_cutoff)
confusion$cellranger_cutoff <- as.factor(confusion$cellranger_cutoff)
confusion$alevin_cutoff <- as.factor(confusion$alevin_cutoff)
cutoff <- confusionMatrix(data = confusion$cellranger_cutoff,
                          reference = confusion$alevin_cutoff)

# Compare confusion matrices
miQC[[2]]
cutoff[[2]]

# Compare "accuracy" (concordance) across tools for miQC and cutoffs
miQC_concordance <- miQC[[3]][1]
cutoff_concordance <- cutoff[[3]][1]