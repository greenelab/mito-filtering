## Since we don't have ground truth labels for our data or the public data,
## we have to use more qualitative metrics to assess if our method is "better"
## than a simple cutoff at percent mitochondrial reads. This script generates
## a UMAP plot of the unfiltered data, what the same mapping would like with
## the two filtering methods, and what they would look like with their own
## UMAP embedding. It also provides a two-way table of number of cells
## included and excluded by each filtering method.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
})

source("config.R")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="16030X2"),
  make_option(c("-p", "--percent"), type="numeric", default=20)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load data
location <- paste(path_to_mito_filtering, "/sce_objects/pre_filtering/",
                  opt$file, ".rds", sep="")
sce <- readRDS(location)


# Get original embedding
set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

sce$cutoff_keep <- sce[,sce$subsets_mito_percent < opt$percent]
png(paste(path_to_mito_filtering, "/plots", opt$file, "_UMAP_cutoff_",
          opt$percent, ".png", sep=""),width=800, height=800)
plotUMAP(sce, colour_by="cutoff_keep")
dev.off()


# Get cutoff-kept-only embedding
set.seed(617)
cutoff <- sce[,sce$cutoff_keep==T]
cutoff <- runUMAP(cutoff,
                  BNPARAM = BiocNeighbors::AnnoyParam(),
                  BPPARAM = BiocParallel::MulticoreParam(),
                  min_dist = 0.5,  repulsion_strength = 0.25,
                  spread = 0.7,
                  n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots", opt$file, "_UMAP_cutoff_",
          opt$percent, "_only.png", sep=""),width=800, height=800)
plotUMAP(cutoff)
dev.off()
