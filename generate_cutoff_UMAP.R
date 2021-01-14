## Since we don't have ground truth labels for our data or the public data,
## we have to use more qualitative metrics to assess if our method is "better"
## than a simple cutoff at percent mitochondrial reads. This script generates
## a UMAP plot of the unfiltered data, to see if the cells being thrown out
## cluster together or are spread throughout the UMAP space. Then, it generates
## a new UMAP embedding using only the kept data and plots that. Ideally, we
## would be able to see some qualitative improvement in the "cleanness" of the
## model UMAP compared to the simple cutoff UMAP.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(scater)
  library(scran)
})

source("config.R")

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "16030X4"),
  make_option(c("-p", "--percent"), type = "numeric", default = 20)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# load data
location <- paste(path_to_mito_filtering, "/sce_objects/pre_filtering/",
                  opt$file, ".rds", sep = "")
sce <- readRDS(location)


# Get original embedding
set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

sce$cutoff_keep <- sce$subsets_mito_percent < opt$percent
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_UMAP_cutoff_",
          opt$percent, ".png", sep = ""), width = 800, height = 800)
plotUMAP(sce, colour_by = "cutoff_keep") +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  ggtitle("Cutoff at 10%") + theme(text = element_text(size = 20))
dev.off()



# Get cutoff-kept-only embedding
set.seed(617)
cutoff <- sce[, sce$cutoff_keep == T]
cutoff <- runUMAP(cutoff,
                  BNPARAM = BiocNeighbors::AnnoyParam(),
                  BPPARAM = BiocParallel::MulticoreParam(),
                  min_dist = 0.5,  repulsion_strength = 0.25,
                  spread = 0.7,
                  n_neighbors = min(15, ncol(cutoff)))

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_UMAP_cutoff_",
          opt$percent, "_only.png", sep = ""), width = 800, height = 800)
plotUMAP(cutoff)
dev.off()

set.seed(617)
g <- buildSNNGraph(cutoff,
                   k = 30,   # higher = bigger clusters
                   BNPARAM = BiocNeighbors::AnnoyParam(),
                   BPPARAM = BiocParallel::MulticoreParam())
clusters <- as.factor(igraph::cluster_louvain(g)$membership)
cutoff$clusters <- clusters
saveRDS(cutoff, paste(path_to_mito_filtering, "/sce_objects/", opt$file,
                      "_clustering_cutoff_", opt$percent, ".rds", sep = ""))
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_clustering_cutoff_",
          opt$percent, ".png", sep = ""), width = 800, height = 800)
plotUMAP(cutoff, colour_by = "clusters")
dev.off()
