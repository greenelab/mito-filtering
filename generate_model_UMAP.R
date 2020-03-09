## Since we don't have ground truth labels for our data or the public data,
## we have to use more qualitative metrics to assess if our method is "better"
## than a simple cutoff at percent mitochondrial reads. This script generates
## a UMAP plot of the unfiltered data and shows which cells would be kept by
## our model, and what they would look like in a UMAP embedding on their own.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(scater)
})

source("config.R")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="16030X2"),
  make_option(c("-m", "--model"), type="character", default="linear")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load data
location <- paste(path_to_mito_filtering, "/sce_objects/post_filtering/",
                  opt$file, "_", opt$model, ".rds", sep="")
sce <- readRDS(location)


# Plot full UMAP embedding
set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/", opt$file,
          "_UMAP.png",sep=""),width=800, height=800)
plotUMAP(sce, colour_by="subsets_mito_percent")
dev.off()


# Plot full UMAP colored by which cells kept
justkeep <- sce[,sce$keep==T]
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_UMAP_", opt$model,
          ".png", sep=""),width=800, height=800)
plotUMAP(sce, colour_by="keep")
dev.off()


# Plot UMAP embedding with only kept cells
set.seed(617)
justkeep <- sce[,sce$keep==T]
justkeep <- runUMAP(justkeep,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_UMAP_", opt$model,
          "_only.png", sep=""), width=800, height=800)
plotUMAP(justkeep)
dev.off()

