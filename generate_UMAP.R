## Since we don't have ground truth

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(splines)
  library(scater)
})

source("config.R")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="16030X2"),
  make_option(c("-m", "--model"), type="character", default="linear"), 
  make_option(c("-p", "--percent"), type="numeric", default=20)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# load data
location <- paste(path_to_mito_filtering, "/sce_objects/post_filtering/",
                  opt$file, "_", opt$model, ".rds", sep="")
sce <- readRDS(location)

sce <- logNormCounts(sce)

set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model,"_UMAP_1_full.png",sep=""),width=800, height=800)
plotUMAP(sce, colour_by="subsets_mito_percent")
dev.off()

justkeep <- sce[,sce$keep==T]
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model,"_UMAP_2_filtered.png",sep=""),width=800, height=800)
plotUMAP(justkeep)
dev.off()

cutoff <- sce[,sce$subsets_mito_percent<=opt$percent]
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model,"_UMAP_3_cutoff_",opt$cutoff,".png",sep="")
          ,width=800, height=800)
plotUMAP(cutoff)
dev.off()

set.seed(617)
justkeep <- runUMAP(justkeep,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model,"_UMAP_4_filtered_only.png",sep=""),width=800, height=800)
plotUMAP(justkeep)
dev.off()

set.seed(617)
cutoff <- runUMAP(cutoff,
                    BNPARAM = BiocNeighbors::AnnoyParam(),
                    BPPARAM = BiocParallel::MulticoreParam(),
                    min_dist = 0.5,  repulsion_strength = 0.25,
                    spread = 0.7,
                    n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model,"_UMAP_5_cutoff_",opt$cutoff,"_only.png",
          sep=""),width=800, height=800)
plotUMAP(cutoff)
dev.off()
