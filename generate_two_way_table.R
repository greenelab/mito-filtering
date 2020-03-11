## This is a simple script to show how many cells are kept or lost using our
## filtering method as opposed to a simple mito percentage cutoff. 

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
sce$cutoff <- sce$subsets_mito_percent < opt$percent

colData(sce)$category<-factor("both_keep",levels=c("both_keep","model_keep","cutoff_keep","both_remove"))

if(sum(colData(sce)$keep==T & colData(sce)$cutoff==T)>0){
  print("yay!")
  colData(sce)[colData(sce)$keep==T & colData(sce)$cutoff==T,]$category<-"both_keep"
} 
if(sum(colData(sce)$keep==T & colData(sce)$cutoff==F)>0){
  print("yes!")
  colData(sce)[colData(sce)$keep==T & colData(sce)$cutoff==F,]$category<-"model_keep"
}
if(sum(colData(sce)$keep==F & colData(sce)$cutoff==T)>0){
  print("hao!")
  colData(sce)[colData(sce)$keep==F & colData(sce)$cutoff==T,]$category<-"cutoff_keep"
}
if(sum(colData(sce)$keep==F & colData(sce)$cutoff==F)>0){
  print("molto bene!")
  colData(sce)[colData(sce)$keep==F & colData(sce)$cutoff==F,]$category<-"both_remove"
}

# Plot full UMAP embedding
set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

png(paste(path_to_mito_filtering, "/plots/UMAP_", opt$file,
          "_", opt$model, "_cutoff_", opt$percent, 
          ".png",sep=""),width=800, height=800)
plotUMAP(sce,colour_by="category")
dev.off()

# add pseudocounts to each category. It's crude but it allows you to still make
# a two-way table even if all cells go one way by either metric.
sce_short <- subset(colData(sce),select=c("keep","cutoff"))
pseudocounts <- data.frame(keep=c("FALSE","FALSE","TRUE","TRUE"),cutoff=c("TRUE","FALSE","TRUE","FALSE"))
sce_short <- rbind(sce_short,pseudocounts)

# create two-way table
x <- table(sce_short$cutoff,sce_short$keep)
rownames(x) <- c("Cutoff removes", "Cutoff keeps")
colnames(x) <- c("\tModel removes", "Model keeps")

# remove pseudocounts
x <- x - 1

# save to file
outfile <- paste(path_to_mito_filtering, "/tables/", opt$file,
                  "_", opt$model, "_cutoff_", opt$percent, ".txt", sep="")
write.table(x, file=outfile, quote=F, sep="\t")

