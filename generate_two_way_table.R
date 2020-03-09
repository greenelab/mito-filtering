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

# create two-way table
sce$cutoff <- sce$subsets_mito_percent < opt$percent
x <- table(sce$cutoff,sce$keep)
rownames(x) <- c("Cutoff removes", "Cutoff keeps")
colnames(x) <- c("\tModel removes", "Model keeps")

# save to file
outfile <- paste(path_to_mito_filtering, "/tables/", opt$file,
                  "_", opt$model, "_cutoff_", opt$percent, ".txt", sep="")
write.table(x, file=outfile, quote=F, sep="\t")
