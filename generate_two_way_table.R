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

