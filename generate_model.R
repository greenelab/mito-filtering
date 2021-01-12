## In single-cell RNA-seq, a lot of cells will die during the dissocation
## process or have been dead already, meaning their RNA is degraded and
## shouldn't be included in downstream analysis. One sign of this is if
## most of the RNA reads are coming from mtDNA, which is usually protected by
## the mitochondrial membrane and thus degrades slower. People usually set some
## rule like "exclude all cells with >15% mitochondrial reads," but what cutoff
## works varies widely by experiment and we'd like a more data-driven way to
## make the decision of if a cell is dead and to exclude it. We hypothesize any
## experiment will have two distributions of total genes expressed and percent
## mitochondrial reads, one of healthy cells and one of dead cells. This script
## creates a mixture model to estimate these distributions and select cells to
## remove based on posterior probability of coming from the dead distribution.

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(flexmix)
  library(splines)
})

source("config.R")

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "16030X2"),
  make_option(c("-m", "--model"), type = "character", default = "linear"),
  make_option(c("-c", "--cutoff"), type = "numeric", default = 0.75)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# load data
location <- paste(path_to_mito_filtering,
                  "/sce_objects/pre_filtering/", opt$file, ".rds", sep = "")
sce <- readRDS(location)
metrics <- as.data.frame(colData(sce))

# visualize data distribution
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_loess.png", sep = ""),
    width = 1200, height = 800)
ggplot(metrics, aes(x = detected, y = subsets_mito_percent)) +
  geom_point(colour = "#33adff") +
  stat_smooth(method = "loess", formula = y~x, size = 2, se = F,
              colour = "black") +
  labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
  ggtitle(opt$file) + theme(text = element_text(size = 20))
dev.off()

# run flexmix to get mixture model
if (opt$model == "linear") {
  model <- flexmix(subsets_mito_percent~detected, data = metrics, k = 2)
} else if (opt$model == "spline") {
  model <- flexmix(subsets_mito_percent~bs(detected), data = metrics, k = 2)
} else if (opt$model == "polynomial") {
  model <- flexmix(subsets_mito_percent~poly(detected, degree = 2), data = metrics, k = 2)
}

# plot model parameters on data
fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_fit.png", sep = ""), width = 1200, height = 800)
ggplot(metrics, aes(x = detected, y = subsets_mito_percent)) +
  geom_point(colour = "#33adff") +
  geom_line(data = fitted_models, inherit.aes = F, aes(x = V1, y = Comp.1), lwd = 2) +
  geom_line(data = fitted_models, inherit.aes = F, aes(x = V1, y = Comp.2), lwd = 2) +
  labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
  ggtitle(opt$file) + theme(text = element_text(size = 20)) +
  ylim(0, max(metrics$subsets_mito_percent))
dev.off()

# Determine which model component is the dead cell distribution,
# since flexmix doesn't number the distributions deterministically.
# Dead cell distribution should always have a higher y-intercept than
# the alive cell distribution, no matter the kind of model used.
intercept1 <- parameters(model, component = 1)[1]
intercept2 <- parameters(model, component = 2)[1]
if (intercept1 > intercept2) {
  dead_cell_dist <- 1
} else {
  dead_cell_dist <- 2
}

# calculate posterior likelihood for each cell
post <- posterior(model)
metrics$probability_dead <- post[, dead_cell_dist]

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_posterior.png", sep = ""), width = 1200, height = 800)
ggplot(metrics,
  aes(x = detected, y = subsets_mito_percent, colour = probability_dead)) +
  labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
  geom_point() +
  ggtitle(opt$file) + theme(text = element_text(size = 20)) +
  geom_line(data = fitted_models, inherit.aes = F, aes(x = V1, y = Comp.1), lwd = 2) +
  geom_line(data = fitted_models, inherit.aes = F, aes(x = V1, y = Comp.2), lwd = 2) +
  ylim(0, max(metrics$subsets_mito_percent))
dev.off()

# determine which cells to keep
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_keep.png", sep = ""), width = 1200, height = 800)
metrics$keep <- metrics$probability_dead <= opt$cutoff
ggplot(metrics, aes(x = detected, y = subsets_mito_percent, colour = keep)) +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
  geom_point() + ggtitle(opt$file)  + theme(text = element_text(size = 20))
table(metrics$keep)
dev.off()

colData(sce)$probability_dead <- metrics$probability_dead
colData(sce)$keep <- metrics$keep

# save object with filtering data
filename <- paste(opt$file, "_", opt$model, ".rds", sep = "")
filepath <- paste(path_to_mito_filtering, "sce_objects",
                  "post_filtering", filename, sep = "/")
saveRDS(sce, file = filepath)
