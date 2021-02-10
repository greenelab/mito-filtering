suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(flexmix)
  library(cowplot)
  library(patchwork)
  library(scater)
  library(scran)
  library(mbkmeans)
  library(reshape2)
  library(viridis)
  library(splines)
  library(graphics)
  library(gridExtra)
  library(grid)
  library(caret)
  library(ggpubr)
  library(gtable)
  library(BiocParallel)
})

source("config.R")
source("utils.R")

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = "16030X4"),
  make_option(c("-m", "--model"), type = "character", default = "linear"),
  make_option(c("-p", "--posterior"), type = "numeric", default = 0.75),
  make_option(c("-c", "--cutoff"), type = "numeric", default = 10)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
rm(opt_parser, option_list)

# load data
location <- paste(path_to_mito_filtering,
                  "/sce_objects/pre_filtering/", opt$file, ".rds", sep = "")
sce <- readRDS(location)
metrics <- as.data.frame(colData(sce))


### FIGURE 1

# visualize data distribution
qc_spike2 <- isOutlier(metrics$subsets_mito_percent, type = "higher")
mad <- attr(qc_spike2, "thresholds")[2]
p1 <- plotNoModel(metrics) +
  geom_hline(yintercept = 5, lwd = line_width) +
  geom_hline(yintercept = mad, linetype = 2, lwd = line_width)

# run flexmix to get mixture model
model <- flexmix(subsets_mito_percent~detected, data = metrics, k = 2)

# plot model parameters on data
fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))

p2 <- plotDistributions(metrics, fitted_models)

metrics <- calculatePosterior(metrics, model)
p3 <- plotPosterior(metrics, fitted_models)

metrics$keep <- metrics$probability_dead <= opt$posterior
p4 <- plotKeepOrToss(metrics)

png("~/Documents/Figures/SmartQC Paper/Figure 1.png",
    width = 5200, height = 3200, res = 300)
p1 + p2 + p3 + p4 + plot_annotation(tag_levels = "A")
dev.off()


### SUPPLEMENTAL FIGURE 1

metrics <- as.data.frame(colData(sce))
model <- flexmix(subsets_mito_percent~bs(detected), data = metrics, k = 2)

# plot model parameters on data
fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))
k2 <- plotDistributions(metrics, fitted_models)

metrics <- calculatePosterior(metrics, model)
k3 <- plotPosterior(metrics, fitted_models)

metrics$keep <- metrics$probability_dead <= opt$posterior
k4 <- plotKeepOrToss(metrics)

png("~/Documents/Figures/SmartQC Paper/Supp Figure 1.png",
    width = 5200, height = 3200, res = 300)
k2 + k3 + k4 + plot_spacer() + plot_annotation(tag_levels = "A")
dev.off()


### FIGURE 2

run_model <- function(file, title) {
  location <- paste(path_to_mito_filtering,
                    "/sce_objects/pre_filtering/", file, ".rds", sep = "")
  sce <- readRDS(location)
  metrics <- as.data.frame(colData(sce))
  model <- flexmix(subsets_mito_percent~detected, data = metrics, k = 2)
  fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))

  metrics <- calculatePosterior(metrics, model)
  p <- plotPosterior(metrics, fitted_models, title)

  p
}

p1 <- run_model("Macosko", "Macosko et al")
p2 <- run_model("Shekhar", "Shekhar et al")
p3 <- run_model("Richard", "Richard et al")
p4 <- run_model("Zeisel", "Zeisel et al")
p5 <- run_model("Lawlor", "Lawlor et al")
p6 <- run_model("Wang", "Wang et al")

png("~/Documents/Figures/SmartQC Paper/Figure 2.png",
    width = 7200, height = 3200, res = 300)
p1 + p2 + p3 + p4 + p5 + p6 + plot_annotation(tag_levels = "A")
dev.off()


### FIGURE 3

r1 <- run_model("16030X2", "")
r2 <- run_model("16030X3", "")
r3 <- run_model("H094_pOme_alevin", "")
r4 <- run_model("H095_pOvaR_alevin", "")
r5 <- run_model("H117_pTum1_alevin", "")
r6 <- run_model("H122_pTum1_alevin", "")

png("~/Documents/Figures/SmartQC Paper/Figure 3.png",
    width = 7200, height = 3200, res = 300)
r1 + r2 + r3 + r4 + r5 + r6 + plot_annotation(tag_levels = "A")
dev.off()


### FIGURE 4

# Plot data as quantified by cellranger
location <- paste(path_to_mito_filtering, "/sce_objects/post_filtering/",
                  "H117_pTum1_cellranger_linear.rds", sep = "")
cellranger <- readRDS(location)
metrics <- as.data.frame(colData(cellranger))
model <- flexmix(subsets_mito_percent~detected, data = metrics, k = 2)
fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))

t1 <- plotDistributions(metrics, fitted_models, title = "Cellranger")
t5 <- plotNoModel(metrics, title = "Cellranger") +
  geom_hline(yintercept = opt$cutoff, lwd = line_width) + ylim(1, 100)

# Plot data as quantified by alevin
location <- paste(path_to_mito_filtering, "/sce_objects/post_filtering/",
                  "H117_pTum1_alevin_linear.rds", sep = "")
alevin <- readRDS(location)
metrics <- as.data.frame(colData(alevin))
t2 <- plotDistributions(metrics, fitted_models, title = "Alevin")
t6 <- plotNoModel(metrics, title = "Alevin") +
  geom_hline(yintercept = opt$cutoff, lwd = line_width) + ylim(1, 100)

# Standardize barcodes
cellranger$barcode <- gsub("-[0-9]", "", colnames(cellranger))
cellranger_metrics <- as.data.frame(subset(colData(cellranger),
                                           select = c("barcode", "keep")))
setnames(cellranger_metrics, c("barcode", "cellranger_status"))
alevin$barcode <- gsub(".*_", "", colnames(alevin))
alevin_metrics <- as.data.frame(subset(colData(alevin),
                                       select = c("barcode", "keep")))
setnames(alevin_metrics, c("barcode", "alevin_status"))

make_confusion_matrix <- function(tbl, title) {
  title <- text_grob(title, just = "center", size = 20, face = "bold")
  cellranger <- text_grob("Keep in Cellranger", just = "center", size = 18)
  alevin <- text_grob("Keep in Alevin", just = "center", size = 18)
  T23 <- text_grob("FALSE", just = "center", size = 16)
  T24 <- text_grob("TRUE", just = "center", size = 16)
  T32 <- text_grob("FALSE", just = "center", size = 16)
  T33 <- text_grob(tbl[1, 1], just = "center", size = 16)
  T34 <- text_grob(tbl[1, 2], just = "center", size = 16)
  T42 <- text_grob("TRUE", just = "center", size = 16)
  T43 <- text_grob(tbl[2, 1], just = "center", size = 16)
  T44 <- text_grob(tbl[2, 2], just = "center", size = 16)

  # Construct a gtable
  leg <- gtable(width = unit(c(4.6, 2.8, 2.8, 2.8), "cm"),
               height = unit(c(2, 1.5, 1, 1, 1), "cm"))

  # Place the grobs into the table
  leg <- gtable_add_grob(leg, cellranger, t = 2, l = 3, r = 4)
  leg <- gtable_add_grob(leg, alevin, t = 4, b = 5, l = 1)
  leg <- gtable_add_grob(leg, T23, t = 3, l = 3)
  leg <- gtable_add_grob(leg, T24, t = 3, l = 4)
  leg <- gtable_add_grob(leg, T32, t = 4, l = 2)
  leg <- gtable_add_grob(leg, T33, t = 4, l = 3)
  leg <- gtable_add_grob(leg, T34, t = 4, l = 4)
  leg <- gtable_add_grob(leg, T42, t = 5, l = 2)
  leg <- gtable_add_grob(leg, T43, t = 5, l = 3)
  leg <- gtable_add_grob(leg, T44, t = 5, l = 4)
  leg <- gtable_add_grob(leg, title, t = 1, l = 1, r = 4)

  leg
}

# Make confusion matrix, where true positive = cell is not filtered out by
# miQC in either the cellranger or the alevin quantifications
twoway <- inner_join(alevin_metrics, cellranger_metrics)
twoway$cellranger_status <- as.factor(twoway$cellranger_status)
twoway$alevin_status <- as.factor(twoway$alevin_status)
miQC <- confusionMatrix(data = twoway$cellranger_status,
                        reference = twoway$alevin_status,
                        dnn = c("Kept in Cellranger",
                              "Kept in Alevin"))

t3 <- make_confusion_matrix(miQC[[2]], "miQC filtering")

# Make basic cutoffs
cellranger$cutoff <- cellranger$subsets_mito_percent > opt$cutoff
alevin$cutoff <- alevin$subsets_mito_percent > opt$cutoff

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
                          reference = confusion$alevin_cutoff,
                          dnn = c("Kept in Cellranger",
                                "Kept in Alevin"))
t4 <- make_confusion_matrix(cutoff[[2]], "10% cutoff")

png("~/Documents/Figures/SmartQC Paper/Figure 4.png",
    width = 6400, height = 3200, res = 300)
t1 + t2 + t3 + t5 + t6 + t4
dev.off()


### SUPPLEMENTAL FIGURE 2

set.seed(621)
sce <- scater::runPCA(sce, ncomponents = 50,
                      ntop = 1000,
                      scale = TRUE,
                      BSPARAM = BiocSingular::RandomParam(),
                      BPPARAM = MulticoreParam(6))

ks <- seq(2, 10)
res <- lapply(ks, function(k) {
  mbkmeans(sce, clusters = k,
           reduceMethod = "PCA",
           calc_wcss = TRUE, num_init = 10)
})

wcss <- sapply(res, function(x) sum(x$WCSS_per_cluster))
elbow <- as.data.frame(cbind(ks, wcss))

png("~/Documents/Figures/SmartQC Paper/Supp Figure 2.png",
    width = 2600, height = 1800, res = 300)
ggplot(elbow, aes(x = ks, y = wcss)) +
  labs(x = "Number of clusters", y = "Within cluster sum of squares") +
  geom_line() + geom_point() +
  theme_bw() +
  theme(text = element_text(size = text_size))
dev.off()

### FIGURE 5

set.seed(617)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

sce$`% Mito` <- sce$subsets_mito_percent

q1 <- plotUMAP(sce, colour_by = "% Mito") +
  theme(text = element_text(size = text_size),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
q2 <- plotUMAP(sce, colour_by = "Keep") +
  scale_fill_manual(name = "Keep", values = c("#999999", "#E69F00")) +
  theme(text = element_text(size = text_size),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

sce$`Keep at 10% cutoff` <- sce$subsets_mito_percent < 10

q3 <- plotUMAP(sce, colour_by = "Keep at 10% cutoff") +
  scale_fill_manual(name = "Keep", values = c("#999999", "#E69F00")) +
  theme(text = element_text(size = text_size),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))


# Based on supplemental figure 2, run mbkmeans with 6 clusters
set.seed(621)
six <-     mbkmeans(sce, clusters = 6,
                    batch_size = 500,
                    reduceMethod = "PCA",
                    calc_wcss = TRUE
)
sce$Clusters <- as.factor(six$Clusters)
q4 <- plotUMAP(sce, colour_by = "Clusters") +
  theme(text = element_text(size = text_size),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

everything <- colData(sce)
A <- table(everything$Clusters)

B <- subset(everything, everything$Keep == T)
B <- table(B$Clusters)

C <- subset(everything, everything$`Keep at 10% cutoff` == T)
C <- table(C$Clusters)

df <- as.data.frame(cbind(A, B, C))
df$cluster <- rownames(df)
df <- melt(df)
df$ind <- paste(df$cluster, df$variable, sep = "")

q5tmp <- ggplot(data = df, aes(x = cluster, y = value, fill = ind)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cluster", y = "Number of cells") + theme_bw() +
  theme(legend.position = "none", text = element_text(size = text_size)) +
  ggtitle("Cluster membership by filtering method") +
  scale_fill_manual(values = c("#7396bf", "#96b1cf", "#b9cbdf", "#e6914d",
                               "#ecad79", "#f2c8a6", "#7dc072", "#9dd095",
                               "#bedfb9", "#ce6964", "#da8f8b", "#edc7c5",
                               "#9b81b1", "#b4a0c5", "#cdc0d8", "#b08a82",
                               "#c3a6a2", "#d7c3c1"))

# Construct the six grobs - three symbols and three labels
L1 <- rectGrob(height = .6, width = .6, gp = gpar(fill = "#7396bf", col = NA))
L2 <- rectGrob(height = .6, width = .6, gp = gpar(fill = "#96b1cf", col = NA))
L3 <- rectGrob(height = .6, width = .6, gp = gpar(fill = "#b9cbdf", col = NA))
T1 <- text_grob("All cells", x = .1, just = "left", size = 16)
T2 <- text_grob("miQC filtering", x = .1, just = "left", size = 16)
T3 <- text_grob("10% cutoff", x = .1, just = "left", size = 16)

# Construct a gtable - 2 columns X 3 rows
leg <- gtable(width = unit(c(1, 3), "cm"), height = unit(c(1, 1, 1), "cm"))

# Place the six groba into the table
leg <- gtable_add_grob(leg, L1, t = 1, l = 1)
leg <- gtable_add_grob(leg, L2, t = 2, l = 1)
leg <- gtable_add_grob(leg, L3, t = 3, l = 1)
leg <- gtable_add_grob(leg, T1, t = 1, l = 2)
leg <- gtable_add_grob(leg, T2, t = 2, l = 2)
leg <- gtable_add_grob(leg, T3, t = 3, l = 2)

g <- ggplotGrob(q5tmp)

# Get the position of the panel,
# add a column to the right of the panel,
# put the legend into that column,
# and then add another spacing column
pos <- g$layout[grepl("panel", g$layout$name), c("t", "l")]
g <- gtable_add_cols(g, sum(leg$widths), pos$l)
g <- gtable_add_grob(g, leg, t = pos$t, l = pos$l + 1)
q5 <- gtable_add_cols(g, unit(20, "pt"), pos$l)
                     
png("~/Documents/Figures/SmartQC Paper/Figure 5.png",
    width = 7400, height = 3200, res = 300)
q1 + q4 + q5 + q2 + q3 + plot_spacer()
dev.off()