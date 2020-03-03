
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tximport)
  library(scater)
  library(scran)
  library(ggplot2)
  library(flexmix)
  library(splines)
})

load_data <- function(tumorname){
  filename <- paste("../../data/tumors/", tumorname, "/alevin_output/alevin/quants_mat.gz",sep="")
  txia <- suppressMessages(tximport(filename, type = "alevin"))
  sce <- SingleCellExperiment(assays = list(counts = txia$counts))
  
  sce
}

x2 <- load_data("16030X2")


gene_map <- fread("../../data/index/gene2symbol.tsv", header = F)
setnames(gene_map, c("ensembl_id", "gene_symbol"))

name_swap <- function(sce){
  rowData(sce)$ensembl_id <- rownames(sce)
  rowData(sce) <- left_join(as.data.frame(rowData(sce)), gene_map)
  rownames(sce) <- rowData(sce)$gene_symbol
  
  rowData(sce)$gene_symbol <- NULL
  rowData(sce)$ensembl_id <- NULL
  
  sce
}

x2 <- name_swap(x2)


mt_genes <- grepl("^MT-",  rownames(x2))

feature_ctrls <- list(mito = rownames(x2)[mt_genes])


run_scater <- function(sce){
  sce <- addPerCellQC(sce, subsets = feature_ctrls,
                      BPPARAM = BiocParallel::MulticoreParam())
  
  sce <- addPerFeatureQC(sce, BPPARAM = BiocParallel::MulticoreParam())
  rowData(sce)$Gene <- rownames(rowData(sce))
  sce
}
x2 <- run_scater(x2)




remove_worst_cells <- function(sce){
  filt <- sce$sum > 500 & sce$detected > 100
  length(which(!filt))
  sce<-sce[,filt]
  sce
}

x2 <- remove_worst_cells(x2)


metrics <- as.data.frame(colData(x2))

p <- ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point()
p + stat_smooth(method="loess",formula=y~x,size=1,se=F,colour="blue")

options(scipen = 5)
set.seed(1010)
model<-flexmix(subsets_mito_percent~bs(detected,Boundary.knots = c(0,6000)),data=metrics,k=2)
model@components

fitted_models <- as.data.frame(cbind(metrics$detected,fitted(model)))

ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point() + 
  geom_line(data=fitted_models, aes(x=V1, y=Comp.1, color="1")) +
  geom_line(data=fitted_models, aes(x=V1, y=Comp.2, color="2"))+ ggtitle("16030X2 spline")


post <- posterior(model)
metrics$posterior_dist1 <- post[,1]
metrics$posterior_dist2 <- post[,2]


ggplot(metrics, aes(x=detected,y=subsets_mito_percent,colour=posterior_dist2)) + geom_point() + ggtitle("16030X2 spline")
# + geom_line(data=fitted_models, aes(x=V1, y=Comp.1)) +
#  geom_line(data=fitted_models, aes(x=V1, y=Comp.2))


metrics$keep<-metrics$posterior_dist2>=0.25
ggplot(metrics, aes(x=detected,y=subsets_mito_percent,colour=keep)) + geom_point() + ggtitle("16030X2 spline")
table(metrics$keep)

x2 <- x2[,metrics$keep]
dim(x2)