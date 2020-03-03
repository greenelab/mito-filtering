suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tximport)
  library(scater)
  library(scran)
  library(ggplot2)
  library(flexmix)
  library(scRNAseq)
  library(biomaRt)
})

#load data
dataname<-"MessmerESCData"
sce <- MessmerESCData()

#name swap
gene_map <- fread("../../data/index/gene2symbol.tsv", header = F)
setnames(gene_map, c("ensembl_id", "gene_symbol"))
gene_map$ensembl_id <- gsub(pattern = "\\.[0-9]+", replacement = "", gene_map$ensembl_id)

rowData(sce)$ensembl_id <- rownames(sce)
sce <- sce[rownames(sce) %in% gene_map$ensembl_id,]
rowData(sce) <- inner_join(as.data.frame(rowData(sce)), gene_map)
rownames(sce) <- rowData(sce)$gene_symbol

#get mito genes
mt_genes <- grepl("^MT-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sum(mt_genes)

#run scater
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

sce <- addPerFeatureQC(sce, BPPARAM = BiocParallel::MulticoreParam())
rowData(sce)$Gene <- rownames(rowData(sce))

#remove worst cells
filt <- sce$sum > 500 & sce$detected > 100
length(which(!filt))
sce<-sce[,filt]


metrics <- as.data.frame(colData(sce))

p <- ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point()
p + stat_smooth(method="loess",formula=y~x,size=1,se=F,colour="blue")  +
  ggtitle(dataname)

p <- ggplot(metrics, aes(x=total_counts,y=subsets_mito_percent)) + geom_point()
p <- ggplot(metrics, aes(x=detected,y=log10_total_counts)) + geom_point()

options(scipen = 5)
set.seed(1010)
model<-flexmix(subsets_mito_percent~detected,data=metrics,k=2)
model@components

slope_1=model@components$Comp.1[[1]]@parameters$coef[2]
intercept_1=model@components$Comp.1[[1]]@parameters$coef[1]

slope_2=model@components$Comp.2[[1]]@parameters$coef[2]
intercept_2=model@components$Comp.2[[1]]@parameters$coef[1]

ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point() +
  geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) +
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2)) +
  ggtitle(dataname)

post <- posterior(model)
metrics$posterior_dist1 <- post[,1]
metrics$posterior_dist2 <- post[,2]

ggplot(metrics, aes(x=detected,y=subsets_mito_percent,colour=posterior_dist1)) +
  geom_point() + geom_abline(data=metrics,aes(slope=slope_1,intercept=intercept_1)) +
  geom_abline(data=metrics,aes(slope=slope_2,intercept=intercept_2)) +
  ggtitle(dataname)

metrics$keep<-metrics$posterior_dist1<=0.75
ggplot(metrics, aes(x=detected,y=subsets_mito_percent,colour=keep)) + geom_point()
table(metrics$keep)
