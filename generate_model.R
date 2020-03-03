suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(flexmix)
  library(splines)
})

source("config.R")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="16030X2"),
  make_option(c("-m", "--model"), type="character", default="linear")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# load data
location <- paste(path_to_mito_filtering,
                  "/sce_objects/pre_filtering/", opt$file, ".rds", sep="")
sce <- readRDS(location)
metrics <- as.data.frame(colData(sce))

# visualize data distribution
png(paste(path_to_mito_filtering,"/plots/",opt$file,"_loess.png",sep=""),width=1200, height=800)
ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point() +
  stat_smooth(method="loess",formula=y~x,size=1,se=F,colour="blue") +
  labs(x="Unique genes found", y="Percent reads mitochondrial") +
  ggtitle(opt$file)
dev.off()

# run flexmix to get mixture model
if(opt$model=="linear"){
  model<-flexmix(subsets_mito_percent~detected, data=metrics, k=2)
} else if(opt$model=="spline"){
  model<-flexmix(subsets_mito_percent~bs(detected), data=metrics, k=2)
} else if(opt$model=="polynomial"){
  model<-flexmix(subsets_mito_percent~poly(detected, degree=2), data=metrics, k=2)
}

# plot model parameters on data
fitted_models <- as.data.frame(cbind(metrics$detected,fitted(model)))

png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_fit.png", sep=""), width=1200, height=800)
ggplot(metrics, aes(x=detected,y=subsets_mito_percent)) + geom_point() +
  geom_line(data=fitted_models, aes(x=V1, y=Comp.1)) +
  geom_line(data=fitted_models, aes(x=V1, y=Comp.2)) +
  labs(x="Unique genes found", y="Percent reads mitochondrial") +
  ggtitle(opt$file)
dev.off()

# determine which model component is the dead cell distribution
intercept1 <- parameters(model, component=1)[1]
intercept2 <- parameters(model, component=2)[1]
if(intercept1 > intercept2){
  dead_cell_dist <- 1
} else {
  dead_cell_dist <- 2
}

# calculate posterior likelihood for each cell
post <- posterior(model)
metrics$probability_dead <- post[,dead_cell_dist]

#TODO: decide if it's important to display the lines in this plot
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_posterior.png", sep=""), width=1200, height=800)
ggplot(metrics, aes(x=detected, y=subsets_mito_percent, colour=probability_dead)) +
  labs(x="Unique genes found", y="Percent reads mitochondrial") +
  geom_point() + ggtitle(opt$file) 
dev.off()

# determine which cells to keep
png(paste(path_to_mito_filtering, "/plots/", opt$file, "_",
          opt$model, "_keep.png", sep=""), width=1200, height=800)
metrics$keep<-metrics$probability_dead <= 0.75
ggplot(metrics, aes(x=detected,y=subsets_mito_percent,colour=keep)) +
  labs(x="Unique genes found", y="Percent reads mitochondrial") +
  geom_point() + ggtitle(opt$file) 
table(metrics$keep)
dev.off()

colData(sce)$probability_dead <- metrics$probability_dead
colData(sce)$keep <- metrics$keep

# save object with filtering data
filename <- paste(opt$file, "_", opt$model, ".rds", sep="")
filepath <- paste(path_to_mito_filtering, "sce_objects",
                  "post_filtering", filename, sep="/")
saveRDS(sce, file=filepath)