library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)
library(tidyverse)
library(data.table)
library(gplots)
library(RColorBrewer)
library(irlba)
library(Rtsne)
library(gridExtra)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase","limma","genefilter","edge","qvalue"))



con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
ls()

class(bottomly.eset)
save(bottomly.eset, file="bottomly.Rdata")
pdata=pData(bottomly.eset)
dim(pdata)

edata=as.matrix(exprs(bottomly.eset))
dim(edata)

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]
png(("bottomly_heatmap_clustered_column_scale.png"),  height = 700, width = 700)
my_palette <- colorRampPalette(c("blue", "white", "orange"))(n=299)
heatmap.2(edata,
          main = "Bottomly et al. Clustered",
          notecol = "black",
          density.info = "none",
          margins = c(12,9), 
          col = my_palette,
          dendrogram ="none",
          scale = "column",
          Colv = FALSE)

edata <- t(scale(t(edata), scale=FALSE, center=TRUE))
svd.out <- svd(edata)
names(svd.out)
## [1] "d" "u" "v"
print(paste("Dimension of left singular vectors:", dim(svd.out$u)))
## [1] "Dimension of left singular vectors: 1049"
## [2] "Dimension of left singular vectors: 21"
print(paste("Length of singular values:",length(svd.out$d)))
## [1] "Length of singular values: 21"
print(paste("Dimension of right singular vectors:",dim(svd.out$v)))
## [1] "Dimension of right singular vectors: 21"
## [2] "Dimension of right singular vectors: 21"
plot(1:ncol(edata), svd.out$v[,1],pch=20)

PC = data.table(svd.out$v,pData(bottomly.eset))

ggplot(PC) + geom_point(aes(x=V2, y=V6, col=as.factor(strain)))

ggplot(PC) + geom_boxplot(aes(x=as.factor(strain), y=V1))
png(("Mostowski-problem3.png"),  height = 700, width = 700)
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V1))
dev.off()
PCu = data.table(svd.out$u,pData(bottomly.eset)) ## hom 3

png(("plot.png"),  height = 700, width = 700)
ggplot(PCu) + geom_point(aes(x=V1, y=V2, col=as.factor(strain)))

plot1 = ggplot(PCu) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) #hom4
plot2 = ggplot(PCu) + geom_violin(aes(x=as.factor(strain), y=V2),draw_quantiles = c(0.25, 0.5, 0.75)) #hom4
plot3 = ggplot(PCu) + geom_violin(aes(x=as.factor(strain), y=V3),draw_quantiles = c(0.25, 0.5, 0.75)) #hom4
plot4 = ggplot(PCu) + geom_violin(aes(x=as.factor(strain), y=V4),draw_quantiles = c(0.25, 0.5, 0.75)) #hom4
plot5 = ggplot(PCu) + geom_violin(aes(x=as.factor(strain), y=V5),draw_quantiles = c(0.25, 0.5, 0.75)) #hom4
grid.arrange(plot1, plot2, plot3, plot4, plot5,  ncol=5)
dev.off()
# Set a seed for reproducible results
set.seed(1)
# complexity is a hyperparameter needed for this algorithm. 30 is a default
#tsne_out <- Rtsne(edata,pca=TRUE,perplexity=30)
#tsne_out = data.table(tsne_out$Y)
#ggplot(tsne_out) + geom_point(aes(x=V1, y=V2))

#edata=as.matrix(exprs(bottomly.eset))
#dim(edata)

#ggplot(edata, col = cl$cluster)
#dim(cl)
tsne_out <- Rtsne(edata,pca=TRUE,perplexity=30)
#tsne_out = data.table(tsne_out$Y)
tsne_out = as.data.frame(tsne_out$Y)
fit_cluster_kmeans = kmeans(scale(tsne_out), 5)
tsne_out$cl_kmeans = factor(fit_cluster_kmeans$cluster)
png(("test.png"),  height = 700, width = 700)
ggplot(tsne_out) + geom_point(aes(x=V1, y=V2, colour =  cl_kmeans))
dev.off()





## hom 4





