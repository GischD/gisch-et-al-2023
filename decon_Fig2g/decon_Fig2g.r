if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("granulator")
library(granulator)
library(ggplot2)

bulkCutRun <- as.matrix(read.table(file="/Users/michelle/Desktop/Rauchman\ Paper/celltype.decon.test/CutRunData.Top10Percent.txt", header=TRUE, row.names=1))
colnames(bulkCutRun) <- c("DNAme", "H3K4me3", "H3K4me1", "H3K27Ac", "H3K27me3")

bulkCutRun.act <- bulkCutRun[,-c(1,5)]

scRNA <- as.matrix(read.table(file="/Users/michelle/Desktop/Rauchman\ Paper/celltype.decon.test/KidneyAtlasRNA.Top10Percent.txt", header=TRUE, row.names=1))
scRNA <- scRNA[,-10]
scRNA <- scRNA[,-4]
scRNA <- scRNA[,-c(9, 8, 10)]
scRNA <- scRNA[,-5]
scRNA <- scRNA[,-1] # remove DTL because too similar to TAL
scRNA <- scRNA[,-3] # remove EC - not on original request


plot_similarity(sigMatrix=scRNA)
#k=4

decon <- deconvolute(m = bulkCutRun.act, sigMatrix = scRNA)
prop <- decon$proportions$dtangle_sig1[1:3, 1:5]

plot_proportions(deconvoluted = decon, method = 'dtangle', signature = 'sig1') + FontSize(14)

plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = TRUE)
