library(ggplot2)

#bedtools fisher -a a.bed -b b.bed -g t.genome
data <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Fisher.test/data.txt", header=TRUE)

data$Sample1 <- as.character(data$Sample1)
#Then turn it back into a factor with the levels in the correct order
data$Sample1 <- factor(data$Sample1, levels=unique(data$Sample1))

data$Sample2 <- as.character(data$Sample2)
#Then turn it back into a factor with the levels in the correct order
data$Sample2 <- factor(data$Sample2, levels=unique(data$Sample2))

level_order <- c("H3K27me3", "H3K27Ac", "H3K4me1", "H3K4me3", "Hyper.DNAme.TI", "Hypo.DNAme.TI", "Hyper.DNAme.Glom", "Hypo.DNAme.Glom", "Bulk.ATAC")

ggplot(data = data,aes(x = as.factor(Sample1), y = as.factor(Sample2))) + 
  geom_point(aes(size=as.numeric(Nlog10pval), color=as.numeric(log10ratio))) +
  xlab('') + ylab('') + scale_size_continuous(name="-log10 p-value") + scale_color_continuous(name="log10 odds Ratio") +
  scale_color_gradient2(low="darkblue", high="darkred", mid="gray", midpoint=0, name="log10 odds ratio") +
  theme_classic() +
  scale_y_discrete(limits = level_order) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) +
  theme(axis.text.y = element_text(size=12))
