#regional mRNA, DNAme, and CutRun Correlation PLUS Atlas RNA
#Top 500 most DE genes from TI vs Glom
data <- read.table(file="gene.data.top500DE.txt", header=TRUE)

data.matrix <- data[,c(6,5,7,8,2,4,3)]
row.names(data.matrix) <- data[,1]
colnames(data.matrix) <- c("H3K4me3 Pro", "H3K4me1 Pro", "H3K27Ac Pro", "H3K27me3 Pro", "TI regional mRNA", "PTS1/S2 snRNA", "C-TAL snRNA")

corrplot(cor(data.matrix, method="spearman"),
         method = "square", type="lower",      
         addrect = 2,              # If order = "hclust", number of cluster rectangles
         rect.col = 3,             # Color of the rectangles
         rect.lwd = 3,             # Line width of the rectangles
         col=colorRampPalette(c("blue4", "white", "red4"))(100),
         tl.cex = 1.7, 
         tl.col = "black") 
