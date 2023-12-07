#install.packages("UpSetR")
library("UpSetR")

dips <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Annotation/bulk.ATAC.peaks/bulk.ATAC.Map/upset.plot/dips.txt")
H3K4me3 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Annotation/bulk.ATAC.peaks/bulk.ATAC.Map/upset.plot/H3K4me3.txt")
H3K4me1 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Annotation/bulk.ATAC.peaks/bulk.ATAC.Map/upset.plot/H3K4me1.txt")
H3K27ac <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Annotation/bulk.ATAC.peaks/bulk.ATAC.Map/upset.plot/H3K27ac.txt")
H3K27me3 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/Annotation/bulk.ATAC.peaks/bulk.ATAC.Map/upset.plot/H3K27me3.txt")


listInput <- list(TI.DNAm.dip = dips[,1], H3K4me3 = H3K4me3[,1], H3K4me1 = H3K4me1[,1], H3K27ac = H3K27ac[,1], H3K27me3 = H3K27me3[,1])
upset(fromList(listInput), group.by = "sets", keep.order=TRUE)


upset(fromList(listInput), order.by = "freq", query.legend = "bottom", text.scale = 2, mainbar.y.label = "Peak Intersections", sets.x.label = "Total ATAC Peaks",
      queries = list(list(query = intersects, params = list("TI.DNAm.dip", "H3K4me3", "H3K27ac"), color = "darkgreen", active = T, query.name = "Active Promoter"), 
                                          list(query = intersects, params = list("TI.DNAm.dip", "H3K4me3", "H3K27ac", "H3K4me1"), color = "darkgreen", active = T), 
                                          list(query = intersects, params = list("TI.DNAm.dip", "H3K27ac", "H3K4me1"), color = "blue", active = T, query.name = "Predicted Enhancer"),
                                          list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3"), color = "darkred", active = T, query.name = "Repressed Region"),
                                          list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3", "H3K4me3"), color = "darkred", active = T),
                                          list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3", "H3K4me3", "H3K4me1"), color = "darkred", active = T),
                                          list(query = intersects, params = list("H3K27me3"), color = "darkred", active = T)))


dips <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/ATAC.Map/PT12a.ATAC.DNAme.bed")
H3K4me3 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/ATAC.Map/PT12a.ATAC.H3K4me3.bed")
H3K4me1 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/ATAC.Map/PT12a.ATAC.H3K4me1.bed")
H3K27me3 <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/ATAC.Map/PT12a.ATAC.H3K27me3.bed")
H3K27ac <- read.table(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/ATAC.Map/PT12a.ATAC.H3K27ac.bed")

listInput <- list(TI.DNAm.dip = dips[,4], H3K4me3 = H3K4me3[,4], H3K4me1 = H3K4me1[,4], H3K27ac = H3K27ac[,4], H3K27me3 = H3K27me3[,4])

upset(fromList(listInput), order.by = "freq", query.legend = "bottom", text.scale = 2, mainbar.y.label = "Peak Intersections", sets.x.label = "Total ATAC Peaks",
      queries = list(list(query = intersects, params = list("TI.DNAm.dip", "H3K4me3", "H3K27ac"), color = "darkgreen", active = T, query.name = "Active Promoter"), 
                     list(query = intersects, params = list("TI.DNAm.dip", "H3K4me3", "H3K27ac", "H3K4me1"), color = "darkgreen", active = T), 
                     list(query = intersects, params = list("TI.DNAm.dip", "H3K27ac", "H3K4me1"), color = "blue", active = T, query.name = "Predicted Enhancer"),
                     list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3"), color = "darkred", active = T, query.name = "Repressed Region"),
                     list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3", "H3K4me3"), color = "darkred", active = T),
                     list(query = intersects, params = list("TI.DNAm.dip", "H3K27me3", "H3K4me3", "H3K4me1"), color = "darkred", active = T),
                     list(query = intersects, params = list("H3K27me3"), color = "darkred", active = T)))
