library(Signac)
library(Seurat)
library(SeuratDisk)
library(stringr)
library(genomation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(future)
library(parallel)
detectCores()
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 32 * 1024 ^ 3)


#Call Peaks
setwd("/Users/michelle/Desktop/Rauchman\ Paper/KPMP_Multiomic")
#load full object
b.atac <- LoadH5Seurat('Kidney_KPMP_Multiome_ATAC-RNA_Seurat_12152021-002.h5Seurat')

#update fragment paths
sample.folder <- list.files(pattern='10X_Dual') 

frags <- Fragments(b.atac)  # get list of fragment objects
Fragments(b.atac) <- NULL  # remove fragment information from assay

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = paste0(sample.folder[[i]],'/atac_fragments.tsv.gz')) # update path
}
Fragments(b.atac) <- frags # assign updated list back to the object
# Traceback
Fragments(b.atac)

#Call peaks in PT-S1, PT-S2, and aPT
peaks.PT12a <- CallPeaks(
  object = b.atac,
  group.by = "subclass.l3",
  idents = c("PT-S1", "PT-S2", "aPT"),
  macs2.path = "/Users/michellepherson/opt/anaconda3/bin/macs2"
)
peaks.PT12a.data <- as.data.frame(peaks.PT12a)

peaks.PT12a.bed <- peaks.PT12a.data[,c(1,2,3)]
write.table(peaks.PT12a.bed, file="peaks.PT12a.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#Call peaks in C-TAL and aTAL1/2
peaks.cTAL.aTAL12 <- CallPeaks(
  object = b.atac,
  group.by = "subclass.l3",
  idents = c("C-TAL", "aTAL1", "aTAL2"),
  macs2.path = "/Users/michelle/opt/anaconda3/envs/macs2/bin/macs2"
)
peaks.cTAL.data <- as.data.frame(peaks.cTAL.aTAL12)

peaks.cTAL.bed <- peaks.cTAL.data[,c(1,2,3)]
write.table(peaks.cTAL.bed, file="peaks.cTAL.aTAL12.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(peaks.cTAL.data, file="peaks.cTAL.aTAL12.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


######PT Object
#load RNA assay only
b.rna <- LoadH5Seurat('Kidney_KPMP_Multiome_ATAC-RNA_Seurat_12152021-002.h5Seurat', assay="RNA")

#subset RNA object for only cell types of interest - PT
b.rna.PT <- subset(b.rna, idents=c("PT-S1", "PT-S2", "aPT"))
colnames(b.rna.PT)

#load called peaks from ATAC data - removed nonstandard chromosome information
peaks.PT12a <- readBed(file="/Volumes/Mobius/Rauchman/Rauchman.KPMP.final.files/KPMP_Multiomic/Analysis/PT12a/peaks.PT12a.bed")

#extract counts for PT specific peaks in PT cells (S1, S2, aPT) only - takes awhile to run
macs2_counts.PT12a <- FeatureMatrix(
  fragments = Fragments(b.atac),
  features = peaks.PT12a,
  cells = colnames(b.rna.PT)
)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#add new PT specific peak assay to RNA PT object
b.rna.PT[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.PT12a,
  fragments = Fragments(b.atac),
  annotation = annotation
)

#normalize PT specific peak counts
DefaultAssay(b.rna.PT) <- "peaks"
b.rna.PT <- FindTopFeatures(b.rna.PT, min.cutoff = 5)
b.rna.PT <- RunTFIDF(b.rna.PT)
b.rna.PT <- RunSVD(b.rna.PT)

DefaultAssay(b.rna.PT) <- "peaks"

#compute the GC content for each peak
b.rna.PT <- RegionStats(b.rna.PT, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to all genes - takes a long time to run 
b.rna.PT <- LinkPeaks(
  object = b.rna.PT,
  peak.assay = "peaks",
  expression.assay = "RNA",
)

linkages <- as.data.frame(b.rna.PT@assays[["peaks"]]@links)
write.table(linkages, file="linkages.RNA.PI.PT12aPeaks.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#saveRDS(b.rna.PT, file="multiome.RNA.PT.PT12aPeaks.linkages.rds")
#b.rna.PT <- readRDS(file="/Users/michelle/Desktop/Rauchman\ Paper/KPMP_Multiomic/multiome.RNA.PT.PT12aPeaks.linkages.rds")

DimPlot(b.rna.PT, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()



######TAL Object
#subset RNA object for only cell types of interest - TAL
b.rna.TAL <- subset(b.rna, idents=c("C-TAL", "aTAL1", "aTAL2"))
colnames(b.rna.TAL)

#load called peaks from ATAC data  removed nonstandard chromosome information
peaks.cTAL.aTAL12 <- readBed(file="/Users/michelle/Desktop/Rauchman\ Paper/KPMP_Multiomic/TAL/peaks.cTAL.aTAL12.bed")

#extract counts for TAL specific peaks in TAL cells (c-TAL, aTAL1, aTAL2) only - takes awhile to run
macs2_counts.TAL <- FeatureMatrix(
  fragments = Fragments(b.atac),
  features = peaks.cTAL.aTAL12,
  cells = colnames(b.rna.TAL)
)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#add new TAL specific peak assay to RNA TAL object
b.rna.TAL[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.TAL,
  fragments = Fragments(b.atac),
  annotation = annotation
)

#normalize TAL specific peak counts
DefaultAssay(b.rna.TAL) <- "peaks"
b.rna.TAL <- FindTopFeatures(b.rna.TAL, min.cutoff = 5)
b.rna.TAL <- RunTFIDF(b.rna.TAL)
b.rna.TAL <- RunSVD(b.rna.TAL)

DefaultAssay(b.rna.TAL) <- "peaks"

#compute the GC content for each peak
b.rna.TAL <- RegionStats(b.rna.TAL, genome = BSgenome.Hsapiens.UCSC.hg38)

saveRDS(b.rna.TAL, file="multiome.RNA.TAL.Peaks.rds")
