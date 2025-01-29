## @author: Marta E. Camarena
## @date: 21st Jan 2025

##---------------- Load libraries ----------------##
library(ChIPseeker)
library(clusterProfiler)

library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)

wd="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/ATACSeq_GSE120508"
indir="/datasets/marta/scTestis_GSE120508/GSE120507_ATAC"

file = file.path(indir,"GSE120507_Combined.narrowPeak")

########### - Process the narrowPeak File - ################
peak <- readPeakFile(file)
covplot(peak, weightCol="V5")

## X chromosome
covplot(peak, weightCol="V5", chrs=c("chrX"), xlim=c(4.5e7, 5e7))

########### - Annotate the File - ################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_annotation <- annotatePeak(
  peak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),  # Define upstream/downstream of TSS
  level = "gene",  # Annotate at the gene level
  annoDb = "org.Hs.eg.db"  # For mapping to gene symbols
)
