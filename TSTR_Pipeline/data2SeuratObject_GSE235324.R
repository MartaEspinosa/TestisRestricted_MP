# Author: Marta E. Camarena
# email: mespinosa@researchmar.net
# Date: 21st Sept 2025
# GSE235324

library(dplyr)
library(Seurat)
library(tidyr)

outdir="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/scMultiomicsTestis_GSE235324"

## NOA1
NOA1 = read.table("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_exp.log2TPM.NOA1.tsv")
metadata_NOA1 = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_Info.NOA1.tsv", header=T, sep="\t")
rownames(metadata_NOA1) = metadata_NOA1$Cell
metadata_NOA1$Cell = NULL
print(unique(metadata_NOA1$Stage))
names(metadata_NOA1)[2] = "Sample"
names(metadata_NOA1)[8] = "CellType"
## Seurat
NOA1_Seurat <- CreateSeuratObject(counts = NOA1, meta.data = metadata_NOA1)
saveRDS(NOA1_Seurat, file = file.path(outdir,"NOA1SeuratObject.rds"))

## NOA2
NOA2 = read.table("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_exp.log2TPM.NOA2.tsv")
metadata_NOA2 = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_Info.NOA2.tsv", header=T, sep="\t")
rownames(metadata_NOA2) = metadata_NOA2$Cell
metadata_NOA2$Cell = NULL
print(unique(metadata_NOA2$Stage))
names(metadata_NOA2)[2] = "Sample"
names(metadata_NOA2)[8] = "CellType"
## Seurat
NOA2_Seurat <- CreateSeuratObject(counts = NOA2, meta.data = metadata_NOA2)
saveRDS(NOA2_Seurat, file = file.path(outdir,"NOA2SeuratObject.rds"))

## Normal
Normal = read.table("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_exp.log2TPM.Normal.tsv")
metadata_Normal = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/GSE235321_Info.Normal.tsv", header=T, sep="\t")
rownames(metadata_Normal) = metadata_Normal$Cell
metadata_Normal$Cell = NULL
print(unique(metadata_Normal$Stage))
names(metadata_Normal)[2] = "Sample"
names(metadata_Normal)[8] = "CellType"
## Seurat
Normal_Seurat <- CreateSeuratObject(counts = Normal, meta.data = metadata_Normal)
saveRDS(Normal_Seurat, file = file.path(outdir,"NormalSeuratObject.rds"))
