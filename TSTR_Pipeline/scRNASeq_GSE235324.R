## @author: Marta E. Camarena
## @date: 21st Jan 2025

##---------------- Load libraries ----------------##
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
# library(patchwork)

wd="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/scMultiomicsTestis_GSE235324"


##-----------0. Importing Seurat Object -----------##
## Created by data2SeuratObj.R
seurat.obj = readRDS(file = file.path(wd,"NormalSeuratObject.rds"))
## how many cells per cell type and patient?
table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$Sample)
seurat.obj$CellType <- factor(seurat.obj$CellType, levels = c("Undiff.SPG-1", "Undiff.SPG-2","Diff.ing SPG", 
                                                              "preL", 
                                                              "L1/2", "L3", 
                                                              "Z", 
                                                              "preP", 
                                                              "P", 
                                                              "D", 
                                                              "SPC7", 
                                                              "S1", "S2", "S3", 
                                                              "ST", "T", "LD"))

##--------------1. Quality Control ---------------##
dir.create(file.path(wd,"plots/01_QC"), recursive = T)
## - Filter out cells with the highest and lowest number of 1% 
## - Filter out cells with more than 10% of mitochondrial genes expression.

seurat.obj[["percent.mt"]]  <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(wd,"plots/01_QC/Vln_feat_counts_mt.png"), width=7.20, height=8.03)

VlnPlot(seurat.obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, group.by = "CellType", pt.size = 0) +
  geom_hline(yintercept=5000)
ggsave(file.path(wd,"plots/01_QC/Vln_feat_counts_mt.celltype.png"), width=12.46, height=9.5)

VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "CellType", pt.size = 0) + geom_hline(yintercept=6000)
ggsave(file.path(wd,"plots/01_QC/Vln_counts_threshold6000.celltype.png"), width=7.20, height=8.03)

VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample", pt.size = 0)
ggsave(file.path(wd,"plots/01_QC/Vln_feat_counts_mt.sample.png"), width=7.20, height=8.03)

VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "Sample", pt.size = 0) + geom_hline(yintercept=6000)
ggsave(file.path(wd,"plots/01_QC/Vln_counts_threshold6000.samples.png"), width=7.20, height=8.03)

## subset based on literature thresholds
# seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 200 & nCount_RNA < 6000 & percent.mt < 10)

##-----------3. Analysis -----------##
## Normalization - already normalized
# seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.obj), 10)
LabelPoints(plot = VariableFeaturePlot(seurat.obj), points = top10, repel = TRUE)
ggsave(file.path(wd,"plots/01_QC/VariableFeaturePlot.png"), width=7.20, height=8.03)

## prior to PCA
all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)
## PCA
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))

## Clustering
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:15)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)

## UMAP
dir.create(file.path(wd,"plots/02_UMAPs"))
seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)
DimPlot(seurat.obj, reduction = "umap") + ggtitle("UMAP per Identity")
ggsave(file.path(wd,"plots/02_UMAPs/umaps_identity.png"), width=7.20, height=8.03)

DimPlot(seurat.obj, reduction = "umap", group.by = "CellType") + ggtitle("UMAP per Cell Type")
ggsave(file.path(wd,"plots/02_UMAPs/umaps_celltypes.png"), width=7.20, height=8.03)

DimPlot(seurat.obj, reduction = "umap", group.by = "CellType", split.by = "Sample") + ggtitle("UMAP per Cell Type and Patient")
ggsave(file.path(wd,"plots/02_UMAPs/umaps_celltypes_patient.png"), width=15.15, height=8.03)

FeaturePlot(seurat.obj, features = c("TNP1","PRM2"), label=T) # sperm


##-------------6. Tumor-specific genes--------------##
seurat.obj <- subset(seurat.obj, subset = CellType != "NA")

candidates = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes/TSTR_candidatesORFs_fullcharacterized.csv")
candidates$gene_name = gsub("ENSG00000287861","LOC107984132", candidates$gene_name)
CTx = candidates %>% subset(coding_noncoding_chr == "CT-X") %>% select(gene_name) %>% unique() %>% pull(gene_name)
CTnonx = candidates %>% subset(coding_noncoding_chr == "CT-nonX") %>% select(gene_name) %>% unique() %>% pull(gene_name)
NCnonX = candidates %>% subset(coding_noncoding_chr == "Noncoding-nonX") %>% select(gene_name) %>% unique() %>% pull(gene_name)

all_TSTR = candidates %>% subset(coding_noncoding_chr != "altORFs") %>% select(gene_name) %>% unique() %>% pull(gene_name)

dir.create(file.path(wd,"plots/05_TSA"))

## CT-X
DotPlot(seurat.obj, features = CTx, group.by = "CellType") + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="CT-X TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/CTx.png"))

## CT-nonX
DotPlot(seurat.obj, features = CTnonx, group.by = "CellType") + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="CT-nonX TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
  
ggsave(file.path(wd,"plots/05_TSA/CTnonX.png"))

## Noncoding-X
DotPlot(seurat.obj, features = NCnonX, group.by = "CellType") + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="Noncoding-nonX TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/Noncoding.png"))


## markers
markers = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/celltype_markers.csv")
DotPlot(seurat.obj, features = markers$gene_name, group.by = "CellType") + 
  RotatedAxis() + 
  labs(x="Marker",
       y="Cell Type",
       title="Stage Markers") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/markers.png"))

DotPlot(seurat.obj, features =  c("TNP1","PRM2","SHBG", "SRD5A2"), group.by = "CellType") + 
  RotatedAxis() + 
  labs(x="Marker",
       y="Cell Type",
       title="Stage Markers") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")

# FeaturePlot(seurat.obj, reduction = "umap", features = NCnonX) 

### 
markers = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/celltype_markers.csv")
DoHeatmap(
  seurat.obj,
  features = markers$gene_name,
  group.by = "CellType"
) +
  scale_fill_gradientn(colors = c("#4682B4", "#FFFFE0", "#CA3433")) +
  NoLegend()

DoHeatmap(
  seurat.obj,
  features = all_TSTR,
  group.by = "CellType"
) +
  scale_fill_gradientn(colors = c("#4682B4", "#FFFFE0", "#CA3433")) +
  NoLegend()

