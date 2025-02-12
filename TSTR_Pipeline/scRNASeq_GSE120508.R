## @author: Marta E. Camarena
## @date: 21st Jan 2025

##---------------- Load libraries ----------------##
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
# library(patchwork)

wd="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/scRNASeq_GSE120508"

## create seurat Object
# counts = read.csv("/datasets/marta/scTestis_GSE120508/GSE112013_scRNA/GSE112013_Combined_UMI_table.txt", sep="\t")
# rownames(counts) = counts$Gene
# counts$Gene = NULL
# seurat.obj = CreateSeuratObject(counts = counts, project = "GSE120508", min.cells = 3, min.features = 200)
# saveRDS(seurat.obj, file = file.path(wd,"SeuratObject.rds"))

#-----------0. Importing Seurat Object -----------##
## Created by data2SeuratObj.R
seurat.obj = readRDS(file = file.path(wd,"SeuratObject.rds"))

##--------------1. Quality Control ---------------##
dir.create(file.path(wd,"plots/01_QC"), recursive = T)
## - Filter out cells with the highest and lowest number of 1% 
## - Filter out cells with more than 10% of mitochondrial genes expression.

seurat.obj[["percent.mt"]]  <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(wd,"plots/01_QC/Vln_feat_counts_mt.png"), width=7.20, height=8.03)

## subset based on literature thresholds
# https://www.nature.com/articles/s41422-018-0099-2#Sec19
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & percent.mt < 20)

##-----------3. Analysis -----------##
## Normalize
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(seurat.obj), 10)
## plot variable features with and without labels
# plot1 <- VariableFeaturePlot(seurat.obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

## Scaling
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
DimPlot(seurat.obj, reduction = "umap", label = T) + ggtitle("UMAP per Identity")
ggsave(file.path(wd,"plots/02_UMAPs/umaps_identity.png"), width=7.20, height=8.03)

##---------------- Clustering ----------------------##
seurat.obj.markers <- FindAllMarkers(seurat.obj, only.pos = TRUE)

seurat.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


seurat.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(seurat.obj, features = top10$gene) + NoLegend()


## Reassigning clusters
FeaturePlot(seurat.obj, features = c("DAZL","MAGEA4","UTF1","ID4","FGFR3","KIT","DMRT1","DMRTB1","STRA8","SYCP3","SPO11","MLH3","ZPBP","TNP1","PRM2"))#"VIM","CD14","VWF","ACTA2","DLK1")
FeaturePlot(seurat.obj,features =  c("TNP1","PRM2","SHBG", "SRD5A2"), label=T) # sperm
FeaturePlot(seurat.obj, features = c("SPO11","MLH3"), label=T)#"VIM","CD14","VWF","ACTA2","DLK1")
FeaturePlot(seurat.obj, features = c("UTF1","ID4","FGFR3"), label=T)#"VIM","CD14","VWF","ACTA2","DLK1")
FeaturePlot(seurat.obj, features = c("KIT","DMRT1"), label=T)#"VIM","CD14","VWF","ACTA2","DLK1")
FeaturePlot(seurat.obj, features = c("ACRV1","ZPBP","TSSK6","PRM3"), label=T)#"VIM","CD14","VWF","ACTA2","DLK1")
FeaturePlot(seurat.obj, features = c("STRA8","SYCP3","SPO11","MLH3","ZPBP","TNP1","PRM2"), label=T)#"VIM","CD14","VWF","ACTA2","DLK1")

new.cluster.ids <- c("Sperm", #0
                     "Sperm", #1
                     "Sperm", #2
                     "Leydig Cells", #3
                     "Early Spermatocytes", #4
                     "Elongated Spermatids", #5
                     "Diff. S'gonia", #6
                     "Round Spermatids", #7
                     "Endothelial Cells", #8
                     "Macrophages", #9
                     "Myoid Cells", #10
                     "SSC", #11
                     "Late Spermatocytes", #12
                     "Round Spermatids", #13
                     "Myoid Cells") #14
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj.renamed <- RenameIdents(seurat.obj, new.cluster.ids)
DimPlot(seurat.obj.renamed, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# DoHeatmap(seurat.obj.renamed, features = c("DAZL","MAGEA4","UTF1","ID4","FGFR3","KIT","DMRT1","DMRTB1","STRA8","SYCP3","SPO11","MLH3","ZPBP","TNP1","PRM2")) + NoLegend()

##-------------6. Tumor-specific genes--------------##
candidates = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes/onlyStep1/TSTR_candidatesORFs_fullcharacterized.csv")
candidates$gene_name = gsub("ENSG00000287861","LOC107984132", candidates$gene_name)
CTx = candidates %>% subset(coding_noncoding_chr == "CT-X") %>% select(gene_name) %>% unique() %>% pull(gene_name)
CTnonx = candidates %>% subset(coding_noncoding_chr == "CT-nonX") %>% select(gene_name) %>% unique() %>% pull(gene_name)
NCnonX = candidates %>% subset(coding_noncoding_chr == "Noncoding-nonX") %>% select(gene_name) %>% unique() %>% pull(gene_name)

all_TSTR = candidates %>% subset(coding_noncoding_chr != "altORFs") %>% select(gene_name) %>% unique() %>% pull(gene_name)

dir.create(file.path(wd,"plots/05_TSA"))

# reorder clusters
seurat.obj.renamed@active.ident <- factor(seurat.obj.renamed@active.ident, 
                            levels=c("SSC", 
                                     "Diff. S'gonia",
                                     "Early Spermatocytes",
                                     "Late Spermatocytes",
                                     "Round Spermatids",
                                     "Elongated Spermatids",
                                     "Sperm",
                                     "Macrophages","Endothelial Cells","Myoid Cells", "Leydig Cells"))

## CT-X
DotPlot(seurat.obj.renamed, features = CTx) + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="CT-X TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/CTx.png"))

## CT-nonX
DotPlot(seurat.obj.renamed, features = CTnonx) + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="CT-nonX TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")

ggsave(file.path(wd,"plots/05_TSA/CTnonX.png"))

## Noncoding-X
DotPlot(seurat.obj.renamed, features = NCnonX) + 
  RotatedAxis() + 
  labs(x="Candidate TSAntigen",
       y="Cell Type",
       title="Noncoding-nonX TSTR by testicular cell types") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/Noncoding.png"))


## markers
markers = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/celltype_markers.csv")
DotPlot(seurat.obj.renamed, features = markers$gene_name) + 
  RotatedAxis() + 
  labs(x="Marker",
       y="Cell Type",
       title="Stage Markers") +
  scale_colour_gradient2(low = "#FDE333", mid = "#009796", high = "#4B0055")
ggsave(file.path(wd,"plots/05_TSA/markers.png"))


# FeaturePlot(seurat.obj, reduction = "umap", features = NCnonX) 

### 
markers = read.csv("/datasets/marta/scMultiomicsTestis_GSE235324/celltype_markers.csv")
DoHeatmap(
  seurat.obj.renamed,
  features = markers$gene_name
) +
  scale_fill_gradientn(colors = c("#4682B4", "#FFFFE0", "#CA3433")) +
  NoLegend()

DoHeatmap(
  seurat.obj.renamed,
  features = all_TSTR) +
  scale_fill_gradientn(colors = c("#4682B4", "#FFFFE0", "#CA3433")) +
  NoLegend()

