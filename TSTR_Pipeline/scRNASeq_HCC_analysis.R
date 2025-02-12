## Semi-Automatized single-cell analysis of Hepatocellular Carcinoma datasets
## @author: Marta Espinosa
## @date: 27th May 2024

##---------------- Load libraries ----------------##
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
# library(patchwork)
wd = "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/scRNASeq/HCC"


##-----------0. Importing Seurat Object -----------##
## Created by data2SeuratObj.R
datasets_directories = list.dirs("/projects_eg/projects/marta/HCC_scRNASeq", recursive = F)[-1]
datasets = gsub("..*/","",datasets_directories)

for(i in length(datasets_directories)) {
  dataset = datasets[i]
  print(dataset)
  seurat.obj = readRDS(file = file.path(datasets_directories[i],"seurat_object.rds"))
  ## how many cells per cell type and patient?
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$patient)
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$site)
  # table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$site, seurat.obj@meta.data$patient)
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$patient, seurat.obj@meta.data$site)
  
  
  ##--------------1. Quality Control ---------------##
  dir.create(file.path(wd,dataset,"plots/01_QC"), recursive = T)
  ## - Filter out cells with the highest and lowest number of 1% 
  ## - Filter out cells with more than 10% of mitochondrial genes expression.
  
  seurat.obj[["percent.mt"]]  <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_feat_counts_mt.png"), width=7.20, height=8.03)
  
  VlnPlot(seurat.obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, group.by = "CellType", pt.size = 0) +
    geom_hline(yintercept=5000)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_feat_counts_mt.celltype.png"), width=7.20, height=8.03)
  
  VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "CellType", pt.size = 0) + geom_hline(yintercept=6000)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_counts_threshold6000.celltype.png"), width=7.20, height=8.03)
  
  VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample", pt.size = 0)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_feat_counts_mt.sample.png"), width=7.20, height=8.03)
  
  VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "Sample", pt.size = 0) + geom_hline(yintercept=6000)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_counts_threshold6000.samples.png"), width=7.20, height=8.03)
  
  VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "site", pt.size = 0) + geom_hline(yintercept=6000)
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_counts_threshold6000.site.png"), width=7.20, height=8.03)
  
  ## subset based on literature thresholds
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 200 & nCount_RNA < 6000 & percent.mt < 10)
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$patient)
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$site)
  # table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$site, seurat.obj@meta.data$patient)
  table(seurat.obj@meta.data$CellType, seurat.obj@meta.data$patient, seurat.obj@meta.data$site)
  ## Filtered | 200 > nFeat > 2500 | 200 > nCount > 6000 | percent.mt < 10%
  VlnPlot(seurat.obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, group.by = "CellType", pt.size = 0) 
  ggsave(file.path(wd,dataset,"plots/01_QC/Vln_feat_counts_mt.filtered.celltype.png"), width=7.20, height=8.03)
  
  ##-----------3. Analysis -----------##
  ## Normalization
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  
  ## Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat.obj), 10)
  LabelPoints(plot = VariableFeaturePlot(seurat.obj), points = top10, repel = TRUE)
  ggsave(file.path(wd,dataset,"plots/01_QC/VariableFeaturePlot.png"), width=7.20, height=8.03)
  
  ## prior to PCA
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  ## PCA
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  
  ## num clusters
  ElbowPlot(seurat.obj)
  eigs<-seurat.obj[["pca"]]@stdev^2
  seurat.obj[["pca"]]@misc$proportion<-eigs/sum(eigs)
  seurat.obj[["pca"]]@misc$cumulative<-cumsum(eigs)/sum(eigs)
  seurat.obj[["pca"]]@misc$cumulative #to check which one explains >90% variation
  
  SelPC<-min(which(seurat.obj[["pca"]]@misc$cumulative>0.9))
  print(SelPC)
  
  ## Clustering
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:SelPC)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
  
  ## UMAP
  dir.create(file.path(wd,dataset,"plots/02_UMAPs"))
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:SelPC)
  DimPlot(seurat.obj, reduction = "umap") + ggtitle(paste0("UMAP per Identity | ",dataset))
  ggsave(file.path(wd,dataset,"plots/02_UMAPs/umaps_identity.png"), width=7.20, height=8.03)
  
  DimPlot(seurat.obj, reduction = "umap", group.by = "CellType") + ggtitle(paste0("UMAP per Cell Type | ",dataset))
  ggsave(file.path(wd,dataset,"plots/02_UMAPs/umaps_celltypes.png"), width=7.20, height=8.03)
  
  DimPlot(seurat.obj, reduction = "umap", group.by = "CellType", split.by = "patient") + ggtitle(paste0("UMAP per Cell Type and Patient | ",dataset))
  ggsave(file.path(wd,dataset,"plots/02_UMAPs/umaps_celltypes_patient.png"), width=15.15, height=8.03)
  
  if (dataset == "GSE149614" | dataset == "GSE151530") {
    DimPlot(subset(seurat.obj, subset = CellType == "Hepatocyte"), reduction = "umap", group.by = "CellType", split.by = "patient") + ggtitle(paste0("UMAP per Cell Type and Patient | ",dataset))
    ggsave(file.path(wd,dataset,"plots/02_UMAPs/umaps_malignantHep_patient.png"), width=15.15, height=8.03)
  } else {
    DimPlot(subset(seurat.obj, subset = CellType == "Malignant cell"), reduction = "umap", group.by = "CellType", split.by = "patient") + ggtitle(paste0("UMAP per Cell Type and Patient | ",dataset))
    ggsave(file.path(wd,dataset,"plots/02_UMAPs/umaps_malignantcells_patient.png"), width=15.15, height=8.03)
  }
  
  ##------------------4. Integration???-----------------##
  ## If integration is needed, perform it. Otherwise, comment until 5. Clusters reassignment
  dir.create(file.path(wd,dataset,"plots/03_Integration"))
  
  ## reload seurat object
  toIntegrate_seurat.obj = readRDS(file =  file.path(datasets_directories[i],"seurat_object.rds"))
  toIntegrate_seurat.obj[["percent.mt"]]  <- PercentageFeatureSet(toIntegrate_seurat.obj, pattern = "^MT-")
  toIntegrate_seurat.obj <- subset(toIntegrate_seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 200 & nCount_RNA < 6000 & percent.mt < 10)
  table(toIntegrate_seurat.obj@meta.data$CellType, toIntegrate_seurat.obj@meta.data$patient)
  
  ## Only for the patients with enough malignant cells
  if (dataset == "GSE149614") {
    malignantcells = subset(toIntegrate_seurat.obj, subset = CellType == "Hepatocyte")
  } else {
    malignantcells = subset(toIntegrate_seurat.obj, subset = CellType == "Malignant cell")
  }
  table(malignantcells@meta.data$CellType, malignantcells@meta.data$patient)
  patients = as.data.frame(table(malignantcells@meta.data$patient))
  patients = patients %>% subset(Freq > 100)
  patients_to_keep = as.vector(patients$Var1)
  
  seurat.obj.patients = subset(toIntegrate_seurat.obj, subset = patient %in% patients_to_keep)
  
  # split the dataset into a list of seurat objects (per patient)
  ifnb.list <- SplitObject(seurat.obj.patients, split.by = "patient")
  
  # normalize and identify variable features for each dataset independently
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) # this command creates an 'integrated' data assay
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  
  ## Now we can run a single integrated analysis on all cells!
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(immune.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:SelPC)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:SelPC)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)
  
  DimPlot(immune.combined, reduction = "umap", group.by = "patient")
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAPs_identity.png"), width=7.20, height=8.03)
  
  DimPlot(immune.combined, reduction = "umap", group.by = "CellType", label = TRUE,
          repel = TRUE)
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAPs_CellType.png"), width=7.20, height=8.03)
  
  DimPlot(immune.combined, reduction = "umap", group.by = "CellType") + ggtitle(paste0("Integrated UMAP | ",dataset))
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAP.png"), width=8.56, height=4.59)
  
  
  DimPlot(immune.combined, reduction = "umap", split.by = "patient", group.by = "CellType") + ggtitle(paste0("Integrated UMAP per Patient | ",dataset))
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAP_CellType_Patient.png"), width=8.56, height=4.59)
  
  ## Select malignant cells
  
  if (dataset == "GSE149614") {
    immune.combined_malignantcells = subset(immune.combined, subset = CellType == "Hepatocyte")
  } else {
    immune.combined_malignantcells = subset(immune.combined, subset = CellType == "Malignant cell")
  }
  
  DimPlot(immune.combined_malignantcells, reduction = "umap", group.by = "patient") + ggtitle(paste0("Integrated malignant cells UMAP per Patient | ",dataset))
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAP_Patient_mixed.png"), width=8.56, height=4.59)
  
  DimPlot(immune.combined_malignantcells, reduction = "umap", group.by = "CellType", split.by = "patient") + ggtitle(paste0("Integrated malignant cells UMAP per Patient | ",dataset))
  ggsave(file.path(wd,dataset,"plots/03_Integration/IntegratedUMAP_MalignantCells_Patient.png"), width=8.56, height=4.59)
  
  
  ##-------------5. Clusters reassignment--------------##
  ## If we integrated:
  seurat.obj = immune.combined
  DefaultAssay(seurat.obj) <- "RNA"
  ## If we did not integrate:
  # seurat.obj = seurat.obj
  # DefaultAssay(seurat.obj) <- "RNA"
  
  markers = read.csv(file.path("/projects_eg/projects/marta/HCC_scRNASeq/cell_markers_HCC_Xie.csv"))
  HCC_markers = markers %>% subset(CellType == "HCC")
  dir.create(file.path(wd,dataset,"plots/04_Clustering"))
  
  ## defining new celltypes wrt file
  Bcells = markers %>% subset(CellType == "B") %>% pull(CellMarker) %>% unique
  Tcells = markers %>% subset(CellType == "CD8.T" | CellType == "NK" | CellType == "NKT" | CellType == "Naive.T") %>% pull(CellMarker) %>% unique
  En = markers %>% subset(CellType == "En") %>% pull(CellMarker) %>% unique
  Ep = markers %>% subset(CellType == "Ep") %>% pull(CellMarker) %>% unique
  Neutrophils = markers %>% subset(CellType == "Neutrophils") %>% pull(CellMarker)
  M = markers %>% subset(CellType == "Mac") %>% pull(CellMarker) %>% unique
  pDC = markers %>% subset(CellType == "pDC" | CellType == "MDDC") %>% pull(CellMarker) %>% unique
  HCC = markers %>% subset(CellType == "HCC") %>% pull(CellMarker) %>% unique
  ILC = markers %>% subset(CellType == "ILC") %>% pull(CellMarker) %>% unique
  Fibro = markers %>% subset(CellType == "Fibroblasts") %>% pull(CellMarker) %>% unique
  
  markers_list = list("B-cell" = Bcells,
                      "T-cell" = Tcells,
                      "Endothelials" = En,
                      "Epithelials" = Ep,
                      "pDC" = pDC,
                      "Neutrophils" = Neutrophils,
                      "Macrophages" = M,
                      "Fibroblasts" = Fibro,
                      "HCC" = HCC,
                      "ILC" = ILC)
  remove_items <- function(lst) {
    if ("EPCAM" %in% lst) {
      lst <- lst[lst != "EPCAM"]
    } 
    if ("PTPRC" %in% lst) {
      lst <- lst[lst != "PTPRC"]
    } 
    return(lst)
  }
  
  # Apply the function to each list in markers_list
  markers_list <- lapply(markers_list, remove_items)
  dotplot_marker = DotPlot(object = seurat.obj,
                           features = markers_list,
                           col.min=0,
                           col.max = 1,
                           cluster.idents = T) +
    theme(axis.text.x.bottom = element_text(angle=90),
          axis.text.x = element_text(angle=90)) +
    labs(title="DotPlot reassigning cell types based on Xie et al.")
  ggsave(file.path(wd,dataset,"plots/04_Clustering/dotplot_XieMarkers_identity.png"), width=15.15, height=8.03)
  
  dotplot_marker = DotPlot(object = seurat.obj,
                           features = markers_list,
                           col.min=0,
                           col.max = 1,
                           cluster.idents = T, group.by = "CellType") +
    theme(axis.text.x.bottom = element_text(angle=90),
          axis.text.x = element_text(angle=90)) +
    labs(title="DotPlot Assigned CellTypes vs cell types based on Xie et al.")
  ggsave(file.path(wd,dataset,"plots/04_Clustering/dotplot_XieMarkers_CellTypes.png"), width=15.15, height=8.03)
  
  ##-------------6. TSTR--------------##
  candidates = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes/onlyStep1/TSTR_candidatesORFs_fullcharacterized.csv")
  genes_candidates = unique(candidates$gene_name)
  dir.create(file.path(wd,dataset,"plots/05_TSA"))
  
  DotPlot(seurat.obj, features = genes_candidates, group.by = "CellType") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("TSTR HCC Candidates expression by cell types | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CellTypes.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = genes_candidates, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("TSTR HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_Patient.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = genes_candidates, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("TSTR HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_Site.png"), width=11.56, height=4.59)
  
  ### CT-X
  CTX = candidates %>% subset(coding_noncoding_chr == "CT-X") %>% pull(gene_name)
  CTX = unique(CTX)
  DotPlot(seurat.obj, features = CTX, group.by = "CellType") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTX HCC Candidates expression by cell types | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTX_CellTypes.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = CTX, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTX HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTX_Patient.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = CTX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTX HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTX_Site.png"), width=11.56, height=4.59)

  ### CT-nonX
  CTnonX = candidates %>% subset(coding_noncoding_chr == "CT-nonX") %>% pull(gene_name)
  CTnonX = unique(CTnonX)
  DotPlot(seurat.obj, features = CTnonX, group.by = "CellType") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTnonX HCC Candidates expression by cell types | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTnonX_CellTypes.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = CTnonX, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTnonX HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTnonX_Patient.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = CTnonX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTnonX HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTnonX_Site.png"), width=11.56, height=4.59)
  
  ### Noncoding-nonX
  NCnonX = candidates %>% subset(coding_noncoding_chr == "Noncoding-nonX") %>% pull(gene_name)
  NCnonX = unique(NCnonX)
  DotPlot(seurat.obj, features = NCnonX, group.by = "CellType") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("Noncoding-nonX HCC Candidates expression by cell types | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_NCnonX_CellTypes.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = NCnonX, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("Noncoding-nonX HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_NCnonX_Patient.png"), width=11.56, height=4.59)
  
  DotPlot(seurat.obj, features = NCnonX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("Noncoding-nonX HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_NCnonX_Site.png"), width=11.56, height=4.59)
  
  
  # if(dataset == "GSE189903") {
  #   P1 = subset(seurat.obj, subset = patient == "P1")
  #   DotPlot(P1, features = genes_candidates, group.by = "CellType") + 
  #     RotatedAxis() + 
  #     labs(x="Candidate TSAntigen",
  #          y="Cell Type",
  #          title=paste0("TSA HCC Candidates expression by cell types | ", dataset))
  #   ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSA_CellTypes.P1.png"), width=11.56, height=4.59)
  #   
  #   DotPlot(P1, features = genes_candidates, group.by = "site") + 
  #     RotatedAxis() + 
  #     labs(x="Candidate TSAntigen",
  #          y="Cell Type",
  #          title=paste0("TSA HCC Candidates expression by site | ", dataset))
  #   ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSA_Site.P1.png"), width=11.56, height=4.59)
  # }
  
  ## only malignant cells - dotplot
  if (dataset == "GSE149614") {
    malignantcells = subset(seurat.obj, subset = CellType == "Hepatocyte")
  } else {
    malignantcells = subset(seurat.obj, subset = CellType == "Malignant cell")
  }
  
  DotPlot(malignantcells, features = genes_candidates, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("TSA HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_Patient_MalignantCells.png"), width=11.56, height=4.59)

  ## CTX
    DotPlot(malignantcells, features = CTx, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTx HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTx_Patient_MalignantCells.png"), width=11.56, height=4.59)
  ## CTnonX
  DotPlot(malignantcells, features = CTnonX, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTnonX HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTnonX_Patient_MalignantCells.png"), width=11.56, height=4.59)
  ## NoncodingnonX
  DotPlot(malignantcells, features = NCnonX, group.by = "patient") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("NCnonX HCC Candidates expression by patients | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_NCnonX_Patient_MalignantCells.png"), width=11.56, height=4.59)
  
  
  DotPlot(malignantcells, features = genes_candidates, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("TSTR HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_Site_MalignantCells.png"), width=11.56, height=3.59)
  ## CTX
  DotPlot(malignantcells, features = CTX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTx HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTX_Site_MalignantCells.png"), width=11.56, height=3.59)
  ## CTnonX
  DotPlot(malignantcells, features = CTnonX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("CTnonX HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_CTnonX_Site_MalignantCells.png"), width=11.56, height=3.59)
  ## NoncodingnonX
  DotPlot(malignantcells, features = NCnonX, group.by = "site") + 
    RotatedAxis() + 
    labs(x="Candidate TSAntigen",
         y="Cell Type",
         title=paste0("NCnonX HCC Candidates expression by site | ", dataset))
  ggsave(file.path(wd,dataset,"plots/05_TSA/dotplot_TSTR_NCnonX_Site_MalignantCells.png"), width=11.56, height=3.59)

  ## Plot expression level
  VlnPlot(malignantcells, features = CTX, split.by = "patient", group.by = "patient",
          pt.size = 0) 
  ggsave(file.path(wd,dataset,"plots/05_TSA/ExpressionLvl_MalignantCells_CTX.png"), width=15.15, height=15.03)
  VlnPlot(malignantcells, features = CTnonX, split.by = "patient", group.by = "patient",
          pt.size = 0) 
  ggsave(file.path(wd,dataset,"plots/05_TSA/ExpressionLvl_MalignantCells_CTnonX.png"), width=15.15, height=15.03)
  VlnPlot(malignantcells, features = NCnonX, split.by = "patient", group.by = "patient",
          pt.size = 0) 
  ggsave(file.path(wd,dataset,"plots/05_TSA/ExpressionLvl_MalignantCells_NCnonX.png"), width=15.15, height=15.03)
  
  
  if (dataset == "GSE151530") {
    DefaultAssay(malignantcells) <- "RNA"
    
    FeaturePlot(malignantcells, reduction = "umap", features = CTX) 
    ## Patient H70
    ## CTX
    FeaturePlot(subset(malignantcells, subset = patient == "H70"), reduction = "umap", features = CTX)
    ggsave(file.path(wd,dataset,"plots/05_TSA/FeaturePlot_MalignantCells_CTx_H70.png"), width=15.15, height=8.03)
    ## CTnonX
    FeaturePlot(subset(malignantcells, subset = patient == "H70"), reduction = "umap", features = CTnonX)
    ggsave(file.path(wd,dataset,"plots/05_TSA/FeaturePlot_MalignantCells_CTnonX_H70.png"), width=15.15, height=8.03)
    ## NCnonX
    FeaturePlot(subset(malignantcells, subset = patient == "H70"), reduction = "umap", features = NCnonX)
    ggsave(file.path(wd,dataset,"plots/05_TSA/FeaturePlot_MalignantCells_NCnonX_H70.png"), width=15.15, height=8.03)
    
    } else if (dataset == "GSE189903") {
    DefaultAssay(malignantcells) <- "RNA"
    ## too few malignan cells

  } else if (dataset == "GSE149614") {
    DefaultAssay(malignantcells) <- "RNA"

  }
  
  ## save counts
  # write.table(as.matrix(GetAssayData(object = seurat.obj, slot = "counts")), 
  #             '/your/path/to/store/csv/counts.csv', 
  #             sep = ',', row.names = T, col.names = T, quote = F)
    
  }
  
