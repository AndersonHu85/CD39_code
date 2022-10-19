library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(ggsci)
options(stringsAsFactors = F)


mito.genes <- grep(
  pattern = "^mt-",
  x = rownames(x = merge@assays$RNA@data),
  value = TRUE)

percent.mito <- Matrix::colSums(merge@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(merge@assays$RNA@counts)

merge <- AddMetaData(
  object = merge,
  metadata = percent.mito,
  col.name = "percent.mito")
VlnPlot(merge, c("nCount_RNA","nFeature_RNA","percent.mito"), group.by = "sample",pt.size = 0)


merge <- subset(merge, subset = percent.mito < 0.1)
merge <- subset(merge, subset = nFeature_RNA < 7500)


merge <- NormalizeData(merge)

merge <- FindVariableFeatures(merge)
VariableFeaturePlot(merge)
merge <- ScaleData(merge, verbose = T)
merge <- RunPCA(merge, npcs = 100, verbose = FALSE)
PCAPlot(merge)
ElbowPlot(merge, ndims = 100)

merge <- RunHarmony(merge, c("sample"))
merge <- FindNeighbors(merge, reduction = "harmony", dims = 1:50, do.plot = T)
merge <- FindClusters(merge, resolution = 0.8)

merge <- RunTSNE(merge, reduction = "harmony", dims = 1:50)

merge <- RunUMAP(merge, reduction = "harmony", dims = 1:50)

celltype <- c(brewer.pal(11,"Set3")[c(1,3:8,10:11)],pal_aaas(alpha = 0.7)(10)[-2],brewer.pal(8,"Set2"),brewer.pal(11,"Paired"))#[13:1]

DimPlot(merge, reduction = "umap", group.by = "celltype2", label = F, repel = TRUE, cols = celltype, pt.size = 0.5
        #,split.by = "sample"
        )

