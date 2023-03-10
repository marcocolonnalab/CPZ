library(Seurat)
library(dplyr)
library(Matrix)
library(readr)

## GATHERING DATA TOGETHER
    
dataFolder <- "~/align3"

dataFolders <- list()
dataFolders[[1]] <- paste(dataFolder, "/A1/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[2]] <- paste(dataFolder, "/A2/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[3]] <- paste(dataFolder, "/A3/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[4]] <- paste(dataFolder, "/A4/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[5]] <- paste(dataFolder, "/A5/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[6]] <- paste(dataFolder, "/B7/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[7]] <- paste(dataFolder, "/B8/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[8]] <- paste(dataFolder, "/B9/outs/filtered_feature_bc_matrix", sep="")

fdata <- list()
fdata[[1]] <- Read10X(dataFolders[[1]])
fdata[[2]] <- Read10X(dataFolders[[2]])
fdata[[3]] <- Read10X(dataFolders[[3]])
fdata[[4]] <- Read10X(dataFolders[[4]])
fdata[[5]] <- Read10X(dataFolders[[5]])
fdata[[6]] <- Read10X(dataFolders[[6]])
fdata[[7]] <- Read10X(dataFolders[[7]])
fdata[[8]] <- Read10X(dataFolders[[8]])

whole <- list()
whole[[1]] <- CreateSeuratObject(counts=fdata[[1]], project = "A1")
whole[[2]] <- CreateSeuratObject(counts=fdata[[2]], project = "A2")
whole[[3]] <- CreateSeuratObject(counts=fdata[[3]], project = "A3")
whole[[4]] <- CreateSeuratObject(counts=fdata[[4]], project = "A4")
whole[[5]] <- CreateSeuratObject(counts=fdata[[5]], project = "A5")
whole[[6]] <- CreateSeuratObject(counts=fdata[[6]], project = "B7")
whole[[7]] <- CreateSeuratObject(counts=fdata[[7]], project = "B8")
whole[[8]] <- CreateSeuratObject(counts=fdata[[8]], project = "B9")

#merge all samples
batch1 <- merge(x=whole[[1]], y=list(whole[[2]], whole[[3]], whole[[4]], whole[[5]]), add.cell.ids = c("A1" , "A2", "A3", "A4", "A5" ))

batch2 <- merge(x=whole[[6]], y=list(whole[[7]], whole[[8]]), add.cell.ids = c("B7", "B8", "B9" ))

## Number of cells before
#cells.before <- length(colnames(x= whole))

cells_plot <- read_tsv("~/dataReCluster.tsv")
cells_subset <- filter(cells_plot, !(Clustering0.9 %in% c(17, 24, 26)))$X1
batch1 <- subset(x = batch1, cells = cells_subset) 
batch2 <- subset(x = batch2, cells = cells_subset)

mito.genes <- c(grep("^MT-", rownames(x = batch1), value = T),
                grep("^mt-", rownames(x = batch1), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = batch1, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = batch1, slot = 'counts'))
batch1[['percent.mito']] <- percent.mito

batch1 <- subset(x = batch1, subset = percent.mito <= 0.05 & nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000)

mito.genes <- c(grep("^MT-", rownames(x = batch2), value = T),
                grep("^mt-", rownames(x = batch2), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = batch2, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = batch2, slot = 'counts'))
batch2[['percent.mito']] <- percent.mito

batch2 <- subset(x = batch2, subset = percent.mito <= 0.05 & nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000) 


# combine datasets together
whole.list <- list("batch1" = batch1, "batch2" = batch2)

for (i in 1:length(whole.list)) {
    whole.list[[i]] <- NormalizeData(whole.list[[i]], verbose = FALSE)
    whole.list[[i]] <- FindVariableFeatures(whole.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}

reference.list <- whole.list[c("batch1", "batch2")]
whole.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
whole <- IntegrateData(anchorset = whole.anchors, dims = 1:30)


DefaultAssay(whole) <- "integrated"


whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mito"))
gc()

## PCA
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                verbose=FALSE)


## TSNE
whole <- RunTSNE(object = whole, dims = 1:20, check_duplicates = FALSE)
whole <- RunUMAP(object = whole, dims = 1:20)

## CLUSTERING
whole <- FindNeighbors(object = whole, dims = 1:20)
whole <- FindClusters(object = whole, resolution = 0.3)

DimPlot(object = whole, reduction="tsne")

## SAVING DATA FOR EXPLORER
expData <- GetAssayData(object = whole, slot = 'data')
save(expData, file="expData.Rda")

dataForPlot <- as.data.frame(whole@reductions$tsne@cell.embeddings)
dataForPlot$Sample <- whole@meta.data$orig.ident
dataForPlot$Cluster <-  Idents(object = whole)
dataForPlot$nUmi <- whole@meta.data$nCount_RNA
dataForPlot$nGene <- whole@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(whole@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(whole@meta.data$nFeature_RNA)

write.table(dataForPlot, "data_for_plot.tsv", sep="\t", quote=F)

DefaultAssay(whole) <- "RNA"

## FINDING ANS SAVING MARKERS
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                logfc.threshold = 0.10)
write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

## SAVING
save(whole, file = "whole_object.Robj")

## Number of cells after
cells.after <- length(colnames(x = whole))
print(paste0("cells.before:",cells.before))
print(paste0("cells.after:",cells.after))
print(paste0("cell.diff:", cells.before-cells.after))
