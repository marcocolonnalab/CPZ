library(Seurat)
library(dplyr)
library(Matrix)
library(readr)

## GATHERING DATA TOGETHER
    
dataFolder <- "~/align"

dataFolders <- list()
dataFolders[[1]] <- paste(dataFolder, "/AS037-1/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[2]] <- paste(dataFolder, "/AS037-2/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[3]] <- paste(dataFolder, "/AS037-3/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[4]] <- paste(dataFolder, "/AS037-4/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[5]] <- paste(dataFolder, "/AS037-5/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[6]] <- paste(dataFolder, "/AS037-6/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[7]] <- paste(dataFolder, "/AS037-7/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[8]] <- paste(dataFolder, "/AS037-8/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[9]] <- paste(dataFolder, "/AS037-9/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[10]] <- paste(dataFolder, "/AS037-10/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[11]] <- paste(dataFolder, "/AS037-11/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[12]] <- paste(dataFolder, "/AS037-12/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[13]] <- paste(dataFolder, "/AS037-13/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[14]] <- paste(dataFolder, "/AS037-14/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[15]] <- paste(dataFolder, "/AS037-15/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[16]] <- paste(dataFolder, "/AS037-16/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[17]] <- paste(dataFolder, "/AS037-17/outs/filtered_feature_bc_matrix", sep="")
dataFolders[[18]] <- paste(dataFolder, "/AS037-18/outs/filtered_feature_bc_matrix", sep="")

fdata <- list()
fdata[[1]] <- Read10X(dataFolders[[1]])
fdata[[2]] <- Read10X(dataFolders[[2]])
fdata[[3]] <- Read10X(dataFolders[[3]])
fdata[[4]] <- Read10X(dataFolders[[4]])
fdata[[5]] <- Read10X(dataFolders[[5]])
fdata[[6]] <- Read10X(dataFolders[[6]])
fdata[[7]] <- Read10X(dataFolders[[7]])
fdata[[8]] <- Read10X(dataFolders[[8]])
fdata[[9]] <- Read10X(dataFolders[[9]])
fdata[[10]] <- Read10X(dataFolders[[10]])
fdata[[11]] <- Read10X(dataFolders[[11]])
fdata[[12]] <- Read10X(dataFolders[[12]])
fdata[[13]] <- Read10X(dataFolders[[13]])
fdata[[14]] <- Read10X(dataFolders[[14]])
fdata[[15]] <- Read10X(dataFolders[[15]])
fdata[[16]] <- Read10X(dataFolders[[16]])
fdata[[17]] <- Read10X(dataFolders[[17]])
fdata[[18]] <- Read10X(dataFolders[[18]])

whole <- list()
whole[[1]] <- CreateSeuratObject(counts=fdata[[1]], project = "A1")
whole[[2]] <- CreateSeuratObject(counts=fdata[[2]], project = "A2")
whole[[3]] <- CreateSeuratObject(counts=fdata[[3]], project = "A3")
whole[[4]] <- CreateSeuratObject(counts=fdata[[4]], project = "A4")
whole[[5]] <- CreateSeuratObject(counts=fdata[[5]], project = "A5")
whole[[6]] <- CreateSeuratObject(counts=fdata[[6]], project = "A6")
whole[[7]] <- CreateSeuratObject(counts=fdata[[7]], project = "A7")
whole[[8]] <- CreateSeuratObject(counts=fdata[[8]], project = "A8")
whole[[9]] <- CreateSeuratObject(counts=fdata[[9]], project = "A9")
whole[[10]] <- CreateSeuratObject(counts=fdata[[10]], project = "A10")
whole[[11]] <- CreateSeuratObject(counts=fdata[[11]], project = "A11")
whole[[12]] <- CreateSeuratObject(counts=fdata[[12]], project = "A12")
whole[[13]] <- CreateSeuratObject(counts=fdata[[13]], project = "A13")
whole[[14]] <- CreateSeuratObject(counts=fdata[[14]], project = "A14")
whole[[15]] <- CreateSeuratObject(counts=fdata[[15]], project = "A15")
whole[[16]] <- CreateSeuratObject(counts=fdata[[16]], project = "A16")
whole[[17]] <- CreateSeuratObject(counts=fdata[[17]], project = "A17")
whole[[18]] <- CreateSeuratObject(counts=fdata[[18]], project = "A18")

#merge all samples
batch1 <- merge(x=whole[[1]], y=list(whole[[2]], whole[[3]], whole[[4]], whole[[5]], whole[[6]], whole[[7]], whole[[8]], whole[[9]], whole[[10]]), add.cell.ids = c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10" ))

batch2 <- merge(x=whole[[11]], y=list(whole[[12]], whole[[13]], whole[[14]], whole[[15]], whole[[16]], whole[[17]], whole[[18]]), add.cell.ids = c("A11" , "A12", "A13", "A14", "A15", "A16", "A17", "A18" ))

## Number of cells before
#cells.before <- length(colnames(x= whole))

cells_plot <- read_tsv("~/dataReCluster.tsv")
cells_subset <- filter(cells_plot, !(Clustering0.8 %in% c(13, 31, 33, 36, 37, 38, 41, 42, 43, 44, 45, 46)))$X1
batch1 <- subset(x = batch1, cells = cells_subset)
batch2 <- subset(x = batch2, cells = cells_subset)


## NORMALIZATION
mito.genes <- c(grep("^MT-", rownames(x = batch1), value = T),
                grep("^mt-", rownames(x = batch1), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = batch1, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = batch1, slot = 'counts'))
batch1[['percent.mito']] <- percent.mito

batch1 <- subset(x = batch1, subset = percent.mito <= 0.05 & nFeature_RNA > 400 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 15000)

mito.genes <- c(grep("^MT-", rownames(x = batch2), value = T),
                grep("^mt-", rownames(x = batch2), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = batch2, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = batch2, slot = 'counts'))
batch2[['percent.mito']] <- percent.mito

batch2 <- subset(x = batch2, subset = percent.mito <= 0.05  & nFeature_RNA > 400 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 15000)


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
whole <- FindClusters(object = whole, resolution = 0.9)

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
