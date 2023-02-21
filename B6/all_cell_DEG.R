#work directory check and set

setwd(dir = "~/all_celltype")

library(tidyverse)
library(Seurat)
library(cowplot)
library(slingshot)
library(tradeSeq)
library(metap)
library(DropletUtils)
library(RColorBrewer)

load("~/whole_object.Robj")

median <- whole@meta.data %>%
  select(nCount_RNA,nFeature_RNA) %>%
  summarise(nUMI=median(nCount_RNA), nGene=median(nFeature_RNA))

whole@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800)

Idents(whole) <- "orig.ident"
whole$condition <- NA
whole$condition[which(str_detect(whole$orig.ident, "A1|A2"))] <- "Normal"
whole$condition[which(str_detect(whole$orig.ident, "A3|A4|A5"))] <- "Demye"
whole$condition[which(str_detect(whole$orig.ident, "B7|B8|B9"))] <- "Remye"
DimPlot(whole, label = T, split.by = "condition")


whole$celltype <- NA
whole$celltype[which(whole$integrated_snn_res.0.3 %in% c(0,1,5,12,15,16))] <- "Ex_neuron"
whole$celltype[which(whole$integrated_snn_res.0.3 %in% c(2,6,8,9,11))] <- "In_neuron"
whole$celltype[which(whole$integrated_snn_res.0.3 %in% c(3,17))] <- "Oligo"
whole$celltype[which(whole$integrated_snn_res.0.3 == 7)] <- "Astro"
whole$celltype[which(whole$integrated_snn_res.0.3 == 4)] <- "Micro"
whole$celltype[which(whole$integrated_snn_res.0.3 == 10)] <- "OPC"
whole$celltype[which(whole$integrated_snn_res.0.3 %in% c(13,14,18))] <- "Vascular_cells"
DimPlot(object = whole, reduction = "umap", label = T, group.by = "celltype")


Idents(whole) <- "celltype"
whole$celltype.condition <- paste(Idents(whole), whole$condition, sep = "_")
Idents(whole) <- "celltype.condition"
#check all the idents
levels(whole)

# DEG: Demye vs Nor
Idents(whole) <- "celltype.condition"
levels(whole)
de <- data.frame()
DefaultAssay(whole) <- "RNA"
for (i in unique(whole@meta.data$celltype)){
  demye <- paste(i, "Demye", sep = "_")
  normal <- paste(i, "Normal", sep = "_")
  de.genes <- FindMarkers(whole, ident.1 = demye, ident.2 = normal)
  de.genes$cluster <- i
  de.genes$gene <- rownames(de.genes)
  de <- rbind(de, de.genes)
}
write.table(de,"~/demye_vs_normal.txt", sep= "\t")
d_n <- read.table("~/demye_vs_normal.txt", sep= "\t")

# DEG: Remye vs Nor
de <- data.frame()
DefaultAssay(whole) <- "RNA"
for (i in unique(whole@meta.data$celltype)){
  remye <- paste(i, "Remye", sep = "_")
  normal <- paste(i, "Normal", sep = "_")
  de.genes <- FindMarkers(whole, ident.1 = remye, ident.2 = normal)
  de.genes$cluster <- i
  de.genes$gene <- rownames(de.genes)
  de <- rbind(de, de.genes)
}
write.table(de,"~/remye_vs_normal.txt", sep= "\t")
r_n <- read.table("~/remye_vs_normal.txt", sep= "\t")

# DEG: Remye vs Demye
de <- data.frame()
DefaultAssay(whole) <- "RNA"
for (i in unique(whole@meta.data$celltype)){
  remye <- paste(i, "Remye", sep = "_")
  demye <- paste(i, "Demye", sep = "_")
  de.genes <- FindMarkers(whole, ident.1 = remye, ident.2 = demye)
  de.genes$cluster <- i
  de.genes$gene <- rownames(de.genes)
  de <- rbind(de, de.genes)
}
write.table(de,"~/remye_vs_demye.txt", sep= "\t")
r_d <- read.table("~/remye_vs_demye.txt", sep= "\t")



Idents(whole) <- "integrated_snn_res.0.3"
cluster.averages <- AverageExpression(whole, add.ident = "orig.ident")

cluster.averages.rna <- cluster.averages$RNA
write.table(cluster.averages.rna,"~/avg_exp_rna.tsv", sep= "\t")



