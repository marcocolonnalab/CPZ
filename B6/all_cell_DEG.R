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



## volcano plots
d_n <- read.table("~/demye_vs_normal.txt", sep= "\t")
#remove "mm10---"
d_n$gene <- gsub("mm10---","", d_n$gene)
# set padj, log2FC threshold value (check the csv colume name)
up <- d_n %>% 
  dplyr::filter(cluster == "Micro") %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5)
down <- d_n %>% 
  dplyr::filter(cluster == "Micro") %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < (-0.5))
sig <- d_n %>% 
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  top_n(wt = abs(avg_log2FC), n = 100)
anno_list <- sig$gene
# anno_list <- c(filter(de.vs.nor p_val_adj < 0.05, avg_log2FC > 1)$gene, filter(de.vs.nor, p_val_adj < 0.05, avg_log2FC < (-0.5))$gene)
# add diffexpressed colum to the table
d_n$diffexpressed <- "NO"
d_n$diffexpressed[d_n$avg_log2FC > 0.5 & d_n$p_val_adj < 0.05] <- "UP"
d_n$diffexpressed[d_n$avg_log2FC < (-0.5) & d_n$p_val_adj < 0.05] <- "DOWN"
# plot volcano
library(ggrepel)
d_n  %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj)), color = diffexpressed) +
  geom_point(size = 0.5) +
  geom_point(data = up, color = "#e42625") +
  geom_point(data = down, color = "#3484bc") +
  geom_text_repel(data = dplyr::filter(de.vs.nor, gene %in% anno_list),
                  aes(label=gene),box.padding = 0.4, color = "black", max.overlaps = 50) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("demye vs. normal") +
  ylab("-log10(adj. p-value)") +
  xlab("log2(Fold Change)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# save pdf
ggsave(filename = "~micro_d_n.pdf",device = "pdf", height = 9, width = 9)

