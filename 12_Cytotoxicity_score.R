# Calculate cytotoxicity function of CD8+ T cells.

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(ggthemes)
library(ggbeeswarm)
library(scales)
library(forcats)

kill_geneset <- list(c("PRF1","IFNG","GZMA","GZMB","GZMK","GZMH","GNLY","NKG7","KLRK1","KLRB1","KLRD1","CTSW","CST7"))

good_samples <- c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17")
cd8t_clusters <- c("CD8_C1","CD8_C2","CD8_C3","CD8_C4","CD8_C5","CD8_C6","CD8_C7")

setwd("/mnt/sdc/singlecell/data/")
file1 <- readRDS("T_total_3_26.rds")
file1 <- subset(file1, subset = orig.ident %in% good_samples)
file1 <- subset(file1, subset = Cluster %in% cd8t_clusters)
file1 <- AddModuleScore(object=file1, features=kill_geneset, name="Killing_score")

file1[["group"]] <- ifelse(file1$orig.ident %in% c("T1","T5","T6","T8","T10","T12","T13","T14","T16","T17"),"LUAD","AISMIA")
head(file1@meta.data)

ggplot(data=file1@meta.data, aes(x=group, y=Killing_score1, fill=group)) +
  geom_boxplot() +
  labs(x="Group", y = "Cytotoxicity_score", fill="Group") +
  #scale_fill_brewer(palette="Blues") + 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  stat_compare_means()+
  scale_fill_brewer(palette="Set1")
ggsave("/mnt/sdc/singlecell/analysis/Killing_score.pdf", plot=last_plot(), width=6, height=8, dpi=600)
