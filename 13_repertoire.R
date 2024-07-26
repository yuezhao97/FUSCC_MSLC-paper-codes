########### Serial analyses using scRepertoire. ###########
########### Created by Yue Zhao on 2023-03-15. ###########
########### Last modified by Yue Zhao on 2024-07-24. ###########

library(scRepertoire)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)

setwd("/Users/zhaoy2/Desktop/sc_project/scRepertoire")

s1 <- read.csv("input/FD1_LC3_filtered_contig_annotations.csv",skipNul = TRUE)
s2 <- read.csv("input/FD1_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s3 <- read.csv("input/FD2_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s4 <- read.csv("input/FD2_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s5 <- read.csv("input/FD4_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s6 <- read.csv("input/FD4_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s7 <- read.csv("input/FD5_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s8 <- read.csv("input/FD5_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s9 <- read.csv("input/FD8_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s10 <- read.csv("input/FD8_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s11 <- read.csv("input/FD9_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s12 <- read.csv("input/FD9_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s13 <- read.csv("input/FD14_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s14 <- read.csv("input/FD14_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
s15 <- read.csv("input/FD16_LC1_filtered_contig_annotations.csv",skipNul = TRUE)
s16 <- read.csv("input/FD16_LC2_filtered_contig_annotations.csv",skipNul = TRUE)
c_list <- list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16)
head(c_list[[10]])
s10[1:5,1:5]
colnames(c_list[[10]])
colnames(c_list[[2]])
combined <- combineTCR(c_list, 
	samples = c("FD1_LUAD","FD1_MIA","FD2_LUAD","FD2_MIA","FD4_LUAD","FD4_MIA","FD5_LUAD","FD5_MIA","FD8_LUAD1","FD8_LUAD2","FD9_LUAD","FD9_MIA","FD14_LUAD","FD14_AIS","FD16_LUAD1","FD16_LUAD2"),
	ID = c("LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD"),
	cells ="T-AB")
example <- addVariable(combined, name="group",
                       variables=c("LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD"))
example[[1]][1:5,ncol(example[[1]])]
#subset <- subsetContig(combined, name = "sample", variables = c("FD1", "FD2"))

# Quantify clonotypes:
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = T, exportTable = T)
quantContig_output
quantContig_output$group <- c("LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD","LUAD","AISMIA","LUAD","AISMIA","LUAD","LUAD")
quantContig_output
ggsave("quantcontig_genent.pdf", plot=last_plot(), width=16, height=6, dpi=300)
quantContig(combined, cloneCall="gene+nt", scale = F)
ggsave("quantcontig_genent_unscaled.pdf", plot=last_plot(), width=16, height=6, dpi=300)
quantContig_output
write.table(quantContig_output, file="5b.txt", sep="\t", quote=F, row.names=F)
median(quantContig_output[quantContig_output$group=="AISMIA",]$scaled)
ggplot(data=quantContig_output, aes(x=group, y=scaled, color=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  labs(x="Group", y = "Normalized number of TCR clonotypes", color="Group") +
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
  scale_color_brewer(palette="Set1")
ggsave("boxplot_scaled.pdf", plot=last_plot(), width=6, height=8, dpi=300)

quantContig_output
ggplot(data=quantContig_output, aes(x=group, y=contigs, fill=group)) +
  geom_boxplot() +
  labs(x="Group", y = "Number of TCR clonotypes", fill="Group") +
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
ggsave("boxplot_unscaled.pdf", plot=last_plot(), width=6, height=8, dpi=300)

# Export a table:
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = T, exportTable = T)
quantContig_output
abundanceContig(combined, cloneCall = "gene", scale = F)
ggsave("abundance.pdf", plot=last_plot(), width=8, height=6, dpi=300)
lengthContig(combined, cloneCall="nt", chain = "TRA")
# Compare Clonotypes: - run this for all samples.
clonotypeplot <- function(data, sample1, sample2, samplename, height=8, width=6){
  plot <- compareClonotypes(data, 
                    numbers = 10, 
                    samples = c(sample1, sample2), 
                    cloneCall="aa", 
                    graph = "alluvial")
  ggsave(paste0("Clonotype_comparison/",samplename,".pdf"), plot=plot, height=height, width=width, dpi=300)
}
clonotypeplot(combined, "FD1_LUAD","FD1_AISMIA", "FD1")
clonotypeplot(combined, "FD2_LUAD","FD2_AISMIA", "FD2")
clonotypeplot(combined, "FD4_LUAD","FD4_AISMIA", "FD4", 8, 7)
clonotypeplot(combined, "FD5_LUAD","FD5_AISMIA", "FD5")
clonotypeplot(combined, "FD8_LUAD1","FD8_LUAD2", "FD8")
clonotypeplot(combined, "FD9_LUAD","FD9_AISMIA", "FD9")
clonotypeplot(combined, "FD14_LUAD","FD14_AISMIA", "FD14", 8, 7)
clonotypeplot(combined, "FD16_LUAD1","FD16_LUAD2", "FD16")

# Visualize gene usage:
vizGenes(combined, gene="V",chain="TRB",plot="bar",order="variance",scale=T)
vizGenes(combined, gene="V",chain="TRB",plot="heatmap",order="variance",scale=T)

# Clonal space homeostasis:
tb1 <- clonalHomeostasis(combined, cloneCall="gene",
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1),
                  exportTable=T)
ggsave("ClonalSpace.pdf", plot=last_plot(), width=16, height=6, dpi=300)
head(tb1)
tb1 <- as.data.frame(tb1)
for (i in 1:nrow(tb1)){
  tb1[i,"large_and_hyperexpanded"] <- sum(tb1[i,4:5])
}
colnames(tb1)
tb1$group <- c("LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","LUAD","LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","LUAD")
write.table(tb1, file="5d.txt", sep="\t", row.names=T, quote=F)
ggplot(data=tb1, aes(x=group, y=large_and_hyperexpanded, color=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  labs(x="Group", y = "Large and hyperexpanded TCRs", color="Group") +
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position = "none")+
  stat_compare_means()+
  scale_color_brewer(palette="Set1")
ggsave("large_and_hyperexpanded_tcrs.pdf", plot=last_plot(), width=4, height=6, dpi=300)

# Clonal proportion:
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 
ggsave("ClonalProportion.pdf", plot=last_plot(), width=16, height=6, dpi=300)

# Overlap analysis:
clonalOverlap(combined, 
              cloneCall = "gene+nt", 
              method = "morisita")
ggsave("ClonalOverlap.pdf", plot=last_plot(), width=16, height=6, dpi=300)
clonesizeDistribution(combined, 
                      cloneCall = "gene+nt", 
                      method="ward.D2")

# Diversity analysis:

clonalDiversity(example, 
                cloneCall = "gene+nt", 
                group.by = "sample", 
                x.axis = "group", 
                exportTable = F,
                n.boots = 100)+
  stat_compare_means()+
  scale_color_manual(values=c("#377eb8","#e41a1c","#377eb8","#e41a1c","#377eb8","#e41a1c","#377eb8","#e41a1c","#377eb8","#377eb8","#377eb8","#e41a1c","#e41a1c","#377eb8","#377eb8","#377eb8"))
ggsave("ClonalDiversity_group_genent.pdf", plot=last_plot(), width=12, height=6, dpi=300)
#write.table(table_shannon, file="5c.txt", sep="\t", row.names=F, quote=F)
