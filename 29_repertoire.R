########### Serial analyses using scRepertoire. ###########
########### Created by Yue Zhao on 2023-03-15. ###########
########### Last modified by Yue Zhao on 2023-03-15. ###########

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
	samples = c("FD1","FD1","FD2","FD2","FD4","FD4","FD5","FD5","FD8","FD8","FD9","FD9","FD14","FD14","FD16","FD16"),
	ID = c("LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD","AISMIA","LUAD1","LUAD2","LUAD","AISMIA","LUAD","AISMIA","LUAD1","LUAD2"),
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

ggplot(data=quantContig_output, aes(x=group, y=scaled, color=group)) +
  geom_boxplot() +
  geom_jitter(size=1, width=0.5)+
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
clonalHomeostasis(combined, cloneCall="gene",
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))
ggsave("ClonalSpace.pdf", plot=last_plot(), width=16, height=6, dpi=300)
clonalHomeostasis(combined, cloneCall="aa")

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
                n.boots = 100)+
  stat_compare_means()
ggsave("ClonalDiversity_group_genent.pdf", plot=last_plot(), width=12, height=6, dpi=300)

