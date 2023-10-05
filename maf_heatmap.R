setwd("/Users/zhaoy2/Desktop/sc_project/maf_total")

s <- "FD2"
snp1 <- read.delim(paste0(s,"_LC1_SNP.maf"), sep="\t", header=T)
snp2 <- read.delim(paste0(s,"_LC2_SNP.maf"), sep="\t", header=T)
indel1 <- read.delim(paste0(s,"_LC1_Indel.maf"), sep="\t", header=T)
indel2 <- read.delim(paste0(s,"_LC2_Indel.maf"), sep="\t", header=T)
dim(indel2)
colnames(snp1)
colnames(indel1)
colnames(snp1) == colnames(indel1)
unique(indel2$Variant_Classification)
unique(snp2$Variant_Classification)

var_types <- c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Del","In_Frame_Ins","Splice_Site","Fusion")
file_fd2 <- rbind(snp1,snp2,indel1,indel2)
file_fd2 <- file_fd2[file_fd2$Variant_Classification %in% var_types,]
nrow(file_fd2)
head(file_fd2)
file_fd2_brief <- file_fd2[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
head(file_fd2_brief)
library(reshape2)
file_fd2_brief_m <- reshape(file_fd2_brief, direction="wide", idvar="Hugo_Symbol", timevar="Tumor_Sample_Barcode")
head(file_fd2_brief_m)

######




samples <- c("FD2","FD4","FD5","FD8","FD9","FD14","FD16") # Run FD1 separately.
for (s in samples){
  snp1 <- read.delim(paste0(s,"_LC1_SNP.maf"), sep="\t", header=T)
  snp2 <- read.delim(paste0(s,"_LC2_SNP.maf"), sep="\t", header=T)
  indel1 <- read.delim(paste0(s,"_LC1_Indel.maf"), sep="\t", header=T)
  indel2 <- read.delim(paste0(s,"_LC2_Indel.maf"), sep="\t", header=T)
  file_out <- rbind(snp1, snp2, indel1, indel2)
  file_out <- file_out[file_out$Variant_Classification %in% var_types,]
  file_out_brief <- file_out[,c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","aaChange")]
  file_out_brief[,"Joint"] <- paste0(file_out_brief$Hugo_Symbol,".",file_out_brief$aaChange)
  write.table(file_out_brief, file=paste0(s,"_all.txt"), sep="\t", quote=F, row.names=F)
}
s <- "FD1"
snp1 <- read.delim(paste0(s,"_LC3_SNP.maf"), sep="\t", header=T)
snp2 <- read.delim(paste0(s,"_LC2_SNP.maf"), sep="\t", header=T)
indel1 <- read.delim(paste0(s,"_LC3_Indel.maf"), sep="\t", header=T)
indel2 <- read.delim(paste0(s,"_LC2_Indel.maf"), sep="\t", header=T)
file_out <- rbind(snp1, snp2, indel1, indel2)
file_out <- file_out[file_out$Variant_Classification %in% var_types,]
file_out_brief <- file_out[,c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","aaChange")]
file_out_brief[,"Joint"] <- paste0(file_out_brief$Hugo_Symbol,".",file_out_brief$aaChange)
write.table(file_out_brief, file=paste0(s,"_all.txt"), sep="\t", quote=F, row.names=F)

library(pheatmap)
library(RColorBrewer)

fd2_heatmap <- read.delim("FD2_for_heatmap.txt", sep="\t", header=T, row.names=1)
fd2_heatmap_t <- t(fd2_heatmap)

pheatmap(fd2_heatmap_t, cluster_rows = F, cluster_cols = F,
         filename = "FD2_heatmap.pdf",width = 8,height = 0.5,
         color = c("white","navyblue","maroon"),
         show_colnames = F, legend=F,
         #border_color = NA
         #color = colorRampPalette(brewer.pal(9,"PuRd"))(100)
         )

samples_for_heatmap <- c("FD1","FD2","FD4","FD5","FD8","FD9","FD14","FD16")
for (s in samples_for_heatmap){
  print(s)
  file1 <- read.delim(paste0(s,"_for_heatmap.txt"), sep="\t", header=T, row.names=1)
  file1_t <- t(file1)
  pheatmap(file1_t, cluster_rows = F, cluster_cols = F,
           filename = paste0(s,"_heatmap.pdf"),width = 8,height = 0.5,
           color = c("white","navyblue","maroon"),
           show_colnames = F, legend=F,
           #border_color = NA
           #color = colorRampPalette(brewer.pal(9,"PuRd"))(100)
  )
}
file1 <- read.delim(paste0(s,"_for_heatmap.txt"), sep="\t", header=T)
head(file1)
duplicated(file1$Gene)
file1$Gene

### Re-organize the table for heatmap:
samples_for_heatmap <- c("FD2","FD4","FD5","FD8","FD9","FD14") # run FD1/FD16 separately.
for (s in samples_for_heatmap){
  print(s)
  file1 <- read.delim(paste0(s,"_all.txt"), sep="\t", header=T)
  genes <- unique(file1$Hugo_Symbol)
  out <- data.frame(Gene=genes)
  s1 <- file1[file1$Tumor_Sample_Barcode == paste0(s,"_LC1"),]
  s2 <- file1[file1$Tumor_Sample_Barcode == paste0(s,"_LC2"),]
  common_file <- s1[s1$Joint %in% s2$Joint,] # common
  parallel_file <- s1[s1$Hugo_Symbol %in% s2$Hugo_Symbol & !(s1$Joint %in% s2$Joint),] # parallel
  out[out$Gene %in% common_file$Hugo_Symbol, "Class"] <- "Common"
  out[out$Gene %in% parallel_file$Hugo_Symbol, "Class"] <- "Parallel_evolution"
  out[!(out$Class %in% c("Common","Parallel_evolution")), "Class"] <- "Others"

  out[out$Gene %in% s1$Hugo_Symbol & out$Class == "Others", "LC1"] <- 1
  out[out$Class == "Parallel_evolution", "LC1"] <- 2
  out[out$Class == "Common", "LC1"] <- 3
  out[!(out$LC1 %in% c(1,2,3)), "LC1"] <- 0
  
  out[out$Gene %in% s2$Hugo_Symbol & out$Class == "Others", "LC2"] <- 1
  out[out$Class == "Parallel_evolution", "LC2"] <- 2
  out[out$Class == "Common", "LC2"] <- 3
  out[!(out$LC2 %in% c(1,2,3)), "LC2"] <- 0
  
  write.table(out, file=paste0(s,"_for_heatmap_new.txt"), sep="\t", quote=F, row.names=F)
}

for (s in samples_for_heatmap){
  out <- read.delim(paste0(s,"_for_heatmap_new.txt"), sep="\t", header=T,row.names = 1)
  out <- out[,-1]
  n <- nrow(out)
  out_t <- t(out)
  out_t <- as.data.frame(out_t)
  pheatmap(out_t, cluster_rows = F, cluster_cols = F,
           filename = paste0(s,"_heatmap.pdf"),width = 2.75,height = 0.5,
           color = c("grey","navyblue","orange","maroon"),
           show_colnames = F, legend=F,
           #border_color = NA
           #color = colorRampPalette(brewer.pal(9,"PuRd"))(100)
  )
}

# FD1:
s <- "FD1"
file1 <- read.delim(paste0(s,"_all.txt"), sep="\t", header=T)
genes <- unique(file1$Hugo_Symbol)
out <- data.frame(Gene=genes)
s1 <- file1[file1$Tumor_Sample_Barcode == paste0(s,"_LC3"),]
s2 <- file1[file1$Tumor_Sample_Barcode == paste0(s,"_LC2"),]
common_file <- s1[s1$Joint %in% s2$Joint,] # common
parallel_file <- s1[s1$Hugo_Symbol %in% s2$Hugo_Symbol & !(s1$Joint %in% s2$Joint),] # parallel
out[out$Gene %in% common_file$Hugo_Symbol, "Class"] <- "Common"
out[out$Gene %in% parallel_file$Hugo_Symbol, "Class"] <- "Parallel_evolution"
out[!(out$Class %in% c("Common","Parallel_evolution")), "Class"] <- "Others"

out[out$Gene %in% s1$Hugo_Symbol & out$Class == "Others", "LC3"] <- 1
out[out$Class == "Parallel_evolution", "LC3"] <- 2
out[out$Class == "Common", "LC3"] <- 3
out[!(out$LC3 %in% c(1,2,3)), "LC3"] <- 0

out[out$Gene %in% s2$Hugo_Symbol & out$Class == "Others", "LC2"] <- 1
out[out$Class == "Parallel_evolution", "LC2"] <- 2
out[out$Class == "Common", "LC2"] <- 3
out[!(out$LC2 %in% c(1,2,3)), "LC2"] <- 0
out
write.table(out, file=paste0(s,"_for_heatmap_new.txt"), sep="\t", quote=F, row.names=F)

out <- read.delim(paste0(s,"_for_heatmap_new.txt"), sep="\t", header=T,row.names = 1)
out <- out[,-1]
n <- nrow(out)
out_t <- t(out)
out_t <- as.data.frame(out_t)
pheatmap(out_t, cluster_rows = F, cluster_cols = F,
         filename = paste0(s,"_heatmap.pdf"),width = 2.75,height = 0.5,
         color = c("grey","navyblue","orange","maroon"),
         show_colnames = F, legend=F,
         #border_color = NA
         #color = colorRampPalette(brewer.pal(9,"PuRd"))(100)
)

# FD16:
s <- "FD16"
out <- read.delim(paste0(s,"_for_heatmap_new.txt"), sep="\t", header=T,row.names = 1)
out <- out[,-1]
n <- nrow(out)
out_t <- t(out)
out_t <- as.data.frame(out_t)
pheatmap(out_t, cluster_rows = F, cluster_cols = F,
         filename = paste0(s,"_heatmap.pdf"),width = 2.75,height = 0.5,
         color = c("grey","navyblue","orange","maroon"),
         show_colnames = F, legend=F,
         #border_color = NA
         #color = colorRampPalette(brewer.pal(9,"PuRd"))(100)
)
