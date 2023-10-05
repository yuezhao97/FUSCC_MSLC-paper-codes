########### R script for comparing normalized CNV load between groups. ###########
########### Created by Yue Zhao on 2023-04-27. ###########
########### Last modified by Yue Zhao on 2023-04-27. ###########

library(stringr)
library(reshape2)

chr_file <- read.delim("/mnt/sdc/singlecell/scripts/chr_pq.new.txt", sep="\t", header=T)
out <- chr_file[,c("chr","len")]
out <- out[!duplicated(out$chr),]
out$chr <- paste0("chr",out$chr)
write.table(out, file="/mnt/sdc/singlecell/scripts/full_chr_length.txt", sep="\t", quote=F, row.names=F)

samples <- c(1,2,4:17)
wgiis <- c()
for (i in samples){
	print(i)
	cnv_file <- read.delim(paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNVgroup_and_CNVtype_in_T",i,"_complete.txt"), sep="\t", header=T)
	cnv_file_new <- merge(cnv_file, out, by="chr", all.x=F, all.y=F)
	print(head(cnv_file_new))
	cnv_file_new$variant_length <- cnv_file_new$end - cnv_file_new$start
	cnv_file_new$normalized_cnv_load <- cnv_file_new$variant_length / cnv_file_new$len
	wgii <- mean(cnv_file_new$normalized_cnv_load)
	wgiis <- rbind(wgiis, wgii)
}
out_final <- cbind(paste0("T",samples), wgiis)
out_final <- as.data.frame(out_final)
colnames(out_final) <- c("Sample","wgii")
luad <- c("T1","T5","T6","T8","T10","T12","T13","T14","T16","T17")
aismia <- c("T2","T4","T7","T9","T11","T15")
out_final[out_final$Sample %in% luad, "group"] <- "LUAD"
out_final[out_final$Sample %in% aismia, "group"] <- "AIS/MIA"
out_final$wgii <- as.numeric(as.character(out_final$wgii))
write.table(out_final, file="/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/wgii.txt", sep="\t", quote=F, row.names=F)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

ggplot(data=out_final, aes(x=group, y=wgii, fill=group)) +
  geom_boxplot() +
  labs(x="", y = "Weighted Genomic Instability Index\n(wGII)", fill="Group") +
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
ggsave("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/wgii_comparison.pdf", width=6, height=8, dpi=600)
