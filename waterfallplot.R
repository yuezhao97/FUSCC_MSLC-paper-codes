setwd("/Users/zhaoy2/Desktop/sc_project/maf_total")
library(ggplot2)
driver_file <- read.csv("/Users/zhaoy2/Desktop/1000_genomics/20220203_lung_drivers.csv")
driver_genes <- driver_file$Gene_Symbol
head(driver_file)

var_types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Splice_Site","In_Frame_Del","Frame_Shift_Ins","In_Frame_Ins","Nonstop_Mutation","Fusion")
samples <- c("FD2","FD4","FD5","FD8","FD9","FD14","FD16") # Run FD1 separately.
for (s in samples){
  snp1 <- read.delim(paste0(s,"_LC1_SNP.maf"), sep="\t", header=T)
  snp2 <- read.delim(paste0(s,"_LC2_SNP.maf"), sep="\t", header=T)
  indel1 <- read.delim(paste0(s,"_LC1_Indel.maf"), sep="\t", header=T)
  indel2 <- read.delim(paste0(s,"_LC2_Indel.maf"), sep="\t", header=T)
  file_out <- rbind(snp1, snp2, indel1, indel2)
  file_out <- file_out[file_out$Variant_Classification %in% var_types,]
  file_out[,"Joint"] <- paste0(file_out$Hugo_Symbol,".",file_out$aaChange)
  write.table(file_out, file=paste0(s,"_all.txt"), sep="\t", quote=F, row.names=F)
}
s <- "FD1"
snp1 <- read.delim(paste0(s,"_LC3_SNP.maf"), sep="\t", header=T)
snp2 <- read.delim(paste0(s,"_LC2_SNP.maf"), sep="\t", header=T)
indel1 <- read.delim(paste0(s,"_LC3_Indel.maf"), sep="\t", header=T)
indel2 <- read.delim(paste0(s,"_LC2_Indel.maf"), sep="\t", header=T)
file_out <- rbind(snp1, snp2, indel1, indel2)
file_out <- file_out[file_out$Variant_Classification %in% var_types,]
file_out[,"Joint"] <- paste0(file_out$Hugo_Symbol,".",file_out$aaChange)
write.table(file_out, file=paste0(s,"_all.txt"), sep="\t", quote=F, row.names=F)


files_all <- list.files("./", pattern="all.txt")
files_all
combined <- c()
for (f in files_all){
  file1 <- read.delim(f, sep="\t", header=T)
  combined <- rbind(combined, file1)
}

library(maftools)
genes_to_plot <- unique(all.maf$Hugo_Symbol)[unique(all.maf$Hugo_Symbol) %in% driver_genes]
genes_to_plot <- append(genes_to_plot, "ALK")
genes_to_plot
all.maf <- combined
all.clin <- read.delim("../sample_info_noname.txt",sep="\t",header=T)
all.clin <- all.clin[all.clin$Tumor_ID_dedup != "T3",]
head(all.clin)
unique(all.maf$Tumor_Sample_Barcode)
all.maf$Tumor_Sample_Barcode <- factor(all.maf$Tumor_Sample_Barcode, levels=c("FD1_LC3","FD1_LC2","FD2_LC1","FD2_LC2","FD4_LC1","FD4_LC2","FD5_LC1","FD5_LC2","FD8_LC1","FD8_LC2","FD9_LC1","FD9_LC2","FD14_LC1","FD14_LC2","FD16_LC1","FD16_LC2"),
                                       labels=c("FD1_LUAD","FD1_MIA","FD2_LUAD","FD2_MIA","FD4_LUAD","FD4_MIA","FD5_LUAD","FD5_MIA","FD8_LUAD1","FD8_LUAD2","FD9_LUAD","FD9_MIA","FD14_LUAD","FD14_AIS","FD16_LUAD1","FD16_LUAD2"))
all.clin$Tumor_Sample_Barcode <- factor(all.clin$Tumor_Sample_Barcode, levels=c("FD1_LC3","FD1_LC2","FD2_LC1","FD2_LC2","FD4_LC1","FD4_LC2","FD5_LC1","FD5_LC2","FD8_LC1","FD8_LC2","FD9_LC1","FD9_LC2","FD14_LC1","FD14_LC2","FD16_LC1","FD16_LC2"),
                                        labels=c("FD1_LUAD","FD1_MIA","FD2_LUAD","FD2_MIA","FD4_LUAD","FD4_MIA","FD5_LUAD","FD5_MIA","FD8_LUAD1","FD8_LUAD2","FD9_LUAD","FD9_MIA","FD14_LUAD","FD14_AIS","FD16_LUAD1","FD16_LUAD2"))
all.clin$Pathology <- factor(all.clin$Pathology, levels=c("AIS","MIA","LUAD"))
head(all.maf)
all_file <- read.maf(maf=all.maf, clinicalData = all.clin, vc_nonSyn = var_types)

plotmafSummary(maf=all_file, rmOutlier = T, addStat = "median", dashboard=T, titvRaw = F)
col=RColorBrewer::brewer.pal(n=9,name="Paired")
names(col)=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Splice_Site","In_Frame_Del","Frame_Shift_Ins","In_Frame_Ins","Nonstop_Mutation","Fusion")
fabcolors=RColorBrewer::brewer.pal(n=3,name="Set1")
names(fabcolors)=c("AIS","MIA","LUAD")
GGOcolors <- RColorBrewer::brewer.pal(n=3,name="Accent")
names(GGOcolors) <- c("GGO","Subsolid","Solid")

print(GGOcolors)
print(fabcolors)
fabcolors_anno=list(Pathology=fabcolors,Radiological_group=GGOcolors)

colnames(all.clin)
head(all.clin)
anno <- all.clin[,c("Tumor_Sample_Barcode","Pathology")]
anno$Pathology <- factor(anno$Pathology, levels=c("AIS","MIA","LUAD"))

pdf("waterfall.pdf", width=6, height=8)
oncoplot(maf = all_file, 
        colors=col,
        annotationDat = anno,
        clinicalFeatures = "Pathology",
        sortByAnnotation = F,
        annotationOrder = c("AIS","MIA","LUAD"),
        annotationColor = fabcolors_anno,
        genes = genes_to_plot,
        writeMatrix=T,
        removeNonMutated = F,
        showTumorSampleBarcodes=F,
        fontSize=1,top=12,
        keepGeneOrder = F,GeneOrderSort = F)
dev.off()
  #PlotOncogenicPathways(maf=all_file, pathways="RTK-RAS")
