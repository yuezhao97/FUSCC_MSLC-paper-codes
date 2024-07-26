setwd("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/")

library(tidyverse)

cell_group=read.table("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/T1/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",header = T,sep = "\t",stringsAsFactors = F)
cell_group=cell_group[!str_detect(cell_group$cell_group_name,"references"),]
cell_group$cell_group_name=str_replace(cell_group$cell_group_name,"all.*observations\\.","")
head(cell_group)

group_cellcount=as.data.frame(table(cell_group$cell_group_name))
colnames(group_cellcount)=c("cell_group_name","cellcount")
group_cellcount$cellratio=group_cellcount$cellcount / sum(group_cellcount$cellcount)
head(group_cellcount)

group_cnvtype=read.table("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNVgroup_and_CNVtype_in_T1.txt",header = T,sep = "\t",stringsAsFactors = F)
head(group_cnvtype)
group_cnvtype=group_cnvtype%>%inner_join(group_cellcount,by="cell_group_name")

alltype=c()
for (i in 1:22) {
  for (j in c("p","q")) {
    for (k in c("gain","loss")) {
      alltype=append(alltype,paste("chr",i,j,"_",k,sep = ""))
    }
  }
}
head(alltype)

cellpercent=c()
for (i in alltype) {
  if(i %in% unique(group_cnvtype$cnv_type)){
    tmp=sum(group_cnvtype[group_cnvtype$cnv_type == i,"cellratio"])
    cellpercent=append(cellpercent,tmp)
  }else{
    cellpercent=append(cellpercent,0)
  }
}
names(cellpercent)=alltype
sample0.stat=as.data.frame(cellpercent)

some.sample.stat <- sample0.stat
for (j in c(2,4:17)){
	print(paste0("T",j,"..."))
	cell_group=read.table(paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/T",j,"/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"),header = T,sep = "\t",stringsAsFactors = F)
	cell_group=cell_group[!str_detect(cell_group$cell_group_name,"references"),]
	cell_group$cell_group_name=str_replace(cell_group$cell_group_name,"all.*observations\\.","")
	group_cellcount=as.data.frame(table(cell_group$cell_group_name))
	colnames(group_cellcount)=c("cell_group_name","cellcount")
	group_cellcount$cellratio=group_cellcount$cellcount / sum(group_cellcount$cellcount)
	group_cnvtype=read.table(paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNVgroup_and_CNVtype_in_T",j,".txt"),header = T,sep = "\t",stringsAsFactors = F)
	group_cnvtype=group_cnvtype%>%inner_join(group_cellcount,by="cell_group_name")
	alltype=c()
	for (i in 1:22) {
  		for (s in c("p","q")) {
    			for (k in c("gain","loss")) {
      				alltype=append(alltype,paste("chr",i,s,"_",k,sep = ""))
    			}
  		}
	}
	cellpercent=c()
	for (i in alltype) {
  		if(i %in% unique(group_cnvtype$cnv_type)){
    			tmp=sum(group_cnvtype[group_cnvtype$cnv_type == i,"cellratio"])
    			cellpercent=append(cellpercent,tmp)
  		}else{
    			cellpercent=append(cellpercent,0)
  		}
	}
	names(cellpercent)=alltype
	one.sample.stat=as.data.frame(cellpercent)
	some.sample.stat <- cbind(some.sample.stat, one.sample.stat)
}
some.sample.stat <- as.data.frame(some.sample.stat)
head(some.sample.stat)
colnames(some.sample.stat) <- c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17")
colnames(some.sample.stat) <- c("FD2_LUAD","FD2_MIA","FD1_MIA","FD1_LUAD","FD4_LUAD","FD4_MIA","FD5_LUAD","FD5_MIA","FD8_LUAD1","FD8_LUAD2","FD9_LUAD","FD9_MIA","FD14_LUAD","FD14_AIS","FD16_LUAD1","FD16_LUAD2")
#some.sample.stat <- some.sample.stat[,c("T2","T4","T7","T9","T11","T15","T1","T5","T6","T8","T10","T12","T13","T14","T16","T17")]
write.table(some.sample.stat, file="/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/final_frequency_all.txt", sep="\t", col.names=T, quote=F)

library(RColorBrewer)
library(scales)
library(pheatmap)

samples <- colnames(some.sample.stat)
group <- c("LUAD","AIS/MIA","AIS/MIA","LUAD","LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","LUAD","LUAD","AIS/MIA","LUAD","AIS/MIA","LUAD","LUAD")
anno_col <- as.data.frame(group)
rownames(anno_col) <- samples
anno_col

pheatmap(some.sample.stat,
         cluster_rows = F,cluster_cols = T,
         color = colorRampPalette(brewer.pal(9,"YlGnBu"))(100),
	annotation_col = anno_col,
         #rev(brewer.pal(n = 7, name ="RdYlBu"))
         #brewer.pal(9,"PuRd")
         #brewer.pal(9,"RdPu"),
         border_color = "grey",
         filename = "/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNV.heatmap.pdf",width = 6,height = 12
)
