setwd("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/")

library(tidyverse)
chr_pq=read.table("/mnt/sdc/singlecell/scripts/chr_pq.new.txt",header = T,sep = "\t",stringsAsFactors = F) 
chr_pq$chr=paste("chr",chr_pq$chr,sep = "")
head(chr_pq)

for (j in c(1:17)){
	cnv_regions=read.table(paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/T",j,"/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"),header = T,sep = "\t",stringsAsFactors = F)
	cnv_regions=cnv_regions%>%filter(state!=3)
	cnv_regions$cell_group_name=str_replace(cnv_regions$cell_group_name,"all.*observations\\.","")
	cnv_regions$cnv_type=""
	for (i in 1:nrow(cnv_regions)) {
  		tmp_chr_pq=chr_pq%>%filter(chr==cnv_regions[i,"chr"])
  		if(cnv_regions[i,"start"] <= tmp_chr_pq[1,"cutoff"] & cnv_regions[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
    		cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
  		}else if (cnv_regions[i,"start"] >= tmp_chr_pq[2,"cutoff"] & cnv_regions[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
    		cnv_regions[i,"cnv_type"]=paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
  		}else {
    		cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"chr"],"p,",cnv_regions[i,"chr"],"q",sep = "")
  		}
  		if (cnv_regions[i,"state"] < 3) {
    		cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_loss",sep = "")
  		}
  		if (cnv_regions[i,"state"] > 3) {
    		cnv_regions[i,"cnv_type"]=paste(cnv_regions[i,"cnv_type"],"_gain",sep = "")
  		}
  		if (str_detect(cnv_regions[i,"cnv_type"],",")) {
    			tmp1=strsplit(cnv_regions[i,"cnv_type"],",")[[1]][1]
    			tmp2=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
    			tmp3=strsplit(strsplit(cnv_regions[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
    			cnv_regions[i,"cnv_type"]=paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
    			rm(list = c("tmp1","tmp2","tmp3"))
  		}
	}
	write.table(cnv_regions,file = paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNVgroup_and_CNVtype_in_T",j,"_complete.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
	cnv_regions_copy=cnv_regions
	cnv_regions_part1=cnv_regions[str_detect(cnv_regions$cnv_type,","),]
	cnv_regions_part2=cnv_regions[!str_detect(cnv_regions$cnv_type,","),]
	cnv_regions_part1_new=as.data.frame(matrix(NA,ncol = ncol(cnv_regions_part1), nrow = nrow(cnv_regions_part1)*2 ))
	colnames(cnv_regions_part1_new)=colnames(cnv_regions_part1)
	for (i in 1:dim(cnv_regions_part1)[1]) {
  		cnv_regions_part1_new[2*i-1,]=cnv_regions_part1[i,]
  		cnv_regions_part1_new[2*i,]  =cnv_regions_part1[i,]

  		tmp_chr_pq=chr_pq[chr_pq$chr == cnv_regions_part1[i,"chr"],]
  		cnv_regions_part1_new[2*i-1,"end"] = tmp_chr_pq[1,"cutoff"]
  		cnv_regions_part1_new[2*i,"start"] = tmp_chr_pq[2,"cutoff"]

  		cnv_regions_part1_new[2*i-1,"cnv_type"] = strsplit(cnv_regions_part1_new[2*i-1,"cnv_type"],",")[[1]][1]
  		cnv_regions_part1_new[2*i,"cnv_type"]   = strsplit(cnv_regions_part1_new[2*i,  "cnv_type"],",")[[1]][2]
	}
	cnv_regions=cnv_regions_part2 %>% rbind(cnv_regions_part1_new)
	cnv_regions$chr=factor(cnv_regions$chr,levels = paste("chr",1:22,sep = ""))
	cnv_regions=cnv_regions%>%arrange(cell_group_name,chr,start)
	cnv_regions$region_len=cnv_regions$end-cnv_regions$start+1
	cnv_regions=cnv_regions[,c("cell_group_name","cnv_type","region_len")]
	cnv_regions=as.data.frame(cnv_regions%>%group_by(cell_group_name,cnv_type)%>%dplyr::summarize(all_region_len=sum(region_len)))
	cnv_regions$final_label=""
	ratio=0.1
	for (i in 1:dim(cnv_regions)[1]) {
  		arm=strsplit(cnv_regions[i,"cnv_type"],"_")[[1]][1]
  		ref_len=chr_pq[paste(chr_pq$chr,chr_pq$arm,sep = "") == arm,"arm_len"]

  		if(cnv_regions[i,"all_region_len"] >= ref_len*ratio){
    			cnv_regions[i,"final_label"]=cnv_regions[i,"cnv_type"]
  		}else{
    			cnv_regions[i,"final_label"]=""
  		}
	}
	cnv_regions$cnv_type=NULL
	colnames(cnv_regions)[3]="cnv_type"
	cnv_regions=cnv_regions%>%filter(cnv_type != "")
	write.table(cnv_regions[,c("cell_group_name","cnv_type")],file = paste0("/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/frequency_heatmap/CNVgroup_and_CNVtype_in_T",j,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
}

