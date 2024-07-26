########### R script for running inferCNV for individual patients. ###########
########### Created by Yue Zhao on 2022-11-10. ###########
########### Last modified by Yue Zhao on 2022-11-10. ###########

# Load packages:
library(Seurat)
library(infercnv)

setwd("/mnt/sdc/singlecell/data/pre_infercnv")

# Load the samplesheet:
samplesheet <- read.delim("infercnv_samplesheet.txt", sep="\t", header=T)

# Initial settings:
data_path <- "/mnt/sdc/singlecell/data/pre_infercnv"
out_path <- "/mnt/sdc/singlecell/inferCNV/individual/"

# Run inverCNV:
geneFile="geneFile2.txt"
for (i in 1:nrow(samplesheet)){
	expFile <- paste0(samplesheet[i,2],".txt")
	groupFiles <- paste0(samplesheet[i,3],".txt")
	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("Fib","Endo"))

	infercnv_obj = infercnv::run(infercnv_obj,
        	                     cutoff=0.1,
                	             out_dir=paste0(out_path,samplesheet[i,1]),
                        	     cluster_by_groups=TRUE,
                             		denoise=TRUE,
                                	output_format="pdf",
                            		 HMM=TRUE)
}
#expFile="expfile_FD2.txt"
#groupFiles="groupfile_FD2.txt"
#geneFile="geneFile2.txt"

#infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
#                                   annotations_file=groupFiles,
#                                    delim="\t",
#                                    gene_order_file= geneFile,
#                                    ref_group_names=c("Fib","Endo"))

#infercnv_obj = infercnv::run(infercnv_obj,
#                             cutoff=0.1,
#                             out_dir=out_path,
#                             cluster_by_groups=TRUE,
#                             denoise=TRUE,
#				output_format="pdf",
 #                            HMM=TRUE)


