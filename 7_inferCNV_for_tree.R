########### R script for running inferCNV for individual tumours to generate trees. ###########
########### Created by Yue Zhao on 2022-11-14. ###########
########### Last modified by Yue Zhao on 2022-11-14. ###########

# Load packages:
library(Seurat)
library(infercnv)

setwd("/mnt/sdc/singlecell/data/pre_infercnv")

# Load the samplesheet:
samplesheet <- read.delim("infercnv_for_tree_samplesheet_2ndhalf.txt", sep="\t", header=T)

# Initial settings:
data_path <- "/mnt/sdc/singlecell/data/pre_infercnv"
out_path <- "/mnt/sdc/singlecell/inferCNV/phylogenetic_tree/"

# Run inverCNV:
geneFile="geneFile2.txt"
for (i in 1:nrow(samplesheet)){
	dir_to_make <- paste0(out_path,samplesheet[i,1])
	if (!dir.exists(dir_to_make)){
		mkdir0 <- paste0("mkdir ",dir_to_make)
		system(mkdir0)
	}
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
					tumor_subcluster_partition_method="random_trees",
                        	     cluster_by_groups=FALSE,
					plot_steps=T,
					scale_data=T,
					analysis_mode="subclusters",
					HMM_type="i6",
                             		denoise=TRUE,
                                	output_format="pdf",
                            		 HMM=TRUE)
}
