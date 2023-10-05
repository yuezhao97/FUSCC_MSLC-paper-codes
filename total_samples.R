#setwd("C:\\Users\\lenovo\\Desktop\\GSE138709_RAW")  ##注意工作目录
#BiocManager::install("celldex")

library(celldex)
library(pheatmap)
library(dplyr)
library(patchwork)
library(readr)
library(devtools)
library(Seurat)
library(SeuratData)
library(patchwork)
library(rtracklayer)
library(ggplot2)
library(stringr)
library(paletteer)
library(scales)
library(SingleR)
library(clustree)

data_total<-list()





new_counts1 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD2_LC1_5/filtered_feature_bc_matrix")
new_counts2 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD2_LC2_5/filtered_feature_bc_matrix")
new_counts3 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC1_5/filtered_feature_bc_matrix")
new_counts4 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC2_5/filtered_feature_bc_matrix")
new_counts5 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD1_LC3_5/filtered_feature_bc_matrix")
new_counts6 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD4_LC1_5/filtered_feature_bc_matrix")
new_counts7 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD4_LC2_5/filtered_feature_bc_matrix")
new_counts8 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD5_LC1_5/filtered_feature_bc_matrix")
new_counts9 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD5_LC2_5/filtered_feature_bc_matrix")
new_counts10 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD9_lc1/filtered_feature_bc_matrix")
new_counts11 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD9_lc2/filtered_feature_bc_matrix")
new_counts12 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD8_LC1_5/filtered_feature_bc_matrix")
new_counts13 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD8_LC2_5/filtered_feature_bc_matrix")
new_counts14 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD14_lc1/filtered_feature_bc_matrix")
new_counts15 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD14_lc2/filtered_feature_bc_matrix")

new_counts16 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD16_lc1/filtered_feature_bc_matrix")
new_counts17 <- Read10X(data.dir = "/home/shuishan/data_share/sc_project_data/10XRNAseq/FD16_lc2/filtered_feature_bc_matrix")


data_total[[1]] <- CreateSeuratObject(counts = new_counts1, min.cells = 3,min.features = 200, project = "T1")
data_total[[2]] <- CreateSeuratObject(counts = new_counts2, min.cells = 3,min.features = 200, project = "T2")
data_total[[3]] <- CreateSeuratObject(counts = new_counts3, min.cells = 3,min.features = 200, project = "T3")
data_total[[4]] <- CreateSeuratObject(counts = new_counts4, min.cells = 3,min.features = 200, project = "T4")
data_total[[5]] <- CreateSeuratObject(counts = new_counts5, min.cells = 3,min.features = 200, project = "T5")
data_total[[6]] <- CreateSeuratObject(counts = new_counts6, min.cells = 3,min.features = 200, project = "T6")
data_total[[7]] <- CreateSeuratObject(counts = new_counts7, min.cells = 3,min.features = 200, project = "T7")
data_total[[8]] <- CreateSeuratObject(counts = new_counts8, min.cells = 3,min.features = 200, project = "T8")
data_total[[9]] <- CreateSeuratObject(counts = new_counts9, min.cells = 3,min.features = 200, project = "T9")
data_total[[10]] <- CreateSeuratObject(counts = new_counts10, min.cells = 3,min.features = 200, project = "T10")
data_total[[11]] <- CreateSeuratObject(counts = new_counts11, min.cells = 3,min.features = 200, project = "T11")
data_total[[12]] <- CreateSeuratObject(counts = new_counts12, min.cells = 3,min.features = 200, project = "T12")
data_total[[13]] <- CreateSeuratObject(counts = new_counts13, min.cells = 3,min.features = 200, project = "T13")
data_total[[14]] <- CreateSeuratObject(counts = new_counts14, min.cells = 3,min.features = 200, project = "T14")
data_total[[15]] <- CreateSeuratObject(counts = new_counts15, min.cells = 3,min.features = 200, project = "T15")

data_total[[16]] <- CreateSeuratObject(counts = new_counts16, min.cells = 3,min.features = 200, project = "T16")
data_total[[17]] <- CreateSeuratObject(counts = new_counts17, min.cells = 3,min.features = 200, project = "T17")



#三、质量控???
#Homo_sapiens.GRCh38.103.chr.gtf.gz
#3.1 载入注释
gtf <- import("/home/shuishan/software/genome/Homo_sapiens.GRCh38.104.chr.gtf.gz")
#gtf <- import("D:\\Homo_sapiens.GRCh38.104.chr.gtf.gz")
gtf <- gtf[!is.na(gtf$gene_name)]
gtf <- gtf[gtf$gene_name!=""]

## protein coding
protein <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                    "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                    "nonsense_mediated_decay", "nontranslating_CDS", 
                    "non_stop_decay", "polymorphic_pseudogene", 
                    "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                    "TR_J_gene", "TR_V_gene")]

## mitochondrial genes
mito <- gtf$gene_name[as.character(seqnames(gtf)) %in% "MT"]

## long noncoding
lincRNA <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("3prime_overlapping_ncrna", "ambiguous_orf", 
                    "antisense_RNA", "lincRNA", "ncrna_host", "non_coding", 
                    "processed_transcript", "retained_intron", 
                    "sense_intronic", "sense_overlapping")]

## short noncoding
sncRNA <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("miRNA", "miRNA_pseudogene", "misc_RNA", 
                    "misc_RNA_pseudogene", "Mt_rRNA", "Mt_tRNA", 
                    "Mt_tRNA_pseudogene", "ncRNA", "pre_miRNA", 
                    "RNase_MRP_RNA", "RNase_P_RNA", "rRNA", "rRNA_pseudogene", 
                    "scRNA_pseudogene", "snlRNA", "snoRNA", 
                    "snRNA_pseudogene", "SRP_RNA", "tmRNA", "tRNA",
                    "tRNA_pseudogene", "ribozyme", "scaRNA", "sRNA")]

## pseudogene
pseudogene <- 
  gtf$gene_name[gtf$transcript_biotype %in% 
                  c("disrupted_domain", "IG_C_pseudogene", "IG_J_pseudogene", 
                    "IG_pseudogene", "IG_V_pseudogene", "processed_pseudogene", 
                    "pseudogene", "transcribed_processed_pseudogene",
                    "transcribed_unprocessed_pseudogene", 
                    "translated_processed_pseudogene", 
                    "translated_unprocessed_pseudogene", "TR_J_pseudogene", 
                    "TR_V_pseudogene", "unitary_pseudogene", 
                    "unprocessed_pseudogene")]

annotations <- list(protein=unique(protein), 
                    mito=unique(mito),
                    lincRNA=unique(lincRNA),
                    sncRNA=unique(sncRNA),
                    pseudogene=unique(pseudogene))


names(annotations )
annotations <- annotations[lengths(annotations)>0]

#3.2质控及其细胞选取
# Seurat会计算基因数以及UMI??? (nFeature and nCount)
# Seurat会将原始数据保存在RNA slot中，
# 每一行对应一个基因，每一列对应一个细???.
#DefaultAssay(pbmc) <- "integrated"
#slotNames(pbmc)
#names(pbmc@assays) ## Assays(scRNA)
#head(Idents(pbmc)) ## 查看一下当前的细胞分组

# 在计算比例时，使用目标基因中的数值除以总的数值???
# `[[`运算符可以给metadata加column.



for (i in 1:17){

  data_total[[i]][["percent.mt"]] <- PercentageFeatureSet(data_total[[i]], pattern = "^MT-")
  #作图查看"nFeature_RNA", "nCount_RNA", "percent.mt"三项比例
  # Visualize QC metrics as a violin plot
  plot=VlnPlot(data_total[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #plot
  ggsave(paste0("./result_both/",i,"_nFeature_RNA & nCount_RNA & percent.mt.png"), plot, width=20 ,height=10)



  percent <- lapply(annotations, function(.ele){
    PercentageFeatureSet(data_total[[i]], features = rownames(data_total[[i]])[rownames(data_total[[i]]) %in% .ele])
  })

  # AddMetaData 会在object@meta.data中加一列。这些信息都会在QC中使用到???
  for(j in seq_along(percent)){
    data_total[[i]][[paste0("percent.", names(percent)[j])]] <- percent[[j]]
  }

  #【一维可视化】对上面scRNA 所有注释的信息，进行可视化
  plot1 <- VlnPlot(object = data_total[[i]], 
                   features = c("nFeature_RNA", "nCount_RNA",
                                paste0("percent.", names(percent))), 
                   ncol = 4)
  #plot1
  ggsave(paste0("./result_both/",i,"_annotations_information.png"), plot1, width=20 ,height=20)




  #VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #【二维可视化】选取两个特征进行画图 ??? FeatureScatter
  # FeatureScatter是一个用来画点图的工具，可以比较不同的metadata???
  plot2 <- FeatureScatter(data_total[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- FeatureScatter(data_total[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot4 <- plot2 + plot3
  #plot4
  ggsave(paste0("./result_both/",i,"_correlation.png"), plot4, width=20 ,height=10)

  #开始真正的细胞过滤啦~

  ## 开始过滤
  #计算对应的动态阈值
  Count1<-median(data_total[[i]]$nCount_RNA)-5*mad(data_total[[i]]$nCount_RNA)
  Count2<-median(data_total[[i]]$nCount_RNA)+5*mad(data_total[[i]]$nCount_RNA)
  mt<-median(data_total[[i]]$percent.mt)+6*mad(data_total[[i]]$percent.mt)


  data_total[[i]] <- subset(x = data_total[[i]], 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &percent.protein > 80 & 
                 percent.mt < mt & nCount_RNA>Count1 & nCount_RNA<Count2)
}



#数据整合之前要对每个样本的seurat对象进行数据标准化和选择高变基因
for (i in 1:length(data_total)) {
  data_total[[i]] <- NormalizeData(data_total[[i]])
  data_total[[i]] <- FindVariableFeatures(data_total[[i]], selection.method = "vst")
}


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data_total)
##以VariableFeatures为基础寻找锚点，运行时间较???
scRNA.anchors <- FindIntegrationAnchors(object.list = data_total , anchor.features = features)
##利用锚点整合数据，运行时间较???
data_total <- IntegrateData(anchorset = scRNA.anchors)

save(data_total,file="./result_both/sample_8_data.Rdata")
#load(file="./result_both/sample_data.Rdata")

#DefaultAssay(data_total) <- "RNA"
pbmc<-data_total



#切换回原来的数据集。
DefaultAssay(pbmc) <- "integrated"



#数据降维
all.genes <- rownames(pbmc)
length(all.genes)
#两千基因全部scale进行之后的PCA降维。
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = all.genes)
#如果使用RNA的assay可以考虑使用以下代码进行降维
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))





# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#VizDimLoadings(pbmc, dims = 1:3, reduction = "pca",ncol = 2)

plot4=DimPlot(pbmc,dims = c(1,2), reduction = "pca")

#plot4
ggsave("./result_both/PC1-PC2.png", plot4, width=10 ,height=10)


#plot_list<-list()
#for (i in 1:10) {
#  p=DimHeatmap(pbmc, dims = i, cells = 500, balanced = TRUE)
#  p
#}

plot5=DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE)


ggsave("./result_both/dims1-10.png", plot5, width=10 ,height=10)


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#pbmc <- JackStraw(pbmc, num.replicate = 100)

#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#plot6=JackStrawPlot(pbmc, dims = 1:40)
#plot6
#ggsave("./result_both/JackStrawPlot.png", plot6, width=10 ,height=10)

plot7=ElbowPlot(pbmc)
#plot7
ggsave("./result_both/ElbowPlot.png", plot7, width=10 ,height=10)



pbmc <- FindNeighbors(pbmc, dims = 1:30)
Idents(pbmc)
pbmc <- FindClusters(pbmc, resolution = 0.1)





#确定配色
#colorlist <- hue_pal()(10) #这个是原本的默认配色
#colorlist[c(1,3,5)] <- "gray"  #将几个不想分析的类改成灰???
#UMAP
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)#, cols =colorlist)
p3 <- p1 + p2
#p3
ggsave("./result_both/umap.png", p2, width=10 ,height=10)

#p4 <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident", cols =colorlist)
#p4
#ggsave("0.08resolution_cluster_final/umap_cluster.png", p3, width=12 ,height=7)
#ggsave("0.08resolution_cluster_final/umap_cluster1.png", p4, width=10 ,height=6)

#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20)
p5 <- DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident")
p6 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE)
#p6
p7 <- p5 + p6
#p7
ggsave("./result_both/tsne.png", p6, width=10 ,height=10)

#p8 <- DimPlot(scRNA, reduction = "tsne", split.by = "orig.ident", cols =colorlist)

#ggsave("0.08resolution_cluster_final/tsne_cluster.png", p7, width=20 ,height=10)
#ggsave("0.08resolution_cluster_final/tsne_cluster1.png", p8, width=20 ,height=10)

#保存一下umapandtsne
p9 <- p2 + p6
ggsave("./result_both/umapandtsne_cluster.png", p9, width=12 ,height=6)


saveRDS(pbmc, file = "./result_both/pbmc_final.rds")





# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

#BiocManager::install("cli")

#cluster0
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)




pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
p13=DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#p13
ggsave("./result_both/headmap-top.png", p13, width=10 ,height=12)


write.table(top10,file = './result_both/cluster-difgene.txt',sep = '\t\t')



celltype_marker=c(
  "CD79A","IGKC",#B CELL
  "CD3D","CD3E", #T CELL
  "LYZ","CD68",  #myeloid cells
  "DCN","C1R","COL1A1",#fibroblasts
  "PECAM1","RAMP2","CLDN5",#Endothelial
  "TPSB2","CPA3","MS4A2",#mast cell
  "EPCAM",#cancer cell
  "NKG7",#NK cell
  "AGER",#AT1
  "SFTPA1",#AT2
  "SCGBlAl","CP",#club cell
  "KRTl7","KRT5",#basal airway epithelial
  "TPPP3",#ciliated airway epithelial
  "CTSD","C1QA","APOC1",#macrophage
  "LST1","S100A9","S100A8",#Granulocyte
  "LA-DQA1","LA-DRB1","LA-DQB1"#,DC
)

DefaultAssay(pbmc) <- "RNA"

p13=VlnPlot(pbmc,features = celltype_marker,pt.size = 0,ncol = 6)


ggsave("./result_both/define_cluster.png", p13, width=30 ,height=30)

P14 <- DotPlot(pbmc, features = celltype_marker, cols = c("green","red")) + RotatedAxis()

ggsave("./result_both/Dotplot.png", P14, width=20 ,height=6)




meta=pbmc@meta.data
hpca.se=HumanPrimaryCellAtlasData() ##第一次载入会下载数据集，可能会慢一些，后面在用时就不用下载了
Blue.se=BlueprintEncodeData() 
Immune.se=DatabaseImmuneCellExpressionData()
Nover.se=NovershternHematopoieticData()
MonacoIm.se=MonacoImmuneData()
#hpca.se
#进行singleR注释
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data") ##获取标准化矩阵

pbmcpbmc3.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BLUE=Blue.se, HPCA=hpca.se,Immune=Immune.se,Nover=Nover.se,MonacoIm=MonacoIm.se), 
                          labels = list(Blue.se$label.main, hpca.se$label.main, Immune.se$label.main, Nover.se$label.main, MonacoIm.se$label.main)) 


table(pbmcpbmc3.hesc$labels,meta$seurat_clusters)
m<-table(pbmcpbmc3.hesc$labels,meta$seurat_clusters)
write.table(m,file="./result_both/SingIeR3.txt",sep="\t")
pbmc3<-pbmc
pbmc3@meta.data$labels <-pbmcpbmc3.hesc$labels

p15=DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap")
ggsave("./result_both/SingleR.png", p15, width=30 ,height=30)
saveRDS(pbmc, file = "./result_both/pbmc_SingleR.rds")








meta=pbmc@meta.data
hpca.se <- HumanPrimaryCellAtlasData()

hpca.se
#进行singleR注释
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data") ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc

#seurat 和 singleR的table表
table(pbmc.hesc$labels,meta$seurat_clusters)

pbmc@meta.data$labels <-pbmc.hesc$labels
print(DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap"))

pbmc3<-pbmc

bpe.se<-BlueprintEncodeData()

pbmcpbmc3.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=bpe.se, HPCA=hpca.se), 
                          labels = list(bpe.se$label.main, hpca.se$label.main)) 

table(pbmcpbmc3.hesc$labels,meta$seurat_clusters)

pbmc3@meta.data$labels <-pbmcpbmc3.hesc$labels

print(DimPlot(pbmc3, group.by = c("seurat_clusters", "labels"),reduction = "umap"))
print(plotScoreHeatmap(pbmc.hesc))

#plotDeltaDistribution(pbmc.hesc, ncol = 4)

tab <- table(label = pbmc.hesc$labels,cluster = meta$seurat_clusters)

#pheatmap(log10(tab + 10))


#####=======查看并提取各个类的细???========#######

table(Idents(pbmc))
# fibroblast                      T CELL            endothelial cell                 plasma cell                  macrophage                      B cell               myofibroblast 
# 7300                        7026                        3118                        3014                        1501                        1371                        1101 
# epithelial cell                           8                           9                          10                   mastocyte Plasmacytoid dendritic cell                 melanophore 
# 867                         295                         260                         182                         159                         138                          66 
prop.table(table(Idents(pbmc)))
# fibroblast                      T CELL            endothelial cell                 plasma cell                  macrophage                      B cell               myofibroblast 
# 0.276536101                 0.266156527                 0.118115009                 0.114175316                 0.056860368                 0.051935753                 0.041707705 
# epithelial cell                           8                           9                          10                   mastocyte Plasmacytoid dendritic cell                 melanophore 
# 0.032843397                 0.011175089                 0.009849231                 0.006894462                 0.006023184                 0.005227669                 0.002500189 




#另外一种自动注释软件
#devtools::install_github("ZJUFanLab/scCATCH")
library(scCATCH)

clu_markers <- findmarkergenes(object = pbmc,
                               species = "Human",
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = "Lung",
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)

clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Lung")

new.cluster.ids <- clu_ann$cell_type
names(new.cluster.ids) <- clu_ann$cluster
pbmc <- RenameIdents(pbmc, new.cluster.ids)
p12<-DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("./result_both/scCATCH_cluster.png", p12, width=30 ,height=30)


