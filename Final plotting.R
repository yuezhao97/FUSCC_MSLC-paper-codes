

#绘制TCR的top10相关图
library(ggplot2)
library(reshape2)
setwd("/Users/wangjun/Desktop/最终方案/figure/TCR/")


data<-read.table(file="./FD1_LUAD_MIA1_top_10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))

ggsave(p1,file="./FD1_LUAD_MIA.png",dpi = 1000,height = 4.52,width = 14.9)





data<-read.table(file="./FD1_LUAD_MIA2_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD1_LUAD_MIA2.png",dpi = 1000)
ggsave(p1,file="./FD1_LUAD_MIA2.png",dpi = 1000,width = 6.86,height =4.09)









data<-read.table(file="./FD2_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD2_LUAD_MIA.png",dpi = 1000)
ggsave(p1,file="./FD2_LUAD_MIA.png",dpi = 1000,width = 7.62,height =4.05)







data<-read.table(file="./FD4_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD4_LUAD_MIA.png",dpi = 1000)
ggsave(p1,file="./FD4_LUAD_MIA.png",dpi = 1000,width = 6.36,height =4.03)





data<-read.table(file="./FD5_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD5_LUAD_MIA.png",dpi = 100)
ggsave(p1,file="./FD5_LUAD_MIA.png",dpi = 1000,width = 6.5,height =3.55)





data<-read.table(file="./FD8_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD8_LUAD_MIA.png",dpi = 100)
ggsave(p1,file="./FD8_LUAD_MIA.png",dpi = 1000,width = 4.92,height =3.42)








data<-read.table(file="./FD9_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD9_LUAD_MIA.png",dpi = 100)
ggsave(p1,file="./FD9_LUAD_MIA.png",dpi = 1000,width = 6.75,height =3.36)





data<-read.table(file="./FD14_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD14_LUAD_MIA.png",dpi = 100)
ggsave(p1,file="./FD14_LUAD_MIA.png",dpi = 1000,width = 5.18,height =3.45)







data<-read.table(file="./FD16_LUAD_MIA_top10.csv",header = T,sep=",")
data<-melt(data)
colnames(data)<-c("Sequence","Tissue","Proportion")
p1<-ggplot(data, aes(
  x = factor(Sequence,levels = unique(Sequence)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Proportion, -Proportion),  # 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(-0.2,0.2, 0.01),
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./FD16_LUAD_MIA.png",dpi = 100)
ggsave(p1,file="./FD16_LUAD_MIA.png",dpi = 1000,width = 6.71,height =3.7)







#探究纤维细胞在两组的分布是否存在差异
library(Seurat)
setwd("/Users/wangjun/Desktop/最终方案/figure/")

fib<-readRDS(file="./Fibroblast.rds")

table(fib$orig.ident)

fib_mia<-c(0.055989583,0.065934066,0.164503245,0.052105423,0.01369484,0.039802484,0.094070696)
fib_luad<-c(0.112368024,0.112368024,0.246222633,0.032636739,0.031098825,0.057657658,0.022685789)

wilcox.test(fib_luad,fib_mia,paired = T)



#绘制scMetabolism的图

tumor<-readRDS(file="./tumor_infercnv.rds")
Idents(tumor)<-"orig.ident"

LUAD<-colnames(subset(tumor,idents = c("T1","T5","T6","T8","T10","T12","T13","T14","T16","T17")))
MIA<-colnames(subset(tumor,idents = c("T2","T3","T4","T7","T9","T11","T15")))

score<-read.table(file="./scMetabolism_tumor_KEGG.csv",header = T,row.names = 1,sep = ",")
dim(score)

dim(tumor)

MIA<-gsub("-",".",MIA)
LUAD<-gsub("-",".",LUAD)

#计算effect size

library(pwr)
pwr.t.test(n=16008,d=,sig.level = 0.05,alternative = "two.sided")

pathway<-rownames(score)
score<-as.matrix(score)
score[1:4,1:4]
num<-c()
for (i in 1:81){
  d_score=abs(mean(score[i,LUAD])-mean(score[i,MIA]))/sd(score)
  power<-pwr.t2n.test(n1=5495 , n2= 10513,d=d_score,power = NULL, sig.level =0.01, alternative = "two.sided" )
  #power<-pwr.t.test(n=16008,d=d_score,power = NULL, sig.level =0.01, alternative = "two.sided" )
  if (power$power>0.99){
    print(pathway[i])
    num<-c(num,i)
    
  }
}
mean(score[1,LUAD])





power<-pwr.t2n.test(n1=5495 , n2= 10513, d=0.05,power = NULL,sig.level =0.05, alternative = "two.sided" )

power$power



score[1:4,1:4]
LUAD[1:4]

#作图
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


score_t<-as.data.frame(t(score))
score_t$Tissue<-"MIA"
score_t[LUAD,]$Tissue<-"LUAD"
table(score_t$Tissue)

score_t$Tissue<-factor(score_t$Tissue,levels = c("MIA","LUAD"))

for (i in 1:81){
  d_score=abs(mean(score[i,LUAD])-mean(score[i,MIA]))/sd(score)
  power<-pwr.t2n.test(n1=5495 , n2= 10513,d=d_score,power = NULL, sig.level =0.01, alternative = "two.sided" )
  #power<-pwr.t.test(n=16008,d=d_score,power = NULL, sig.level =0.01, alternative = "two.sided" )
  if (power$power>0.99){
    print(pathway[i])
    p<-ggplot(score_t,aes(x=`Tissue`,y=score_t[,i],color=`Tissue`))+geom_boxplot()
    
    p1<-p+theme_bw()+theme(panel.grid = element_blank())+scale_color_manual(values=c("#4DBBD5FF","#E64B35FF"))+stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('MIA','LUAD')),method="wilcox.test")
    p2<-p1+theme(axis.text.x=element_text(hjust = 0.5,size=12,face = "bold",color = "black"),
                 axis.text.y=element_text(size=12,face = "bold",color = "black"),
                 axis.title.y=element_text(size=12,face = "bold",color = "black"),
                 axis.title.x = element_text(size=12,face = "bold",color = "black"),
                 plot.title=element_text(hjust=0.5, size=12,face = "bold",color = "black"))+    
      labs(y=colnames(score_t)[i],x="Tissue",title="Boxplot")+   #添加标签
      #ylim(0,0.061)+
      theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))#去除网格
    p2
    filenames=paste("./scMetabolism/total/",colnames(score_t)[i],".png",sep="")
    print(filenames)
    filenames=gsub(" / ","_",filenames)
    ggsave(p2,file=filenames,dpi=1000,width = 2.97 ,height =  6.25)
    
  }
}

show_col(pal[c(2,9)])
#4DBBD5FF,#E64B35FF



#
MIA_score<-rowMeans(score[,MIA])
LUAD_score<-rowMeans(score[,LUAD])

MIA_score
length(MIA_score)
mean_data<-data.frame(matrix(NA,81,2))
colnames(mean_data)<-c("MIA","LUAD")

mean_data$MIA<-MIA_score
mean_data$LUAD<-LUAD_score

rownames(mean_data)<-pathway

mean_data[1:4,]

mean_data<-mean_data[num,]



data<-melt(t(mean_data))

data[1:4,]

min(data$Score)
data$Score<-data$Score+0.15




colnames(data)<-c("Tissue","Pathway","Score")
data$Score<-((data$Score-min(data$Score))/(max(data$Score)-min(data$Score)))
data$Tissue<-factor(data$Tissue,levels = c("LUAD","MIA"))
data$Tissue
p1<-ggplot(data, aes(
  x = factor(Pathway,levels = unique(Pathway)),             # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
  y = ifelse(Tissue == "LUAD", Score, -Score),# 判断分组情况，将两个柱子画在0的两侧
  fill = Tissue)) +
  geom_bar(stat = 'identity')+                                # 画柱形图
  coord_flip()+                                               # x轴与y轴互换位置                                                   # 标签大小
  scale_y_continuous(                                        # 调整y轴
    labels = abs,                                             # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
p1
ggsave(p1,file="./Score.png",dpi = 100)
ggsave(p1,file="./Score.png",dpi = 1000,width = 7,height =10)



mean_data[c("Citrate cycle (TCA cycle)"),]
pathway

data[1:4,]



dim(mean_data)

mean_data$MIA<-(mean_data$MIA-min(mean_data))/(max(mean_data)-min(mean_data))
mean_data$LUAD<-(mean_data$LUAD-min(mean_data))/(max(mean_data)-min(mean_data))

new_data<-mean_data[(mean_data$LUAD >0.1 | mean_data$MIA>0.1),]

new_data$FC<-new_data$LUAD/new_data$MIA

rownames(new_data)

new_data$col<-ifelse(new_data$FC>1,"#E64B35FF","#4DBBD5FF")

table(new_data$col)
new_data$col<-factor(new_data$col,levels = c("#4DBBD5FF","#E64B35FF"))

p2<-ggplot(new_data, aes(x = factor(rownames(new_data),levels = unique(rownames(new_data))), y = FC))+
  geom_point(size=ifelse(new_data$FC > (1/new_data$FC),new_data$FC*3+3,(1/new_data$FC)*3+3),color=new_data$col)+
  coord_flip()+
  scale_y_continuous(                                        # 调整y轴
    breaks=seq(0,2.5, 0.5),                                          # 刻度设置为绝对值
    expand = expansion(mult = c(0.1, 0.1))) +
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))

ggsave(p2,file="./Score_plot.png",dpi = 1000,width = 6,height =8)

pathway


write.table(mean_data,file="./Score_plot_data.csv",col.names = T,row.names = T,sep = ",")

?write.table








library(monocle)

load(file="/Users/wangjun/Desktop/monocle/CDS_1500.Rdata")
#load(file="/Users/wangjun/Desktop/monocle/CDS_2000.Rdata")
load(file="/Users/wangjun/Desktop/monocle/CDS_1000.Rdata")


load(file="/Users/wangjun/Desktop/monocle/CDS_2000.Rdata")





table(CDS$orig.ident)

library(paletteer) 
library(scales)
pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal_new<-(colorRampPalette(pal)(13))
show(pal)


p3<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+facet_wrap("~Type", nrow = 1)+
  scale_color_manual(values=pal[c(1,5,2)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p3


ggsave(p3,file="./monocle_Type.png")
ggsave(p3,file="/Users/wangjun/Desktop//monocle_Type.png",dpi = 1000,width = 5.02,height =3.43)



pal[c(9,2,7)]

pal[1:10]
pal[c(9,2,7)]
#00A087FF #4DBBD5FF #B09C85FF #91D1C2FF #3C5488FF #8491B4FF #F39B7FFF #7E6148FF #E64B35FF 

plot_cell_trajectory(CDS, color_by = "Type", cell_size = 1)+scale_color_manual(values=c("#0F99B2","#E36666","#F2B77C","#66CCFF"))+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

CDS$Type

p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p2<-plot_cell_trajectory(CDS, color_by = "Pseudotime", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

p1
p2

table(CDS$orig.ident)


CDS


CDS <- orderCells(CDS)

?orderCells



#查看各个marker的表达
pData(CDS)$ADH1C = log2(exprs(CDS)['ADH1C',]+1)
p2=plot_cell_trajectory(CDS, color_by = "ADH1C" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2



pData(CDS)$CAPN13 = log2(exprs(CDS)['CAPN13',]+1)
p2=plot_cell_trajectory(CDS, color_by = "CAPN13" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2


pData(CDS)$DDIT4 = log2(exprs(CDS)['DDIT4',]+1)
p2=plot_cell_trajectory(CDS, color_by = "DDIT4" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2



pData(CDS)$LDHA = log2(exprs(CDS)['LDHA',]+1)
p2=plot_cell_trajectory(CDS, color_by = "LDHA" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2


pData(CDS)$MGP = log2(exprs(CDS)['MGP',]+1)
p2=plot_cell_trajectory(CDS, color_by = "MGP" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2

# 将特定特征值作为细胞的属性
CDS_order <- order(-CDS$MGP)

# 设置排序后的细胞顺序
CDS <- CDS[CDS_order, ]


CDS$cell_order


# 绘制拟时序结果的坐标并对细胞排序
plot_cell_trajectory(CDS, color_by = "MGP", reorder_cells = TRUE)

CDS_test <- orderCells(CDS, gene = "MGP")

library(Seurat)
CDS_test<-orderCells(monocle_result, cell_ids = CDS$cell_order)
?reducedDim

??reducedDim



# 绘制轨迹图，使用 plot_cells 函数
# 将 plot_cells 函数中的 cells 参数设置为排序后的细胞顺序，以控制图层顺序
plot_cells(CDS, cells = CDS_order)



plot_cell_trajectory(CDS, cells = CDS_order,color_by = "MGP")


p2=plot_cell_trajectory(CDS,cells = CDS_order, color_by = "MGP" ,cell_size = 1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2

?plot_cell_trajectory



# 按照特定特征值对细胞进行排序
cell_order <- order(CDS$MGP)

# 将单细胞转录组数据转换为DataFrame
sce_df <- as.data.frame(t(assays(CDS)$counts))

sce_df$cell_feature <- cell_feature
sce_df$cell_order <- cell_order

# 使用ggplot2包绘制单细胞转录组轨迹图
ggplot(sce_df, aes(x = cell_order, y = CellName)) +
  geom_point(aes(size = cell_feature), color = "blue", alpha = 0.7) +
  scale_size_continuous(range = c(1, 5)) +
  theme_minimal() +
  labs(x = "Cell Order", y = "Cell Name") +
  coord_cartesian(clip = "off")










CDS$cell_order <- CDS_order
plot_cell_trajectory(CDS, color_by = "MGP",cell_order_attribute = "cell_order")





#提取CDS的坐标重新绘制
library(monocle)
library(ggplot2)

# 假设你已经运行了 Monocle 的主成分分析 (PCA) 或者 t-distributed stochastic neighbor embedding (t-SNE) 等降维算法，并将结果保存在 Monocle 对象中

# 提取降维后的坐标
reduced_dim <- reducedDim(CDS, "DDRTree")  # 假设使用了主成分分析 (PCA) 降维
coordinates <- reduced_dim$cell_embeddings  # 提取降维后的细胞坐标



summary(CDS@reducedDimS)
summary(CDS@reducedDimK)

dim(CDS@reducedDimS)

CDS@reducedDimS[1:2,1:2]




Gene_list<-c("ADH1C", "CAPN13", "DDIT4", "LDHA", "MGP", "MUC5B", "PCP4L1")


for (i in 1:length(Gene_list)){
#for (i in 1:1){
  pData(CDS)$score = log2(exprs(CDS)[Gene_list[i],]+1)
  score<-CDS$score
  data <- data.frame(x = CDS@reducedDimS[1,], y = CDS@reducedDimS[2,], score = score)
  data<-data[order(data$score), ]
  colnames(data)<-c("Component_1","Component_2","score")
  p1<-ggplot(data, aes(x = Component_1, y = Component_2,color=score)) +
    geom_point(size = 1) +  # 绘制细胞点
    #geom_point(position = position_jitter(height = 0.3)) +
    scale_color_gradient(paste0(Gene_list[i],"_Exp"),low = "#F8DBDA", high = "#D31E18")+   #修改颜色
    #scale_color_viridis()+#自定义颜色
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.box.background=element_rect(colour = "white",fill = "white"))
  
  ggsave(p1,file=paste0("/Users/wangjun/Desktop/图片补充/",Gene_list[i],"_Exp.pdf"),dpi = 500)
  
}














pData(CDS)$MUC5B = log2(exprs(CDS)['MUC5B',]+1)
p2=plot_cell_trajectory(CDS, color_by = "MUC5B" ,cell_size = 1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2



pData(CDS)$PCP4L1 = log2(exprs(CDS)['PCP4L1',]+1)
p2=plot_cell_trajectory(CDS, color_by = "PCP4L1" ,cell_size = .1)  + scale_colour_gradient(low = "#6666CC", high = "#FF9900")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
p2












p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+
  scale_color_manual(values=pal[c(1,5,2)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p2<-plot_cell_trajectory(CDS, color_by = "Pseudotime", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

p4<-p1+p2
p4
ggsave(p4,file="/Users/wangjun/Desktop/monocle_Pseudotime.png")
ggsave(p4,file="/Users/wangjun/Desktop/monocle_Pseudotime_4.png",dpi = 1000,width = 6.24,height =3.04)

p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+
  scale_color_manual(values=pal[c(1,5,2)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p2<-plot_cell_trajectory(CDS, color_by = "Pseudotime", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

p3<-p1+p2
p3
ggsave(p4,file="/Users/wangjun/Desktop/monocle_Pseudotime.png")
ggsave(p3,file="/Users/wangjun/Desktop/monocle_Pseudotime_3.png",dpi = 1000,width = 6.24,height =3.04)






p3<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5,show_branch_points = F)+facet_wrap("~orig.ident", nrow = 1)+
  scale_color_manual(values=pal[c(1,5,2)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title.x=element_text(size=5,face = "bold"),axis.title.y=element_text(size=5,face = "bold"),axis.text.x=element_text(size=5,face = "bold",color = "black"),axis.text.y=element_text(size=5,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 5))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p3

ggsave(p3,file="./monocle_Pseudotime_orig.ident.png",dpi = 1000,width = 10.5,height =2.18)



table(CDS$Type[CDS$orig.ident!="T3" & CDS$State==3 ])

CDS$State


table(tumor)

colnames(CDS)[CDS$State==1]

tumor<-AddMetaData(tumor,CDS$State,col.name = "CDS_result")

Idents(tumor)<-"CDS_result"

table(Idents(tumor))

Idents(tumor)<-"orig.ident"
table(Idents(tumor))

tumor<-subset(tumor,idents="T3",invert=T)

table(tumor$CDS_result)

Idents(tumor)<-"CDS_result"


new.cluster.ids <- c("1","other","other")

levels(tumor)
names(new.cluster.ids) <- levels(tumor)
new.cluster.ids
tumor <- RenameIdents(tumor, new.cluster.ids)



table(Idents(pbmc))




library(dplyr)

scRNA.markers <- FindAllMarkers(object = tumor,
                                min.pct = 0.1, thresh.use = 0)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
marker_gene_2 <- scRNA.markers %>% filter(p_val_adj <= 0.05 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

dim(marker_gene)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/tumor_CDS_DEG_gene_update_3_17_new.csv")

?FindAllMarkers


#通路富集结果调整
#读取数据
MIA_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/MIA_gene.csv",sep = "\t")
MIA_data[1:4,]

LUAD_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/LUAD_gene.csv",sep = "\t")
LUAD_data[1:4,]


#拟时序通路富集更新
data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/tumor_CDS_DEG_gene_update_3_17_new.csv",sep = ",",header = T)
data[1:4,]

DEG=data[data$cluster%in%1,]$gene

DEG=data[data$cluster%in%"other",]$gene


library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(ggplot2)
library(org.Hs.eg.db)
library(BiocManager)


geneList1<-bitr(unique(MIA_data$V1), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)

geneList1<-bitr(unique(DEG), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)


head(geneList1)




GO2<-enrichGO( geneList1$ENTREZID,#GO富集分析
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               #pvalueCutoff = 0.03,#设定p值阈值
               qvalueCutoff = 0.05,#设定q值阈值
               pAdjustMethod = "BH",
               readable = TRUE,minGSSize=24,maxGSSize=290)

GO2@result$Description


kk2<-enrichKEGG( geneList1$ENTREZID,#GO富集分析
                 organism = "hsa",
                 keyType = "kegg",#设定读取的gene ID类型
                 #ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                 pvalueCutoff = 0.05,#设定p值阈值
                 #qvalueCutoff = 0.05,#设定q值阈值
                 pAdjustMethod = "BH")



df<-GO2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/monocle_up_Go.csv",sep = ",",col.names = T,row.names = T)


#肿瘤金标准驱动基因
data_golden<-read.table(file = "/Users/wangjun/Desktop/IntOGen-DriverGenes_LUAD.tsv",sep = "\t",header = T)

interact(data_golden$Symbol,DEG)      
            
intersect(data_golden$Symbol,DEG)




df<-kk2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/MIA_Kegg.csv",sep = ",",col.names = T,row.names = T)








geneList1<-bitr(unique(LUAD_data$V1), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)

head(geneList1)




GO2<-enrichGO( geneList1$ENTREZID,#GO富集分析
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 0.05,#设定p值阈值
               #qvalueCutoff = 0.05,#设定q值阈值
               pAdjustMethod = "BH",
               readable = TRUE)

kk2<-enrichKEGG( geneList1$ENTREZID,#GO富集分析
                 organism = "hsa",
                 keyType = "kegg",#设定读取的gene ID类型
                 #ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                 pvalueCutoff = 0.05,#设定p值阈值
                 #qvalueCutoff = 0.05,#设定q值阈值
                 pAdjustMethod = "BH")



df<-GO2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/LUAD_Go.csv",sep = ",",col.names = T,row.names = T)



df<-kk2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/LUAD_Kegg.csv",sep = ",",col.names = T,row.names = T)







#拟时序分支1的差异基因富集

fate_1_gene<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/tumor_CDS_DEG_gene.csv",sep = ",",header = T)
fate_1_gene[1:4,]

gene<-marker_gene[marker_gene$cluster==1,]$gene



geneList1<-bitr(unique(gene), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)

head(geneList1)




GO2<-enrichGO( geneList1$ENTREZID,#GO富集分析
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 0.05,#设定p值阈值
               #qvalueCutoff = 0.05,#设定q值阈值
               pAdjustMethod = "BH",
               readable = TRUE)

kk2<-enrichKEGG( geneList1$ENTREZID,#GO富集分析
                 organism = "hsa",
                 keyType = "kegg",#设定读取的gene ID类型
                 #ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                 pvalueCutoff = 0.05,#设定p值阈值
                 #qvalueCutoff = 0.05,#设定q值阈值
                 pAdjustMethod = "BH")







df<-GO2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/fate1_Go.csv",sep = ",",col.names = T,row.names = T)



df<-kk2@result
df[1:4,]
write.table(df,file = "/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/fate1_Kegg.csv",sep = ",",col.names = T,row.names = T)






#计算TCR的多样性指数
library(vegan)


#读取FD_LC1的TCR数据
FD1_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD1_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD1_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD1_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD1_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD1_LC1_data$Frequency, index = "inv")
Inverse.Simpson


FD1_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD1_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD1_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD1_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD1_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD1_LC2_data$Frequency, index = "inv")
Inverse.Simpson


FD1_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD1_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD1_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD1_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD1_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD1_LC2_data$Frequency, index = "inv")
Inverse.Simpson


FD1_LC3_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD1_LC3/Clonotype_Frequency.csv",header = T,sep = ",")
FD1_LC3_data[1:4,]

Shannon.Wiener <- diversity(FD1_LC3_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD1_LC3_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD1_LC3_data$Frequency, index = "inv")
Inverse.Simpson


FD2_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD2_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD2_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD2_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD2_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD2_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD2_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD2_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD2_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD2_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD2_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD2_LC2_data$Frequency, index = "inv")
Inverse.Simpson







FD4_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD4_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD4_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD4_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD4_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD4_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD4_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD4_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD4_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD4_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD4_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD4_LC2_data$Frequency, index = "inv")
Inverse.Simpson







FD5_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD5_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD5_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD5_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD5_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD5_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD5_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD5_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD5_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD5_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD5_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD5_LC2_data$Frequency, index = "inv")
Inverse.Simpson







FD8_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD8_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD8_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD8_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD8_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD8_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD8_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD8_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD8_LC2_data[1:4,]

sum(FD8_LC2_data$Frequency)


Shannon.Wiener <- diversity(FD8_LC2_data$Frequency/2215, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD8_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD8_LC2_data$Frequency, index = "inv")
Inverse.Simpson


d50.index <-function(x,type='raw'){
  if(type=='raw'){
    myfreqs <- table(x)/length(x) 
    myvec <- sort(as.numeric(myfreqs),decreasing = T)
  }else{
    myvec=sort(x,decreasing = T)
  }
  len=length(myvec)
  state=cumsum(myvec)>sum(myvec)/2  
  (len -sum(state))/len # 
}

d50.index(FD8_LC2_data$Frequency)

100*0.01*log2(0.01)

?diversity










FD9_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD9_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD9_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD9_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD9_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD9_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD9_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD9_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD9_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD9_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD9_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD9_LC2_data$Frequency, index = "inv")
Inverse.Simpson







FD14_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD14_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD14_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD14_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD14_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD14_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD14_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD14_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD14_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD14_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD14_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD14_LC2_data$Frequency, index = "inv")
Inverse.Simpson









FD16_LC1_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD16_LC1/Clonotype_Frequency.csv",header = T,sep = ",")
FD16_LC1_data[1:4,]

Shannon.Wiener <- diversity(FD16_LC1_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD16_LC1_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD16_LC1_data$Frequency, index = "inv")
Inverse.Simpson




FD16_LC2_data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR/FD16_LC2/Clonotype_Frequency.csv",header = T,sep = ",")
FD16_LC2_data[1:4,]

Shannon.Wiener <- diversity(FD16_LC2_data$Frequency, index = "shannon")
Shannon.Wiener
Simpson <- diversity(FD16_LC2_data$Frequency, index = "simpson")
Simpson
Inverse.Simpson <- diversity(FD16_LC2_data$Frequency, index = "inv")
Inverse.Simpson



TCR_LUAD<-c(7.127372,6.890004,7.587322,6.719225,6.294129,7.403371)
TCR_MIA<-c(5.214451,6.173589,7.173777,7.066926,6.343318,6.730217)

wilcox.test(TCR_MIA,TCR_LUAD,paired=TRUE)


TCR<-data.frame(matrix(NA,14,2))
colnames(TCR)<-c("Samples","Index")
TCR[1:7,]$Samples="MIA"
TCR[8:14,]$Samples="LUAD"


TCR[1:7,]$Index=TCR_MIA
TCR[8:14,]$Index=TCR_LUAD

TCR$Samples<-factor(TCR$Samples,levels = c("MIA","LUAD"))






TCR_LUAD<-c(7.127372,6.890004,7.587322,6.719225,6.294129,7.403371,7.236271,6.342968,7.62522,7.255371)
TCR_MIA<-c(5.214451,6.173589,7.173777,7.066926,6.343318,6.730217)

wilcox.test(TCR_MIA,TCR_LUAD)


TCR<-data.frame(matrix(NA,17,2))
colnames(TCR)<-c("Samples","Index")
TCR[1:7,]$Samples="MIA"
TCR[8:17,]$Samples="LUAD"


TCR[1:7,]$Index=TCR_MIA
TCR[8:17,]$Index=TCR_LUAD

TCR$Samples<-factor(TCR$Samples,levels = c("MIA","LUAD"))





TCR_new<-TCR
p<-ggplot(data=TCR_new)+
  geom_point(aes(x=factor(Samples),y=Index,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Shannon_Wiener Index")+
  ylim(5,8)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/TCR/Shannon_Wiener/total.png",dpi = 1000,width = 1.87,height = 3.41)







#进行T细胞的绘图

#载入数据

setwd("/Users/wangjun/Desktop/最终方案/figure/")

T_total<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/T_total_update.rds")


T_total_test<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/T_total_12_10.rds")
T_total<-T_total_test

table(Idents(T_total_test))

table(subset(T_total_test,idents = c("CD4_C2"))$orig.ident)

Idents(T_total)<-"orig.ident"
table(Idents(T_total))

T_total<-subset(T_total,idents="T3",invert=T)
Idents(T_total)<-"Cluster"



T_total=FindVariableFeatures(object = T_total)
f=VariableFeatures(object = T_total)
T_total <- ScaleData(T_total, features = f)
T_total <- RunPCA(T_total, features = f)
T_total <- FindNeighbors(T_total, dims = 1:30)
T_total <- FindClusters(T_total, resolution = 0.2)



library(paletteer) 
library(scales)
pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal[c(2,3,7,10,5,8,4,1,9)]

pal_T<-pal[c(2,3,5,1,7,4,6)]
pal_T[5]<-"#E97123"

p5<-DimPlot(T_total,reduction = 'tsne', label.box = T,  label = F,repel = T,cols = pal_T)
ggsave(p5,file="/Users/wangjun/Desktop/T_dimplot.png",dpi = 1000,width = 7.45,height = 4.67)










show_col(colorRampPalette(c("#E64B35FF","#EFC000FF"))(4))
show_col("#EFC000FF")

#pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,5,8,7,1,9)] 
pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,10,5,8,4,1,9)] 

pal_new
pal_new<-(colorRampPalette(pal)(16))

pal_new[14]<-"#B5553D"
pal_new[15]<-"#DA4B3A"

show_col(pal_new)

#UMAP
T_total <- RunUMAP(T_total, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T_total, reduction = "umap", group.by = "Cluster",cols = pal_new)
p3 <- DimPlot(T_total, reduction = "umap", label = TRUE, repel = TRUE)#, cols =colorlist)
p1


Idents(T_total)<-"Cluster"

#鉴定marker
marker<-c(
  "IL7R","CD4",#CD4_T
  "CD8A","CD8B",#CD8_T
  "TCF7","SELL","LEF1","CCR7",#naive_t
  "LAG3","TIGIT","PDCD1","CTLA4","KLRG1",#exhausted
  "IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7",#cytotoxic
  "IL2RA","FOXP3","IKZF2","TGFB1","TGFB3","TGFBI","TGFBR1",#treg
  "MAF","CXCR5","CXCL13", # and PDCD1 #tfh
  "IRF4","CREM","NR4A2",#th17
  "STAT4","IL12RB2", # and IFNG  #th1
  "GATA3","STAT6","IL4",#th2
  "TRDC","TRGC2","TRGC1","CD38","MKI67",#gd_t
  #"ITGAE",#Trm
  "CXCR4"#Tem
  #"CXCR6","PTGER4" #CD4 Trm
  "IL15RA","CD44","MBD2","BLIMP1"#memory
  )

Idents(T_total)<-"Cluster"

DotPlot(T_total, features = marker, cols = c("#4DBBD5FF","#E64B35FF")) + RotatedAxis()



table(Idents(T_total))

#亚群鉴定

p1<-DotPlot(T_total, features = marker, cols = c("#4DBBD5FF","#E64B35FF")) + RotatedAxis()


ggsave(p1,file="/Users/wangjun/Desktop/T_dotplot.png",dpi=1000,width = 13.4,height = 4.31)


# CD8_C1 CD8+ GZMK+ T cells 炎性衰老的保守标志,CD8+ Tem
# CD8_C2 CD8+ GZMK+ T cells
# CD8_C3 Exhausted T cells
# CD8_C4 ? CD8+ GZMB+ T cells
# CD8_C5 CD8+ GZMB+ T cells
# CD8_C6 CD8+ GZMB+ T cells
# CD8_C7 CD8+ GZMB+ T cells

# CD4_C1 Naïve T cells
# CD4_C2 Treg cells
# CD4_C3 T helper
# CD4_C4 Naïve T cells
# CD4_C5 memory T cell 
# CD4_C6 memory T cell 
# CD4_C7 memory T cell 


table(Idents(T_total))

#鉴定出相应的细胞类群后
new.cluster.ids <- c("CD8+ GZMK+ T cells","CD8+ GZMK+ T cells","Exhausted T cells","CD8+ GZMB+ T cells","CD8+ GZMB+ T cells","CD8+ GZMB+ T cells","CD8+ GZMB+ T cells","Naïve T cells","Th cells","Naïve T cells","Memory T cells","Memory T cells","Memory T cells","γδT_C1","γδT_C2","Treg cells")
names(new.cluster.ids) <- levels(T_total)
T_total <- RenameIdents(T_total, new.cluster.ids)

table(Idents(T_total))
Idents(T_total)<-factor(Idents(T_total),levels=c("Naïve T cells","Treg cells","Memory T cells","Th cells","CD8+ GZMK+ T cells","CD8+ GZMB+ T cells","Exhausted T cells","γδT_C1","γδT_C2"))


p1 <- DimPlot(T_total, reduction = "umap")
p3 <- DimPlot(T_total, reduction = "umap", label = F, repel = TRUE,cols = pal_T)#, cols =colorlist)
p3


#TNSE
T_total <- RunTSNE(T_total, reduction = "pca", dims = 1:20)
p5 <- DimPlot(T_total, reduction = "tsne",group.by = "Cluster")#,cols = pal_new)
  

table(Idents(T_total))
length(Idents(T_total))


#计算比例

df = data.frame(clu=names(table(Idents(T_total))),
                per=sprintf("%1.2f%%", 100*table(Idents(T_total))/length(Idents(T_total))))
df
T_total$per = df[match(Idents(T_total),df$clu),2]

T_total$new = paste0(Idents(T_total),"(",T_total$per,")")

table(T_total$new )

Idents(T_total)<-"new"

Idents(T_total)<- factor(Idents(T_total),levels = c("Naïve T cells(66.57%)","Treg cells(10.13%)","CD8+ GZMK+ T cells(10.63%)","CD8+ GZMB+ T cells(10.14%)","Exhausted T cells(0.70%)","γδT_C1(1.21%)","γδT_C2(0.63%)"))

Idents(T_total)<- factor(Idents(T_total),levels = c("Naïve T cells(52.17%)","Treg cells(10.13%)","Th cells(7.62%)","Memory T cells(6.78%)","CD8+ GZMK+ T cells(10.63%)","CD8+ GZMB+ T cells(10.14%)","Exhausted T cells(0.70%)","γδT_C1(1.21%)","γδT_C2(0.63%)"))

DimPlot(T_total,reduction = 'tsne', label.box = T,  label = F,repel = T,cols = pal_T)

pal_T<-pal[c(3,4,9,7,5,1,8,6,2)]
p6<-DimPlot(T_total,reduction = 'tsne', label.box = T,  label = F,repel = T,cols = pal_T)
p6
ggsave(p6,file="/Users/wangjun/Desktop/T_dimplot.png",dpi=1000,width = 7.62,height = 4.7)


table(Idents(T_total))

T_total<-AddMetaData(T_total,Idents(T_total),col.name="Annotation")

?AddMetaData

table(T_total$Cluster)


saveRDS(T_total,file="/Users/wangjun/Desktop/T_total_3_31.rds")


pal

#TNSE
T_total <- RunTSNE(T_total, reduction = "pca", dims = 1:20)
p5 <- DimPlot(T_total, reduction = "tsne",group.by = "Cluster")#,cols = pal_new)

p5

ggsave(p5,file="/Users/wangjun/Desktop/figure2_1.png",dpi = 1000,width = 5.42,height = 4.62)



library(Seurat)

#计算比例在样本中的分布
saveRDS(T_total,file="/Users/wangjun/Desktop/T_total_3_26.rds")
T_total<-readRDS(file="/Users/wangjun/Desktop/T_total_3_26.rds")

table(Idents(T_total))

table(subset(T_total,idents = "Naïve T cells(66.57%)")$orig.ident)
table(subset(T_total,idents = "Treg cells(10.13%)")$orig.ident)
table(subset(T_total,idents = "CD8+ GZMK+ T cells(10.63%)")$orig.ident)
table(subset(T_total,idents = "CD8+ GZMB+ T cells(10.14%)")$orig.ident)
table(subset(T_total,idents = "Exhausted T cells(0.70%)")$orig.ident)
table(subset(T_total,idents = "γδT_C1(1.21%)")$orig.ident)
table(subset(T_total,idents = "γδT_C2(0.63%)")$orig.ident)


table(subset(T_total,idents = "Naïve T cells(52.17%)")$orig.ident)
table(subset(T_total,idents = "Treg cells(10.13%)")$orig.ident)
table(subset(T_total,idents = "Th cells(7.62%)")$orig.ident)
table(subset(T_total,idents = "Memory T cells(6.78%)")$orig.ident)
table(subset(T_total,idents = "CD8+ GZMK+ T cells(10.63%)")$orig.ident)

table(subset(T_total,idents = "CD8+ GZMB+ T cells(10.14%)")$orig.ident)
table(subset(T_total,idents = "Exhausted T cells(0.70%)")$orig.ident)


table(subset(T_total,idents = "γδT_C1(1.21%)")$orig.ident)
table(subset(T_total,idents = "γδT_C2(0.63%)")$orig.ident)







#计算注释后细胞比例
T_cells<-read.table(file="/Users/wangjun/Desktop/T细胞比例.csv",header = T,sep = ",")

T_cells


for (i in 1:length(T_cells$names)){
#for (i in 3:4){
  MIA<-as.numeric(T_cells[i,c(3,4,7,9,11,15)])
  LUAD<-as.numeric(T_cells[i,c(2,5,6,8,10,12,13,14,16,17)])
  p=wilcox.test(MIA,LUAD)
  print(paste0(T_cells$names[i]," : ",p$p.value))
}


for (i in 1:length(T_cells$names)){
  #for (i in 3:4){
  MIA<-as.numeric(T_cells[i,c(3,4,7,9,11,15)])
  LUAD<-as.numeric(T_cells[i,c(2,5,6,8,10,14)])
  p=wilcox.test(MIA,LUAD,alternative = "greater",paired = T)
  print(paste0(T_cells$names[i]," : ",p$p.value))
}

GZMB_MIA<-T_cells[6,c(3,4,7,9,11,15)]
GZMB_LUAD<-T_cells[6,c(2,5,6,8,10,12,13,14,16,17)]


GZMB<-data.frame(matrix(NA,16,2))
colnames(GZMB)<-c("Samples","Proportion")
GZMB[1:6,]$Samples="MIA"
GZMB[7:16,]$Samples="LUAD"


GZMB[1:6,]$Proportion=as.numeric(GZMB_MIA)
GZMB[7:16,]$Proportion=as.numeric(GZMB_LUAD)

GZMB$Samples<-factor(GZMB$Samples,levels = c("MIA","LUAD"))




GZMB<-data.frame(matrix(NA,12,2))
colnames(GZMB)<-c("Samples","Proportion")
GZMB[1:6,]$Samples="MIA"
GZMB[7:12,]$Samples="LUAD"


GZMB[1:6,]$Proportion=as.numeric(GZMB_MIA)
GZMB[7:12,]$Proportion=as.numeric(GZMB_LUAD[c(1:5,8)])

GZMB$Samples<-factor(GZMB$Samples,levels = c("MIA","LUAD"))





GZMB_new<-GZMB[c(6,12),]
#GZMB_new<-GZMB
p<-ggplot(data=GZMB_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/GZMB/6.png",dpi = 1000,width = 2.12,height = 3.44)




















new_marker<-c(
  "TCF7","SELL","LEF1","CCR7",#Naïve markers
  "CD44", "IL7R", "MBD2",#Memory T
  "IL2RA","FOXP3","IKZF2",#Treg  markers
  "TRDC","TRGC2","TRGC1","CD38","MKI67",#gd_T
  "IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7",#Cytokines and effector molecules
  "CD28","TNFRSF14","ICOS","TNFRSF9",#Co-stimulatory molecules
  "EOMES","HOPX","TBX21",'ZEB2','ZNF683','HIF1A','ID2','TOX',#Transcription factors
  "LAG3","TIGIT","PDCD1","HAVCR2","CTLA4"#Exhausted
)

Idents(T_total)<-"Cluster"


cd8_1<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C1"))[new_marker,])
cd8_2<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C2"))[new_marker,])
cd8_3<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C3"))[new_marker,])
cd8_4<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C4"))[new_marker,])
cd8_5<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C5"))[new_marker,])
cd8_6<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C6"))[new_marker,])
cd8_7<-rowMeans(GetAssayData(subset(T_total,idents="CD8_C7"))[new_marker,])


cd4_1<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C1"))[new_marker,])
cd4_2<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C2"))[new_marker,])
cd4_3<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C3"))[new_marker,])
cd4_4<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C4"))[new_marker,])
cd4_5<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C5"))[new_marker,])
cd4_6<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C6"))[new_marker,])
cd4_7<-rowMeans(GetAssayData(subset(T_total,idents="CD4_C7"))[new_marker,])


gd_1<-rowMeans(GetAssayData(subset(T_total,idents="γδT_C1"))[new_marker,])
gd_2<-rowMeans(GetAssayData(subset(T_total,idents="γδT_C2"))[new_marker,])






length(new_marker)
#进行矩阵拼接
data_matrix<-data.frame(matrix(NA,40,16))



colnames(data_matrix)<-c("CD8-C1","CD8-C2","CD8-C3","CD8-C4","CD8-C5","CD8-C6","CD8-C7","CD4-C1","CD4-C2","CD4-C3","CD4-C4","CD4-C5","CD4-C6","CD4-C7","γδ_C1","γδ_C2")
rownames(data_matrix)<-new_marker

data_matrix[,1]<-as.numeric(cd8_1)
data_matrix[,2]<-as.numeric(cd8_2)
data_matrix[,3]<-as.numeric(cd8_3)
data_matrix[,4]<-as.numeric(cd8_4)
data_matrix[,5]<-as.numeric(cd8_5)
data_matrix[,6]<-as.numeric(cd8_6)
data_matrix[,7]<-as.numeric(cd8_7)

data_matrix[,8]<-as.numeric(cd4_1)
data_matrix[,9]<-as.numeric(cd4_2)
data_matrix[,10]<-as.numeric(cd4_3)
data_matrix[,11]<-as.numeric(cd4_4)
data_matrix[,12]<-as.numeric(cd4_5)
data_matrix[,13]<-as.numeric(cd4_6)
data_matrix[,14]<-as.numeric(cd4_7)


data_matrix[,15]<-as.numeric(gd_1)
data_matrix[,16]<-as.numeric(gd_2)


data_zscore<-data_matrix
# zscore
for (i in 1:length(rownames(data_matrix))) {
  data_zscore[i,] = (data_matrix[i,] - mean(as.numeric(data_matrix[i,]))) / sd(data_matrix[i,])
}


data_zscore_new<-data_zscore[,c("CD8-C1","CD8-C2","CD8-C4","CD8-C5","CD8-C6","CD8-C7","CD8-C3","CD4-C1","CD4-C4","CD4-C3","CD4-C5","CD4-C6","CD4-C7","CD4-C2","γδ_C1","γδ_C2")]

library(pheatmap)
p2<-pheatmap(data_zscore_new,cluster_cols=F,cluster_rows = F,color = colorRampPalette(colors = c("#4DBBD5FF","white","#E64B35FF"))(100))
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/figure2_2.png",dpi = 1000,width = 3.84,height = 5.87)

p2<-pheatmap(data_zscore_new,cluster_cols=F,cluster_rows = F)
p2

ggsave(p2,file="/Users/wangjun/Desktop/figure2_2(2).png",dpi = 1000,width = 3.84,height = 5.87)

write.table(data_matrix,file="/Users/wangjun/Desktop/最终方案/figure/T细胞marker表达量.csv",col.names = T,row.names = T,sep = ",")


 
#CD8_C3分析

CD8_C3_MIA<-c(0.004591368,0,0.001901593,0.000248324,0,0.002728513)
CD8_C3_LUAD<-c(0.002240896,0.011021619,0.025987006,0.000688705,0,0.004106776)

wilcox.test(CD8_C3_LUAD,CD8_C3_MIA,paired = TRUE)



CD8_C3_MIA<-c(0.004591368,0,0.001901593,0.000248324,0,0.002728513)
CD8_C3_LUAD<-c(0.002240896,0.011021619,0.025987006,0.000688705,0,0.004106776,0.001286174,0.002327747,0.001393405,0.036308623)

p=wilcox.test(CD8_C3_LUAD,CD8_C3_MIA)




CD8_C3<-data.frame(matrix(NA,16,2))
colnames(CD8_C3)<-c("Samples","Proportion")
CD8_C3[1:6,]$Samples="MIA"
CD8_C3[7:16,]$Samples="LUAD"


CD8_C3[1:6,]$Proportion=CD8_C3_MIA
CD8_C3[7:16,]$Proportion=CD8_C3_LUAD

CD8_C3$Samples<-factor(CD8_C3$Samples,levels = c("MIA","LUAD"))




CD8_C3


CD8_C3_new<-CD8_C3[c(1,8),]
CD8_C3_new<-CD8_C3
p<-ggplot(data=CD8_C3_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.04)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/CD8_C3/total.png",dpi = 1000,width = 2.23,height = 3.44)



#gdT_C1

gdT_C1_MIA<-c(0.022956841,0.008361204,0.004278583,0.0104296,0.013888889,0.006821282)
gdT_C1_LUAD<-c(0.002240896,0.005086901,0.011494253,0.012396694,0.024282561,0.0275154)

wilcox.test(gdT_C1_LUAD,gdT_C1_MIA,paired = T)
shapiro.test(gdT_C1_MIA)


gdT_C1_MIA<-c(0.022956841,0.008361204,0.004278583,0.0104296,0.013888889,0.006821282)
gdT_C1_LUAD<-c(0.002240896,0.005086901,0.011494253,0.012396694,0.024282561,0.0275154,0.008360129,0.011173184,0.015791918,0.024710035)

wilcox.test(gdT_C1_LUAD,gdT_C1_MIA)



gdT_C1<-data.frame(matrix(NA,16,2))
colnames(gdT_C1)<-c("Samples","Proportion")
gdT_C1[1:6,]$Samples="MIA"
gdT_C1[7:16,]$Samples="LUAD"


gdT_C1[1:6,]$Proportion=gdT_C1_MIA
gdT_C1[7:16,]$Proportion=gdT_C1_LUAD

gdT_C1$Samples<-factor(gdT_C1$Samples,levels = c("MIA","LUAD"))




gdT_C1_new<-gdT_C1[c(7,14),]
gdT_C1_new<-gdT_C1
p<-ggplot(data=gdT_C1_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.03)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/gdT_C1/total.png",dpi = 1000,width = 2.12,height = 3.44)



#gd_C2
gdT_C2_MIA<-c(0.014692378,0.010033445,0.007844069,0.017879315,0.003472222,0)
gdT_C2_LUAD<-c(0.004481793,0.002119542,0.004497751,0.003443526,0.004415011,0.001232033)
wilcox.test(gdT_C2_MIA,gdT_C2_LUAD,paired = TRUE)






gdT_C2_MIA<-c(0.014692378,0.010033445,0.007844069,0.017879315,0.003472222,0)
gdT_C2_LUAD<-c(0.004481793,0.002119542,0.004497751,0.003443526,0.004415011,0.001232033,0,0.000465549,0.013469577,0.004034291)
wilcox.test(gdT_C2_MIA,gdT_C2_LUAD)



gdT_C2<-data.frame(matrix(NA,16,2))
colnames(gdT_C2)<-c("Samples","Proportion")
gdT_C2[1:6,]$Samples="MIA"
gdT_C2[7:16,]$Samples="LUAD"


gdT_C2[1:6,]$Proportion=gdT_C2_MIA
gdT_C2[7:16,]$Proportion=gdT_C2_LUAD

gdT_C2$Samples<-factor(gdT_C2$Samples,levels = c("MIA","LUAD"))

dim(gdT_C2)
gdT_C2_new<-gdT_C2[c(7,14),]
gdT_C2_new<-gdT_C2
p<-ggplot(data=gdT_C2_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.02)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/gdT_C2/FD14_LUAD_MIA.png",dpi = 1000,width = 2.12,height = 3.44)



































#CD4_C2
CD4_C2_MIA<-c(0.214334009,0.19464145,0.053422371,0.079669031,0.072682324,0.036395147,0.066846725)
CD4_C2_LUAD<-c(0.171933583,0.152392947,0.13188929,0.076923077,0.043956044,0.174190971,0.124614436,0.182407407,0.098383372,0.242022582)
wilcox.test(CD4_C2_MIA,CD4_C2_LUAD)


table(subset(T_total,idents="CD4_C2")$orig.ident)


CD4_C2<-data.frame(matrix(NA,17,2))
colnames(CD4_C2)<-c("Samples","Proportion")
CD4_C2[1:7,]$Samples="MIA"
CD4_C2[8:17,]$Samples="LUAD"


CD4_C2[1:7,]$Proportion=CD4_C2_MIA
CD4_C2[8:17,]$Proportion=CD4_C2_LUAD

CD4_C2$Samples<-factor(CD4_C2$Samples,levels = c("MIA","LUAD"))


#CD4_C2_new<-CD4_C2[c(7,14),]
CD4_C2_new<-CD4_C2
p<-ggplot(data=CD4_C2_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  #ylim(0,0.05)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/CD4_C2/FD14_LUAD_MIA.png",dpi = 1000,width = 2.23,height = 3.44)



CD4_C2_seurat<-subset(T_total,idents="CD4_C2")



CD4_C2_seurat <- FindNeighbors(CD4_C2_seurat, dims = 1:30)
#设置resolution
CD4_C2_seurat <- FindClusters(CD4_C2_seurat, resolution = 0.1)

DotPlot(CD4_C2_seurat, features = new_marker, cols = c("green","red")) + RotatedAxis()

table(subset(CD4_C2_seurat,idents=0)$orig.ident)

CD4_C2_new<-subset(CD4_C2_seurat,idents=0)



#CD4_C2
CD4_C2_MIA<-c(0.061524334,0.051839465,0.074637509,0.068537373,0.034722222,0.057298772)
CD4_C2_LUAD<-c(0.133893557,0.144128868,0.09870065,0.070247934,0.221381745,0.093358105,0.039735099,0.151129363,0.177839851,0.087459807)
wilcox.test(CD4_C2_MIA,CD4_C2_LUAD)




CD4_C2_MIA<-c(0.061524334,0.051839465,0.074637509,0.068537373,0.034722222,0.057298772)
CD4_C2_LUAD<-c(0.133893557,0.144128868,0.09870065,0.070247934,0.221381745,0.093358105)
wilcox.test(CD4_C2_MIA,CD4_C2_LUAD,paired = TRUE)




CD4_C2<-data.frame(matrix(NA,14,2))
colnames(CD4_C2)<-c("Samples","Proportion")
CD4_C2[1:7,]$Samples="MIA"
CD4_C2[8:14,]$Samples="LUAD"


CD4_C2[1:7,]$Proportion=CD4_C2_MIA
CD4_C2[8:14,]$Proportion=CD4_C2_LUAD

CD4_C2$Samples<-factor(CD4_C2$Samples,levels = c("MIA","LUAD"))


CD4_C2_new<-CD4_C2[c(7,14),]
#CD4_C2_new<-CD4_C2
p<-ggplot(data=CD4_C2_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.25)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/CD4_C2/FD14_LUAD_MIA.png",dpi = 1000,width = 2.23,height = 3.44)




table(CDS[CDS$State==1]$Type)

table(CDS$Type[CDS$State==1])
table(CDS$Type[CDS$State==2])
table(CDS$Type[CDS$State==3])



T_total_part<-subset(T_total,idents="CD4_C2",invert=T)

table(Idents(T_total_part))

table(Idents)


new.cluster.ids <- c("CD4_C2")
levels(CD4_C2_new)
names(new.cluster.ids) <- levels(CD4_C2_new)
new.cluster.ids
CD4_C2_new <- RenameIdents(CD4_C2_new, new.cluster.ids)

T_total_new<-merge(T_total_part,y=CD4_C2_new,add.cell.ids = NULL)

table(Idents(T_total_new))

saveRDS(T_total_new,file="./T_total_12_10.rds")

table(Idents(T_total))



#绘制NK细胞的聚类图

NK_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/NK_cells_update.rds")
table(Idents(NK_cells))

NK_cells<-AddMetaData(NK_cells,Idents(NK_cells),col.name = "Cluster")

Idents(NK_cells)<-"orig.ident"
table(Idents(NK_cells))
NK_cells<-subset(NK_cells,idents="T3",invert=T)

Idents(NK_cells)<-"Cluster"


df = data.frame(clu=names(table(NK_cells$Cluster)),
                per=sprintf("%1.2f%%", 100*table(NK_cells$Cluster)/length(NK_cells$Cluster)))
NK_cells$per = df[match(NK_cells$Cluster,df$clu),2]
NK_cells$new = paste0(NK_cells$Cluster,"(",NK_cells$per,")")
table(NK_cells$new)
Idents(NK_cells)<-"new"

Idents(NK_cells)<-factor(Idents(NK_cells),levels = c("CD56dimCD16+NK cells(90.49%)","CD56brightCD16-NK cells(9.51%)"))






#UMAP
NK_cells <- RunUMAP(NK_cells, reduction = "pca", dims = 1:20)
p1 <- DimPlot(NK_cells, reduction = "umap",cols = c("#E64B35FF","#4DBBD5FF"))
p3 <- DimPlot(NK_cells, reduction = "umap", label = TRUE, repel = TRUE)#, cols =colorlist)
p1


#TNSE
NK_cells <- RunTSNE(NK_cells, reduction = "pca", dims = 1:10)
p5 <- DimPlot(NK_cells, reduction = "tsne",cols = c("#E64B35FF","#4DBBD5FF"))
p5


p5 <- DimPlot(NK_cells, reduction = "tsne",cols = pal[c(1,5,3,2,7,4)])
p5

ggsave(p5,file="/Users/wangjun/Desktop/NK_final_new.png",dpi = 1000,width = 6.73,height =4.27 )


p5 <- DimPlot(NK_cells, reduction = "tsne",split.by = "orig.ident",cols = c("#E64B35FF","#4DBBD5FF"))+
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title.x=element_text(size=5,face = "bold"),axis.title.y=element_text(size=5,face = "bold"),axis.text.x=element_text(size=5,face = "bold",color = "black"),axis.text.y=element_text(size=5,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 5))

p5
ggsave(p5,file="/Users/wangjun/Desktop/NK_orig_ident_new.png",dpi = 1000,width = 17.2,height =1.94)






pal

NK_markers<-c(
  "GZMB","GZMH","IFNG", #效应器功能
  "KIR",#杀手细胞免疫球蛋白样受体
  "S1PR5","S1PR1PR1","CXCR1","CXCR2","CX3CR1",#组织循环
  "CD27","TNFSF4","TNFS8","TNFSF13","TNFSF13B",#免疫调节
  "CCR7","SELL","CXCR3","CCR5"#组织归巢
)


markers<-c(
  "NCAM1",
  "FCGR3A","FCGR3B",
  "NCR1","NCR2","NCR3"
)

NK_total_markers<-c(
  "NCAM1",
  "FCGR3A","FCGR3B",
  "GZMB","GZMH","IFNG", #效应器功能
  "KIR",#杀手细胞免疫球蛋白样受体
  "S1PR5","S1PR1PR1","CXCR1","CXCR2","CX3CR1",#组织循环
  "CD27","TNFSF4","TNFS8","TNFSF13","TNFSF13B",#免疫调节
  "CCR7","SELL","CXCR3","CCR5"#组织归巢
)

p5<-DotPlot(NK_cells,features = NK_total_markers, cols = c("#4DBBD5FF","#E64B35FF")) + RotatedAxis()
  theme(legend.key.size = unit(8, "pt"))
p5  

ggsave(p5,file="./NK_dotplot.png",dpi = 1000,width = 8.78,height =3.7 )




NK_cells <- FindNeighbors(NK_cells, reduction = "pca", dims = 1:30) 
NK_cells <- FindClusters(NK_cells, resolution = 0.2)


#计算CD56dimCD16+NK cells的比例
table(NK_cells$orig.ident)

table(subset(NK_cells,idents="CD56dimCD16+NK cells")$orig.ident)

table(subset(NK_cells,idents="CD56dimCD16+NK cells")$orig.ident)


table(Idents(NK_cells))




#CD4_C2
NK_plus_MIA<-c(0.908256881,0.966412214,0.920560748,0.965313653,0.93129771,0.926345609)
NK_plus_LUAD<-c(0.173913043,0.175,0.907120743,0.628571429,0.8375,0.972098214)
wilcox.test(NK_plus_MIA,NK_plus_LUAD,paired = TRUE)


NK_plus_MIA<-c(0.908256881,0.966412214,0.920560748,0.965313653,0.93129771,0.926345609)
NK_plus_LUAD<-c(0.173913043,0.175,0.907120743,0.628571429,0.8375,0.972098214,0.948542024,0.450819672,0.804469274,0.885981308)
wilcox.test(NK_plus_MIA,NK_plus_LUAD)



NK_plus<-data.frame(matrix(NA,14,2))
colnames(NK_plus)<-c("Samples","Proportion")
NK_plus[1:7,]$Samples="MIA"
NK_plus[8:14,]$Samples="LUAD"


NK_plus[1:7,]$Proportion=NK_plus_MIA
NK_plus[8:14,]$Proportion=NK_plus_LUAD

NK_plus$Samples<-factor(NK_plus$Samples,levels = c("MIA","LUAD"))

NK_plus_new<-NK_plus[c(7,14),]
#NK_plus_new<-NK_plus
p<-ggplot(data=NK_plus_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,1)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/NK_plus/FD14_LUAD_MIA.png",dpi = 1000,width = 2.12,height = 3.44)


test








#髓系细胞

Myeloid_cells<-readRDS(file="./macrophages_update.rds")
table(Idents(Myeloid_cells))

levels(Myeloid_cells)

Idents(Myeloid_cells)<-factor(Idents(Myeloid_cells),levels =c("Alveolar resident Macrophage","Perivascular resident Macrophage","Anti-inflammatory macrophages","Proliferating Macrophage","Monocytes Derived DC","Migratory Conventional Dendritic cell","Conventional Dendritic cell type 1","Conventional Dendritic cell type 2","Classical CD14+ Monocyte","Non-classical CD16+ Monocyte") )



#TNSE
Myeloid_cells <- RunTSNE(Myeloid_cells, reduction = "pca", dims = 1:20)
p5 <- DimPlot(Myeloid_cells, reduction = "tsne",cols = pal[1:10])
p5


ggsave(p5,file="./Myeloid_final.png",dpi = 1000,width = 6.97,height =3.99 )

p5 <- DimPlot(Myeloid_cells, reduction = "tsne",split.by = "orig.ident",cols = pal[1:10])+
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title.x=element_text(size=5,face = "bold"),axis.title.y=element_text(size=5,face = "bold"),axis.text.x=element_text(size=5,face = "bold",color = "black"),axis.text.y=element_text(size=5,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 5))

p5
ggsave(p5,file="./Myeloid_orig_ident.png",dpi = 1000,width = 17.2,height =1.94)





#肿瘤细胞通路富集结果
MIA_kegg<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/MIA_Kegg.csv",sep = ",",header  = T,row.names = 1)

MIA_kegg[1:4,1:4]

data<-MIA_kegg[MIA_kegg$select==1,c("Description","GeneRatio","BgRatio","p.adjust","Count")]


MIA_GO<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/MIA_Go.csv",sep = ",",header  = T,row.names = 1)

MIA_GO[1:4,1:4]
colnames(MIA_GO)
data<-MIA_GO[MIA_GO$select==1,c("Description","GeneRatio","BgRatio","p.adjust","Count")]


LUAD_Kegg<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/LUAD_Kegg.csv",sep = ",",header  = T,row.names = 1)

LUAD_Kegg[1:4,1:4]
colnames(LUAD_Kegg)
data<-LUAD_Kegg[LUAD_Kegg$select==1,c("Description","GeneRatio","BgRatio","p.adjust","Count")]






colnames(data)
data

#MIA_Kegg
(6/99)/(41/8209)
(6/99)/(88/8209)
(9/99)/(205/8209)
(7/99)/(136/8209)

data$`Enrich factor` <- c(12.13452,5.653581,3.640355,4.2679)


#MIA_Go
(5/170)/(41/18723)
(13/170)/(369/18723)
(13/170)/(371/18723)
(13/170)/(446/18723)
(10/170)/(284/18723)
(3/170)/(17/18723)

data$`Enrich factor` <- c(13.43113,3.880105,3.859188,3.210222,3.878003,19.43564)



#LUAD_Kegg
(14/177)/(78/8209)
(15/177)/(223/8209)
(12/177)/(157/8209)
(8/177)/(108/8209)
(7/177)/(92/8209)

data$`Enrich factor` <- c(8.324352,3.119632,3.544856,3.435447,3.528801)





data

data$`-log10(p.adjust)`=-log10(data$p.adjust)


data_new<-data[,c("Description","Count","Enrich factor","-log10(p.adjust)")]
colnames(data_new)<-c("Pathway","Count","Enrich_Factor","-log10(P_value)")

data_new
data_new<-data_new[order(-data_new$Enrich_Factor),]




p1<-ggplot() +  geom_bar(data = data_new, 
                         aes(y = Enrich_Factor, x = Pathway,fill=`-log10(P_value)`),
                         stat = "identity",
                         width = 0.7, 
                         position = position_dodge(width = 0.9)) + 
  scale_fill_gradient2(high="#E64B35FF",low ="white",midpoint = 1 )+ theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()+
  ylim(0,20)+
  scale_x_discrete(limits=factor(data_new$Pathway))+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))

p1



ggsave(p1,file="./MIA_KEGG.png",dpi = 1000,width = 6.72,height = 2.15)

ggsave(p1,file="./MIA_GO.png",dpi = 1000,width = 6.72,height = 2.87)

ggsave(p1,file="./LUAD_KEGG.png",dpi = 1000,width = 6.72,height = 2.42)














LUAD_GO<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG/LUAD_Go2.csv",sep = ",",header  = T,row.names = 1,quote ="")

dim(LUAD_GO)

LUAD_GO[1:4,1:4]
colnames(LUAD_GO)
data<-LUAD_GO[LUAD_GO$select==1,c("Description","GeneRatio","BgRatio","p.adjust","Count")]

colnames(data)
data

#LUAD_Go1
(18/274)/(277/18723)
(9/274)/(95/18723)
(8/274)/(87/18723)
(10/274)/(141/18723)
(11/274)/(189/18723)

data$`Enrich factor` <- c(4.440354,6.473569,6.283413,4.846249,3.977002)




data$`-log10(p.adjust)`=-log10(data$p.adjust)


data_new<-data[,c("Description","Count","Enrich factor","-log10(p.adjust)")]
colnames(data_new)<-c("Pathway","Count","Enrich_Factor","-log10(P_value)")

data_new
data_new<-data_new[order(-data_new$Enrich_Factor),]


p1<-ggplot() +  geom_bar(data = data_new, 
                         aes(y = Enrich_Factor, x = Pathway,fill=`-log10(P_value)`),
                         stat = "identity",
                         width = 0.7, 
                         position = position_dodge(width = 0.9)) + 
  scale_fill_gradient2(high="#E64B35FF",low ="white",midpoint = 1 )+ theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()+
  ylim(0,20)+
  scale_x_discrete(limits=factor(data_new$Pathway))+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))

p1

ggsave(p1,file="./LUAD_GO_1.png",dpi = 1000,width = 6.72,height = 2.42)









data<-LUAD_GO[LUAD_GO$select==2,c("Description","GeneRatio","BgRatio","p.adjust","Count")]

colnames(data)
data

#LUAD_Go1
(12/274)/(62/18723)
(9/274)/(34/18723)
(15/274)/(216/18723)


data$`Enrich factor` <- c(13.22557,18.08791,4.745286)




data$`-log10(p.adjust)`=-log10(data$p.adjust)


data_new<-data[,c("Description","Count","Enrich factor","-log10(p.adjust)")]
colnames(data_new)<-c("Pathway","Count","Enrich_Factor","-log10(P_value)")

data_new
data_new<-data_new[order(-data_new$Enrich_Factor),]


p1<-ggplot() +  geom_bar(data = data_new, 
                         aes(y = Enrich_Factor, x = Pathway,fill=`-log10(P_value)`),
                         stat = "identity",
                         width = 0.7, 
                         position = position_dodge(width = 0.9)) + 
  scale_fill_gradient2(high="#E64B35FF",low ="white",midpoint = 1 )+ theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()+
  ylim(0,20)+
  scale_x_discrete(limits=factor(data_new$Pathway))+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))

p1

ggsave(p1,file="./LUAD_GO_2.png",dpi = 1000,width = 6.72,height = 1.65)



pbmc<-readRDS(file="/Users/wangjun/Desktop/pbmc_final.rds")

table(pbmc$orig.ident)










data<-LUAD_GO[LUAD_GO$select==3,c("Description","GeneRatio","BgRatio","p.adjust","Count")]

colnames(data)
data

#LUAD_Go1
(21/274)/(437/18723)
(18/274)/(467/18723)
(14/274)/(309/18723)
(14/274)/(310/18723)
(8/274)/(123/18723)
(16/274)/(420/18723)
(10/274)/(221/18723)
(6/274)/(96/18723)



data$`Enrich factor` <- c(3.283694,2.633786,3.095954,3.085967,4.444365,2.603128,3.091951,4.270757)




data$`-log10(p.adjust)`=-log10(data$p.adjust)


data_new<-data[,c("Description","Count","Enrich factor","-log10(p.adjust)")]
colnames(data_new)<-c("Pathway","Count","Enrich_Factor","-log10(P_value)")

data_new
data_new<-data_new[order(-data_new$Enrich_Factor),]


p1<-ggplot() +  geom_bar(data = data_new, 
                         aes(y = Enrich_Factor, x = Pathway,fill=`-log10(P_value)`),
                         stat = "identity",
                         width = 0.7, 
                         position = position_dodge(width = 0.9)) + 
  scale_fill_gradient2(high="#E64B35FF",low ="white",midpoint = 1 )+ theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()+
  ylim(0,20)+
  scale_x_discrete(limits=factor(data_new$Pathway))+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))

p1

ggsave(p1,file="./LUAD_GO_3.png",dpi = 1000,width = 6.72,height = 2.99)











data<-LUAD_GO[LUAD_GO$select==4,c("Description","GeneRatio","BgRatio","p.adjust","Count")]

colnames(data)
data

#LUAD_Go1
(11/274)/(126/18723)
(9/274)/(98/18723)
(5/274)/(36/18723)
(7/274)/(93/18723)
(4/274)/(25/18723)
(11/274)/(229/18723)
(7/274)/(99/18723)




data$`Enrich factor` <- c(5.965502,6.275398,9.490572,5.143278,10.93314,3.282329,4.831564)



data$`-log10(p.adjust)`=-log10(data$p.adjust)


data_new<-data[,c("Description","Count","Enrich factor","-log10(p.adjust)")]
colnames(data_new)<-c("Pathway","Count","Enrich_Factor","-log10(P_value)")

data_new
data_new<-data_new[order(-data_new$Enrich_Factor),]


p1<-ggplot() +  geom_bar(data = data_new, 
                         aes(y = Enrich_Factor, x = Pathway,fill=`-log10(P_value)`),
                         stat = "identity",
                         width = 0.7, 
                         position = position_dodge(width = 0.9)) + 
  scale_fill_gradient2(high="#E64B35FF",low ="white",midpoint = 1 )+ theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()+
  ylim(0,20)+
  scale_x_discrete(limits=factor(data_new$Pathway))+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))

p1

ggsave(p1,file="./LUAD_GO_4.png",dpi = 1000,width = 6.72,height = 2.99)



#NK cells

NK_MIA<-c(0.129415257,0.163586414,0.06505548,0.187907364,0.039215686,0.067276539)
NK_LUAD<-c(0.005903491,0.013770695,0.042505593,0.01820704,0.014424811,0.11786372,0.146372081,0.030333168,0.02328001,0.083934735)

wilcox.test(NK_MIA,NK_LUAD)







library(dplyr)
library(reshape2)
#计算大类的显著性
total_data<-read.table(file="/Users/wangjun/Desktop/最新版结果/大类比例计算显著性.csv",sep = ",",header = T,row.names = 1)
total_data[1:4,1:4]

data_melt<-melt(t(total_data))

head(data_melt)

colnames(data_melt)<-c("Sample","Cells","Proportion")


Tissue<-c("LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"LUAD",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"LUAD")
#LUAD	MIA	LUAD	MIA	LUAD	MIA	LUAD	MIA	LUAD	LUAD	LUAD	MIA	LUAD	MIA	LUAD	LUAD


Tissue_total<-rep(Tissue,10)
dim(data_melt)

data_melt$Tissue <- Tissue

unique(data_melt$Cells)

#总样本显著性检验
for (i in unique(data_melt$Cells)){
  new_data<-data_melt[data_melt$Cells==i,]
  p<-wilcox.test(new_data[new_data$Tissue=="MIA",]$Proportion,new_data[new_data$Tissue=="LUAD",]$Proportion)
  print(paste0(i,":",p$p.value))
}
data_melt[data_melt$Cells=="%B cells",]

data_melt[data_melt$Cells=="%tumor_inferCNV",]





colnames(total_data_new)

total_data_new<-total_data[,c(1,2,3,4,5,6,7,8,11,12,13,14)]

data_melt<-melt(t(total_data_new))

colnames(data_melt)<-c("Sample","Cells","Proportion")


Tissue<-c("LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA",	"LUAD",	"MIA")
#LUAD	MIA	LUAD	MIA	LUAD	MIA	LUAD	MIA	LUAD	LUAD	LUAD	MIA	LUAD	MIA	LUAD	LUAD

Tissue_total<-rep(Tissue,10)
dim(data_melt)

data_melt$Tissue <- Tissue


#配对样本样本显著性检验
for (i in unique(data_melt$Cells)){
  new_data<-data_melt[data_melt$Cells==i,]
  p<-wilcox.test(new_data[new_data$Tissue=="MIA",]$Proportion,new_data[new_data$Tissue=="LUAD",]$Proportion,paired = TRUE)
  print(paste0(i,":",p$p.value))
}
data_melt[data_melt$Cells=="%tumor_inferCNV",]



B_cells_data<-data_melt[data_melt$Cells=="%B cells",]

B_cells_data$Tissue<-factor(B_cells_data$Tissue,levels = c("MIA","LUAD"))


B_cells_new<-B_cells_data
p<-ggplot(data=B_cells_new)+
  geom_point(aes(x=factor(Tissue),y=Proportion,col=factor(Tissue)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.1)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/B_cells/total.png",dpi = 1000,width = 2.12,height = 3.44)

#data_melt[data_melt$Cells=="%tumor_inferCNV",]

#准备绘制NK细胞大类的比例图
NK_cells_data<-data_melt[data_melt$Cells=="%NK cells",]

NK_cells_data$Tissue<-factor(NK_cells_data$Tissue,levels = c("MIA","LUAD"))


NK_cells_new<-NK_cells_data[c(11,12),]
p<-ggplot(data=NK_cells_new)+
  geom_point(aes(x=factor(Tissue),y=Proportion,col=factor(Tissue)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.2)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/NK_cells/FD14.png",dpi = 1000,width = 2.12,height = 3.44)

1

NK_cells_data$


#出NK细胞和B细胞大类的箱线图
p<-ggplot(NK_cells_data,aes(x=`Tissue`,y=`Proportion`,color=`Tissue`))+geom_boxplot()

p1<-p+theme_bw()+theme(panel.grid = element_blank())+scale_color_manual(values=c("#4DBBD5FF","#E64B35FF"))+stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('MIA','LUAD')),method="wilcox.test")
p2<-p1+theme(axis.text.x=element_text(hjust = 0.5,size=12,face = "bold",color = "black"),
             axis.text.y=element_text(size=12,face = "bold",color = "black"),
             axis.title.y=element_text(size=12,face = "bold",color = "black"),
             axis.title.x = element_text(size=12,face = "bold",color = "black"),
             plot.title=element_text(hjust=0.5, size=12,face = "bold",color = "black"))+    
  labs(y="NK cells",x="Tissue",title="Boxplot")+   #添加标签
  #ylim(0,0.061)+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))#去除网格
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/NK_cells.png",dpi = 1000,width = 3,height = 4.54)






p<-ggplot(B_cells_data,aes(x=`Tissue`,y=`Proportion`,color=`Tissue`))+geom_boxplot()

p1<-p+theme_bw()+theme(panel.grid = element_blank())+scale_color_manual(values=c("#4DBBD5FF","#E64B35FF"))+stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('MIA','LUAD')),method="wilcox.test")
p2<-p1+theme(axis.text.x=element_text(hjust = 0.5,size=12,face = "bold",color = "black"),
             axis.text.y=element_text(size=12,face = "bold",color = "black"),
             axis.title.y=element_text(size=12,face = "bold",color = "black"),
             axis.title.x = element_text(size=12,face = "bold",color = "black"),
             plot.title=element_text(hjust=0.5, size=12,face = "bold",color = "black"))+    
  labs(y="B cells",x="Tissue",title="Boxplot")+   #添加标签
  #ylim(0,0.061)+
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black",size = 1))#去除网格
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/B_cells.png",dpi = 1000,width = 3,height = 4.54)













#T细胞拟时序分析
library(monocle)

load(file="/Users/wangjun/Desktop/monocle_t/CDS_1000.Rdata")
#load(file="/Users/wangjun/Desktop/monocle/CDS_2000.Rdata")

CDS$Cluster
p3<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+facet_wrap("~Cluster", nrow = 1)+
  #scale_color_manual(values=pal[c(1,2,5)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 6),axis.title.x=element_text(size=6,face = "bold"),axis.title.y=element_text(size=6,face = "bold"),axis.text.x=element_text(size=6,face = "bold",color = "black"),axis.text.y=element_text(size=6,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 6))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p3

ggsave(p3,file="./monocle_Type.png")
ggsave(p3,file="./monocle_Type.png",dpi = 1000,width = 5.02,height =3.43)


#T细胞的拟时序没有得到有价值的结果，舍去



#计算TCR最高丰度

TCR_MIA<-c(0.153027823,0.046867959,0.027217269,0.018162948,0.038358608,0.009251472)
TCR_LUAD<-c(0.009866512,0.046786922,0.01793722,0.027895182,0.02237927,0.007578628,0.012586277,0.041598248,0.008250825,0.010383747)

wilcox.test(TCR_LUAD,TCR_MIA)



TCR_MIA<-c(0.153027823,0.046867959,0.027217269,0.018162948,0.038358608,0.009251472)
TCR_LUAD<-c(0.009866512,0.046786922,0.01793722,0.027895182,0.02237927,0.007578628)

wilcox.test(TCR_LUAD,TCR_MIA,paired=TRUE)


TCR<-data.frame(matrix(NA,16,2))
colnames(TCR)<-c("Samples","Index")
TCR[1:6,]$Samples="MIA"
TCR[7:16,]$Samples="LUAD"


TCR[1:6,]$Index=TCR_MIA
TCR[7:16,]$Index=TCR_LUAD

TCR$Samples<-factor(TCR$Samples,levels = c("MIA","LUAD"))


TCR_new<-TCR
p<-ggplot(data=TCR_new)+
  geom_point(aes(x=factor(Samples),y=Index,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  ylim(0,0.08)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/TCR/最高丰度/FD14.png",dpi = 1000,width = 1.87,height = 3.41)








a<-c(6.343318,7.62522,7.236271,6.890004,7.587322)
b<-c(6.294129,7.255371,6.342968,6.173589,7.173777)


a_updata<-c(5.214451,6.890004,6.890004,7.587322,7.066926,7.62522,6.343318,7.403371,7.236271)
b_updata<-c(7.127372,5.950745,6.173589,7.173777,6.719225,7.255371,6.294129,6.730217,6.342968)


tumor_a<-c(1.796703,1.990221,2.037571,1.689583,1.876247,1.888009,1.674485,1.933822)
tumor_b<-c(1.859713,1.763369,1.825259,1.605069,1.660304,1.448939,1.600487,1.656552)


shapiro.test(tumor_a)


m=wilcox.test(tumor_a,tumor_b,paired = TRUE)


t.test(tumor_a,tumor_b,paired = TRUE)


TCR<-data.frame(matrix(NA,16,2))
colnames(TCR)<-c("Samples","Index")
TCR[1:8,]$Samples="High"
TCR[9:16,]$Samples="Low"


TCR[1:8,]$Index=tumor_a
TCR[9:16,]$Index=tumor_b

TCR$Samples<-factor(TCR$Samples,levels = c("Low","High"))

TCR


TCR_new<-TCR[c(8,16),]
p<-ggplot(data=TCR_new)+
  geom_point(aes(x=factor(Samples),y=Index,col=factor(Samples)),size=2)+
  xlab("TCR diversity")+
  ylab("Tumor heterogeneity")+
  ylim(1.4,2.2)+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/TCR/肿瘤异质性相关/8_updata.png",dpi = 1000,width = 2.12,height = 3.41)




#肿瘤异质性计算
T1<-c(18.88782359,18.12080537,15.86768936,12.27229147,15.77181208,7.046979866,5.800575264,6.232023011)

T2<-c(1.350521793,26.21240025,13.44383057,0.368324125,23.38858195,14.48741559,0.061387354,10.62001228,6.077348066)

T3<-c(1.516070346,18.07155852,16.19163129,21.10369921,17.04063069,12.49241965,5.518496058)
T4<-c(22.61846605,16.41426478,19.73619932,16.75622863,5.178309722,0.683927699,7.52320469,7.425500733)

T5<-c(13.51616063,5.093046033,15.08325171,21.59647405,26.19980411,5.82761998,8.178256611)


T6<-c(21.78304787,16.07378129,6.587615283,9.617918314,0.966183575,10.8476065,7.158541941,10.32059728,16.64470795)

T7<-c(22.57231405,13.68801653,28.20247934,10.95041322,10.12396694,5.681818182,1.033057851,5.268595041)


T8<-c(36.31653321,8.101742817,8.054639661,20.63118229,8.431464908,12.19971738)

T9<-c(18.25264329,20.53422371,15.91541458,5.564830273,21.48024485,7.011686144)



T10<-c(49.52647439,10.69737409,9.384416703,9.448988377,8.22212656,6.823073612)

T11<-c(12.57823541,9.428629114,7.651928124,8.015344236,0.444175247,18.61498082,15.12214819,23.44033919)



T12<-c(21.4532872,18.12283737,15.1816609,6.314878893,14.57612457,7.525951557,14.10034602)


T13<-c(38.03977273,26.10795455,6.761363636,5.340909091,0.511363636,6.448863636,8.4375,6.590909091)



T14<-c(28.56137883,21.10446711,9.004572635,10.83362645,16.95392191,6.718255364)

T15<-c(36.84798808,14.67958271,20.71535022,5.923994039,11.10283159,7.414307004)


T16<-c(23.80405664,9.911978569,12.05510907,13.89207807,0.727133563,11.40451588,13.54764638,11.97856869)

T17<-c(15.94344222,13.01803998,0.048756704,28.03510483,7.069722087,22.9156509,5.655777669)


diversity(T3, index = "shannon")






#拟时序更新

library(monocle)

#load(file="/Users/wangjun/Desktop/monocle_1.7/CDS_1300.Rdata")
load(file="/Users/wangjun/Desktop/monocle/CDS_2000.Rdata")

#table(Idents(CDS))

CDS

p3<-plot_cell_trajectory(CDS, color_by = "Type", cell_size = 0.5)+facet_wrap("~Cell_type", nrow = 1)+#scale_color_manual(values=pal[c(1,2,5)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 6),axis.title.x=element_text(size=6,face = "bold"),axis.title.y=element_text(size=6,face = "bold"),axis.text.x=element_text(size=6,face = "bold",color = "black"),axis.text.y=element_text(size=6,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 6))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p3

table(CDS)


table(CDS$Type[CDS$State==2])


p4<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+facet_wrap("~Cell_type", nrow = 1)+#scale_color_manual(values=pal[c(1,2,5)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 6),axis.title.x=element_text(size=6,face = "bold"),axis.title.y=element_text(size=6,face = "bold"),axis.text.x=element_text(size=6,face = "bold",color = "black"),axis.text.y=element_text(size=6,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 6))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p4

table(CDS$State[CDS$Cell_type=="Tumor cells"])


CDS_tumor<-CDS[CDS$Cell_type=="Tumor cells"]

table(CDS_tumor$Cell_type)

table(CDS_tumor$Cell_type[CDS_tumor$Cell_type=="Tumor cells"])


CDS@reducedDimS[,1:4]

plotdf2=as.data.frame(t(CDS@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=CDS$Pseudotime
plotdf2$State=CDS$State
plotdf2$Cell_type=CDS$Cell_type
plotdf2$Type=CDS$Type


plotdf2_new<-plotdf2[plotdf2$Cell_type=="Tumor cells",]
dim(plotdf2_new)

table(plotdf2_new$State)
plotdf2_new_luad<-plotdf2_new[plotdf2_new$Type=="LUAD",]
plotdf2_new_mia<-plotdf2_new[plotdf2_new$Type=="MIA",]

dim(plotdf2_new)


library(dplyr)
plotdf2_new_mia%>%ggplot(aes(component1,component2,color=State))+
  geom_point()+
  xlim(-10,5)+
  ylim(-10,5)+
  theme_minimal()

table(plotdf2_new_luad$State)














pal[c(9,2,7)]

pal[1:10]
pal[c(9,2,7)]
#00A087FF #4DBBD5FF #B09C85FF #91D1C2FF #3C5488FF #8491B4FF #F39B7FFF #7E6148FF #E64B35FF 

plot_cell_trajectory(CDS, color_by = "Type", cell_size = 1)+scale_color_manual(values=c("#0F99B2","#E36666","#F2B77C","#66CCFF"))+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

CDS$Type

table(CDS$orig.ident)




p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p2<-plot_cell_trajectory(CDS, color_by = "Pseudotime", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)



p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5)+
  scale_color_manual(values=pal[c(1,2,5)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p2<-plot_cell_trajectory(CDS, color_by = "Pseudotime", cell_size = 0.5)+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title.x=element_text(size=14,face = "bold"),axis.title.y=element_text(size=14,face = "bold"),axis.text.x=element_text(size=14,face = "bold",color = "black"),axis.text.y=element_text(size=14,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 20))+NoLegend()
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)

p4<-p1+p2
p4
ggsave(p4,file="./monocle_Pseudotime.png")
ggsave(p4,file="./monocle_Pseudotime_4.png",dpi = 1000,width = 6.24,height =3.04)


p3<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 0.5,show_branch_points = F)+facet_wrap("~orig.ident", nrow = 1)+
  scale_color_manual(values=pal[c(1,5,2)])+ 
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title.x=element_text(size=5,face = "bold"),axis.title.y=element_text(size=5,face = "bold"),axis.text.x=element_text(size=5,face = "bold",color = "black"),axis.text.y=element_text(size=5,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 5))
#p1<-plot_cell_trajectory(CDS, color_by = "State", cell_size = 1)
p3

ggsave(p3,file="./monocle_Pseudotime_orig.ident.png",dpi = 1000,width = 10.5,height =2.18)



table(CDS$Type[CDS$orig.ident!="T3" & CDS$State==3 ])



#肿瘤细胞比上正常细胞
table(Idents(Epi_total))

Epi_TN<-subset(Epi_total,idents=c("AT1 cells","AT2 cells","Tumor cells"))

new.cluster.ids <- c("Normal cells","Normal cells","Tumor cells")
names(new.cluster.ids) <- levels(Epi_TN)
Epi_TN <- RenameIdents(Epi_TN, new.cluster.ids)

table(Idents(Epi_TN))

T_N<-Idents(Epi_TN)
Epi_TN<-AddMetaData(Epi_TN,T_N,col.name = "T_N")

Idents(Epi_TN)<-"orig.ident"

library(dplyr)

T1<-subset(Epi_TN,idents="T1")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T1_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T2")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T2_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T3")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T3_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T4")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T4_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T5")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T5_DEG_gene.csv")


T1<-subset(Epi_TN,idents="T6")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T6_DEG_gene.csv")


T1<-subset(Epi_TN,idents="T7")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T7_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T8")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T8_DEG_gene.csv")




T1<-subset(Epi_TN,idents="T9")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T9_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T10")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T10_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T11")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T11_DEG_gene.csv")


T1<-subset(Epi_TN,idents="T12")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T12_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T13")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T13_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T14")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T14_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T15")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T15_DEG_gene.csv")


T1<-subset(Epi_TN,idents="T16")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T16_DEG_gene.csv")



T1<-subset(Epi_TN,idents="T17")
Idents(T1)<-"T_N"
scRNA.markers <- FindAllMarkers(object = T1, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)
marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=1) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)
write.csv(marker_gene,"/Users/wangjun/Desktop/最终方案/figure/病理分期差异基因/DEG_new/T17_DEG_gene.csv")









#
macrophage<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/macrophages_update.rds")

table(Idents(macrophage))

table(subset(macrophage,idents="Monocytes Derived DC")$orig.ident)
table(subset(macrophage,idents="Alveolar resident Macrophage")$orig.ident)
table(subset(macrophage,idents="Perivascular resident Macrophage")$orig.ident)
table(subset(macrophage,idents="Classical CD14+ Monocyte")$orig.ident)
table(subset(macrophage,idents="Anti-inflammatory macrophages")$orig.ident)
table(subset(macrophage,idents="Non-classical CD16+ Monocyte")$orig.ident)
table(subset(macrophage,idents="Conventional Dendritic cell type 2")$orig.ident)
table(subset(macrophage,idents="Proliferating Macrophage")$orig.ident)
table(subset(macrophage,idents="Migratory Conventional Dendritic cell")$orig.ident)
table(subset(macrophage,idents="Conventional Dendritic cell type 1")$orig.ident)

table(macrophage$orig.ident)


macrophage_indenty<-read.table(file="/Users/wangjun/Desktop/macrophage_percentage.csv",sep = ",",header = T,row.names = 1)

macrophage_indenty[1:4,1:4]

for (i in 1:length(rownames(macrophage_indenty))){
#for (i in 6:6){
  LUAD=macrophage_indenty[i,c("T1","T5","T6","T8","T10","T12","T13","T14","T16","T17")]
  MIA=macrophage_indenty[i,c("T2","T4","T7","T9","T11","T15")]
  LUAD_pair=macrophage_indenty[i,c("T1","T5","T6","T8","T10","T14")]
  P=wilcox.test(as.numeric(LUAD),as.numeric(MIA))$p.value
  print(paste0(rownames(macrophage_indenty)[i],":",P))
  #P=wilcox.test(as.numeric(LUAD_pair),as.numeric(MIA),paired = T)$p.value
  #print(paste0(rownames(macrophage_indenty)[i],":",P))
}





MIA
LUAD

wilcox.test(LUAD,MIA)



#不同分期细胞差异基因计算
#导入数据
library(Seurat)
library(dplyr)
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/B_cells_11_15.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/B_cells_DEG_gene.csv")








#计算T细胞的差异基因
T_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/T_total_12_10.rds")



table(Idents(T_cells))

B_cells<-subset(T_cells,idents = c("CD4_C1","CD4_C2","CD4_C3","CD4_C4","CD4_C5","CD4_C6","CD4_C7"))

B_cells<-subset(T_cells,idents = c("CD8_C1","CD8_C2","CD8_C3","CD8_C4","CD8_C5","CD8_C6","CD8_C7"))



table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/CD8_T_cells_DEG_gene.csv")







#计算肥大细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/Mast_cells.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/Mast_cells_DEG_gene.csv")










#计算成纤维细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/Fibroblast.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/Fibroblast_cells_DEG_gene.csv")








#计算内皮细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/Endothelial.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/Endothelial_cells_DEG_gene.csv")







#计算NK细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/NK_cells_update.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/NK_cells_DEG_gene.csv")










#计算Plasma细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/Plasma_cells.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/Plasma_cells_DEG_gene.csv")










#计算AT_cluster细胞差异
B_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/AT_cluster.rds")

table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/AT_cells_DEG_gene.csv")






#计算macrophages细胞差异
macrophages_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/macrophages_update.rds")

table(Idents(macrophages_cells))


B_cells<-subset(macrophages_cells,idents = c("Alveolar resident Macrophage","Perivascular resident Macrophage","Anti-inflammatory macrophages","Proliferating Macrophage"))


B_cells<-subset(macrophages_cells,idents = c("Monocytes Derived DC","Conventional Dendritic cell type 2","Migratory Conventional Dendritic cell","Conventional Dendritic cell type 1"))

B_cells<-subset(macrophages_cells,idents = c("Classical CD14+ Monocyte","Non-classical CD16+ Monocyte"))




table(Idents(B_cells))
table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/DC_cells_DEG_gene.csv")









#计算T细胞的差异基因
T_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/T_total_12_10.rds")



table(Idents(T_cells))

B_cells<-subset(T_cells,idents = c("CD4_C2"))



table(B_cells$orig.ident)

Idents(B_cells)<-"orig.ident"

B_cells_ture<-subset(B_cells,idents=c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17"))

table(Idents(B_cells_ture))

B_cells_ture$Type = "LUAD"
B_cells_ture$Type[B_cells_ture$orig.ident %in% c("T2","T4","T7","T9","T11","T15")] = "MIA"
table(B_cells_ture$Type)


Idents(B_cells_ture)<-"Type"

scRNA.markers <- FindAllMarkers(object = B_cells_ture, only.pos = TRUE, 
                                min.pct = 0.25, thresh.use = 0.25)

marker_gene <- scRNA.markers %>% filter(p_val_adj <= 0.01 & avg_log2FC>=0.5) %>% group_by(cluster) %>% top_n(n = 1000, wt = p_val_adj)

write.csv(scRNA.markers,"/Users/wangjun/Desktop/差异基因计算/CD4_C2_T_cells_DEG_gene.csv")

table(Idents(B_cells_ture))







#计算拟时序LUAD分支的通路富集
data<-read.table(file="/Users/wangjun/Desktop/最新版结果/肿瘤细胞拟时序/tumor_CDS_DEG_gene_update.csv",header = T,sep = ",")

data[1:4,]

DEG<-data[data$cluster==2,]$gene

#差异基因
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(ggplot2)
library(org.Hs.eg.db)


geneList1<-bitr(unique(DEG), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)

head(geneList1)


GO2<-enrichGO( geneList1$ENTREZID,#GO富集分析
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 0.05,#设定p值阈值
               #qvalueCutoff = 0.05,#设定q值阈值
               pAdjustMethod = "BH",
               readable = TRUE)

kk2<-enrichKEGG( geneList1$ENTREZID,#GO富集分析
                 organism = "hsa",
                 keyType = "kegg",#设定读取的gene ID类型
                 #ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                 pvalueCutoff = 0.05,#设定p值阈值
                 #qvalueCutoff = 0.05,#设定q值阈值
                 pAdjustMethod = "BH")


GO2@result[1:8]

kk2

write.table(GO2@result,file="/Users/wangjun/Desktop/最新版结果/肿瘤细胞拟时序/GO2_down.csv",sep = ",",col.names = T)



#WES筛选出的Silent基因表达量
library(Seurat)
#导入肿瘤细胞数据
data<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/tumor_infercnv.rds")

data<-readRDS(file="/Users/wangjun/Desktop/pbmc_final.rds")
DefaultAssay(data)<-"RNA"


Idents(data)<-"orig.ident"
table(Idents(data))



FD1_LC2_silent=c("CEP131","OR11H12","SHISA2","CHST5","IDH3A","RNF213","CROCC")
FD1_LC2_data<-subset(data,idents = "T4")


data_FD1_LC2_silent<-rowMeans(data.frame(GetAssayData(FD1_LC2_data))[FD1_LC2_silent,])
mean(data_FD1_LC2_silent[!is.na(data_FD1_LC2_silent)])



FD1_LC2_silent_no=c("OR2L3","GGT2","AZIN2","MAPKBP1","KAT2A","ALPP","SLAMF9","DRICH1","ANKLE1","COLQ","MEX3D","ATAD3B","RBM47","MUC12","RPL21","CDC25B","IL17D","ZXDA","PDSS1","UNC93B1","SHROOM4","GNAQ","WIZ","SLC34A3","DUSP5","CCDC61","RIOK1")
FD1_LC2_data<-subset(data,idents = "T4")

data_FD1_LC2_silent<-rowMeans(data.frame(GetAssayData(FD1_LC2_data))[FD1_LC2_silent_no,])
mean(data_FD1_LC2_silent[!is.na(data_FD1_LC2_silent)])







FD1_LC3_silent=c("GPD2","DDX11","NBPF10","TENT4A","GSR","RASIP1","GON4L","PLEC","WASHC1","HOXC13","FRG2C","CASP10","OGFR","APC","SSC5D","CRY1","SQLE","FLVCR1")
FD1_LC3_data<-subset(data,idents = "T5")

data_FD1_LC3_silent<-rowMeans(data.frame(GetAssayData(FD1_LC3_data))[FD1_LC3_silent,])
mean(data_FD1_LC3_silent[!is.na(data_FD1_LC3_silent)])



FD1_LC3_silent_no=c("EPS15","SEC24C","CAPN14","NEUROD2","ACSF2","FBXO41","RBM10","COL15A1","NEFH","MAP7","BRPF1","MUC12","PNMT","C16orf90","ZBTB48","FRG2C","TRIOBP","MS4A14","CYP2C9","OR10K2","TNPO3","NLRX1","COL16A1","SEMA6A","PPIL2","TEKT3","ARMCX4","TRIM25","GPR62","ALPP","MID1","WIZ","AQP11","CACNA2D2","ZDHHC19","FLNA","ZFHX4","OR5F1","ZXDA","BMP6","B4GALT5","SYMPK","TRO","MME","RAB11FIP4","USHBP1","BAP1","RGR","CRY1","CDH7","JAML","KLHL40","LCT","EGFR","MAP3K5","ZNF846")
FD1_LC3_data<-subset(data,idents = "T5")

data_FD1_LC3_silent<-rowMeans(data.frame(GetAssayData(FD1_LC3_data))[FD1_LC3_silent_no,])
mean(data_FD1_LC3_silent[!is.na(data_FD1_LC3_silent)])



#FD2

FD2_LC1_silent=c("RTN1","PPM1H","SYT6","TRIM65","PCDHGA3","AR","LEPROTL1","MRGBP","ABCA1","ZNF845","C7orf26","ZPBP2","GRIN2A","BIRC5","TEX13C","KANK1","DDX11","ATAD2","SLC13A3","NXNL1","MRPL38","TBXAS1","RHOT1","LOC400499","KMT5A","BOP1","NBPF20","ASAH2")
FD2_LC1_data<-subset(data,idents = "T1")

data_FD2_LC1_silent<-rowMeans(data.frame(GetAssayData(FD2_LC1_data))[FD2_LC1_silent,])
mean(data_FD2_LC1_silent[!is.na(data_FD2_LC1_silent)])


FD2_LC1_silent_no=c("SVEP1","PAPSS2","ATP13A1","RENBP","WIPF3","ZNF8","ZNF257","LAMA1","ANKRD36","ARHGAP10","SHANK2","NKX2-3","MAMLD1","MAST3","POTEH","POU3F2","KDM3B","EPC2","TRPV1","PI4KA","CNTLN","CCSER2","DEK","ADAM30","ITCH","SULT1A2","NAV2","EHD1","SECTM1","REEP1","PTX4","GIN1","KMT5A","SBNO1","GNAS","VSIG4","MED12L","MEGF10","ZNF329","KCNH4","NAP1L1","ADCY8","GRIN2A","MDM2","BRPF1","NLRP9","FOXJ2","POLL","KIF13A","CELF2","NCOR2","LRRN1","CSAD","ANP32C","HCN1","VEGFC","ZDHHC8","NXNL1","FCRL3","AARD","MROH2A","PIEZO1","AGTPBP1","APBB1","FKTN","DSG1","MYBBP1A","USH2A","LMTK3","GFPT2","APC2","TDRD1","HSPA12A","NBR1","MAP3K15","CNGB1","OR8H2","HAPLN4","ACSL4","DACT1","MRPS34","DOLPP1","ARID3B")
FD2_LC1_data<-subset(data,idents = "T1")

data_FD2_LC1_silent<-rowMeans(data.frame(GetAssayData(FD2_LC1_data))[FD2_LC1_silent_no,])
mean(data_FD2_LC1_silent[!is.na(data_FD2_LC1_silent)])






FD2_LC2_silent=c("SLC12A2","PLIN4","ZNF845","SLC2A3","STAB1","FOXC1","KCNA6","PCDHB9","HSPG2","ALDH3B1","NUFIP1","MLLT3","SOX1")
FD2_LC2_data<-subset(data,idents = "T2")

data_FD2_LC2_silent<-rowMeans(data.frame(GetAssayData(FD2_LC2_data))[FD2_LC2_silent,])
mean(data_FD2_LC2_silent[!is.na(data_FD2_LC2_silent)])



FD2_LC2_silent_no=c("TJP3","FOXD4","GSC2","RSPH6A","LRRCC1","MTCL1","MYL7","ANKRD36","C1QL1","MTHFD2L","KCNMA1","NACAD","TENM2","ZC3H3","ZNF831","CSRNP1","KLF12","SYNE3","TRIM67","KMT5A","ITGA8","THSD7A","PRDM9","DEAF1","UTF1","SP8")
FD2_LC2_data<-subset(data,idents = "T2")

data_FD2_LC2_silent<-rowMeans(data.frame(GetAssayData(FD2_LC2_data))[FD2_LC2_silent_no,])
mean(data_FD2_LC2_silent[!is.na(data_FD2_LC2_silent)])



#FD4


FD4_LC1_silent=c("EGFR","SLC28A1","C1QTNF1","NWD2","ZNF587","FPR2","GJD2","SLC43A1","ADAM19","SPTA1","SMOX","TLX2","PKD1L3","RBPJL","APLNR","H3C10","C8B","ARAP3","ALDH3B1","INPP5A","XIRP1","NIPAL4","ERN1","NDUFS2","FUT7","D2HGDH","DST","AQP12B","APOB","GPR87","CDH18","ASXL3","PPFIA3","HMCN2","RETN","MAGEB5","RAX2","PHYHD1","CDRT1","KIRREL2","ZNF467","KRTAP24-1","KMT2C","RERG","SPATA31D1","ANK2","TUT4","TMEM132B","ATP8B4","RIPOR1","EML5","LRCOL1","SMC3","BPIFB1","SLC9B1","IL20","RUNDC3B","RANBP2","SIMC1","SHD","NUTM1","GRIA3","FREM3","NFIX","MIR1-1HG","ACTL8","ARVCF","COL7A1","CDH22","POTEH","NALCN","TRMT6","IL36B","KCND2","TSGA10IP","TOMM40L","FADS6","UTRN","BANK1","ISLR2","PAX2","OR5M9","PEAR1","RBM24","KCNH2","ITGA11","CYP4F3","CDX2","FAM186A","OLFM1","TBC1D22A","OBSCN","FOXD4L6","SLC4A2","CDK5RAP3","CNTN6","MAP3K1","CD248","MYO7B","SLC45A1","PCDHAC1","COL20A1","XAB2","RBM20","UBAP2L","SV2C","MED12","KSR1","TMEM63A","OR8H2","SH2D7","KLHL1","UNC80","ADAMTS4","EPHA10","KMT2B","TMEM14B","TGM5","AOC2","CACNB1","TRMT10C","HOXA1","SRRM2","CYBB","ARMC3","MYLK","RYR1","ELOA2","OR1A1","TANGO6","PMPCA","VCAN","TNRC18","SLC7A10","SLC6A7","ZC3H11B","ACAN","CTSV","HK3","SALL4","HMGB3","ADARB2","HTT","LCP2","C3orf20","EDA2R","CAMK2B","GNAQ","SHPK","SPATA22","CCDC191","EPS8L3","PLCB1","LOC100996413","API5","OR9I1","FREM2","DDX11","ADPGK","TCEAL3","C11orf95","CASR")
FD4_LC1_data<-subset(data,idents = "T6")

data_FD4_LC1_silent<-rowMeans(data.frame(GetAssayData(FD4_LC1_data))[FD4_LC1_silent,])
mean(data_FD4_LC1_silent[!is.na(data_FD4_LC1_silent)])


FD4_LC1_silent_no=c("GCOM1","EIF4G3","MTTP","KIF1A","ASCL4","FBN2","GRIN2A","RBM33","SCP2","LRFN5","ARHGEF15","ITSN1","NRP1","TPST2","PRUNE2","OR51S1","HSPB2","ALDH1A2","KLRG1","CX3CR1","TTC14","IGDCC3","DZIP1L","GIMAP6","HPD","TF","PCDH15","PKP2","CCNB3","MUC7","NTN1","ZNF599","MCM3AP","KRTAP22-1","MPL","DNAH5","ZNF512B","PCDHA7","SPATA31D1","IFT140","PCDH7","ZNF75A","KCNG2","SLFN11","TP53AIP1","PITPNM3","C10orf90","POU5F2","SPATA31C1","MAGEA11","PROX1","TCHHL1","CACNA2D3","RIOK1","PDPR","NCAM1","OR4D10","PCDHGB5","SHQ1","CACNA2D1","NPAP1","GFRA1","SLC26A1","PYGM","PLPPR4","ZDHHC11","FRMPD1","ANKRD17","IGSF9","NSUN2","CYLC2","IL25","MUC17","MAP1A","HS6ST1","CGB2","PCDHB3","BTNL3","TLL2","DRD3","NRK","GABRG2","LMAN2","ANKDD1A","PIK3CA","ZNF423","AR","DEK","DNER","CCDC8","CCDC57","ZNF76","NUMA1","ALDH4A1","CR1","SLC7A4","LNPEP","SIGLECL1","PCED1A","TDRD9","PSD2","CDH11","MATK","KIF27","CD160","APCS","KALRN","NINL","MCIDAS","TSHZ3","BUB1B-PAK6","NCKAP5L","CYFIP2","NLRP1","RPS12","MAGI1","CHRNA3","KRTAP4-7","GH2","AHSG","ZNF713","TBC1D4","SHC4","DNAH8","TNFAIP2","CDH4","MMP8","HES6","KCNB1","KRT6A","MN1","WDR6","NLRP3","MUC12","MBLAC2","COL6A6","FUT6","ANK1","ARID5B","GTF2E1","CHRM5","LRRIQ3","RERG","PMS2","CROCC","PAWR","ALMS1","MRC1","ADGRB2","KCP","FCRLA","SOAT2","ARID3A","CNOT8","FBXO30","CHD7","MICAL3","GATA5","EPHA7","UBE3A","MYCL","CDH20","MYH2","GSC","HIRA","GOLGA6B","ANO2","CDA","RSPH9","ELSPBP1","TMPRSS7","PCDHB9","ZNF415","CYP4F22","RIMBP2","MCAT","CNN1","KRAS","RBM10","NPPB","IL7R","ZNF485","DST","ZDHHC14","CAPN10","POTEA","CCDC171","ARFGEF2","TRIML2","KIF19","STYXL2","GFAP","SORCS1","CLIP1","COL20A1","TTC17","PLEKHH1","PRR23C","P3H2","ADCYAP1R1","PPP1R37","RBP2","FAM221B","ZFHX3","UBFD1","KCNK2","LGALS8","GMNC","PNPLA6","FLNB","PRDM8","LCE3D","NLRP8","DNAH6","SNX15","SLC28A3","ARHGEF9","MAP3K19","TSC22D1","RYR2","WFDC5","CYP2C18","SH3GLB1","OR51L1","FAT3","PRODH2","LRRK2","OR2C1","ARHGAP6","TSPAN16","PTGS1","GAREM2","KCNIP1","H6PD","SPTBN4","CHIA","RAF1","ALOX5AP","CEACAM6","PLA2R1","CCDC168","ADGRA1","PARG","PABPC3","CWF19L2","KRTAP15-1","KCNH5","SMAD9","OR5L1","C4orf54","DCLK1","FAM71B","FAM171A2","VAT1L","TRPM6","C16orf54","NUP210","THSD7B","TAS2R1","CYP11B1","USP29","LYN","UTF1","ARMC8","PCDHGB2","HEPH","KCNH6","KIAA1671","CELF6","KHDRBS1","ZNF148","SYT6","ITIH1","ZFHX4","UGT1A4","PGLYRP3","TMEM117","PIF1","DNTTIP2","PCSK7","CPXM2","WWP1","IMPG1","FAM166C","LCMT2","KRTAP10-12","ADGRG4","KCNH7","NAV3","TEPSIN","USP26","PUS1","PAX4","GCSAML","PCDHA10","MAP3K1","KIAA0040","OR52E4","ADAMTSL1","ZNF483","AHNAK","TBX6","ALX4","COL5A3","KCNAB1","FAM186A","PCDHB7","FOXN4","STAG1","SERPINB12","PARP8","SEMA3E","NBPF11","TNRC6B","ARHGEF4","LMLN2","MYLK4","HOXA1","ALPG","UNC5CL","GCN1","POT1","PELI2","RPF2","KCNN1","C2CD4A","NLRP12","FOXD4L3","PTPRN2","DACH1","ZC3H13","GJA8","CDX1","SLC17A3","GLDC","FERD3L","OR2T8","CAMSAP2","SEC31A","NAP1L3","OR10A7","MAEL","ADGRF4","PTPRB","PCDH11X","C3orf56","MPPED1","UPK1A","GDF2","CCBE1","DOCK9","OR10G7","INPP5A","FRMD5","ECEL1","SAGE1","TBC1D16","OR1A1","LOXHD1","FAM43B","FBN1","FMO3","MASP1","SETBP1","DHX33","ARHGEF38","GPR101","SYNE1","AEBP1","OR2W5","GPR78","SRC","EFCAB5","BTAF1","BTBD16","MAGEB2","POLQ","GPR87","TSGA10IP","TEX35","FMNL3","EYS","PCDHB6","TG","PTPN23","FXYD5","DCSTAMP","RIN3","PDILT","GAB4","PALB2","GLI3","KRT6B","DNAH3","ABCB1","VEZT","ZNF280D","WDR97","MYBPC3","CCSER1","DAB1","ZFPM2","CLRN2","CHD5","PABPC1","CUBN","AJAP1","DMP1","TNKS2","MYO18B","DNAL1","DCAF13","FEM1B","LRRC53","ITGA7","DGKB","COL3A1","ZEB2","CLTB","SLC12A1","OTOGL","TNRC18","FAM151A","CDH6","APC2","CD163L1","ATXN1L","EPPK1","OTOG","FNDC1","PKD2L1","KBTBD3","SYPL1","KCTD6","ITGA4","ESR1","QARS1","SPPL3","COL22A1","COPB2","PIGT","CACNA1I","GABRG3","RBM20","AHI1","SETD2","MMRN1","MAS1","DMRTB1","OR8K5","KRT76","WDFY4","OR5M1","TM9SF3","XPNPEP1","USH2A","SLC9A9","UNC80","IZUMO3","OLFML2A","COL1A2","PPRC1","EXOC3L2","CBLC","DBNL","NRXN1","GLI2","LOC100509620","COL4A1","PDE10A","MXRA5","OR4C6","MED23","ZNF226","SETD1B","CSMD2","KRTAP5-10","MROH9","ANKS3","NOLC1","ANKRD36","GRIP1","BSPRY","LINGO2","ABCA4","CKMT1B","EXPH5","CCDC40","MMP2","RYR3","TMEM14B","SLC52A3","EMILIN1","OR6K6","OR7E24","NAALADL2","RIMS1","DOCK2","SLC5A12","RAI1","NWD2","PRC1","ZNF675","RTKN2","IQCF5","CACUL1","CNBD1","ZNF114","MST1R","OR10A2","PCDHB4","SLFN13","KCNH1","ELMO1","F5","FREM1","P2RX7","CNNM1","CACNA1E","NBPF9","PCDH12","ODAM","LHX9","WIZ","SMCHD1","CHGB","SHISA9","ZNF514","THRA","MFSD4A","ZNF236","CBY2","ZSCAN32","PRKDC","MEGF8","NR2F1","SF3B4","OR2AK2")
FD4_LC1_data<-subset(data,idents = "T6")

data_FD4_LC1_silent<-rowMeans(data.frame(GetAssayData(FD4_LC1_data))[FD4_LC1_silent_no,])
mean(data_FD4_LC1_silent[!is.na(data_FD4_LC1_silent)])





FD4_LC2_silent=c("OR2T33","STRIP2","VWA5B2","ITPRID1","NKX6-1","OR52A1","OR8G2P","ADCY7","CDH23","NOTUM","GAPVD1","ACTR1B","PIAS2","OR5AR1","UBC","MYH1","KLHL3","YIPF7","HOXA9","AKAP6","GIMAP1","AJAP1","ISM1","LOC101928841","ANKRD52","STX7","TRPC7","FER1L6","RASAL3","RP1L1","SIGLEC8","DRD5","TAS1R3","HSPA12A","SMC3","WASHC2A","FAM98B","CLSTN2","ZNF606","CCDC89","YBX1","RRAS","ZNF711","PML","VSTM4","KMT5A","RANBP2","NCLN","DNAH5","RAB3IL1","NYAP1","TMEM255A","EPHB3","LRRC38","ADAMTS20","NMUR1","GFPT2","PNLIP","FMOD","CLK4","GNAQ","TCEAL3","NBR1","PIEZO2","OR56A3","LGALS12","CHMP1A","RGS9","RHOBTB1","ANKK1","COL3A1","ANKRD35","HNRNPA0","UNC45B","PYHIN1","OR8G1","WNT10B","PCDH18","PCSK5","CCDC60","CPS1","IPO9","MROH9","POM121L2","AFF2","TNFAIP2","MYH4","FGF13","HGF","PCDHA6","PCDH15","EXOC3","TM7SF3","STAR","NR4A2","SPEN","BTBD11","ZNF544","BHLHE22","MAGI2","P2RX2","TBX15","OBSCN","BTNL3","HOXB1","VPS33A","MGAT5B","CILP","KCNJ12","ADRB3","C1orf167","CD163L1","ALG13","DOK3","PITX1","NUFIP1","BAHCC1","SLC14A2","SLC5A10","POP4","L2HGDH","SKI","KIAA1109","MXRA5","KMT2C","ROBO3","CUX1","HOXD10","NSD1","SIRPD","GLT1D1","KLK11","AMBRA1","UMOD","CFAP46","KRT35","FCRL6","STAB2","SYT1","TTC9B","NSUN5","ASXL2","TRIM52","NALCN","FGB","UGT1A4","HRNR","IGFN1","GTSE1","DSG3","ZNF81","MAP1B","ATP6AP1L","DISP3","OR5AS1","PCNT","ITGBL1","IFI16","CAPN10","PPP4R3B","CYP11B1","POTEJ","TRPM4","SEMA3F","SLC12A5","STARD8","EXO1")
FD4_LC2_data<-subset(data,idents = "T7")

data_FD4_LC2_silent<-rowMeans(data.frame(GetAssayData(FD4_LC2_data))[FD4_LC2_silent,])
mean(data_FD4_LC2_silent[!is.na(data_FD4_LC2_silent)])



FD4_LC2_silent_no<-c("ABCB4","KCNG1","SGPL1","CPSF3","ENPP1","IRGC","HTT","ADNP2","SERPINA4","GPR162","BSN","TRAK1","ILKAP","KIF3C","QPCT","OR14C36","KIAA0040","LDB2","FBN2","LPAR4","BCAT1","ODF2","BEND4","MOGAT3","MYO7A","IGFALS","UTS2R","RRBP1","KRT1","REPIN1","IRAK2","ANKRD52","CWF19L2","IGFN1","FAT1","MID1IP1","ATP6V1B1","SLC8A1","KL","HJV","ZNF835","KERA","SCD5","TECRL","XKR5","BEST3","EXOSC4","C22orf42","MKI67","CDA","AMPH","NRXN2","MSGN1","SLC25A18","HCCS","PLPPR3","GABRB1","PGLS","CHRNA3","LCE1F","SLIT3","CAD","THSD7B","CPB1","ATG4C","CLSTN2","CATSPER2","ARAP2","ADAMTS16","DTX3L","TMEM236","HAS1","ICOS","FBXO21","LRP2","POM121L2","TCF23","CDR1","TMEM132D","PLEKHM3","OTOGL","HELZ2","CCDC88C","C3AR1","ABHD17A","KLHL24","STUM","KRT3","CHST3","TTC37","SORL1","NAV3","ZNF735","SLC6A15","DSCAM","SCN10A","SVEP1","TANGO6","ZNF615","HOXB3","IK","GASK1B","KIF1A","ASH1L","ZNF729","INHBC","CNPY3-GNMT","PATL1","PSD2","PCDHB2","TM6SF2","RELN","MAPK4","REC114","KCNH8","SDK1","ZNF451","KIAA0100","TRIM49","RFX6","DOCK10","CDK18","PDZRN3","NPHP4","SPTB","CCDC187","OR1C1","SIGLEC7","PTPRN","PARP14","NPIPB12","PPIP5K2","ELMO1","PSG1","DOCK2","TTC30B","CDH12","IQGAP3","CENPE","PRR14L","SALL2","KMT5A","PIKFYVE","AGRN","CRYBA4","CALHM3","DNAH11","ESPNL","XPC","C2CD4C","PYGB","IFRD2","KIF4B","ZIC1","CDH10","HPSE2","SFRP1","MFNG","CORO2A","VCAN","SLC24A3","ENOX1","GZMA","CLGN","RPGRIP1L","EFCAB1","PABPC1","DMD","SYTL4","USP6","MAP4K2","GPA33","PCLO","BOLL","WDFY4","ZMYM2","NLRP14","MAP7D2","ADGRL3","CRTC2","SP8","VSIG8","KRTAP13-4","GNAS","SLC1A3","ZNF662","DENR","ACSL6","SEMA5A","SCN11A","ENTPD6","SNX25","UNC13A","CFAP46","SYNE2","BRINP3","URB1","FRAS1","ANKHD1","CDH18","DIAPH3","SUN3","MYO1G","MYH13","TCF3","ZAN","FERMT2","CDC45","SETX","NPIPB11","IL36G","PPP1R3G","MMP2","SPTA1","ZNF80","AHNAK","RYR3","EFHB","AGMO","STRN","PACSIN2","ALMS1","NUP210L","KDM4C","SPATA31C2","NALCN","QSER1","IFNA4","RFX4","TLR6","CD1E","PDHA2","LIG1","SRRM4","ANKRD36C","MEGF8","ZNF804A","CHRNA9","OR4A16","RRP15","E2F1","C6","BRPF3","GPSM1","SCLY","GLRA1","PCDHGB2","RAG2","CNTNAP4","PLIN4","ACTN2","PDGFRB","KCNA5","SGCD","TMEM14B","EHHADH","TPCN2","KCNA1","CFAP65","EBF2","PLXNA1","TRIM58","PWWP3A","DYNC2H1","PLBD1","DKK1","FA2H","DTHD1","STXBP5L","MUC12","EIF5B","SYT6","TRPM3","POTEH","CADPS","ZNF616","PKHD1L1","PPFIA2","DEK","NYNRIN","CEP295NL","KALRN","EGR4","ANKRD30A","PCDHA10","ELMOD1","XIRP2","BEND5","RBFOX3","SPDYE21","NTF3","ZNF853","ESRRB","ZNF765","ADRA2A","ITGA11","CACNA1D","ZNF415","F5","RBM14-RBM4","ADAMTS6","TAF1L","WIZ","DHRS9","MICAL2","OR14I1","TEX13C","KCNU1","TMEM86A","CSMD1","ZFYVE1","PRDM4","PITPNB","CRACD","SPATS2","DDX11","LAMA4","SPON1","CIAO1","ASB3","ARHGAP5","B3GNT4","OR2T12","EPB41L3","PKHD1","FOXD4L3","PLEKHG1","CALHM5","SAMD4A","KRT32","PRAC2","FIGNL1","PDE2A","ADGRG4","ZFHX4","DIDO1","C10orf67","C4orf50","APC","EPPK1","SERPINB3","SATB1","OR2L13","TH","GAB4","NCAPG","ANK2","PPP1R26","KCNJ4","PRAMEF12","LTBP4","OTOF","PCDH10","C4orf54","ADAMTS7","LAMC3","PPP1R16B","VHLL","IRF2BPL","SCN3A","RIT2","DLK2","INTS14","LYAR","CDH19","EPG5","PRPF40A","MARCHF4","HNF1A","KCNJ12","CCER1","UNC93A","JPH3","FAM205A","GRM7","GPAM","MUC5B","GRM5","TENM3","COL6A3","GSTA5","SDK2","TRIP12","SLC26A8","CFAP47","AGO1","GRB14","ZP1","BRSK2","CTNND2","P2RX6","IL20RA","NUP35","STARD8","PCDHA9","GABRD","SYTL2","CLDN20","OR6C68","ADGRF2","ZRANB3","DCTN2","GOLGA6C","ATP6AP1L","CCDC148","SLC22A6","COL21A1","SH3BP4","OR2G6","PRSS23","AJM1","FAM227A","MCMBP","CCDC40","SIRPA","KRTAP7-1","OR6F1","ATP11A","UBOX5","ZFHX3","L3MBTL3","KRAS","KAZN","MS4A8","RBM10","COL24A1","ASB7","PERM1","NID2","SAFB2","COL18A1","SLC22A7","HMCN2","SHANK1","EBNA1BP2","COL6A5","TRIM49B","UBAP2L","TCERG1L","OLFML2B","SYN2","MUC16","EPHA3","REEP4","XPO5","IL36A","POSTN","NBPF9","SLC12A6","CASR","TTN","RFPL4AL1","PTF1A","DENND1A","TICRR","ERCC6","SLC9A4","CROCC2","TSPOAP1","ANKRD12","GPR101","ANKK1","STRA8","FOLH1","KRT36","PLCG2","OR9K2","PCDHB9","SLITRK3","NPR1","NFATC1","ZEB2","TRPC4","MYCT1","PTPRQ","PXDNL","RNF146","PERP","C17orf80","TAF3","GLRB","TENM2","USP19","TCF7L2","OR56A5","RABGAP1L","TPR","ADAMTS8","ZNF667","CCDC60","ZDHHC8","FAT4","CYP11B1","CD6","IRF2BP2","KCTD1","PATJ","ZNF284","ZNF469","LPAR3","KCNA4","INPP5E","KCNB2","PIK3C2B","SIGLEC1","ACACB","LRP1B","KLHL32","SAP25","KCTD16","DRC7","TUBB4B","TRAPPC8","DAW1","OR2W5","HTR1F","ZNF365","PCDH11X","ZFAND1","HSPA8","PCDHB12","AOC1","GNRHR","PCDHA1","FOXD1","MYO5B","TBC1D16","CEACAM20","CCDC170","ANKMY1","LYSMD4","VWDE")
FD4_LC2_data<-subset(data,idents = "T7")

data_FD4_LC2_silent<-rowMeans(data.frame(GetAssayData(FD4_LC2_data))[FD4_LC2_silent_no,])
mean(data_FD4_LC2_silent[!is.na(data_FD4_LC2_silent)])


#FD5

FD5_LC1_silent=c("FRG2C","TMEM151A","DENND3","OR4Q2","POU3F2","MUC16","NLRP5","CYP2W1","ZAN","AMIGO1","EEF1A1","EXOSC6","PCDHA2","ELOVL7","OTOP3","CTNNBL1")
FD5_LC1_data<-subset(data,idents = "T8")

data_FD5_LC1_silent<-rowMeans(data.frame(GetAssayData(FD5_LC1_data))[FD5_LC1_silent,])
mean(data_FD5_LC1_silent[!is.na(data_FD5_LC1_silent)])


FD5_LC1_silent_no=c("KMT2A","TMEM150B","ANKRD63","PLEKHA7","CHD4","TXNDC5","CLOCK","ANKRD36C","PCDHB7","MUC16","CCER1","ASH1L","DDX11","BRD9","TPRN","DEK","NDUFAF6","ESPN","GRIK5","VWCE","ASH2L","CERCAM","GPR180","ZAN","USP32","PARVG","LCE1E","RPS6KA2","CYP1A2","SRR","CHD1","SMPD1","SOWAHA","PLXNB3","ARAP1","FOXA2","EGFR","OBSCN","FOXE3","PRAMEF15","POLR3H","RIOK1","DST","KDM4A","ATAD2","PABPC1","SORBS1","RBM10")
FD5_LC1_data<-subset(data,idents = "T8")

data_FD5_LC1_silent<-rowMeans(data.frame(GetAssayData(FD5_LC1_data))[FD5_LC1_silent_no,])
mean(data_FD5_LC1_silent[!is.na(data_FD5_LC1_silent)])





FD5_LC2_silent=c("TUBB8B","HPD","AKR1B10","PKNOX1","FRG2C","TMEM200C","OR5AS1","KMT2D","FAM81A","SCAF1","DLGAP4","C20orf96","EP400","AKNA","FAT1","DTNB","NOP14","EXOSC6","EPHA2")
FD5_LC2_data<-subset(data,idents = "T9")

data_FD5_LC2_silent<-rowMeans(data.frame(GetAssayData(FD5_LC2_data))[FD5_LC2_silent,])
mean(data_FD5_LC2_silent[!is.na(data_FD5_LC2_silent)])



FD5_LC2_silent_no<-c("ITPR3","UNC80","RAET1L","HAUS7","CACNA1A","TYW1","GK2","NR4A3","BMP2K","SPDYE21","ARAP3","CRYBG2","RRP7A","CDH23","NLRP6","PDE4C","PIK3C3","PRR18","MAML3","RAVER2","SAMD11","TEX15","ENOX1","EGFR","ARNTL","CADM4","NANOG","KBTBD13","NDUFB7","CHD3","NEB","DCUN1D4","LIMK1","RIOK1","AGBL5","FLACC1")
FD5_LC2_data<-subset(data,idents = "T9")

data_FD5_LC2_silent<-rowMeans(data.frame(GetAssayData(FD5_LC2_data))[FD5_LC2_silent_no,])
mean(data_FD5_LC2_silent[!is.na(data_FD5_LC2_silent)])





#FD8

FD8_LC1_silent=c("GPD1L","CNIH2","MAST1","MYBBP1A","H3C1","SRGAP2","ASH1L","CTNNA2","OR4K17","P3H1","CDIPT","RUFY4","MAB21L2","AP2A1","OR10P1","SUGP2","DYRK2","ACIN1","BICD1","NANOGP8","PLPPR2","ALKBH6","FAT1","GTF2IRD1","SLC12A2","ST8SIA5","KCNU1","CEACAM1","FANCI","HSPA6","UROC1","TMEM63B","ABHD11","OR5F1","SIN3A","PDE1B","TRHR","SYNE1","ALPK2","ATP1A3","IDI2","SMG7","NELL2","RHOBTB1","NDUFA7","SYMPK","B3GNT4","PACSIN1","CRNKL1","RAD54L","RLF","ZNF285","OR52B6","AVIL","NR2F6","NCOA6","DENND4A","HSPG2","PCDHGA3","ADGRB1","ZNF845","C6orf118","NLRP12","PRR35","TXNDC16","PCDHB15","CA4","TEX45","RASSF5")
FD8_LC1_data<-subset(data,idents = "T12")

data_FD8_LC1_silent<-rowMeans(data.frame(GetAssayData(FD8_LC1_data))[FD8_LC1_silent,])
mean(data_FD8_LC1_silent[!is.na(data_FD8_LC1_silent)])


FD8_LC1_silent_no=c("PREX2","SEL1L2","HSPA5","HRNR","RPS6KL1","TMCO4","CPNE8","CCDC63","INTS1","PRICKLE1","PLEC","SLC26A5","GABRB2","CLASP1","PYGO1","CCNYL1","GSTK1","TTC13","MYH7B","CRYBG3","AKAP13","KCNK13","ABI3BP","FAM120C","GTPBP1","PIGR","MAEL","SLC33A1","SYP","ABI1","DPYS","WDR97","TAF6L","TPX2","ATP8B1","VWF","DPP9","ASCL1","SLC8A2","PCDHGC3","SEPTIN4","TLL2","POLQ","MAP3K13","RYR3","EHBP1L1","PLA2G15","ANP32E","SNX17","CIART","PKHD1L1","AGO2","ANKRD50","IGFBP7","VPS45","BTNL8","PATZ1","EMSY","SLC35A2","POMT2","SYNPO2","TTN","ZNHIT1","ESPL1","CRNKL1","UCN3","SNX27","CDH16","CEBPA","CENPE","TMEM131","MRTFB","HORMAD1","PDE3B","CAP1","CDH10","OSR2","CSMD3","RTP1","POFUT2","SRD5A1","TNS1","CCDC105","HGFAC","SIRT1","USH2A","MRGPRX1","ADGRB2","EGFR","NCAM1","NRBP2","KIFAP3","PCDHB12","WDR43","ITIH6","ZSWIM8","UPRT","FAM193A","KAZN","DGKH","GOLGA6L4","RAI2","SLC25A31","GAREM2","LRRC37B","MRPS31","FBRS","ABCC3","MYO9B","GCKR","WDR87","RPGRIP1","ZNF284","PPRC1","RASSF8","DNAH9","PCYT1A","CTLA4","ELFN1","ZNF672","BRD9","CCDC17","KIF4B","ADAR","A2ML1","DSG3","GTF3C1","SRSF1","TUBB2B","PRRC2B","BOP1","DEXI","IKBIP","IGFL3","VPS18","ZNF239","WDR90","ZNF446","HEATR3","ZFYVE16","SOX17","SLC15A1","LY9","PIK3CA","SMG7","ZNF71","AGAP9","TRIM49","USE1","ERBB2","KIF14","CLMP","TLX3","SHH","ZNF605","ATF7IP","LRP8","NLRP5","TEX15","NT5E")
FD8_LC1_data<-subset(data,idents = "T12")

data_FD8_LC1_silent<-rowMeans(data.frame(GetAssayData(FD8_LC1_data))[FD8_LC1_silent_no,])
mean(data_FD8_LC1_silent[!is.na(data_FD8_LC1_silent)])





FD8_LC2_silent=c("SLCO4A1","PRKD3","PCDHGB5","FUS","PDHA1","DENND4B","NTN5","CFAP53","SRGAP2","YIF1B","GOLGA6L7","DNAH5","FNDC1","SRRM2","HGFAC","CLEC17A","DUSP2","TNRC18","PSMC3","FOLH1","RANBP2","KLF17","NANOGP8")
FD8_LC2_data<-subset(data,idents = "T13")

data_FD8_LC2_silent<-rowMeans(data.frame(GetAssayData(FD8_LC2_data))[FD8_LC2_silent,])
mean(data_FD8_LC2_silent[!is.na(data_FD8_LC2_silent)])



FD8_LC2_silent_no<-c("MUC16","DOCK1","NBR1","TBP","PRKD2","ARHGAP6","KCNT2","ARID1B","HOOK3","C16orf92","SORL1","GLMP","F2RL1","EPPK1","KIAA0754","RYR1","EOMES","CLSTN2","TCF12","NET1","KLHDC3","MUC3A","OR8G5","CD276","MBD3","HLTF","TRIM49","CD200R1L","SORBS1","TDG","AHNAK2","PCDHA7","RECQL4","MPND","EGFR","CRHR2","GABRG1","BTNL3","SFTPB","ANKRD36","SEMA6A")
FD8_LC2_data<-subset(data,idents = "T13")

data_FD8_LC2_silent<-rowMeans(data.frame(GetAssayData(FD8_LC2_data))[FD8_LC2_silent_no,])
mean(data_FD8_LC2_silent[!is.na(data_FD8_LC2_silent)])






#FD9


FD9_LC1_silent=c("CAPNS1","ARID3A","PHLDB1","SCN8A","TIAM1","GET1","ZNF316","MYOF","OR5H1","LRFN3","IFNL2","FEZ2","NAV1","XIRP2","DSE","MYT1L","OR1J4","OR51M1")
FD9_LC1_data<-subset(data,idents = "T10")

data_FD9_LC1_silent<-rowMeans(data.frame(GetAssayData(FD9_LC1_data))[FD9_LC1_silent,])
mean(data_FD9_LC1_silent[!is.na(data_FD9_LC1_silent)])


FD9_LC1_silent_no=c("FBXO30","PDE6G","SIMC1","LRRC66","LOC101928120","TNRC18","RPS6KC1","ADAMTS18","ANKRD24","LRIF1","OR5H1","KIF18B","ZNF875","PPP1R13L","MCM7","CREB3L1","FGFBP3","TEX45","POU5F1B","NOL8","ZFP92","WDR62","FLG2","SLC44A5","KIF19","PTP4A2","TP53","EPHA7","SMARCAD1","TNIP1","IL27","FOXN4")
FD9_LC1_data<-subset(data,idents = "T10")

data_FD9_LC1_silent<-rowMeans(data.frame(GetAssayData(FD9_LC1_data))[FD9_LC1_silent_no,])
mean(data_FD9_LC1_silent[!is.na(data_FD9_LC1_silent)])





FD9_LC2_silent=c("MMP2","KMT5A","NPTXR","OR5H1","FOXJ2","ASTN2")
FD9_LC2_data<-subset(data,idents = "T11")

data_FD9_LC2_silent<-rowMeans(data.frame(GetAssayData(FD9_LC2_data))[FD9_LC2_silent,])
mean(data_FD9_LC2_silent[!is.na(data_FD9_LC2_silent)])



FD9_LC2_silent_no<-c("PSMB5","HSPG2","ERBB2","VWA8","EPB41L5","OSBPL3","OR5H1","LZTS2","NTN5","PCDHGA8","ANKRD36C","CD160")
FD9_LC2_data<-subset(data,idents = "T11")

data_FD9_LC2_silent<-rowMeans(data.frame(GetAssayData(FD9_LC2_data))[FD9_LC2_silent_no,])
mean(data_FD9_LC2_silent[!is.na(data_FD9_LC2_silent)])





#FD14

FD14_LC1_silent=c("RBM6","HADHA","DNMBP","TDG","MAMDC4","ANKS6","ELK1","CSF3R","ARSD","C6","CARTPT","KCNG3","HOXD13","FAM53C","GSK3B","POFUT1","KMT2D","GSN","ZDHHC11","SLC7A1","PAX5","ST3GAL3","MON1A","ADAM11","TMEM14B","CPXM2","MED30","ARID1B","OR2L2","CBLN1","OGFOD2","GRIK3","C20orf173","GRIN2B","FOXF2","ILF3","PGLYRP4","IRS2","GABRR2","MICAL1","PLEKHA4","GEMIN5","SHANK1","KLHL32","PERM1","ABCA7","STYXL2","RABGGTB","ZNF587","NUFIP1","TEX9","RBM23","FOXO3B","SP8","HTT","RYR3","TRMT5","SLC24A4","SLC38A10","ANKRD27","SPPL2B","CCDC62","FAM221B","FOXO4","F2R","WDR25","ARID5A","MICAL3","LMAN1L","XAF1","TFAP2C","FOXC1","BAHCC1","SALL1","CHIC1","PLXNC1","SLC12A4","ADAMTS6","MRPL51","GTF3C3","PC","NETO2","ZNF534","SLC7A11","YBX1","PRSS33","FDPS","PLPPR3","PDZD7","TEKT5","SYTL1","BRINP1","DNASE1L2","JAK3","MTARC1")
FD14_LC1_data<-subset(data,idents = "T14")

data_FD14_LC1_silent<-rowMeans(data.frame(GetAssayData(FD14_LC1_data))[FD14_LC1_silent,])
mean(data_FD14_LC1_silent[!is.na(data_FD14_LC1_silent)])


FD14_LC1_silent_no=c("LRIT3","IKBKE","MFSD14B","LTBP2","SENP3","OR51L1","ZNF75D","GSTM2","TMEM14B","CFAP65","ZNF646","HOXB4","COBL","ABCA4","DZIP1L","SNX17","TDG","AKAP8L","TNPO2","NCOA4","KIF3C","SLIT1","IFNA2","STX2","ARFGAP1","RDH8","GALNT14","SMG9","TOGARAM1","DYNC2I1","PRR12","MYO1G","EGFR","MFAP3L","AAAS","UTRN","BCL11A","ACOXL","IL2RA","HYAL4","CECR2","LAMA1","SOCS2","FOXJ1","PEG3","PDGFRL","VAV2","MYH15","ANO8","CDH11","NANOS1","KLHDC7A","PROSER3","TTN","MLX","TLL1","FMO4","YARS1","ECM1","C2CD3","NCF4","MRGPRX3","MMP25","CD53","PTPRB","ZNF586","POGK","MAGED1","PLEKHG2","KRTCAP3","ZFHX4","POLG","NOTCH2","LUC7L3","LOXHD1","ANKRD17","VWF","DOCK6","ARHGEF17","FOXE1","PTPRO","ILF3","RYR2","SHISA8","EEF1E1","PARD6G","HNRNPU","TRAF7","OR51A2","FMN2","NDRG4","IGDCC4","SMYD3","CACNA1G","CENPV","IRX3","AURKB","BCAS3","SLC2A3","CELSR3","ESAM","NENF","NPIPB11","ZMIZ2","PCDHAC1","FASN","CHD1","POLD1","ALOX15B","HTRA3","ANO2","SPTA1","ZDHHC11","GUCY2C","RARA","ANAPC1","SPTLC3","ZBTB47","UNC5C","SCTR","MAPK8IP2","WDFY4","VPS39","PLD4","TMEM199","OR13C5","RTL6","SIGIRR","POM121L12","XIRP1","TICAM1","ANKRD11","LARP6","GAS8","SRP68","AMDHD2","GATM","SEZ6L","ITPR3","TTLL4","C11orf95","SLC5A4","SNX29","VSTM2L","ATXN7","MS4A10","SOS1","NBPF12","APC2","HYAL2","SLC34A2","ZNF784","KCNQ4","CHI3L1","PLCD3","LNP1","PIK3R2","PIAS3","FAAP20","PTPRG","GPC2","BEST2","KDM7A","CCDC130","BIRC6","PRAMEF17","GRIN2D","QRICH1","NEDD9","C1QTNF1","SLC6A5","NFAT5","PRUNE2","RAP1GAP","NEUROD2","SGCG","LRRC71","ENTPD7","AUTS2","COL27A1","GPKOW","CKB","MIGA2","SAFB2","PTPRM","TEKT5","MCF2L","JPH2","PKN1","SIMC1","NCOR2","IFNW1","PER1","SLX4","AHNAK2","SP9","MDK","CHPF","TYW1","PPFIA2","MUC16","UNC93A","CYRIA","BRF1","KIRREL1","FGFR4","GRK2","EXT1","CASQ1","BICRA","DHX8","ABHD17C","CASKIN2","DNMT3B","OFD1","BTBD3","C17orf75","IQCA1L","SERGEF","PPP1R3F","LCT")
FD14_LC1_data<-subset(data,idents = "T14")

data_FD14_LC1_silent<-rowMeans(data.frame(GetAssayData(FD14_LC1_data))[FD14_LC1_silent_no,])
mean(data_FD14_LC1_silent[!is.na(data_FD14_LC1_silent)])





FD14_LC2_silent=c("ASS1","PARP12","GPA33","CLEC18B","UBAP2","CCDC102A","SPTB","ALPK3","SERPINC1","INCA1","IRGC","OR2D2","TRPM3","PLEKHG4","PLCE1","GTF2IRD2","SLFN14","NUP88","RNF185","MYPN","ANKRD9","EPHA2","ZNF407","TMEM14B","PLEKHM2","ACLY","IRF2BPL","RGS7","PLCG2","ARAP1","USP36","C3orf18","ZNF587","UBASH3A","CHD5","GHRL","UBB","CDH2","TNFAIP8L1","CDCP1","MBD6","APOL3","CCDC62","LRP8","UPF3A","KANSL3","MTHFR","ZNF316","MUC17","PHKB","WNK1","ZFP92","ARHGEF11")
FD14_LC2_data<-subset(data,idents = "T15")

data_FD14_LC2_silent<-rowMeans(data.frame(GetAssayData(FD14_LC2_data))[FD14_LC2_silent,])
mean(data_FD14_LC2_silent[!is.na(data_FD14_LC2_silent)])



FD14_LC2_silent_no<-c("GSE1","DEAF1","SDK1","PIK3CG","FBLN7","LLPH","PTPRJ","MROH1","FXR2","PTGER1","SLC5A12","TDG","THAP4","C2CD4C","CPAMD8","FAM126B","GIGYF2","SEC24B","FGFRL1","DOCK3","SLC23A3","SLC6A8","DCAF1","U2AF2","CBARP","ZP1","UBR5","JMJD6","APC2","QRICH2","CAMK2A","THG1L","STX11","BAZ2A","MSS51","ARC","STARD7","CASP10","ZNF773","IKZF4","COLGALT2","TEKT3","POM121","QPCTL","FSCB","ISLR2","CHD3","AFF4","PRAMEF17","TP63","ANAPC1","NR1H4","FBXO11","CDAN1","SH2D2A","AKR7A2","SLC2A3","KIAA0100","POLR2M","TTC12","XIRP2","DDRGK1","MYEOV","SPATA31D1","AUTS2","ARID5B","LCT","AHNAK","TANC2","PEBP1","MUC5B","LHCGR","EPB41L2","MYH14","NCK2","SLC35G4","URB2","CDC25A","DENND3","ZFYVE1","CDC42BPB","LOXL4","OR13C5","PRICKLE1","HTRA1","PHLDB1","MAGEA6","ASCL1","ADAM33","TRIL","TXNDC5","NBPF11","RNF150","MEX3C","GFER","CASTOR1","ZFPM1","CELSR1","PCM1","CACHD1","SLFN5","DAAM2","VAV1","RABEP1","CMA1","VCPKMT","FAM234A","CACNA2D2","AP3M1","TMEM14B","NHLRC3","IFITM1","SLC26A1","EHF","CLEC18B","FMO1","CD14","NBPF9","SYNPO2","YTHDF3","TMEM163","MOB1A","CDH15","TRPM2","COL7A1","MTFMT","DNAH11","RNF227","CRHR2","RARA","MEP1A","OBSCN","DMBX1","RBMXL3","DDX31","SEMA5A","MIEF1","EIF4A3","TEX14","GALNT16")
FD14_LC2_data<-subset(data,idents = "T15")

data_FD14_LC2_silent<-rowMeans(data.frame(GetAssayData(FD14_LC2_data))[FD14_LC2_silent_no,])
mean(data_FD14_LC2_silent[!is.na(data_FD14_LC2_silent)])




#FD16

FD16_LC1_silent=c("RTP2","STUB1","ONECUT3","POLR2I","KLF4","EIF2S3B","PLD5","STK32C","CXCR2","MYBPH","ABCD1","NPAS4","SHROOM4","RHOBTB2","EPB41L3","RNF213","PUM1","PORCN","VWA5B2","MYH7","DDX23","SETDB1","DNPEP","ZMYM4","SLC22A23","ANK1","KCNJ3","PRUNE2","DPP3","NHSL2","SHANK1","KPNB1","TAOK3","TNKS1BP1","COL1A1","PRDM9","TRIM2","FAM171A2","KCNQ1","FAM135B","SMARCD3","OPN4","TRIL","CIT","SUPT16H","MEGF8","SLC24A2","DNMT3A","EML3","SLC25A42","UAP1","NTN1","XRCC6","ABCA2","HS6ST3","CHD8","TMEM245","E2F3","BMP5","UBC","NEB","HMCN2","CCR10","CEP19","PHKG2","NCKIPSD","ANO4","ASAP2","SRD5A2","MTUS2","RFK","ACTR1A","GLI4","SEC23B","THBS1","RANBP10","MED12","NOTCH3","PCDH19","POU6F1","ZNF865","AACS","ALDH1L2","UNC93B1","TAMALIN","KDM3A","SPTAN1","KDELR3","CIAPIN1","DCAF7","AHNAK2","PCDH20","LOXL2","LRFN2")
FD16_LC1_data<-subset(data,idents = "T16")

data_FD16_LC1_silent<-rowMeans(data.frame(GetAssayData(FD16_LC1_data))[FD16_LC1_silent,])
mean(data_FD16_LC1_silent[!is.na(data_FD16_LC1_silent)])


FD16_LC1_silent_no=c("TOGARAM2","ATP13A1","RBMS3","CALR","MYH9","FOS","LDLRAD3","KDM3B","AHNAK","TBC1D30","CXCR2","SLC5A3","GLI4","ADGRB3","PABPC3","MYH7","EGFR","RNF113B","DTL","KRT10","PRR16","NUDT16L1","IQSEC2","KXD1","DNPEP","CTNNB1","HNRNPA2B1","SRSF9","PCDH20","SUN2","H1-4","KCNN4","COL1A1","GBP4","IFI30","POTEJ","DDX5","ROGDI","NCL","SLC12A3","AAK1","VAT1","ARMC9","VWA5B2","ZNF292","SPTAN1","CD3E","SSR2","SRD5A2","BRD4","ACP2","BTBD17","CT45A2","NKG7","ZFP36","PCDH19")
FD16_LC1_data<-subset(data,idents = "T16")

data_FD16_LC1_silent<-rowMeans(data.frame(GetAssayData(FD16_LC1_data))[FD16_LC1_silent_no,])
mean(data_FD16_LC1_silent[!is.na(data_FD16_LC1_silent)])





FD16_LC2_silent=c("SOX9","SPTBN4","ATP1A2","ESRRB","HYAL2","ARPC1B","TEAD4","LLGL2","MGAT3","PLPP3","ZNF703","PCDHGA10","RAP1GAP","FLAD1","KCNN2","MASP2","GRIPAP1","GSX1","TRIOBP","GAL3ST4","MKI67","SLC2A1","FOXL2","PCOLCE2","GALNT10","ADGRB2","RAP2A","TFE3","ZFHX3","PIK3R5","SON","HID1","MYO7A","NRXN2","MYCL","PPP1R3F","UBC","CTDNEP1","PDZD8","MYRFL","AHNAK2","HIC1","MEN1","COL4A5","ATG16L2","HDAC11","DIP2C","TIMM29","GRIA2","NOTCH2","ZNF574","KCNH1","FBN3","KIF1C","BAP1","ANKRD52","EIF3B","VIT","KCNA3","KANK4","ABCA1","NELFB","ACTN2","PDE7B","CELSR3","EVL","DLGAP3","CACNA1H","FBXO39","HERC2","MEX3B","CAMTA1","CLDN9","IQSEC2","ITPRID2","RHOV","SMARCD3","MAN2C1","CLEC16A","KCNS1","CAD","OXTR","DCLK1","ELAVL4","ZNF281","CD276","PPP2R5E","BTBD11","NRTN","MEIS1","HK2","NUMA1","CCDC171","GGT7","KLF2","CNNM2","WNT4","GARNL3","PIK3R2","FGD1","SGSM2","COL4A6")
FD16_LC2_data<-subset(data,idents = "T17")

data_FD16_LC2_silent<-rowMeans(data.frame(GetAssayData(FD16_LC2_data))[FD16_LC2_silent,])
mean(data_FD16_LC2_silent[!is.na(data_FD16_LC2_silent)])



FD16_LC2_silent_no<-c("CT45A1","TMEM202","E2F3","TATDN1","OXTR","TIMM29","AIPL1","MYL6B","VPS13C","RPRD2","KRT6B","CALR","NBEAL1","SCUBE3","FLNA","DDX5","ISM2","PCDHGA10","MAN2C1","TNRC18","BCL9L","EML5","FCGR2A","HYAL2","BCL2L12","FN1","CLEC4M","COL4A6","PRDM9","UPRT","CLEC16A","EGFR","MYLIP","MYCL","NBPF11","SLK","WNT2B","MKI67","AHNAK2","KCTD1")
FD16_LC2_data<-subset(data,idents = "T17")

data_FD16_LC2_silent<-rowMeans(data.frame(GetAssayData(FD16_LC2_data))[FD16_LC2_silent_no,])
mean(data_FD16_LC2_silent[!is.na(data_FD16_LC2_silent)])





#另一样本的表达量

FD1_LC2_silent=c("CEP131","OR11H12","SHISA2","CHST5","IDH3A","RNF213","CROCC")
FD1_LC2_data<-subset(data,idents = "T4")

data_FD1_LC2_silent<-rowMeans(data.frame(GetAssayData(FD1_LC3_data))[FD1_LC2_silent,])
mean(data_FD1_LC2_silent[!is.na(data_FD1_LC2_silent)])



FD1_LC2_silent_no=c("OR2L3","GGT2","AZIN2","MAPKBP1","KAT2A","ALPP","SLAMF9","DRICH1","ANKLE1","COLQ","MEX3D","ATAD3B","RBM47","MUC12","RPL21","CDC25B","IL17D","ZXDA","PDSS1","UNC93B1","SHROOM4","GNAQ","WIZ","SLC34A3","DUSP5","CCDC61","RIOK1")
FD1_LC2_data<-subset(data,idents = "T4")

data_FD1_LC2_silent<-rowMeans(data.frame(GetAssayData(FD1_LC3_data))[FD1_LC2_silent_no,])
mean(data_FD1_LC2_silent[!is.na(data_FD1_LC2_silent)])




FD1_LC3_silent=c("GPD2","DDX11","NBPF10","TENT4A","GSR","RASIP1","GON4L","PLEC","WASHC1","HOXC13","FRG2C","CASP10","OGFR","APC","SSC5D","CRY1","SQLE","FLVCR1")
FD1_LC3_data<-subset(data,idents = "T5")

data_FD1_LC3_silent<-rowMeans(data.frame(GetAssayData(FD1_LC2_data))[FD1_LC3_silent,])
mean(data_FD1_LC3_silent[!is.na(data_FD1_LC3_silent)])



FD1_LC3_silent_no=c("EPS15","SEC24C","CAPN14","NEUROD2","ACSF2","FBXO41","RBM10","COL15A1","NEFH","MAP7","BRPF1","MUC12","PNMT","C16orf90","ZBTB48","FRG2C","TRIOBP","MS4A14","CYP2C9","OR10K2","TNPO3","NLRX1","COL16A1","SEMA6A","PPIL2","TEKT3","ARMCX4","TRIM25","GPR62","ALPP","MID1","WIZ","AQP11","CACNA2D2","ZDHHC19","FLNA","ZFHX4","OR5F1","ZXDA","BMP6","B4GALT5","SYMPK","TRO","MME","RAB11FIP4","USHBP1","BAP1","RGR","CRY1","CDH7","JAML","KLHL40","LCT","EGFR","MAP3K5","ZNF846")
FD1_LC3_data<-subset(data,idents = "T5")

data_FD1_LC3_silent<-rowMeans(data.frame(GetAssayData(FD1_LC2_data))[FD1_LC3_silent_no,])
mean(data_FD1_LC3_silent[!is.na(data_FD1_LC3_silent)])












#FD2

FD2_LC1_silent=c("RTN1","PPM1H","SYT6","TRIM65","PCDHGA3","AR","LEPROTL1","MRGBP","ABCA1","ZNF845","C7orf26","ZPBP2","GRIN2A","BIRC5","TEX13C","KANK1","DDX11","ATAD2","SLC13A3","NXNL1","MRPL38","TBXAS1","RHOT1","LOC400499","KMT5A","BOP1","NBPF20","ASAH2")
FD2_LC1_data<-subset(data,idents = "T1")

data_FD2_LC1_silent<-rowMeans(data.frame(GetAssayData(FD2_LC2_data))[FD2_LC1_silent,])
mean(data_FD2_LC1_silent[!is.na(data_FD2_LC1_silent)])


FD2_LC1_silent_no=c("SVEP1","PAPSS2","ATP13A1","RENBP","WIPF3","ZNF8","ZNF257","LAMA1","ANKRD36","ARHGAP10","SHANK2","NKX2-3","MAMLD1","MAST3","POTEH","POU3F2","KDM3B","EPC2","TRPV1","PI4KA","CNTLN","CCSER2","DEK","ADAM30","ITCH","SULT1A2","NAV2","EHD1","SECTM1","REEP1","PTX4","GIN1","KMT5A","SBNO1","GNAS","VSIG4","MED12L","MEGF10","ZNF329","KCNH4","NAP1L1","ADCY8","GRIN2A","MDM2","BRPF1","NLRP9","FOXJ2","POLL","KIF13A","CELF2","NCOR2","LRRN1","CSAD","ANP32C","HCN1","VEGFC","ZDHHC8","NXNL1","FCRL3","AARD","MROH2A","PIEZO1","AGTPBP1","APBB1","FKTN","DSG1","MYBBP1A","USH2A","LMTK3","GFPT2","APC2","TDRD1","HSPA12A","NBR1","MAP3K15","CNGB1","OR8H2","HAPLN4","ACSL4","DACT1","MRPS34","DOLPP1","ARID3B")
FD2_LC1_data<-subset(data,idents = "T1")

data_FD2_LC1_silent<-rowMeans(data.frame(GetAssayData(FD2_LC2_data))[FD2_LC1_silent_no,])
mean(data_FD2_LC1_silent[!is.na(data_FD2_LC1_silent)])






FD2_LC2_silent=c("SLC12A2","PLIN4","ZNF845","SLC2A3","STAB1","FOXC1","KCNA6","PCDHB9","HSPG2","ALDH3B1","NUFIP1","MLLT3","SOX1")
FD2_LC2_data<-subset(data,idents = "T2")

data_FD2_LC2_silent<-rowMeans(data.frame(GetAssayData(FD2_LC1_data))[FD2_LC2_silent,])
mean(data_FD2_LC2_silent[!is.na(data_FD2_LC2_silent)])



FD2_LC2_silent_no=c("TJP3","FOXD4","GSC2","RSPH6A","LRRCC1","MTCL1","MYL7","ANKRD36","C1QL1","MTHFD2L","KCNMA1","NACAD","TENM2","ZC3H3","ZNF831","CSRNP1","KLF12","SYNE3","TRIM67","KMT5A","ITGA8","THSD7A","PRDM9","DEAF1","UTF1","SP8")
FD2_LC2_data<-subset(data,idents = "T2")

data_FD2_LC2_silent<-rowMeans(data.frame(GetAssayData(FD2_LC1_data))[FD2_LC2_silent_no,])
mean(data_FD2_LC2_silent[!is.na(data_FD2_LC2_silent)])



#FD4


FD4_LC1_silent=c("EGFR","SLC28A1","C1QTNF1","NWD2","ZNF587","FPR2","GJD2","SLC43A1","ADAM19","SPTA1","SMOX","TLX2","PKD1L3","RBPJL","APLNR","H3C10","C8B","ARAP3","ALDH3B1","INPP5A","XIRP1","NIPAL4","ERN1","NDUFS2","FUT7","D2HGDH","DST","AQP12B","APOB","GPR87","CDH18","ASXL3","PPFIA3","HMCN2","RETN","MAGEB5","RAX2","PHYHD1","CDRT1","KIRREL2","ZNF467","KRTAP24-1","KMT2C","RERG","SPATA31D1","ANK2","TUT4","TMEM132B","ATP8B4","RIPOR1","EML5","LRCOL1","SMC3","BPIFB1","SLC9B1","IL20","RUNDC3B","RANBP2","SIMC1","SHD","NUTM1","GRIA3","FREM3","NFIX","MIR1-1HG","ACTL8","ARVCF","COL7A1","CDH22","POTEH","NALCN","TRMT6","IL36B","KCND2","TSGA10IP","TOMM40L","FADS6","UTRN","BANK1","ISLR2","PAX2","OR5M9","PEAR1","RBM24","KCNH2","ITGA11","CYP4F3","CDX2","FAM186A","OLFM1","TBC1D22A","OBSCN","FOXD4L6","SLC4A2","CDK5RAP3","CNTN6","MAP3K1","CD248","MYO7B","SLC45A1","PCDHAC1","COL20A1","XAB2","RBM20","UBAP2L","SV2C","MED12","KSR1","TMEM63A","OR8H2","SH2D7","KLHL1","UNC80","ADAMTS4","EPHA10","KMT2B","TMEM14B","TGM5","AOC2","CACNB1","TRMT10C","HOXA1","SRRM2","CYBB","ARMC3","MYLK","RYR1","ELOA2","OR1A1","TANGO6","PMPCA","VCAN","TNRC18","SLC7A10","SLC6A7","ZC3H11B","ACAN","CTSV","HK3","SALL4","HMGB3","ADARB2","HTT","LCP2","C3orf20","EDA2R","CAMK2B","GNAQ","SHPK","SPATA22","CCDC191","EPS8L3","PLCB1","LOC100996413","API5","OR9I1","FREM2","DDX11","ADPGK","TCEAL3","C11orf95","CASR")
FD4_LC1_data<-subset(data,idents = "T6")

data_FD4_LC1_silent<-rowMeans(data.frame(GetAssayData(FD4_LC2_data))[FD4_LC1_silent,])
mean(data_FD4_LC1_silent[!is.na(data_FD4_LC1_silent)])


FD4_LC1_silent_no=c("GCOM1","EIF4G3","MTTP","KIF1A","ASCL4","FBN2","GRIN2A","RBM33","SCP2","LRFN5","ARHGEF15","ITSN1","NRP1","TPST2","PRUNE2","OR51S1","HSPB2","ALDH1A2","KLRG1","CX3CR1","TTC14","IGDCC3","DZIP1L","GIMAP6","HPD","TF","PCDH15","PKP2","CCNB3","MUC7","NTN1","ZNF599","MCM3AP","KRTAP22-1","MPL","DNAH5","ZNF512B","PCDHA7","SPATA31D1","IFT140","PCDH7","ZNF75A","KCNG2","SLFN11","TP53AIP1","PITPNM3","C10orf90","POU5F2","SPATA31C1","MAGEA11","PROX1","TCHHL1","CACNA2D3","RIOK1","PDPR","NCAM1","OR4D10","PCDHGB5","SHQ1","CACNA2D1","NPAP1","GFRA1","SLC26A1","PYGM","PLPPR4","ZDHHC11","FRMPD1","ANKRD17","IGSF9","NSUN2","CYLC2","IL25","MUC17","MAP1A","HS6ST1","CGB2","PCDHB3","BTNL3","TLL2","DRD3","NRK","GABRG2","LMAN2","ANKDD1A","PIK3CA","ZNF423","AR","DEK","DNER","CCDC8","CCDC57","ZNF76","NUMA1","ALDH4A1","CR1","SLC7A4","LNPEP","SIGLECL1","PCED1A","TDRD9","PSD2","CDH11","MATK","KIF27","CD160","APCS","KALRN","NINL","MCIDAS","TSHZ3","BUB1B-PAK6","NCKAP5L","CYFIP2","NLRP1","RPS12","MAGI1","CHRNA3","KRTAP4-7","GH2","AHSG","ZNF713","TBC1D4","SHC4","DNAH8","TNFAIP2","CDH4","MMP8","HES6","KCNB1","KRT6A","MN1","WDR6","NLRP3","MUC12","MBLAC2","COL6A6","FUT6","ANK1","ARID5B","GTF2E1","CHRM5","LRRIQ3","RERG","PMS2","CROCC","PAWR","ALMS1","MRC1","ADGRB2","KCP","FCRLA","SOAT2","ARID3A","CNOT8","FBXO30","CHD7","MICAL3","GATA5","EPHA7","UBE3A","MYCL","CDH20","MYH2","GSC","HIRA","GOLGA6B","ANO2","CDA","RSPH9","ELSPBP1","TMPRSS7","PCDHB9","ZNF415","CYP4F22","RIMBP2","MCAT","CNN1","KRAS","RBM10","NPPB","IL7R","ZNF485","DST","ZDHHC14","CAPN10","POTEA","CCDC171","ARFGEF2","TRIML2","KIF19","STYXL2","GFAP","SORCS1","CLIP1","COL20A1","TTC17","PLEKHH1","PRR23C","P3H2","ADCYAP1R1","PPP1R37","RBP2","FAM221B","ZFHX3","UBFD1","KCNK2","LGALS8","GMNC","PNPLA6","FLNB","PRDM8","LCE3D","NLRP8","DNAH6","SNX15","SLC28A3","ARHGEF9","MAP3K19","TSC22D1","RYR2","WFDC5","CYP2C18","SH3GLB1","OR51L1","FAT3","PRODH2","LRRK2","OR2C1","ARHGAP6","TSPAN16","PTGS1","GAREM2","KCNIP1","H6PD","SPTBN4","CHIA","RAF1","ALOX5AP","CEACAM6","PLA2R1","CCDC168","ADGRA1","PARG","PABPC3","CWF19L2","KRTAP15-1","KCNH5","SMAD9","OR5L1","C4orf54","DCLK1","FAM71B","FAM171A2","VAT1L","TRPM6","C16orf54","NUP210","THSD7B","TAS2R1","CYP11B1","USP29","LYN","UTF1","ARMC8","PCDHGB2","HEPH","KCNH6","KIAA1671","CELF6","KHDRBS1","ZNF148","SYT6","ITIH1","ZFHX4","UGT1A4","PGLYRP3","TMEM117","PIF1","DNTTIP2","PCSK7","CPXM2","WWP1","IMPG1","FAM166C","LCMT2","KRTAP10-12","ADGRG4","KCNH7","NAV3","TEPSIN","USP26","PUS1","PAX4","GCSAML","PCDHA10","MAP3K1","KIAA0040","OR52E4","ADAMTSL1","ZNF483","AHNAK","TBX6","ALX4","COL5A3","KCNAB1","FAM186A","PCDHB7","FOXN4","STAG1","SERPINB12","PARP8","SEMA3E","NBPF11","TNRC6B","ARHGEF4","LMLN2","MYLK4","HOXA1","ALPG","UNC5CL","GCN1","POT1","PELI2","RPF2","KCNN1","C2CD4A","NLRP12","FOXD4L3","PTPRN2","DACH1","ZC3H13","GJA8","CDX1","SLC17A3","GLDC","FERD3L","OR2T8","CAMSAP2","SEC31A","NAP1L3","OR10A7","MAEL","ADGRF4","PTPRB","PCDH11X","C3orf56","MPPED1","UPK1A","GDF2","CCBE1","DOCK9","OR10G7","INPP5A","FRMD5","ECEL1","SAGE1","TBC1D16","OR1A1","LOXHD1","FAM43B","FBN1","FMO3","MASP1","SETBP1","DHX33","ARHGEF38","GPR101","SYNE1","AEBP1","OR2W5","GPR78","SRC","EFCAB5","BTAF1","BTBD16","MAGEB2","POLQ","GPR87","TSGA10IP","TEX35","FMNL3","EYS","PCDHB6","TG","PTPN23","FXYD5","DCSTAMP","RIN3","PDILT","GAB4","PALB2","GLI3","KRT6B","DNAH3","ABCB1","VEZT","ZNF280D","WDR97","MYBPC3","CCSER1","DAB1","ZFPM2","CLRN2","CHD5","PABPC1","CUBN","AJAP1","DMP1","TNKS2","MYO18B","DNAL1","DCAF13","FEM1B","LRRC53","ITGA7","DGKB","COL3A1","ZEB2","CLTB","SLC12A1","OTOGL","TNRC18","FAM151A","CDH6","APC2","CD163L1","ATXN1L","EPPK1","OTOG","FNDC1","PKD2L1","KBTBD3","SYPL1","KCTD6","ITGA4","ESR1","QARS1","SPPL3","COL22A1","COPB2","PIGT","CACNA1I","GABRG3","RBM20","AHI1","SETD2","MMRN1","MAS1","DMRTB1","OR8K5","KRT76","WDFY4","OR5M1","TM9SF3","XPNPEP1","USH2A","SLC9A9","UNC80","IZUMO3","OLFML2A","COL1A2","PPRC1","EXOC3L2","CBLC","DBNL","NRXN1","GLI2","LOC100509620","COL4A1","PDE10A","MXRA5","OR4C6","MED23","ZNF226","SETD1B","CSMD2","KRTAP5-10","MROH9","ANKS3","NOLC1","ANKRD36","GRIP1","BSPRY","LINGO2","ABCA4","CKMT1B","EXPH5","CCDC40","MMP2","RYR3","TMEM14B","SLC52A3","EMILIN1","OR6K6","OR7E24","NAALADL2","RIMS1","DOCK2","SLC5A12","RAI1","NWD2","PRC1","ZNF675","RTKN2","IQCF5","CACUL1","CNBD1","ZNF114","MST1R","OR10A2","PCDHB4","SLFN13","KCNH1","ELMO1","F5","FREM1","P2RX7","CNNM1","CACNA1E","NBPF9","PCDH12","ODAM","LHX9","WIZ","SMCHD1","CHGB","SHISA9","ZNF514","THRA","MFSD4A","ZNF236","CBY2","ZSCAN32","PRKDC","MEGF8","NR2F1","SF3B4","OR2AK2")
FD4_LC1_data<-subset(data,idents = "T6")


FD4<-read.table(file="/Users/wangjun/Desktop/FD4.csv",sep=",",header = T)
FD4$FD4_LC1
FD4_LC1_silent_no=FD4$FD4_LC1
FD4_LC2_silent_no=FD4$FD4_LC2



data_FD4_LC1_silent<-rowMeans(data.frame(GetAssayData(FD4_LC2_data))[FD4_LC1_silent_no,])
mean(data_FD4_LC1_silent[!is.na(data_FD4_LC1_silent)])





FD4_LC2_silent=c("OR2T33","STRIP2","VWA5B2","ITPRID1","NKX6-1","OR52A1","OR8G2P","ADCY7","CDH23","NOTUM","GAPVD1","ACTR1B","PIAS2","OR5AR1","UBC","MYH1","KLHL3","YIPF7","HOXA9","AKAP6","GIMAP1","AJAP1","ISM1","LOC101928841","ANKRD52","STX7","TRPC7","FER1L6","RASAL3","RP1L1","SIGLEC8","DRD5","TAS1R3","HSPA12A","SMC3","WASHC2A","FAM98B","CLSTN2","ZNF606","CCDC89","YBX1","RRAS","ZNF711","PML","VSTM4","KMT5A","RANBP2","NCLN","DNAH5","RAB3IL1","NYAP1","TMEM255A","EPHB3","LRRC38","ADAMTS20","NMUR1","GFPT2","PNLIP","FMOD","CLK4","GNAQ","TCEAL3","NBR1","PIEZO2","OR56A3","LGALS12","CHMP1A","RGS9","RHOBTB1","ANKK1","COL3A1","ANKRD35","HNRNPA0","UNC45B","PYHIN1","OR8G1","WNT10B","PCDH18","PCSK5","CCDC60","CPS1","IPO9","MROH9","POM121L2","AFF2","TNFAIP2","MYH4","FGF13","HGF","PCDHA6","PCDH15","EXOC3","TM7SF3","STAR","NR4A2","SPEN","BTBD11","ZNF544","BHLHE22","MAGI2","P2RX2","TBX15","OBSCN","BTNL3","HOXB1","VPS33A","MGAT5B","CILP","KCNJ12","ADRB3","C1orf167","CD163L1","ALG13","DOK3","PITX1","NUFIP1","BAHCC1","SLC14A2","SLC5A10","POP4","L2HGDH","SKI","KIAA1109","MXRA5","KMT2C","ROBO3","CUX1","HOXD10","NSD1","SIRPD","GLT1D1","KLK11","AMBRA1","UMOD","CFAP46","KRT35","FCRL6","STAB2","SYT1","TTC9B","NSUN5","ASXL2","TRIM52","NALCN","FGB","UGT1A4","HRNR","IGFN1","GTSE1","DSG3","ZNF81","MAP1B","ATP6AP1L","DISP3","OR5AS1","PCNT","ITGBL1","IFI16","CAPN10","PPP4R3B","CYP11B1","POTEJ","TRPM4","SEMA3F","SLC12A5","STARD8","EXO1")
FD4_LC2_data<-subset(data,idents = "T7")

data_FD4_LC2_silent<-rowMeans(data.frame(GetAssayData(FD4_LC1_data))[FD4_LC2_silent,])
mean(data_FD4_LC2_silent[!is.na(data_FD4_LC2_silent)])



FD4_LC2_silent_no<-c("ABCB4","KCNG1","SGPL1","CPSF3","ENPP1","IRGC","HTT","ADNP2","SERPINA4","GPR162","BSN","TRAK1","ILKAP","KIF3C","QPCT","OR14C36","KIAA0040","LDB2","FBN2","LPAR4","BCAT1","ODF2","BEND4","MOGAT3","MYO7A","IGFALS","UTS2R","RRBP1","KRT1","REPIN1","IRAK2","ANKRD52","CWF19L2","IGFN1","FAT1","MID1IP1","ATP6V1B1","SLC8A1","KL","HJV","ZNF835","KERA","SCD5","TECRL","XKR5","BEST3","EXOSC4","C22orf42","MKI67","CDA","AMPH","NRXN2","MSGN1","SLC25A18","HCCS","PLPPR3","GABRB1","PGLS","CHRNA3","LCE1F","SLIT3","CAD","THSD7B","CPB1","ATG4C","CLSTN2","CATSPER2","ARAP2","ADAMTS16","DTX3L","TMEM236","HAS1","ICOS","FBXO21","LRP2","POM121L2","TCF23","CDR1","TMEM132D","PLEKHM3","OTOGL","HELZ2","CCDC88C","C3AR1","ABHD17A","KLHL24","STUM","KRT3","CHST3","TTC37","SORL1","NAV3","ZNF735","SLC6A15","DSCAM","SCN10A","SVEP1","TANGO6","ZNF615","HOXB3","IK","GASK1B","KIF1A","ASH1L","ZNF729","INHBC","CNPY3-GNMT","PATL1","PSD2","PCDHB2","TM6SF2","RELN","MAPK4","REC114","KCNH8","SDK1","ZNF451","KIAA0100","TRIM49","RFX6","DOCK10","CDK18","PDZRN3","NPHP4","SPTB","CCDC187","OR1C1","SIGLEC7","PTPRN","PARP14","NPIPB12","PPIP5K2","ELMO1","PSG1","DOCK2","TTC30B","CDH12","IQGAP3","CENPE","PRR14L","SALL2","KMT5A","PIKFYVE","AGRN","CRYBA4","CALHM3","DNAH11","ESPNL","XPC","C2CD4C","PYGB","IFRD2","KIF4B","ZIC1","CDH10","HPSE2","SFRP1","MFNG","CORO2A","VCAN","SLC24A3","ENOX1","GZMA","CLGN","RPGRIP1L","EFCAB1","PABPC1","DMD","SYTL4","USP6","MAP4K2","GPA33","PCLO","BOLL","WDFY4","ZMYM2","NLRP14","MAP7D2","ADGRL3","CRTC2","SP8","VSIG8","KRTAP13-4","GNAS","SLC1A3","ZNF662","DENR","ACSL6","SEMA5A","SCN11A","ENTPD6","SNX25","UNC13A","CFAP46","SYNE2","BRINP3","URB1","FRAS1","ANKHD1","CDH18","DIAPH3","SUN3","MYO1G","MYH13","TCF3","ZAN","FERMT2","CDC45","SETX","NPIPB11","IL36G","PPP1R3G","MMP2","SPTA1","ZNF80","AHNAK","RYR3","EFHB","AGMO","STRN","PACSIN2","ALMS1","NUP210L","KDM4C","SPATA31C2","NALCN","QSER1","IFNA4","RFX4","TLR6","CD1E","PDHA2","LIG1","SRRM4","ANKRD36C","MEGF8","ZNF804A","CHRNA9","OR4A16","RRP15","E2F1","C6","BRPF3","GPSM1","SCLY","GLRA1","PCDHGB2","RAG2","CNTNAP4","PLIN4","ACTN2","PDGFRB","KCNA5","SGCD","TMEM14B","EHHADH","TPCN2","KCNA1","CFAP65","EBF2","PLXNA1","TRIM58","PWWP3A","DYNC2H1","PLBD1","DKK1","FA2H","DTHD1","STXBP5L","MUC12","EIF5B","SYT6","TRPM3","POTEH","CADPS","ZNF616","PKHD1L1","PPFIA2","DEK","NYNRIN","CEP295NL","KALRN","EGR4","ANKRD30A","PCDHA10","ELMOD1","XIRP2","BEND5","RBFOX3","SPDYE21","NTF3","ZNF853","ESRRB","ZNF765","ADRA2A","ITGA11","CACNA1D","ZNF415","F5","RBM14-RBM4","ADAMTS6","TAF1L","WIZ","DHRS9","MICAL2","OR14I1","TEX13C","KCNU1","TMEM86A","CSMD1","ZFYVE1","PRDM4","PITPNB","CRACD","SPATS2","DDX11","LAMA4","SPON1","CIAO1","ASB3","ARHGAP5","B3GNT4","OR2T12","EPB41L3","PKHD1","FOXD4L3","PLEKHG1","CALHM5","SAMD4A","KRT32","PRAC2","FIGNL1","PDE2A","ADGRG4","ZFHX4","DIDO1","C10orf67","C4orf50","APC","EPPK1","SERPINB3","SATB1","OR2L13","TH","GAB4","NCAPG","ANK2","PPP1R26","KCNJ4","PRAMEF12","LTBP4","OTOF","PCDH10","C4orf54","ADAMTS7","LAMC3","PPP1R16B","VHLL","IRF2BPL","SCN3A","RIT2","DLK2","INTS14","LYAR","CDH19","EPG5","PRPF40A","MARCHF4","HNF1A","KCNJ12","CCER1","UNC93A","JPH3","FAM205A","GRM7","GPAM","MUC5B","GRM5","TENM3","COL6A3","GSTA5","SDK2","TRIP12","SLC26A8","CFAP47","AGO1","GRB14","ZP1","BRSK2","CTNND2","P2RX6","IL20RA","NUP35","STARD8","PCDHA9","GABRD","SYTL2","CLDN20","OR6C68","ADGRF2","ZRANB3","DCTN2","GOLGA6C","ATP6AP1L","CCDC148","SLC22A6","COL21A1","SH3BP4","OR2G6","PRSS23","AJM1","FAM227A","MCMBP","CCDC40","SIRPA","KRTAP7-1","OR6F1","ATP11A","UBOX5","ZFHX3","L3MBTL3","KRAS","KAZN","MS4A8","RBM10","COL24A1","ASB7","PERM1","NID2","SAFB2","COL18A1","SLC22A7","HMCN2","SHANK1","EBNA1BP2","COL6A5","TRIM49B","UBAP2L","TCERG1L","OLFML2B","SYN2","MUC16","EPHA3","REEP4","XPO5","IL36A","POSTN","NBPF9","SLC12A6","CASR","TTN","RFPL4AL1","PTF1A","DENND1A","TICRR","ERCC6","SLC9A4","CROCC2","TSPOAP1","ANKRD12","GPR101","ANKK1","STRA8","FOLH1","KRT36","PLCG2","OR9K2","PCDHB9","SLITRK3","NPR1","NFATC1","ZEB2","TRPC4","MYCT1","PTPRQ","PXDNL","RNF146","PERP","C17orf80","TAF3","GLRB","TENM2","USP19","TCF7L2","OR56A5","RABGAP1L","TPR","ADAMTS8","ZNF667","CCDC60","ZDHHC8","FAT4","CYP11B1","CD6","IRF2BP2","KCTD1","PATJ","ZNF284","ZNF469","LPAR3","KCNA4","INPP5E","KCNB2","PIK3C2B","SIGLEC1","ACACB","LRP1B","KLHL32","SAP25","KCTD16","DRC7","TUBB4B","TRAPPC8","DAW1","OR2W5","HTR1F","ZNF365","PCDH11X","ZFAND1","HSPA8","PCDHB12","AOC1","GNRHR","PCDHA1","FOXD1","MYO5B","TBC1D16","CEACAM20","CCDC170","ANKMY1","LYSMD4","VWDE")
FD4_LC2_data<-subset(data,idents = "T7")

data_FD4_LC2_silent<-rowMeans(data.frame(GetAssayData(FD4_LC1_data))[FD4_LC2_silent_no,])
mean(data_FD4_LC2_silent[!is.na(data_FD4_LC2_silent)])


#FD5

FD5_LC1_silent=c("FRG2C","TMEM151A","DENND3","OR4Q2","POU3F2","MUC16","NLRP5","CYP2W1","ZAN","AMIGO1","EEF1A1","EXOSC6","PCDHA2","ELOVL7","OTOP3","CTNNBL1")
FD5_LC1_data<-subset(data,idents = "T8")

data_FD5_LC1_silent<-rowMeans(data.frame(GetAssayData(FD5_LC2_data))[FD5_LC1_silent,])
mean(data_FD5_LC1_silent[!is.na(data_FD5_LC1_silent)])


FD5_LC1_silent_no=c("KMT2A","TMEM150B","ANKRD63","PLEKHA7","CHD4","TXNDC5","CLOCK","ANKRD36C","PCDHB7","MUC16","CCER1","ASH1L","DDX11","BRD9","TPRN","DEK","NDUFAF6","ESPN","GRIK5","VWCE","ASH2L","CERCAM","GPR180","ZAN","USP32","PARVG","LCE1E","RPS6KA2","CYP1A2","SRR","CHD1","SMPD1","SOWAHA","PLXNB3","ARAP1","FOXA2","EGFR","OBSCN","FOXE3","PRAMEF15","POLR3H","RIOK1","DST","KDM4A","ATAD2","PABPC1","SORBS1","RBM10")
FD5_LC1_data<-subset(data,idents = "T8")

data_FD5_LC1_silent<-rowMeans(data.frame(GetAssayData(FD5_LC2_data))[FD5_LC1_silent_no,])
mean(data_FD5_LC1_silent[!is.na(data_FD5_LC1_silent)])





FD5_LC2_silent=c("TUBB8B","HPD","AKR1B10","PKNOX1","FRG2C","TMEM200C","OR5AS1","KMT2D","FAM81A","SCAF1","DLGAP4","C20orf96","EP400","AKNA","FAT1","DTNB","NOP14","EXOSC6","EPHA2")
FD5_LC2_data<-subset(data,idents = "T9")

data_FD5_LC2_silent<-rowMeans(data.frame(GetAssayData(FD5_LC1_data))[FD5_LC2_silent,])
mean(data_FD5_LC2_silent[!is.na(data_FD5_LC2_silent)])



FD5_LC2_silent_no<-c("ITPR3","UNC80","RAET1L","HAUS7","CACNA1A","TYW1","GK2","NR4A3","BMP2K","SPDYE21","ARAP3","CRYBG2","RRP7A","CDH23","NLRP6","PDE4C","PIK3C3","PRR18","MAML3","RAVER2","SAMD11","TEX15","ENOX1","EGFR","ARNTL","CADM4","NANOG","KBTBD13","NDUFB7","CHD3","NEB","DCUN1D4","LIMK1","RIOK1","AGBL5","FLACC1")
FD5_LC2_data<-subset(data,idents = "T9")

data_FD5_LC2_silent<-rowMeans(data.frame(GetAssayData(FD5_LC1_data))[FD5_LC2_silent_no,])
mean(data_FD5_LC2_silent[!is.na(data_FD5_LC2_silent)])





#FD8

FD8_LC1_silent=c("GPD1L","CNIH2","MAST1","MYBBP1A","H3C1","SRGAP2","ASH1L","CTNNA2","OR4K17","P3H1","CDIPT","RUFY4","MAB21L2","AP2A1","OR10P1","SUGP2","DYRK2","ACIN1","BICD1","NANOGP8","PLPPR2","ALKBH6","FAT1","GTF2IRD1","SLC12A2","ST8SIA5","KCNU1","CEACAM1","FANCI","HSPA6","UROC1","TMEM63B","ABHD11","OR5F1","SIN3A","PDE1B","TRHR","SYNE1","ALPK2","ATP1A3","IDI2","SMG7","NELL2","RHOBTB1","NDUFA7","SYMPK","B3GNT4","PACSIN1","CRNKL1","RAD54L","RLF","ZNF285","OR52B6","AVIL","NR2F6","NCOA6","DENND4A","HSPG2","PCDHGA3","ADGRB1","ZNF845","C6orf118","NLRP12","PRR35","TXNDC16","PCDHB15","CA4","TEX45","RASSF5")
FD8_LC1_data<-subset(data,idents = "T12")

data_FD8_LC1_silent<-rowMeans(data.frame(GetAssayData(FD8_LC2_data))[FD8_LC1_silent,])
mean(data_FD8_LC1_silent[!is.na(data_FD8_LC1_silent)])


FD8_LC1_silent_no=c("PREX2","SEL1L2","HSPA5","HRNR","RPS6KL1","TMCO4","CPNE8","CCDC63","INTS1","PRICKLE1","PLEC","SLC26A5","GABRB2","CLASP1","PYGO1","CCNYL1","GSTK1","TTC13","MYH7B","CRYBG3","AKAP13","KCNK13","ABI3BP","FAM120C","GTPBP1","PIGR","MAEL","SLC33A1","SYP","ABI1","DPYS","WDR97","TAF6L","TPX2","ATP8B1","VWF","DPP9","ASCL1","SLC8A2","PCDHGC3","SEPTIN4","TLL2","POLQ","MAP3K13","RYR3","EHBP1L1","PLA2G15","ANP32E","SNX17","CIART","PKHD1L1","AGO2","ANKRD50","IGFBP7","VPS45","BTNL8","PATZ1","EMSY","SLC35A2","POMT2","SYNPO2","TTN","ZNHIT1","ESPL1","CRNKL1","UCN3","SNX27","CDH16","CEBPA","CENPE","TMEM131","MRTFB","HORMAD1","PDE3B","CAP1","CDH10","OSR2","CSMD3","RTP1","POFUT2","SRD5A1","TNS1","CCDC105","HGFAC","SIRT1","USH2A","MRGPRX1","ADGRB2","EGFR","NCAM1","NRBP2","KIFAP3","PCDHB12","WDR43","ITIH6","ZSWIM8","UPRT","FAM193A","KAZN","DGKH","GOLGA6L4","RAI2","SLC25A31","GAREM2","LRRC37B","MRPS31","FBRS","ABCC3","MYO9B","GCKR","WDR87","RPGRIP1","ZNF284","PPRC1","RASSF8","DNAH9","PCYT1A","CTLA4","ELFN1","ZNF672","BRD9","CCDC17","KIF4B","ADAR","A2ML1","DSG3","GTF3C1","SRSF1","TUBB2B","PRRC2B","BOP1","DEXI","IKBIP","IGFL3","VPS18","ZNF239","WDR90","ZNF446","HEATR3","ZFYVE16","SOX17","SLC15A1","LY9","PIK3CA","SMG7","ZNF71","AGAP9","TRIM49","USE1","ERBB2","KIF14","CLMP","TLX3","SHH","ZNF605","ATF7IP","LRP8","NLRP5","TEX15","NT5E")
FD8_LC1_data<-subset(data,idents = "T12")

data_FD8_LC1_silent<-rowMeans(data.frame(GetAssayData(FD8_LC2_data))[FD8_LC1_silent_no,])
mean(data_FD8_LC1_silent[!is.na(data_FD8_LC1_silent)])





FD8_LC2_silent=c("SLCO4A1","PRKD3","PCDHGB5","FUS","PDHA1","DENND4B","NTN5","CFAP53","SRGAP2","YIF1B","GOLGA6L7","DNAH5","FNDC1","SRRM2","HGFAC","CLEC17A","DUSP2","TNRC18","PSMC3","FOLH1","RANBP2","KLF17","NANOGP8")
FD8_LC2_data<-subset(data,idents = "T13")

data_FD8_LC2_silent<-rowMeans(data.frame(GetAssayData(FD8_LC1_data))[FD8_LC2_silent,])
mean(data_FD8_LC2_silent[!is.na(data_FD8_LC2_silent)])



FD8_LC2_silent_no<-c("MUC16","DOCK1","NBR1","TBP","PRKD2","ARHGAP6","KCNT2","ARID1B","HOOK3","C16orf92","SORL1","GLMP","F2RL1","EPPK1","KIAA0754","RYR1","EOMES","CLSTN2","TCF12","NET1","KLHDC3","MUC3A","OR8G5","CD276","MBD3","HLTF","TRIM49","CD200R1L","SORBS1","TDG","AHNAK2","PCDHA7","RECQL4","MPND","EGFR","CRHR2","GABRG1","BTNL3","SFTPB","ANKRD36","SEMA6A")
FD8_LC2_data<-subset(data,idents = "T13")

data_FD8_LC2_silent<-rowMeans(data.frame(GetAssayData(FD8_LC1_data))[FD8_LC2_silent_no,])
mean(data_FD8_LC2_silent[!is.na(data_FD8_LC2_silent)])






#FD9


FD9_LC1_silent=c("CAPNS1","ARID3A","PHLDB1","SCN8A","TIAM1","GET1","ZNF316","MYOF","OR5H1","LRFN3","IFNL2","FEZ2","NAV1","XIRP2","DSE","MYT1L","OR1J4","OR51M1")
FD9_LC1_data<-subset(data,idents = "T10")

data_FD9_LC1_silent<-rowMeans(data.frame(GetAssayData(FD9_LC2_data))[FD9_LC1_silent,])
mean(data_FD9_LC1_silent[!is.na(data_FD9_LC1_silent)])


FD9_LC1_silent_no=c("FBXO30","PDE6G","SIMC1","LRRC66","LOC101928120","TNRC18","RPS6KC1","ADAMTS18","ANKRD24","LRIF1","OR5H1","KIF18B","ZNF875","PPP1R13L","MCM7","CREB3L1","FGFBP3","TEX45","POU5F1B","NOL8","ZFP92","WDR62","FLG2","SLC44A5","KIF19","PTP4A2","TP53","EPHA7","SMARCAD1","TNIP1","IL27","FOXN4")
FD9_LC1_data<-subset(data,idents = "T10")

data_FD9_LC1_silent<-rowMeans(data.frame(GetAssayData(FD9_LC2_data))[FD9_LC1_silent_no,])
mean(data_FD9_LC1_silent[!is.na(data_FD9_LC1_silent)])





FD9_LC2_silent=c("MMP2","KMT5A","NPTXR","OR5H1","FOXJ2","ASTN2")
FD9_LC2_data<-subset(data,idents = "T11")

data_FD9_LC2_silent<-rowMeans(data.frame(GetAssayData(FD9_LC1_data))[FD9_LC2_silent,])
mean(data_FD9_LC2_silent[!is.na(data_FD9_LC2_silent)])



FD9_LC2_silent_no<-c("PSMB5","HSPG2","ERBB2","VWA8","EPB41L5","OSBPL3","OR5H1","LZTS2","NTN5","PCDHGA8","ANKRD36C","CD160")
FD9_LC2_data<-subset(data,idents = "T11")

data_FD9_LC2_silent<-rowMeans(data.frame(GetAssayData(FD9_LC1_data))[FD9_LC2_silent_no,])
mean(data_FD9_LC2_silent[!is.na(data_FD9_LC2_silent)])





#FD14

FD14_LC1_silent=c("RBM6","HADHA","DNMBP","TDG","MAMDC4","ANKS6","ELK1","CSF3R","ARSD","C6","CARTPT","KCNG3","HOXD13","FAM53C","GSK3B","POFUT1","KMT2D","GSN","ZDHHC11","SLC7A1","PAX5","ST3GAL3","MON1A","ADAM11","TMEM14B","CPXM2","MED30","ARID1B","OR2L2","CBLN1","OGFOD2","GRIK3","C20orf173","GRIN2B","FOXF2","ILF3","PGLYRP4","IRS2","GABRR2","MICAL1","PLEKHA4","GEMIN5","SHANK1","KLHL32","PERM1","ABCA7","STYXL2","RABGGTB","ZNF587","NUFIP1","TEX9","RBM23","FOXO3B","SP8","HTT","RYR3","TRMT5","SLC24A4","SLC38A10","ANKRD27","SPPL2B","CCDC62","FAM221B","FOXO4","F2R","WDR25","ARID5A","MICAL3","LMAN1L","XAF1","TFAP2C","FOXC1","BAHCC1","SALL1","CHIC1","PLXNC1","SLC12A4","ADAMTS6","MRPL51","GTF3C3","PC","NETO2","ZNF534","SLC7A11","YBX1","PRSS33","FDPS","PLPPR3","PDZD7","TEKT5","SYTL1","BRINP1","DNASE1L2","JAK3","MTARC1")
FD14_LC1_data<-subset(data,idents = "T14")

data_FD14_LC1_silent<-rowMeans(data.frame(GetAssayData(FD14_LC2_data))[FD14_LC1_silent,])
mean(data_FD14_LC1_silent[!is.na(data_FD14_LC1_silent)])


FD14_LC1_silent_no=c("LRIT3","IKBKE","MFSD14B","LTBP2","SENP3","OR51L1","ZNF75D","GSTM2","TMEM14B","CFAP65","ZNF646","HOXB4","COBL","ABCA4","DZIP1L","SNX17","TDG","AKAP8L","TNPO2","NCOA4","KIF3C","SLIT1","IFNA2","STX2","ARFGAP1","RDH8","GALNT14","SMG9","TOGARAM1","DYNC2I1","PRR12","MYO1G","EGFR","MFAP3L","AAAS","UTRN","BCL11A","ACOXL","IL2RA","HYAL4","CECR2","LAMA1","SOCS2","FOXJ1","PEG3","PDGFRL","VAV2","MYH15","ANO8","CDH11","NANOS1","KLHDC7A","PROSER3","TTN","MLX","TLL1","FMO4","YARS1","ECM1","C2CD3","NCF4","MRGPRX3","MMP25","CD53","PTPRB","ZNF586","POGK","MAGED1","PLEKHG2","KRTCAP3","ZFHX4","POLG","NOTCH2","LUC7L3","LOXHD1","ANKRD17","VWF","DOCK6","ARHGEF17","FOXE1","PTPRO","ILF3","RYR2","SHISA8","EEF1E1","PARD6G","HNRNPU","TRAF7","OR51A2","FMN2","NDRG4","IGDCC4","SMYD3","CACNA1G","CENPV","IRX3","AURKB","BCAS3","SLC2A3","CELSR3","ESAM","NENF","NPIPB11","ZMIZ2","PCDHAC1","FASN","CHD1","POLD1","ALOX15B","HTRA3","ANO2","SPTA1","ZDHHC11","GUCY2C","RARA","ANAPC1","SPTLC3","ZBTB47","UNC5C","SCTR","MAPK8IP2","WDFY4","VPS39","PLD4","TMEM199","OR13C5","RTL6","SIGIRR","POM121L12","XIRP1","TICAM1","ANKRD11","LARP6","GAS8","SRP68","AMDHD2","GATM","SEZ6L","ITPR3","TTLL4","C11orf95","SLC5A4","SNX29","VSTM2L","ATXN7","MS4A10","SOS1","NBPF12","APC2","HYAL2","SLC34A2","ZNF784","KCNQ4","CHI3L1","PLCD3","LNP1","PIK3R2","PIAS3","FAAP20","PTPRG","GPC2","BEST2","KDM7A","CCDC130","BIRC6","PRAMEF17","GRIN2D","QRICH1","NEDD9","C1QTNF1","SLC6A5","NFAT5","PRUNE2","RAP1GAP","NEUROD2","SGCG","LRRC71","ENTPD7","AUTS2","COL27A1","GPKOW","CKB","MIGA2","SAFB2","PTPRM","TEKT5","MCF2L","JPH2","PKN1","SIMC1","NCOR2","IFNW1","PER1","SLX4","AHNAK2","SP9","MDK","CHPF","TYW1","PPFIA2","MUC16","UNC93A","CYRIA","BRF1","KIRREL1","FGFR4","GRK2","EXT1","CASQ1","BICRA","DHX8","ABHD17C","CASKIN2","DNMT3B","OFD1","BTBD3","C17orf75","IQCA1L","SERGEF","PPP1R3F","LCT")
FD14_LC1_data<-subset(data,idents = "T14")

data_FD14_LC1_silent<-rowMeans(data.frame(GetAssayData(FD14_LC2_data))[FD14_LC1_silent_no,])
mean(data_FD14_LC1_silent[!is.na(data_FD14_LC1_silent)])





FD14_LC2_silent=c("ASS1","PARP12","GPA33","CLEC18B","UBAP2","CCDC102A","SPTB","ALPK3","SERPINC1","INCA1","IRGC","OR2D2","TRPM3","PLEKHG4","PLCE1","GTF2IRD2","SLFN14","NUP88","RNF185","MYPN","ANKRD9","EPHA2","ZNF407","TMEM14B","PLEKHM2","ACLY","IRF2BPL","RGS7","PLCG2","ARAP1","USP36","C3orf18","ZNF587","UBASH3A","CHD5","GHRL","UBB","CDH2","TNFAIP8L1","CDCP1","MBD6","APOL3","CCDC62","LRP8","UPF3A","KANSL3","MTHFR","ZNF316","MUC17","PHKB","WNK1","ZFP92","ARHGEF11")
FD14_LC2_data<-subset(data,idents = "T15")

data_FD14_LC2_silent<-rowMeans(data.frame(GetAssayData(FD14_LC1_data))[FD14_LC2_silent,])
mean(data_FD14_LC2_silent[!is.na(data_FD14_LC2_silent)])



FD14_LC2_silent_no<-c("GSE1","DEAF1","SDK1","PIK3CG","FBLN7","LLPH","PTPRJ","MROH1","FXR2","PTGER1","SLC5A12","TDG","THAP4","C2CD4C","CPAMD8","FAM126B","GIGYF2","SEC24B","FGFRL1","DOCK3","SLC23A3","SLC6A8","DCAF1","U2AF2","CBARP","ZP1","UBR5","JMJD6","APC2","QRICH2","CAMK2A","THG1L","STX11","BAZ2A","MSS51","ARC","STARD7","CASP10","ZNF773","IKZF4","COLGALT2","TEKT3","POM121","QPCTL","FSCB","ISLR2","CHD3","AFF4","PRAMEF17","TP63","ANAPC1","NR1H4","FBXO11","CDAN1","SH2D2A","AKR7A2","SLC2A3","KIAA0100","POLR2M","TTC12","XIRP2","DDRGK1","MYEOV","SPATA31D1","AUTS2","ARID5B","LCT","AHNAK","TANC2","PEBP1","MUC5B","LHCGR","EPB41L2","MYH14","NCK2","SLC35G4","URB2","CDC25A","DENND3","ZFYVE1","CDC42BPB","LOXL4","OR13C5","PRICKLE1","HTRA1","PHLDB1","MAGEA6","ASCL1","ADAM33","TRIL","TXNDC5","NBPF11","RNF150","MEX3C","GFER","CASTOR1","ZFPM1","CELSR1","PCM1","CACHD1","SLFN5","DAAM2","VAV1","RABEP1","CMA1","VCPKMT","FAM234A","CACNA2D2","AP3M1","TMEM14B","NHLRC3","IFITM1","SLC26A1","EHF","CLEC18B","FMO1","CD14","NBPF9","SYNPO2","YTHDF3","TMEM163","MOB1A","CDH15","TRPM2","COL7A1","MTFMT","DNAH11","RNF227","CRHR2","RARA","MEP1A","OBSCN","DMBX1","RBMXL3","DDX31","SEMA5A","MIEF1","EIF4A3","TEX14","GALNT16")
FD14_LC2_data<-subset(data,idents = "T15")

data_FD14_LC2_silent<-rowMeans(data.frame(GetAssayData(FD14_LC1_data))[FD14_LC2_silent_no,])
mean(data_FD14_LC2_silent[!is.na(data_FD14_LC2_silent)])




#FD16

FD16_LC1_silent=c("RTP2","STUB1","ONECUT3","POLR2I","KLF4","EIF2S3B","PLD5","STK32C","CXCR2","MYBPH","ABCD1","NPAS4","SHROOM4","RHOBTB2","EPB41L3","RNF213","PUM1","PORCN","VWA5B2","MYH7","DDX23","SETDB1","DNPEP","ZMYM4","SLC22A23","ANK1","KCNJ3","PRUNE2","DPP3","NHSL2","SHANK1","KPNB1","TAOK3","TNKS1BP1","COL1A1","PRDM9","TRIM2","FAM171A2","KCNQ1","FAM135B","SMARCD3","OPN4","TRIL","CIT","SUPT16H","MEGF8","SLC24A2","DNMT3A","EML3","SLC25A42","UAP1","NTN1","XRCC6","ABCA2","HS6ST3","CHD8","TMEM245","E2F3","BMP5","UBC","NEB","HMCN2","CCR10","CEP19","PHKG2","NCKIPSD","ANO4","ASAP2","SRD5A2","MTUS2","RFK","ACTR1A","GLI4","SEC23B","THBS1","RANBP10","MED12","NOTCH3","PCDH19","POU6F1","ZNF865","AACS","ALDH1L2","UNC93B1","TAMALIN","KDM3A","SPTAN1","KDELR3","CIAPIN1","DCAF7","AHNAK2","PCDH20","LOXL2","LRFN2")
FD16_LC1_data<-subset(data,idents = "T16")

data_FD16_LC1_silent<-rowMeans(data.frame(GetAssayData(FD16_LC2_data))[FD16_LC1_silent,])
mean(data_FD16_LC1_silent[!is.na(data_FD16_LC1_silent)])


FD16_LC1_silent_no=c("TOGARAM2","ATP13A1","RBMS3","CALR","MYH9","FOS","LDLRAD3","KDM3B","AHNAK","TBC1D30","CXCR2","SLC5A3","GLI4","ADGRB3","PABPC3","MYH7","EGFR","RNF113B","DTL","KRT10","PRR16","NUDT16L1","IQSEC2","KXD1","DNPEP","CTNNB1","HNRNPA2B1","SRSF9","PCDH20","SUN2","H1-4","KCNN4","COL1A1","GBP4","IFI30","POTEJ","DDX5","ROGDI","NCL","SLC12A3","AAK1","VAT1","ARMC9","VWA5B2","ZNF292","SPTAN1","CD3E","SSR2","SRD5A2","BRD4","ACP2","BTBD17","CT45A2","NKG7","ZFP36","PCDH19")
FD16_LC1_data<-subset(data,idents = "T16")

data_FD16_LC1_silent<-rowMeans(data.frame(GetAssayData(FD16_LC2_data))[FD16_LC1_silent_no,])
mean(data_FD16_LC1_silent[!is.na(data_FD16_LC1_silent)])





FD16_LC2_silent=c("SOX9","SPTBN4","ATP1A2","ESRRB","HYAL2","ARPC1B","TEAD4","LLGL2","MGAT3","PLPP3","ZNF703","PCDHGA10","RAP1GAP","FLAD1","KCNN2","MASP2","GRIPAP1","GSX1","TRIOBP","GAL3ST4","MKI67","SLC2A1","FOXL2","PCOLCE2","GALNT10","ADGRB2","RAP2A","TFE3","ZFHX3","PIK3R5","SON","HID1","MYO7A","NRXN2","MYCL","PPP1R3F","UBC","CTDNEP1","PDZD8","MYRFL","AHNAK2","HIC1","MEN1","COL4A5","ATG16L2","HDAC11","DIP2C","TIMM29","GRIA2","NOTCH2","ZNF574","KCNH1","FBN3","KIF1C","BAP1","ANKRD52","EIF3B","VIT","KCNA3","KANK4","ABCA1","NELFB","ACTN2","PDE7B","CELSR3","EVL","DLGAP3","CACNA1H","FBXO39","HERC2","MEX3B","CAMTA1","CLDN9","IQSEC2","ITPRID2","RHOV","SMARCD3","MAN2C1","CLEC16A","KCNS1","CAD","OXTR","DCLK1","ELAVL4","ZNF281","CD276","PPP2R5E","BTBD11","NRTN","MEIS1","HK2","NUMA1","CCDC171","GGT7","KLF2","CNNM2","WNT4","GARNL3","PIK3R2","FGD1","SGSM2","COL4A6")
FD16_LC2_data<-subset(data,idents = "T17")

data_FD16_LC2_silent<-rowMeans(data.frame(GetAssayData(FD16_LC1_data))[FD16_LC2_silent,])
mean(data_FD16_LC2_silent[!is.na(data_FD16_LC2_silent)])



FD16_LC2_silent_no<-c("CT45A1","TMEM202","E2F3","TATDN1","OXTR","TIMM29","AIPL1","MYL6B","VPS13C","RPRD2","KRT6B","CALR","NBEAL1","SCUBE3","FLNA","DDX5","ISM2","PCDHGA10","MAN2C1","TNRC18","BCL9L","EML5","FCGR2A","HYAL2","BCL2L12","FN1","CLEC4M","COL4A6","PRDM9","UPRT","CLEC16A","EGFR","MYLIP","MYCL","NBPF11","SLK","WNT2B","MKI67","AHNAK2","KCTD1")
FD16_LC2_data<-subset(data,idents = "T17")

data_FD16_LC2_silent<-rowMeans(data.frame(GetAssayData(FD16_LC1_data))[FD16_LC2_silent_no,])
mean(data_FD16_LC2_silent[!is.na(data_FD16_LC2_silent)])




#秩和检验

silent<-c(1.209121654,0.988264041,0.90131913,0.748170745,0.943735919,1.171395969,1.010454331,0.879134664,0.960056584,0.958864341,0.937850281,0.741663028,0.893168901,1.156047338,0.978205971,1.090418042)
non_silent<-c(1.023227259,0.953452563,0.909143393,0.627257208,0.925974667,1.108750346,1.260663493,0.798760982,0.941984614,1.126805002,0.811522695,0.655313108,0.800924862,1.146946534,0.935742424,1.038854908)
wilcox.test(silent,non_silent,paired = T,alternative = "greater")



mean(silent)
mean(non_silent)

median(silent)
median(non_silent)

silent<-c(0.924554974,0.942367507,1.182219641,0.805693694,0.803229099,1.133272342,1.138031278,0.657156341,0.862025114,1.223457362,1.086931403,0.770452687,0.836282014,1.17115365,1.089363796,0.980439259)
non_silent<-c(1.083200946,0.935641235,1.31090988,0.576816221,0.881240169,1.085718415,1.333901311,0.82400616,0.926880822,0.892498718,0.845361281,0.595146002,0.810525534,1.117294712,1.047000984,0.85156179)
wilcox.test(silent,non_silent,paired = T,alternative = "greater")





WES<-data.frame(matrix(NA,32,2))
colnames(WES)<-c("Type","Index")
WES[1:16,]$Type="silent"
WES[17:32,]$Type="non_silent"


WES[1:16,]$Index=silent
WES[17:32,]$Index=non_silent

WES$Type<-factor(WES$Type,levels = c("silent","non_silent"))

library(ggplot2)
WES_new<-WES[c(1:32),]
p<-ggplot(data=WES_new)+
  geom_point(aes(x=factor(Type),y=Index,col=factor(Type)),size=2)+
  ylab("Gene expression rate")+
  ylim(0.6,1.3)+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p

ggsave(p,file="/Users/wangjun/Desktop/最新版结果/WES表达量绘图/WES_total.png",dpi=1000,width = 2.12,height = 3.44)

for (i in 1:16){
  WES_new<-WES[c(i,16+i),]
  p<-ggplot(data=WES_new)+
    geom_point(aes(x=factor(Type),y=Index,col=factor(Type)),size=2)+
    ylab("Gene expression rate")+
    ylim(0.6,1.3)+
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
    scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
  ggsave(p,file=paste0("/Users/wangjun/Desktop/最新版结果/WES表达量绘图/WES_",i,".png"),dpi=1000,width = 2.12,height = 3.44)
  
}



silent<-c(0.1757919,0.1428508,0.06837187,0.1236825,0.0608689,0.08851123,0.354885,0.08166332,0.05994096,0.1784226,0.08780944,0.03867367,0.1007641,0.1276341,0.1356129,0.1430261)
non_silent<-c(0.1301674,0.07532649,0.1004688,0.05262555,0.07416145,0.07159976,0.1717978,0.0830725,0.08665259,0.09296912,0.05880304,0.09935569,0.08297958,0.1195959,0.4492534,0.2023653)
wilcox.test(silent,non_silent,paired = T,alternative = "greater")







#T细胞差异亚群计算

CD8_C1<-c(0.025356577,0.039485767	,0.031372549,	0.039423485	,0.045150502,	0.083958021	,0.104112194,	0.010674931	,0.015396076,	0.044881493	,0.033906177,	0.066225166,	0.036458333	,0.12936345	,0.202592087,	0.086126629	,0.10096463)


#配对分析结果

wilcox.test(CD8_C1[c(3,4,6,8,12,14)],CD8_C1[c(2,5,7,9,13,15)],paired = T,alternative = "less")

#分组分析结果
wilcox.test(CD8_C1[c(3,4,6,8,12,14,10,11,16,17)],CD8_C1[c(2,5,7,9,13,15)],alternative = "less")




CD8_C2<-c(0.216323296,	0.060606061	,0.019047619	,0.008902077	,0.128762542	,0.074212894	,0.057047778	,0.00137741	,0.015147753	,0.099344428	,0.080352996	,0.034216336	,0.026041667	,0.033264887,	0.003410641	,0.012104283,	0.004501608)
#配对分析结果

wilcox.test(CD8_C2[c(3,4,6,8,12,14)],CD8_C2[c(2,5,7,9,13,15)],paired = T,alternative = "less")

#分组分析结果
wilcox.test(CD8_C2[c(3,4,6,8,12,14,10,11,16,17)],CD8_C2[c(2,5,7,9,13,15)],alternative = "less")


CD8_C4<-c(0.050713154,	0.053259871	,0.058263305	,0.024162781,	0.078595318,	0.049975012,	0.12051343,	0.097107438,	0.11522225,	0.042864347,	0.044588946,	0.087196468	,0.040798611,	0.029568789	,0.016371078,	0.053072626,	0.04437299)
#配对分析结果

wilcox.test(CD8_C4[c(3,4,6,8,12,14)],CD8_C4[c(2,5,7,9,13,15)],paired = T,alternative = "less")

#分组分析结果
wilcox.test(CD8_C4[c(3,4,6,8,12,14,10,11,16,17)],CD8_C4[c(2,5,7,9,13,15)],alternative = "less")



#*
CD8_C5<-c(0.064976228,	0.282828283,	0.003921569,	0.001695634,	0.087792642,	0.00149925,	0.002376991,	0.001033058,	0.019369258,	0.025214322,	0.010682768,	0.001103753	,0.006944444	,0.026694045,	0.003410641,	0.000465549	,0)
#配对分析结果

wilcox.test(CD8_C5[c(3,4,6,8,12,14)],CD8_C5[c(2,5,7,9,13,15)],paired = T,alternative = "less")
#分组分析结果
wilcox.test(CD8_C5[c(3,4,6,8,12,14,10,11,16,17)],CD8_C5[c(2,5,7,9,13,15)],alternative = "less")



#*
CD8_C6<-c(0.000792393,	0.013774105	,0.001680672,	0.005934718,	0.036789298,0.00049975	,0.001663894,	0.005853994,	0.034020363,	0	,0.000464468	,0.001103753,	0.000868056	,0.058316222	,0.036152797,	0	,0)
#配对分析结果

wilcox.test(CD8_C6[c(3,4,6,8,12,14)],CD8_C6[c(2,5,7,9,13,15)],paired = T,alternative = "less")
#分组分析结果
wilcox.test(CD8_C6[c(3,4,6,8,12,14,10,11,16,17)],CD8_C6[c(2,5,7,9,13,15)],alternative = "less")






CD4_C1<-c(0.250396197	,0.19651056,	0.393277311,	0.148791861,	0.425585284	,0.55822089	,0.585928215,	0.080922865	,0.136329774,	0.461422088	,0.671620994	,0.636865342	,0.731770833	,0.525256674	,0.645975443	,0.640130354	,0.719614148)
#配对分析结果

wilcox.test(CD4_C1[c(3,4,6,8,12,14)],CD4_C1[c(2,5,7,9,13,15)],paired = T,alternative = "less")
#分组分析结果
wilcox.test(CD4_C1[c(3,4,6,8,12,14,10,11,16,17)],CD4_C1[c(2,5,7,9,13,15)],alternative = "less")




CD4_C3<-c(0.021394612	,0.056932966	,0.080112045,	0.030097499,	0.063545151,	0.051474263,	0.021155217,	0.387052342,	0.162900422,	0.033282905,	0.021365536,	0.05187638,	0.039930556,	0.003285421,	0.002046385,	0.011173184,	0.017363344)
#配对分析结果

wilcox.test(CD4_C3[c(3,4,6,8,12,14)],CD4_C3[c(2,5,7,9,13,15)],paired = T,alternative = "greater")
#分组分析结果
wilcox.test(CD4_C3[c(3,4,6,8,12,14,10,11,16,17)],CD4_C3[c(2,5,7,9,13,15)],alternative = "greater")





CD4_C4<-c(0.003169572,	0.050505051	,0.219047619,0.328529038	,0.022575251,	0.004997501,	0.004040884	,0.136707989	,0.231437795	,0.002017146,	0	,0.017660044	,0.029513889	,0.002053388	,0.008185539	,0,	0.003858521)
#配对分析结果

wilcox.test(CD4_C4[c(3,4,6,8,12,14)],CD4_C4[c(2,5,7,9,13,15)],paired = T,alternative = "greater")
#分组分析结果
wilcox.test(CD4_C4[c(3,4,6,8,12,14,10,11,16,17)],CD4_C4[c(2,5,7,9,13,15)],alternative = "less")



CD4_C5<-c(0.000792393,0.004591368,0.002240896,0.001695634,0.00083612,0.027486257,0.001663894,0.005165289,0.001489943,0.002017146,0.00510915,0.013245033,0.024305556,0.003285421,0.008185539,0.004189944,0.002572347)
#配对分析结果

wilcox.test(CD4_C5[c(3,4,6,8,12,14)],CD4_C5[c(2,5,7,9,13,15)],paired = T,alternative = "greater")
#分组分析结果
wilcox.test(CD4_C5[c(3,4,6,8,12,14,10,11,16,17)],CD4_C5[c(2,5,7,9,13,15)],alternative = "less")


CD4_C6<-c(0,0.001836547,0.010644258,0.002543451,0.016722408,0.005747126,0.009032565,0,0,0.002521432,0.007431491,0.014348786,0.006076389,0.002053388,0.000682128,0.000931099,0.008360129)
#配对分析结果

wilcox.test(CD4_C6[c(3,4,6,8,12,14)],CD4_C6[c(2,5,7,9,13,15)],paired = T,alternative = "greater")
#分组分析结果
wilcox.test(CD4_C6[c(3,4,6,8,12,14,10,11,16,17)],CD4_C6[c(2,5,7,9,13,15)],alternative = "less")


CD4_C7<-c(0.003961965,0.024793388,0.03697479,0.245866893,0.022575251,0.001249375,0.003803185,0.187327824,0.171591756,0,0.000464468,0.007726269,0.005208333,0.002874743,0.006139154,0,0.001286174)
#配对分析结果

wilcox.test(CD4_C7[c(3,4,6,8,12,14)],CD4_C7[c(2,5,7,9,13,15)],paired = T,alternative = "greater")
#分组分析结果
wilcox.test(CD4_C7[c(3,4,6,8,12,14,10,11,16,17)],CD4_C7[c(2,5,7,9,13,15)],alternative = "less")






#读取T细胞数据
T_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/T_total_12_10.rds")
table(Idents(T_cells))

#鉴定显著的T细胞亚群


#髓系细胞
library(Seurat)
Myeloid_cells<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/macrophages_update.rds")
table(Idents(Myeloid_cells))


levels(Myeloid_cells)

Idents(Myeloid_cells)<-factor(Idents(Myeloid_cells),levels =c("Alveolar resident Macrophage","Perivascular resident Macrophage","Anti-inflammatory macrophages","Proliferating Macrophage","Monocytes Derived DC","Migratory Conventional Dendritic cell","Conventional Dendritic cell type 1","Conventional Dendritic cell type 2","Classical CD14+ Monocyte","Non-classical CD16+ Monocyte") )

library(paletteer) 
library(scales)
pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal[c(2,3,7,10,5,8,4,1,9)]

#pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,5,8,7,1,9)] 
pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,10,5,8,4,1,9)] 

pal_new
pal_new<-(colorRampPalette(pal)(16))



df = data.frame(clu=names(table(Idents(Myeloid_cells))),
                per=sprintf("%1.2f%%", 100*table((Idents(Myeloid_cells)))/length(Idents(Myeloid_cells))))

Myeloid_cells$per = df[match(Idents(Myeloid_cells),df$clu),2]
Myeloid_cells$new = paste0(Idents(Myeloid_cells),"(",Myeloid_cells$per,")")
table(Myeloid_cells$new)


DimPlot(sce,reduction = 'umap',  
        group.by = 'new',
        label.box = T,  label = F,repel = T)


Myeloid_cells$new<-factor(Myeloid_cells$new,levels = c("Monocytes Derived DC(30.48%)","Classical CD14+ Monocyte(17.72%)","Alveolar resident Macrophage(12.91%)","Perivascular resident Macrophage(12.06%)","Anti-inflammatory macrophages(8.75%)","Non-classical CD16+ Monocyte(8.65%)","Conventional Dendritic cell type 2(3.45%)","Proliferating Macrophage(2.64%)","Migratory Conventional Dendritic cell(1.69%)","Conventional Dendritic cell type 1(1.66%)"))


Idents(Myeloid_cells)<-"orig.ident"

Myeloid_cells<-subset(Myeloid_cells,idents = "T3",invert=T)


library(ggplot2)
#TNSE
Myeloid_cells <- RunTSNE(Myeloid_cells, reduction = "pca", dims = 1:20)
p5 <- DimPlot(Myeloid_cells, reduction = "tsne",cols = pal[1:10])
p5

Idents(Myeloid_cells)<-"new"

ggsave(p5,file="./Myeloid_final.png",dpi = 1000,width = 6.97,height =3.99 )

p5 <- DimPlot(Myeloid_cells, reduction = "tsne",split.by = "orig.ident",cols = pal[1:10])+
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title.x=element_text(size=5,face = "bold"),axis.title.y=element_text(size=5,face = "bold"),axis.text.x=element_text(size=5,face = "bold",color = "black"),axis.text.y=element_text(size=5,face = "bold",color = "black"),legend.box.background=element_rect(colour = "white",fill = "white",size = 5))

p5
ggsave(p5,file="/Users/wangjun/Desktop/最新版结果/Myeloid_orig_ident_update.png",dpi = 1000,width = 23,height =2)



p5 <- DimPlot(Myeloid_cells, reduction = "tsne",group.by = "new",cols = pal[1:10])

p5
ggsave(p5,file="/Users/wangjun/Desktop/最新版结果/Myeloid_update.png",dpi = 1000,width = 8.62,height =4.72)











#绘制CD8-C5和CD8-C6的分组图

#CD8_C5

CD8_C5_MIA<-c(0.282828283,0.087792642,0.002376991,0.019369258,0.006944444,0.003410641)
CD8_C5_LUAD<-c(0.003921569,0.001695634,0.00149925,0.001033058,0.025214322,0.010682768,0.001103753,0.026694045,0.000465549,0)
wilcox.test(CD8_C5_MIA,CD8_C5_LUAD)

length(CD8_C5_LUAD)


#CD8_C6
CD8_C6_MIA<-c(0.013774105,0.036789298,0.001663894,0.034020363,0.000868056,0.036152797)
CD8_C6_LUAD<-c(0.001680672,0.005934718,0.00049975,0.005853994,0,0.000464468,0.001103753,0.058316222,0,0)
wilcox.test(CD8_C6_MIA,CD8_C6_LUAD)






CD8_C5<-data.frame(matrix(NA,17,2))
colnames(CD8_C5)<-c("Samples","Proportion")
CD8_C5[1:7,]$Samples="MIA"
CD8_C5[8:17,]$Samples="LUAD"


CD8_C5[1:7,]$Proportion=CD8_C5_MIA
CD8_C5[8:17,]$Proportion=CD8_C5_LUAD

CD8_C5$Samples<-factor(CD8_C5$Samples,levels = c("MIA","LUAD"))


#CD8_C5_new<-CD8_C5[c(7,14),]
CD8_C5_new<-CD8_C5
p<-ggplot(data=CD8_C5_new)+
  geom_point(aes(x=factor(Samples),y=Proportion,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Proportion")+
  #ylim(0,0.05)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/最终方案/figure/CD8_C5/FD14_LUAD_MIA.png",dpi = 1000,width = 2.23,height = 3.44)






#CDS下调基因富集展示
#开始做GO-KEGG的通路富集展示图
data_go<-read.table(file="/Users/wangjun/Desktop/monocle_down_Go_select.csv",header = T,sep=",")


data_go<-data_go[data_go$select%in%1,]

data_go



#开始作图
colnames(data_go)

#df1<-data_kegg[,c(2,3,4,6,9)]
df1<-data_go[,c(3,4,5,7,10)]


df1<-df1[order(-df1$Count),]


df1$`-log10(P_value)`=-log10(as.numeric(df1$p.adjust))

df1$GeneRatio
df1$BgRatio

(7/37)/(480/18723)
(5/37)/(336/18723)
(5/37)/(487/18723)
(4/37)/(112/18723)
(4/37)/(188/18723)
(4/37)/(284/18723)
(3/37)/(109/18723)
(3/37)/(142/18723)
(2/37)/(24/18723)
(2/37)/(27/18723)
(2/37)/(36/18723)


df1$Count




df1$enrich_factor=c(7.38,7.53,5.20,18.07,10.77,7.13,13.93,10.69,42.17,37.48,28.11)


head(df1)
df2<-df1[,c(1,5,6,7)]

colnames(df2)<-c("Pathway","Count","-log10(P_value)","Enrich_Factor")
df2$Count
df2$Count<-as.numeric(df2$Count)


p1<-ggplot(df2,aes(`Enrich_Factor`,Pathway) )+ 
  geom_point(aes(color=`-log10(P_value)`,size= `Count`) ) +
  scale_size_continuous(range = c(4,8))+
  scale_color_gradient2(high="red",mid = "#CC3300",low ="blue",midpoint = 1.75 )+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 


p1


ggsave(p1,file="/Users/wangjun/Desktop/图片合集材料包/monocle_down_go.png",dpi=1000,width = 7.15,height = 4.22 )




library(BiocManager)
BiocManager::install("Rgraphviz")
install.packages("Rgraphviz")


library(Rgraphviz)
plotGOgraph(GO2,firstSigNodes = 1,useFullNames = TRUE)#GO-BP功能网络图



GO2


#通路热图，小样本适合，大样本基因太多
enrichplot::heatplot(GO2@result,showCategory = 100,label_format = 10)#基因-通路关联热图

?enrichplot::heatplot

install.packages("ggnewscale")
library("ggnewscale")

p1<-enrichplot::emapplot(pairwise_termsim(GO2),
                     showCategory = 30, 
                     color = "p.adjust", 
                     layout = "kk",
                     nWords = 1,
                     cex_label_group = 4,
                     label_format = 50)#通路间关联网络图

ggsave(p1,file="/Users/wangjun/Desktop//emapplot_go_down.png",dpi=1000,height = 12,width = 12)


pairwise_termsim(GO2@result)

GO2@result$Description



subset(GO2)




enrichplot::emapplot(pairwise_termsim(GO2))

?enrichplot::emapplot


?enrichplot::emapplot

#可以选择特定的通路在KEGG中打开
browseKEGG(kk,"hsa00360")#选择其中的hsa05166通路进行展示

#绘制GO三个板块的功能网，但是并不是很实用。
GO3<-enrichGO( names(geneList),#GO富集分析
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",#设定读取的gene ID类型
               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 0.05,#设定p值阈值
               #qvalueCutoff = 0.05,#设定q值阈值
               pAdjustMethod = "BH",
               readable = TRUE)

library(Rgraphviz)
plotGOgraph(GO3,firstSigNodes = 1,useFullNames = TRUE)#GO-BP功能网络图




#计算特异性macro亚群的MIA和LUAD的差异基因。



#读取抗原亲和力较高的TCR对应的barcode
FD1_LC2_barcode_low<-read.table(file="/Users/wangjun/Desktop/sc新抗原分析结果v1/neo_matrix_panpep/snv/FD1_LC2_barcode_low.csv",sep = ",")

FD1_LC2_barcode_low[1:4,]

dim(FD1_LC2_barcode_low)

T_total<-readRDS(file="/Users/wangjun/Desktop/T_total_3_31.rds")

table(T_total$Annotation)

Idents(T_total)<-"orig.ident"

T_FD1_LC2<-subset(T_total,idents = "T4")

anno<-data.frame(T_FD1_LC2$Annotation)

my_list <- strsplit(rownames(anno), "_")

for (i in 1:length(rownames(anno))){
  print(i)
  rownames(anno)[i]<-my_list[[i]][3]
}
anno$barcode<-rownames(anno)


table(anno[anno$barcode%in%FD1_LC2_barcode_low$V1,]$T_FD1_LC2.Annotation)
dim(FD1_LC2_barcode_low)

nchar("CTTACCGAGCTAAGAT-1")



colnames(T_total)


anno<-data.frame(T_total$orig.ident)

my_list <- strsplit(rownames(anno), "_")

table(nchar(anno$barcode))

tail(anno$barcode)

anno$barcode<-rownames(anno)
anno[1:3,]
dim(anno)

for (i in 1:length(rownames(anno))){
  print(i)
  if (nchar(anno$barcode[i])==18){
    anno$barcode[i]<-anno$barcode[i]
  }
  else{
    anno$barcode[i]<-my_list[[i]][3]
  }
}

length(anno$barcode)
dim(T_total)


anno$barcode<-rownames(anno)

length(unique(anno$barcode))
length(anno$barcode)

barcode<-unique(anno$barcode)

write.table(barcode,file="/Users/wangjun/Desktop/barcode.csv",sep = ",")




#新抗原统计
T_total<-readRDS(file="/Users/wangjun/Desktop/T_total_3_31.rds")
data<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_snv_index_result_new_best.csv",sep = ",",header = T)

data[1:4,1:4]

MIA<-c("FD1_LC2","FD2_LC2","FD4_LC2","FD5_LC2","FD9_LC2","FD14_LC2")
LUAD<-c("FD1_LC3","FD2_LC1","FD4_LC1","FD5_LC1","FD8_LC1","FD8_LC2","FD9_LC1","FD14_LC1","FD16_LC1","FD16_LC2")

MIA_neo<-table(data[data$SampleID%in%MIA,]$SampleID)

LUAD_neo<-table(data[data$SampleID%in%LUAD,]$SampleID)





MIA_neo<-c(3,29,59,7,2,0)
LUAD_neo<-c(8,36,5,2,5,83,16,30,11,6)

LUAD_neo<-c(8,36,83,16,6,5)

wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)



#强结合
table(data[data$SampleID%in%MIA & data$tumor_Score_.rank<0.5,]$SampleID)
table(data[data$SampleID%in%LUAD & data$tumor_Score_.rank<0.5,]$SampleID)
wilcox.test(c(9,15,1),c(9,18,3),alternative = "less",paired = T)
#V = 0, p-value = 0.1855


#弱结合
table(data[data$SampleID%in%MIA & data$tumor_Score_.rank>0.5,]$SampleID)
table(data[data$SampleID%in%LUAD & data$tumor_Score_.rank>0.5,]$SampleID)
wilcox.test(c(3,20,44,6,2),c(5,27,65,13,4),alternative = "less",paired = T)
#V = 0, p-value = 0.02838

Idents(T_total)

#读取和TCR对应的barcode
sample="FD2_LC1"
sample_T="T1"

sample_list=c("FD1_LC2","FD1_LC3","FD2_LC1","FD4_LC1","FD4_LC2","FD5_LC1","FD5_LC2","FD8_LC1","FD8_LC2","FD9_LC1","FD9_LC2","FD14_LC1","FD14_LC2","FD16_LC1","FD16_LC2")
sample_T_list=c("T4","T5","T1","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")

dim(anno)

library(Seurat)
Idents(T_total)<-"orig.ident"


anno[1:4,]


data[1:4,]


for (i in 1:length(sample_list)){
#for (i in 1:1){
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=c("LAG3","TIGIT","PDCD1","CTLA4","KLRG1")
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  col_n<-in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode
  print(length(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3)))
  print(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3))
  T_total_set<-subset(T_total,idents = sample_T)
  print(length(rownames(anno[anno$barcode%in%col_n & anno$T_total.orig.ident%in%sample_T,])))
  if (length(rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))!=1) {
  score<-(colMeans(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set)))))
  print(paste0(sample,"-",length(score),"-",mean(score)))}
  else{
    score<-(mean(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set)))))
    print(paste0(sample,"-","1","-",score))
    }
  }


print(in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR)


#存下每个细胞的耗竭得分
score_list<-list()

for (i in 1:length(sample_list)){
  #for (i in 1:1){
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=c("LAG3","TIGIT","PDCD1","CTLA4","KLRG1")
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  col_n<-in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode
  print(length(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3)))
  print(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3))
  T_total_set<-subset(T_total,idents = sample_T)
  if (length(rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))!=1) {
    score_list[[i]]<-colMeans(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set))))
    print(paste0(sample,"-",length(score_list[[i]]),"-",score_list[[i]]))}
  else{
    score_list[[i]]<-colMeans(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set))))
    print(paste0(sample,"-","1","-",score_list[[i]]))
  }
}











mean(colMeans(data.frame(GetAssayData(T_total_set))))



data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))]

dim(anno)

rownames(anno[anno$barcode%in%col_n ,])

in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode


col_n<-in_put3[in_put3$TRB_CDR3%in%"CASSYSGLAGETQYF",]$barcode

rownames(anno[anno$barcode%in%col_n,])
anno[anno$barcode%in%"ATAGACCAGCCGTCGT-1",]

anno$barcode[1:4]


MIA_score<-c(3.59758146888463,1.96071428301553,2.52069350339907,1.70044438570224,3.82881163796315)
LUAD_score<-c(3.0829921708233,6.76177368116464,6.24942286457179,3.18344116377111,5.57894295254999,1.17585726928332,2.05387325108412,5.73325784688997,3.3671847456927,4.87664944379975)

LUAD_score<-c(3.0829921708233,6.76177368116464,6.24942286457179,3.18344116377111,5.57894295254999,1.17585726928332,2.05387325108412,5.73325784688997)

wilcox.test(MIA_score,LUAD_score,alternative = "less",paired = T)
















#新抗原INDEL统计
#新抗原统计
data<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_indel_index_result_new_best.csv",sep = ",",header = T)

table(data$SampleID)

MIA<-c("FD1_LC2","FD4_LC2")
LUAD<-c("FD2_LC1","FD5_LC1","FD8_LC2","FD14_LC1","FD16_LC1")

MIA_neo<-table(data[data$SampleID%in%MIA,]$SampleID)

LUAD_neo<-table(data[data$SampleID%in%LUAD,]$SampleID)

LUAD_neo<-c(8,36,83,16,6)

wilcox.test(MIA_neo,LUAD_neo,alternative = "less")

#强结合
table(data[data$SampleID%in%MIA & data$tumor_Score_.rank<0.5,]$SampleID)
table(data[data$SampleID%in%LUAD & data$tumor_Score_.rank<0.5,]$SampleID)
wilcox.test(c(9,15,1),c(9,18,3),alternative = "less",paired = T)
#V = 0, p-value = 0.1855


#弱结合
table(data[data$SampleID%in%MIA & data$tumor_Score_.rank>0.5,]$SampleID)
table(data[data$SampleID%in%LUAD & data$tumor_Score_.rank>0.5,]$SampleID)
wilcox.test(c(3,20,44,6,2),c(5,27,65,13,4),alternative = "less",paired = T)
#V = 0, p-value = 0.02838




MIA<-c("FD1_LC2","FD4_LC2")
LUAD<-c("FD2_LC1","FD5_LC1","FD8_LC2","FD14_LC1","FD16_LC1")

#存下每个细胞的耗竭得分
score_list2<-list()

for (i in 1:length(sample_list)){
  #for (i in 1:1){
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=c("LAG3","TIGIT","PDCD1","CTLA4","HAVCR2","ENTPD1")
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  col_n<-in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode

  print(length(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3)))
  print(unique(in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$TRB_CDR3))
  T_total_set<-subset(T_total,idents = sample_T)
  
  print(length(rownames(anno[anno$barcode%in%col_n & anno$T_total.orig.ident%in%sample_T,])))
    if (length(rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))!=1) {
    score_list2[[i]]<-colMeans(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set))))
    print(paste0(sample,"-",length(score_list2[[i]]),"-",score_list2[[i]]))}
  else{
    score_list2[[i]]<-mean(data.frame(GetAssayData(T_total_set))[gene_list,gsub("-",".",rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))])/mean(colMeans(data.frame(GetAssayData(T_total_set))))
    print(paste0(sample,"-","1","-",score_list2[[i]]))
  }
}



MIA<-c(1,5,7,11,13)
LUAD<-c(2,3,4,6,8,9,10,12,14,15)





MIA_score<-c(score_list[[1]],score_list[[5]],score_list[[7]],score_list[[11]],score_list[[13]],
             score_list2[[1]],score_list2[[5]],score_list2[[7]],score_list2[[11]],score_list2[[13]])

LUAD_score<-c(score_list[[2]],score_list[[3]],score_list[[4]],score_list[[6]],score_list[[8]],score_list[[9]],score_list[[10]],score_list[[12]],score_list[[14]],score_list[[15]],
              score_list2[[2]],score_list2[[3]],score_list2[[4]],score_list2[[6]],score_list2[[8]],score_list2[[9]],score_list2[[10]],score_list2[[12]],score_list2[[14]],score_list2[[15]])


wilcox.test(MIA_score,LUAD_score,alternative = "less")


library(lsr)
cohensD(MIA_score,LUAD_score)




length(MIA_score)+length(LUAD_score)



data<-data.frame(matrix(NA,802,2))
colnames(data)<-c("Exhausted_Score","Type")

data[1:243,]$Exhausted_Score<-MIA_score
data[1:243,]$Type<-"MIA"
data[244:802,]$Exhausted_Score<-LUAD_score
data[244:802,]$Type<-"LUAD"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))



p<-ggplot(data, aes(x=Type, y=Exhausted_Score, fill=Type)) + 
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(title = "Violin_Plot", x = "Tissue", y = "Exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/volin_plot.png",dpi = 1000,height = 6.37,width = 5.14)



length(score_list)
length(score_list2)

for (i in 1:15){
  score<-c(score_list[[i]],score_list2[[i]])
  cat(mean(score))
  cat("\n")
}




#TCR克隆绘图
MIA_neo<-c(3.582632,0,2.045445,2.520694,1.700444,3.828812)
LUAD_neo<-c(3.082992,6.295444,6.249423,3.126719,5.578943,0.9877201,2.053873,5.500828,3.120805,4.876649)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less")



MIA_neo<-c(3.582632,0,2.045445,2.520694,1.700444,3.828812)
LUAD_neo<-c(3.082992,6.295444,6.249423,3.126719,2.053873,5.500828)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)



Neoantigen<-data.frame(matrix(NA,16,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:16,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:16,]$Neoantigen_count=as.numeric(LUAD_neo)

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))




Neoantigen<-data.frame(matrix(NA,12,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:12,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:12,]$Neoantigen_count=as.numeric(LUAD_neo[c(1:5,8)])

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))



library(ggplot2)

Neoantigen_new<-Neoantigen[c(6,12),]
#Neoantigen_new<-Neoantigen
p<-ggplot(data=Neoantigen_new)+
  geom_point(aes(x=factor(Samples),y=Neoantigen_count,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("TCR_clone_count")+
  ylim(0,8)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/Exhausted_score/6.png",dpi = 1000,width = 2.12,height = 3.44)

#




















MIA_neo<-c(7,0,61,7,2,29)
LUAD_neo<-c(8,7,83,17,6,37)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)



MIA_neo<-c(7,0,61,7,2,29)
LUAD_neo<-c(8,7,83,17,6,37,30,12,6,2)


MIA_neo<-c(7,0,61,7,2,29)
LUAD_neo<-c(8,7,83,17,30,12,6,37,6,2)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less")





#新抗原绘图

Neoantigen<-data.frame(matrix(NA,16,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:16,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:16,]$Neoantigen_count=as.numeric(LUAD_neo)

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))




Neoantigen<-data.frame(matrix(NA,12,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:12,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:12,]$Neoantigen_count=as.numeric(LUAD_neo[c(1:5,8)])

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))



library(ggplot2)

Neoantigen_new<-Neoantigen[c(1,7),]
#Neoantigen_new<-Neoantigen
p<-ggplot(data=Neoantigen_new)+
  geom_point(aes(x=factor(Samples),y=Neoantigen_count,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("Neoantigen_count")+
  ylim(0,100)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/Neoantigen_count/1.png",dpi = 1000,width = 2.12,height = 3.44)











#TCR克隆绘图

MIA_neo<-c(8,0,119,33,10,73)
LUAD_neo<-c(48,29,154,166,7,25,10,74,41,5)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less")


MIA_neo<-c(8,0,119,33,10,73)
LUAD_neo<-c(48,29,154,166,10,74)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)





MIA_neo<-c(0.006546645,0,0.027921164,0.008562532,0.008920607,0.061396131)
LUAD_neo<-c(0.027858387,0.016347238,0.03836572,0.070160609,0.002310231,0.011286682,0.011778563,0.028040925,0.016646366,0.002736727)


MIA_neo<-c(0.006546645,0,0.027921164,0.008562532,0.008920607,0.061396131)
LUAD_neo<-c(0.027858387,0.016347238,0.03836572,0.070160609,0.011778563,0.028040925)

wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)






Neoantigen<-data.frame(matrix(NA,16,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:16,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:16,]$Neoantigen_count=as.numeric(LUAD_neo)

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))




Neoantigen<-data.frame(matrix(NA,12,2))
colnames(Neoantigen)<-c("Samples","Neoantigen_count")
Neoantigen[1:6,]$Samples="MIA"
Neoantigen[7:12,]$Samples="LUAD"


Neoantigen[1:6,]$Neoantigen_count=as.numeric(MIA_neo)
Neoantigen[7:12,]$Neoantigen_count=as.numeric(LUAD_neo[c(1:5,8)])

Neoantigen$Samples<-factor(Neoantigen$Samples,levels = c("MIA","LUAD"))



library(ggplot2)
#Neoantigen_new<-Neoantigen[c(1,7),]
Neoantigen_new<-Neoantigen
p<-ggplot(data=Neoantigen_new)+
  geom_point(aes(x=factor(Samples),y=Neoantigen_count,col=factor(Samples)),size=2)+
  xlab("Tissue")+
  ylab("TCR_clone_count")+
  ylim(0,200)+
  theme(plot.title = element_text(hjust = 0.5,size = 14),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))+
  scale_discrete_manual(values=c("#4DBBD5FF","#E64B35FF"),aesthetics = 'colour')+NoLegend()
p
ggsave(p,file="/Users/wangjun/Desktop/TCR_clone/total_new.png",dpi = 1000,width = 2.12,height = 3.44)





library(Seurat)
#读取大类文件
pbmc<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/pbmc_merge_new.rds")
table(Idents(pbmc))





pbmc=FindVariableFeatures(object = pbmc)
f=VariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc, features = f)
pbmc <- RunPCA(pbmc, features = f)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.06)


pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.06)



table(Idents(pbmc))

pbmc<-AddMetaData(pbmc,Idents(pbmc),col.name = "Annotation")


Idents(pbmc)<-"orig.ident"

pbmc<-subset(pbmc,idents="T3",invert=T)

#UMAP
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Annotation")
p3 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)#, cols =colorlist)
p1

p6




#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:30)
p5 <- DimPlot(pbmc, reduction = "tsne")
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5


pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:25)
p5 <- DimPlot(pbmc, reduction = "tsne")
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5



#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20)
p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5



#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:18)
p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5







#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:15)
p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5


#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:10)
p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5

for (i in 10:30){
  pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:i)
  p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
  #p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
  ggsave(p5,file=paste0("/Users/wangjun/Desktop/",i,"_new.png"))
}


#TNSE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:5)
p5 <- DimPlot(pbmc, reduction = "tsne")#,cols = pal)
#p6 <- DimPlot(pbmc, reduction = "tsne", repel = TRUE)
p5




prop.table(table(pbmc$Annotation))


table(Idents(pbmc))

Idents(pbmc)<-"Annotation"

Idents(pbmc)<-"orig.ident"

pbmc<-subset(pbmc,idents="T3",invert=T)

Idents(pbmc)<-"new"
table(Idents(pbmc))

#鉴定出相应的细胞类群后
new.cluster.ids <- c("T cells","Epithelial cells","Myeloid cells","NK cells","Fibroblasts","Endothelial cells","B cells","Mast cells","Plasma cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

table(pbmc$)





Myeloid_cells<-pbmc
df = data.frame(clu=names(table(Idents(Myeloid_cells))),
                per=sprintf("%1.2f%%", 100*table((Idents(Myeloid_cells)))/length(Idents(Myeloid_cells))))

Myeloid_cells$per = df[match(Idents(Myeloid_cells),df$clu),2]
Myeloid_cells$new = paste0(Idents(Myeloid_cells),"(",Myeloid_cells$per,")")
table(Myeloid_cells$new)


table(pbmc$orig.ident)



Myeloid_cells$new<-factor(Myeloid_cells$new,levels = c("T cells(38.92%)","Epithelial cells(22.19%)","Myeloid cells(11.60%)","NK cells(6.92%)","Fibroblasts(6.46%)","Endothelial cells(6.39%)","B cells(4.33%)","Mast cells(2.02%)","Plasma cells(1.16%)"))

Myeloid_cells$new<-factor(Myeloid_cells$new,levels = c("T cells(37.47%)","Epithelial cells(22.14%)","Myeloid cells(12.18%)","NK cells(6.99%)","Fibroblasts(6.58%)","Endothelial cells(6.35%)","B cells(5.06%)","Mast cells(2.02%)","Plasma cells(1.21%)"))

Myeloid_cells$new<-factor(Myeloid_cells$new,levels = c("T cells(38.43%)","Epithelial cells(22.37%)","Myeloid cells(11.69%)","NK cells(6.98%)","Fibroblasts(6.51%)","Endothelial cells(6.45%)","B cells(4.36%)","Mast cells(2.04%)","Plasma cells(1.17%)"))


Idents(Myeloid_cells)<-"new"







#TNSE
Myeloid_cells <- RunTSNE(Myeloid_cells, reduction = "pca", dims = 1:16)
p5 <- DimPlot(Myeloid_cells,reduction = "tsne",cols = pal[c(1,10,6,5,7,2,3,4,9)])
p5




pbmc<-Myeloid_cells
Idents(pbmc)<-"new"
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:25)
p5 <- DimPlot(Myeloid_cells,reduction = "tsne",cols = pal[c(1,10,6,5,7,2,3,4,9)])
p5

table(pbmc$orig.ident)



ggsave(p5,file="/Users/wangjun/Desktop/Total_new.png",dpi = 1000,width = 7.76,height = 5.44)

ggsave(p5,file="/Users/wangjun/Desktop/Total_new.pdf",dpi = 1000,width = 7.76,height = 5.44)

table(pbmc$new)

Idents(pbmc)<-"Annotation"

markers<-c("CD3D","CD3E","NKG7","CD79A","IGKC","LYZ","CD68","TPSB2","CPA3","MS4A2","EPCAM","DCN","C1R","COL1A1","PECAM1","CLDN5")

p1<-DotPlot(pbmc, features = markers, cols = c("#4DBBD5FF","#E64B35FF")) + RotatedAxis()

table(Idents(pbmc))
Idents(pbmc)<-factor(Idents(pbmc),levels = c("T cells","NK cells","B cells","Plasma cells","Myeloid cells","Mast cells","Epithelial cells","Fibroblasts","Endothelial cells"))

library(ggplot2)
ggsave(p1,file="/Users/wangjun/Desktop/更新结果/dotplot.pdf",width = 9.08,height = 3.76)



for (i in 1:length(markers)){
  p<-FeaturePlot(pbmc,reduction = "tsne", features = markers[i],cols = c("grey","red"),pt.size = 0.5, repel = TRUE,order = T)+ RotatedAxis()
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/",markers[i],"_no_label_new.pdf"),width = 6.22,height = 5.93)
  
}


for (i in 1:length(markers)){
  p<-FeaturePlot(pbmc,reduction = "tsne", features = markers[i],cols = c("grey","red"),label = TRUE,label.size =3,pt.size = 0.5, repel = TRUE,order = T)+ RotatedAxis()
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/",markers[i],"_label_new.pdf"),width = 6.22,height = 5.93)
  
}



p<-FeaturePlot(pbmc,reduction = "tsne", features = "CCR2",cols = c("grey","red"),pt.size = 1,label = TRUE,label.size =3, repel = TRUE,order = T)+ RotatedAxis()

p<-FeaturePlot(pbmc,reduction = "tsne", features = "CCR2",cols = c("grey","red"),pt.size = 0.5, repel = TRUE,order = T)+ RotatedAxis()

p








pal

col<-pal[c(1,10,6,5,7,2,3,4,9)]
col[4]<-"#FFE79F"
col[4]<-"#ffe79f"




library(paletteer) 
library(scales)
pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal[c(2,3,7,10,5,8,4,1,9)]

#pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,5,8,7,1,9)] 
pal <- paletteer_d("ggsci::nrc_npg")[c(2,3,7,10,5,8,4,1,9)] 





table(Idents(Myeloid_cells))


Idents(pbmc)<-"orig.ident"

#鉴定出相应的细胞类群后
new.cluster.ids <- c("LUAD","MIA","MIA","LUAD","LUAD","MIA","LUAD","MIA","LUAD","MIA","LUAD","LUAD","LUAD","MIA","LUAD","LUAD")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)





pbmc<-AddMetaData(pbmc,Idents(pbmc),col.name="Type")

table(Idents(pbmc))


Idents(pbmc)<-"orig.ident"

p5<-DimPlot(pbmc,group.by = "orig.ident",reduction = "tsne")


p5<-DimPlot(pbmc,split.by = "Type",reduction = "tsne")


p5<-DimPlot(Myeloid_cells,split.by = "Type",reduction = "tsne")


p5<-DimPlot(pbmc,group.by = "orig.ident",reduction = "tsne")

p5<-DimPlot(pbmc,group.by = "Type",reduction = "tsne")


ggsave(p5,file="/Users/wangjun/Desktop/更新结果/sample_orig.ident.pdf",width = 6.22,height = 5.93)
ggsave(p5,file="/Users/wangjun/Desktop/更新结果/sample_type.pdf",width = 6.22,height = 5.93)
ggsave(p5,file="/Users/wangjun/Desktop/更新结果/sample_type_split.pdf",width = 11,height = 5.93)


cols = pal[c(1,10,6,5,7,2,3,4,9)]









#读取基因数据
sample_list=c("FD1_LC2","FD1_LC3","FD2_LC1","FD4_LC1","FD4_LC2","FD5_LC1","FD5_LC2","FD8_LC1","FD8_LC2","FD9_LC1","FD9_LC2","FD14_LC1","FD14_LC2","FD16_LC1","FD16_LC2")
sample_T_list=c("T4","T5","T1","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")



data<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_snv_index_result.csv",sep = ",",header = T)
data2<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_indel_index_result.csv",sep = ",",header = T)


data_total<-rbind(data[,c(1,3)],data2[,c(1,3)])



#读取T细胞数据
library(Seurat)
Tumor<-readRDS(file="/Users/wangjun/Desktop/最终方案/figure/tumor_infercnv.rds")
Idents(Tumor)<-"orig.ident"



unique(data2[data2$SampleID%in%"FD2_LC2",]$Gene)
unique(data2[data2$SampleID%in%"FD16_LC2",]$Gene)

names(Tumor_set@meta.data)


score_list<-list()

for (i in 1:length(sample_list)){
  #for (i in 1:1){
  print(i)
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=unique(data[data$SampleID%in%sample,]$Gene)
  print(gene_list)
  Tumor_set<-subset(Tumor,idents = sample_T)
  
  Tumor_set <- AddModuleScore(object=Tumor_set, features=gene_list, name="to_compare")
  
  score_list[[i]]<-rowMeans(Tumor_set@meta.data[,13:(12+length(gene_list))])
  }


score_list2<-list()
for (i in 1:length(sample_list)){
  #for (i in 1:1){
  print(i)
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=unique(data2[data2$SampleID%in%sample,]$Gene)
  print(gene_list)
  if (length(gene_list)!=0){
    Tumor_set<-subset(Tumor,idents = sample_T)
    
    Tumor_set <- AddModuleScore(object=Tumor_set, features=gene_list, name="to_compare")
    
    if (length(gene_list)!=1){
      score_list2[[i]]<-rowMeans(Tumor_set@meta.data[,13:(12+length(gene_list))])
    }
    else{
      score_list2[[i]]<-Tumor_set@meta.data[,13:(12+length(gene_list))]
    }
  }
  if (length(gene_list)==0){
    score_list2[[i]]<-NULL
  }
}

score_list2[[15]]<-NULL


score_list2[[14]]


length(score_list2)

table(Idents(Tumor))


MIA_score<-c(score_list[[1]],score_list[[5]],score_list[[7]],score_list[[11]],score_list[[13]],
             score_list2[[1]],score_list2[[5]],score_list2[[7]],score_list2[[11]],score_list2[[13]])

LUAD_score<-c(score_list[[2]],score_list[[3]],score_list[[4]],score_list[[6]],score_list[[8]],score_list[[9]],score_list[[10]],score_list[[12]],score_list[[14]],score_list[[15]],
              score_list2[[2]],score_list2[[3]],score_list2[[4]],score_list2[[6]],score_list2[[8]],score_list2[[9]],score_list2[[10]],score_list2[[12]],score_list2[[14]],score_list2[[15]])





score_total_list<-list()

for (i in 1:length(sample_list)){
  #for (i in 1:1){
  print(i)
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  gene_list=unique(data_total[data_total$SampleID%in%sample,]$Gene)
  print(gene_list)
  Tumor_set<-subset(Tumor,idents = sample_T)
  
  Tumor_set <- AddModuleScore(object=Tumor_set, features=gene_list, name="to_compare")
  
  score_total_list[[i]]<-rowMeans(Tumor_set@meta.data[,13:(12+length(gene_list))])
}

names(Tumor_set@meta.data)

for (i in 1:15){
  s=mean(score_total_list[[i]])
  print(paste0(sample_list[i],":",s))
  
}











#绘制一下小提琴图

wilcox.test(MIA_score,LUAD_score)


length(MIA_score)+length(LUAD_score)



data<-data.frame(matrix(NA,22340,2))
colnames(data)<-c("Neoantigen_Score","Type")

data[1:6200,]$Neoantigen_Score<-MIA_score
data[1:6200,]$Type<-"MIA"
data[6201:22340,]$Neoantigen_Score<-LUAD_score
data[6201:22340,]$Type<-"LUAD"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))



mean(MIA_score)
mean(LUAD_score)



library(ggplot2)
p<-ggplot(data, aes(x=Type, y=Neoantigen_Score, fill=Type)) + 
  geom_violin(trim = FALSE, size = 0.7, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5,size=0.01) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(title = "Violin_Plot", x = "Tissue", y = "Neoantigen_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

p

t.test(MIA_score,LUAD_score,alternative = "greater")

library(lsr)
d=cohensD(MIA_score,LUAD_score)


p<-ggplot(data, aes(x=Type, y=Neoantigen_Score, fill=Type)) + 
  ggtitle(paste0("P = 6.902e-08","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_Neoantigen_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p


ggsave(p,file="/Users/wangjun/Desktop/violin_plot_Neoantigen_Score_new.pdf",height = 6.37,width = 5.14)





t.test(MIA_score,LUAD_score)

library(lsr)
cohensD(MIA_score,LUAD_score)


ggsave(p,file="/Users/wangjun/Desktop/volin_plot.png",dpi = 1000,height = 6.37,width = 5.14)




for (i in 1:15){
  s=mean(c(score_list[[i]],score_list2[[i]]))
  print(paste0(sample_list[i],":",s))
  
}


mean(c(score_list[[3]],score_list2[[1]]))

length(score_list)

table(Idents(Tumor))


#计算每个样本的新抗原相关基因的表达量检验
MIA_neo_exp<-c(0.00133980243522386,-0.00778282057081626,-0.00900144442922002,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.0119153716769999,0.104798733518385,-0.0140605621295198,-0.0522850686551056,-0.0433967387739034,-0.00213491823983392,-0.047735717398372,-0.049092356836148,0.03447649555235,0.0014202156078874)

wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater")


MIA_neo_exp<-c(0.00133980243522386,-0.00778282057081626,-0.00900144442922002,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.0119153716769999,-0.0140605621295198,-0.0522850686551056,-0.0433967387739034,-0.049092356836148)
wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater",paired = T)



MIA_neo_exp<-c(-0.012113429326935,-0.0167801548100993,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.00361242206404616,-0.0219726189882752,-0.097666121143288,-0.0897094438146464,-0.0960410431710384)
wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater",paired = T)


MIA_neo_exp<-c(-0.012113429326935,-0.0167801548100993,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.00361242206404616,-0.0219726189882752,-0.097666121143288,-0.0897094438146464,-0.0960410431710384,0.148025326015436,-0.0433967387739034,-0.00363452727604194,0.00509901536902821,0.00142336464471198)
wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater")





length(MIA_neo_exp)+length(LUAD_neo_exp)


data<-data.frame(matrix(NA,10,2))
colnames(data)<-c("Exhausted_Score","Type")

data[1:8120,]$Exhausted_Score<-MIA_score
data[1:8120,]$Type<-"MIA"
data[8121:25609,]$Exhausted_Score<-LUAD_score
data[8121:25609,]$Type<-"LUAD"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))





p<-ggplot(data, aes(x=Type, y=Exhausted_Score, fill=Type)) + 
  geom_violin(trim = FALSE, size = 0.7, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5,size=0.01) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(title = "Violin_Plot", x = "Tissue", y = "Exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

p


p<-ggplot(data, aes(x=Type, y=Exhausted_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0003","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_MANA_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p





#获取含有新抗原的CD8+T 细胞的barcode
sample_list=c("FD1_LC2","FD1_LC3","FD2_LC1","FD4_LC1","FD4_LC2","FD5_LC1","FD5_LC2","FD8_LC1","FD8_LC2","FD9_LC1","FD9_LC2","FD14_LC1","FD14_LC2","FD16_LC1","FD16_LC2")
sample_T_list=c("T4","T5","T1","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")

library(Seurat)
T_total<-readRDS(file = "/Users/wangjun/Desktop/T_total_3_26.rds")

Idents(T_total)<-"orig.ident"

table(Idents(T_total))


anno<-data.frame(T_total$orig.ident)

anno$barcode<-rownames(anno)
anno[1:3,]
dim(anno)

my_list <- strsplit(rownames(anno), "_")


for (i in 1:length(rownames(anno))){
  print(i)
  if (nchar(anno$barcode[i])==18){
    anno$barcode[i]<-anno$barcode[i]
  }
  else{
    anno$barcode[i]<-my_list[[i]][3]
  }
}

length(anno$barcode)
dim(T_total)



length(unique(anno$barcode))
length(anno$barcode)

barcode<-unique(anno$barcode)


data<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_snv_index_result_new_best.csv",sep = ",",header = T)

data[1:4,]

barcode_list<-list()
for (i in 1:length(sample_list)){
  #for (i in 1:1){
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  print(i)
  #gene_list=c("LAG3","TIGIT","PDCD1","CTLA4","KLRG1")
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  col_n<-in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode
  T_total_set<-subset(T_total,idents = sample_T)
  print(length(rownames(anno[anno$barcode%in%col_n & anno$T_total.orig.ident%in%sample_T,])))
  print(rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,]))
  barcode_list[[i]]<-rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,])
  }


final_convert<-list()
for (i in 1:length(sample_list)){
  sample=sample_list[i]
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  anno_barcode<-anno[rownames(anno)%in%barcode_list[[i]],]
  anno_barcode$barcode_sample<-rownames(anno_barcode)
  TCR_barcode<-in_put3[in_put3$barcode%in%anno_barcode$barcode,]
  best_TCR<-data[data$SampleID%in%sample,c(1,2,3,4,16)]
  print(length(rownames(anno_barcode)))
  print(length(rownames(TCR_barcode)))
  print(length(rownames(best_TCR)))
  merge_1<-merge(TCR_barcode,anno_barcode,by.x = "barcode",by.y = "barcode")
  #print(paste0("merge_1:",length(rownames(merge_1))))
  merge_2<-merge(merge_1,best_TCR,by.x = "TRB_CDR3",by.y = "BestMatchTCR")
  #print(paste0("merge_2:",length(rownames(merge_2))))
  final_csv<-merge_2[,c("SampleID","HLA","Gene","Base_site","TRA_CDR3","TRB_CDR3","CDR3","Frequency","barcode","barcode_sample","T_total.orig.ident")]
  #print(paste0("final_csv:",length(rownames(final_csv))))
  final_convert[[i]]<-final_csv
}

final_frame<-data.frame()
for (i in 1:length(sample_list)){
  final_frame<-rbind(final_frame,final_convert[[i]])
  #print(dim(final_convert[[i]]))
}
write.table(final_frame,file="/Users/wangjun/Desktop/barcode_TCR_snv_mutation.csv",col.names = T,row.names = F,sep = ",")








#Indel
data<-read.table(file="/Users/wangjun/Desktop/matrix_for_panpep_indel_index_result_new_best.csv",sep = ",",header = T)
colnames(data)

#存下每个细胞的耗竭得分
barcode_list2<-list()
for (i in 1:length(sample_list)){
  #for (i in 1:1){
  sample=sample_list[i]
  sample_T=sample_T_list[i]
  print(i)
  #gene_list=c("LAG3","TIGIT","PDCD1","CTLA4","KLRG1")
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  col_n<-in_put3[in_put3$TRB_CDR3%in%data[data$SampleID%in%sample,]$BestMatchTCR,]$barcode
  T_total_set<-subset(T_total,idents = sample_T)
  print(length(rownames(anno[anno$barcode%in%col_n & anno$T_total.orig.ident%in%sample_T,])))
  barcode_list2[[i]]<-rownames(anno[anno$barcode%in%col_n &anno$T_total.orig.ident%in%sample_T,])
}


final_convert<-list()
for (i in 1:length(sample_list)){
  sample=sample_list[i]
  in_put3<-read.table(file=paste0("/Users/wangjun/Desktop/TCR_update/",sample,"/Clonotype_Frequency_barcode.csv"),sep = ",",header = T)
  anno_barcode<-anno[rownames(anno)%in%barcode_list2[[i]],]
  anno_barcode$barcode_sample<-rownames(anno_barcode)
  TCR_barcode<-in_put3[in_put3$barcode%in%anno_barcode$barcode,]
  best_TCR<-data[data$SampleID%in%sample,c(1,2,3,4,13)]
  print(length(rownames(anno_barcode)))
  print(length(rownames(TCR_barcode)))
  print(length(rownames(best_TCR)))
  merge_1<-merge(TCR_barcode,anno_barcode,by.x = "barcode",by.y = "barcode")
  #print(paste0("merge_1:",length(rownames(merge_1))))
  merge_2<-merge(merge_1,best_TCR,by.x = "TRB_CDR3",by.y = "BestMatchTCR")
  #print(paste0("merge_2:",length(rownames(merge_2))))
  final_csv<-merge_2[,c("SampleID","HLA","Gene","Base_site","TRA_CDR3","TRB_CDR3","CDR3","Frequency","barcode","barcode_sample","T_total.orig.ident")]
  #print(paste0("final_csv:",length(rownames(final_csv))))
  final_convert[[i]]<-final_csv
}

data

final_frame<-data.frame()
for (i in 1:length(sample_list)){
  final_frame<-rbind(final_frame,final_convert[[i]])
  #print(dim(final_convert[[i]]))
}
write.table(final_frame,file="/Users/wangjun/Desktop/barcode_TCR_indel_mutation.csv",col.names = T,row.names = F,sep = ",")







MIA_barcode<-c(barcode_list[[1]],barcode_list[[5]],barcode_list[[7]],barcode_list[[11]],barcode_list[[13]],
             barcode_list2[[1]],barcode_list2[[5]],barcode_list2[[7]],barcode_list2[[11]],barcode_list2[[13]])

LUAD_barcode<-c(barcode_list[[2]],barcode_list[[3]],barcode_list[[4]],barcode_list[[6]],barcode_list[[8]],barcode_list[[9]],barcode_list[[10]],barcode_list[[12]],barcode_list[[14]],barcode_list[[15]],
              barcode_list2[[2]],barcode_list2[[3]],barcode_list2[[4]],barcode_list2[[6]],barcode_list2[[8]],barcode_list2[[9]],barcode_list2[[10]],barcode_list2[[12]],barcode_list2[[14]],barcode_list2[[15]])

LUAD_barcode












T_total_copy<-T_total


names(T_total@meta.data)

Idents(T_total)<-"Cluster"

table(Idents(T_total))


T_total<-subset(T_total,idents = c("CD8_C1","CD8_C2","CD8_C3","CD8_C4","CD8_C5","CD8_C6","CD8_C7"))

Idents(T_total)<-"orig.ident"

table(Idents(T_total))


new.cluster.ids=c("LUAD","MIA","MIA","LUAD","LUAD","MIA","LUAD","MIA","LUAD","MIA","LUAD","LUAD","LUAD","MIA","LUAD","LUAD")

names(new.cluster.ids) <- levels(T_total)
T_total <- RenameIdents(T_total, new.cluster.ids)

T_total<-AddMetaData(T_total,Idents(T_total),col.name = "Type")








MANA_score<-c("CXCL13","HLA-DRA","HLA-DRB5","HLA-DQA1","HLA-DRB1","HLA-DQB1","HLA-DPA1","HLA-DPB1")


T_total <- AddModuleScore(object=T_total, features=MANA_score, name="MANA")


#score_list2[[i]]<-rowMeans(Tumor_set@meta.data[,13:(12+length(gene_list))])
names(T_total@meta.data)

T_total@meta.data[,24:(23+length(MANA_score))][1:4,]

MANA_score_exp<-rowMeans(T_total@meta.data[,24:(23+length(MANA_score))])

MIA_MANA_score_exp<-rowMeans(subset(T_total,idents = "MIA")@meta.data[,24:(23+length(MANA_score))])
LUAD_MANA_score_exp<-rowMeans(subset(T_total,idents = "LUAD")@meta.data[,24:(23+length(MANA_score))])

length(MIA_MANA_score_exp)


barcode_list<-c(MIA_barcode,LUAD_barcode)




t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

library(lsr)
cohensD(MANA_score_exp[barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])




p=t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!(names(MANA_score_exp)%in%barcode_list)])


length(MANA_score_exp[names(MANA_score_exp)%in%barcode])


#绘图
data<-data.frame(matrix(NA,7590,2))
colnames(data)<-c("MANA_Score","Type")

data[1:339,]$MANA_Score<-MANA_score_exp[names(MANA_score_exp)%in%barcode_list]
data[1:339,]$Type<-"TCR"
data[340:7590,]$MANA_Score<-MANA_score_exp[!names(MANA_score_exp)%in%barcode_list]
data[340:7590,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0003","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_MANA_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MANA_Score_new.pdf",height = 6.37,width = 5.14)









t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])+length(MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


#绘图
data<-data.frame(matrix(NA,3597,2))
colnames(data)<-c("MIA_MANA_Score","Type")

data[1:117,]$MIA_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[1:117,]$Type<-"MIA_TCR"
data[118:3597,]$MIA_MANA_Score<-MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode]
data[118:3597,]$Type<-"MIA_nonTCR"

data$Type<-factor(data$Type,levels = c("MIA_nonTCR","MIA_TCR"))


p=t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")
d=cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])
d
p<-ggplot(data, aes(x=Type, y=MIA_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = ",round(p$p.value,4),"   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs( x = "Tissue", y = "MIA_MANA_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MIA_MANA_Score_new.pdf",height = 6.37,width = 5.14)





LUAD_MANA_score_exp

t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")

library(lsr)
cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,3993,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD_TCR"
data[223:3993,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[223:3993,]$Type<-"LUAD_nonTCR"

data$Type<-factor(data$Type,levels = c("LUAD_nonTCR","LUAD_TCR"))

data[1:4,]

p=t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")
d=cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = ",round(p$p.value,4),"   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "LUAD_MANA_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_LUAD_MANA_Score_new.pdf",height = 6.37,width = 5.14)





t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")


library(lsr)
cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])



t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(unique(names(MANA_score_exp)))

p=t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "less")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])
length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,339,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD"
data[223:339,]$LUAD_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[223:339,]$Type<-"MIA"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.2898","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "TCR_MANA_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MANA_TCR_Score_new.pdf",height = 6.37,width = 5.14)








#换一种MANA的基因集

MANA_score<-c("CXCL13", "HLA-DRA", "HLA-DRB5", "HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "CCL3", "GZMA", "GEM", "ENTPD1", "HLA-DPA1", "TNS3", "MIR4435-2HG", "IFNG", "HLA-DPB1")

GetAssayData(T_total_set)[MANA_score,1:4]


T_total <- AddModuleScore(object=T_total, features=MANA_score, name="MANA")
names(T_total@meta.data)



MIA_MANA_score_exp<-rowMeans(subset(T_total,idents = "MIA")@meta.data[,24:(23+length(MANA_score))])
LUAD_MANA_score_exp<-rowMeans(subset(T_total,idents = "LUAD")@meta.data[,24:(23+length(MANA_score))])


MANA_score_exp<-rowMeans(T_total@meta.data[,24:(23+length(MANA_score))])

barcode_list<-c(MIA_barcode,LUAD_barcode)


t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])


length(MANA_score_exp[names(MANA_score_exp)%in%barcode])

p=t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])
p$p.value
#绘图
data<-data.frame(matrix(NA,7590,2))
colnames(data)<-c("MANA_Score","Type")

data[1:339,]$MANA_Score<-MANA_score_exp[names(MANA_score_exp)%in%barcode_list]
data[1:339,]$Type<-"TCR"
data[340:7590,]$MANA_Score<-MANA_score_exp[!names(MANA_score_exp)%in%barcode_list]
data[340:7590,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))



p<-ggplot(data, aes(x=Type, y=MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 5.336e-06","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_MANA_Score_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MANA_Score_2_new.pdf",height = 6.37,width = 5.14)








#

t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


length(MIA_MANA_score_exp)
length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])

p=t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

d=cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])

#绘图
data<-data.frame(matrix(NA,3597,2))
colnames(data)<-c("MIA_MANA_Score","Type")

data[1:117,]$MIA_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[1:117,]$Type<-"TCR"
data[118:3597,]$MIA_MANA_Score<-MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode]
data[118:3597,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))



p<-ggplot(data, aes(x=Type, y=MIA_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0003","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "MIA_MANA_Score_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MIA_MANA_Score_2_new.pdf",height = 6.37,width = 5.14)











t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")

library(lsr)
cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])




p=t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")

d=cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,3993,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"TCR"
data[223:3993,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[223:3993,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))

p

p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0011","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "LUAD_MANA_Score_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_LUAD_MANA_Score_2_new.pdf",height = 6.37,width = 5.14)












#比较同样有新抗原TCR的CD8+T 细胞，MIA和LUAD那组高


t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(unique(names(MANA_score_exp)))

p=t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "less")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])
length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,339,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD"
data[223:339,]$LUAD_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[223:339,]$Type<-"MIA"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.3273","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "TCR_MANA_Score_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MANA_TCR_Score_2_new.pdf",height = 6.37,width = 5.14)












#查看virus score的表达情况，这里MANA没有修改

#换一种MANA的基因集

MANA_score<-c("YPEL5", "CLDND1", "CCND3", "KLF3", "NUDT11", "LYAR", "PIK3R1", "RARA", "AUTS2", "PTGER2", "ARL4A", "GPR15", "KLRC1", "GPR183", "IL7R")

GetAssayData(T_total_set)[MANA_score,1:4]


T_total <- AddModuleScore(object=T_total, features=MANA_score, name="MANA")
names(T_total@meta.data)



MIA_MANA_score_exp<-rowMeans(subset(T_total,idents = "MIA")@meta.data[,24:(23+length(MANA_score))])
LUAD_MANA_score_exp<-rowMeans(subset(T_total,idents = "LUAD")@meta.data[,24:(23+length(MANA_score))])


MANA_score_exp<-rowMeans(T_total@meta.data[,24:(23+length(MANA_score))])

barcode_list<-c(MIA_barcode,LUAD_barcode)


t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "less")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])


length(MANA_score_exp[names(MANA_score_exp)%in%barcode])



p=t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "less")

d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])

#绘图
data<-data.frame(matrix(NA,7590,2))
colnames(data)<-c("MANA_Score","Type")

data[1:339,]$MANA_Score<-MANA_score_exp[names(MANA_score_exp)%in%barcode_list]
data[1:339,]$Type<-"TCR"
data[340:7590,]$MANA_Score<-MANA_score_exp[!names(MANA_score_exp)%in%barcode_list]
data[340:7590,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))

p

p<-ggplot(data, aes(x=Type, y=MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0109","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_virus_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_virus_Score_new.pdf",height = 6.37,width = 5.14)








#

t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


length(MIA_MANA_score_exp)
length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])




p=t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

d=cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


#绘图
data<-data.frame(matrix(NA,3597,2))
colnames(data)<-c("MIA_MANA_Score","Type")

data[1:117,]$MIA_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[1:117,]$Type<-"TCR"
data[118:3597,]$MIA_MANA_Score<-MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode]
data[118:3597,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))



p<-ggplot(data, aes(x=Type, y=MIA_MANA_Score, fill=Type)) +
  ggtitle(paste0("P = 0.2908","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "MIA_virus_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MIA_virus_Score_new.pdf",height = 6.37,width = 5.14)











t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "less")

library(lsr)
cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])


p=t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "less")

d=cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])



#绘图
data<-data.frame(matrix(NA,3993,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"TCR"
data[223:3993,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[223:3993,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))



p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0006","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "LUAD_virus_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_LUAD_virus_Score_new.pdf",height = 6.37,width = 5.14)




#比较同样有新抗原TCR的CD8+T 细胞，MIA和LUAD那组高


t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "less")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[!names(MANA_score_exp)%in%MIA_barcode])


#断

p=t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "less")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])
length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,339,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD"
data[223:339,]$LUAD_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[223:339,]$Type<-"MIA"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.01644","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "TCR_virus_Score_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_virus_TCR_Score_new.pdf",height = 6.37,width = 5.14)








#正经的耗竭得分来咯
MANA_score<-c("HAVCR2","CXCL13","CCL3","SIRPG","IFNG","TIGIT","GZMB","PDCD1","PARK7","TNFRSF9","ACP5","CTLA4","RBPJ","CXCR6","CD27","FKBP1A","BST2","TPI1","MIR155HG","PTTG1","CD63","SAMSN1","RGS1","ITGAE","HLA-DRA","IGFLR1","KRT86","ENTPD1","DUSP4","SIT1","TOX","PHLDA1","CCND2","GPR25","LAYN","PRDX5","SARDH","FASLG","ANXA5","CTSD","PDIA6","RANBP1","FKBP1A","SDCBP2","COTL1","TNFRSF1B","IDH2","CD38","CD82","LAG3","MIR497HG","APOBEC3C","ITM2A","COX5A","IFI35","NDFIP2","TNFRSF18","KRT81","DNPH1","RGS2","HMGN1","DYNLL1","SNRPB","SYNGR2","RAB27A","PSMC3","GALM","FABP5","UBE2L6","MYO7A","PRDX3","DDIT4","STMN1","CDK2AP2","VCAM1","SNAP47","PSMB3","ISG15","HLA-DRB5","CKS2","TNIP3","CD7","PSMD4","ATP6V1C2","PSMD8")

GetAssayData(T_total_set)[MANA_score,1:4]


T_total <- AddModuleScore(object=T_total, features=MANA_score, name="MANA")
names(T_total@meta.data)



MIA_MANA_score_exp<-rowMeans(subset(T_total,idents = "MIA")@meta.data[,c(24:108)])
LUAD_MANA_score_exp<-rowMeans(subset(T_total,idents = "LUAD")@meta.data[,c(24:108)])


MANA_score_exp<-rowMeans(T_total@meta.data[,c(24:108)])

barcode_list<-c(MIA_barcode,LUAD_barcode)


t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

library(lsr)
cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])


length(MANA_score_exp[names(MANA_score_exp)%in%barcode])


p=t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])


#绘图
data<-data.frame(matrix(NA,7590,2))
colnames(data)<-c("MANA_Score","Type")

data[1:339,]$MANA_Score<-MANA_score_exp[names(MANA_score_exp)%in%barcode_list]
data[1:339,]$Type<-"TCR"
data[340:7590,]$MANA_Score<-MANA_score_exp[!names(MANA_score_exp)%in%barcode_list]
data[340:7590,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 4.509e-13","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_exhausted_Score_new.pdf",height = 6.37,width = 5.14)








#

t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


length(MIA_MANA_score_exp)
length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])


p=t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

d=cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])

#绘图
data<-data.frame(matrix(NA,3597,2))
colnames(data)<-c("MIA_MANA_Score","Type")

data[1:117,]$MIA_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[1:117,]$Type<-"TCR"
data[118:3597,]$MIA_MANA_Score<-MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode]
data[118:3597,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))



p<-ggplot(data, aes(x=Type, y=MIA_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.0003","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "MIA_exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MIA_exhausted_Score_new.pdf",height = 6.37,width = 5.14)











p=t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")

library(lsr)
d=cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,3993,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"TCR"
data[223:3993,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[223:3993,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 3.656e-10","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "LUAD_exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_LUAD_exhasted_Score_new.pdf",height = 6.37,width = 5.14)




#比较同样有新抗原TCR的CD8+T 细胞，MIA和LUAD那组高


p=t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])
length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,339,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD"
data[223:339,]$LUAD_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[223:339,]$Type<-"MIA"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))

d
p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 3.194e-06","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "TCR_Exhausted_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_exhasted_TCR_Score_new.pdf",height = 6.37,width = 5.14)












#免疫检查点得分
MANA_score<-c("CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1")
GetAssayData(T_total_set)[MANA_score,1:4]


T_total <- AddModuleScore(object=T_total, features=MANA_score, name="MANA")
names(T_total@meta.data)


MIA_MANA_score_exp<-rowMeans(subset(T_total,idents = "MIA")@meta.data[,24:(23+length(MANA_score))])
LUAD_MANA_score_exp<-rowMeans(subset(T_total,idents = "LUAD")@meta.data[,24:(23+length(MANA_score))])


MANA_score_exp<-rowMeans(T_total@meta.data[,24:(23+length(MANA_score))])

barcode_list<-c(MIA_barcode,LUAD_barcode)


p=t.test(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list],alternative = "greater")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%barcode_list],MANA_score_exp[!names(MANA_score_exp)%in%barcode_list])


length(MANA_score_exp[names(MANA_score_exp)%in%barcode])


#绘图
data<-data.frame(matrix(NA,7590,2))
colnames(data)<-c("MANA_Score","Type")

data[1:339,]$MANA_Score<-MANA_score_exp[names(MANA_score_exp)%in%barcode_list]
data[1:339,]$Type<-"TCR"
data[340:7590,]$MANA_Score<-MANA_score_exp[!names(MANA_score_exp)%in%barcode_list]
data[340:7590,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 1.791e-07","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "Total_immune_checkpoint_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_immune_checkpoint_Score_new.pdf",height = 6.37,width = 5.14)








#

p=t.test(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
d=cohensD(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode],MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode])


length(MIA_MANA_score_exp)
length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])


#绘图
data<-data.frame(matrix(NA,3597,2))
colnames(data)<-c("MIA_MANA_Score","Type")

data[1:117,]$MIA_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[1:117,]$Type<-"TCR"
data[118:3597,]$MIA_MANA_Score<-MIA_MANA_score_exp[!names(MIA_MANA_score_exp)%in%MIA_barcode]
data[118:3597,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=MIA_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 0.1967","   Cohen’s d = ",round(d,2)))+
    geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "MIA_immune_checkpoint_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_MIA_immune_checkpoint_Score_new.pdf",height = 6.37,width = 5.14)











p=t.test(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode],alternative = "greater")

library(lsr)
d=cohensD(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode],LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode])


length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,3993,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"TCR"
data[223:3993,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[!names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[223:3993,]$Type<-"nonTCR"

data$Type<-factor(data$Type,levels = c("nonTCR","TCR"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 1.054e-06","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "LUAD_immune_checkpoint_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_LUAD_immune_checkpoint_Score_new.pdf",height = 6.37,width = 5.14)




#比较同样有新抗原TCR的CD8+T 细胞，MIA和LUAD那组高


p=t.test(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode],alternative = "greater")

library(lsr)
d=cohensD(MANA_score_exp[names(MANA_score_exp)%in%LUAD_barcode],MANA_score_exp[names(MANA_score_exp)%in%MIA_barcode])

length(MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode])
length(LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode])

#绘图
data<-data.frame(matrix(NA,339,2))
colnames(data)<-c("LUAD_MANA_Score","Type")

data[1:222,]$LUAD_MANA_Score<-LUAD_MANA_score_exp[names(LUAD_MANA_score_exp)%in%LUAD_barcode]
data[1:222,]$Type<-"LUAD"
data[223:339,]$LUAD_MANA_Score<-MIA_MANA_score_exp[names(MIA_MANA_score_exp)%in%MIA_barcode]
data[223:339,]$Type<-"MIA"

data$Type<-factor(data$Type,levels = c("MIA","LUAD"))


p
p<-ggplot(data, aes(x=Type, y=LUAD_MANA_Score, fill=Type)) + 
  ggtitle(paste0("P = 2.639e-08","   Cohen’s d = ",round(d,2)))+
  geom_violin(trim = FALSE, size = 1.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(x = "Tissue", y = "immune_checkpoint_Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
p
ggsave(p,file="/Users/wangjun/Desktop/更新结果/violin_plot_immune_checkpoint_TCR_Score_new.pdf",height = 6.37,width = 5.14)


table(pbmc$orig.ident)


#重新检查合并情况
data<-readRDS("/Users/wangjun/Desktop/最终方案/figure/T_total_update.rds")

table(data$orig.ident)

setwd("/Users/wangjun/Desktop/最终方案/figure/")

#载入EPCAM+
epcam<-readRDS(file="./epcam.rds")
#载入T细胞
T_cells<-readRDS(file="/Users/wangjun/Desktop/T_total_3_26.rds")
#载入NK细胞
NK_cells<-readRDS(file="./NK_cells_update.rds")
#载入B细胞
B_cells<-readRDS(file="./B_cells_11_15.rds")
#载入髓系细胞
Myeloid_cells<-readRDS(file="./macrophages_update.rds")
#载入肥大细胞
Mast_cells<-readRDS(file="./Mast_cells.rds")
#载入内皮细胞
Endothelial<-readRDS(file="./Endothelial.rds")
#载入成纤维细胞
Fibroblast<-readRDS(file="./Fibroblast.rds")

#载入浆细胞
Plasma_cells<-readRDS(file="./Plasma_cells.rds")



#修改各个亚群的Idents
table(Idents(epcam))
new.cluster.ids <- c("Epithelial cells","Epithelial cells")
names(new.cluster.ids) <- levels(epcam)
epcam <- RenameIdents(epcam, new.cluster.ids)




table(Idents(T_cells))
new.cluster.ids <- c("T cells","T cells","T cells","T cells","T cells","T cells","T cells")
names(new.cluster.ids) <- levels(T_cells)
T_cells <- RenameIdents(T_cells, new.cluster.ids)


table(Idents(NK_cells))
new.cluster.ids <- c("NK cells","NK cells")
names(new.cluster.ids) <- levels(NK_cells)
NK_cells <- RenameIdents(NK_cells, new.cluster.ids)


table(Idents(B_cells))
new.cluster.ids <- c("B cells","B cells","B cells","B cells","B cells","B cells")
names(new.cluster.ids) <- levels(B_cells)
B_cells <- RenameIdents(B_cells, new.cluster.ids)


table(Idents(Myeloid_cells))
new.cluster.ids <- c("Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells")
names(new.cluster.ids) <- levels(Myeloid_cells)
Myeloid_cells <- RenameIdents(Myeloid_cells, new.cluster.ids)

table(Myeloid_cells$orig.ident)






table(Idents(Mast_cells))
new.cluster.ids <- c("Mast cells","Mast cells","Mast cells","Mast cells","Mast cells","Mast cells","Mast cells","Mast cells","Mast cells")
names(new.cluster.ids) <- levels(Mast_cells)
Mast_cells <- RenameIdents(Mast_cells, new.cluster.ids)


table(Idents(Endothelial))
new.cluster.ids <- c("Endothelial cells","Endothelial cells","Endothelial cells","Endothelial cells")
names(new.cluster.ids) <- levels(Endothelial)
Endothelial <- RenameIdents(Endothelial, new.cluster.ids)


table(Idents(Fibroblast))
new.cluster.ids <- c("Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts")
names(new.cluster.ids) <- levels(Fibroblast)
Fibroblast <- RenameIdents(Fibroblast, new.cluster.ids)



table(Idents(Plasma_cells))
new.cluster.ids <- c("Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells")
names(new.cluster.ids) <- levels(Plasma_cells)
Plasma_cells <- RenameIdents(Plasma_cells, new.cluster.ids)



pbmc<-merge(epcam,y=c(T_cells,NK_cells,Myeloid_cells,Mast_cells,B_cells,Plasma_cells,Endothelial,Fibroblast),add.cell.ids = NULL)

saveRDS(pbmc,file="./pbmc_merge_new.rds")

table(Idents(pbmc))

table(pbmc$orig.ident)


table(T_cells$orig.ident)









#圆圈折线图绘制
percent_data<-read.table("/Users/wangjun/Desktop/更新结果/macrophage//macrophage_percentage.csv",sep = ",",header = T,row.names = 1)
percent_data[1:2,]

names<-rownames(percent_data)

p=wilcox.test(as.numeric(percent_data[1,c(2,4,6,8,10,14)]),as.numeric(percent_data[1,c(1,3,5,7,9,13,11,12,15,16)]))
p$p.value

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]))
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T,alternative = "less")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  
  # 创建示例数据
  data <- data.frame(group = c("MIA", "MIA", "MIA", "MIA","MIA","MIA","MIA","MIA","MIA","MIA"),
                     value = c(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),NA,NA,NA,NA),
                     connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                     connect_value = as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]))
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  # 创建分组箱形图
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    ggtitle(paste0("P-value : ",round(p$p.value,4)))+    
    xlab("Tissue")+
    ylab(rownames(percent_data)[i])+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  p <- p + scale_x_discrete(limits = c("MIA", "LUAD"))
  p
  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/macrophage/",names[i],"_greater.pdf"),width = 2.12,height = 3.44)
  
}



percent_data<-read.table("/Users/wangjun/Desktop/更新结果/cluster/大类比例.csv",sep = ",",header = T,row.names = 1)
percent_data[,1:3]


names<-c("Epithelial cells","T cells","Myeloid cells","NK cells","B cells","Mast cells","Plasma cells","Endothelial cells","Fibroblasts")


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T,alternative = "less")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}




#for (i in 1:length(rownames(percent_data))){
for (i in 5:5){
  # 创建示例数据
  data <- data.frame(group = c("MIA", "MIA", "MIA", "MIA","MIA","MIA","MIA","MIA","MIA","MIA"),
                     value = c(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),NA,NA,NA,NA),
                     connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                     connect_value = as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]))
  #p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  
  # 创建分组箱形图
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    ggtitle(paste0("P-value : ",round(p$p.value,4)))+
    xlab("Tissue")+
    #ylab(rownames(percent_data)[i])+
    ylab("B cells(4.36%)")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  p <- p + scale_x_discrete(limits = c("MIA", "LUAD"))
  p
  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/cluster/",names[i],"_paired.pdf"),width = 2.12,height = 3.44)
  
}





#开始统计T细胞的比例数据
table(T_cells$new)
Idents(T_cells)<-"new"

data<-data.frame(matrix(NA,7,16))
data[1,]<-table(subset(T_cells,idents="Naïve T cells(66.57%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[2,]<-table(subset(T_cells,idents="CD8+ GZMK+ T cells(10.63%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[3,]<-table(subset(T_cells,idents="CD8+ GZMB+ T cells(10.14%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[4,]<-table(subset(T_cells,idents="Treg cells(10.13%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[5,]<-table(subset(T_cells,idents="Exhausted T cells(0.70%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[6,]<-table(subset(T_cells,idents="γδT_C1(1.21%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]
data[7,]<-table(subset(T_cells,idents="γδT_C2(0.63%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]


write.table(data,file = "/Users/wangjun/Desktop/更新结果/T.csv",sep=",")

data

table(subset(T_cells,idents="Naïve T cells(66.57%)")$orig.ident)

table(subset(T_cells,idents="Naïve T cells(66.57%)")$orig.ident)[c("T5","T4","T1","T2","T6","T7","T8","T9","T12","T13","T10","T11","T14","T15","T16","T17")]

table(subset(T_cells,idents="γδT_C2(0.63%)")$orig.ident)

table(T_cells$orig.ident)



percent_data<-read.table("/Users/wangjun/Desktop/更新结果/T_cells/T.csv",sep = ",",header = T,row.names = 1)
percent_data[1:2,]

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T,alternative = "less")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}

paste0(strsplit(rownames(percent_data)[2],split = " ")[[1]][1:2])
strsplit()

name<-c("Naïve T cells","CD8+ GZMK+ T cells","CD8+ GZMB+ T cells","Treg cells","Exhausted T cells","gdT_C1","gdT_C2")

i=1
#for (i in 1:length(rownames(percent_data))){
for (i in 7:7){ 
  # 创建示例数据
  data <- data.frame(group = c("MIA", "MIA", "MIA", "MIA","MIA","MIA","MIA","MIA","MIA","MIA"),
                     value = c(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),NA,NA,NA,NA),
                     connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                     connect_value = as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]))
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  #p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  # 创建分组箱形图
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    ggtitle(paste0("P-value : ",round(p$p.value,4)))+
    xlab("Tissue")+
    #ylab(rownames(percent_data)[i])+
    ylab("gdT_C2(0.63%)")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  p <- p + scale_x_discrete(limits = c("MIA", "LUAD"))
  p
  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/T_cells/",name[i],"_greater.pdf"),width = 2.12,height = 3.44)
  
}






percent_data<-read.table("/Users/wangjun/Desktop/更新结果/NK_cells/NK_cells.csv",sep = ",",header = T,row.names = 1)
percent_data[1:2,]

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "greater")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}

for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


for (i in 1:length(rownames(percent_data))){
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T,alternative = "less")
  
  print(paste0(rownames(percent_data)[i]," : ",p$p.value))
}


name<-c("CD56dimCD16+NK cells")


for (i in 1:length(rownames(percent_data))){
  
  # 创建示例数据
  data <- data.frame(group = c("MIA", "MIA", "MIA", "MIA","MIA","MIA","MIA","MIA","MIA","MIA"),
                     value = c(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),NA,NA,NA,NA),
                     connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                     connect_value = as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]))
  
  p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")
  #p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13)]),paired = T)
  # 创建分组箱形图
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    ggtitle(paste0("P-value : ",round(p$p.value,4)))+
    xlab("Tissue")+
    ylab("CD56dimCD16+K cells(90.49%)")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  p <- p + scale_x_discrete(limits = c("MIA", "LUAD"))
  p
  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/更新结果/NK_cells/",name[i],"_less.pdf"),width = 2.12,height = 3.44)
  
}






#可视化MIA和LUAD两组种的silent和non_silent突变情况

MIA_silent<-c(1.209121654,0.748170745,1.171395969,0.879134664,0.741663028,1.156047338)

MIA_non_silent<-c(1.023227259,0.627257208,1.108750346,0.798760982,0.655313108,1.146946534)

wilcox.test(MIA_silent,MIA_non_silent,paired = T,alternative = "greater")



LUAD_silent<-c(0.988264041,0.90131913,0.943735919,1.010454331,0.960056584,0.958864341,0.937850281,0.893168901,0.978205971,1.090418042)

LUAD_non_silent<-c(0.953452563,0.909143393,0.925974667,1.260663493,0.941984614,1.126805002,0.811522695,0.800924862,0.935742424,1.038854908)

wilcox.test(LUAD_silent,LUAD_non_silent,paired = T,alternative = "greater")

wilcox.test(LUAD_silent,MIA_silent,alternative = "greater")

wilcox.test(LUAD_non_silent,MIA_non_silent,alternative = "greater")


library(ggplot2)
for (i in 1:length(rownames(percent_data))){
  
  # 创建示例数据
  data <- data.frame(group = c("MIA_silent","MIA_silent","MIA_silent","MIA_silent","MIA_silent","MIA_silent"),
                     value = c(as.numeric(MIA_silent)),
                     connect_group = c("MIA_non_silent","MIA_non_silent","MIA_non_silent","MIA_non_silent","MIA_non_silent","MIA_non_silent"),
                     connect_value = as.numeric(MIA_non_silent))
  
  #data <- data.frame(group = c("LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent"),
  #                   value = c(as.numeric(LUAD_silent)),
  #                   connect_group = c("LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent"),
  #                   connect_value = as.numeric(LUAD_non_silent))
  
  
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    #ggtitle(paste0("P-value : 0.01563"))+
    ggtitle(paste0("P-value : 0.01563"))+
    xlab("Mutation type")+
    ylab("Gene expression rate")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  #p <- p + scale_x_discrete(limits = c("MIA_silent", "MIA_non_silent"))
  p <- p + scale_x_discrete(limits = c("MIA_silent", "MIA_non_silent"))
  
  p <- p+theme(axis.text.x = element_text(angle = 60,hjust = 0.8,vjust = 0.8)) 

  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/MIA_silent.pdf"),width = 2.12,height = 4.5)
  
}



MIA_neo_exp<-c("0.00133980243522386","-0.00778282057081626",-0.00900144442922002,"0.00143224302738983","0.009240162196046")
LUAD_neo_exp<-c("-0.0119153716769999",-0.0140605621295198,"-0.0522850686551056","-0.0433967387739034","-0.00213491823983392",-0.047735717398372,"-0.049092356836148","0.03447649555235","0.0014202156078874","0.104798733518385")

MIA_neo_exp<-c(-0.0119153716769999,-0.00778282057081626,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(0.00133980243522386,-0.0219726189882752,-0.0522850686551056,-0.0897094438146464,-0.049092356836148,0.104798733518385,-0.0433967387739034,-0.00213491823983392,0.03447649555235,0.0014202156078874)


MIA_neo_exp<-c(-0.0119153716769999,-0.00778282057081626,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(0.00133980243522386,-0.0219726189882752,-0.0522850686551056,-0.0897094438146464,-0.049092356836148)


wilcox.test(MIA_neo_exp,LUAD_neo_exp,paired = T,alternative = "greater")

wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater")



#修订之后
MIA_neo_exp<-c(-0.012113429326935,-0.0167801548100993,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.00361242206404616,-0.0219726189882752,-0.097666121143288,-0.0897094438146464,-0.0960410431710384,0.148025326015436,-0.0433967387739034,-0.00363452727604194,0.00509901536902821,0.00142336464471198)

MIA_neo_exp<-c(-0.012113429326935,-0.0167801548100993,-0.0926577713344851,0.00143224302738983,0.009240162196046)
LUAD_neo_exp<-c(-0.00361242206404616,-0.0219726189882752,-0.097666121143288,-0.0897094438146464,-0.0960410431710384)
wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater",paired = T)


wilcox.test(MIA_neo_exp,LUAD_neo_exp,alternative = "greater")


unique(data_total[data_total$SampleID%in%"FD2_LC1",]$Gene)


library(ggplot2)
for (i in 1:length(rownames(percent_data))){
  
  # 创建示例数据
  data <- data.frame(group = c("MIA_neo_exp","MIA_neo_exp","MIA_neo_exp","MIA_neo_exp","MIA_neo_exp"),
                     value = c(as.numeric(MIA_neo_exp),NA,NA,NA,NA,NA),
                     connect_group = c("LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp","LUAD_neo_exp"),
                     connect_value = as.numeric(LUAD_neo_exp))
  
  #data <- data.frame(group = c("LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent","LUAD_silent"),
  #                   value = c(as.numeric(LUAD_silent)),
  #                   connect_group = c("LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent","LUAD_non_silent"),
  #                   connect_value = as.numeric(LUAD_non_silent))
  
  
  p <- ggplot(data, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    #ggtitle(paste0("P-value : 0.01563"))+
    ggtitle(paste0("P-value : 0.1562"))+
    xlab("Tiisue type")+
    ylab("Neoantigen expression")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=14),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
  #p <- p + scale_x_discrete(limits = c("MIA_silent", "MIA_non_silent"))
  p <- p + scale_x_discrete(limits = c("MIA_neo_exp", "LUAD_neo_exp"))
  
  p <- p+theme(axis.text.x = element_text(angle = 60,hjust = 0.8,vjust = 0.8)) 
  
  ggsave(p,file=paste0("/Users/wangjun/Desktop/Neoantigen_exp_paired.pdf"),width = 2.12,height = 4.5)

}










#绘制figure_2d的图
library(reshape2)
library(ggplot2)
data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/Fig_2D_细胞比例计算.csv",header=T,row.names = 1,sep = ",")

data
data[1:4,1:4]


colSums(data)

for (i in 1:length(colnames(data))){
  data_new[,i]<-data[,i]/colSums(data)[i]
}

length(colnames[data])

colSums(data_new)
(colSums(data)[8])

rownames(data_new)

data_new<-data_new[c("Endothelial cells","Fibroblasts","Epithelial cells","mast cells","Myeloid cells","Plasma cells","B cells","NK cells","T cells"),]

data_melt<-melt(t(data_new))

data_melt[1:4,]

colnames(data_melt)<-c("Sample_type","Cell_type","Percent")


P1_D<-ggplot(data_melt, aes(x= Sample_type, y = Percent, fill = Cell_type))+ ## 使用ggplot2语法
  geom_bar(stat = "identity", position = "stack",width = 0.7)+ ## 添加柱子
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),#ggplot2绘图，有点粗糙，见谅
        axis.text.x.bottom = element_text(angle = 60, hjust = 1,size = 10),
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values = pal[c(2,7,10,4,6,9,3,5,1)])

ggsave(P1_D,file="/Users/wangjun/Desktop/最终方案/figure/Figure1_D.pdf",height = 4.23,width = 8.43)


library(paletteer) 
library(scales)

pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal_new<-(colorRampPalette(pal)(13))
show_col(pal_new)

pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 
pal[c(2,3,7,10,5,8,4,1,9)]

pal[c(1,10,6,5,7,2,3,4,9)]


pal[c(1:13)]


pal[c(2,7,10,4,6,9,3,5,1)]

pal_new[c()]


show_col("#FFB2A8")


#绘制figure4B
library(Seurat)


T_cell<-readRDS(file="/Users/wangjun/Desktop/T_total_3_31.rds")



#鉴定marker
marker<-c(
  "IL7R","CD4",#CD4_T
  "CD8A","CD8B",#CD8_T
  "TCF7","SELL","LEF1","CCR7",#naive_t
  "LAG3","TIGIT","PDCD1","CTLA4","KLRG1",#exhausted
  "IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7",#cytotoxic
  "IL2RA","FOXP3","IKZF2","TGFB1","TGFB3","TGFBI","TGFBR1",#treg
  "MAF","CXCR5","CXCL13", # and PDCD1 #tfh
  "IRF4","CREM","NR4A2",#th17
  "STAT4","IL12RB2", # and IFNG  #th1
  "GATA3","STAT6","IL4",#th2
  "TRDC","TRGC2","TRGC1","CD38","MKI67",#gd_t
  #"ITGAE",#Trm
  "CXCR4"#Tem
  #"CXCR6","PTGER4" #CD4 Trm
  #"IL15RA","CD44","MBD2","BLIMP1"#memory
)

table(Idents(T_cell))
table(T_cell$Annotation)


marker<-c(
  "IL7R","CD4",#CD4_T
  "CD8A","CD8B",#CD8_T
  "TCF7","SELL","LEF1","CCR7",#naive_t
  #"LAG3","TIGIT","PDCD1","CTLA4","KLRG1",#exhausted
  "LAG3","TIGIT","PDCD1","CTLA4", "HAVCR2", "ENTPD1",#exhausted
  "IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7",#cytotoxic
  "IL2RA","FOXP3","IKZF2","TGFB1","TGFB3","TGFBI","TGFBR1",#treg
  "MAF","CXCR5","CXCL13", # and PDCD1 #tfh
  "IRF4","CREM","NR4A2",#th17
  "STAT4","IL12RB2", # and IFNG  #th1
  "GATA3","STAT6","IL4",#th2
  "TRDC","TRGC2","TRGC1","CD38","MKI67",#gd_t
  #"ITGAE",#Trm
  "CXCR4"#Tem
  #"CXCR6","PTGER4" #CD4 Trm
  #"IL15RA","CD44","MBD2"#memory
)




new.cluster.ids <- c("Naïve T cells","Treg cells","Th cells","Memory T cells","CD8+ GZMK+ T cells","CD8+ GZMB+ T cells","Exhausted T cells","γδT_C1","γδT_C2")
names(new.cluster.ids) <- levels(T_cell)
T_cell <- RenameIdents(T_cell, new.cluster.ids)



install.packages("extrafont")
library(extrafont)



T_dotplot<-DotPlot(T_cell, features = marker, cols = c("#4DBBD5FF","#E64B35FF")) + RotatedAxis()

ggsave(T_dotplot,file="/Users/wangjun/Desktop/最终方案/figure/Figure4_B.pdf",height = 3.95,width = 14)









table(Idents(T_total))
cairo_pdf(filename="/Users/wangjun/Desktop/最终方案/figure/Figure4_B.pdf",height = 3.95,width = 14)

T_dotplot
dev.off()


p5<-DimPlot(T_cell,reduction = 'tsne', label.box = T,  label = F,repel = T,cols = )

ggsave(p5,file="/Users/wangjun/Desktop/T_dimplot.png",dpi = 1000,width = 7.45,height = 4.67)





#TNSE
T_total <- RunTSNE(T_cell, reduction = "pca", dims = 1:20)
p5 <- DimPlot(T_cell, reduction = "tsne",cols = pal[c(1:9)])


ggsave(file="/Users/wangjun/Desktop/最终方案/figure/Figure4_A.png",p5,height = 4,24,width = 6.98)


pal <- paletteer_d("ggsci::nrc_npg")[c(1:10)] 

dev.off()







#更新TCR多样性的图片


data<-read.table(file="/Users/wangjun/Desktop/最终方案/figure/TCR_tumor.csv",header = T,sep = ",")


data_melt<-melt(t(data))[,c(1,3)]

colnames(data_melt)<-c("group","value")

data_new <- data.frame(group = c("Low","Low","Low","Low","Low","Low","Low","Low"), 
                   value = c(as.numeric(data$Low)),
                   connect_group = c("High","High","High","High","High","High","High","High"),
                   connect_value = as.numeric(data$High))
#p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")


p=wilcox.test(as.numeric(data[,c(2)]),as.numeric(data[,c(1)]),paired = T)
# 创建分组箱形图
p

p1 <- ggplot(data_new, aes(x = group, y = value, group = group)) +
    #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
    ggtitle(paste0("P-value : ",round(p$p.value,4)))+    
    xlab("TCR diversity")+
    ylab("Tumor heterogeneity")+
    # 绘制连接线
    geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
                 color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
    # 绘制数据点
    geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
    geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
    theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=12))+
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
p2 <- p1 + scale_x_discrete(limits = c("Low", "High"))
p2
  
ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/Figure5_A.pdf",width = 2.12,height = 3.44)
  







MIA_neo<-c(7,0,61,7,2,29)
LUAD_neo<-c(8,7,83,17,6,37,30,12,6,2)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less")


data_new <- data.frame(group = c("MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA"), 
                       value = c(as.numeric(MIA_neo),NA,NA,NA,NA),
                       connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                       connect_value = as.numeric(LUAD_neo))
#p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")


p=wilcox.test(as.numeric(MIA_neo),as.numeric(LUAD_neo[1:6]),paired = T,alternative = "less")
# 创建分组箱形图
p

p1 <- ggplot(data_new, aes(x = group, y = value, group = group)) +
  #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
  ggtitle(paste0("P-value : ",round(p$p.value,4)))+    
  xlab("Tissue")+
  ylab("Neoantigen_count")+
  # 绘制连接线
  geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
               color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
  # 绘制数据点
  geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
  geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
p2 <- p1 + scale_x_discrete(limits = c("MIA", "LUAD"))
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/Figure5_E.pdf",width = 2.12,height = 3.44)




#

MIA_neo<-c(8,0,119,33,10,73)
LUAD_neo<-c(48,29,154,166,10,74)
wilcox.test(MIA_neo,LUAD_neo,alternative = "less",paired = T)


MIA_neo<-c(8,0,119,33,10,73)
LUAD_neo<-c(48,29,154,166,10,74,7,25,41,5)


data_new <- data.frame(group = c("MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA","MIA"), 
                       value = c(as.numeric(MIA_neo),NA,NA,NA,NA),
                       connect_group = c("LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD","LUAD"),
                       connect_value = as.numeric(LUAD_neo))
#p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")


p=wilcox.test(as.numeric(MIA_neo),as.numeric(LUAD_neo[1:6]),paired = T,alternative = "less")
# 创建分组箱形图
p

p1 <- ggplot(data_new, aes(x = group, y = value, group = group)) +
  #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
  ggtitle(paste0("P-value : ",round(p$p.value,4)))+    
  xlab("Tissue")+
  ylab("Neoantigen_count")+
  # 绘制连接线
  geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
               color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
  # 绘制数据点
  geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
  geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
p2 <- p1 + scale_x_discrete(limits = c("MIA", "LUAD"))
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/Figure5_F.pdf",width = 2.12,height = 3.44)





#Gene expression rate
silent<-c(1.209121654,0.988264041,0.90131913,0.748170745,0.943735919,1.171395969,1.010454331,0.879134664,0.960056584,0.958864341,0.937850281,0.741663028,0.893168901,1.156047338,0.978205971,1.090418042)
non_silent<-c(1.023227259,0.953452563,0.909143393,0.627257208,0.925974667,1.108750346,1.260663493,0.798760982,0.941984614,1.126805002,0.811522695,0.655313108,0.800924862,1.146946534,0.935742424,1.038854908)


data_new <- data.frame(group = c("Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent","Silent"), 
                       value = c(as.numeric(silent)),
                       connect_group = c("Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent","Non_silent"),
                       connect_value = as.numeric(non_silent))
#p=wilcox.test(as.numeric(percent_data[i,c(2,4,6,8,10,14)]),as.numeric(percent_data[i,c(1,3,5,7,9,13,11,12,15,16)]),alternative = "less")


p=wilcox.test(as.numeric(non_silent),as.numeric(silent),paired = T,alternative = "less")
# 创建分组箱形图
p

p1 <- ggplot(data_new, aes(x = group, y = value, group = group)) +
  #geom_boxplot(fill = "lightgray", width = 0.5, outlier.shape = NA) +
  ggtitle(paste0("P-value : ",round(p$p.value,4)))+    
  xlab("Mutation type")+
  ylab("Gene expression rate")+
  # 绘制连接线
  geom_segment(aes(x = group, xend = connect_group, y = value, yend = connect_value),
               color = "grey50", size = 1, linetype = "solid",alpha=0.8)+
  # 绘制数据点
  geom_point(aes(x = group, y = value), size = 2.5, color = "#4DBBD5FF") +
  geom_point(aes(x = connect_group, y = connect_value), size = 2.5, color = "#E64B35FF") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=12))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.8,0.3),legend.box.background=element_rect(colour = "white",fill = "white"))
p2 <- p1 + scale_x_discrete(limits = c("Silent", "Non_silent"))
p2

ggsave(p2,file="/Users/wangjun/Desktop/最终方案/figure/Figure5_J.pdf",width = 2.12,height = 3.44)











