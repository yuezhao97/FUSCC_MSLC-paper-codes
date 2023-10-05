setwd("/Users/zhaoy2/Desktop/sc_project/inferCNV")

file1 <- read.delim("final_frequency_all.txt", sep="\t", header=T,row.names=1)
head(file1)
file1[file1 >= 0.1] <- 1
file1[file1 < 0.1] <- 0
file1 <- file1[,-2]
head(file1)

all <- c("T1","T2","T4","T5","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17")
FD1 <- c("T4","T5")
FD2 <- c("T1","T2")
FD4 <- c("T6","T7")
FD5 <- c("T8","T9")
FD8 <- c("T12","T13")
FD9 <- c("T10","T11")
FD14 <- c("T14","T15")
FD16 <- c("T16","T17")

dist_all <- c()
for (i in all){
  col1 <- as.data.frame(file1[,i])
  dist_row <- c()
  for (j in all){
    col2 <- as.data.frame(file1[,j])
    col_to_compare <- cbind(col1, col2)
    colnames(col_to_compare) <- c(i,j)
    rownames(col_to_compare) <- rownames(file1)
    #print(head(col_to_compare))
    distance2 <- 0
    for (r in 1:nrow(col_to_compare)){
      distance2 <- distance2 + (col_to_compare[r,1] - col_to_compare[r,2])^2
    }
    distance <- sqrt(distance2)
    dist_row <- rbind(dist_row, distance)
  }
  dist_all <- cbind(dist_all, dist_row)
}

colnames(dist_all) <- all
rownames(dist_all) <- all
head(dist_all)
write.table(dist_all, file="genomic_distance/genomic_distance_infercnv.txt", sep="\t", quote=F)

self_dists <- c(dist_all["T1","T2"], dist_all["T4","T5"], dist_all["T6","T7"], dist_all["T8","T9"], dist_all["T12","T13"], dist_all["T10","T11"], dist_all["T14","T15"], dist_all["T16","T17"])

class(self_dists)
self_dists <- as.data.frame(self_dists)
colnames(self_dists) <- "Genomic_distance"
self_dists$Group <- "Self"
other_dists <- c(dist_all["T1","T4"], dist_all["T1","T7"], dist_all["T1","T9"],dist_all["T1","T13"],dist_all["T1","T11"],dist_all["T1","T15"],dist_all["T1","T17"],
                 dist_all["T5","T2"],dist_all["T5","T7"],dist_all["T5","T9"],dist_all["T5","T13"],dist_all["T5","T11"],dist_all["T5","T15"],dist_all["T5","T17"],
                 dist_all["T6","T4"],dist_all["T6","T2"],dist_all["T6","T9"],dist_all["T6","T13"],dist_all["T6","T11"],dist_all["T6","T15"],dist_all["T6","T17"],
                 dist_all["T8","T4"],dist_all["T8","T2"],dist_all["T8","T7"],dist_all["T8","T13"],dist_all["T8","T11"],dist_all["T8","T15"],dist_all["T8","T17"],
                 dist_all["T12","T4"],dist_all["T12","T2"],dist_all["T12","T7"],dist_all["T12","T9"],dist_all["T12","T11"],dist_all["T12","T15"],dist_all["T12","T17"],
                 dist_all["T10","T4"],dist_all["T10","T2"],dist_all["T10","T7"],dist_all["T10","T9"],dist_all["T10","T13"],dist_all["T10","T15"],dist_all["T10","T17"],
                 dist_all["T14","T4"],dist_all["T14","T2"],dist_all["T14","T7"],dist_all["T14","T9"],dist_all["T14","T13"],dist_all["T14","T11"],dist_all["T14","T17"],
                 dist_all["T16","T4"],dist_all["T16","T2"],dist_all["T16","T7"],dist_all["T16","T9"],dist_all["T16","T13"],dist_all["T16","T11"],dist_all["T16","T15"])
other_dists <- as.data.frame(other_dists)
colnames(other_dists) <- "Genomic_distance"
other_dists$Group <- "Other"
out_dists <- rbind(self_dists, other_dists)
out_dists <- out_dists[,c(2,1)]
head(out_dists)
out_dists$Group <- factor(out_dists$Group, levels=c("Self","Other"))

library(ggplot2)
library(ggpubr)
ggplot(data=out_dists, aes(x=Group, y=Genomic_distance, fill=Group)) +
  geom_boxplot() +
  labs(x="Group", y = "Genomic distance", fill="Group") +
  #scale_fill_brewer(palette="Blues") + 
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
  scale_fill_manual(values = c("pink","lightblue"))
  #scale_fill_brewer(palette="Set1")
ggsave("genomic_distance/Genomic_distance_self_vs_other.pdf", plot=last_plot(), width=6, height=8, dpi=600)
