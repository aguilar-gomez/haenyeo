library(ggplot2)
library(tidyverse)
library(ggrepel)

#Bed format of GrCh38 assembly annotation
annotation_grch38 <- read.delim("annotation_grch38.bed")

k=0
outquantile<-.999
filename<-paste0("small_scan_k",k,".txt")
scan.haenyeo <- read.delim(filename)
scan.haenyeo$CHROM<-paste0("chr",scan.haenyeo$chr)
#Find outliers
cutoff<-quantile(scan.haenyeo $lle.ratio,outquantile)
#cutoff<-cuts[k+1]
hits<-scan.haenyeo[scan.haenyeo $lle.ratio>cutoff,]
anno<-merge(hits,annotation_grch38,by="CHROM")

annotated<-anno[(anno$position>=anno$start) & (anno$position<=anno$end),]

anno_reduced<-annotated %>% group_by(gname) %>% top_n(1, lle.ratio)
anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, start)
anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, end)
anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, position)
anno_reduced<-unique(anno_reduced[c("CHROM","position","lle.ratio","gname")])
anno_red<-anno_reduced[!startsWith(anno_reduced$gname,"LOC"),]

write.table(anno_reduced,paste("outliers_selscan",outquantile,".txt"),quote = F,append = F,row.names = F)

#colors=c('#008000','#FF0000','#7FCC12','#F0E442', '#0000FF','#FFA500','#8B2121')
colors<-c("#56B4E9","#7FCC12", '#8B2121','#008000','#F0E442', 
          '#0000FF','#FF0000','#FFA500')


half_anno<-merge(half_data,anno_red,all.x = T)

newdata <- half_anno[order(parse_number(half_anno$CHROM)),] 

n_of_chr<-length(unique(newdata$CHROM))
newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))


ggplot(newdata,aes(x=SNP,y=lle.ratio,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+
  ggtitle(paste0("selection on k=",k))+theme_classic(base_size=14)+
  ylim(c(0,NA))  +scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(colors[k+1], "gray"), unique(length(axis_set$CHROM))))+
  geom_label_repel(size=2,max.overlaps = 60,fontface = 'bold',color="black") +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
  #ylim( c(0, max(alldata$pbs1) ))


ggsave(filename=paste0("outlier",outquantile,"_k",k,"selscan_anno.png"),width=10,height=4)


