library(ggplot2)
library(tidyverse)
library(ggrepel)

annotation_grch38 <- read.delim("annotation_grch38.bed")

#Jeju island component is 0
k=0
outquantile<-.999
filename<-paste0("small_scan_k",k,".txt")
scan.haenyeo <- read.delim(filename)
scan.haenyeo$CHROM<-paste0("chr",scan.haenyeo$chr)
#Number1: Find outliers, top 0.1%
cutoff<-quantile(scan.haenyeo $lle.ratio,outquantile)
hits<-scan.haenyeo[scan.haenyeo $lle.ratio>cutoff,]

#Number 2:the SNP is not the only one within ± 100kb (200kb window) that is among the most extreme 0.1% likelihood ratios.
window_200k<-hits %>% group_by(win_start=position-250000,win_end=position+250000)
regions<-window_200k[c("CHROM","win_start","win_end")]
temp<-merge(regions,hits,by="CHROM")
hits_windows<-temp[(temp$position>=temp$win_start) & (temp$position<=temp$win_end),]
#At least another 2 SNPs that are also outliers
temp2<-hits_windows[c("CHROM","win_start","win_end","position","lle.ratio","ID")]
hitsPerWindow2<-temp2%>% group_by(CHROM,win_start)%>%filter(n() > 1)
ld_hits2<-unique(hitsPerWindow2[c("CHROM","position","lle.ratio","ID")])

#3 The SNP does not have any neighbor (+-100kb) with higher likelihood ratio
top_neighbor<-hitsPerWindow2%>%group_by(CHROM,win_start)%>%top_n(1,lle.ratio)
top_neighbor2<-top_neighbor%>%group_by(CHROM,win_start)%>%top_n(1,position)
filtered_hits<-unique(top_neighbor2[c("CHROM","position","lle.ratio","ID")])



window_200k<-filtered_hits %>% group_by(win_start=position-250000,win_end=position+250000)
regions<-window_200k[c("CHROM","win_start","win_end")]
temp<-merge(regions,filtered_hits,by="CHROM")
hits_windows_2<-temp[(temp$position>=temp$win_start) & (temp$position<=temp$win_end),]
top_neighbor3<-hits_windows_2%>%group_by(CHROM,win_start)%>%top_n(1,lle.ratio)
filtered_hits2<-unique(top_neighbor3[c("CHROM","position","lle.ratio","ID")])



sorted_hits<-filtered_hits2[order(-filtered_hits2$lle.ratio),]
top10SNPs<-head(sorted_hits,10)

write.table(top10SNPs,paste0("top10SNPs_",outquantile,".txt"),quote = F,append = F,row.names = F)

write.table(sorted_hits,paste0("topSNPs_",outquantile,".txt"),quote = F,append = F,row.names = F)



#4 Annotate
anno<-merge(filtered_hits,annotation_grch38,by="CHROM")
annotated<-anno[(anno$position>=anno$start) & (anno$position<=anno$end),]

#Top SNP of each annotation
anno_reduced<-annotated %>% group_by(gname) %>% top_n(1, lle.ratio)
anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, position)
anno_reduced<-unique(anno_reduced[c("CHROM","position","lle.ratio","gname")])




colors<-c("#56B4E9","#7FCC12", '#8B2121','#008000','#F0E442', 
          '#0000FF','#FF0000','#FFA500')

anno_collapsed<-anno_reduced %>%
  group_by(CHROM,position,lle.ratio) %>%
  summarize(gname = str_c(gname, collapse = ", "))
noloc<-anno_collapsed[!startsWith(anno_collapsed$gname,"LOC"),]
half_data<-scan.haenyeo[scan.haenyeo$lle.ratio>quantile(scan.haenyeo$lle.ratio,.5),]
half_anno<-merge(half_data,noloc,all.x = T)

newdata <- half_anno[order(parse_number(half_anno$CHROM)),] 

n_of_chr<-length(unique(newdata$CHROM))
newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))


ggplot(newdata,aes(x=SNPs,y=lle.ratio,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+
  ggtitle(paste0("Selection on Jeju component"))+theme_classic(base_size=14)+
  ylim(c(0,NA))  +scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(colors[k+1], "gray"), unique(length(axis_set$CHROM))))+
  geom_label_repel(size=3,max.overlaps = 60,fontface = 'bold') +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))


ggsave(filename=paste0("outlier_ldhits_noloc_min3",outquantile,"_k",k,"selscan_anno.png"),width=10,height=4)

