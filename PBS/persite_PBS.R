setwd("~/Documents/Berkeley/haenyeo/pbsScanv11/")
library(ggplot2)
library(stringr)
library(ggrepel)
library(tidyverse)
library(latex2exp)

custom_title<-"Jeju Island (Haenyeo+Jeju) vs Seoul vs Chinese"
pbs1<-"HAE_JEJU"
outquantile<-.999
color1<-"blue"

persitefst <- read.delim("hsc4pbs", header=TRUE)
persitefst <- read.delim("hsc4pbsMC", header=TRUE)

#Nothing!
pos2check<-persitefst[persitefst$WEIR_AND_COCKERHAM_FST==1,]
pos2check<-persitefst[persitefst$WEIR_AND_COCKERHAM_FST.1==1,]
pos2check<-persitefst[persitefst$WEIR_AND_COCKERHAM_FST.2==1,]

#Convert 1 to .999
persitefst$WEIR_AND_COCKERHAM_FST.1[persitefst$WEIR_AND_COCKERHAM_FST.1==1]<-.99
persitefst$WEIR_AND_COCKERHAM_FST.2[persitefst$WEIR_AND_COCKERHAM_FST.2==1]<-.99
persitefst$WEIR_AND_COCKERHAM_FST[persitefst$WEIR_AND_COCKERHAM_FST==1]<-.99

persitefst$fst12 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST))
persitefst$fst13 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST.1))
persitefst$fst23 <- (-log(1-persitefst$WEIR_AND_COCKERHAM_FST.2))

persitefst$pbs1 = (fst12 + fst13 - fst23) / 2
persitefst$pbs2 = (fst12 + fst23 - fst13) / 2
persitefst$pbs3 = (fst13 + fst23 - fst12) / 2

#quantile(persitefst$pbs1)
#0%          25%          50%          75%         100% 
#-2.095085120 -0.006527248 -0.002356547  0.008298388  1.641828180 


pbs_fil<-persitefst[persitefst$pbs1>0,]

#Optional, run if you did not convert ones to .999
gezero<-persitefst[persitefst$pbs1>0,]
pbs_fil<-gezero[gezero$pbs1<10,]

#quantile(pbs_fil$pbs1)
#0%          25%          50%          75%         100% 
#1.690610e-09 4.545145e-03 1.184028e-02 2.532697e-02 8.771309e-01 

toplot<-quantile(pbs_fil$pbs1,.75)
pbs_fil2plot<-pbs_fil[pbs_fil$pbs1>toplot,]

################chr6#################
chr6_only<-pbs_fil2plot[pbs_fil2plot$CHROM=="6",]
cutoff<-quantile(chr6_only$pbs1,.99)
chr6hits<-chr6_only[chr6_only$pbs1>=.11,]
chr6hits$CHROM<-paste0("chr",chr6hits$CHROM)
write.table(chr6hits[c("CHROM","POS","POS","pbs1")],
            paste0(pbs1,"_pbs_outliers_top_persite_chr6",outquantile,".bed"),
            quote = F,sep = "\t",row.names = F)
annotatedchr6<- read.delim(paste0("HAE_JEJU_pbs_anno25_chr6.bed"), header=FALSE)
colnames(annotatedchr6)<-c("CHROM","BIN_START","BIN_END","pbs1","gene")

anno_reduced<-annotatedchr6 %>% group_by(gene) %>% top_n(1, pbs1)%>% top_n(1, BIN_START)
anno_commas<-anno_reduced %>%
  group_by(CHROM,BIN_START,BIN_END,pbs1)%>%
  summarize(gene= str_c(gene, collapse = "\n "), .groups = 'drop')

chr6_only$CHROM<-paste0("chr",chr6_only$CHROM)
temp<-merge(chr6_only,anno_commas,by.x=c("CHROM","POS","pbs1"),by.y=c("CHROM","BIN_END","pbs1"),all.x = T)


newdata <- temp[order(parse_number(temp$CHROM)),] 


newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))

h<-ggplot(newdata, aes(x=SNPs,y=pbs1,label=gene,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "gray"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))+
  xlab("Chromosome")+ylab(TeX("$PBS_{HAE+JEJU}$"))    +
  #geom_text_repel(size=3,max.overlaps = 50,fontface = 'bold',ylim = c(.23,.45))
  geom_label_repel(size=3,max.overlaps = 60,fontface = 'bold',ylim = c(.5,1))


ggsave(paste0("anno_",outquantile,pbs1,"_pbs1_persite_chr6.png"),width=10,height=5)


nolabel<-ggplot(newdata, aes(x=SNPs,y=pbs1,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "gray"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))
  xlab("Chromosome")+ylab(TeX("$PBS_{HAE+JEJU}$"))  


ggsave(paste0(outquantile,pbs1,"_pbs1_persite.png"),width=10,height=5)

###########Find outliers PBS########
cutoff<-quantile(pbs_fil2plot$pbs1,.9999)
hits<-pbs_fil2plot[pbs_fil2plot$pbs1>=cutoff,]
#hits<-pbs_fil2plot[pbs_fil2plot$pbs1>=.23,]
outquantile<-"pbs1_23"
outquantile<-"pbs2_2"
hits$CHROM<-paste0("chr",hits$CHROM)
write.table(hits[c("CHROM","POS","POS","pbs1")],
            paste0(pbs1,"_pbs_outliers_top_persite",outquantile,".bed"),
            quote = F,sep = "\t",row.names = F)

annotated<- read.delim(paste0("HAE_JEJU_pbs_anno25.bed"), header=FALSE)
annotated<- read.delim(paste0("HAE_JEJU_pbs_anno23.bed"), header=FALSE)
annotated<- read.delim(paste0("HAE_JEJUpbs2_anno2.bed"), header=FALSE)

#annotated<- read.delim(paste0("HAE_JEJU_pbs_anno3.48.bed"), header=FALSE)
colnames(annotated)<-c("CHROM","BIN_START","BIN_END","pbs1","gene")

anno_reduced<-annotated %>% group_by(gene) %>% top_n(1, pbs1)%>% top_n(1, BIN_START)
anno_commas<-anno_reduced %>%
  group_by(CHROM,BIN_START,BIN_END,pbs1)%>%
  summarize(gene= str_c(gene, collapse = "\n "), .groups = 'drop')

pbs_fil2plot$CHROM<-paste0("chr",pbs_fil2plot$CHROM)
temp<-merge(pbs_fil2plot,anno_commas,by.x=c("CHROM","POS","pbs1"),by.y=c("CHROM","BIN_END","pbs1"),all.x = T)


newdata <- temp[order(parse_number(temp$CHROM)),] 

n_of_chr<-length(unique(newdata$CHROM))
newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))

newdata$gene[newdata$gene=="LOC105372869"]<-"AVPR1B-DT"
newdata$gene[newdata$gene=="LOC105374216 TNIK"]<-"TNIK"

h<-ggplot(newdata, aes(x=SNPs,y=pbs1,label=gene,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "black"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))+
  xlab("Chromosome")+ylab(TeX("$PBS_{HAE+JEJU}$"))    +
  #geom_text_repel(size=3,max.overlaps = 50,fontface = 'bold',ylim = c(.23,.45))
  geom_label_repel(size=3,max.overlaps = 60,fontface = 'bold')


ggsave(paste0("anno_",outquantile,pbs1,"_pbs1_persite_missing.png"),width=10,height=5)


nolabel<-ggplot(newdata, aes(x=SNPs,y=pbs1,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "black"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))
  xlab("Chromosome")+ylab(TeX("$PBS_{HAE+JEJU}$"))  


ggsave(paste0(outquantile,pbs1,"_pbs1_persite_missing.png"),width=10,height=5)

########### Find outliers FST#########
cutoff<-quantile(pbs_fil2plot$WEIR_AND_COCKERHAM_FST,.999)
hits<-pbs_fil2plot[pbs_fil2plot$WEIR_AND_COCKERHAM_FST>=cutoff,]

hits<-pbs_fil2plot[pbs_fil2plot$WEIR_AND_COCKERHAM_FST>=.25,]

outquantile<-"fst01_25"
write.table(hits[c("CHROM","POS","POS","WEIR_AND_COCKERHAM_FST")],
            paste0(pbs1,"-SEOULfst_outliers_top_persite",outquantile,".bed"),
            quote = F,sep = "\t",row.names = F)


annotated<- read.delim(paste0("HAE_JEJU-SEOUL_fst_anno25.bed"), header=FALSE)
colnames(annotated)<-c("CHROM","BIN_START","BIN_END","pbs1","gene")

anno_reduced<-annotated %>% group_by(gene) %>% top_n(1, pbs1)%>% top_n(1, BIN_START)
anno_commas<-anno_reduced %>%
  group_by(CHROM,BIN_START,BIN_END,pbs1)%>%
  summarize(gene= str_c(gene, collapse = "\n "), .groups = 'drop')

colnames(anno_commas)<-c("CHROM","POS","BIN_END","WEIR_AND_COCKERHAM_FST","gene")
temp<-merge(pbs_fil2plot,anno_commas,by=c("CHROM","POS","WEIR_AND_COCKERHAM_FST"),all.x = T)

newdata <- temp[order(parse_number(temp$CHROM)),] 

n_of_chr<-length(unique(newdata$CHROM))
newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))



custom_title<-"Jeju Island (Haenyeo+Jeju) vs Seoul "
h<-ggplot(newdata, aes(x=SNPs,y=WEIR_AND_COCKERHAM_FST,label=gene,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "hotpink2"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))+
  xlab("Chromosome")+ylab(TeX("$FST_{HAE+JEJU-SEOUL}$"))    +
  #geom_text_repel(size=3,max.overlaps = 50,fontface = 'bold',ylim = c(.23,.45))
  geom_label_repel(size=3,max.overlaps = 60,fontface = 'bold')

h


ggsave(paste0("anno_",outquantile,pbs1,"-SEOUL_fst_persite.png"),width=10,height=5)


nolabel<-ggplot(newdata, aes(x=SNPs,y=WEIR_AND_COCKERHAM_FST,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+theme_classic(base_size=14)+ggtitle(custom_title)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c(color1, "hotpink2"), unique(length(axis_set$CHROM))))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  #ylim( c(0, max(alldata$pbs1) ))
  xlab("Chromosome")+ylab(TeX("$PBS_{HAE+JEJU}$"))  

nolabel
ggsave(paste0(outquantile,pbs1,"-SEOUL_fst_persite.png"),width=10,height=5)

