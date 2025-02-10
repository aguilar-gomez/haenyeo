library(ggplot2)
library(tidyverse)
library(ggrepel)
library(data.table)
setwd("~/Documents/Berkeley/Haenyeo/selscan_outCHS")

annotation_grch38 <- read.delim("annotation_grch38.bed")

k=1

scan.haenyeo <- read.delim("selscan_jejuComp_CHS_small.txt")
scan.haenyeo$CHROM<-paste0("chr",scan.haenyeo$chr)
onlyLLE <- read.csv("onlyLLE", sep="")


#Find outliers
cutoff<-quantile(scan.haenyeo$lle.ratio,.99)
#cutoff<-cuts[k+1]
hits<-scan.haenyeo[scan.haenyeo $lle.ratio>cutoff,]

#Number 2:the SNP is not the only one within ± 100kb (200kb window) that is among the most extreme 1% likelihood ratios.
window_200k<-hits %>% group_by(win_start=pos-100000,win_end=pos+100000)
regions<-window_200k[c("CHROM","win_start","win_end")]

# Create an empty list to store the results for each chromosome
results_list <- list()
# Define the range of chromosomes (1 to 22)
chromosomes <- 1:22

# Loop through each chromosome
for (k in chromosomes) {
  chrom=paste0("chr",k)
  print(k)
  print(chrom)
  # Subset regions and hits for the current chromosome
  regions_chr <- regions[regions$CHROM == chrom, ]
  hits_chr <- hits[hits$CHROM == chrom, ]
  
  # Perform the inner join for the current chromosome
  joined_data <- inner_join(regions_chr, hits_chr, by = "CHROM", suffix = c("_regions", "_hits"),relationship = "many-to-many")
  
  # Filter the joined data for the current chromosome
  hits_windows_chr <- joined_data %>%
    filter(pos >= win_start & pos <= win_end)
  
  temp2<-hits_windows_chr[c("CHROM","win_start","win_end","pos","lle.ratio","rs")]
  #remove filter to have at least two
  hitsPerWindow<-temp2%>% group_by(CHROM,win_start)%>%filter(n() > 1)
  
  #3 The SNP does not have any neighbor (+-100kb) with higher likelihood ratio
  top_neighbor<-hitsPerWindow%>%group_by(CHROM,win_start)%>%top_n(1,lle.ratio)
  df<-unique(top_neighbor[c("CHROM","pos","lle.ratio","rs")])
  
  filtered_hits <- df %>%
    group_by(group_id = cumsum(c(0, diff(pos) > 100000))) %>%
    filter(row_number() == which.max(lle.ratio)) %>%
    ungroup() %>%
    select(-group_id)
  # Store the results in the list
  results_list[[as.character(k)]] <-filtered_hits
  #results_list[[as.character(k)]] <-top_neighbor
}


merged_results <- bind_rows(results_list)

#Do top 50 because many are going to fall in the same gene
top50_rs_numbers<- head(merged_results[order(merged_results$lle.ratio, decreasing = TRUE), ] ,100)
#write.table(top50_rs_numbers,paste0("top100_August29",outquantile,".txt"),quote = F,append = F,row.names = F)


anno<-merge(top50_rs_numbers,annotation_grch38,by="CHROM",all.x = T)
annotated<-anno[(anno$pos>=anno$start) & (anno$pos<=anno$end),]


anno_reduced<-annotated %>% group_by(gname) %>% top_n(1, lle.ratio)
#anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, start)
#anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, end)
#anno_reduced<-anno_reduced %>% group_by(gname) %>% top_n(1, position)
anno_reduced<-unique(anno_reduced[c("CHROM","pos","lle.ratio","gname")])
#anno_red<-anno_reduced[!startsWith(anno_reduced$gname,"LOC"),]

anno_collapsed<-anno_reduced %>%
  group_by(CHROM,pos,lle.ratio) %>%
  summarize(gname = str_c(gname, collapse = ", ")) 

sum(onlyLLE$lle.ratio==0)

#anno_collapsed$gname[anno_collapsed$gname=="LOC105373062, PRR5-ARHGAP8, ARHGAP8"]<-"ARHGAP8"

#Merge and fill out unique missing values 
merged<-merge(anno_collapsed,top50_rs_numbers,all = T)
merged<- head(merged[order(merged$lle.ratio, decreasing = TRUE), ] ,200)
missing<-sum(is.na(merged$gname))
merged$gname[is.na(merged$gname)]<-paste0("missing",1:missing)

anno_collapsed<-anno_collapsed%>% group_by(gname) %>% top_n(1, pos)
hola<-anno_collapsed%>%
  arrange(desc(lle.ratio))  

toplot<-hola[1:25,]

toplot[c(1,3),]

tested<-merged[c(1,7,11,12,13,14,15,17,18,19),]

#tested$gname[tested$pos==161704797]<-"FCRLB"

#tested$gname[startsWith(tested$gname,"missing")]<-"*"
tested$gname[startsWith(tested$gname,"missing")]<-""
#window_results<-merged %>% group_by(win_start=pos-100000,win_end=pos+100000)
#regions<-window_results[c("CHROM","win_start","win_end")]
#temp<-merge(regions,merged,by="CHROM")
#temp2<-temp[(temp$pos>=temp$win_start) & (temp$pos<=temp$win_end),]
#top_SNP<-temp2%>%group_by(CHROM,win_start)%>%top_n(1,lle.ratio)
#top_SNP2<-top_SNP%>%group_by(CHROM,win_start)%>%top_n(1,pos)
#top_SNP3<-unique(top_SNP2[c("CHROM","pos","lle.ratio","rs","gname")])
#top_SNP3$gname[is.na(top_SNP3$gname)]<-paste0("missing",1:60)

#merged<-top_SNP3 %>% group_by(gname) %>% top_n(1, pos)


write.table(merged,paste0("top100_annotated_JEJKORCHS_selscanv2.txt"),quote = F,append = F,row.names = F)


hola<-merge(merged,scan.haenyeo)%>%
  arrange(desc(lle.ratio)) 
write.table(hola,paste0("SNPS_allelefreq_3popsCHS_noflip.txt"),quote = F,append = F,row.names = F)


#colors=c('#008000','#FF0000','#7FCC12','#F0E442', '#0000FF','#FFA500','#8B2121')
colors<-c("#00b300", '#8B2121','#008000',"#56B4E9", 
          '#0000FF','#FF0000','#FFA500')

noloc<-anno_collapsed[!startsWith(anno_collapsed$gname,"LOC"),]
half_data<-scan.haenyeo[scan.haenyeo$lle.ratio>quantile(scan.haenyeo$lle.ratio,.5),]
half_anno<-merge(half_data, tested,all.x = T)

newdata <- half_anno[order(parse_number(half_anno$CHROM)),] 

n_of_chr<-length(unique(newdata$CHROM))
newdata$SNPs<-1:nrow(newdata)
axis_set <- newdata %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(SNPs))


# Degrees of freedom (typically 1 for LRT comparing two models)
degrees_of_freedom <- 1

# Calculate p-values for each LRT value
newdata$p_values <- 1 - pchisq(newdata$lle.ratio, df = 1 * degrees_of_freedom)

newdata$logpval<- (-log10(newdata$p_values))

cbPaletteP5 <- c( "#DBD2C8","#FD767D","#74491F","#314F70","#a0c5d3",
                  "#a5021e","#6c6e00","#F4C250", "#f29900","#3F2B13","#463A68")


ggplot(newdata,aes(x=SNPs,y=lle.ratio,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1)+
  #ggtitle(paste0("Selection on Jeju component"))+
  theme_classic(base_size=14)+
  ylim(c(0,NA))  +scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#314F70","#a0c5d3"), unique(length(axis_set$CHROM))))+
  geom_text_repel(size=3,max.overlaps = 60,fontface = 'bold',nudge_y = 6,nudge_x = -3000) +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))+
  geom_point(data=newdata%>% filter(!is.na(gname)),
             pch=21,
             size=4) +ylab("Likelihood ratio")+
  geom_hline(yintercept=33, linetype="dashed", color = "red")
  #ylim( c(0, max(alldata$pbs1) ))
ggsave(filename=paste0("JEJKORCHS_onlyTestedLabeled_thresholdLLE.png"),width=10,height=3)

ggplot(newdata,aes(x=SNPs,y=logpval,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1)+theme_classic(base_size=14)+
 scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#314F70","#a0c5d3"), unique(length(axis_set$CHROM))))+
  geom_text_repel(size=3,max.overlaps = 60,
                  fontface = 'bold',nudge_y = 1,nudge_x = -2500) +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))+
  geom_point(data=newdata%>% filter(!is.na(gname)),
             pch=21,
             size=4) +ylab(bquote(~-log[10](pvalue)))+
  geom_hline(yintercept=8, linetype="dashed", color = "red")
ggsave(filename=paste0("JEJKORCHS_threshold_log10pval.png"),width=10,height=3)


ggsave(filename=paste0("JEJKORCHS_onlyTestedLabeled_nobox.pdf"),width=10,height=3)

df <- top50_rs_numbers %>%
  group_by(CHROM) %>%
  mutate(diff = ifelse(CHROM == lag(CHROM) & lle.ratio != lag(lle.ratio), abs(pos- lag(pos)), NA)) %>%
  ungroup()


onlychr8<-scan.haenyeo[scan.haenyeo$CHROM=="chr8",]
chr8Top<-results_list[[8]]


anno8<-merge(chr8Top,annotation_grch38,by="CHROM",all.x = T)
annotated8<-anno8[(anno8$pos>=anno8$start) & (anno8$pos<=anno8$end),]


anno_reduced<-annotated8 %>% group_by(gname) %>% top_n(1, lle.ratio)
anno_reduced<-unique(anno_reduced[c("CHROM","pos","lle.ratio","gname")])
anno_collapsed<-anno_reduced %>%
  group_by(CHROM,pos,lle.ratio) %>%
  summarize(gname = str_c(gname, collapse = ", ")) 

top_rsnumbers<- head(anno_collapsed[order(anno_collapsed$lle.ratio, decreasing = TRUE), ] ,20)

onlychr8_2plot<-merge(onlychr8,top_rsnumbers,all.x = T)

ggplot(onlychr8_2plot,aes(x=pos,y=lle.ratio,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1.5)+
  ggtitle(paste0("Selection on Jeju component"))+theme_classic(base_size=14)+
  ylim(c(0,NA))  +
  scale_color_manual(values = rep(c(colors[k+1], "gray"), unique(length(axis_set$CHROM))))+
  geom_label_repel(size=3,max.overlaps = 60,fontface = 'bold') +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+xlim(4e7,5e7)
ggsave(filename=paste0("chr8_selscan_outCHS_anno_region.png"),width=10,height=4)





merged<-merge(anno_collapsed,chr8Top,all = T)
adios<-merge(merged,scan.haenyeo)%>%
  arrange(desc(lle.ratio)) 
write.table(adios,paste0("chr8_SNPS_allelefreq_3popsCHS_noflip.txt"),quote = F,append = F,row.names = F)


library(reshape2)
alleleFreqOutliersCHSchr8.frq <- read.csv("~/Documents/Berkeley/haenyeo/Analysis2023/Selscan/selscan3pop/outCHS/alleleFreqOutliersCHSchr8.frq.strat", sep="")
freqs<-dcast(alleleFreqOutliersCHSchr8.frq,SNP~CLST,value.var="MAF")


complettable<-merge(adios, freqs,by.x = "rs",by.y="SNP")%>%
  arrange(desc(lle.ratio)) 

write.table(complettable,paste0("chr8_SNPS_allelefreq_3popsCHS_noflip_perPOP.tsv"),quote = F,append = F,row.names = F,sep = "\t")



bye<-newdata[!is.na(newdata$gname),c("rs","pos","lle.ratio","chr","gname","p_values","logpval")]
write.table(bye,"labeledSNPS_pval.txt",quote = F,append = F,row.names = F,sep = "\t")


###########Inflation factor########################################
onlyLLE <- read.csv("onlyLLE", sep="")
simpleQQPlot <- function(observedPValues, min = 0, max = 1, color, font_size = 1.5) {
  x = -log10(sort(runif(n = length(observedPValues), min = min, max = max)))
  y = -log10(sort(observedPValues))
  all1 = c(x, y)
  range = c(min(all1), max(all1))
  
  # Increase font size for labels and axes
  plot(x, y,
       xlab = "-log10(expectedPValues)",
       ylab = "-log10(observedPValues)",
       xlim = range, ylim = range, col = color,
       cex.lab = font_size,   # Font size for labels
       cex.axis = font_size   # Font size for axis ticks
  )
  
  abline(0, 1, col = "red")  # Red line at y = x
}

ggQQPlot <- function(observedPValues, min = 0, max = 1, color = "blue", font_size = 1.5) {
  # Generate expected p-values (from uniform distribution)
  expectedPValues <- -log10(sort(runif(n = length(observedPValues), min = min, max = max)))
  
  # Calculate the observed p-values
  observedPValues <- -log10(sort(observedPValues))
  
  # Create a data frame for ggplot
  data <- data.frame(
    expected = expectedPValues,
    observed = observedPValues
  )
  
  # Create the QQ plot using ggplot
  ggplot(data, aes(x = expected, y = observed)) +
    geom_point(color = color) +   # Scatter plot with specified color
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Reference line y = x
    xlab("-log10(expected P-values)") +
    ylab("-log10(observed P-values)") +
    theme_minimal() +
    theme(
      text = element_text(size = font_size * 5),  # Adjust font size based on input
      axis.title = element_text(size = font_size * 5)
    ) +
    coord_equal()  # Ensure equal scaling for x and y axes
}
#Get 50% of the values
half<-as.data.frame(apply(onlyLLE, 2, FUN = function(x) sort(x, decreasing = T)[1:floor(length(x)/2)]))
half$p_values <- 1 - pchisq(half$lle.ratio, df = 1 )

png("qqplot_pval_halfDATA.png", width=2000, height=2000,res=300)
simpleQQPlot(half$p_values,min=0,max=1,"black")
dev.off()


ggQQPlot(half$p_values,min=0,max=1,"black")
ggsave(filename="QQPlot_BeforeCorrection.png",width=10,height=10)

onlyLLE<-half
# Step 1: Convert p-values to chi-square statistics
onlyLLE$chi_square_stats <- qchisq(1 - onlyLLE$p_values, df = 1)
# Step 2: Calculate the median of the chi-square statistics
median_chisq <- quantile(onlyLLE$chi_square_stats,.99)
# Step 3: Calculate the inflation factor
expected_median <- qchisq(0.99, df = 1)  # Expected median for chi-square with 1 df is ~0.455
inflation_factor <- median_chisq / expected_median
# Step 4: Adjust the chi-square statistics using the inflation factor
onlyLLE$adjusted_chi_square <- onlyLLE$chi_square_stats / inflation_factor
# Step 5: Convert the adjusted chi-square statistics back to p-values
onlyLLE$adjusted_p_values <- 1 - pchisq(onlyLLE$adjusted_chi_square, df = 1)
onlyLLE$adjusted_p_values[onlyLLE$adjusted_p_values == 0] <- 1.110223e-16

png("qqplot_pval_halfCorrected99th.png", width=2000, height=2000,res=300)
simpleQQPlot(onlyLLE$adjusted_p_values ,min=0,max=1,"black")
dev.off()




newdata$chi_square <- qchisq(1-newdata$p_values, df = 1)
newdata$chi_adjust<-newdata$chi_square/inflation_factor
newdata$pval_adjust <- 1-pchisq(newdata$chi_adjust, df = 1)
newdata$logpvalA<- (-log10(newdata$pval_adjust))
ggplot(newdata,aes(x=SNPs,y=logpvalA,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1)+theme_classic(base_size=14)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#314F70","#a0c5d3"), unique(length(axis_set$CHROM))))+
  geom_text_repel(size=3,max.overlaps = 60,
                  fontface = 'bold',nudge_y = 1,nudge_x = -2500) +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))+
  geom_point(data=newdata%>% filter(!is.na(gname)),
             pch=21,
             size=4) +ylab(bquote(~-log[10](pvalue)))


ggsave(filename=paste0("JEJKORCHS_onlyTestedLabeled_Pval10Corrected_IF99th.png"),width=10,height=3)

tested$pval <- 1-pchisq(tested$lle.ratio, df = 1)
tested$chi_square <- qchisq(1-tested$pval, df = 1)
tested$chi_adjust<-tested$chi_square/inflation_factor
tested$pval_adjust <- 1-pchisq(tested$chi_adjust, df = 1)

newdata$chi_square <- qchisq(1-newdata$p_values, df = 1)
newdata$chi_adjust<-newdata$chi_square/inflation_factor
newdata$pval_adjust <- 1-pchisq(newdata$chi_adjust, df = 1)
newdata$logpvalA<- (-log10(newdata$pval_adjust))
ggplot(newdata,aes(x=SNPs,y=chi_adjust,label=gname,color=as.factor(parse_number(CHROM))))+
  geom_point(size=1)+theme_classic(base_size=14)+
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#314F70","#a0c5d3"), unique(length(axis_set$CHROM))))+
  geom_text_repel(size=3,max.overlaps = 60,
                  fontface = 'bold',nudge_y = 1,nudge_x = -2500) +
  xlab("Chromosome")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))+
  geom_point(data=newdata%>% filter(!is.na(gname)),
             pch=21,
             size=4) +ylab("Chi-Squared Adjusted Values")

ggsave(filename=paste0("JEJKORCHS_onlyTestedLabeled_Chi2Corrected.png"),width=10,height=3)
