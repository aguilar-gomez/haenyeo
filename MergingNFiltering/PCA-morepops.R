library(ggplot2)
library(ggrepel)
library(gplots)


eigvectors <- read.table("name.eigenvec")
eigenval <- read.table("name.eigenval")


vecs<-eigvectors[c(3:22)]
haenyeo_meta <- read.delim("haenyeo_meta.tsv")
bamnames<-eigvectors$V2[c(1:75)]
Number<-as.numeric(sapply(strsplit(bamnames,"\\_"), `[`, 1))
samples<-data.frame(bamnames,Number)

#1000 genomes
Female_Phase3_list <- read.delim2("Female_Phase3_list", header=FALSE)
meta_out<-Female_Phase3_list[,c(1,4,5,6)]
colnames(meta_out)<-c("bamnames","pop","group","continent")

AsianFemale_Phase3_list <- read.delim2("AsianFemale_Phase3_list", header=FALSE)
meta_out2<-AsianFemale_Phase3_list[,c(1,4,5,6)]
colnames(meta_out2)<-c("bamnames","pop","group","continent")

CHB_Phase3_list <- read.delim2("CHBlist", header=FALSE)
meta_out3<-CHB_Phase3_list[,c(1,4,5,6)]
colnames(meta_out3)<-c("bamnames","pop","group","continent")


meta<-rbind(meta_out,meta_out2,meta_out3)
metadata<-merge(samples,haenyeo_meta)
colnames(eigvectors)<-c("sample","bamnames",paste0("PC",c(1:20)))

alldata<-merge(eigvectors,metadata,all.x = T)

finaldata<-merge(alldata,meta,all.x = T)

# Eigenvalues
eigvalue <- eigenval$V1/sum(eigenval$V1);
cat(signif(eigvalue , digits=3)*100,"\n");

# Plot

col1 <- c(finaldata$Group[c(1:75)],finaldata$group[c(76:223)],
          rep("KOREA1K",75),finaldata$group[c(299:459)])
alldata$Population <- ifelse(col1==0, "Jeju",
               ifelse(col1==1, "Haenyeo",
                      ifelse(col1==2, "Seoul",
                             ifelse(col1=="KOREA1K", "Korea1K",
                                    ifelse(col1=="Kinh Vietnamese", "Vietnamese",
                                           ifelse(col1=="Kinh,Kinh Vietnamese", "Vietnamese",
                                                  ifelse(col1=="Southern Han Chinese", "Southern_Chinese",
                                                         ifelse(col1=="Han Chinese", "Chinese",
                                    col1))))))))


col1 <- c(finaldata$Group[c(1:75)],finaldata$pop[c(76:223)],
          rep("KOREA1K",75),finaldata$pop[c(299:459)])
alldata$Pop3 <- ifelse(col1==0, "JEJ",
                             ifelse(col1==1, "HAE",
                                    ifelse(col1==2, "SEO",
                                           ifelse(col1=="KOREA1K", "KOR",
                                col1))))

famfile<-data.frame(alldata2$Population,alldata2$bamnames,0,0,0,-9)
write.table(famfile,file="pcapops.allchr.maf5.fam", quote = F,row.names = F,col.names = F,sep = "\t")
famfilepop3<-data.frame(alldata2$Pop3,alldata2$bamnames,0,0,0,-9)
write.table(famfilepop3,file="threeletterpop.fam", quote = F,row.names = F,col.names = F,sep = "\t")

comp<-c(2,3)
title <- paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="",collapse="")
xlabel = paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)",sep="")
ylabel = paste("PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")
rownames(alldata)<-eigvectors$V2

cbPalette <- c( "blue","hotpink2","red","orange", "#56B4E9",  "#F0E442", "black","#00b300","gray","purple","limegreen","green")
col2hex(cbPalette)
forpong <- c("firebrick2","hotpink2", "blue","purple","#F0E442", "gray","blue","#56B4E9")
#[1] "#0000FF" "#FF0000" "#56B4E9" "#F0E442" "#000000" "#EE6AA7" "#A020F0" "#BEBEBE"
ggplot() + geom_point(data=alldata, aes_string(x=x_axis, y=y_axis, color="Population"),alpha = .7)+
  scale_colour_manual(values=cbPalette)+ theme_bw(base_size = 18) +xlab(xlabel)+ylab(ylabel)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename=paste0("PCA_","PC",comp[1],"PC",comp[2],".png"),width=7,height=4)

barplot(cumsum(eigvalue),xlab = "Number of Components",ylab="Variance (%)",main="PCA Haenyeo")


barplot(eigvalue,xlab = "Number of Components",ylab="Variance (%)",main="PCA Haenyeo")


#Plot with labels
comp<-c(1,2)
title <- paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="",collapse="")
xlabel = paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)",sep="")
ylabel = paste("PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")
rownames(alldata)<-eigvectors$V2
subsetvietnam<-alldata[alldata$Pop3=="KHV",]
ggplot() + geom_point(data=alldata, aes_string(x=x_axis, y=y_axis, color="Population"),size=3,alpha = .8)+
  scale_colour_manual(values=cbPalette)+ theme_bw(base_size = 18) +xlab(xlabel)+ylab(ylabel)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel(data=subsetvietnam, aes(x=PC1, y=PC2,label=rownames(subsetvietnam), color="Population"))


ggsave(filename=paste0("PCA_labels","PC",comp[1],"PC",comp[2],".png"),width=9,height=4)
