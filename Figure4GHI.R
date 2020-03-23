library(reshape2)
library(plyr)
library(dplyr)

#import data
SZmix_stats_human <- read.table("SZ-RG_S7_L001_human_out_gene_exon_tagged.dge.summary.txt",header=T,stringsAsFactors = F)
SZmix_stats_mouse <- read.table("SZ-RG_S7_L001_mouse_out_gene_exon_tagged.dge.summary.txt",header=T,stringsAsFactors = F)

SZmix_data_human <- read.table("SZ-RG_S7_L001_human_out_gene_exon_tagged.dge.txt",header=T,stringsAsFactors = F)
SZmix_data_mouse <- read.table("SZ-RG_S7_L001_human_out_gene_exon_tagged.dge.txt",header=T,stringsAsFactors = F)

#colnames(SZmix_data_human) <- c('gene_id','CB2','CB1','CB4','CB3')
colnames(SZmix_data_human) <- c('gene_id','G3cell','G5cell','R1cell','R3cell')
#colnames(SZmix_data_mouse) <- c('gene_id','CB2','CB1','CB4','CB3')
colnames(SZmix_data_mouse) <- c('gene_id','G3cell','G5cell','R1cell','R3cell')
#SZmix_data_mix <- cbind(SZmix_data_human[ ,c('gene_id','G3cell','G5cell')],SZmix_data_mouse[ ,c('R1cell','R3cell')])
SZmix_data_merge<-merge(SZmix_data_human,SZmix_data_mouse,by="gene_id")
SZmix_data_submix <- SZmix_data_merge[ ,c(1:3,8,9)]
colnames(SZmix_data_submix) <- c('gene_id','G3cell','G5cell','R1cell','R3cell')

#melt so can do some basic calculations
SZmix_data_submix_melt <- melt(SZmix_data_submix, id.vars='gene_id')
SZmix_data_submix_sum<- ddply(SZmix_data_submix_melt, c('gene_id'), summarise,
                   sum = sum(value), sd = sd(value),
                   sem = sd(value)/sqrt(length(value)))
#Add the column of sum and sd to the original gene matrix table
SZmix_data_submix_merge <- merge(SZmix_data_submix,SZmix_data_submix_sum,by='gene_id')
#Subset data based on if sum>50 and sem>50, aim for sub to be ~100 genes
SZmix_sub <- subset(SZmix_data_submix_merge, c(sem > 30))#320 genes

SZmix_sub_sub <-SZmix_sub[ ,1:5]
SZmix_sub_melt <- melt(SZmix_sub_sub, id.vars='gene_id')
SZmix_sub_melt$value <- as.numeric(as.character(SZmix_sub_melt$value))

#make the heatmap
library(ggplot2)
library(viridis)
b <- c("R1cell"="firebrick1","R3cell"="firebrick3","G3cell"="chartreuse3","G5cell"="forestgreen")
png(filename='SZ_variablegenes_hmap.png', width=800, height=1400)
ggplot(SZmix_sub_melt,aes(variable,gene_id,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_x_discrete(limits=c("R1cell","R3cell","G3cell","G5cell"),
                   labels = c("B16-1 cell","B16-3 cells","U87-3cells","U87-5 cells")) +
  scale_fill_viridis(name="Gene \nExpression",option ="D") +
  labs(x="Cell Type/Number", y="Gene") +
  theme(axis.title.x = element_text(colour="black",size = 36,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(axis.title.y = element_text(colour="black",size = 36,margin = margin(t=0,r=0,b=0,l=0))) +
  theme(legend.title = element_text(colour="black", size=20, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 18)) +
  theme(axis.text.y = element_text(colour="black", size = 10)) +
  theme(axis.text.x = element_text(size = 22,colour = b)) +
  theme(axis.title = element_text(colour="black", size = 22)) +
  theme(plot.margin = margin(1, 0.5, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 0.1)) +
  theme(legend.text = element_text(colour="black", size = 28)) +
  theme(legend.title = element_text(colour="black", size = 30)) 
graphics.off()


pdf(file="SZ_heatmap.pdf",width=18, height=18)
ggplot(SZmix_sub_melt,aes(variable,gene_id,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_x_discrete(limits=c("R1cell","R3cell","G3cell","G5cell"),
                   labels = c("B16-1 cell","B16-3 cells","U87-3cells","U87-5 cells")) +
  scale_fill_viridis(name="Gene \nExpression",option ="D") +
  labs(x="Cell Type/Number", y="Gene") +
  theme(axis.title.x = element_text(colour="black",size = 48,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(axis.title.y = element_text(colour="black",size = 48,margin = margin(t=0,r=0,b=0,l=0))) +
  theme(legend.title = element_text(colour="black", size=32, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 28)) +
  theme(axis.text.y = element_text(colour="black", size = 7)) +
  theme(axis.text.x = element_text(size = 38,colour = b)) +
  theme(axis.title = element_text(colour="black", size = 22)) +
  theme(plot.margin = margin(1, 0.5, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 0.1)) +
  theme(legend.text = element_text(colour="black", size = 38,vjust = 2)) +
  theme(legend.title = element_text(colour="black", size = 34)) +
  theme(legend.background = element_rect(size=8)) +
guides(fill = guide_colourbar(barwidth = 3, barheight = 20))
dev.off()


write.table(SZmix_data_submix,"SZ_geneMatrix.txt",sep="\t",row.names = F)
write.table(SZmix_sub,"SZ_geneMatrix_SUB.txt",sep="\t",row.names = F)

#make density plots
#Rdata_density <- subset(SZmix_sub_melt, variable %in% barcodes)
png(filename='SZ_densityplots.png', width=800, height=600)
ggplot(SZmix_sub_melt,aes(value,colour=variable,fill=variable))+
  labs(x="Gene Expression (counts)", y="Density") +
  scale_fill_manual(values = c("forestgreen","chartreuse2","firebrick4","firebrick1"),
                      labels = c(" U87-5 cells", " U87-3 cells", " B16-3 cells"," B16-1 cell")) +
  geom_density(colour=NA,alpha=0.6) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 26)) +
  theme(legend.title.align = 0.8) +
  theme(legend.position = c(0.8, 0.85)) +
  theme(axis.text.y = element_text(colour="black", size = 28,angle = 90,vjust=0, hjust=0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 32)) +
  theme(axis.title.y = element_text(colour="black",size = 36,margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.title.x = element_text(colour="black",size = 36,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(plot.margin = margin(2, 1, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +  
  scale_size(range = c(2,12)) +
  theme(panel.grid.major = element_line(colour="grey", size=0.8, linetype="dashed")) +
  theme(panel.grid.minor = element_line(colour="grey", size=0.8, linetype="dashed")) +
  guides(fill = guide_legend(override.aes = list(size=10),reverse=T)) 
graphics.off()

pdf(file="SZ_density.pdf",width=14, height=10)
ggplot(SZmix_sub_melt,aes(value,colour=variable,fill=variable))+
  labs(x="Gene Expression (counts)", y="Density") +
  scale_fill_manual(values = c("forestgreen","chartreuse2","firebrick4","firebrick1"),
                    labels = c(" U87-5 cells", " U87-3cells", " B16-3 cells"," B16-1 cell")) +
  geom_density(colour=NA,alpha=0.6) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 30)) +
  theme(legend.title.align = 0.8) +
  theme(legend.position = c(0.8, 0.85)) +
  theme(axis.text.y = element_text(colour="black", size = 32,angle = 90,vjust=0, hjust=0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 34)) +
  theme(axis.title.y = element_text(colour="black",size = 38,margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.title.x = element_text(colour="black",size = 38,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(plot.margin = margin(2, 1, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +  
  scale_size(range = c(2,12)) +
  theme(panel.grid.major = element_line(colour="grey", size=0.8, linetype="dashed")) +
  theme(panel.grid.minor = element_line(colour="grey", size=0.8, linetype="dashed")) +
  guides(fill = guide_legend(override.aes = list(size=10),reverse=T)) 
dev.off()

#make the genes vs read plot
setwd('/Users/erica/Desktop/SZ_RNAseq')
stats<- c()
for (x in list.files(pattern="*.dge.summary.txt")) {
  u<-read.table(x, header=T,stringsAsFactors = F)
  u$Label = factor(x)
  stats <- rbind(stats, u)
  cat(x, "\n ")
}

library(reshape2)
library(plyr)
library(dplyr)
library(data.table)
slim <-t(data.frame(strsplit(as.character(stats$Label),"_",fixed = T)))
SZ_stats <- cbind(stats[ ,1:4],slim[ ,c(1,4)])

library(mgsub)
SZ_stats_names <- cbind(SZ_stats,mgsub(SZ_stats$CELL_BARCODE,c("GCTNGTCGTTTG","GCANGAACCACA","AGANATGAGGGT","ACGNACCTGTTC"),
                                       c("R3cell","R1cell","G5cell","G3cell")))
SZ_stats_names$Label <- paste(SZ_stats_names$'1',SZ_stats_names$'2')
SZ_stats_names <- SZ_stats_names[c(9,10,15,16) ,c(2:4,7,8)]
colnames(SZ_stats_names)<- c("Gene_reads","Trans_num","Gene_number","Cell_number","Exp")
write.table(SZ_stats_names,"SZstats.txt")

library(ggplot2)
pdf(file="RG_geneVreads_trans.pdf",width=14, height=10)
ggplot(data=SZ_stats_names) + geom_point(aes(x=Gene_reads,y=Gene_number,colour=Cell_number,size=8)) +
  guides(size=F) +
  scale_colour_manual("Sample",values = c("forestgreen","chartreuse2","firebrick4","firebrick1"),
                      labels = c(" U87-5 cells", " U87-3cells", " B16-3 cells"," B16-1 cell")) +
  ylab("Number of Genes Detected") + xlab("Number of Reads") +
  guides(size=F) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 30)) +
  theme(legend.title.align = 0.8) +
  theme(legend.position = c(0.15, 0.85)) +
  theme(axis.text.y = element_text(colour="black", size = 32,angle = 90,vjust=0, hjust=0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 34)) +
  theme(axis.title.y = element_text(colour="black",size = 38,margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.title.x = element_text(colour="black",size = 38,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(plot.margin = margin(1, 2, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +  
  scale_size(range = c(2,12)) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = T)) +
  theme(panel.grid.major = element_line(colour="grey", size=0.8, linetype="dashed")) +
  theme(panel.grid.minor = element_line(colour="grey", size=0.8, linetype="dashed"))
dev.off()

library(umap)
#separate data columns from labels
SZmix_data_submix_mat <- SZmix_data_submix[ ,2:5]
SZmix_data_submix_row.names <- SZmix_data_submix$gene_id


SZmix_data_submix_mat <- SZmix_data_submix[ ,2:5]
SZmix_data_submix_row.names <- SZmix_data_submix$gene_id
SZmix_data_submix_mat_t <- t(SZmix_data_submix_mat)
SZmix_umap <- umap(SZmix_data_submix_mat_t,n_epochs=500,n_neighbors=4)

#futz with data to prepare for plotting
SZmix_umap_lay <- as.data.frame(SZmix_umap$layout) 
colnames(SZmix_umap_lay)<- c('Dim1','Dim2')
SZmix_umap_lay$sample <- rownames(SZmix_umap_lay)

library(ggplot2)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Paired")
#save as pdf
pdf(file="SZ_umap.pdf",width=14, height=10)
ggplot(SZmix_umap_lay, aes(Dim1,Dim2,colour=sample)) +
  geom_point(aes(size=8)) +
  scale_colour_manual("Sample",values = c("forestgreen","chartreuse2","firebrick4","firebrick1"),
                      labels = c(" U87-5 cells", " U87-3cells", " B16-3 cells"," B16-1 cell")) +
  guides(size=F) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 30)) +
  theme(legend.title.align = 0.8) +
  theme(legend.position = c(0.15, 0.85)) +
  theme(axis.text.y = element_text(colour="black", size = 32,angle = 90,vjust=0, hjust=0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 34)) +
  theme(axis.title.y = element_text(colour="black",size = 38,margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.title.x = element_text(colour="black",size = 38,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +  
  scale_size(range = c(2,12)) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = T)) +
  theme(panel.grid.major = element_line(colour="grey", size=0.8, linetype="dashed")) +
  theme(panel.grid.minor = element_line(colour="grey", size=0.8, linetype="dashed"))
dev.off()

library(Rtsne)
SZ_normal <- normalize_input(SZmix_data_submix_mat_t)
colMeans(SZ_normal)
range(SZ_normal)

#Barnes-Hut t-Distributed Stochastic Neighbor Embedding
# Set a seed if you want reproducible results
set.seed(42)
tsne_out <- Rtsne(SZ_normal,pca=FALSE,perplexity=1,theta=0.0,pca_center=F,normalize=F) 
# Run TSNE
#Exact t-SNE can be computed by setting theta=0.0
SZ_tsne <-as.data.frame(tsne_out$Y) 
colnames(SZ_tsne) <- c('TSNE.1','TSNE.2')
SZ_tsne$sample <- c("U87-5 cells", "U87-3cells", "B16-3 cells","B16-1 cell")

pdf(file="SZ_tsne.pdf",width=14, height=10)
ggplot(SZ_tsne, aes(TSNE.1,TSNE.2,colour=sample)) +
  geom_point(aes(size=8)) +
  scale_colour_manual(values = c("forestgreen","chartreuse2","firebrick4","firebrick1"),
                      labels = c(" U87-5 cells", " U87-3 cells", " B16-3 cells"," B16-1 cell")) +
  scale_x_continuous(breaks=c(seq(-1200,1200,300))) + 
  scale_y_continuous(breaks=c(seq(-2000,2000,500))) +
  guides(size=F) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 30)) +
  theme(legend.title.align = 0.8) +
  theme(legend.position = c(0.15, 0.85)) +
  theme(axis.text.y = element_text(colour="black", size = 34,angle = 90,vjust=0, hjust=0.5)) +
  theme(axis.text.x = element_text(colour="black", size = 32)) +
  theme(axis.title.y = element_text(colour="black",size = 36,margin = margin(t=0,r=20,b=0,l=0))) +
  theme(axis.title.x = element_text(colour="black",size = 36,margin = margin(t=20,r=0,b=0,l=0))) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +  
  scale_size(range = c(2,12)) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse=T)) +
  theme(panel.grid.major = element_line(colour="grey", size=1, linetype="dashed")) +
  theme(panel.grid.minor = element_line(colour="grey", size=1, linetype="dashed"))
dev.off()
