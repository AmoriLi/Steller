library(dplyr)
library(tidyr)
library(Seurat) #v4.3
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(Matrix)
library(stringr)
library(reshape2)
library(grDevices)
library(pheatmap)
library(cluster)
library(png)
library(reshape2)
library(purrr)
library(broom)
library(aplot)
library(data.table)
library(tibble)
library(cowplot)
library(VennDiagram)
library(colorspace)
library(plotly)
library(igraph)
source("~/NPC_project/BMK/pipeline/amplicon/function.R")
source("~/NPC_project/BMK/pipeline/spatial_plot.r")
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")
exp_col<-c('#034f84','#92a8d1','#d6d4e0','#f4a688',"#ED797B","#d64161","#c94c4c")

wkdir<-"~/NPC_project/BMK/tissue/brain/P5/Merged_2512/3_clone"
setwd(wkdir)
#wkdir<-"~/NPC_project/BMK/tissue/brain/P5/Merge/integrate"#args[1] 
if(!dir.exists(wkdir)){
    dir.create(wkdir,recursive = TRUE)
}

### read, cell, clonebc preprocess stat
df1<-read.csv("../../S3000/amplicon_all/tmp/rd.stat.txt",sep = '\t',header = FALSE)
#df1<-df1[c(1,2,3,4,8,9),]
colnames(df1)<-c('step','reads')
df1$step<-factor(df1$step,levels = unique(df1$step))
get_col_plot(df1,'step','reads')
ggsave('S3000.read.barplot.png',width = 5,height = 5,dpi = 400)

df<-read.csv("./S3000/amplicon_all/tmp/cell.stat.txt",sep = '\t',row.names = 1)
colnames(df)<-c('step','reads')
df$step<-factor(df$step,levels = unique(df$step))
get_col_plot(df,'step','reads')
ggsave('S3000.cell.barplot.png',width = 5,height = 5,dpi = 400)

df2<-read.csv("./S3000/amplicon_all/tmp/clonebc.stat.txt",sep = '\t',row.names = 1)
colnames(df2)<-c('step','reads')
df2$step<-factor(df2$step,levels = unique(df2$step))
get_col_plot(df2,'step','reads')
ggsave('S3000.clonebc.barplot.png',width = 5,height = 5,dpi = 400)

### merge two section amplicon and rna dataset
indir<-"~/NPC_project/BMK/tissue/brain/P5"
sample<-c("S3000","S3000_50cs3")
df_list<-list()
for(i in sample){
    df<-read.csv(file.path(indir,i,"amplicon_all/tmp/qc3_cell_clonebc_umi.tsv"),row.names =1,sep = '\t')[,c(1,2,3,4,5)]
    colnames(df)<-c('cell','clonebc','cell_x','cell_y','umi')
    df$cell<-paste(i,df$cell,sep = "_")
    df_list[[i]]<-df
}
n_distinct(df_list[['S3000']]$clonebc)
n_distinct(df_list[['S3000_50cs3']]$clonebc)
df<-bind_rows(df_list)
da<-create_object(df)
da$sample<-str_split(colnames(da),"_cell",simplify = TRUE)[,1]
table(da$sample)
da<-get_bc_qc_plot(da,saveplot = FALSE,outdir = file.path(wkdir,"2_qc"))
da<-da[,da$nCount_CloneBC<=73]
da<-da[rowSums(da@assays$CloneBC@counts)>0,]
table(da$sample)
saveRDS(da,file.path(wkdir,'clonebc.cell.rds'))

data<-data.frame(sample=c('S3000','S3000_50cs3'),
                num=c(2003,2589))
pdf(file.path(wkdir,'2603.CloneBC_diversity.barplot.comparison.pdf'),width = 3.2,height = 3)
options(repr.plot.width=3.2,repr.plot.height=3)
ggplot(data, aes(x = sample, y=num,fill=sample)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = get_specific_group_color('ident')) + 
        theme_minimal() + theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
        ylab(label = "CloneBC diversity")+
        labs(fill = "")
dev.off()

df<-da@meta.data
options(repr.plot.width = 3, repr.plot.height = 3)
p1 <- ggplot(df, aes(y = log2(nCount_CloneBC), x=sample,fill=sample)) + 
    geom_violin(color = "#e9ecef") + 
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +  
    #stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red", color = "black") +  
    scale_fill_manual(values = get_specific_group_color('ident')) + 
    theme_minimal() + theme(axis.text = element_text(size = 10), 
    axis.title = element_text(size = 15), legend.position = "none") + 
    labs(fill = "")
plot(p1)
#ggsave('2603.nCount_CloneBC.byIdent.Violin_boxplot.png',width = 3,height = 3,dpi = 400)

p1 <- ggplot(df, aes(x = nFeature_CloneBC)) + 
    geom_histogram(color = "#e9ecef",fill= "#ada397",
        binwidth = 1) + theme_minimal() + theme(axis.text = element_text(size = 10), 
    axis.title = element_text(size = 15), legend.position = "none") + 
    labs(title = 'Total',fill = "")
p2 <- ggplot(df %>% filter(sample=='S3000'), aes(x = nFeature_CloneBC)) + 
    geom_histogram(color = "#e9ecef", fill=get_specific_group_color('ident')['S3000'], #position = "identity", 
        binwidth = 1) + theme_minimal() + 
    coord_cartesian(ylim = c(0, 3100)) +  
    theme(axis.text = element_text(size = 10), 
    axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
    labs(title = 'S3000',fill = "")
p3 <- ggplot(df %>% filter(sample=='S3000_50cs3'), aes(x = nFeature_CloneBC)) + 
    geom_histogram(color = "#e9ecef", fill=get_specific_group_color('ident')['S3000_50cs3'], #position = "identity", 
        binwidth = 1) + 
    scale_fill_manual(values = get_specific_group_color('ident')) + 
    coord_cartesian(ylim = c(0, 3100)) + 
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
    axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
    labs(title = 'S3000_50cs3',fill = "")
pdf(file.path(wkdir,'2603.nFeature_CloneBC.Total_Ctrl_CS3.histogram.pdf'),width = 6,height = 4)
options(repr.plot.width = 6, repr.plot.height = 4)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
#ggsave('2603.nFeature_CloneBC.Total_Ctrl_CS3.histogram.png',width = 6,height = 3,dpi = 400)

summary(da@meta.data[,'nCount_CloneBC'])
summary(da@meta.data[da$sample=='S3000','nCount_CloneBC'])
summary(da@meta.data[da$sample=='S3000_50cs3','nCount_CloneBC'])
summary(da@meta.data[,'nFeature_CloneBC'])
summary(da@meta.data[da$sample=='S3000','nFeature_CloneBC'])
summary(da@meta.data[da$sample=='S3000_50cs3','nFeature_CloneBC'])

df2<-df[df$cell %in% colnames(da) & df$clonebc %in% rownames(da),]
write.csv(df2,file.path(wkdir,'cell.clonebc.umi_table.txt'))

### Merge with rna
da_rna<-readRDS(file.path(indir,"Merged_2512/1_cluster/object.integrated.spatial.clustered.rds"))
tcell1<-colnames(da_rna)[da_rna$orig.ident=='S3000']
tcell2<-colnames(da_rna)[da_rna$orig.ident=='S3000_50cs3']
apcell1<-colnames(da)[da$sample=='S3000']
apcell2<-colnames(da)[da$sample=='S3000_50cs3']
itcell1<-intersect(tcell1,apcell1)
itcell2<-intersect(tcell2,apcell2)
length(itcell1)/length(tcell1)
length(itcell2)/length(tcell2)
length(itcell1)
length(itcell2)
data<-data.frame(sample=c('S3000','S3000_50cs3'),
                percent=c(2.7,3.2))
pdf(file.path(wkdir,'2603.CloneBC_cell_percent.barplot.comparison.pdf'),width = 3,height = 3)
options(repr.plot.width=3,repr.plot.height=3)
ggplot(data, aes(x = sample, y=percent,fill=sample)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = col_ident) + 
        theme_minimal() + theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
        ylab(label = "CloneBC+ cells (%)")+
        labs(fill = "")
dev.off()

da_merge<-get_clone_rna_object(da,da_rna)
da_merge<-da_merge[,colnames(da_merge@assays$Spatial@counts)]

## clonebc and mcherry expression correlation
da_merge$mCherry<-ifelse(da_merge@assays$Spatial@counts["mcherry",]>0,"mCherry","none")
da_merge$CloneID<-ifelse(colnames(da_merge) %in% colnames(da_merge@assays$CloneBC),'CloneID','none')

options(repr.plot.width=10,repr.plot.height=3)
FeaturePlot(da_merge,features = c('mcherry','nCount_CloneBC','nFeature_CloneBC'),max.cutoff = 30,min.cutoff = 0.5,
            pt.size = 1,ncol = 3,order = TRUE)+
    coord_fixed()
ggsave('UMAP.mcherry.CloneBC_count_feature.png',width = 10,height = 3,dpi = 500)

sparse.gbm.T<-as(da@assays$CloneBC@counts,"dgTMatrix") # convert to coo_matrix
writeMM(sparse.gbm.T,file=file.path(wkdir,"clonebc.cell.matrix.mtx"))
#system("gzip 2107_450.matrix.mtx.gz")
write.table(rownames(da),file=file.path(wkdir,"clonebc.tsv"),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(da@meta.data,file=file.path(wkdir,"cloenbc_cell.meta.tsv"),sep='\t')

###Run in shell: python clonecalling.py $(wkdir)
clones<-read.csv(file.path(wkdir,"tmp/clones.csv"),header = TRUE)
rownames(clones)<-clones$cell
colnames(clones)<-c('clone.id','cell.barcode','edge.den')

clones$section<-gsub('_cell_[0-9]*','',clones$cell.barcode)
table(clones$section)
n_distinct(clones$clone.id[clones$section=='S3000'])
n_distinct(clones$clone.id[clones$section=='S3000_50cs3'])

data<-data.frame(sample=c('S3000','S3000_50cs3'),
                percent=c(1158,1302))
pdf(file.path(wkdir,'Clone_cell_number.barplot.comparison.pdf'),width = 3,height = 3)
options(repr.plot.width=3,repr.plot.height=3)
ggplot(data, aes(x = sample, y=percent,fill=sample)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = col_ident) + 
        theme_minimal() + theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
        ylab(label = "Clone+ cells (%)")+
        labs(fill = "")
dev.off()

data<-data.frame(sample=c('S3000','S3000_50cs3'),
                percent=c(582,659))
pdf(file.path(wkdir,'Clone_number.barplot.comparison.pdf'),width = 3,height = 3)
options(repr.plot.width=3,repr.plot.height=3)
ggplot(data, aes(x = sample, y=percent,fill=sample)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = col_ident) + 
        theme_minimal() + theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), legend.text = element_text(size = 10)) + 
        ylab(label = "#Clone")+
        labs(fill = "")
dev.off()

da_merge<-clone_metric_plot(da_merge,clones,saveplot = TRUE,outdir = wkdir)
saveRDS(da_merge,file.path(wkdir,"RNA.allcloneidcell.merged.rds"))

da<-da_merge[,colnames(da_merge@assays$CloneBC@counts)]
saveRDS(da,file.path(wkdir,"all.cloneidcell.rna.object.rds"))


