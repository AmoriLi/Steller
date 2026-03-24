library(Seurat)
library(tidyverse)
library(ggplot2)
library(DoubletFinder)
source("/public/home/zxli_gibh/NPC_project/BMK/pipeline/RNA/Helper.R")
#source("/public/home/zxli_gibh/NPC_project/BMK/pipeline/spatial_plot.r")

wd<-"~/NPC_project/SeekSpace/B29_left/section1/RNA/Steller"
if(!dir.exists(wd)){
    dir.create(wd,recursive = TRUE)
}
setwd(wd)
getwd()

da <- Read10X('../Outs/B29_left11_filtered_feature_bc_matrix')
da <- CreateSeuratObject(counts=da,project='B29_left1')

m.s.genes<-c('Exo1','Msh2','Mcm4','Chaf1b','Rrm2','Cenpu','Mrpl36','Gmnn','Hells','Cdc6','Gins2','Uhrf1',
             'Cdc45','Slbp','Ubr7','Fen1','Rad51ap1','Mcm5','Nasp','Cdca7','Blm','Usp1','Ung','Prim1',
             'Clspn','Mcm6','Dtl','Pola1','Dscc1','Tipin','Wdr76','Casp8ap2','Tyms','Ccne2','Rrm1','Polr1b',
             'Rfc2','Rad51','E2f8','Mcm7','Pcna')
#m.g2m.genes<-convertHumanGeneList(cc$g2m.genes)
m.g2m.genes<-c('Dlgap5','Ctcf','Smc4','Kif20b','Cdc25c','Gtse1','Tpx2','Hmgb2','Cks1brt','Cdca2',
               'Top2a','Cks2','Cdca3','G2e3','Ttk','Ncapd2','Lbr','Anp32e','Ckap2','Tacc3','Nek2',
               'Cenpe','Kif11','Anln','Hjurp','Aurkb','Rangap1','Cks1b','Hmmr','Ckap5','Cdc20',
               'Psrc1','Kif23','Ect2','Kif2c','Ndc80','Nuf2','Cdca8','Birc5','Cenpf','Ube2c','Jpt1',
               'Pimreg','Nusap1','Mki67','Tubb4b','Bub1','Cenpa','Ccnb2','Aurka','Ckap2l')
Hb_gene<-rownames(da)[grepl("^Hb",rownames(da))]

rps_genes<-rownames(da@assays$RNA@counts)[grep("^Rps",rownames(da@assays$RNA@counts))]

rpl_genes<-rownames(da@assays$RNA@counts)[grep("^Rpl",rownames(da@assays$RNA@counts))]

mt_genes<-rownames(da@assays$RNA@counts)[grep("^mt-",rownames(da@assays$RNA@counts))]
da<-PercentageFeatureSet(da,"^mt-",col.name = "percent_mt")
da<-PercentageFeatureSet(da,"^Rp[s,l]",col.name = "percent_rp")
da<-PercentageFeatureSet(da,"^Hb",col.name = "percent_hbb")

da <- NormalizeData(da, normalization.method = "LogNormalize", scale.factor = 10000) %>%
FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
ScaleData() %>%
RunPCA() %>%
FindNeighbors(dims = 1:15) %>%
FindClusters(resolution = seq(0.2, 1.0, 0.2)) %>%
RunUMAP(dims = 1:15)

spatial_df <- read.table('../Outs/B29_left11_filtered_feature_bc_matrix/cell_locations.tsv.gz', row.names = 1, sep = '\t',header = T)
colnames(spatial_df) <- c("spatial_1","spatial_2")
spatial_matrix <- as.matrix(spatial_df)
spatial_matrix_sorted <- spatial_matrix[match(row.names(da@meta.data),row.names(spatial_matrix)), ]
da@reductions$spatial <- CreateDimReducObject(embeddings = spatial_matrix_sorted, key='spatial_', assay='RNA')

samplename = 'B29_left1'
size_x = 55050
size_y = 19906
da@misc$info[[`samplename`]]$size_x = as.integer(size_x)
da@misc$info[[`samplename`]]$size_y = as.integer(size_y)

img = '../Outs/B29_left11_aligned_DAPI.png'  

#base64
img_64 = base64enc::dataURI(file = img)
da@misc$info[[`samplename`]]$img = img_64
#png
img_gg <- png::readPNG(img)
img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1, "npc"))
da@misc$info[[`samplename`]]$img_gg = img_grob

img = '../Outs/B29_left11_aligned_HE_TIMG.png'  

#base64
img_64 = base64enc::dataURI(file = img)
da@misc$info[[`samplename`]]$img_he = img_64
#png
img_gg <- png::readPNG(img)
img_grob <- grid::rasterGrob(img_gg, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1, "npc"))
da@misc$info[[`samplename`]]$img_he_gg = img_grob

saveRDS(da,"B29_left1_raw.rds")

##gender
sum(da@assays$RNA@counts['Sry',])

da$mcherry<-da@assays$RNA@counts['mcherry',]
options(repr.plot.height=5, repr.plot.width=5)
FeaturePlot(da, reduction = 'umap', 
            features=c('mcherry'), 
            order = TRUE,
            pt.size = 0.1)
table(da@assays$RNA@counts['mcherry',]>0)

#dimplot umap
options(repr.plot.height=7, repr.plot.width=7)
DimPlot(da, reduction = 'umap',group.by = "RNA_snn_res.0.4",label = TRUE)+
    scale_color_manual(values = get_specific_group_color("div_col"))

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res <- paramSweep(da, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- da$RNA_snn_res.0.2
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(da@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
da <- doubletFinder(da, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#da <- doubletFinder(da, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_276", sct = FALSE)
head(da@meta.data)
dbnm<-grep("DF.",colnames(da@meta.data),value = TRUE)
table(da@meta.data[,dbnm])

options(repr.plot.width=5,repr.plot.height=5)
DimPlot(da,reduction="umap",group.by = dbnm)

da_s<-da[,(da@meta.data[,dbnm]=="Singlet")]
da_s

conf_cells<-rownames(da_s@meta.data %>% filter(nCount_RNA>800, nCount_RNA<20000, nFeature_RNA>600, nFeature_RNA<8000, 
                                              percent_mt<10, percent_rp<2))
head(conf_cells)
length(conf_cells)

da_s$good_quality<-ifelse(colnames(da_s) %in% conf_cells,TRUE,FALSE)
table(da_s$good_quality)

tmp<-da_s@meta.data %>% group_by(RNA_snn_res.0.4) %>% 
    summarise(count = sum(!good_quality),
    proportion = count / n())
tmp

#remove clusters with >70% low quality cells
da_s<-da_s[,!da_s$RNA_snn_res.0.4 %in% c(0,8)]
da_s

table(da_s$mcherry>0)

da_s <- NormalizeData(da_s, normalization.method = "LogNormalize", scale.factor = 10000) %>%
FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
ScaleData() %>%
RunPCA() %>%
FindNeighbors(dims = 1:20) %>%
FindClusters(resolution = seq(0.2, 1.0, 0.2)) %>%
RunUMAP(dims = 1:20)

#dimplot umap
options(repr.plot.height=7, repr.plot.width=7)
DimPlot(da_s, reduction = 'umap',group.by = "RNA_snn_res.0.2",label = TRUE)+
    scale_color_manual(values = get_specific_group_color("div_col"))

ct_map<-c("EN","IN","IN","Astro","IN_Lhx6","OPC","EN","EN","Olig","EN_TH","VLMC","Olig","Immune","CA_DG","EN")
names(ct_map)<-c(0:14)
ct_map

da_s2<-da_s[,da_s$RNA_snn_res.0.2!=15]
da_s2$celltype<-as.character(ct_map[match(da_s2$RNA_snn_res.0.2,names(ct_map))])
table(da_s2$celltype)

ct_col<-get_specific_group_color("div_col")[1:length(unique(da_s2$celltype))]
names(ct_col)<-unique(da_s2$celltype)
#dimplot umap
options(repr.plot.height=7, repr.plot.width=7)
DimPlot(da_s2, reduction = 'umap',group.by = "celltype",label = TRUE,cols=ct_col)

# DAPI
options(repr.plot.height=5, repr.plot.width=15)
ImageSpacePlot(obj=da_s2, group_by = "celltype",type="DAPI",size=0.01)+
    scale_color_manual(values = ct_col)

options(repr.plot.height=5, repr.plot.width=15)
ImageSpacePlot(obj=da_s2[,da_s2$mcherry>0], group_by = "celltype",type="DAPI",size=0.5)+
    scale_color_manual(values = ct_col)

#dimplot umap
da_s2$mcherry_bi<-ifelse(da_s2$mcherry>0,TRUE,FALSE)
table(da_s2$mcherry_bi)
options(repr.plot.height=7, repr.plot.width=7)
DimPlot(da_s2, reduction = 'umap',group.by = "mcherry_bi",order=TRUE,cols=c("grey","#c1502e"),pt.size = 0.1)

options(repr.plot.height=5, repr.plot.width=15)
ImageSpacePlot(obj=da_s2, group_by = "celltype",type="HE",size=0.01)+
    scale_color_manual(values = ct_col)

options(repr.plot.height=5, repr.plot.width=10)
FeaturePlot(obj=da_s, feature=c("Meg3","Rbfox3"),n=2)

options(repr.plot.height=5, repr.plot.width=10)
FeaturePlot(obj=da_s, feature=c("C1qb","Cx3cr1"),order=TRUE,n=2)

options(repr.plot.height=5, repr.plot.width=10)
FeaturePlot(obj=da_s, feature=c("Gfap","Aldoc"),order=TRUE,
            n=2)

options(repr.plot.height=5, repr.plot.width=15)
FeaturePlot(obj=da_s, feature=c("Plp1","Mobp","Pdgfra"),n=3)

options(repr.plot.height=5, repr.plot.width=5)
FeaturePlot(obj=da_s, feature=c("Vtn"),order=TRUE,n=1)

options(repr.plot.height=5, repr.plot.width=5)
FeaturePlot(obj=da_s, feature=c("Prox1"),#order=TRUE,
            n=1)

options(repr.plot.height=5, repr.plot.width=20)
FeaturePlot(obj=da_s, feature=c("Slc17a7","Slc17a6","Satb2","Fezf2"),n=4)

options(repr.plot.height=5, repr.plot.width=15)
FeaturePlot(obj=da_s, feature=c("Lhx6","Gad1","Slc32a1"),n=3)

saveRDS(da_s2,"B29_left_slide1_filtered_celltype.rds")
write.table(colnames(da_s2),"clean.cellbc.txt",quote = FALSE,row.names = FALSE,col.names=FALSE)
