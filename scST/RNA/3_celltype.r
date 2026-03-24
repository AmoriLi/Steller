library(Matrix)
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(patchwork)
library(dplyr)
source("/public/home/zxli_gibh/NPC_project/BMK/pipeline/RNA/Helper.R")
#source("/public/home/zxli_gibh/NPC_project/BMK/pipeline/spatial_plot.r")
source("/public/home/zxli_gibh/NPC_project/BMK/pipeline/amplicon/function.R")
exp_col<-c('#034f84','#92a8d1','#d6d4e0','#f4a688',"#ED797B","#d64161","#c94c4c")

###2025.11
#Test6a: Integrate Bandler_2022, Ratz_2022 and our B29-P10 SeekSpace as well as P5 scRNAseq
#define final level 1/2/3 celltypes

###2025-12-01
# add B29 s1 data into integration
#adjust L2 celltype annotation due to the weired glioblast/neuroblast distribution in BMK RCTD result
wd<-"~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/Merged_2512/3_clone"
setwd(wd)
ref1<-readRDS("~/NPC_project/SeekSpace/B29_left/section1/RNA/Steller/B29_left_slide1_filtered_celltype.rds")
ref2<-readRDS("~/NPC_project/SeekSpace/B29_left/section2/RNA/Steller/B29_left_slide2_filtered_celltype.rds")
da<-readRDS("all.cloneidcell.rna.object.rds")
rm_cells<-colnames(da@assays$CloneBC@counts)[!(colnames(da@assays$CloneBC@counts) %in% colnames(da@assays$Spatial@counts))]
da<-da[,!colnames(da) %in% rm_cells]
da[['RNA']]<-CreateAssayObject(counts = da@assays$Spatial@counts)
da2<-readRDS("~/NPC_project/BMK/tissue/brain/P5/Merged_2512/3_clone/all.cloneidcell.rna.object.rds")
rm_cells<-colnames(da2@assays$CloneBC@counts)[!(colnames(da2@assays$CloneBC@counts) %in% colnames(da2@assays$Spatial@counts))]
da2<-da2[,!colnames(da2) %in% rm_cells]
da2[['RNA']]<-CreateAssayObject(counts = da2@assays$Spatial@counts)

ref<-readRDS("~/NPC_project/scRNAseq/all_sc_integrate/1_celltype_integrated/test6a-v2.fn_ct3_tmp.only_snST.rds")
ref$lineage_tmp<-"RGC"
ref$lineage_tmp[ref$final_ct3_tmp %in% c('Macro_Igf1','Micro_Igf1','Micro_Itgam','Div_Micro_Itgam')]<-"Myeloid"
ref$lineage_tmp[ref$final_ct3_tmp %in% c('Vascular','Epithelia','Endothelia')]<-"Other"
ifnb.list <- list(ref,da,da2)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    DefaultAssay(x)<-'RNA'
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
#ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#    x <- ScaleData(x, features = features, verbose = FALSE)
#    x <- RunPCA(x, features = features, verbose = FALSE)
#})
#ifnb.list<-PrepSCTIntegration(ifnb.list)
anchors <- FindIntegrationAnchors(object.list = ifnb.list,anchor.features = features, reduction = "cca")
# this command creates an 'integrated' data assay
da.combined <- IntegrateData(anchorset = anchors)

# Run the standard workflow for visualization and clustering
da.combined <- ScaleData(da.combined, verbose = FALSE)
da.combined <- RunPCA(da.combined, npcs = 30, verbose = FALSE)
da.combined <- RunUMAP(da.combined, reduction = "pca", dims = 1:20)
da.combined <- FindNeighbors(da.combined, reduction = "pca", dims = 1:20)
da.combined <- FindClusters(da.combined, resolution = 0.2)
saveRDS(da.combined,"2602.test7.integration.object.rds")

da.combined<-FindClusters(da.combined,resolution = 0.25)
options(repr.plot.width=10,repr.plot.height=9)
DimPlot(da.combined,pt.size = 0.01,label.size = 5,cols=get_specific_group_color('div_col'),label = TRUE,repel = TRUE
       )

show_celltype_marker_umap(da.combined,major_ct = 'EN')
show_celltype_marker_umap(da.combined,major_ct = 'IN')
show_celltype_marker_umap(da.combined,major_ct = 'precursor')
show_celltype_marker_umap(da.combined,major_ct = 'glia')
show_celltype_marker_umap(da.combined,major_ct = 'other')
da.combined$integrate_ct1_tmp<-"Undefined"
da.combined$integrate_ct1_tmp[da.combined$seurat_clusters %in% c(0,4,6,7,10,12,13,18)]<-"EN"
da.combined$integrate_ct1_tmp[da.combined$seurat_clusters %in% c(2,3,8,14,16)]<-"IN"
da.combined$integrate_ct1_tmp[da.combined$seurat_clusters %in% c(10,15,17)]<-"Other"
da.combined$integrate_ct1_tmp[da.combined$seurat_clusters %in% c(1,5,9,11)]<-"Glia"
options(repr.plot.width=10,repr.plot.height=9)
DimPlot(da.combined,group.by = 'integrate_ct1_tmp',
        pt.size = 0.01,cols=get_specific_group_color('div_col')#,label = TRUE,repel = TRUE
       )

saveRDS(da.combined,"test7.Integrate_ct1_tmp.rds")

#### adjust image and cell position
da.combined<-readRDS("test7.Integrate_ct1_tmp.rds")
da.combined
image <- Read10X_Image(image.dir = "~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/RNA/BST/05.AllheStat/heAuto_level_matrix/subdata/L9_Spatial",image.name = "tissue_lowres_image.png")
DefaultAssay(image)<-'Spatial'
da.combined@images$B33_s1<-image
g <- rasterGrob(da.combined@images$B33_s1@image, interpolate = TRUE)
ggplot() +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void() # 

img<-da.combined@images$B33_s1@image
img_rot <- img
img_w <- ncol(img_rot)
img_h <- nrow(img_rot)
img_grob <- rasterGrob(img_rot, interpolate = TRUE)

meta_tmp<-da.combined@meta.data[da.combined$orig.ident=='B33_s1',]
dim(meta_tmp)
plot_data <- align_spots(meta_tmp, img_h = img_h,img_w = img_w,
                         scale_x = 0.825, scale_y=0.9,
                         off_x = 168, off_y = 35, rotate = 0,flip_y = TRUE)
head(plot_data)
ggplot(plot_data, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void() +
  labs(subtitle = "")

da.combined$x_px<-NA
da.combined$y_px<-NA
da.combined@meta.data[rownames(plot_data),c('x_px','y_px')]<-plot_data[,c('x_px','y_px')]
da.combined@images$B33_s1@image<-img_rot

#S3000
img1<-da.combined@images$S3000@image
img_rot1 <- aperm(img1, c(2, 1, 3))[,rev(1:nrow(img1)), ]
img_w1 <- ncol(img_rot1)
img_h1 <- nrow(img_rot1)
img_grob1 <- rasterGrob(img_rot1, interpolate = TRUE)
meta_tmp<-da.combined@meta.data[da.combined$orig.ident=='S3000',]
dim(meta_tmp)
plot_data1 <- align_spots(meta_tmp, img_h = img_h1,img_w = img_w1,
                         scale_x = 0.935, scale_y=0.83,
              ne           off_x = 65, off_y = 90, rotate = 270,flip_y = TRUE)
head(plot_data1)
ggplot(plot_data1, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob1, xmin = 0, xmax = img_w1, ymin = 0, ymax = img_h1) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w1), ylim = c(0, img_h1)) +
  theme_void() +
  labs(subtitle = "")

da.combined@meta.data[rownames(plot_data1),c('x_px','y_px')]<-plot_data1[,c('x_px','y_px')]
da.combined@images$S3000@image<-img_rot1

#S3000_50cs3
img2<-da.combined@images[['S3000_50cs3']]@image
img_rot2 <- aperm(img2, c(2, 1, 3))[rev(1:ncol(img2)),, ]
img_w2 <- ncol(img_rot2)
img_h2 <- nrow(img_rot2)
img_grob2 <- rasterGrob(img_rot2, interpolate = TRUE)
meta_tmp<-da.combined@meta.data[da.combined$orig.ident=='S3000_50cs3',]
dim(meta_tmp)
plot_data2 <- align_spots(meta_tmp, img_h = img_h2,img_w = img_w2,
                         scale_x = 1, scale_y=0.87,
                         off_x = 0, off_y = 6, rotate = 90,flip_y = TRUE)
head(plot_data2)
ggplot(plot_data2, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob2, xmin = 0, xmax = img_w2, ymin = 0, ymax = img_h2) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w2), ylim = c(0, img_h2)) +
  theme_void() +
  labs(subtitle = "")

da.combined@meta.data[rownames(plot_data2),c('x_px','y_px')]<-plot_data2[,c('x_px','y_px')]
da.combined@images$S3000_50cs3@image<-img_rot2

### Align A13 two sections
plot_data1$slice<-'S3000'
plot_data2$slice<-'S3000_50cs3'
# Manually adjust
shift_x <- 135    # Left right
shift_y <- 95  # up down
angle_deg <- -2  # rotate
scale_val <- 0.95    # scale
angle_rad <- angle_deg * (pi / 180)
mov_dots_aligned <- plot_data2 %>%
  mutate(
    x_new = scale_val * (x_px * cos(angle_rad) - y_px * sin(angle_rad)),
    y_new = scale_val * (x_px * sin(angle_rad) + y_px * cos(angle_rad)),
    x_aligned = x_new + shift_x,
    y_aligned = y_new + shift_y
  )

plot_data <- rbind(
  plot_data1 %>% select(x=x_px, y=y_px, slice),
  mov_dots_aligned %>% select(x = x_aligned, y = y_aligned, slice)
)
ggplot(plot_data, aes(x = x, y = y, color = slice)) +
  geom_point(alpha = 0.5, size = 0.25) +
  scale_color_manual(values = c("S3000" = "black", "S3000_50cs3" = "red")) +
  coord_fixed() + 
  theme_minimal() +
  labs(title = paste0("Alignment Check (Angle: ", angle_deg, "°, Offset: ", shift_x, ", ", shift_y, ")"),
       subtitle = "Black: Reference, Red: Aligned Moving")

pos_ad<-rbind(da.combined@meta.data %>% filter(orig.ident=='B33_s1') %>% 
              select(x_ad=x_px,y_ad=y_px,slice=orig.ident) %>% mutate(x_ad=x_ad+1200),
             plot_data %>% select(x_ad=x,y_ad=y,slice))
head(pos_ad)
ggplot(pos_ad, aes(x = x_ad, y = y_ad, color = slice)) +
  geom_point(alpha = 0.5, size = 0.1) +
  scale_color_manual(values = c("S3000" = "black", "S3000_50cs3" = "red",'B33_s1'='green')) +
  coord_fixed() + 
  theme_minimal() 

write.csv(pos_ad,'A13_aligned.B33_s1.adjusted.clonecell.position.csv')
da.combined@meta.data[rownames(pos_ad),c('x_ad','y_ad')]<-pos_ad[,c('x_ad','y_ad')]
saveRDS(da.combined,"test7.adjusted.Integrate_ct1_tmp.rds")

da.combined$lineage_tmp<-'RGC'
da.combined$lineage_tmp[da.combined$seurat_clusters%in%c(12,21) ]<-'Other'
da.combined$lineage_tmp[da.combined$seurat_clusters%in%c(19) ]<-'Myeloid'

###### Refine RGC ralated cell types
ref_s2<-ref[,ref$lineage_tmp=='RGC' & ref$final_ct3_tmp!='Remove']
DefaultAssay(ref_s2)<-"RNA"
ref_s2<- FindVariableFeatures(ref_s2,nfeatures = 3000) %>% 
        ScaleData() %>% 
        RunPCA(npcs = 20, verbose = FALSE) %>% 
        RunUMAP(reduction = "pca", dims = 1:20) %>%
        FindNeighbors(reduction = "pca", dims = 1:20) %>%
        FindClusters(resolution = 0.6)
ref_s2<-FindClusters(ref_s2,resolution = 0.8)
ct3_map<-c('Astro_1','OPC','Unknown_1','CPN_L23_1','SPN_1','CThPN','PN_L4','Olig_1','Astro_2','GLUT_TH_1',
           'CPN_L5','Unknown_2','Olig_2','CPN_L6','DG','IN_Lhx6','IN_Six3','IN_Lhx6_Sst','IN_Cck','Olig_3',
          'SPN_2','SPN_3','PN_CA3','PN_L6_Nr4a2','CPN_L5_2','PN_CA1','Olig_4','GLUT_TH_2','SCPN_L5','CPN_L23_2',
          'IN_Lhx6_PAL')
names(ct3_map)<-as.character(c(0:30))
ref_s2$final_ct3_tmp2<-as.character(ct3_map[match(ref_s2$seurat_clusters,names(ct3_map))])
unique(ref_s2$final_ct3_tmp2)
ct32_map<-c('DG','Olig','Unknown','Astro','GLUT_TH','Olig','GLUT_TH','SPN','OPC','Astro','IN_Six3','CPN',
            'Unknown','PN_L4','Olig','SPN','PN_L6b','Olig','CPN','SPN','CPN','IN_Lhx6','IN_Lhx6',
            'IN_Cck','CPN','CPN','SCPN','CThPN','CA1','CA3','IN_Lhx6')
names(ct32_map)<-unique(ref_s2$final_ct3_tmp2)
ref_s2$final_ct2_tmp2<-as.character(ct32_map[match(ref_s2$final_ct3_tmp2,names(ct32_map))])
unique(ref_s2$final_ct2_tmp2)
saveRDS(ref_s2,"test7a.snST.RGC.final_ct2.rds")

da2<-da.combined[,da.combined$lineage_tmp=='RGC']#da.combined[,da.combined$lg_type_by_score=='RGC1']
cells<-intersect(colnames(ref_s2),colnames(da2))
da2$final_ct3_tmp2<-NA
da2$final_ct2_tmp2<-NA
da2@meta.data[cells,c('final_ct3_tmp2','final_ct2_tmp2')]<-ref_s2@meta.data[cells,c('final_ct3_tmp2','final_ct2_tmp2')]
da2<-FindVariableFeatures(da2, nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.6) %>%
    RunUMAP(dims = 1:20)
da2<-FindClusters(da2,resolution = 0.65)
ct_map<-c('Astro','Unknown','SPN','CPN_L23','OPC',
           'CThPN','PN_L4','GLUT','CPN_L5','Olig',
          'Olig','IN_Lhx6','DG','CPN_L6','Mixed',
          'SPN','CPN_L23?','IN_Six3','IN_Lhx6','PN_CA3',
          'IN_Vip','PN_CA1','SCPN','PN_L6b','CPN_L5',
          'Olig','GLUT')
names(ct_map)<-as.character(c(0:26))
ct_map
ct_map2<-c('Astro','Unknown','SPN','CPN','OPC',
           'CThPN','PN_L4','GLUT','CPN','Olig',
          'Olig','IN','DG','CPN','Mixed',
          'SPN','CPN?','IN','IN','CA',
          'IN','CA','SCPN','PN_L6b','CPN',
          'Olig','GLUT')
names(ct_map2)<-as.character(c(0:26))
ct_map2
da2$final_ct3<-as.character(ct_map[match(as.character(da2$seurat_clusters),names(ct_map))])
unique(da2$final_ct3)
da2$final_ct2<-as.character(ct_map2[match(as.character(da2$seurat_clusters),names(ct_map2))])
unique(da2$final_ct2)
saveRDS(da2,'test7a.RGC.snST.scST_clonecells.final_ct23.rds')

###########
da.combined<-readRDS("test7.Integrate_ct1_tmp.rds")
da<-readRDS('test7a.RGC.snST.scST_clonecells.final_ct23.rds')
da.combined$final_ct3<-da.combined$lineage_tmp
da.combined$final_ct2<-da.combined$lineage_tmp
da.combined@meta.data[colnames(da),c('final_ct3','final_ct2')]<-da@meta.data[,c('final_ct3','final_ct2')]
unique(da.combined$final_ct3)
unique(da.combined$final_ct2)

da2<-da.combined[,!is.na(da.combined$final_ct3)]
da2$final_ct3<-ifelse(da2$final_ct3=='CPN_L23?','CPN_L23',da2$final_ct3)
da2$final_ct2<-ifelse(da2$final_ct2=='CPN?','CPN',da2$final_ct2)

#pdf("test7a.snST.scST.integrated.final_ct3.UMAP.pdf",width = 11,height = 9)
options(repr.plot.width=11,repr.plot.height=9)
DimPlot(da2,group.by = 'final_ct2',cols = col_ct2,raster = TRUE,
        pt.size = 0.75,label = TRUE,repel = TRUE,
       )+coord_fixed()
#dev.off()

### ct3 umap splited by ident
da2$final_ct3_hl<-ifelse(is.na(da2$nCount_CloneBC),'background',da2$final_ct3)
da2$final_ct2_hl<-ifelse(is.na(da2$nCount_CloneBC),'background',da2$final_ct2)
col_ct2_hl<-get_specific_group_color('ct2')
col_ct2_hl[18]<-'#D3D3D366'
names(col_ct2_hl)[18]<-'background'
col_ct2_hl
da2$final_ct2_hl<-factor(da2$final_ct2_hl,levels=c('background',unique(da2$final_ct2)))
unique(da2$final_ct2_hl)
#pdf("test7a.snST.scST.integrated.final_ct3.only_S3000.highlighted_UMAP.pdf",width = 11,height = 9)
options(repr.plot.width=10,repr.plot.height=9)
DimPlot(da2[,!da2$orig.ident %in% c('S3000_50cs3','B33_s1')],group.by = 'final_ct2_hl',cols = col_ct2_hl,
        raster = FALSE,order = TRUE,
        pt.size = 0.25#,label = TRUE,repel = TRUE,
       )+coord_fixed()
ggsave("test7a.snST.scST.integrated.final_ct3.only_S3000.highlighted_UMAP.png",width = 11,height = 9,dpi = 400)
#dev.off()
#pdf("test7a.snST.scST.integrated.final_ct3.only_S3000_50cs3.highlighted_UMAP.pdf",width = 11,height = 9)
options(repr.plot.width=10,repr.plot.height=9)
DimPlot(da2[,!da2$orig.ident %in% c('S3000','B33_s1')],group.by = 'final_ct2_hl',cols = col_ct2_hl,raster = FALSE,order = TRUE,
        pt.size = 0.25#,label = TRUE,repel = TRUE,
       )+coord_fixed()
ggsave("test7a.snST.scST.integrated.final_ct3.only_S3000_50cs3.highlighted_UMAP.png",width = 11,height = 9,dpi = 400)
#dev.off()


pdf("test7a.snST.scST.integrated.final_ct3.UMAP.pdf",width = 11,height = 9)
options(repr.plot.width=11,repr.plot.height=9)
DimPlot(da2,group.by = 'final_ct3',cols = col_ct3,raster = TRUE,
        pt.size = 0.75,label = TRUE,repel = TRUE,
       )+coord_fixed()
dev.off()
pdf("test7a.snST.scST.integrated.final_ct2.UMAP.pdf",width = 11,height = 9)
options(repr.plot.width=11,repr.plot.height=9)
DimPlot(da2,group.by = 'final_ct2',cols = col_ct2,raster = TRUE,
        pt.size = 0.75,label = TRUE,repel = TRUE
       )+coord_fixed()
dev.off()

da2$sample<-'snST'
#da2$sample[da2$orig.ident=='B33_s1']<-'B33_fixed_cs3'
da2$sample[da2$orig.ident=='S3000']<-'A13'
da2$sample[da2$orig.ident=='S3000_50cs3']<-'A13_cs3'
da2$mouse<-'B29'
#da2$mouse[da2$orig.ident=='B33_s1']<-'B33'
da2$mouse[da2$orig.ident %in% c('S3000','S3000_50cs3')]<-'A13'

pdf('test7a.snST.final_ct3.spatial.pdf',width = 30,height=20)
show_cluster_in_b29_spatial(da2,ref1,ref2,group = 'final_ct3',all = TRUE,col = get_specific_group_color('ct3'),
                            saveplot = TRUE,filename = "test7a.snST.final_ct3"
                           )
#ggsave('test7a.snST.final_ct3.spatial.png',width = 30,height=20,dpi = 300)
dev.off()
pdf('test7a.snST.final_ct2.spatial.pdf',width = 30,height=20)
show_cluster_in_b29_spatial(da2,ref1,ref2,group = 'final_ct2',all = TRUE,col = col_ct2,saveplot = TRUE,filename = "test7a.snST.final_ct2")
dev.off()


mk<-c('Slc17a7',#neuron and excitotary neuron pan-markers
      'Tle4','Pcp4',#CThPN, SCPN
      'Nr4a2',#PN_L6b
      'Satb2',#CPN pan-marker
      'Plxnd1','Cux1',#CPN subtype(CPN_L6 marker?)
      'Rorb','Wfs1','Cpne4','Prox1','Slc17a6',#EN subtypes
      'Slc32a1','Gad1',#IN pan-marker
      'Lhx6','Vip','Six3','Foxp1',#IN sutypes
      'Slc1a3','Pdgfra','Plp1',
      'Ttr','Vtn','Cldn5',#Other
      'Ly86'#Myeloid
     )
mk_ct<-c('CThPN','SCPN','PN_L6b','CPN_L6','CPN_L5','CPN_L23',#'CPN_L23?',
         'PN_L4','PN_CA1','PN_CA3','DG','GLUT','IN_Lhx6','IN_Vip','IN_Six3','SPN',
        'Astro','OPC','Olig',
         'Other','Myeloid',
         'Mixed','Unknown')

### mk_col<-ct_col[mk_ct]
#names(mk_col)<-mk
#mk_col
da2$final_ct3<-factor(da2$final_ct3,levels = mk_ct)

pdf("test7a_final_ct2.snST_scST_A13.canonical_marker_stackViolinPlot.pdf",width = 7,height = 8)
options(repr.plot.width=8,repr.plot.height=8)
VlnPlot(da2[,da2$mouse!='B33'],
        flip = TRUE,assay = 'RNA',
        cols = colorRampPalette(c("#A8C17B", "#B5EAE0", "#F8F9FA", "#B89A9A", "#EBB9B9"),space = "Lab")(27),
        features = paste0('rna_',mk),#combine = TRUE,
        pt.size = 0,group.by = 'final_ct3',stack=TRUE,#sort=TRUE
       )+
    theme(legend.position = 'none')
#ggsave("test7a_final_ct3.canonical_marker_stackViolinPlot.png",width = 8,height = 8,dpi = 400)
dev.off()

write.csv(da2@meta.data,'test7a.snST.scST.integrated.final_ct23.csv')
saveRDS(da2,"test7a.snST.scST.integrated.final_ct23.rds")

###Add Spatial Region
da_adj<-readRDS("./test7.adjusted.Integrate_ct1_tmp.rds")
cell_region<-read.table("../2_region/L9/cell_region.2603.txt",sep = '\t')
colnames(cell_region)<-c("barcode","L1_region","L2_region")
head(cell_region)
cells<-colnames(da_adj)[da_adj$orig.ident %in% c('S3000','S3000_50cs3')]
cell_meta<-da_adj@meta.data[cells,]
cell_meta$barcode<-rownames(cell_meta)
cell_region<-cell_region[cell_region$barcode %in% cells,]
cell_meta<-merge(cell_meta,cell_region,by = 'barcode',all = FALSE)
dim(cell_meta)

df<-cell_meta
for (i in 1:nrow(df)) {
    if (df$L1_region[i]=="Undefined") {
        df$L1_region[i] <- get_most_frequent_value(df, i,"L1_region")
        df$L2_region[i] <- get_most_frequent_value(df, i,"L2_region")  
        #df$L3_region[i] <- get_most_frequent_value(df, i,"L3_region")  
  }    
} 
rownames(df)<-df$barcode
table(df$L1_region)
table(df$L2_region)
da_adj@meta.data[,c("L1_region","L2_region")]<-df[colnames(da_adj),c("L1_region","L2_region")]

plot_SpatialDim_across_sections(da_adj,group = 'final_ct3',col = get_specific_group_color('ct3'),pt_size = 1.2,
                                saveplot = TRUE,filepath = getwd()
                               )
plot_SpatialDim_across_sections(da_adj,group = 'final_ct2',col = get_specific_group_color('ct2'),pt_size = 1,
                                saveplot = TRUE,filepath = getwd()
                               )

plot_SpatialDim_across_sections(da_adj,group = 'L1_region',col = get_specific_group_color('L1_region'),pt_size = 1.2,
                                saveplot = TRUE,filepath = getwd()
                               )
plot_SpatialDim_across_sections(da_adj,group = 'L2_region',col = get_specific_group_color('L2_region'),pt_size = 1.2,
                                saveplot = TRUE,filepath = getwd()
                               )

cells<-colnames(da_adj)[da_adj$orig.ident %in% c(#'B33_s1',
                                                 'S3000','S3000_50cs3')]

saveRDS(da_adj,'test7a.adjusted.final_ct23.region_L12.2603.rds')
write.csv(da_adj@meta.data,'test7a.adjusted.final_ct23.region_L12.2603.csv')

#### relative abundance dotplot and barplot
da_adj<-readRDS('test7a.adjusted.final_ct23.region_L12.2603.rds')
table(colnames(da2) %in% colnames(da_adj))
da2@meta.data[colnames(da_adj),c("L1_region","L2_region",#"L3_region",
                                 'x_ad','y_ad')]<-da_adj@meta.data[,c("L1_region","L2_region",#"L3_region",
                                                                      'x_ad','y_ad')]

tmp_df<-da2@meta.data %>% filter(orig.ident %in% c(#'B33_s1',
                                                   'S3000','S3000_50cs3'))
tmp_df<-tmp_df[!is.na(tmp_df$L1_region),]
tmp_df$final_ct3<-factor(tmp_df$final_ct3,
                         levels = c('CPN_L23',#'CPN_L23?',
                                    'PN_L4','CPN_L5','CPN_L6','PN_L6b','SCPN','CThPN','PN_CA1',
                                    'PN_CA3','DG','IN_Lhx6','IN_Vip','SPN','GLUT','IN_Six3',
                                    'Astro','OPC','Olig','Myeloid',
                                    'Other','Mixed','Unknown'))
tmp_df$L2_region<-factor(tmp_df$L2_region,
                        levels=c('Ctx_L1','Ctx_L23','Ctx_L4','Ctx_L5','Ctx_L6','OLF',
                                 'CA1','CA3','DG','WM?',
                                 'STR','CTXsp','PAL',
                                 'HY','Hb','LD','MD','AV','AM','VA/VL/VM','RE','ZI',
                                 'RT','VZ',
                                 'Other'
                                )
                        )
get_domain_celltype_bubble_bar_plot_v2(tmp_df,group1 = 'final_ct3',group2 = 'L2_region',w=9,h=20,mulp = 4,
                                       col1 = col_ct3,col2 = col_r2,saveplot = TRUE,outdir=getwd()
                                      )

saveRDS(da2,'test7a.snST.scST.integrated.final_ct23.L12_region.2603.rds')

### Fate
da2$fate<-"Undefined"
da2$fate[da2$final_ct2 %in% c('CPN','SCPN','PN_L4','PN_L6b','CThPN')]<-'PN'
da2$fate[da2$final_ct2 %in% c('CA')]<-'CA'
da2$fate[da2$final_ct2 %in% c('DG')]<-'DG'
da2$fate[da2$final_ct3 %in% c('IN_Lhx6','IN_Vip')]<-'IN'
da2$fate[da2$final_ct3 %in% c('IN_Six3')]<-'GABA'
da2$fate[da2$final_ct2 %in% c('SPN')]<-'SPN'
da2$fate[da2$final_ct2 %in% c('Olig','Astro','OPC')]<-'Glia'
da2$fate[da2$final_ct2 %in% c('GLUT')]<-'GLUT'
da2@meta.data[rownames(da2@reductions$umap@cell.embeddings),c('UMAP_1','UMAP_2')]<-da2@reductions$umap@cell.embeddings[,c('UMAP_1','UMAP_2')]
write.csv(da2@meta.data,'test7a.snST.scST.integrated.final_ct23.L123_region.umap.csv')


















