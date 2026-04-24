### based on L9 expression matrix
library(purrr)
library(cowplot)
library(png)
library(grid)
#library(pdftools)
library(cluster)
library(ggpubr)
library(stringr)
library(plotly)
library(patchwork)
library(RANN)
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")
#source("~/NPC_project/BMK/pipeline/spatial_plot.r")

wd<-"~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/Merged_2512/2_region/L9"
setwd(wd)
getwd()
if(!dir.exists(wd)){dir.create(wd,recursive=TRUE)}
script.dir <- "~/NPC_project/BMK/pipeline/RNA" #dirname(script.name) 

obj1<-readRDS("object.spatial.clustered.rds")
obj2<-readRDS("~/NPC_project/BMK/tissue/brain/P5/Merged_2512/2_region/L9/P5_S3000_object.spatial.clustered.rds")
obj3<-readRDS("~/NPC_project/BMK/tissue/brain/P5/Merged_2512/2_region/L9/P5_S3000_50cs3_object.spatial.clustered.rds")

ifnb.list <- list(obj1,obj2,obj3)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    DefaultAssay(x)<-'Spatial'
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
saveRDS(da.combined,"test3.integration.spatial.clustered.region.rds")

da.combined <- FindClusters(da.combined, resolution = 0.25)

#res=0.15
options(repr.plot.height=10,repr.plot.width=30)
tmp_col<-get_specific_group_color('div_col')[1:n_distinct(da.combined$seurat_clusters)]
names(tmp_col)<-unique(da.combined$seurat_clusters)
SpatialDimPlot(da.combined,cols = tmp_col,pt.size.factor = 10,crop = FALSE,label = TRUE)

rg1_map<-c('HY','Ctx','subCtx','TH','HIP','Ctx','Other','subCtx','TH','Tract','TH','Other')
names(rg1_map)<-c(0:11)
rg1_map
da.combined$L1_region_tmp<-as.character(rg1_map[match(da.combined$integrated_snn_res.0.2,names(rg1_map))])
unique(da.combined$L1_region_tmp)

options(repr.plot.height=10,repr.plot.width=30)
tmp_col<-get_specific_group_color('div_col')[1:n_distinct(da.combined$L1_region_tmp)]
names(tmp_col)<-unique(da.combined$L1_region_tmp)
SpatialDimPlot(da.combined,group.by = 'L1_region_tmp',cols = tmp_col,pt.size.factor = 10,crop = FALSE)

saveRDS(da.combined[,da.combined$L1_region_tmp=='Ctx'],"test3.Integrate_Ctx_tmp.rds")
saveRDS(da.combined[,da.combined$L1_region_tmp=='subCtx'],"test3.Integrate_subCtx_tmp.rds")
saveRDS(da.combined[,da.combined$L1_region_tmp=='HIP'],"test3.Integrate_HIP_tmp.rds")
saveRDS(da.combined[,da.combined$L1_region_tmp=='TH'],"test3.Integrate_TH_tmp.rds")
saveRDS(da.combined[,da.combined$L1_region_tmp=='HY'],"test3.Integrate_HY_tmp.rds")
saveRDS(da.combined[,da.combined$L1_region_tmp%in%c('Other')],"test3.Integrate_Other_tmp.rds")

###After brain region annotation of each structure
da_ctx<-readRDS('test3.Integrate_Ctx_L3_region_tmp.rds')
da_subctx<-readRDS('test3.Integrate_subCtx_L3_region_tmp.rds')
da_hy<-readRDS('test3.Integrate_HY_L3_region_tmp.rds')
da_th<-readRDS('test3.Integrate_TH_L3_region_tmp.rds')#da_hy<-readRDS('test3.Integrate_HY_L3_region_tmp.rds')
da_hip<-readRDS('test3.Integrate_HIP_L3_region_tmp.rds')
da_other<-readRDS("test3.Integrate_Other_L3_region_tmp.rds")

cells<-c(colnames(da_ctx),colnames(da_subctx),colnames(da_hip),colnames(da_th),colnames(da_hy),colnames(da_other))
ct3_tmp<-as.character(c(da_ctx$L3_region_tmp,da_subctx$L3_region_tmp,da_hip$L3_region_tmp,
                        da_th$L3_region_tmp,da_hy$L3_region_tmp,da_other$L3_region_tmp))
head(cells)
unique(ct3_tmp)

da.combined$L3_region_tmp<-'Other'
da.combined@meta.data[cells,'L3_region_tmp']<-ct3_tmp
unique(da.combined$L3_region_tmp)

##Cluster spatial smoothing
# apply smooth
test_cl<-c()
for(i in unique(da.combined$orig.ident)){
    tmp_da<-da.combined@meta.data[da.combined$orig.ident==i,]
    tmp_me<-smooth_spatial_clusters(tmp_da, "L3_region_tmp", k = 8)
    names(tmp_me)<-rownames(tmp_da)
    test_cl<-c(test_cl,tmp_me)
}

da.combined@meta.data[names(test_cl),'L3_region_sm']<-as.character(test_cl)

#find nearest neighbor's region to undefined cells
df<-da.combined@meta.data
for(i in unique(df$orig.ident)){
    tmp_df<-df %>% filter(orig.ident==i)
    #print(dim(tmp_df))
    for (i in 1:nrow(tmp_df)) {
        if (tmp_df$L3_region_sm[i]=="Undefined") {
            #df$L1_region[i] <- get_most_frequent_value(df, i,"L1_region")
            #df$L2_region[i] <- get_most_frequent_value(df, i,"L2_region")  
            tmp_df$L3_region_sm[i] <- get_most_frequent_value(tmp_df, i,"L3_region_sm")  
        }    
    }
    df[rownames(tmp_df),'L3_region_sm']<-tmp_df[,'L3_region_sm']  
} 

da.combined@meta.data[rownames(df),'L3_region_sm']<-df$L3_region_sm

da.combined$L3_region_sm<-ifelse(da.combined$L3_region_sm=='LVZ','Other',da.combined$L3_region_sm)
table(da.combined$L3_region_sm)

tmp_col<-get_specific_group_color('div_col')[1:length(unique(da.combined$L3_region_sm))]
names(tmp_col)<-unique(da.combined$L3_region_sm)
tmp_col

p1<-plot_spatial_group(da.combined,'S3000','L3_region_sm',tmp_col,pt.size = 0.8)
p2<-plot_spatial_group(da.combined,'S3000_50cs3','L3_region_sm',tmp_col,pt.size = 0.8)
p3<-plot_spatial_group(da.combined,'B33_s1','L3_region_sm',tmp_col,pt.size = 0.8)
wrap_plots(p1,p2,p3)+
    plot_layout(nrow = 1, guides = "collect"
               ) &
    theme(
      legend.position = "right",
      legend.key.height = unit(0.001, "cm"), # 极短长条
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8)
    )

rg3_rg2_map<-c('Other','Ctx_L23','Ctx_L1','Ctx_L4','Ctx_L5','Ctx_L6','Ctx_L6','CA1','Ctx_L23','Ctx_L23','Ctx_L6','Ctx_L4',
          'Ctx_L5','CA3','CA1','CA1','Other','TH_5','CA1','DG','DG','TH_3','CA3','TH_2','TH_6','TH_8','TH_7','STR','TH_1',
              'TH_9','HY_1','HY_6','TH_4','HY_4','PAL_1','HY_2','Ctx_L6','PAL','HY_3','HY_7','PAL','HY_5','Ctx_L6','STR','STR',
              'Ctx_L5','STR','OLF','STR','STR','Ctx_L23','Ctx_L23','VZ3')
names(rg3_rg2_map)<-unique(da.combined$L3_region_sm)
rg3_rg2_map
da.combined$L2_region_tmp<-as.character(rg3_rg2_map[match(da.combined$L3_region_sm,names(rg3_rg2_map))])
unique(da.combined$L2_region_tmp)

rg2_rg1_map<-c('Other','Ctx','Ctx','Ctx','Ctx','Ctx','Hip','Hip','TH','Hip','TH','TH','TH','TH','TH','subCtx','TH',
          'TH','HY','HY','TH','HY','subCtx','HY','subCtx','HY','HY','HY','Ctx','VZ')
names(rg2_rg1_map)<-unique(da.combined$L2_region_tmp)
rg2_rg1_map
da.combined$L1_region<-as.character(rg2_rg1_map[match(da.combined$L2_region_tmp,names(rg2_rg1_map))])
unique(da.combined$L1_region)

saveRDS(da.combined,"test3.Integrate_L123_region.2601.rds")

write.csv(da.combined@meta.data,'test3.Integrate_L123_region.meta.2601.csv')

### adjust image and delineate Left-right
meta_tmp<-da.combined@meta.data[da.combined$orig.ident=='B33_s1',]
dim(meta_tmp)
#img_xmin <- min(meta_tmp$x.axis); img_xmax <- max(meta_tmp$x.axis)
#img_ymin <- min(meta_tmp$y.axis);  img_ymax <- max(meta_tmp$y.axis)
plot_data <- align_spots(meta_tmp, img_h = img_h,img_w = img_w,
                         scale_x = 0.79, scale_y=0.93,
                         off_x = 168, off_y = 5, rotate = 0,flip_y = TRUE)
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
img<-da_all@images$S3000@image
img_rot <- aperm(img, c(2, 1, 3))[,rev(1:nrow(img)), ]
img_w <- ncol(img_rot)
img_h <- nrow(img_rot)
img_grob <- rasterGrob(img_rot, interpolate = TRUE)

meta_tmp<-da_all@meta.data[da_all$orig.ident=='S3000',]
dim(meta_tmp)
plot_data <- align_spots(meta_tmp, img_h = img_h,img_w = img_w,
                         scale_x = 0.935, scale_y=0.85,
                         off_x = 65, off_y = 80, rotate = 270,flip_y = TRUE)
head(plot_data)
ggplot(plot_data, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void() +
  labs(subtitle = "")

da.combined@meta.data[rownames(plot_data),c('x_px','y_px')]<-plot_data[,c('x_px','y_px')]
da.combined@images$S3000@image<-img_rot

#S3000_50cs3
img<-da_all@images$S3000_50cs3@image
img_rot <- aperm(img, c(2, 1, 3))[rev(1:ncol(img)),, ]
img_w <- ncol(img_rot)
img_h <- nrow(img_rot)
img_grob <- rasterGrob(img_rot, interpolate = TRUE)
meta_tmp<-da_all@meta.data[da_all$orig.ident=='S3000_50cs3',]
dim(meta_tmp)
plot_data <- align_spots(meta_tmp, img_h = img_h,img_w = img_w,
                         scale_x = 1, scale_y=0.87,
                         off_x = 0, off_y = 6, rotate = 90,flip_y = TRUE)
head(plot_data)
ggplot(plot_data, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void() +
  labs(subtitle = "")

da.combined@meta.data[rownames(plot_data),c('x_px','y_px')]<-plot_data[,c('x_px','y_px')]
da.combined@images$S3000_50cs3@image<-img_rot

da_all<-readRDS(file.path(wkdir,"test3.Integrate_L123_region.rds"))
da_all
#S3000
img<-da_all@images[['S3000']]@image
img_w <- ncol(img)
img_h <- nrow(img)
img_grob <- rasterGrob(img, interpolate = TRUE)
df1<-da.combined@meta.data[da.combined$orig.ident=='S3000',]
x_anchor <- 630  
y_anchor <- 800  
manual_angle <- 0.75
manual_res <- correct_brain_side_manual(df1, x_anchor, y_anchor, manual_angle)
ggplot(manual_res$data, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  geom_point(aes(color = side), size = 1, alpha = 0.4) +
  annotate("point", x = x_anchor, y = y_anchor, color = "yellow", size = 4, shape = 4) +
  
  scale_color_manual(values = c("L" = "#0072B2", "R" = "#D55E00")) +
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void()

da.combined$LR<-NA
da.combined@meta.data[rownames(manual_res$data),'LR']<-manual_res$data$side
# --------------------------------
#S3000_50cs3
img<-da_all@images[['S3000_50cs3']]@image
img_w <- ncol(img)
img_h <- nrow(img)
img_grob <- rasterGrob(img, interpolate = TRUE)
df1<-da.combined@meta.data[da.combined$orig.ident=='S3000_50cs3',]
x_anchor <- 500  
y_anchor <- 800  
manual_angle <- 2 
manual_res <- correct_brain_side_manual(df1, x_anchor, y_anchor, manual_angle)
ggplot(manual_res$data, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  geom_point(aes(color = side), size = 1, alpha = 0.4) +
  annotate("point", x = x_anchor, y = y_anchor, color = "yellow", size = 4, shape = 4) +
  
  scale_color_manual(values = c("L" = "#0072B2", "R" = "#D55E00")) +
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void()
da.combined@meta.data[rownames(manual_res$data),'LR']<-manual_res$data$side
da.combined@meta.data[da.combined$orig.ident=='B33_s1','LR']<-'L'
saveRDS(da.combined,file.path(wd,'test3.integration.spatial.region.adj_img_pos.LR.rds'))
write.csv(da.combined@meta.data,file.path(wd,"test3.integration.spatial.region.adj_img_pos.LR.meta.csv"))


meta<-read.csv("test3.Integrate_L123_region.meta.2601.csv",row.names = 1,header = TRUE)
head(meta)

df1<-read.csv("S3000_L9_map2_L1.txt",sep = '\t',header = FALSE)
colnames(df1)<-c("L1","L9")
df1<-as.data.frame(apply(df1,2,FUN = function(x){paste0("S3000_",x)}))
head(df1)

df2<-read.csv("S3000_50cs3_L9_map2_L1.txt",sep = '\t',header = FALSE)
colnames(df2)<-c("L1","L9")
df2<-as.data.frame(apply(df2,2,FUN = function(x){paste0("S3000_50cs3_",x)}))
head(df2)

df3<-read.csv("S3000_50cs3_L9_map2_L1.txt",sep = '\t',header = FALSE)
colnames(df3)<-c("L1","L9")
df3<-as.data.frame(apply(df3,2,FUN = function(x){paste0("B33_s1_",x)}))
head(df3)

### 20260311: adjust the L1/L2 brain region annotation
tmp_meta<-da.combined@meta.data
tmp_meta$L2_region_sm2<-tmp_meta$L3_region_sm
tmp_meta$L2_region_sm2[tmp_meta$L2_region_sm2 %in% c('Ctx_L6_2','STR_6','STR_3')]<-'CTXsp'
tmp_meta$L2_region_sm2[tmp_meta$L2_region_sm2 %in% c('OLF','STR_5')]<-'OLF'
tmp_meta$L2_region_sm2[grepl('HY_',tmp_meta$L2_region_sm2)]<-'HY'
tmp_meta$L2_region_sm2[grepl('PAL_',tmp_meta$L2_region_sm2)]<-'PAL'
tmp_meta$L2_region_sm2[grepl('Ctx_L23_',tmp_meta$L2_region_sm2)]<-'Ctx_L23'
tmp_meta$L2_region_sm2[grepl('Ctx_L4_',tmp_meta$L2_region_sm2)]<-'Ctx_L4'
tmp_meta$L2_region_sm2[grepl('Ctx_L5_',tmp_meta$L2_region_sm2)]<-'Ctx_L5'
tmp_meta$L2_region_sm2[grepl('Ctx_L6_',tmp_meta$L2_region_sm2)]<-'Ctx_L6'
tmp_meta$L2_region_sm2[grepl('STR_',tmp_meta$L2_region_sm2)]<-'STR'
tmp_meta$L2_region_sm2[grepl('CA1_',tmp_meta$L2_region_sm2)]<-'CA1'
tmp_meta$L2_region_sm2[grepl('CA3_',tmp_meta$L2_region_sm2)]<-'CA3'
tmp_meta$L2_region_sm2[grepl('DG_',tmp_meta$L2_region_sm2)]<-'DG'
tmp_meta$L2_region_sm2[grepl('VZ*',tmp_meta$L2_region_sm2)]<-'VZ'

map_ad<-paste0('TH_',c(1:9))
names(map_ad)<-c('RT','LD','MD','ZI','Hb','AM','RE','AV','VA/VL/VM')
#map_ad
tmp_ad<-tmp_meta %>% filter(L2_region_sm2 %in% map_ad)

tmp_ad$L2_region_sm2<-as.character(names(map_ad)[match(tmp_ad$L2_region_sm2,map_ad)])

tmp_meta[rownames(tmp_ad),'L2_region_sm2']<-tmp_ad$L2_region_sm2

da.combined$L2_region_sm2<-tmp_meta[colnames(da.combined),'L2_region_sm2']
col_r2<-get_specific_group_color('L2_region')

#png("L3_region_sm.spatial.png",height = 1000,width = 3000,res = 300)
options(repr.plot.height=10,repr.plot.width=20)
p1<-plot_spatial_group(da.combined,'S3000','L2_region_sm2',col_r2,pt.size = 2)
p2<-plot_spatial_group(da.combined,'S3000_50cs3','L2_region_sm2',col_r2,pt.size = 2)
#p3<-plot_spatial_group(da.combined,'B33_s1','L3_region_sm',col_r3,pt.size = 2)
wrap_plots(p1,p2
          )+
    plot_layout(nrow = 1, guides = "collect"
               ) &
    theme(
      legend.position = "none"#"right",
      #legend.key.height = unit(0.001, "cm"),
      #legend.spacing.y = unit(0, "cm"),
      #legend.text = element_text(size = 8)
    )
ggsave("2603.L2_region_sm2.spatial.png",height = 10,width = 20,dpi = 300)

rg2_rg1_map<-c('Other',rep('Ctx',5),rep('Hip',2),'Other','EPI','Hip',rep('TH',5),'CNU',rep('TH',2),'HY','TH','CNU','Ctx','Ctx','VZ')
names(rg2_rg1_map)<-unique(da.combined$L2_region_sm2)
rg2_rg1_map
da.combined$L1_region_sm2<-as.character(rg2_rg1_map[match(da.combined$L2_region_sm2,names(rg2_rg1_map))])
unique(da.combined$L1_region_sm2)

col_r1<-get_specific_group_color('L1_region')
#png("L3_region_sm.spatial.png",height = 1000,width = 3000,res = 300)
options(repr.plot.height=10,repr.plot.width=20)
p1<-plot_spatial_group(da.combined,'S3000','L1_region_sm2',col_r1,pt.size = 2)
p2<-plot_spatial_group(da.combined,'S3000_50cs3','L1_region_sm2',col_r1,pt.size = 2)
#p3<-plot_spatial_group(da.combined,'B33_s1','L3_region_sm',col_r3,pt.size = 2)
wrap_plots(p1,p2
          )+
    plot_layout(nrow = 1, guides = "collect"
               ) &
    theme(
      legend.position = "none"#"right",
      #legend.key.height = unit(0.001, "cm"), 
      #legend.spacing.y = unit(0, "cm"),
      #legend.text = element_text(size = 8)
    )
ggsave("2603.L1_region_sm2.spatial.png",height = 10,width = 20,dpi = 300)

saveRDS(da.combined,"test3.Integrate_L123_region.2603.rds")
write.csv(da.combined@meta.data,'test3.Integrate_L123_region.meta.2603.csv')

meta<-read.csv("test3.Integrate_L123_region.meta.2603.csv",row.names = 1,header = TRUE)
head(meta)

df1<-read.csv("S3000_L9_map2_L1.txt",sep = '\t',header = FALSE)
colnames(df1)<-c("L1","L9")
df1<-as.data.frame(apply(df1,2,FUN = function(x){paste0("S3000_",x)}))
head(df1)

df2<-read.csv("S3000_50cs3_L9_map2_L1.txt",sep = '\t',header = FALSE)
colnames(df2)<-c("L1","L9")
df2<-as.data.frame(apply(df2,2,FUN = function(x){paste0("S3000_50cs3_",x)}))
head(df2)

df<-rbind(df1,df2)
df<-df[df$L9 %in% rownames(meta),]
#df[,c('L1_region','L2_region','L3_region')]<-"Undefined"
df[,c('L1_region','L2_region')]<-meta[df$L9,c('L1_region_sm2','L2_region_sm2')]
head(df)

write.table(df,"L9_L1_region.2603.txt",sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)

cell_df1<-read.csv("~/NPC_project/BMK/tissue/brain/P5/S3000/RNA/BST/07.CellSplit/cell_split_result/all_barcode_num.txt",sep = '\t',header = FALSE)
colnames(cell_df1)<-c("L1","cell")
cell_df1<-as.data.frame(apply(cell_df1,2,FUN = function(x){paste0("S3000_",x)}))
head(cell_df1)

cell_df2<-read.csv("~/NPC_project/BMK/tissue/brain/P5/S3000_50cs3/RNA/BST/07.CellSplit/cell_split_result/all_barcode_num.txt",sep = '\t',header = FALSE)
colnames(cell_df2)<-c("L1","cell")
cell_df2<-as.data.frame(apply(cell_df2,2,FUN = function(x){paste0("S3000_50cs3_",x)}))
head(cell_df2)

cell_df<-rbind(cell_df1,cell_df2)
#cell_df<-cell_df#[!duplicated(cell_df$cell),]
#rownames(cell_df)<-cell_df$L1
dim(cell_df)

write.table(cell_df,file.path(wd,"cell_L1.2603.txt"),sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)
