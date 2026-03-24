library(Matrix)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)
library(purrr)
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")
exp_col<-c('#034f84','#92a8d1','#d6d4e0','#f4a688',"#ED797B","#d64161","#c94c4c")

wd<-"~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/Merged_2512/5_SpatialFate"
if(!dir.exists(wd)){
    dir.create(wd)
}

da<-readRDS("../3_clone/test7a.snST.scST.integrated.final_ct23.L12_region.2603.rds")#test7a.snST.scST.integrated.final_ct23.L123_region.rds")
da_adj<-readRDS("../3_clone/test7a.adjusted.final_ct23.region_L12.2603.rds")#test7a.adjusted.final_ct23.region_L123.rds")

### re-align cells on two sections of A13 by non-rigid registration
#S3000
img<-da_adj@images[['S3000']]@image
img_w <- ncol(img)
img_h <- nrow(img)
img_grob <- rasterGrob(img, interpolate = TRUE)
df1<-da_adj@meta.data[da_adj$orig.ident=='S3000',]
anchor1 <- c(633,765)#c(445,692)#c(435,812)
anchor2 <- c(365,98)#
ggplot(df1, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  annotate("point", x = anchor1[1], y = anchor1[2], color = "yellow", size = 4, shape = 4) +

  annotate("point", x = anchor2[1], y = anchor2[2], color = "red", size = 4, shape = 4) +
  coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
  theme_void() 

#S3000_50cs3
img2<-da_adj@images[['S3000_50cs3']]@image
img_w2 <- ncol(img2)
img_h2 <- nrow(img2)
img_grob2 <- rasterGrob(img2, interpolate = TRUE)
df2<-da_adj@meta.data[da_adj$orig.ident=='S3000_50cs3',]
anchor21 <- c(500,715)#c(292,630)#c(276,742) 
anchor22 <- c(235,38)
ggplot(df2, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob2, xmin = 0, xmax = img_w2, ymin = 0, ymax = img_h2) +
  annotate("point", x = anchor21[1], y = anchor21[2], color = "yellow", size = 4, shape = 4) +

  annotate("point", x = anchor22[1], y = anchor22[2], color = "red", size = 4, shape = 4) +
  coord_fixed(xlim = c(0, img_w2), ylim = c(0, img_h2)) +
  theme_void() 

vec_ref <- c(anchor22[1] - anchor21[1], anchor22[2] - anchor21[2])
vec_mov <- c(anchor2[1] - anchor1[1], anchor2[2] - anchor1[2])

dist_ref <- sqrt(sum(vec_ref^2))
dist_mov <- sqrt(sum(vec_mov^2))
suggest_scale <- dist_ref / dist_mov #1.07#

angle_ref <- atan2(vec_ref[2], vec_ref[1])
angle_mov <- atan2(vec_mov[2], vec_mov[1])
suggest_angle_deg <- (angle_ref - angle_mov) * 180 / pi

cat("suggested scale factor:", suggest_scale, "\n")
cat("suggected rotation angle:", suggest_angle_deg, "\n")
suggest_angle_rad <- angle_ref - angle_mov

m1_x_rotated <- suggest_scale * (anchor1[1] * cos(suggest_angle_rad) - anchor1[2] * sin(suggest_angle_rad))
m1_y_rotated <- suggest_scale * (anchor1[1] * sin(suggest_angle_rad) + anchor1[2] * cos(suggest_angle_rad))
suggest_tx <- anchor21[1] - m1_x_rotated
suggest_ty <- anchor21[2] - m1_y_rotated
cat("suggested x move:", suggest_tx, "\n")
cat("suggested y move:", suggest_ty, "\n")
rigid_angle_rad
rigid_tx
rigid_ty
df1$slice<-'S3000'
df2$slice<-'S3000_50cs3'

shift_x <- -150   
shift_y <- -125  
angle_deg <- 2.8  
scale_val <- 1.07    
angle_rad <- angle_deg * (pi / 180)

mov_dots_aligned <- df1 %>%
  mutate(
    x_new = scale_val * (x_px * cos(angle_rad) - y_px * sin(angle_rad)),
    y_new = scale_val * (x_px * sin(angle_rad) + y_px * cos(angle_rad)),
    x_aligned = x_new + shift_x,
    y_aligned = y_new + shift_y
  )

df_anchor1<-as.data.frame(rbind(anchor1,anchor2))
colnames(df_anchor1)<-c('x_px','y_px')
df_anchor1
df_anchor1_ad<-df_anchor1 %>%
  mutate(
    x_new = scale_val * (x_px * cos(angle_rad) - y_px * sin(angle_rad)),
    y_new = scale_val * (x_px * sin(angle_rad) + y_px * cos(angle_rad)),
    x_aligned = x_new + shift_x,
    y_aligned = y_new + shift_y
  )
df_anchor1_ad

plot_data <- rbind(
  df2 %>% select(x=x_px, y=y_px, slice),
  mov_dots_aligned %>% select(x = x_aligned, y = y_aligned, slice)
)
ggplot(plot_data, aes(x = x, y = y, color = slice)) +
  geom_point(alpha = 0.5, size = 0.25) +
  scale_color_manual(values = c("S3000" = "black", "S3000_50cs3" = "red")) +
  coord_fixed() + 
  theme_minimal() +
  labs(title = paste0("Alignment Check (Angle: ", angle_deg, "°, Offset: ", shift_x, ", ", shift_y, ")"),
       subtitle = "Black: Reference, Red: Aligned Moving")

meta<-da@meta.data %>% filter(mouse %in% c('A13','B33'))
table(is.na(meta$L2_region))
for(i in unique(meta$orig.ident)){
    tmp_df<-meta %>% filter(orig.ident==i)
    #print(dim(tmp_df))
    for (i in 1:nrow(tmp_df)) {
        if (is.na(tmp_df$L2_region[i])) {
            #df$L1_region[i] <- get_most_frequent_value(df, i,"L1_region")
            #df$L2_region[i] <- get_most_frequent_value(df, i,"L2_region") 
            if(!is.null(get_most_frequent_value(tmp_df, i,"L2_region"))){
                tmp_df$L2_region[i] <- get_most_frequent_value(tmp_df, i,"L2_region")  
                tmp_df$L3_region[i] <- get_most_frequent_value(tmp_df, i,"L3_region")  
            }else{
                tmp_df$L2_region[i] <- 'Undefined' 
                tmp_df$L3_region[i] <- 'Undefined'
            }        
        }    
    }
    meta[rownames(tmp_df),'L3_region']<-tmp_df[,'L3_region']  
    meta[rownames(tmp_df),'L2_region']<-tmp_df[,'L2_region'] 
} 

test_df<-meta
test_df[,c('x_ad','y_ad')]<-plot_data[rownames(test_df),c('x','y')]

df_ad<-read.csv('../3_clone/test7a.scST.clone_cell.lineage_final_ct2.error_adjusted.meta.csv',row.names = 1)
rownames(df_ad)<-df_ad$cell
head(df_ad)
test_df<-test_df[rownames(df_ad),] %>% mutate(fate=df_ad$fate,UMAP_1=df_ad$UMAP_1,UMAP_2=df_ad$UMAP_2) %>% filter(mouse=='A13')
write.csv(test_df,'test7a.snST.scST_A13.integrated.final_ct23.L123_region.umap.pos_ad.csv')

### add region background
cell_region<-read.csv("A13.all.cell.region.LR.pos_ad.meta.csv",row.names = 1)
head(cell_region)
meta<-test_df %>% filter(!is.na(clone.id2),clone.id2!='none',mouse=='A13')
#head(meta)
meta$clone<-paste0(meta$mouse,"_",meta$clone.id2)
meta$cell<-rownames(meta)
#meta[,c('x_ad','y_ad')]<-test_df[rownames(meta),c('x_ad','y_ad')]

### Hippocampus neuron spatial orientation
ca_clones<-meta %>% filter(fate%in%c('CA','DG'),mouse=='A13') %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2)
head(ca_clones)
### remove outlier clone pairs by spatial distance
dist_results <- ca_clones %>%
  group_by(clone) %>%
  filter(n() >= 2) %>%
  do({
    coords <- as.matrix(.[, c("x_ad", "y_ad")])
    dists <- as.numeric(dist(coords))
    data.frame(distance = dists)
  }) %>%
  ungroup()
head(dist_results)
options(repr.plot.width=5,repr.plot.height=5)
ggplot(dist_results, aes(x = distance)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#7897AB", alpha = 0.6) +
  geom_density(color = "#D885A3", size = 1) +
  stat_ecdf(geom = "step", color = "#98BA7D", size = 0.8) +
  theme_minimal() +
  labs(title = "Distribution of Clonal Inter-cell Distances",
       x = "Distance (units)", y = "Density / Cumulative Prob")

threshold <- quantile(dist_results$distance, 0.9) 
dist_filtered <- dist_results %>%
  filter(distance <= threshold)
ggplot(dist_filtered, aes(x = distance)) +
  geom_density(fill = "#7897AB", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Filtered Clonal Distances (Excluding Outliers)")
df_edges <- get_clonal_edges(ca_clones)
dist_threshold <- threshold

df_edges_filtered <- df_edges %>%
  filter(distance <= dist_threshold)
cell_sub<-unique(c(df_edges_filtered$cell1,df_edges_filtered$cell2))
length(cell_sub)
clone_df<-ca_clones %>% filter(cell %in% cell_sub) %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2) %>% as.data.frame()
rownames(clone_df)<-clone_df$cell
head(clone_df)
write.csv(clone_df,'CA_DG.clones.outlier_removed.pos_ad.meta.csv')
object<-da_adj
section<-c('S3000_50cs3')
group<-'final_ct3'
pt_size=2
p_list<-list()
df_tmp<-clone_df
for (ident in section){
    img <- object@images[[ident]]@image
    img_w <- ncol(img); img_h <- nrow(img)
    img_grob <- grid::rasterGrob(img, interpolate = TRUE)
    df<-df_tmp #%>% filter(orig.ident==ident)
    p_list[[ident]] <- ggplot(df,aes(x = x_ad,y = y_ad))+ 
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
      geom_path(aes(group = clone), alpha = 0.8, size = 0.3,linewidth = 0.5,color='white') +
      geom_point(size=pt_size,shape=21,stroke=0.1,color='white',aes_string(fill=group))+
      scale_fill_manual(values = col_ct3)+
      xlab(paste0("")) +
      ylab(paste0("")) + 
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
      theme_void() +
      labs(title = ident)+    
      theme(
      legend.key.height = unit(0.001, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8),

    legend.position = "right"
    )
}
options(repr.plot.width=10,repr.plot.height=10)
wrap_plots(p_list)
ggsave('CA_DG_fate.spatial_clone_final_ct3.cs3_plot.png',width = 10,height = 10,dpi = 500)

##### Cortical PN spatial distribution
pn_clones<-meta %>% filter(fate=='PN',mouse=='A13') %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2)
head(pn_clones)
### remove outlier clone pairs by spatial distance
dist_results <- pn_clones %>%
  group_by(clone) %>%
  filter(n() >= 2) %>%
  do({
    coords <- as.matrix(.[, c("x_ad", "y_ad")])
    dists <- as.numeric(dist(coords))
    data.frame(distance = dists)
  }) %>%
  ungroup()
head(dist_results)
options(repr.plot.width=5,repr.plot.height=5)
ggplot(dist_results, aes(x = distance)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#7897AB", alpha = 0.6) +
  geom_density(color = "#D885A3", size = 1) +
  stat_ecdf(geom = "step", color = "#98BA7D", size = 0.8) +
  theme_minimal() +
  labs(title = "Distribution of Clonal Inter-cell Distances",
       x = "Distance (units)", y = "Density / Cumulative Prob")

threshold <- quantile(dist_results$distance, 0.8)  

dist_filtered <- dist_results %>%
  filter(distance <= threshold)

ggplot(dist_filtered, aes(x = distance)) +
  geom_density(fill = "#7897AB", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Filtered Clonal Distances (Excluding Outliers)")
df_edges <- get_clonal_edges(pn_clones)

dist_threshold <- threshold

df_edges_filtered <- df_edges %>%
  filter(distance <= dist_threshold)
cell_sub<-unique(c(df_edges_filtered$cell1,df_edges_filtered$cell2))
length(cell_sub)
clone_df<-pn_clones %>% filter(cell %in% cell_sub) %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2) %>% as.data.frame()
rownames(clone_df)<-clone_df$cell
head(clone_df)
object<-da_adj
section<-c('S3000_50cs3')
group<-'final_ct3'
pt_size=2
p_list<-list()
df_tmp<-clone_df
for (ident in section){
    img <- object@images[[ident]]@image
    img_w <- ncol(img); img_h <- nrow(img)
    img_grob <- grid::rasterGrob(img, interpolate = TRUE)
    df<-df_tmp #%>% filter(orig.ident==ident)
    p_list[[ident]] <- ggplot(df,aes(x = x_ad,y = y_ad))+ 
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
    geom_point(data = cell_region_ad %>% filter(group=='S3000_50cs3',grepl('Ctx_',L2_region) | grepl('OLF',L2_region)), 
             aes(x = x_px, y = y_px, color = L2_region), 
             size = pt_size * 0.5, 
             alpha = 0.2,         
             stroke = 0) +
      scale_color_manual(values = col_r2) + 
      geom_path(aes(group = clone), alpha = 0.8, size = 0.3,linewidth = 0.5,color='white') +
      geom_point(size=pt_size,shape=21,stroke=0.1,color='white',aes_string(fill=group))+
      scale_fill_manual(values = col_ct3)+
      xlab(paste0("")) +
      ylab(paste0("")) + 
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
      theme_void() +
      labs(title = ident)+    
      theme(
      legend.key.height = unit(0.001, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8),
    legend.position = "right"
    )
}
options(repr.plot.width=10,repr.plot.height=10)
wrap_plots(p_list)
ggsave('PN_fate.spatial_clone_final_ct3.cs3_plot.png',width = 10,height = 10,dpi = 500)
write.csv(clone_df,'PN.clones.outlier_removed.pos_ad.meta.csv')
pdf('PN.clone.upset.plot.pdf',width = 5,height = 5.5)
options(repr.plot.width=5,repr.plot.height=5.5)
#glut_clones$clone<-glut_clones$clone.id2
get_upset_plot(clone_df,group = 'final_ct3',only_confident_class = FALSE)
dev.off()

#### thalamus neuron
glut_clones<-meta %>% filter(fate %in% c('GLUT','GABA'),mouse=='A13') %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2)
glut_clones
### remove outlier clone pairs by spatial distance
dist_results <- glut_clones %>%
  group_by(clone) %>%
  filter(n() >= 2) %>%
  do({
    coords <- as.matrix(.[, c("x_ad", "y_ad")])
    dists <- as.numeric(dist(coords))
    data.frame(distance = dists)
  }) %>%
  ungroup()

head(dist_results)
options(repr.plot.width=5,repr.plot.height=5)
ggplot(dist_results, aes(x = distance)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#7897AB", alpha = 0.6) +
  geom_density(color = "#D885A3", size = 1) +
  stat_ecdf(geom = "step", color = "#98BA7D", size = 0.8) +
  theme_minimal() +
  labs(title = "Distribution of Clonal Inter-cell Distances",
       x = "Distance (units)", y = "Density / Cumulative Prob")

threshold <- quantile(dist_results$distance, 0.8) 
dist_filtered <- dist_results %>%
  filter(distance <= threshold)

ggplot(dist_filtered, aes(x = distance)) +
  geom_density(fill = "#7897AB", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Filtered Clonal Distances (Excluding Outliers)")
df_edges <- get_clonal_edges(glut_clones)

dist_threshold <- threshold#quantile(df_edges$distance, 0.95)
df_edges_filtered <- df_edges %>%
  filter(distance <= dist_threshold)

cell_sub<-unique(c(df_edges_filtered$cell1,df_edges_filtered$cell2))
length(cell_sub)
clone_df<-glut_clones %>% filter(cell %in% cell_sub) %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2) %>% as.data.frame()
rownames(clone_df)<-clone_df$cell
head(clone_df)
### adjust inter TH_3 & TH_6 cells on the boundary due to low resolution mapping error
ad_cells<-clone_df %>% filter(L2_region=='TH_3',y_ad<470) %>% pull(cell)
ad_cells
clone_df_ad<-clone_df
clone_df_ad[ad_cells,'L2_region']<-'TH_6'
cell_region_ad<-cell_region
cell_region_ad[ad_cells,'L2_region']<-'TH_6'
## 260317 update thalamus region
map_ad<-paste0('TH_',c(1:9))
names(map_ad)<-c('RT','LD','MD','ZI','Hb','AM','RE','AV','VA/VL/VM')
clone_df_ad$L2_region_sm2<-as.character(names(map_ad)[match(clone_df_ad$L3_region,map_ad)])
cell_region$L2_region_sm2<-as.character(names(map_ad)[match(cell_region$L3_region,map_ad)])
object<-da_adj
section<-c('S3000_50cs3')
group<-'final_ct3'
pt_size=2
p_list<-list()
df_tmp<-clone_df_ad
for (ident in section){
    img <- object@images[[ident]]@image
    img_w <- ncol(img); img_h <- nrow(img)
    img_grob <- grid::rasterGrob(img, interpolate = TRUE)
    df<-df_tmp #%>% filter(orig.ident==ident)
    p_list[[ident]] <- ggplot(df,aes(x = x_ad,y = y_ad))+ 
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
      geom_point(data = cell_region %>% filter(group=='S3000_50cs3',L2_region_sm2 %in% names(map_ad)), 
             aes(x = x_px, y = y_px, color = L2_region_sm2), 
             size = pt_size * 0.5, 
             alpha = 0.2,         
             stroke = 0) +
      scale_color_manual(values = get_specific_group_color('L2_region')) + 
    
      geom_path(aes(group = clone), alpha = 0.8, size = 0.3,linewidth = 0.5,color='white') +
      geom_point(size=pt_size,shape=21,stroke=0.1,color='white',aes_string(fill=group))+
      #scale_color_manual(values = col_r2)+
      scale_fill_manual(values = get_specific_group_color('ct3'))+
      xlab(paste0("")) +
      ylab(paste0("")) + 
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
      theme_void() +
      labs(title = ident)+    
      theme(
      legend.key.height = unit(0.001, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8),
    legend.position = "right"
    )
}
options(repr.plot.width=10,repr.plot.height=10)
wrap_plots(p_list)
ggsave('GLUT_GABA_fate.spatial_clone_final_ct3.cs3_plot.png',width = 10,height = 10,dpi = 500)

pdf('TH_Nucleus.clone.upsetplot.pdf',width = 3.5,height = 4)
#pdf('TH_Nucleus.clone.upsetplot.pdf',width = 3.5,height = 4)
options(repr.plot.width=3.5,repr.plot.height=4)
get_upset_plot(clone_df_ad,group = 'nucleus',only_confident_class = FALSE)
dev.off()

tmp_dorsal<-clone_df_ad %>% filter(nucleus=='Dorsal') %>% group_by(clone) %>% 
    mutate(type=ifelse(n_distinct(L2_region)==1,'Uniq','Shared')) %>% 
    ungroup() %>% select(clone,nucleus,type) %>% distinct() %>% 
      mutate(tot_cl = n_distinct(clone)) %>% group_by(type) %>%
      mutate(cl_num = n_distinct(clone), prop = cl_num / tot_cl) %>% ungroup() %>%
    select(nucleus,type,cl_num,prop) %>% distinct()
tmp_dorsal

tmp_vm<-clone_df_ad %>% filter(nucleus=='Ventro&Medial') %>% group_by(clone) %>% 
    mutate(type=ifelse(n_distinct(L2_region)==1,'Uniq','Shared')) %>% 
    ungroup() %>% select(clone,nucleus,type) %>% distinct() %>% 
      mutate(tot_cl = n_distinct(clone)) %>% group_by(type) %>%
      mutate(cl_num = n_distinct(clone), prop = cl_num / tot_cl) %>% ungroup() %>%
    select(nucleus,type,cl_num,prop) %>% distinct()
tmp_vm

tmp_plot<-rbind(tmp_dorsal,tmp_vm)
tmp_plot
options(repr.plot.width=2.5,repr.plot.height=3)
ggplot(tmp_plot,aes(fill = type, y = prop, x = nucleus)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = get_specific_group_color('div_col')[c(2,5)]) +
  theme_minimal()
ggsave('TH_nucleus_subregion.clone_shared_prop.barplot.png',width = 2.5,height = 3,dpi = 300)

#### SPN
spn_clones<-meta %>% filter(fate=='SPN',mouse=='A13') %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2)
spn_clones
### remove outlier clone pairs by spatial distance
dist_results <- spn_clones %>%
  group_by(clone) %>%
  filter(n() >= 2) %>%
  do({
    coords <- as.matrix(.[, c("x_ad", "y_ad")])
    dists <- as.numeric(dist(coords))
    data.frame(distance = dists)
  }) %>%
  ungroup()

head(dist_results)
options(repr.plot.width=5,repr.plot.height=5)
ggplot(dist_results, aes(x = distance)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#7897AB", alpha = 0.6) +
  geom_density(color = "#D885A3", size = 1) +
  stat_ecdf(geom = "step", color = "#98BA7D", size = 0.8) +
  theme_minimal() +
  labs(title = "Distribution of Clonal Inter-cell Distances",
       x = "Distance (units)", y = "Density / Cumulative Prob")


threshold <- quantile(dist_results$distance, 0.98) 

dist_filtered <- dist_results %>%
  filter(distance <= threshold)

ggplot(dist_filtered, aes(x = distance)) +
  geom_density(fill = "#7897AB", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Filtered Clonal Distances (Excluding Outliers)")
df_edges <- get_clonal_edges(spn_clones)

dist_threshold <- threshold#quantile(df_edges$distance, 0.95)

df_edges_filtered <- df_edges %>%
  filter(distance <= dist_threshold)

cell_sub<-unique(c(df_edges_filtered$cell1,df_edges_filtered$cell2))
length(cell_sub)
clone_df<-spn_clones %>% filter(cell %in% cell_sub,y_ad<500) %>% group_by(clone) %>% 
    mutate(clone_size=n_distinct(cell)) %>% ungroup() %>% filter(clone_size>=2) %>% as.data.frame()
rownames(clone_df)<-clone_df$cell
head(clone_df)

object<-da_adj
section<-c('S3000_50cs3')
group<-'final_ct3'
pt_size=2
p_list<-list()
df_tmp<-clone_df
for (ident in section){
    img <- object@images[[ident]]@image
    img_w <- ncol(img); img_h <- nrow(img)
    img_grob <- grid::rasterGrob(img, interpolate = TRUE)
    df<-df_tmp #%>% filter(orig.ident==ident)
    p_list[[ident]] <- ggplot(df,aes(x = x_ad,y = y_ad))+ 
      # 1. 底层 HE
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
    # alpha 设置低一些，避免遮挡点；size 设置细线
      geom_path(aes(group = clone), alpha = 0.8, size = 0.3,linewidth = 0.5,color='white') +
      geom_point(size=pt_size,shape=21,stroke=0.1,color='white',aes_string(fill=group))+
      scale_fill_manual(values = col_ct3)+
      xlab(paste0("")) +
      ylab(paste0("")) + 
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
      theme_void() +
      labs(title = ident)+    
      theme(
      legend.key.height = unit(0.001, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8),
    legend.position = "right"
    )
}
options(repr.plot.width=10,repr.plot.height=10)
wrap_plots(p_list)
#ggsave('GLUT_GABA_fate.spatial_clone_final_ct3.cs3_plot.png',width = 10,height = 10,dpi = 500)

##Check SPN subtype clonal relationship
da<-FindClusters(da,resolution = 1)
da_spn<-da[,da$seurat_clusters %in% c(4,28)]
da_spn<-FindClusters(da_spn,resolution = 0.1)
da@meta.data[,c('UMAP_1','UMAP_2')]<-da@reductions$umap@cell.embeddings
da_spn<-da[,da$seurat_clusters %in% c(4,21,28) & da$UMAP_2 < -7.5]
da_spn<-FindClusters(da_spn,resolution = 0.25)
da$clone_tmp<-clone_df[colnames(da),'clone']
plot_df_all <- da@meta.data %>%
    filter(final_ct3 == 'SPN', seurat_clusters %in% c(4, 21, 28), UMAP_2 < -7.5)
df_bg <- plot_df_all %>% filter(is.na(clone_tmp))
df_fg <- plot_df_all %>% filter(!is.na(clone_tmp),seurat_clusters!=21)
write.csv(plot_df_all,'SPN.subtype.clone.meta.csv')
spn_col<-get_specific_group_color('div_col')[c(2,3)]
names(spn_col)<-c('4','28')
spn_col
options(repr.plot.width=4.5,repr.plot.height=5.5)
ggplot() +
    geom_point(data = df_bg, aes(x = UMAP_1, y = UMAP_2),
                size = 2, shape = 21, fill = "grey90", color = "white", stroke = 0.1) +
    geom_path(data = df_fg %>% arrange(clone_tmp, UMAP_1),
            aes(x = UMAP_1, y = UMAP_2, group = clone_tmp),
            alpha = 0.8, linewidth = 0.3, color = '#7897AB') +
    geom_point(data = df_fg, aes(x = UMAP_1, y = UMAP_2, fill = as.factor(seurat_clusters)),
            size = 4, shape = 21, stroke = 0.1, color = 'white') +
    scale_fill_manual(values = spn_col, name = "Cluster") +
    coord_fixed() +
    theme_minimal() +
    theme(panel.background = element_blank(),panel.grid = element_blank(),legend.text = element_text(size = 15),
          axis.text = element_text(size = 15),axis.title = element_text(size = 15), 
          axis.line.x.bottom = element_line(linewidth = 0.25),
          axis.line.y.left = element_line(linewidth = 0.25),legend.position = "right")
#ggsave('SPN_subtype.clone_link.UMAP_subset.png',width = 4.5,height = 5.5,dpi = 500)

df_fg$clone<-df_fg$clone_tmp
get_upset_plot(df_fg,group = 'seurat_clusters',only_confident_class = FALSE)
tmp<-df_fg %>% group_by(clone_tmp) %>% 
    mutate(type=ifelse(n_distinct(seurat_clusters)==1,'Uniq','Shared')) %>% ungroup()
tmp_foxp1<-tmp %>% filter(seurat_clusters==4) %>% select(clone_tmp,type) %>% distinct() %>% 
      mutate(tot_cl = n_distinct(clone_tmp),cluster='Foxp1') %>% group_by(type) %>%
      mutate(cl_num = n_distinct(clone_tmp), prop = cl_num / tot_cl) %>% ungroup() %>% select(cluster,type,cl_num,prop) %>% distinct()
tmp_foxp1
tmp_foxp2<-tmp %>% filter(seurat_clusters==28) %>% select(clone_tmp,type) %>% distinct() %>% 
      mutate(tot_cl = n_distinct(clone_tmp),cluster='Foxp2') %>% group_by(type) %>%
      mutate(cl_num = n_distinct(clone_tmp), prop = cl_num / tot_cl) %>% ungroup() %>% select(cluster,type,cl_num,prop) %>% distinct()
tmp_foxp2

tmp_plot<-rbind(tmp_foxp1,tmp_foxp2)
tmp_plot
options(repr.plot.width=2.5,repr.plot.height=3)
ggplot(tmp_plot,aes(fill = type, y = prop, x = cluster)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = get_specific_group_color('div_col')[c(2,5)]) +
  theme_minimal()
ggsave('SPN_subtype.clone_shared_prop.barplot.png',width = 2.5,height = 3,dpi = 300)

tmp_df<-data.frame(x=c(rep(1,2)),type=c('Shared','Uniq'),cl_num=c(5,40))
tmp_df
options(repr.plot.width=2,repr.plot.height=3)
ggplot(tmp_df,aes(fill = type, y = cl_num,x=x)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = get_specific_group_color('div_col')[c(2,5)]) +
  theme_minimal()
ggsave('SPN_subtype.clone_shared_uniq.barplot.png',width = 2,height = 3,dpi = 300)

genePlot(da,gene = c('Foxp1'))
ggsave('Foxp1.expression.alldata.UMAP.png',width = 8,height = 8,dpi = 400)
genePlot(da,gene = c('Foxp2'))
ggsave('Foxp2.expression.alldata.UMAP.png',width = 8,height = 8,dpi = 400)

object<-da_adj
section<-c('S3000_50cs3')
group<-'cluster'
pt_size=2
p_list<-list()
df_tmp<-test %>% filter(cluster %in% c(4,28))
for (ident in section){
    img <- object@images[[ident]]@image
    img_w <- ncol(img); img_h <- nrow(img)
    img_grob <- grid::rasterGrob(img, interpolate = TRUE)
    df<-df_tmp #%>% filter(orig.ident==ident)
    p_list[[ident]] <- ggplot(df,aes(x = x_ad,y = y_ad))+ 
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
      geom_path(aes(group = clone), alpha = 0.8, size = 0.3,linewidth = 0.5,color='white') +
      geom_point(size=pt_size,shape=21,stroke=0.1,color='white',aes_string(fill=group))+
      scale_fill_manual(values = spn_col)+
      xlab(paste0("")) +
      ylab(paste0("")) + 
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
      theme_void() +
      labs(title = ident)+    
      theme(
      legend.key.height = unit(0.001, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = 8),
    legend.position = "right"
    )
}
options(repr.plot.width=10,repr.plot.height=10)
wrap_plots(p_list)
ggsave('SPN_fate.spatial_clone_final_ct3.cs3_plot.png',width = 10,height = 10,dpi = 500)

### confirm by published scLT data
da_sc<-readRDS("~/NPC_project/scRNAseq/all_sc_integrate/1_celltype_integrated/test6a.final_ct123_clone.rds")
da_sc
genes_of_interest <- c("Foxp1", "Foxp2")  
cell_groups <- c("SPN_1", "SPN_2")    

plot_data <- FetchData(
  object = da_sc[,da_sc$study=='Bandler_2022'],
  vars = c(genes_of_interest, "final_ct3"),  
  slot = "scale.data" 
)

plot_data <- plot_data[plot_data$final_ct3 %in% cell_groups, ]
head(plot_data)
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_of_interest),  
    names_to = "Gene",                 
    values_to = "Expression"           
  )

plot_data_long$final_ct3 <- factor(plot_data_long$final_ct3, levels = cell_groups)
plot_data_long$Gene <- factor(plot_data_long$Gene, levels = genes_of_interest)

options(repr.plot.width=3.5,repr.plot.height=4)
ggplot(plot_data_long, aes(x = final_ct3, y = Expression, fill = final_ct3)) +
  geom_violin(trim = FALSE, scale = "width") + 
  facet_wrap(~ Gene, scales = "free_y") +  
  theme_bw() +
  labs(
    x = "Cell Group",
    y = "Log-normalized Expression",
    title = "Gene Expression by Cell Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c('#EBBD7EB2','#BBA1D4B2')) 
#ggsave('Bandler.SPN_subtype.Foxp12.VlnPlot.png',width = 3.5,height = 4,dpi = 400)

tmp_df<-da_sc@meta.data
tmp_df$clone<-tmp_df$final_clone
pdf('Bandler_2022.e12_e13.SPN_subtype.clone_upset.pdf',width = 3.5,height = 4)
options(repr.plot.width=3.5,repr.plot.height=4)
get_upset_plot(tmp_df %>% filter(!is.na(clone),study=='Bandler_2022', final_ct2=='SPN', stage %in% c('e12','e13')
                                ),group = 'final_ct3',only_confident_class = FALSE)
dev.off()

