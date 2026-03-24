library(Matrix)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(grid)
library(purrr)
library(princurve)
library(igraph)
Sys.setenv(V8_NO_INIT = "true")
options(V8.disabled = TRUE)
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")

wd<-"~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/Merged_2512/5_SpatialFate"
if(!dir.exists(wd)){
    dir.create(wd)
}
setwd(wd)
### Hippocampus
da_adj<-readRDS("../3_clone/test7a.adjusted.final_ct23.region_L123.rds")
lr_df<-read.csv("../2_region/L9/test3.integration.spatial.region.adj_img_pos.LR.meta.csv",row.names = 1)
head(lr_df)
map_df<-read.table("../2_region/L9/cell_L1_L9_region.txt",sep = '\t',header = FALSE)
colnames(map_df)<-c("L1","L9",'L1_region','L2_region','L3_region','cell')
map_df<-map_df %>% select('L9','cell') %>% distinct()
#df1<-as.data.frame(apply(df1,2,FUN = function(x){paste0("S3000_",x)}))
head(map_df)
map_df$LR<-lr_df$LR[match(map_df$L9,rownames(lr_df))]
cell_region<-read.table("../2_region/L9/cell_region.txt",sep = '\t')
colnames(cell_region)<-c("barcode","L1_region","L2_region","L3_region")
rownames(cell_region)<-cell_region$barcode
head(cell_region)
pos_df1<-read.table('~/NPC_project/BMK/tissue/brain/P5/S3000/RNA/BST/07.CellSplit/mtx/cells_center.txt',
                    sep = '\t',row.names = 1,col.names = c('x.axis','y.axis'))
rownames(pos_df1)<-paste0('S3000_',rownames(pos_df1))
pos_df1$group<-'S3000'
head(pos_df1)
pos_df2<-read.table('~/NPC_project/BMK/tissue/brain/P5/S3000_50cs3/RNA/BST/07.CellSplit/mtx/cells_center.txt',
                    sep = '\t',row.names = 1,col.names = c('x.axis','y.axis'))
rownames(pos_df2)<-paste0('S3000_50cs3_',rownames(pos_df2))
pos_df2$group<-'S3000_50cs3'
head(pos_df2)
pos_df<-rbind(pos_df1,pos_df2)
head(pos_df)
cell<-intersect(rownames(cell_region),rownames(pos_df))
pos_df<-pos_df[cell,]
pos_df[,c('L2_region','L3_region','LR')]<-cell_region[cell,c('L2_region','L3_region','LR')]
head(pos_df)

#S3000
img1<-da_adj@images$S3000@image
img_w1 <- ncol(img1)
img_h1 <- nrow(img1)
img_grob1 <- rasterGrob(img1, interpolate = TRUE)
plot_data1 <- align_spots(pos_df %>% filter(group=='S3000'), img_h = img_h1,img_w = img_w1,
                         scale_x = 0.935, scale_y=0.83,
                         off_x = 65, off_y = 90, rotate = 270,flip_y = TRUE)
head(plot_data1)
ggplot(plot_data1, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob1, xmin = 0, xmax = img_w1, ymin = 0, ymax = img_h1) +
  geom_point(color = "cyan", size = 0.5, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w1), ylim = c(0, img_h1)) +
  theme_void() +
  labs(subtitle = "")

#S3000_50cs3
img2<-da_adj@images$S3000_50cs3@image
img_w2 <- ncol(img2)
img_h2 <- nrow(img2)
img_grob2 <- rasterGrob(img2, interpolate = TRUE)
plot_data2 <- align_spots(pos_df %>% filter(group=='S3000_50cs3'), img_h = img_h2,img_w = img_w2,
                         scale_x = 1, scale_y=0.87,
                         off_x = 0, off_y = 6, rotate = 90,flip_y = TRUE)
head(plot_data2)
ggplot(plot_data2, aes(x = x_px, y = y_px)) +
  annotation_custom(img_grob2, xmin = 0, xmax = img_w2, ymin = 0, ymax = img_h2) +
  geom_point(color = "cyan", size = 1, alpha = 0.3) + 
  coord_fixed(xlim = c(0, img_w2), ylim = c(0, img_h2)) +
  theme_void() +
  labs(subtitle = "")

plot_data<-rbind(plot_data1,plot_data2)
pos_df[,c('x_px','y_px')]<-plot_data[rownames(pos_df),c('x_px','y_px')]
write.csv(pos_df,'A13.all.cell.region.LR.pos_ad.meta.csv')

ggplot(pos_df %>% filter(group=='S3000',L3_region=='DG_1'), aes(x = x_px, y = y_px,color=L3_region)) +
  annotation_custom(img_grob1, xmin = 0, xmax = img_w1, ymin = 0, ymax = img_h1) +
  geom_point(size = 0.25, alpha = 0.3) + 
  scale_color_manual(values = col_r3) +
  coord_fixed(xlim = c(0, img_w1), ylim = c(0, img_h1)) +
  theme_void() +
  labs(subtitle = "")

ggplot(pos_df %>% filter(group=='S3000_50cs3',L3_region=='DG_1'), aes(x = x_px, y = y_px,color=L3_region)) +
  annotation_custom(img_grob2, xmin = 0, xmax = img_w2, ymin = 0, ymax = img_h2) +
  geom_point(size = 0.25, alpha = 0.3) + 
  scale_color_manual(values = col_r3) +
  coord_fixed(xlim = c(0, img_w2), ylim = c(0, img_h2)) +
  theme_void()

### CA1
pos_df$group2<-paste0(pos_df$group,"-",pos_df$LR)
head(pos_df)
ca1_results <- pos_df %>%
    filter(L3_region == "CA1_3") #%>%
ca1_results_list <-split(ca1_results,ca1_results$group2) %>%
    map(function(sub_df) {
        current_group <- as.character(unique(sub_df$group2))
      
        points_df <- sub_df %>% select(x_px, y_px)
        
        hull_df <- concaveman(as.matrix(points_df), concavity = 3) %>% as.data.frame()
        colnames(hull_df) <- c("x_px", "y_px")
        axis<-get_precise_axis(points_df)
        
        return(list(
        points = points_df,
        hull = hull_df,
        axis = axis
        ))
    })
saveRDS(ca1_results_list,'CA1.LR.mask.axis.rds')
ca_clones<-read.csv('CA_DG.clones.outlier_removed.pos_ad.meta.csv',row.names = 1)
ca_clones$LR<-cell_region$LR[match(rownames(ca_clones),cell_region$barcode)]
ca_clones[,c('x_px','y_px')]<-da_adj@meta.data[rownames(ca_clones),c('x_px','y_px')]
ca_clones$group2<-paste0(ca_clones$orig.ident,'-',ca_clones$LR)
head(ca_clones)
ca1_clone_angle_results_list <- imap(ca1_results_list, function(spatial_info, group_name) {
    
    current_clones <- ca_clones %>%
        filter(group2 == group_name, final_ct3 == 'PN_CA1') %>%
        group_by(clone) %>%
        filter(n() >= 2) %>%
        ungroup()
    
    current_axis <- spatial_info$axis
    
    res_angles <- compute_pair_angles(current_clones, current_axis)
    
    return(data.frame(
        angle = res_angles,
        group = group_name
    ))
})
final_angle_df <- bind_rows(ca1_clone_angle_results_list)
final_angle_df$section<-str_split(final_angle_df$group,'-',simplify = TRUE)[,1]
head(final_angle_df)
unique(final_angle_df$section)
write.csv(final_angle_df,'PN_CA1.clone_cell_pair.angle_to_axis.csv')
plot_data <- final_angle_df %>%
  mutate(angle_capped = ifelse(angle > 45, 47.5, angle)) 
ggplot(plot_data, aes(x = angle_capped)) +
  geom_histogram(
    binwidth = 5, 
    boundary = 0, 
    fill = get_specific_group_color('ct3')['PN_CA1'], 
    color = "white", 
    alpha = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(2.5, 47.5, by = 5), 
    labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40","40-45", ">45")
  ) +
  coord_cartesian(xlim = c(0, 50)) + 
  labs(
    title = "PN_CA1 Clone Orientation",
    subtitle = "5° Bins | Angles >45° Merged",
    x = "Angle Intervals (Degrees)",
    y = "Count of Pairs"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.x = element_line(color = "black")
  )
ggsave('PN_CA1.clone_cell_pair_orientation_to_axis.adjusted.png',width = 5,height = 6,dpi = 300)

### CA3
ca3_results <- pos_df %>%
    filter(L3_region == "CA3_2") #%>%
ca3_results_list <-split(ca3_results,ca3_results$group2) %>%
    map(function(sub_df) {
        current_group <- as.character(unique(sub_df$group2))
        points_df <- sub_df %>% select(x_px, y_px)
        hull_df <- concaveman(as.matrix(points_df), concavity = 3) %>% as.data.frame()
        colnames(hull_df) <- c("x_px", "y_px")
        axis<-get_precise_axis(points_df)
        return(list(
        points = points_df,
        hull = hull_df,
        axis = axis
        ))
    })
saveRDS(ca3_results_list,'CA3.LR.mask.axis.rds')
ca3_clone_angle_results_list <- imap(ca3_results_list, function(spatial_info, group_name) {
    
    current_clones <- ca_clones %>%
        filter(group2 == group_name, final_ct3 == 'PN_CA3') %>%
        group_by(clone) %>%
        filter(n() >= 2) %>%
        ungroup()
    
    current_axis <- spatial_info$axis
    
    res_angles <- compute_pair_angles(current_clones, current_axis)
    
    return(data.frame(
        angle = res_angles,
        group = group_name
    ))
})

final_angle_df2 <- bind_rows(ca3_clone_angle_results_list)
final_angle_df2$section<-str_split(final_angle_df2$group,'-',simplify = TRUE)[,1]
head(final_angle_df2)
unique(final_angle_df2$section)
write.csv(final_angle_df2,'PN_CA3.clone_cell_pair.angle_to_axis.csv')
plot_data <- final_angle_df2 %>%
  filter(!is.na(angle)) %>% 
  mutate(angle_capped = ifelse(angle > 45, 47.5, angle))

ggplot(plot_data, aes(x = angle_capped)) +
  geom_histogram(
    binwidth = 5, 
    boundary = 0, 
    fill = get_specific_group_color('ct3')['PN_CA3'], 
    color = "white", 
    alpha = 0.8
  ) +

  scale_x_continuous(
    breaks = seq(2.5, 47.5, by = 5), 
    labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40","40-45", ">45")
  ) +
  coord_cartesian(xlim = c(0, 50)) + 
  
  labs(
    title = "PN_CA3 Clone Orientation",
    subtitle = paste0("5° Bins | n = ", nrow(plot_data), " pairs | >45° Merged"),
    x = "Angle Intervals (Degrees)",
    y = "Count of Pairs"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.x = element_line(color = "black")
  )
ggsave('PN_CA3.clone_cell_pair_orientation_to_axis.adjusted.png',width = 5,height = 6,dpi = 300)

### DG
get_precise_long_axis_v12 <- function(df) {
    mat <- as.matrix(unique(df[, c("x_px", "y_px")]))
    if(nrow(mat) < 15) return(NULL)
    dists <- as.matrix(dist(mat))
    g <- mst(graph_from_adjacency_matrix(dists, weighted = TRUE, mode = "undirected"))
    v_ends <- farthest_vertices(g)$vertices
    path_nodes <- as.numeric(shortest_paths(g, from = v_ends[1], to = v_ends[2])$vpath[[1]])
    path_mat <- mat[path_nodes, ]
    pull_radius <- 20 
    
    corrected_path <- path_mat
    for(i in 1:nrow(path_mat)) {
        curr_pt <- path_mat[i, ]
        d_to_cells <- sqrt((mat[,1] - curr_pt[1])^2 + (mat[,2] - curr_pt[2])^2)
        neighbors <- mat[d_to_cells <= pull_radius, , drop = FALSE]
        
        if(nrow(neighbors) > 5) {
            corrected_path[i, 1] <- mean(neighbors[, 1])
            corrected_path[i, 2] <- mean(neighbors[, 2])
        }
    }
    n_nodes <- nrow(corrected_path)
    trim_size <- max(2, floor(n_nodes * 0.04))
    corrected_path <- corrected_path[trim_size:(n_nodes - trim_size), ]
    t_seq <- 1:nrow(corrected_path)
    fit_x <- smooth.spline(t_seq, corrected_path[,1], spar = 0.55, df = 20)
    fit_y <- smooth.spline(t_seq, corrected_path[,2], spar = 0.55, df = 20)
    res <- data.frame(
        x_px = predict(fit_x, seq(1, nrow(corrected_path), length.out = 200))$y,
        y_px = predict(fit_y, seq(1, nrow(corrected_path), length.out = 200))$y
    )
    return(res)
}
dg_results <- pos_df %>%
    filter(L3_region == "DG_1") #%>%
dg_results_list <-split(dg_results,dg_results$group2) %>%
    map(function(sub_df) {
        current_group <- as.character(unique(sub_df$group2))
        points_df <- sub_df %>% select(x_px, y_px)
        hull_edges <- get_v_mask(points_df, alpha_val = 8) 
        axis <- get_precise_long_axis_v12(points_df)
        axis$group2 <- current_group
        return(list(
            points = points_df,
            hull = hull_edges, 
            axis = axis
        ))
    })

dg_clone_angle_results_list <- imap(dg_results_list, function(spatial_info, group_name) {
    current_clones <- ca_clones %>%
        filter(group2 == group_name, final_ct3 == 'DG') %>%
        group_by(clone) %>%
        filter(n() >= 2) %>%
        ungroup()
    current_axis <- spatial_info$axis
    res_angles <- compute_pair_angles(current_clones, current_axis)
    return(data.frame(
        angle = res_angles,
        group = group_name
    ))
})
final_angle_df3 <- bind_rows(dg_clone_angle_results_list)
final_angle_df3$section<-str_split(final_angle_df3$group,'-',simplify = TRUE)[,1]
head(final_angle_df3)
unique(final_angle_df3$section)
write.csv(final_angle_df3,'DG.clone_cell_pair.angle_to_axis.csv')

plot_data <- final_angle_df3 %>%
  mutate(angle_capped = ifelse(angle > 45, 47.5, angle)) 
ggplot(plot_data, aes(x = angle_capped)) +
  geom_histogram(
    binwidth = 5, 
    boundary = 0, 
    fill = get_specific_group_color('ct3')['DG'], 
    color = "white", 
    alpha = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(2.5, 47.5, by = 5), 
    labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40","40-45", ">45")
  ) +
  coord_cartesian(xlim = c(0, 50)) + 
  
  labs(
    title = "DG Clone Orientation",
    subtitle = "5° Bins | Angles >45° Merged",
    x = "Angle Intervals (Degrees)",
    y = "Count of Pairs"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.x = element_line(color = "black")
  )
ggsave('DG.clone_cell_pair_orientation_to_axis.adjusted.png',width = 5,height = 6,dpi = 300)

ggplot() +
    annotation_custom(img_grob1, xmin = 0, xmax = img_w1, ymin = 0, ymax = img_h1) +
    geom_polygon(data = ca1_results_list[['S3000-L']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA1_3'], alpha = 0.4, color = 'white', linewidth = 0.25) +
    geom_path(data = ca1_results_list[['S3000-L']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_polygon(data = ca1_results_list[['S3000-R']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA1_3'], alpha = 0.4, color = "white", linewidth = 0.25) +
    geom_path(data = ca1_results_list[['S3000-R']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_segment(data = dg_results_list[['S3000-L']]$hull,
                 aes(x = x1, y = y1, xend = x2, yend = y2), 
                 color = col_r3['DG_1'], linewidth = 0.5) +
    #CA3
    geom_polygon(data = ca3_results_list[['S3000-L']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA3_2'], alpha = 0.4, color = 'white', linewidth = 0.25) +
    geom_path(data = ca3_results_list[['S3000-L']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_polygon(data = ca3_results_list[['S3000-R']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA3_2'], alpha = 0.4, color = "white", linewidth = 0.25) +
    geom_path(data = ca3_results_list[['S3000-R']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_path(data = dg_results_list[['S3000-L']]$axis, 
              aes(x = x_px, y = y_px), 
              color = "white", linewidth = 0.5) +
    geom_segment(data = dg_results_list[['S3000-R']]$hull, 
                 aes(x = x1, y = y1, xend = x2, yend = y2), 
                 color = col_r3['DG_1'], linewidth = 0.5) +
    geom_path(data = dg_results_list[['S3000-R']]$axis, 
              aes(x = x_px, y = y_px), 
              color = "white", linewidth = 0.5) +
    coord_fixed(xlim = c(0, img_w1), ylim = c(0, img_h1), expand = FALSE) +
    theme_void()
ggsave('S3000.CA1_CA3_DG.mask.axis.png',width=10,height = 10,dpi = 500)

ggplot() +
    annotation_custom(img_grob2, xmin = 0, xmax = img_w2, ymin = 0, ymax = img_h2) +
    #CA1
    geom_polygon(data = ca1_results_list[['S3000_50cs3-L']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA1_3'], alpha = 0.4, color = 'white', linewidth = 0.25) +
    geom_path(data = ca1_results_list[['S3000_50cs3-L']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_polygon(data = ca1_results_list[['S3000_50cs3-R']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA1_3'], alpha = 0.4, color = "white", linewidth = 0.25) +
    geom_path(data = ca1_results_list[['S3000_50cs3-R']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_segment(data = dg_results_list[['S3000_50cs3-L']]$hull,
                 aes(x = x1, y = y1, xend = x2, yend = y2), 
                 color = col_r3['DG_1'], linewidth = 0.5) +
    #CA3
    geom_polygon(data = ca3_results_list[['S3000_50cs3-L']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA3_2'], alpha = 0.4, color = 'white', linewidth = 0.25) +
    geom_path(data = ca3_results_list[['S3000_50cs3-L']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_polygon(data = ca3_results_list[['S3000_50cs3-R']]$hull, aes(x = x_px, y = y_px),
        fill = col_r3['CA3_2'], alpha = 0.4, color = "white", linewidth = 0.25) +
    geom_path(data = ca3_results_list[['S3000_50cs3-R']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    geom_path(data = dg_results_list[['S3000_50cs3-L']]$axis, 
              aes(x = x_px, y = y_px), 
              color = "white", linewidth = 0.5) +
    geom_segment(data = dg_results_list[['S3000_50cs3-R']]$hull, 
                 aes(x = x1, y = y1, xend = x2, yend = y2), 
                 color = col_r3['DG_1'], linewidth = 0.5) +
    geom_path(data = dg_results_list[['S3000_50cs3-R']]$axis, 
              aes(x = x_px, y = y_px), 
              color = "white", linewidth = 0.5) +
    coord_fixed(xlim = c(0, img_w2), ylim = c(0, img_h2), expand = FALSE) +
    theme_void()
ggsave('S3000_50cs3.CA1_CA3_DG.mask.axis.png',width=10,height = 10,dpi = 500)

### cortical neurons
pos_df<-read.csv('A13.all.cell.region.LR.pos_ad.meta.csv',header = TRUE,row.names = 1)
head(pos_df)
dim(pos_df)
pos_df2<-pos_df %>% filter(L2_region %in% c('Ctx_L23','OLF','Ctx_L6','Ctx_L4','Ctx_L5','Ctx_L1'))
pos_df2$group2<-paste0(pos_df2$group,"-",pos_df2$LR)
head(pos_df2)

ctx_results_list <- split(pos_df2, pos_df2$group2) %>%
  map(function(sub_df) {
    points_df <- sub_df %>% select(x_px, y_px)
    hull_df <- get_cortex_boundary(points_df)
    pca_res <- prcomp(points_df)
    axis <- get_precise_axis(points_df) 
    
    return(list(
      points = points_df,
      hull = hull_df,
      axis = axis
    ))
  })

saveRDS(ctx_results_list,'Ctx.LR.mask.axis.rds')
ggplot() +
    annotation_custom(img_grob1, xmin = 0, xmax = img_w1, ymin = 0, ymax = img_h1) +
    #geom_polygon(data = ctx_results_list[['S3000-L']]$hull, aes(x = x_px, y = y_px),
    #    fill = get_specific_group_color('L2_region')['Ctx_L6'], alpha = 0.4, color = 'white', linewidth = 0.25) +
    geom_path(data = ctx_results_list[['S3000-L']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    #geom_polygon(data = ctx_results_list[['S3000-R']]$hull, aes(x = x_px, y = y_px),
    #    fill = get_specific_group_color('L2_region')['Ctx_L6'], alpha = 0.4, color = "white", linewidth = 0.25) +
    geom_path(data = ctx_results_list[['S3000-R']]$axis, aes(x = x_px, y = y_px), color = "white", linewidth = 0.5) +
    coord_fixed(xlim = c(0, img_w1), ylim = c(0, img_h1), expand = FALSE) +
    theme_void()
#ggsave('S3000.CA1_3.mask_axis.png',width = 10,height = 10,dpi = 500)

clones<-read.csv('PN.clones.outlier_removed.pos_ad.meta.csv',row.names = 1)
clones$LR<-pos_df$LR[match(rownames(clones),rownames(pos_df))]
clones[,c('x_px','y_px')]<-da_adj@meta.data[rownames(clones),c('x_px','y_px')]
clones$group2<-paste0(clones$orig.ident,'-',clones$LR)
head(clones)

clone_angle_results_list <- imap(ctx_results_list, function(spatial_info, group_name) {

    current_axis <- spatial_info$axis
   
    current_clones <- clones %>%
        filter(group2 == group_name) %>%
        group_by(clone) %>%
        filter(n() >= 2) %>%
        ungroup()

    res_angles <- compute_pair_angles(current_clones, current_axis)
    
    return(data.frame(
        angle = res_angles,
        group = group_name
    ))
})

final_angle_df <- bind_rows(clone_angle_results_list)
options(repr.plot.width=5, repr.plot.height=6)

# 
plot_data <- final_angle_df %>%
  filter(!is.na(angle)) %>% 
  mutate(angle_capped = ifelse(angle > 45, 47.5, angle))

ggplot(plot_data, aes(x = angle_capped)) +
  geom_histogram(
    binwidth = 5, 
    boundary = 0, 
    fill = get_specific_group_color('L1_region')['Ctx'], 
    color = "white", 
    alpha = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(2.5, 47.5, by = 5), 
    labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40","40-45", ">45")
  ) +
  coord_cartesian(xlim = c(0, 50)) + 
  
  labs(
    title = "PN Clone Orientation",
    subtitle = paste0("5° Bins | n = ", nrow(plot_data), " pairs | >45° Merged"),
    x = "Angle Intervals (Degrees)",
    y = "Count of Pairs"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line.x = element_line(color = "black")
  )
ggsave('PN.clone_cell_pair_orientation_to_axis.adjusted.png', width = 5, height = 6, dpi = 300)

