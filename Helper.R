library(ggplot2)
plot_spatial_group<-function(object,ident,group,col,pt.size=0.6){
    img<-object@images[[ident]]@image

    img_w <- ncol(img)
    img_h <- nrow(img)
    img_grob <- rasterGrob(img, interpolate = TRUE)
    meta_tmp<-object@meta.data[object$orig.ident==ident,]
    dim(meta_tmp)
    ggplot(meta_tmp, aes_string(x = "x_px", y = "y_px",fill=group)) +
      annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
      geom_point(size = pt.size, alpha = 0.8,shape=21,stroke = 0) + 
      scale_fill_manual(values = col)+
      coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h)) +
      theme_void()
}

library(RANN)
smooth_spatial_clusters <- function(df, cluster_col, k = 6) {
  coords <- df[, c("x_px", "y_px")]
  clusters <- df[[cluster_col]]

  nn <- RANN::nn2(coords, k = k + 1) 
  
  new_clusters <- sapply(1:nrow(coords), function(i) {
    neighbor_indices <- nn$nn.idx[i, ]
    neighbor_labels <- clusters[neighbor_indices]
    
    tab <- table(neighbor_labels)
    return(names(tab)[which.max(tab)])
  })
  
  return(new_clusters)
}

get_most_frequent_value <- function(df, row, group) {  
  x <- df$x.axis[row]  
  y <- df$y.axis[row]  

  nearby_rows <- df[abs(df$x.axis - x) <= 100 & abs(df$y.axis - y) <= 100 & (df[,group] != "Undefined"),]   
  
  most_frequent <- nearby_rows[,group] %>%  
    na.omit() %>%  
    table() %>%  
    which.max()  
  
  return(names(most_frequent))  
} 

align_spots <- function(meta, img_w=NULL,img_h=NULL,scale_x = 0.8, scale_y = 0.9, 
                        off_x = 0, off_y = 0,
                        rotate=0,
                        flip_y = FALSE#,module='Cluster38'
                       ) {
  df<- data.frame(
    row.names = rownames(meta),
    x_raw = meta$x.axis,
    y_raw = meta$y.axis#,
    #score = meta[,module]
  ) %>%
    mutate(
      #
      x_n = (x_raw - min(x_raw)) / (max(x_raw) - min(x_raw)),
      y_n = (y_raw - min(y_raw)) / (max(y_raw) - min(y_raw)),
    )
  if(flip_y) { df$y_n <- 1 - df$y_n } #
  if(rotate==0){
      df <- df %>%
            mutate(
              x_px = (x_n * img_w * scale_x) + off_x,
              y_px = (y_n * img_h * scale_y) + off_y
            )
  } 
  if(rotate==90){
      df <- df %>%
            mutate(
              x_px = ((1-y_n) * img_w * scale_x) + off_x,
              y_px = (x_n * img_h * scale_y) + off_y
            )
  }  
  if(rotate==270){
      df <- df %>%
            mutate(
              x_px = (y_n * img_w * scale_x) + off_x,
              y_px = ((1 - x_n) * img_h * scale_y) + off_y
            )
  }
  return(df)
}

correct_brain_side_manual <- function(df, x_ann, y_ann, angle_deg) {
  angle_rad <- angle_deg * pi / 180

  df <- df %>%
    mutate(
      x_rel = x_px - x_ann,
      y_rel = y_px - y_ann,
      x_rot = x_rel * cos(angle_rad) + y_rel * sin(angle_rad)
    )
  df$side <- ifelse(df$x_rot < 0, "R", "L")
  line_length <- 1000 
  midline_df <- data.frame(
    x = x_ann + c(-line_length * sin(angle_rad), line_length * sin(angle_rad)),
    y = y_ann + c(line_length * cos(angle_rad), -line_length * cos(angle_rad))
  )
  return(list(data = df, line = midline_df))
}

genePlot<-function(object,gene){
    tmp<-as.data.frame(object@reductions$umap@cell.embeddings)
    if(!is.na(match(gene,rownames(object@assays$integrated@data)))){
        tmp[,gene]<-object@assays$integrated@data[gene,]}
    else{
        if(!is.na(match(gene,rownames(object@assays$RNA@data)))){
            tmp[,gene]<-object@assays$RNA@data[gene,]
        }
        else{
            tmp[,gene]<-object@meta.data[,gene]
        }
        }
    names(tmp) <- gsub("-", ".", names(tmp))
    gene<-gsub("-", ".", gene)
    options(repr.plot.width=10,repr.plot.height=10)
    ggplot(tmp, aes(x = UMAP_1, y = UMAP_2)) + 
            geom_point(aes_string(color = gene), alpha = 0.75,size=0.1,shape = 21) + scale_color_gradientn(colours = exp_col)+
            theme_classic()
}

######################
ImageSpacePlot = function(obj, group_by, type="DAPI", sample=names(obj@misc$info)[1], size=1, alpha=1,color=MYCOLOR){
    MYCOLOR=c(
          "#6394ce", "#2a4c87", "#eed500", "#ed5858",
          "#f6cbc2", "#f5a2a2", "#3ca676", "#6cc9d8",
          "#ef4db0", "#992269", "#bcb34a", "#74acf3",
          "#3e275b", "#fbec7e", "#ec4d3d", "#ee807e",
          "#f7bdb5", "#dbdde6", "#f591e1", "#51678c",
          "#2fbcd3", "#80cfc3", "#fbefd1", "#edb8b5",
          "#5678a8", "#2fb290", "#a6b5cd", "#90d1c1",
          "#a4e0ea", "#837fd3", "#5dce8b", "#c5cdd9",
          "#f9e2d6", "#c64ea4", "#b2dfd6", "#dbdfe7",
          "#dff2ec", "#cce8f3", "#e74d51", "#f7c9c4",
          "#f29c81", "#c9e6e0", "#c1c5de", "#750000"
          )
        
        raster_type <- switch(type,
                          HE = "img_he_gg",
                          DAPI = "img_gg",
                          stop("Invalid type. Must be 'HE' or 'DAPI'.")
                             )

        spatial_coord1 <- as.data.frame(obj[[group_by]])
        colnames(spatial_coord1) <- group_by
        spatial_coord2 <- as.data.frame(obj@reductions$spatial@cell.embeddings)
        spatial_coord <-cbind(spatial_coord2,spatial_coord1)
 
    
        ImageSpacePlot <- ggplot2::ggplot() + ggplot2::annotation_custom(grob = obj@misc$info[[sample]][[raster_type]],
        xmin = 0, xmax = obj@misc$info[[sample]]$size_x, 
        ymin = 0, ymax = obj@misc$info[[sample]]$size_y) +
        ggplot2::geom_point(data = spatial_coord, ggplot2::aes(x = spatial_1,y = spatial_2, color = !!sym(group_by), 
                            fill = !!sym(group_by)), size=size, alpha=alpha)+
        labs(size = group_by) + guides(alpha = "none")+ 
        ggplot2::theme_classic()+
        scale_color_manual(values = color)+ coord_fixed()
    return(ImageSpacePlot)
}

###################
FeatureSpacePlot = function(obj, feature, type="DAPI", sample=names(obj@misc$info)[1], size=1, alpha=c(1,1),color=c("lightgrey","blue")){
    raster_type <- switch(type,
                          HE = "img_he_gg",
                          DAPI = "img_gg",
                          stop("Invalid type. Must be 'HE' or 'DAPI'.")
                             )

    spatial_coord1 <- as.data.frame(obj@reductions$spatial@cell.embeddings)
    spatial_coord2 <- FetchData(obj,feature)
    colnames(spatial_coord2) <- feature
    spatial_coord <-cbind(spatial_coord1,spatial_coord2)

    FeatureSpacePlot <-ggplot2::ggplot() + ggplot2::annotation_custom(grob = obj@misc$info[[sample]][[raster_type]],
        xmin = 0, xmax = obj@misc$info[[sample]]$size_x, ymin = 0, ymax = obj@misc$info[[sample]]$size_y) +
        ggplot2::geom_point(data = spatial_coord, ggplot2::aes(x = spatial_1, y = spatial_2,color = !!sym(feature),alpha = !!sym(feature)), size=size)+
        labs(color = feature)+
        guides(alpha = "none")+
        ggplot2::theme_classic()+
        ggplot2::scale_alpha_continuous(range=alpha)+
        scale_color_gradient(low=color[1],high = color[2])+ coord_fixed()
    return(FeatureSpacePlot)
    }

show_cluster_in_b29_spatial<-function(obj1,obj2,obj3,group='seurat_clusters',type='DAPI',all=FALSE){
    cells1<-colnames(obj1)[obj1$study=='B29_s1']
    cells1a<-gsub("B29_left1_","",cells1)
    cells2<-colnames(obj1)[obj1$study=='B29_s2']
    cells2a<-gsub("B29_left2_","",cells2)
    ref1_tmp<-obj2[,cells1a]
    ref1_tmp$tmp<-as.character(obj1@meta.data[cells1,group])
    ref1_tmp<-ref1_tmp[,!is.na(ref1_tmp$tmp)]
    ref2_tmp<-obj3[,cells2a]
    ref2_tmp$tmp<-as.character(obj1@meta.data[cells2,group])
    ref2_tmp<-ref2_tmp[,!is.na(ref2_tmp$tmp)]
    ct_col<-get_specific_group_color("div_col")[1:length(unique(c(ref1_tmp$tmp,ref2_tmp$tmp)))]
    names(ct_col)<-unique(c(ref1_tmp$tmp,ref2_tmp$tmp))
    if(all==FALSE){
        for(i in unique(c(ref1_tmp$tmp,ref2_tmp$tmp))){
            
        # DAPI
        if(i %in% unique(ref1_tmp$tmp)){
            options(repr.plot.height=5, repr.plot.width=15)
            p1<-ImageSpacePlot(obj=ref1_tmp[,ref1_tmp$tmp==i], group_by = "tmp",type=type,size=0.01)+
                scale_color_manual(values = ct_col)
            #scale_color_manual(values = get_specific_group_color("div_col"))
            plot(p1)
        }
        if(i %in% unique(ref2_tmp$tmp)){
            options(repr.plot.height=5, repr.plot.width=15)
            p2<-ImageSpacePlot(obj=ref2_tmp[,ref2_tmp$tmp==i], group_by = "tmp",type=type,size=0.01)+
                scale_color_manual(values = ct_col)
            #scale_color_manual(values = get_specific_group_color("div_col"))
            plot(p2)
            }
        }
    }
    else{
        options(repr.plot.height=10, repr.plot.width=30)
        p1<-ImageSpacePlot(obj=ref1_tmp, group_by = "tmp",type=type,size=0.01)+
            scale_color_manual(values = ct_col)
            #scale_color_manual(values = get_specific_group_color("div_col"))
        plot(p1)
        options(repr.plot.height=10, repr.plot.width=30)
        p2<-ImageSpacePlot(obj=ref2_tmp, group_by = "tmp",type=type,size=0.01)+
            scale_color_manual(values = ct_col)
            #scale_color_manual(values = get_specific_group_color("div_col"))
        plot(p2)
    }
    
}

show_celltype_marker_umap<-function(obj,major_ct){
    if(major_ct=='EN'){
        genes<-c('Slc17a6','Slc17a7','Satb2','Prox1','Rorb','Pcp4','Tle4','Eomes','Crym','Nr4a2')
        for(i in genes){
            p<-genePlot(obj,i)
            plot(p)
        }
    }
    if(major_ct=='IN'){
        genes<-c('Slc32a1','Gad1','Dlx1','Dlx2','Lhx6','Cck','Meis2','Foxp1','Six3','Foxp2','Nkx2-1','Mef2c','Sst','Pvalb')
        for(i in genes){
            p<-genePlot(obj,i)
            plot(p)
        }
    }
    if(major_ct=='precursor'){
        genes<-c('Sox2','Pax6','Ascl1','Neurog2','Emx1','Eomes','Tubb3','Olig1','Dlx1','Dlx2','Neurod1','Aqp4')
        for(i in genes){
            p<-genePlot(obj,i)
            plot(p)
        }
    }
    if(major_ct=='glia'){
        genes<-c('Sox2','Pax6','Gfap','Aldoc','Slc7a10','S100a6','Plp1','Pdgfra','Tmem212')
        for(i in genes){
            p<-genePlot(obj,i)
            plot(p)
        }
    }
    if(major_ct=='other'){
        genes<-c('Ly86','Mrc1','Vtn','Ttr','Cldn5','Igf1','Itgam','Top2a')
        for(i in genes){
            p<-genePlot(obj,i)
            plot(p)
        }
    }

}

##
generate_lowsat_cols<-function(n=10,
                            h=210,         
                            c=15,          
                            l=seq(25, 95, length.out = n)   
                              ){

    cols <- grDevices::hcl(h = h, c = c, l = l, fixup = TRUE, alpha = TRUE)
    
    print(cols)
    return(cols)
    
    library(ggplot2)
    df <- data.frame(x = 1:n,
                     y = rnorm(n),
                     grp = factor(1:n))
    p<-ggplot(df, aes(x, y, fill = grp)) +
      geom_col(width = .7) +
      scale_fill_manual(values = cols)
    plot(p)
}

plot_SpatialDim_across_sections <- function(object, group, col,pt_size=1,saveplot=FALSE,filepath
                                #add_spot
                               ) {
    section<-c('S3000','S3000_50cs3'#,'B33_s1'
              )
    p_list<-list()
    for (ident in section){
        img <- object@images[[ident]]@image
        img_w <- ncol(img); img_h <- nrow(img)
        img_grob <- grid::rasterGrob(img, interpolate = TRUE)
        df<-object@meta.data[object$orig.ident==ident,]
        p_list[[ident]] <- ggplot(df,aes(x = x_px,y = y_px))+ 
          # HE/dapi
          annotation_custom(img_grob, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
          geom_point(size=pt_size,shape=21,stroke=0.1,aes_string(fill=group))+
          scale_fill_manual(values = col)+
          xlab(paste0("")) +
          ylab(paste0("")) + 
          coord_fixed(xlim = c(0, img_w), ylim = c(0, img_h), clip = 'off') +
          theme_void() +
          labs(title = ident)#+
        #guides(fill = guide_legend(ncol=1,reverse = TRUE,override.aes = list(alpha = 1)),
        #       alpha = 'none', 
        #       color = 'none') +
        theme(
          legend.key.height = unit(0.001, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.text = element_text(size = 8),
        legend.position = "right"
        )
    }
    options(repr.plot.width=20,repr.plot.height=7)
    wrap_plots(p_list)
    if(saveplot){
        ggsave(file.path(filepath,paste0(group,".clonecell.spatialplot.png")),width=20,height = 7,dpi = 400)
    }
    #plot(p_list[['S3000_50cs3']])
    #plot(p_list[['B33_s1']])
}

get_domain_celltype_bubble_bar_plot_v2<-function(df,group1,group2,col1,col2,w=12,h=12,mulp=3,saveplot=FALSE,outdir){
    df<-df[,c(group1,group2)]
    df<-df[!is.na(df[,1]) & !is.na(df[,2]),]
    group1<-colnames(df)[1]
    group2<-colnames(df)[2]
    colnames(df)<-c('group1','group2')
    df2<-df %>% group_by(group1,group2) %>% summarise(num=n()) %>% group_by(group2) %>% mutate(pct=num/sum(num))
    df2$group2<-factor(df2$group2,levels = rev(levels(df$group2)))
    df2$group1<-factor(df2$group1,levels = levels(df$group1))
    #df2$domain<-factor(df2$domain,levels = domain)
    #par(mar=c(8,8,8,8))
    p1 = ggplot(df2, aes(x = group1, y = group2)) + 
      geom_point(aes(size = pct*mulp, fill = group1), alpha = 0.75, shape = 21,stroke=0.25) + 
      scale_size_continuous(limits = c(0.1, mulp), range = c(0.1,mulp), breaks = c(0.1*mulp,0.25*mulp,0.5*mulp,0.75*mulp),labels = c(0.1,0.25,0.5,0.75)) + 
      labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
      theme(legend.key=element_blank(), 
      plot.margin = margin(t = 2, r = 1, b = 1, l = 1, unit = "pt"),
      axis.text.x = element_text(colour = "black", size = 8, face = "bold", angle = 30, vjust = 1, hjust = 1), 
      axis.text.y = element_text(colour = "black", face = "bold", size = 8), 
      legend.text = element_text(size = 8, face ="bold", colour ="black"), 
      legend.title = element_text(size = 8, face = "bold"), 
      panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
      legend.position = "right",
      aspect.ratio = 0.65
           ) +  
      scale_fill_manual(values = col1, guide = FALSE) #+ 
      #scale_y_discrete(limits = rev(levels(df$group2))) 
    options(repr.plot.width=w,repr.plot.height=h)
    plot(p1)
    if(saveplot){
        ggsave(file.path(outdir,paste0(group1,"_",group2,"_toplabel_bubbleplot.png")),width = 6.5,height = 3.5,dpi = 300)
    }
    
    dm_num<-df2 %>% group_by(group2) %>% summarise(num=n())
    #print(dm_num)
    dm_num$group2<-factor(dm_num$group2,levels = levels(df$group2))
    dm_num<-dm_num[order(dm_num$group2),]
    #print(dm_num)
    ct_num<-df2 %>% group_by(group1) %>% summarise(num=n())
    ct_num$group1<-factor(ct_num$group1,levels = levels(df$group1))
    ct_num<-ct_num[order(ct_num$group1),]
    
    # Increase bottom margin
    par(mar=c(8,4,4,4))
    options(repr.plot.width=7.5,repr.plot.height=3.5)
    ymax<-(max(dm_num$num)+1)
    print(ymax)
    p2<-ggplot(dm_num, aes(x = group2, y = num, fill = group2)) +
        geom_bar(stat = "identity", color = NA) +
        #geom_text(aes(label = num), vjust = -1) +
        scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +  # 
        scale_fill_manual(values = col2) +
        #ylim(0, ymax) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = 'none',
              panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
              axis.text.y = element_text(size = 10),      # 
          axis.ticks.y = element_line(color = "black"), # 
          axis.line.y = element_line(color = "black")
             )
    plot(p2)
    if(saveplot){
        ggsave(file.path(outdir,paste0("clonecell_",group2,"_cellbarplot.png")),width = 7.5,height = 3.5,dpi = 300)
    }
    
    par(mar=c(8,4,4,4))
    options(repr.plot.width=8.5,repr.plot.height=3.5)
    ymax<-(max(ct_num$num)+1)
    p3<-ggplot(ct_num, aes(x = group1, y = num, fill = group1)) +
          geom_bar(stat = "identity", color = NA) +
          #geom_text(aes(label = num), vjust = -1) +
          scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +  # 
          scale_fill_manual(values = col1) +
          #ylim(0, ymax) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = 'none',
               panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
                axis.text.y = element_text(size = 10),      # 
          axis.ticks.y = element_line(color = "black"), # 
          axis.line.y = element_line(color = "black")
               )
    plot(p3)
    if(saveplot){
        ggsave(file.path(outdir,paste0("clonecell_",group1,"_numberbarplot.png")),width = 8.5,height = 3.5,dpi = 300)
    }
}   

#202603 defined group colors
col_ident<-c("#ada397","#EAD9C1", "#A3B8C8")
names(col_ident)<-c('all','S3000','S3000_50cs3')

col_r1<-c("#7897AB", 
  "#D885A3", 
  "#98BA7D", 
  "#E2C275",
  "#827397",
  "#FF9F9F",
  "#87AAAA", 
 "#A3C7D6")
names(col_r1)<-c('Other','Ctx','Hip','EPI','TH','CNU','HY','VZ')

col_r2<-c("#D885A3","#7897AB", "#98BA7D", "#FF9F9F", "#827397", "#EBD8C3", "#A2B38B", "#E6BA95",
  "#87AAAA", "#D6AD60", "#B2A4D4", "#F29393", "#7FB5FF", "#B4CFB0", "#E5D3B3", "#94B49F",
  "#FFB2A6", "#867EBA", "#FDCEB9", "#A7D2CB", "#F1D4D4", "#7B8FA1", "#E2C275", "#85A389",
  "#D9A5B3")
names(col_r2)<-c('Other','Ctx_L23','Ctx_L1','Ctx_L4','Ctx_L5','Ctx_L6','CA1','CA3','WM?','Hb','DG','MD','LD','AM','AV','RE','STR','RT','VA/VL/VM','HY','ZI','PAL','CTXsp','OLF','VZ')

col_ct3<-c(
  "#7897AB", "#D885A3", "#98BA7D", "#FF9F9F", "#827397", "#EBD8C3", 
  "#A2B38B", "#E6BA95", "#87AAAA", "#D6AD60", "#B2A4D4", "#F29393", 
  "#7FB5FF", "#B4CFB0", "#E5D3B3", "#94B49F", "#FFB2A6", "#867EBA", 
  "#FDCEB9", "#A7D2CB", "#F1D4D4", "#7B8FA1","#4E6C50"
)
names(col_ct3)<-c('PN_CA1','Olig','Unknown','Astro','GLUT','SPN','Other','DG','OPC','IN_Six3','IN_Lhx6','CPN_L5','CPN_L23',
'CThPN','CPN_L6','PN_L4','IN_Vip','PN_L6b','SCPN','PN_CA3','Myeloid','Mixed'#,'CPN_L23?'
                 )
col_ct2<-c(
  "#7897AB", "#D885A3", "#98BA7D", "#EBD8C3", "#827397", "#FF9F9F", 
  "#A2B38B", "#E6BA95", "#87AAAA", "#D6AD60", "#B2A4D4", "#F29393", 
  "#7FB5FF", "#B4CFB0", "#E5D3B3", "#94B49F", "#FFB2A6", "#867EBA", "#A7D2CB","#4E6C50"
)
names(col_ct2)<-c('CA','Olig','Unknown','Astro','GLUT','SPN','Other','DG','OPC','IN','CPN',
'CThPN','PN_L4','PN_L6b','SCPN','Myeloid','Mixed'#,'CPN?'
                 )

get_specific_group_color<-function(x,lowsat=TRUE){ #x:['major_celltype','minor_celltype','div_col']
    if(x=='major_celltype'){
        return(major_celltype_col)
        }
    if(x=='minor_celltype'){
        return(minor_celltype_col_map)
    }
    if(x=='div_col'){
	if(lowsat){
	col2<-alpha(desaturate(col, 0.4),alpha = 0.7)
	return(col2)
	}else{
	return(col)
	}
    }
    if(x=='ident'){
        return(col_ident)
    }
    if(x=='L2_region'){
        return(col_r2)
    }
    if(x=='L1_region'){
        return(col_r1)
    }
    if(x=='ct2'){
        return(col_ct2)
    }
    if(x=='ct3'){
        return(col_ct3)
    }
}

get_clone_rna_object<-function(object,object_rna){
    object_rna@assays$CloneBC<-object@assays$CloneBC
    object_rna$nCount_CloneBC<-object$nCount_CloneBC
    object_rna$nFeature_CloneBC<-object$nFeature_CloneBC
    object_rna$nCount_CloneBC[is.na(object_rna$nCount_CloneBC)]<-0
    object_rna$nFeature_CloneBC[is.na(object_rna$nFeature_CloneBC)]<-0
    return(object_rna)
}

clone_metric_plot<-function(object,clones,saveplot=TRUE,outdir){
    object$clone.id<-"none"
    meta<-object@meta.data#[colnames(object@assays$CloneBC@counts),]
    meta[clones$cell.barcode,'clone.id']<-clones$clone.id
    meta_sub<-meta[colnames(object@assays$CloneBC@counts),]
    ic_df<-as.data.frame(table(meta_sub$clone.id!='none'))
    colnames(ic_df)<-c("Identify_Clone","Cells")
    #options(repr.plot.width=2.5,repr.plot.height=4)
    p1<-ggplot(ic_df,aes(Identify_Clone,Cells))+
    geom_col()+
    theme_minimal()+
    ggtitle("Cells with CloneID")+
    theme(axis.title = element_text(size = 15),plot.title = element_text(size = 13,face = "bold",hjust = 0.5)) 
    cl_size<-as.data.frame(table(meta_sub[meta_sub$clone.id!='none','clone.id']))
    colnames(cl_size)<-c("clone.id","clonesize")
    cl_size$clone.id<-as.character(cl_size$clone.id)
    cl_size<-cl_size[cl_size$clonesize>1,]
    #meta[!meta$clone.id %in% cl_size$clone.id,'clone.id']<-'none'
    #hist(cl_size$cells,main="Clone size",xlab="Cells")
    p2<-gghistogram(cl_size$clonesize,fill = "grey",title="Clone size",binwidth=1,xlab="cells per clone")+
    theme_minimal()+
    theme(axis.title = element_text(size = 15),plot.title = element_text(size = 14,face = "bold",hjust = 0.5)) 
    options(repr.plot.width=5,repr.plot.height=4)
    par(mfrow=c(1,2),mar = c(3,3,3,3))
    plot(p1+p2)
    if(saveplot){
        ggsave(file.path(outdir,"cloneID_cell_metric.png"),width = 5,height = 4,dpi = 300)
    }
    #reorder the clone id by clone size and assign new id
    cl_size$clone.id2<-as.character(rank(-cl_size$clonesize, ties.method = "first"))
    #cl_map<-c('none',cl_size$clone.id2)
    #names(cl_map)<-c('none',cl_size$clone.id)
    #cl_map
    #meta$clone.id2<-cl_map[match(meta$clone,names(cl_map))]
    meta[,c("clone.id2","clonesize")]<-cl_size[match(meta$clone.id,cl_size$clone.id),c("clone.id2","clonesize")]
    meta$clone.id2[is.na(meta$clone.id2)]<-"none"
    meta$clonesize[is.na(meta$clonesize)]<-0
    #meta$nCount_CloneBC[is.na(meta$nCount_CloneBC)]<-0
    #meta$nFeature_CloneBC[is.na(meta$nFeature_CloneBC)]<-0
    object@meta.data<-meta
    return(object)
}

get_upset_plot<-function(df,group='first_type',only_confident_class=TRUE,min_group_set_size=1){
    if(group=='domain'){
        df[is.na(df$domain),'domain']<-'Undefined'
    }
    if(only_confident_class){
        df_cl<-df %>% filter(clone!='none' & {{ group }} !='Undefined' & spot_class%in%c('singlet','doublet_uncertain')
                    )
        
    }else{
        df_cl<-df %>% filter(clone!='none' & {{ group }} !='Undefined')
    }
    df_cl$clone<-factor(df_cl$clone,levels = unique(df_cl$clone))

    inputlist<-split(df_cl$clone,df_cl[,group])
    nset<-table(lapply(inputlist,length)>min_group_set_size)[['TRUE']]
    upset(fromList(inputlist), 
          text.scale = c(1.3, 1.3, 1, 1, 2, 1.3), #scaling of all axis titles, tick labels, and numbers above the intersection size bars.
          order.by = "freq",
          nsets = nset,
          cutoff = 2,#set a cutoff for the number of intersections per group of sets
          main.bar.color = "#8dabef",
          sets.bar.color = "#78c999"#,empty.intersections = "on"
         )
}

align_spots <- function(meta, img_w=NULL,img_h=NULL,scale_x = 0.8, scale_y = 0.9, 
                        off_x = 0, off_y = 0,
                        rotate=0,
                        flip_y = FALSE#,module='Cluster38'
                       ) {
  df<- data.frame(
    row.names = rownames(meta),
    x_raw = meta$x.axis,
    y_raw = meta$y.axis#,
    #score = meta[,module]
  ) %>%
    mutate(
      x_n = (x_raw - min(x_raw)) / (max(x_raw) - min(x_raw)),
      y_n = (y_raw - min(y_raw)) / (max(y_raw) - min(y_raw)),
    )
  if(flip_y) { df$y_n <- 1 - df$y_n } 
  if(rotate==0){
      df <- df %>%
            mutate(
              x_px = (x_n * img_w * scale_x) + off_x,
              y_px = (y_n * img_h * scale_y) + off_y
            )
  }
  if(rotate==90){
      df <- df %>%
            mutate(
              x_px = ((1-y_n) * img_w * scale_x) + off_x,
              y_px = (x_n * img_h * scale_y) + off_y
            )
  }
  if(rotate==270){
      df <- df %>%
            mutate(
              x_px = (y_n * img_w * scale_x) + off_x,
              y_px = ((1 - x_n) * img_h * scale_y) + off_y
            )
  }
  return(df)
}

get_precise_axis <- function(df,x=8) {
    mat <- as.matrix(df[, c("x_px", "y_px")])
    fit <- principal_curve(mat, 
                           smoother = "smooth_spline", 
                           df = x, 
                           maxit = 150,   
                           stretch = 2)   

    res <- as.data.frame(fit$s[fit$ord, ])
    colnames(res) <- c("x_px", "y_px")
    return(res)
}

calculate_angle <- function(v1, v2) {
    dot_prod <- sum(v1 * v2)
    mag1 <- sqrt(sum(v1^2))
    mag2 <- sqrt(sum(v2^2))
    cos_theta <- abs(dot_prod / (mag1 * mag2)) 
    theta_rad <- acos(pmin(pmax(cos_theta, -1), 1))
    return(theta_rad * 180 / pi)
}

compute_pair_angles <- function(curr_df, axis_ref) {
    idx_pairs <- t(combn(1:nrow(curr_df), 2))

    angles <- apply(idx_pairs, 1, function(p) {
        p1 <- curr_df[p[1], ]
        p2 <- curr_df[p[2], ]
        v_pair <- c(p2$x_px - p1$x_px, p2$y_px - p1$y_px)
        mid_p <- c((p1$x_px + p2$x_px)/2, (p1$y_px + p2$y_px)/2)
        dists <- sqrt((axis_ref$x_px - mid_p[1])^2 + (axis_ref$y_px - mid_p[2])^2)
        near_idx <- which.min(dists)
        ref_idx <- ifelse(near_idx < nrow(axis_ref), near_idx + 1, near_idx - 1)
        v_axis <- c(axis_ref$x_px[ref_idx] - axis_ref$x_px[near_idx],axis_ref$y_px[ref_idx] - axis_ref$y_px[near_idx])
        res_angle <- calculate_angle(v_pair, v_axis)
    })
    return(angles)
}

get_cortex_boundary <- function(points_df) {
  hull_indices <- chull(points_df$x_px, points_df$y_px)
  hull_df <- points_df[hull_indices, ]
  return(hull_df)
}

get_cortex_boundary_flex <- function(points_df, sections = 20) {
  boundary_points <- points_df %>%
    mutate(x_bin = cut(x_px, breaks = sections)) %>%
    group_by(x_bin) %>%
    summarise(
      top = list(filter(cur_data(), y_px == max(y_px))),
      bottom = list(filter(cur_data(), y_px == min(y_px)))
    ) %>%
    tidyr::unnest(cols = c(top, bottom), names_repair = "unique") %>%
    select(x_px = x_px...2, y_px = y_px...3) 
  upper <- boundary_points %>% arrange(x_px)
  lower <- boundary_points %>% arrange(desc(x_px))
  
  hull_df <- bind_rows(upper, lower)
  return(hull_df)
}