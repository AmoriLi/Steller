library(Seurat)
library(png)
library(dplyr)
library(tibble)
library(ggplot2)
library(Matrix)
library(broom)
#library(grid)
#library(gridExtra)
#library(magick)
library(purrr)
library(cowplot)
#library(figpatch)
#library(pdftools)
library(cluster)
library(ggpubr)
library(plotly)
library(patchwork)
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")

args <- commandArgs(trailingOnly = T)
indir = args[1] #"/public/home/zxli_gibh/NPC_project/BMK/tissue/brain/P5" #args[1] #e.g. /public/home/zxli_gibh/NPC_project/BMK/tissue/brain/E15.5/2024.5.7
#slices<-c("S3000-A4","S3000-B2")#list.files(indir,'^S')
script.dir<-"~/NPC_project/BMK/pipeline/RNA"
outdir<-file.path(indir,"Merged_2504/1_cluster")
if(!dir.exists(outdir)) {dir.create(outdir,recursive=TRUE)}
exp_col<-c('#034f84','#92a8d1','#d6d4e0','#f4a688',"#ED797B","#d64161","#c94c4c")

slices<-c("S3000","S3000_50cs3")
object_list<-list()
for(i in 1:length(slices)){
    mtx = file.path(indir,slices[i],"RNA/BST/07.CellSplit/mtx")
    image = file.path(indir,slices[i],"RNA/BST/05.AllheStat/allhe/he_roi_small.png") 
    coord = file.path(indir,slices[i],"RNA/BST/07.CellSplit/mtx/barcodes_pos.tsv.gz") 
    fluo = file.path(indir,slices[i],"RNA/BST/07.CellSplit/images/fluorescence.tif") 
    tmpdir = file.path(indir,slices[i],"RNA/BST/07.CellSplit/spatial") 
    if(!dir.exists(tmpdir)) {dir.create(tmpdir)}
    file.copy(from = paste0(script.dir,"/scalefactors_json.json"),to = tmpdir)
    file.copy(from = image,to = tmpdir)
    he_png <- readPNG(source = image)
    #tissue.positions <- read.table(gzfile(coord,"rt"))
    tissue.positions <- read.csv(gzfile(coord),sep='\t',header=F)

    tissue.positions$tissue <- 1
    names(tissue.positions) <- c("barcodes","row","col","tissue")
    tissue.positions <- tissue.positions[,c("barcodes","tissue","row","col")]

    width <- ncol(he_png)
    height <- nrow(he_png)

    tiff <- paste0("less ",fluo,"| cat")
    zoom = system(tiff,intern = TRUE)
    zoom_scale_width <- width/as.numeric(unlist(strsplit(unlist(strsplit(zoom,split = " "))[3],split = "x"))[1])
    zoom_scale_height <- height/as.numeric(unlist(strsplit(unlist(strsplit(zoom,split = " "))[3],split = "x"))[2])
    tissue.positions$imagerow <- tissue.positions$col*zoom_scale_height
    tissue.positions$imagecol <- tissue.positions$row*zoom_scale_width
    barcode = read.table(paste0(mtx,"/barcodes.tsv.gz"),header = F)
    tissue.positions = tissue.positions[tissue.positions$barcodes%in%barcode$V1,]

    write.table(x = tissue.positions,file = file.path(tmpdir,"tissue_positions_list.csv"),
                quote = F,sep = ",",
                append = F,
                col.names = F,
                row.names = F)

    cmd = paste0("cp ",image," ",paste0(tmpdir,"/tissue_lowres_image.png"))
    system(cmd)
    ##read expression matrix
    data <- Read10X(data.dir = mtx)
    ##create seurat object
    object <- CreateSeuratObject(counts = data,min.cells = 0,min.features = 0,assay = "Spatial")
    ##load image
    image <- Read10X_Image(image.dir = tmpdir,image.name = "tissue_lowres_image.png")
    ##assign spatial assay for image
    DefaultAssay(image) <-"Spatial"
    object[[slices[i]]] <- image

    object_list[[slices[i]]]<-object
}

object<-object_list[[1]]
object$orig.ident<-slices[1]
object$barcode<-paste(object$orig.ident,rownames(object@meta.data),sep = "_")
slices_rn<-gsub("-",".",slices)
object$x.axis<-object@images[[slices_rn[1]]]@coordinates$row
object$y.axis<-object@images[[slices_rn[1]]]@coordinates$col
object$imagerow<-object@images[[slices_rn[1]]]@coordinates$imagerow
object$imagecol<-object@images[[slices_rn[1]]]@coordinates$imagecol
object<-RenameCells(object,new.names = object$barcode)
#head(da@meta.data)
for(i in 2:length(slices)){
    temp<-object_list[[i]]
    temp$orig.ident<-slices[i]
    temp$barcode<-paste(temp$orig.ident,rownames(temp@meta.data),sep = "_")
    temp$x.axis<-temp@images[[slices_rn[i]]]@coordinates$row
    temp$y.axis<-temp@images[[slices_rn[i]]]@coordinates$col
    temp$imagerow<-temp@images[[slices_rn[i]]]@coordinates$imagerow
    temp$imagecol<-temp@images[[slices_rn[i]]]@coordinates$imagecol
    temp<-RenameCells(temp,new.names = temp$barcode)
    #head(da@meta.data)
    object<-merge(object,temp)
}

#modify the coord when multiple samples prevent for overlapping
object$x.coord<-object$x.axis
object$y.coord<-object$y.axis
xmax<-max(object$x.axis[object$orig.ident=='S3000'])
ymax<-max(object$y.axis[object$orig.ident=='S3000'])
object$x.coord[object$orig.ident=="S3000_50cs3"]<-object$x.coord[object$orig.ident=="S3000_50cs3"] + xmax + 50
object$y.coord[object$orig.ident=="S3000_50cs3"]<-object$y.coord[object$orig.ident=="S3000_50cs3"] + ymax + 50
saveRDS(object = object,file = paste0(outdir,"/object.merged.spatial.raw.rds"))

da_s<-da[,da$nCount_Spatial<4000 & da$nCount_Spatial>300 & 
         da$nFeature_Spatial<2000 & da$nFeature_Spatial>150 & 
         da$spots>8 & da$spots<80]
da_s

#modify the coord when multiple samples prevent for overlapping
object<-da_s
object$x.coord<-object$x.axis
object$y.coord<-object$y.axis
xmax<-max(object$x.axis[object$orig.ident=='S3000'])
ymax<-max(object$y.axis[object$orig.ident=='S3000'])

#
right_rotation_matrix <- matrix(c(0, 1, -1, 0), ncol = 2)
right_rotated_points <- as.matrix(object@meta.data[object$orig.ident=="S3000",c("x.coord","y.coord")]) %*% right_rotation_matrix
left_rotation_matrix <- matrix(c(0, -1, 1, 0), ncol = 2)
left_rotated_points <- as.matrix(object@meta.data[object$orig.ident=="S3000_50cs3",c("x.coord","y.coord")]) %*% left_rotation_matrix
cell_meta<-rbind(as.data.frame(right_rotated_points),as.data.frame(left_rotated_points))
head(cell_meta)
meta<-object@meta.data
meta$x.coord<-cell_meta[rownames(meta),"V1"]
meta$y.coord<-cell_meta[rownames(meta),"V2"]
meta$y.coord<-ifelse(meta$orig.ident=="S3000",meta$y.coord+10000,meta$y.coord-10000)
meta$x.coord<-ifelse(meta$orig.ident=="S3000",meta$x.coord+1200,meta$x.coord)
options(repr.plot.width=20,repr.plot.height=9)
get_spatial_dim_plot(meta,group = "orig.ident",x="x.coord",y="y.coord")
object@meta.data<-meta[colnames(object),]

da_s<-object
writeMM(da_s@assays$Spatial@counts,file.path(outdir,"counts.mtx"))
#system("gzip counts.mtx")
write.csv(da_s@meta.data,file.path(outdir,"meta.csv"))
write.csv(rownames(da_s),file.path(outdir,"gene_names.csv"))
saveRDS(da_s,file.path(outdir,"object.merged.spatial.adjusted.rds"))

### integration strategy
da<-readRDS("object.merged.spatial.adjusted.rds")
da1<-da[,da$orig.ident=='S3000']
da2<-da[,da$orig.ident=='S3000_50cs3']
ifnb.list <- list(da1,da2)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list,assay=c('Spatial','Spatial'))
anchors <- FindIntegrationAnchors(object.list = ifnb.list,assay=c('Spatial','Spatial'),
                                  anchor.features = features, reduction = "cca")

# this command creates an 'integrated' data assay
da <- IntegrateData(anchorset = anchors)
# Run the standard workflow for visualization and clustering
da <- ScaleData(da, verbose = FALSE)
da <- RunPCA(da, npcs = 30, verbose = FALSE#,features = tf$`Gene Symbol`
            )
da <- RunUMAP(da, reduction = "pca", dims = 1:30)
da <- FindNeighbors(da, reduction = "pca", dims = 1:30)
da <- FindClusters(da, resolution = 0.2)
saveRDS(da,"object.integrated.spatial.clustered.rds")

da<-readRDS(file.path(outdir,"object.integrated.spatial.clustered.rds"))
col1<-c("#EAD9C1", "#A3B8C8"
        )
names(col1)<-unique(da$orig.ident)
#
plot_df <- data.frame(
    nCount = da@meta.data$nCount_Spatial,
    nFeature = da@meta.data$nFeature_Spatial,
    group = da@meta.data$orig.ident
)

plot_df %>% group_by(group) %>%
  summarise(
    col1_mean = mean(nCount, na.rm = TRUE),
    col1_median = median(nCount, na.rm = TRUE),
    
    col2_mean = mean(nFeature, na.rm = TRUE),
    col2_median = median(nFeature, na.rm = TRUE),
    
    n = n()
)

options(repr.plot.width=2.5,repr.plot.height=3.5)
p <- ggplot(plot_df, aes(x = group, y = nCount, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", color = "black") +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) + 
    geom_hline(yintercept = c(300, 4000), linetype = 2, color = "darkgrey", linewidth = 0.5) +
    scale_fill_manual(values = col1)+
    theme_classic() +
    labs(y = "nCount_Spatial", x = "Sample") +
    NoLegend() 

print(p)
ggsave(file.path(outdir,"nCount_Spatial.byIdents.Vlnplot.png"), p, width = 2.5, height = 3.5)

options(repr.plot.width=2.5,repr.plot.height=3.5)
p <- ggplot(plot_df, aes(x = group, y = nFeature, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", color = "black") +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = c(300, 2000), linetype = 2, color = "darkgrey", linewidth = 0.5) +
    scale_fill_manual(values = col1)+
    theme_classic() +
    labs(y = "nFeature_Spatial", x = "Sample") +
    NoLegend() 

print(p)
ggsave(file.path(outdir,"nFeature_Spatial.byIdents.Vlnplot.png"), p, width = 2.5, height = 3.5)

col_cluster<-get_specific_group_color('div_col')[1:n_distinct(da$seurat_clusters)]
names(col_cluster)<-unique(da$seurat_clusters)
col_cluster

options(repr.plot.width=6,repr.plot.height=5)
DimPlot(da,cols = col_cluster,pt.size = 0.75,raster = TRUE)+
    coord_fixed()
ggsave(file.path(outdir,'cluster.dimplot.png'),width = 6,height = 5,dpi = 400)

options(repr.plot.width=6,repr.plot.height=5)
DimPlot(da,cols = col1,group.by = 'orig.ident',pt.size = 1.5)+coord_fixed()
ggsave(file.path(outdir,"UMAP.byIdent.png"),width=6,height=5,dpi=300)

options(repr.plot.width=20,repr.plot.height=10)
get_spatial_dim_plot(da@meta.data,group = "seurat_clusters",x = "x.coord",y="y.coord",pt.size = 0.5,stroke = 0.01) +
    scale_fill_manual(values = col_cluster) +
    coord_fixed()+
    guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(file.path(outdir,"Spatial.Cluster.DimPlot.png"),width=20,height=10,dpi=300)

