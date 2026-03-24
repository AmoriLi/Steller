library(purrr)
library(cowplot)
library(png)
#library(pdftools)
library(cluster)
library(ggpubr)
library(stringr)
library(plotly)
library(patchwork)
source("~/NPC_project/BMK/pipeline/RNA/Helper.R")

wd<-"~/NPC_project/BMK/tissue/brain/P5/Merged_2512/2_region/L9"
if(!dir.exists(wd)){dir.create(wd,recursive=TRUE)}
indir<-"~/NPC_project/BMK/tissue/brain/P5/S3000/RNA/BST/05.AllheStat"
script.dir <- "~/NPC_project/BMK/pipeline/RNA" 
outdir<-wd
mtx = file.path(indir,"heAuto_level_matrix/subdata/L9_heAuto")
image = file.path(indir,"allhe/he_roi_small.png") 
coord = file.path(indir,"heAuto_level_matrix/subdata/L9_heAuto/barcodes_pos.tsv.gz") 
tmpdir = file.path(indir,"heAuto_level_matrix/subdata/L9_Spatial") 
if(!dir.exists(tmpdir)) {dir.create(tmpdir)}
file.copy(from = paste0(script.dir,"/scalefactors_json.json"),to = tmpdir)
file.copy(from = image,to = tmpdir)
he_png <- readPNG(source = image)
    #tissue.positions <- read.table(gzfile(coord,"rt"))
tissue.positions <- read.csv(gzfile(coord),sep='\t',header=F)

tissue.positions$tissue <- 1
names(tissue.positions) <- c("barcodes","row","col","tissue")
tissue.positions <- tissue.positions[,c("barcodes","tissue","row","col")]

tissue.positions$imagerow <- tissue.positions$col #*zoom_scale_height
tissue.positions$imagecol <- tissue.positions$row #*zoom_scale_width
    ##提取mtx下的barcode,保证一致
barcode = read.table(paste0(mtx,"/barcodes.tsv.gz"),header = F)
tissue.positions = tissue.positions[tissue.positions$barcodes%in%barcode$V1,]

write.table(x = tissue.positions,file = file.path(tmpdir,"tissue_positions_list.csv"),
                quote = F,sep = ",",
                append = F,
                col.names = F,
                row.names = F)

cmd = paste0("cp ",image," ",paste0(tmpdir,"/tissue_lowres_image.png"))
system(cmd)
data <- Read10X(data.dir = mtx)
object <- CreateSeuratObject(counts = data,min.cells = 0,min.features = 0,assay = "Spatial")
image <- Read10X_Image(image.dir = tmpdir,image.name = "tissue_lowres_image.png")
DefaultAssay(image) <-"Spatial"
object[["S3000"]] <- image

object$orig.ident<-'S3000'
object$barcode<-paste(object$orig.ident,rownames(object@meta.data),sep = "_")
slices_rn<-'S3000'#gsub("-",".",slices)
object$x.axis<-object@images[[slices_rn[1]]]@coordinates$row
object$y.axis<-object@images[[slices_rn[1]]]@coordinates$col
object$imagerow<-object@images[[slices_rn[1]]]@coordinates$imagerow
object$imagecol<-object@images[[slices_rn[1]]]@coordinates$imagecol
object<-RenameCells(object,new.names = object$barcode)
#head(da@meta.data)

saveRDS(object = object,file = paste0(outdir,"/P5_S3000_object.spatial.raw.rds"))

m.s.genes<-c('Exo1','Msh2','Mcm4','Chaf1b','Rrm2','Cenpu','Mrpl36','Gmnn','Hells','Cdc6','Gins2','Uhrf1',
             'Cdc45','Slbp','Ubr7','Fen1','Rad51ap1','Mcm5','Nasp','Cdca7','Blm','Usp1','Ung','Prim1',
             'Clspn','Mcm6','Dtl','Pola1','Dscc1','Tipin','Wdr76','Casp8ap2','Tyms','Ccne2','Rrm1','Polr1b',
             'Rfc2','Rad51','E2f8','Mcm7','Pcna')

m.g2m.genes<-c('Dlgap5','Ctcf','Smc4','Kif20b','Cdc25c','Gtse1','Tpx2','Hmgb2','Cks1brt','Cdca2',
            'Top2a','Cks2','Cdca3','G2e3','Ttk','Ncapd2','Lbr','Anp32e','Ckap2','Tacc3','Nek2',
            'Cenpe','Kif11','Anln','Hjurp','Aurkb','Rangap1','Cks1b','Hmmr','Ckap5','Cdc20',
             'Psrc1','Kif23','Ect2','Kif2c','Ndc80','Nuf2','Cdca8','Birc5','Cenpf','Ube2c','Jpt1',
             'Pimreg','Nusap1','Mki67','Tubb4b','Bub1','Cenpa','Ccnb2','Aurka','Ckap2l')
###remove Hb/mt/ribo genes that interfere celltype marker identification
object <- object[!grepl("^Hb",rownames(object)),colnames(object)]
object <- object[!grepl("^mt-",rownames(object)),colnames(object)]#PercentageFeatureSet(da,"^mt-",col.name = "percent_mt")
object <- object[!grepl("^Rp[s,l]",rownames(object)),colnames(object)]#PercentageFeatureSet(da,"^Rp[s,l]",col.name = "percent_ribo")

#object <- CellCycleScoring(object,s.features = m.s.genes,g2m.features = m.g2m.genes)

object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
object <- ScaleData(object)

object <- RunPCA(object=object,pc.genes = VariableFeatures(object))
object <- FindNeighbors(object=object, reduction = "pca",dims = 1:30,k.param = 80,verbose = F)
object <- FindClusters(object=object, verbose = F,resolution = 1.6)
object <- RunUMAP(object=object, reduction = "pca",dims = 1:30)

saveRDS(object,paste0(outdir,"/P5_S3000_object.spatial.clustered.rds"))