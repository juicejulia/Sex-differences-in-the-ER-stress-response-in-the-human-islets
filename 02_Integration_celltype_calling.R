# Integration, v5 workflow
# Integration
## Switch to Seurat v5 workflow for streamlined integration
library(Seurat)
library(SeuratWrappers)
library(SeuratData)
soupx_merged_data <- readRDS("soupx_merged_data.rds")
ERstress.integrated = soupx_merged_data
ERstress.integrated<- JoinLayers(ERstress.integrated,  assay = "RNA")
ERstress.integrated[["RNA"]] <- split(ERstress.integrated[["RNA"]], f = ERstress.integrated$sample)

# analysis without integration
#ERstress.integrated <- SCTransform(ERstress.integrated)
#ERstress.integrated <- RunPCA(ERstress.integrated)
#ERstress.integrated <- RunUMAP(ERstress.integrated, dims = 1:30)
#DimPlot(ERstress.integrated, reduction = "umap", group.by = c("condition", "donor"))
# canonical workflow results in better looking UMAP. 
## also, IntegrateLayers does not work for SCTransformed data in some of the methods

ERstress.integrated <- NormalizeData(ERstress.integrated)
ERstress.integrated <- FindVariableFeatures(ERstress.integrated)
ERstress.integrated <- ScaleData(ERstress.integrated)
ERstress.integrated <- RunPCA(ERstress.integrated)
ERstress.integrated <- FindNeighbors(ERstress.integrated, dims = 1:30, reduction = "pca")
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "unintegrated_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ERstress.integrated, reduction = "umap.unintegrated", group.by = c("condition", "donor","unintegrated_clusters"))
list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.unintegrated"))
}
dev.off()


##cca
ERstress.integrated <- IntegrateLayers(
  object = ERstress.integrated, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
ERstress.integrated <- FindNeighbors(ERstress.integrated, reduction = "integrated.cca", dims = 1:30)
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "cca_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(
  ERstress.integrated,
  reduction = "umap.cca",
  group.by = c("condition", "donor", "cca_clusters"),
  label.size = 2
)
list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.cca"))
}
dev.off()

##rpca
ERstress.integrated <- IntegrateLayers(
  object = ERstress.integrated, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)
ERstress.integrated <- FindNeighbors(ERstress.integrated, reduction = "integrated.rpca", dims = 1:30)
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "rpca_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(
  ERstress.integrated,
  reduction = "umap.rpca",
  group.by = c("condition", "donor", "rpca_clusters"),
  label.size = 2
)
list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.rpca"))
}
dev.off()


###harmony
ERstress.integrated <- IntegrateLayers(
  object = ERstress.integrated, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
ERstress.integrated <- FindNeighbors(ERstress.integrated, reduction = "harmony", dims = 1:30)
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "harmony_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(
  ERstress.integrated,
  reduction = "umap.harmony",
  group.by = c("condition", "donor", "harmony_clusters"),
  label.size = 2
)

list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.harmony"))
}
dev.off()


###mnn
ERstress.integrated <- IntegrateLayers(
  object = ERstress.integrated, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
ERstress.integrated <- FindNeighbors(ERstress.integrated, reduction = "integrated.mnn", dims = 1:30)
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "mnn_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
DimPlot(
  ERstress.integrated,
  reduction = "umap.mnn",
  group.by = c("condition", "donor", "mnn_clusters"),
  label.size = 2
)

list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.mnn"))
}
dev.off()

## sctransform
ERstress.integrated <- SCTransform(ERstress.integrated, vars.to.regress = "percent.mt", verbose = FALSE,method = "glmGamPoi")
ERstress.integrated <- RunPCA(ERstress.integrated,reduction.name="pca_sctransform")
ERstress.integrated <- FindNeighbors(ERstress.integrated, dims = 1:30, reduction = "pca_sctransform")
ERstress.integrated <- FindClusters(ERstress.integrated, resolution = 2, cluster.name = "unintegrated_sctransform_clusters")
ERstress.integrated <- RunUMAP(ERstress.integrated, dims = 1:30, reduction = "pca_sctransform", reduction.name = "umap.unintegrated.sctransform")
DimPlot(ERstress.integrated, reduction = "umap.unintegrated", group.by = c("condition", "donor","unintegrated_sctransform_clusters"))
list <- c("GCG", "INS", "SST", "GHRL","PPY","CFTR","CPA2","SPARC", "VWF","PTPRC","BRCA1")
pdf("test.pdf")
for (i in list)
{
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.unintegrated"))
}
dev.off()

## Choose rpca because it is more conservative. Seurat recommends RPCA where: 
# (a) A substantial fraction of cells in one dateset have no matching type in the other
# (b) Datasets originate from the same platform (i.e. multiple lands of 10x genomics)
# (c) There are a large number of datasets or cells to integrate
# (d) Integration results look good. 
list <- c("PRSS1","INS","GCG","PNLIP","CPA2","CPB1","CTRC","CFTR","PKHD1",
          "SST","DCLK1","ZEB2", "ENG","BRCA1","PTPRC","GHRL","PPY","COL1A1", 
          "COL1A2", "COL6A3", "COL3A1", "TIMP3",  "VSTM2L", "DEFB1","LUM",
          "VWF","LAPTM5", "COL22A1","MEIS2","PRSS23","PLVAP","RGCC","PECAM1",
          "ESM1","SERPINE1","CLDN5","STC1","MMP1","GNG11","IAPP","DLK1","G6PC2",
          "SPP1","KRT7","KRT18","KRT19","SOX9","NKX2-2","CHGA","SPARC")


DefaultAssay(ERstress.integrated)="RNA"
for (i in list)
{
  pdf(paste0(i,".ERstress.integrated.umap.pdf"))
  print (FeaturePlot(ERstress.integrated, features = i, cols = c("lightgrey", "red"), reduction="umap.rpca"))
  dev.off()
}

for (i in list)
{
  pdf(paste0(i,".ERstress.integrated.violin.pdf"),height=5,width=10)
  print (VlnPlot(ERstress.integrated, features = i, group.by="rpca_clusters")+NoLegend())
  dev.off()
}


new.cluster.ids <- c("fibroblast","ductal","alpha","alpha","ductal","alpha","ductal","doublets",
                     "acinar","ductal","doublets","fibroblast","acinar","alpha","ductal","ductal",
                     "endothelial","ductal","endothelial","beta","fibroblast","alpha","beta","acinar","fibroblast",
                     "alpha","doublets","alpha","ductal","doublets","fibroblast","beta","doublets",
                     "alpha","PP","endothelial","doublets","alpha","delta","immune","immune",
                     "acinar","ductal","beta","alpha","acinar","endothelial","doublets","doublets",
                     "doublets","endothelial","immune","doublets","doublets","beta","immune","delta",
                     "doublets","fibroblast","beta","immune"
)

ERstress.integrated$rpca_clusters <- factor(ERstress.integrated$rpca_clusters, 
                                            levels = c(0:60))

names(new.cluster.ids) <- levels(ERstress.integrated$rpca_clusters)
Idents(ERstress.integrated)=ERstress.integrated$rpca_clusters
ERstress.integrated <- RenameIdents(ERstress.integrated, new.cluster.ids)
ERstress.integrated$cell.type <- Idents(ERstress.integrated)

DimPlot(ERstress.integrated, reduction = "umap.rpca", label = TRUE, pt.size = 0.5,raster=FALSE)  #UMAP_1
ERstress.integrated <- JoinLayers(ERstress.integrated)

## Azimuth reference mapping
## Used the annotation of epsilon cells from azimuth
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
InstallData("pancreasref.SeuratData")
ERstress.integrated <- RunAzimuth(ERstress.integrated,reference="pancreasref")
DimPlot(ERstress.integrated, reduction = "umap.rpca", group.by="predicted.annotation.l1",label = TRUE, pt.size = 0.5,raster=FALSE) 
Idents(ERstress.integrated) = ERstress.integrated$"predicted.annotation.l1"
epsilon.cells.azimuth <- colnames(subset(ERstress.integrated,idents="epsilon"))

DimPlot(ERstress.integrated,reduction="umap.rpca",cells.highlight=epsilon.cells.azimuth,cols.highlight=c("firebrick2"),cols="grey") #UMAP_epsilon
ERstress.integrated$cell.type.azimuth <- ifelse(ERstress.integrated$predicted.annotation.l1%in%c("activated_stellate","quiescent_stellate"),"fibroblast",
                                                ifelse(ERstress.integrated$predicted.annotation.l1=="gamma","PP",
                                                       ERstress.integrated$predicted.annotation.l1))

# concordance between the manual supervised annotation and the azimuth annotation, excluding the doublets population
sum(ERstress.integrated$cell.type.azimuth==ERstress.integrated$cell.type)
# 119223
sum(ERstress.integrated$cell.type!="doublets")
# 134330
# 119223/134330=0.8875382

ERstress.integrated$cell.type <- as.character(ERstress.integrated$cell.type)
ERstress.integrated$cell.type.final <- ifelse(ERstress.integrated$predicted.annotation.l1=="epsilon","epsilon",ERstress.integrated$cell.type)                                                       
# Carry over the epsilon cell marker 

Idents(ERstress.integrated)=ERstress.integrated$cell.type.final
Idents(ERstress.integrated)=factor(Idents(ERstress.integrated),levels=c("alpha","beta","delta","epsilon","PP",
                                                                        "ductal","acinar","fibroblast","endothelial",
                                                                        "immune","doublets"))
DimPlot(ERstress.integrated, reduction = "umap.rpca",  label = TRUE, pt.size = 0.5,raster=FALSE) 
DimPlot(ERstress.integrated, reduction = "umap.rpca", cols=c("#88BDFF","#FF8000","#FFFA0A","#FF6347","#10C3BE",
                                                             "#704700","#7CAE00","#AEA200","#C77CFF",
                                                             "#F96DDD","#90EE90"),
        raster=FALSE,label = FALSE)  #UMAP3
DimPlot(ERstress.integrated, reduction = "umap.rpca",  group.by="donor",label = FALSE, pt.size = 0.5,raster=FALSE,alpha=0.4) #UMAP_4
DimPlot(ERstress.integrated, reduction = "umap.rpca",  group.by="condition",label = FALSE, pt.size = 0.5,raster=FALSE,alpha=0.4) #UMAP_5
DimPlot(ERstress.integrated, reduction = "umap.rpca",  group.by="sex",cols= c("#7FC97F", "#BEAED4"), label = FALSE, pt.size = 0.5,raster=FALSE,alpha=0.2) #UMAP_6

#Dotplot
features=c("GCG", "INS", "SST","GHRL","PPY","CFTR","PRSS1","SPARC", "PLVAP","PTPRC")
Idents(ERstress.integrated)=factor(Idents(ERstress.integrated),levels=c("doublets","immune","endothelial","fibroblast","acinar","ductal",
                                                                        "PP","epsilon","delta","beta","alpha"
))

DotPlot(ERstress.integrated, 
        assay="RNA",
        features =features) +   
  RotatedAxis() +
  theme(axis.text.y = element_text(face="bold",size = 15), axis.text.x = element_text(face="bold",size = 15)) +
  xlab("") + ylab("")

#Stacked violin
library(Seurat)
library(ggplot2)
library(cowplot)

list <- c("PRSS1","INS","GCG","PNLIP","CPA2","CPB1","CTRC","CFTR","PKHD1","SST","DCLK1","ZEB2", "ENG","BRCA1","PTPRC","GHRL","PPY","COL1A1", "COL1A2", "COL6A3", "COL3A1", "TIMP3",  "VSTM2L", "DEFB1","LUM","VWF","LAPTM5", "COL22A1","MEIS2","PRSS23","PLVAP","RGCC","PECAM1","ESM1","SERPINE1","CLDN5","STC1","MMP1","GNG11","IAPP","DLK1","G6PC2","SPP1","KRT7","KRT18","KRT19","SOX9","NKX2-2","CHGA","SPARC")
pdf("test.pdf",height=5,width=10)
for (i in list)
{
  print (VlnPlot(ERstress.integrated, features = i,pt.size=0)+
           geom_boxplot(width=0.3,fill="white", outlier.size=0.1) +
           theme(legend.position = 'none'))
}
dev.off()

features=c("GCG", "INS", "SST","GHRL","PPY","CFTR","PRSS1","SPARC", "PLVAP","PTPRC")
VlnPlot(ERstress.integrated, features, stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none") + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
               geom="pointrange", color="black",
               shape = 18, size = 0.75,
               position = position_dodge(width = 0.9)) # adding mean and sd

pdf("test.pdf",height=18,width=10)
VlnPlot(ERstress.integrated, features, stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none") +
  stat_summary(fun=median, 
               geom="point", color="black",
               size = 2,
               position = position_dodge(width = 0.9)) # adding median
dev.off()



saveRDS(ERstress.integrated,"ERstress.integrated.rds")

