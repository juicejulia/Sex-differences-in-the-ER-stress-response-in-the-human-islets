library(Seurat)
library(SeuratData)
library(patchwork)
library(clustree)
library(limma)
library(Biobase)
library(edgeR)
library(plyr)
library(dplyr)
library(SoupX)


# SoupX
# SoupX to remove ambient RNA
## Modified from https://github.com/Gaulton-Lab/HPAP-scRNA-seq/blob/main/HPAP-SoupX.R

samples <- c('HP23166_0h_Tg','HP23166_12h_Tg','SAMN35848421_0h_Tg','SAMN35848421_12h_Tg','SAMN35848421_48h_Tg','SAMN36705973_0h_Tg',
             'SAMN36705973_12h_Tg','SAMN36705973_48h_Tg','SAMN36823227_0h_Tg','SAMN36823227_12h_Tg','SAMN36823227_48h_Tg','SAMN37871873_0h_Tg',
             'SAMN37871873_12h_Tg','SAMN37871873_48h_Tg','SAMN39523303_0h_Tg','SAMN39523303_12h_Tg','SAMN39523303_48h_Tg')



contam_frac_results0 <- NULL
contam_frac_results <- data.frame()
for (x in samples){
  
  data= readRDS(paste0(x,'.rds'))
  DefaultAssay(data)="RNA"
  toc <- GetAssayData(object = data, slot = 'counts') 
  tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
  tod_1 <- tod[rownames(toc),]
  
  metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                     as.data.frame(Idents(data))
  ))
  colnames(metadata) <- c('RD1','RD2','Cluster')
  
  sc <- SoupChannel(tod_1,toc)
  sc <- setDR(sc,metadata[colnames(sc$toc),c('RD1','RD2')])
  sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
  sc <- autoEstCont(sc)
  save <- sc$fit$dd[,c('gene', 'soupExp')]
  colnames(save) <- c('gene', x)
  write.table(save, paste0(x, '_gene_level_soup_exp.tsv'), sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  contamination_fraction <- mean(sc$metaData$rho*100)
  contam_frac_results0$Sample <- x
  contam_frac_results0$ContaminationFraction <- contamination_fraction
  contam_frac_results <- rbind(contam_frac_results,contam_frac_results0) #This will output a table with calculated contamination fraction for each sample
  
  out <- adjustCounts(sc, roundToInt=TRUE) #Stochastically adjusts counts while rounding to maintain overall contamination fraction while outputting integer counts
  
  data2 <- CreateSeuratObject(out)
  data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
  data2 <- SCTransform(data2, vst.flavor = "v2", verbose = FALSE) 
  data2 <- RunPCA(data2, npcs=30) 
  data2 <- RunUMAP(data2,reduction="pca",dims=1:30,verbose=FALSE) 
  data2 <- FindNeighbors(data2,reduction = "pca", dims = 1:30, verbose = FALSE)
  data2 <- FindClusters(data2,resolution=1)
  saveRDS(data2, file = sprintf('%s_SoupX.rds',x))
}    

write.table(contam_frac_results,"contam_frac_results.csv",quote=F,sep=",")

#Create a merged Seurat object from the individual sample post-SoupX Seurat objects
soupx_files <- list.files(pattern='_SoupX.rds')
soupx_data <- list()

for (x in soupx_files){
  sample <- sub('_[^_]*$', '', x)
  tmp <- readRDS(x)
  soupx_data[[sample]] <- tmp
}

soupx_merged_data <- merge(soupx_data[[samples[[1]]]], y=soupx_data[samples[2:length(samples)]], add.cell.ids=samples, project='ERstress')
soupx_merged_data$sample <- sub('_[^_]*$', '',rownames(soupx_merged_data@meta.data))
soupx_merged_data$donor <- sub('\\_.*','',soupx_merged_data$sample)
soupx_merged_data$condition <- sub("^.*?\\_", "",soupx_merged_data$sample)
soupx_merged_data$sex <- ifelse(soupx_merged_data$donor%in%c("HP23166","SAMN36705973","SAMN36823227"),"female","male")
soupx_merged_data[['percent.mt']] <- PercentageFeatureSet(soupx_merged_data, pattern = '^MT-')
saveRDS(soupx_merged_data,"soupx_merged_data.rds")
