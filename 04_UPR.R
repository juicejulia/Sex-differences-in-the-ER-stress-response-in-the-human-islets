# UPR pathway
## UPR pathway modules
## pathways genes were curated from GO terms: ATF6 (GO:0036500), IRE1 (GO:0036498), and PERK (GO:0036499) 

library(Seurat)
library(pheatmap)
library(plyr)
library(dplyr)

PERK_pathway <- c("TMED2","NFE2L2","EIF2S1","DDIT3","EIF2AK3","ATF4","RPAP2","QRICH1")
IRE1_pathway <- c("PARP16","DNAJC10","ERN1","XBP1","PTPN1","VAPB", "HSPA5","HERPUD1")
ATF6_pathway <- c("CREBZF","ATF6B","MBTPS1","MBTPS2","ATF6")
genelist <- c(PERK_pathway, IRE1_pathway, ATF6_pathway)


## Add module scores
## select only alpha, beta, and delta cells
features <- list(PERK_pathway=PERK_pathway, IRE1_pathway=IRE1_pathway, ATF6_pathway=ATF6_pathway)

Endo_ERstress.integrated <- subset(ERstress.integrated, idents=c("alpha","beta","delta"))
temp <- AddModuleScore(object= Endo_ERstress.integrated, 
                       features=features,
                       name='UPR')
names(temp@meta.data)[grep("UPR", names(temp@meta.data))] <- names(features)
meta.data <- temp@meta.data
# Group by condition and sex, and calculate the mean for each score
module_data <- meta.data[, c("donor", "condition", "sex","cell.type.final", "PERK_pathway", "IRE1_pathway", "ATF6_pathway")]

module_average <- module_data %>%
  group_by(cell.type.final,condition, sex) %>%
  dplyr::summarise(
    across(c(PERK_pathway,IRE1_pathway,ATF6_pathway),
           mean, na.rm=TRUE),
    .groups = "drop") %>%
  mutate(cell_type_sex_condition = paste0(cell.type.final,"_", sex, "_",condition))

module_average_1 <-  module_average %>% 
  arrange(cell.type.final, sex, condition)

module_average_combined <- as.data.frame(module_average_1[, c("cell_type_sex_condition",
                                                              "ATF6_pathway", 
                                                              "IRE1_pathway", 
                                                              "PERK_pathway")])
rownames(module_average_combined) = module_average_combined[,1]
module_average_combined = module_average_combined[,-1]
module_average_combined_1 <- as.data.frame(t(module_average_combined))
paletteLength <- 100
myColor <- colorRampPalette(c("grey", "white", "orange"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(module_average_combined_1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(module_average_combined_1)/paletteLength, max(module_average_combined_1), length.out=floor(paletteLength/2)))


pheatmap(module_average_combined_1,cluster_cols=F,cluster_rows=F,
         color=myColor, breaks=myBreaks)


## use pseudobulk gene expression values to examine individual UPR marker genes
pseudobulk_ERstress <- AggregateExpression(ERstress.integrated, assays = "RNA", return.seurat = T, group.by = c("cell.type.final", "sex","condition"))
pseudobulk_ERstress_1 <- pseudobulk_ERstress[genelist,]
Idents(pseudobulk_ERstress_1) = pseudobulk_ERstress_1$cell.type.final

## focusing on beta cells
pseudobulk_ERstress_beta = subset(pseudobulk_ERstress_1, idents="beta")
temp <- pseudobulk_ERstress_beta@assays$RNA$data
beta.meta.data <- pseudobulk_ERstress_beta@meta.data
pathways <- ifelse(rownames(temp)%in%PERK_pathway,"PERK_pathway",ifelse(rownames(temp)%in%IRE1_pathway,"IRE1_pathway",ifelse(rownames(temp)%in%ATF6_pathway,"ATF6_pathway","NONE")))
pathways <- data.frame(pathways)
rownames(pathways) = rownames(temp)


ann_colors = list(
  sex = c(female = "#7FC97F", male = "#BEAED4"),
  condition=c("0h-Tg"="#FFFFFF","12h-Tg"="#A6CEE3","48h-Tg"="#1F78B4"),
  celltype=c(alpha="#88BDFF", beta="#FF8000",delta="#FFFA0A"))


p <- pheatmap(temp,
              main = "beta",
              cluster_cols = FALSE, 	
              clustering_method="ward.D2",
              #clustering_distance_rows = "correlation",
              annotation = beta.meta.data[,c(3,4)],
              annotation_row = pathways,
              scale="row",
              annotation_colors=ann_colors,
              color = colorRampPalette(c("grey", "white", "firebrick"))(100)
)
# add colors to row dendrogram
row_dend <- p$tree_row
# Cut the dendrogram at the second branch level to assign clusters
clusters <- cutree(row_dend, k = 2)  # k = 2 splits into two major branches
# Create a data frame for row annotations based on clusters
row_annotation <- cbind(data.frame(Cluster = factor(clusters)),pathways)
# Define colors for the second branch level clusters
cluster_colors <- c("1" = "#A6761D", "2" = "#E7298A")  # Adjust colors as needed
ann_colors = list(
  sex = c(female = "#7FC97F", male = "#BEAED4"),
  condition=c("0h-Tg"="#FFFFFF","12h-Tg"="#A6CEE3","48h-Tg"="#1F78B4"),
  celltype=c(alpha="#88BDFF", beta="#FF8000",delta="#FFFA0A"),
  Cluster=cluster_colors)


pheatmap(temp,
         main = "beta",
         cluster_cols = FALSE, 	
         clustering_method="ward.D2",
         annotation = beta.meta.data[,c(3,4)],
         annotation_row = row_annotation,
         scale="row",
         annotation_colors=ann_colors,
         color = colorRampPalette(c("grey", "white", "firebrick"))(100)
)

