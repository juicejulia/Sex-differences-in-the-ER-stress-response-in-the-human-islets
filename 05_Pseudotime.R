# pseudotime
## Pseudotime analysis with slingshot
## Adapted from this tutorial:https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
#and this: https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html

library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(scales)
library(viridis)
library(UpSetR)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(knitr)
library(ggplot2)
library(gridExtra)
library(tradeSeq)
library(cowplot)

## focusing on beta cells
## further clean up the beta cell population and only retain those cells with consensus annotation between our own pipeline and azimuth

Idents(ERstress.integrated)=ERstress.integrated$cell.type.final
beta = subset(ERstress.integrated,idents="beta")

Idents(beta) = beta$cell.type.azimuth
beta.sub = subset(beta,idents="beta")

DimPlot(beta, reduction = "umap.rpca",  group.by="sex",label = FALSE, pt.size = 0.5,raster=FALSE,alpha=0.4) 
DimPlot(beta.sub, reduction = "umap.rpca",  group.by="sex",label = FALSE, pt.size = 0.5,raster=FALSE,alpha=0.4) 



# Convert to singleCellExperiment
beta_sce <- as.SingleCellExperiment(beta.sub, assay = "RNA")

library(slingshot)
beta_sce = slingshot(beta_sce,reducedDim='UMAP.RPCA', clusterLabels=colData(beta_sce)$condition,start.clus="0h_Tg", dist.method="mnn",approx_points=100)

# Differential topology test
set.seed(821)
library(condiments)
library(dplyr)
BPPARAM <- BiocParallel::MulticoreParam(workers = 20)

topologyTest(SlingshotDataSet(beta_sce), beta_sce$sex, rep = 100,
             methods = "KS_mean", threshs = .01, paralle=T)
#Generating permuted trajectories
#Running KS-mean test
#   method thresh  statistic    p.value
#1 KS_mean   0.01  0.02511581   0.00542437


library(grDevices)
library(fields)
grad <- viridis::plasma(100, begin = 0, end = 1)
plotcol <- grad[cut(beta_sce$slingPseudotime_1, breaks=100)]
pt <- beta_sce$slingPseudotime_1
plot(reducedDims(beta_sce)$UMAP.RPCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(beta_sce), lwd=2, col='black')
image.plot(legend.only = TRUE, zlim = range(pt, na.rm = TRUE), col = grad,
           legend.lab = "Pseudotime", legend.width = 1.2)



#boxplot comparing the pseudotimes of samples with different Tg treatment
library(ggpmisc)
library(ggpubr)
library(ggplot2)

meta.data.beta <-data.frame(cbind(beta_sce$condition,beta_sce$sex,beta_sce$slingPseudotime_1))
rownames(meta.data.beta)=colnames(beta_sce)
colnames(meta.data.beta)=c("condition","sex","slingPseudotime_1")
meta.data.beta$condition = factor(meta.data.beta$condition)
meta.data.beta$sex = factor(meta.data.beta$sex)
meta.data.beta$slingPseudotime_1 = as.numeric(meta.data.beta$slingPseudotime_1)
my_comparisons <- list(c('0h_Tg','12h_Tg'),c('0h_Tg','48h_Tg'),c('12h_Tg','48h_Tg'))

pdf('test.pdf',height=7,width=3)
ggplot(meta.data.beta, aes(x=condition,y=slingPseudotime_1,fill=condition)) +     
  geom_boxplot(outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,15))+
  #stat_compare_means(comparisons=my_comparisons)+
  theme(text = element_text(size=20)) + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


# Differential progression between different sexes
ggplot(df, aes(x = slingPseudotime_1)) +
  geom_density(alpha = .4, aes(fill = sex), col = "transparent") +
  geom_density(aes(col = sex), fill = "transparent",
               linewidth = 1.5) +
  labs(x = "Pseudotime", fill = "Type") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  facet_wrap(~condition)

ggplot(df, aes(x = slingPseudotime_1)) +
  geom_density(alpha = .4, aes(fill = condition), col = "transparent") +
  geom_density(aes(col = condition), fill = "transparent", linewidth = 1.5) +
  guides(col = "none") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  labs(x = "Pseudotime", fill = "Type") +
  facet_wrap(~sex)

library(tradeSeq)

# fit negative binomial GAM
# filter genes to have minimum shared counts of 20, and minimum expressed in 3 cells, in line with scVelo, and to improve fitGAM execusion speed
gene_cells <- rowSums(counts(beta_sce)!=0)
beta_sce_filtered <- beta_sce[gene_cells > 3,]
gene_sums <- rowSums(counts(beta_sce_filtered))
beta_sce_filtered <-  beta_sce_filtered[gene_sums > 20,]



icMat <- evaluateK(counts=beta_sce_filtered,
                   nGenes = 500,
                   k = 3:10, 
                   conditions=factor(beta_sce$sex),
                   parallel=TRUE)

plot_evalutateK_results(icMat,k=3:10)

beta_sce_filtered <- fitGAM(beta_sce_filtered, nknots=6, conditions=factor(beta_sce_filtered$sex),parallel=TRUE) 
# fit a general additive model 


# Differential expression along pseudotime
ATres <- associationTest(beta_sce_filtered, lineages=TRUE, l2fc = log2(2)) # testing whether the average gene expression is significantly changed along pseudotime
ATres$p.adjust_lineage1_conditionfemale <- p.adjust(ATres$pvalue_lineage1_conditionfemale, "fdr")
ATres$p.adjust_lineage1_conditionmale <- p.adjust(ATres$pvalue_lineage1_conditionmale, "fdr")
write.table(ATres,"genes_associated_with_pseudotime.csv",quote=F,sep=",")


FemaleGenes <-  rownames(ATres)[
  which(ATres$p.adjust_lineage1_conditionfemale <= 0.05)
]
MaleGenes <-  rownames(ATres)[
  which(ATres$p.adjust_lineage1_conditionmale <= 0.05)
]

length(FemaleGenes)
## [1] 4798
length(MaleGenes)
## [1] 7416
UpSetR::upset(fromList(list(Female = FemaleGenes, Male = MaleGenes)))


# Differential expression between two sexes 
condRes <- conditionTest(beta_sce_filtered)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
condRes$padj1 <- p.adjust(condRes$pvalue, "bonferroni") # use a more stringent p value correction method

write.table(condRes,"genes_differential_expression_between_two_sexes_logFC0.csv",quote=F,sep=",")


sum(condRes$padj1 <= 0.01, na.rm = TRUE)
# 4361
conditionGenes <- rownames(condRes)[condRes$padj1 <= 0.01]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]



# Genes in the "unfolded protein response" term that display differential dynamics between females and males
#genelist=c("EIF4A2","ERO1A","EIF4A1","EIF4A3","CEBPG","SDAD1","RRP9","HSP90B1","HERPUD1","CKS1B","RPS14","EXOSC5","SERP1","KIF5B","DNAJB9","BANF1","EIF4E","SEC11A","HSPA9","XBP1","NPM1","HSPA5","IMP3","EEF2","PDIA6","EIF2S1","LSM4","YIF1A","VEGFA","SPCS3","SPCS1","MTHFD2","SRPRA","SRPRB","CALR","ATP6V0D1","ATF4")

# Genes in the UPR term, matching the pseudobulk heatmap
genelist= c("TMED2","ATF6","HERPUD1","PTPN1","NFE2L2","XBP1","DNAJC10","ATF4","DDIT3","HSPA5",
            "ERN1","CREBZF","EIF2S1","RPAP2","QRICH1","EIF2AK3","MBTPS2","PARP16","ATF6B","MBTPS1","VAPB")
scales <- brewer.pal(7,"Accent")[1:2]
pdf("beta_selected_genes_UPR.pdf")
for (i in genelist)
{
  gene = i
  print(plotSmoothers(beta_sce_filtered, assays(beta_sce_filtered)$counts,
                      gene = gene,
                      alpha = 1, border = TRUE, curvesCols=scales)+
          scale_color_manual(values=scales)+
          ggtitle(gene))
}
dev.off()

# visualization of the dynamics of UPR genes of females and males with heatmap
yhatSmooth <- 
  predictSmooth(beta_sce_filtered, gene = genelist, nPoints = 100, tidy = FALSE) 
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))
heatSmooth_female <- pheatmap(yhatSmoothScaled[, 1:100],
                              cluster_cols = FALSE,
                              #cluster_rows = FALSE,
                              show_rownames = TRUE, 
                              show_colnames = FALSE, 
                              main = "female", 	
                              clustering_method="ward.D2",
                              annotation_colors=my_colour,
                              legend = FALSE,
                              silent = TRUE,
                              border_color = NA
)
# add colors to row dendrogram
row_dend <- heatSmooth_female$tree_row
# Cut the dendrogram at the fourth branch level to assign clusters
clusters <- cutree(row_dend, k = 4)  
row_annotation <- data.frame(clusters)
my_colour = list(
  clusters = brewer.pal(7,"Accent")[1:4]
)
heatSmooth_female <- pheatmap(yhatSmoothScaled[, 1:100],
                              cluster_cols = FALSE,
                              #cluster_rows = FALSE,
                              show_rownames = TRUE, 
                              show_colnames = FALSE, 
                              main = "female", 	
                              clustering_method="ward.D2",
                              annotation_colors=my_colour,
                              legend = FALSE,
                              silent = TRUE,
                              border_color = NA,
                              annotation_row = row_annotation,
                              annotation_legend=FALSE,
                              annotation_names_row=FALSE
)

matchingHeatmap_male <- 
  pheatmap(yhatSmoothScaled[heatSmooth_female$tree_row$order, 101:200],
           cluster_cols = FALSE, cluster_rows = FALSE,
           show_rownames = TRUE, show_colnames = FALSE, main = "male",
           legend = TRUE, silent = TRUE, border_color = NA
  )


plot_grid(heatSmooth_female[[4]], matchingHeatmap_male[[4]], ncol = 2)

saveRDS(beta_sce,"beta_sub_sce_filtered_MNN.rds")

