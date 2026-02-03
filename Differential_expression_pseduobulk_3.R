library(Seurat)
# differential expression analysis with pseudobulk approach
ERstress.integrated
## within each sex
## beta
beta = subset(ERstress.integrated,idents="beta")
Idents(beta)=beta$sex
pseudo_beta <- AggregateExpression(beta, assays = "RNA", return.seurat = T, group.by = c("sex","condition",  "donor"))
pseudo_beta$sex.condition <- paste(pseudo_beta$sex, pseudo_beta$condition, sep = "_")
Idents(pseudo_beta) <- "sex.condition"

beta.bulk.de.female.12h.0h <- FindMarkers(object = pseudo_beta, 
                                          ident.1 = "female_12h-Tg",
                                          ident.2 = "female_0h-Tg", 
                                          logfc.threshold=0,
                                          min.pct=0,
                                          test.use = "DESeq2")
beta.bulk.de.female.12h.0h$p_val_adj <- p.adjust(beta.bulk.de.female.12h.0h$p_val, method='fdr')
beta.bulk.de.female.12h.0h <- na.omit(beta.bulk.de.female.12h.0h)

beta.bulk.de.female.48h.0h <- FindMarkers(object = pseudo_beta, 
                                          ident.1 = "female_48h-Tg",
                                          ident.2 = "female_0h-Tg", 
                                          logfc.threshold=0,
                                          min.pct=0,
                                          test.use = "DESeq2",
                                          min.cells.group=2)
beta.bulk.de.female.48h.0h$p_val_adj <- p.adjust(beta.bulk.de.female.48h.0h$p_val, method='fdr')
beta.bulk.de.female.48h.0h <- na.omit(beta.bulk.de.female.48h.0h)

beta.bulk.de.female.48h.12h <- FindMarkers(object = pseudo_beta, 
                                           ident.1 = "female_48h-Tg",
                                           ident.2 = "female_12h-Tg", 
                                           logfc.threshold=0,
                                           min.pct=0,
                                           test.use = "DESeq2",
                                           min.cells.group=2)
beta.bulk.de.female.48h.12h$p_val_adj <- p.adjust(beta.bulk.de.female.48h.12h$p_val, method='fdr')
beta.bulk.de.female.48h.12h <- na.omit(beta.bulk.de.female.48h.12h)

beta.bulk.de.male.12h.0h <- FindMarkers(object = pseudo_beta, 
                                        ident.1 = "male_12h-Tg",
                                        ident.2 = "male_0h-Tg", 
                                        logfc.threshold=0,
                                        min.pct=0,
                                        test.use = "DESeq2")
beta.bulk.de.male.12h.0h$p_val_adj <- p.adjust(beta.bulk.de.male.12h.0h$p_val, method='fdr')
beta.bulk.de.male.12h.0h <- na.omit(beta.bulk.de.male.12h.0h)

beta.bulk.de.male.48h.0h <- FindMarkers(object = pseudo_beta, 
                                        ident.1 = "male_48h-Tg",
                                        ident.2 = "male_0h-Tg", 
                                        logfc.threshold=0,
                                        min.pct=0,
                                        test.use = "DESeq2",
                                        min.cells.group=2)
beta.bulk.de.male.48h.0h$p_val_adj <- p.adjust(beta.bulk.de.male.48h.0h$p_val, method='fdr')
beta.bulk.de.male.48h.0h <- na.omit(beta.bulk.de.male.48h.0h)

beta.bulk.de.male.48h.12h <- FindMarkers(object = pseudo_beta, 
                                         ident.1 = "male_48h-Tg",
                                         ident.2 = "male_12h-Tg", 
                                         logfc.threshold=0,
                                         min.pct=0,
                                         test.use = "DESeq2",
                                         min.cells.group=2)
beta.bulk.de.male.48h.12h$p_val_adj <- p.adjust(beta.bulk.de.male.48h.12h$p_val, method='fdr')
beta.bulk.de.male.48h.12h <- na.omit(beta.bulk.de.male.48h.12h)

write.table(beta.bulk.de.female.12h.0h,"beta.bulk.de.female.12h.0h.csv")
write.table(beta.bulk.de.female.48h.0h,"beta.bulk.de.female.48h.0h.csv")
write.table(beta.bulk.de.female.48h.12h,"beta.bulk.de.female.48h.12h.csv")
write.table(beta.bulk.de.male.12h.0h,"beta.bulk.de.male.12h.0h.csv")
write.table(beta.bulk.de.male.48h.0h,"beta.bulk.de.male.48h.0h.csv")
write.table(beta.bulk.de.male.48h.12h,"beta.bulk.de.male.48h.12h.csv")



# GSEA

library(fgsea)
library(data.table)
library(ggplot2)
library(msigdbr)

hallmarks <- gmtPathways("h.all.v2025.1.Hs.symbols.gmt")
beta.bulk.de.female.12h.0h$rank <- sign(beta.bulk.de.female.12h.0h$avg_log2FC) * (-log10(beta.bulk.de.female.12h.0h$p_val))

rl1 <- beta.bulk.de.female.12h.0h$rank
names(rl1) <- rownames(beta.bulk.de.female.12h.0h)

fgseaRes <- fgsea(
  pathways = hallmarks,
  stats = rl1,
  minSize = 15,
  maxSize = 5000,
  eps = 0.0
)
fwrite(fgseaRes, file="beta.bulk.de.female.12h.0h.GSEA.txt", sep="\t", sep2=c("", " ", ""))


sigPathways <- fgseaRes[padj < 0.05][order(pval), pathway]

plotGseaTable(hallmarks[sigPathways], rl1, fgseaRes, 
              gseaParam=0.5)


## Select the union of significant terms (FDR < 0.05).
female.12h.pathways <- read.table("beta.bulk.de.female.12h.0h.GSEA.txt",header=T,sep="\t") 
female.12h.pathways <- na.omit(female.12h.pathways)
female.12h.pathways.terms <- female.12h.pathways[(female.12h.pathways$padj < 0.05),]$pathway

female.48h.pathways <- read.table("beta.bulk.de.female.48h.0h.GSEA.txt",header=T,sep="\t") 
female.48h.pathways <- na.omit(female.48h.pathways)
female.48h.pathways.terms <- female.48h.pathways[(female.48h.pathways$padj < 0.05),]$pathway

male.12h.pathways <- read.table("beta.bulk.de.male.12h.0h.GSEA.txt",header=T,sep="\t") 
male.12h.pathways <- na.omit(male.12h.pathways)
male.12h.pathways.terms <- male.12h.pathways[(male.12h.pathways$padj < 0.05),]$pathway

male.48h.pathways <- read.table("beta.bulk.de.male.48h.0h.GSEA.txt",header=T,sep="\t") 
male.48h.pathways <- na.omit(male.48h.pathways)
male.48h.pathways.terms <- male.48h.pathways[(male.48h.pathways$padj < 0.05),]$pathway

de.terms = Reduce(union,list(female.12h.pathways.terms,female.48h.pathways.terms,male.12h.pathways.terms,male.48h.pathways.terms))
female.12h.pathways.1 <- female.12h.pathways[female.12h.pathways$pathway%in%de.terms,]
female.12h.pathways.1 <- female.12h.pathways.1[,c("pathway","padj")]
female.48h.pathways.1 <- female.48h.pathways[female.48h.pathways$pathway%in%de.terms,]
female.48h.pathways.1 <- female.48h.pathways.1[,c("pathway","padj")]
male.12h.pathways.1 <- male.12h.pathways[male.12h.pathways$pathway%in%de.terms,]
male.12h.pathways.1 <- male.12h.pathways.1[,c("pathway","padj")]
male.48h.pathways.1 <- male.48h.pathways[male.48h.pathways$pathway%in%de.terms,]
male.48h.pathways.1 <- male.48h.pathways.1[,c("pathway","padj")]

combined.pathways <- merge(female.12h.pathways.1,male.12h.pathways.1,by="pathway",all=TRUE)
combined.pathways <- merge(combined.pathways,female.48h.pathways.1,by="pathway",all=TRUE)
combined.pathways <- merge(combined.pathways,male.48h.pathways.1,by="pathway",all=TRUE)
rownames(combined.pathways) = combined.pathways$pathway
combined.pathways = combined.pathways[,-1]
colnames(combined.pathways) =c("female.12h","male.12h","female.48h","male.48h")
# recoding NA to 1
combined.pathways[is.na(combined.pathways)] <- 1
combined.pathways.1 <- -log10(combined.pathways)
#pheatmap(combined.pathways.1, 
#         clustering_method="ward.D2",
#         cluster_cols = FALSE
#         )
combined.pathways.1 <- combined.pathways.1[,c(1,3,2,4)]
library(gplots)
heatmap.2(data.matrix(combined.pathways.1), 
          dendrogram = "none",
          trace="none",
          density.info="none",
          scale="none", 
          col=colorRampPalette(c("white","red"))(n=10), 
          breaks=c(0,1,2,3,4,5,6,7,8,9,10),
          key=TRUE,
          cexCol = 1, cexRow=1, 
          Rowv=NULL, 
          Colv=NULL,
          sepwidth=c(0.02,0.02),
          sepcolor="black",
          colsep=0:ncol(combined.pathways.1),
          rowsep=0:nrow(combined.pathways.1))


## the same scripts used for other cell types for pathway analysis and visualization

## combining alpha, beta and delta cells 
## top 5 most significant pathways in alpha, beta, and delta cells based on p.adjust < 0.05 and select union of all pathways, creating corresponding dataframe

get_top_pathways <- function(
    celltypes  = c("beta", "alpha", "delta"),
    sexes      = c("female", "male"),
    contrasts  = c("12h.0h", "48h.0h"),
    dir        = ".",                    # directory containing files
    outfile    = "top5_pathways_summary.csv",
    top_n      = 5
) {
  build_filename <- function(cell, sex, contrast) {
    sprintf("%s/%s.bulk.de.%s.%s.GSEA.txt", dir, cell, sex, contrast)
  }
  
  combos <- expand.grid(cell = celltypes, sex = sexes, contrast = contrasts, stringsAsFactors = FALSE)
  
  results_list <- list()   # top 5 per combo
  full_list    <- list()   # FULL pathway+padj+NES per combo (deduplicated)
  summary_rows <- list()   # long-format top-5 summary
  
  for (i in seq_len(nrow(combos))) {
    cell     <- combos$cell[i]
    sex      <- combos$sex[i]
    contrast <- combos$contrast[i]
    fname    <- build_filename(cell, sex, contrast)
    key_name <- paste(cell, sex, contrast, sep = ".")
    
    # Read table safely
    df <- tryCatch(
      read.table(fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) {
        warning(sprintf("Failed to read %s: %s (skipping)", fname, e$message))
        return(NULL)
      }
    )
    if (is.null(df)) next
    
    
    # Ensure padj is numeric, then keep only needed cols and drop NAs
    if (!is.numeric(df$padj)) {
      suppressWarnings(df$padj <- as.numeric(df$padj))
    }
    df <- na.omit(df[, c("pathway", "padj", "NES")])
    if (nrow(df) == 0) next
    
    # Deduplicate same pathway per file using min(padj) (change FUN if needed)
    df <- df[order(df$pathway, df$padj, na.last = NA), ]
    df_all <- df[!duplicated(df$pathway), c("pathway", "padj", "NES")]
    full_list[[key_name]] <- df_all
    
    # Order for top N
    ord_idx <- order(df_all$padj, na.last = NA)
    df_ord  <- df_all[ord_idx, , drop = FALSE]
    top_df  <- head(df_ord[, c("pathway", "padj","NES")], top_n)
    results_list[[key_name]] <- top_df
    
    # Add to summary rows (long format)
    if (nrow(top_df) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        celltype = cell,
        sex      = sex,
        contrast = contrast,
        rank     = seq_len(nrow(top_df)),   # <-- rows, not columns
        pathway  = top_df$pathway,
        padj     = top_df$padj,
        NES		 = top_df$NES,
        stringsAsFactors = FALSE
      )
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        celltype = cell,
        sex      = sex,
        contrast = contrast,
        rank     = integer(0),
        pathway  = character(0),
        padj     = numeric(0),
        NES		 = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Bind and write summary
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df, file = outfile, row.names = FALSE)
  } else {
    summary_df <- data.frame(
      celltype = character(0), sex = character(0), contrast = character(0),
      rank = integer(0), pathway = character(0), padj = numeric(0), NES = numberic(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    top_list  = results_list,
    full_list = full_list,
    summary   = summary_df
  ))
}




# Helper: construct a wide NES matrix for any set of pathway names
NES_matrix_from_full <- function(full_list, pathways = NULL, fill = NA_real_) {
  if (is.null(pathways)) {
    pathways <- unique(unlist(lapply(full_list, function(df) df$pathway)))
  }
  combo_names <- names(full_list)
  
  mat <- matrix(fill, nrow = length(pathways), ncol = length(combo_names),
                dimnames = list(pathways, combo_names))
  
  for (j in seq_along(combo_names)) {
    df <- full_list[[j]]
    if (is.null(df) || nrow(df) == 0) next
    if (!all(c("pathway", "NES") %in% names(df))) next
    
    m <- match(df2$pathway, pathways)
    keep <- !is.na(m)
    mat[m[keep], j] <- df2$NES[keep]
  }
  
  mat
}


padj_matrix_from_full <- function(full_list, pathways = NULL, agg_fun = min, fill = NA_real_) { 
  if (is.null(pathways)) { 
    pathways <- unique(unlist(lapply(full_list, function(df) df$pathway))) 
  }
  combo_names <- names(full_list)
  
  mat <- matrix(fill, nrow = length(pathways), ncol = length(combo_names),
                dimnames = list(pathways, combo_names))
  
  for (j in seq_along(combo_names)) {
    df <- full_list[[j]]
    if (is.null(df) || nrow(df) == 0) next
    if (!all(c("pathway", "padj") %in% names(df))) next
    
    m <- match(df2$pathway, pathways)
    keep <- !is.na(m)
    mat[m[keep], j] <- df2$NES[keep]
  }
  
  mat
}

res <- get_top_pathways(
  celltypes = c("beta", "alpha", "delta"),
  sexes     = c("female", "male"),
  contrasts = c("12h.0h", "48h.0h"),
  dir       = "/gpfs/research/jwanggroup/jwang_group/scRNA-seq/ERstress/integration",
  outfile   = "top5_pathways_summary.csv",
  top_n     = 5
)

# Your union of top pathways (can also be any list you define)
combined.pathways <- unique(res$summary$pathway)

# Build the wide padj and NES matrix across ALL combos using FULL tables
padj_mat <- padj_matrix_from_full(res$full_list, pathways = combined.pathways)
NES_mat <- NES_matrix_from_full(res$full_list, pathways = combined.pathways)

# Save
write.csv(NES_mat, "NES_mat_combined_pathways.csv")

# Plot combined heatmap
combined.pathways <- NES_mat


rownames(combined.pathways) = sub("^.*?\\_", "", rownames(combined.pathways))
combined.pathways = combined.pathways[,c("alpha.female.12h.0h","alpha.female.48h.0h","alpha.male.12h.0h","alpha.male.48h.0h",
                                         "beta.female.12h.0h","beta.female.48h.0h","beta.male.12h.0h","beta.male.48h.0h",
                                         "delta.female.12h.0h","delta.female.48h.0h","delta.male.12h.0h","delta.male.48h.0h")]
my_colour = list(
  sex = c(female = "#7FC97F", male = "#BEAED4"),
  condition=c(Tg_12h="#A6CEE3",Tg_48h="#1F78B4"),
  celltype=c(alpha="#88BDFF", beta="#FF8000",delta="#FFFA0A"))

my_sample_col <- data.frame(sex = rep(c("female","female","male","male"),3),condition=rep(c("Tg_12h","Tg_48h"),6),
                            celltype=c(rep("alpha",4),rep("beta",4),rep("delta",4)))

rownames(my_sample_col)=colnames(combined.pathways)

# recoding NA to 0
combined.pathways[is.na(combined.pathways)] <- 0
#combined.pathways.1 <- -log10(combined.pathways)


# Define colors and breaks
color_palette <- c(colorRampPalette(c("blue","white"))(5),colorRampPalette(c("white","red"))(5))
breaks <- seq(-3, 3, length.out = 11)

library(pheatmap)
# Create heatmap
pheatmap(
  combined.pathways,  
  annotation_col = my_sample_col,
  annotation_colors=my_colour,
  cluster_rows = TRUE,   # Disable row clustering (reordered manually)
  cluster_cols = FALSE,  # Disable column clustering
  clustering_method = "ward.D2",
  color = color_palette,
  breaks = breaks,
  border_color = "black",  # Add separators
  cellwidth = 15,          # Adjust column width
  cellheight = 15,         # Adjust row height
  fontsize_row = 10,       # Adjust row font size
  fontsize_col = 10,       # Adjust column font size
  legend = TRUE,           # Show color key
  legend_breaks = breaks,  # Ensure legend matches breaks
  legend_labels = breaks,  # Label each break
  annotation_names_row = TRUE,  # Show row annotations
  annotation_names_col = TRUE   # Show column annotations
)


# adding blank space between celltypes

# Add blank columns between groups
combined.pathways.1_with_space <- combined.pathways
group_indices <- c(4, 8)  # Indices where groups end (modify as needed)

for (idx in rev(group_indices)) {
  blank_col <- matrix(NA, nrow = nrow(combined.pathways.1_with_space), ncol = 1)
  combined.pathways.1_with_space <- cbind(
    combined.pathways.1_with_space[, 1:idx],
    blank_col,
    combined.pathways.1_with_space[, (idx + 1):ncol(combined.pathways.1_with_space)]
  )
}

# Update annotations to match the new matrix
my_sample_col_with_space <- my_sample_col
for (idx in rev(group_indices)) {
  blank_row <- data.frame(sex = NA, condition = NA, celltype = NA)
  rownames(blank_row) <- paste0("Spacer", idx)  # Assign unique rownames for spacers
  my_sample_col_with_space <- rbind(
    my_sample_col_with_space[1:idx, ],
    blank_row,
    my_sample_col_with_space[(idx + 1):nrow(my_sample_col_with_space), ]
  )
}


# Custom border colors: NA for added spaces, black for the rest
custom_border_colors <- matrix("black", nrow = nrow(combined.pathways.1_with_space), 
                               ncol = ncol(combined.pathways.1_with_space))
blank_cols <- which(apply(is.na(combined.pathways.1_with_space), 2, all))
custom_border_colors[, blank_cols] <- NA

# Create heatmap
pheatmap(
  combined.pathways.1_with_space,  
  annotation_col = my_sample_col_with_space,
  annotation_colors = my_colour,
  cluster_rows = TRUE,   # Enable row clustering
  cluster_cols = FALSE,  # Disable column clustering
  #  clustering_method = "ward.D2",
  color = color_palette,
  breaks = breaks,
  border_color = custom_border_colors,  # Custom grid line colors
  cellwidth = 15,           # Adjust column width
  cellheight = 15,          # Adjust row height
  fontsize_row = 10,        # Adjust row font size
  fontsize_col = 10,        # Adjust column font size
  legend = TRUE,            # Show color key
  legend_breaks = breaks,   # Ensure legend matches breaks
  legend_labels = breaks,   # Label each break
  annotation_names_row = TRUE,  # Show row annotations
  annotation_names_col = TRUE   # Show column annotations
)




