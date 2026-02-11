# Clean up the dataframes
HP23166_0h <- read.table("HP23166_0h_S1_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(HP23166_0h) <- sub(";.*$", "", HP23166_0h$BC)
HP23166_0h_1 <- HP23166_0h[(HP23166_0h$WT.calls!=0)|(HP23166_0h$MUT.calls!=0),]
HP23166_0h_1 <- HP23166_0h_1[,c("WT.calls","MUT.calls")]
write.table(HP23166_0h_1,"XBP1_GoT/HP23166_0h.csv",quote=F,sep=",")

HP23166_12h <- read.table("HP23166_12h_S2_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(HP23166_12h) <- sub(";.*$", "", HP23166_12h$BC)
HP23166_12h_1 <- HP23166_12h[(HP23166_12h$WT.calls!=0)|(HP23166_12h$MUT.calls!=0),]
HP23166_12h_1 <- HP23166_12h_1[,c("WT.calls","MUT.calls")]
write.table(HP23166_12h_1,"XBP1_GoT/HP23166_12h.csv",quote=F,sep=",")

SAMN35848421_0h <- read.table("SAMN35848421_0h_S3_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_0h) <- sub(";.*$", "", SAMN35848421_0h$BC)
SAMN35848421_0h_1 <- SAMN35848421_0h[(SAMN35848421_0h$WT.calls!=0)|(SAMN35848421_0h$MUT.calls!=0),]
SAMN35848421_0h_1 <- SAMN35848421_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_0h_1,"XBP1_GoT/SAMN35848421_0h.csv",quote=F,sep=",")

SAMN35848421_12h <- read.table("SAMN35848421_12h_S4_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_12h) <- sub(";.*$", "", SAMN35848421_12h$BC)
SAMN35848421_12h_1 <- SAMN35848421_12h[(SAMN35848421_12h$WT.calls!=0)|(SAMN35848421_12h$MUT.calls!=0),]
SAMN35848421_12h_1 <- SAMN35848421_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_12h_1,"XBP1_GoT/SAMN35848421_12h.csv",quote=F,sep=",")

SAMN35848421_48h <- read.table("SAMN35848421_48h_S5_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_48h) <- sub(";.*$", "", SAMN35848421_48h$BC)
SAMN35848421_48h_1 <- SAMN35848421_48h[(SAMN35848421_48h$WT.calls!=0)|(SAMN35848421_48h$MUT.calls!=0),]
SAMN35848421_48h_1 <- SAMN35848421_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_48h_1,"XBP1_GoT/SAMN35848421_48h.csv",quote=F,sep=",")

SAMN36705973_0h <- read.table("SAMN36705973_0h_S6_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_0h) <- sub(";.*$", "", SAMN36705973_0h$BC)
SAMN36705973_0h_1 <- SAMN36705973_0h[(SAMN36705973_0h$WT.calls!=0)|(SAMN36705973_0h$MUT.calls!=0),]
SAMN36705973_0h_1 <- SAMN36705973_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_0h_1,"XBP1_GoT/SAMN36705973_0h.csv",quote=F,sep=",")

SAMN36705973_12h <- read.table("SAMN36705973_12h_S7_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_12h) <- sub(";.*$", "", SAMN36705973_12h$BC)
SAMN36705973_12h_1 <- SAMN36705973_12h[(SAMN36705973_12h$WT.calls!=0)|(SAMN36705973_12h$MUT.calls!=0),]
SAMN36705973_12h_1 <- SAMN36705973_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_12h_1,"XBP1_GoT/SAMN36705973_12h.csv",quote=F,sep=",")

SAMN36705973_48h <- read.table("SAMN36705973_48h_S8_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_48h) <- sub(";.*$", "", SAMN36705973_48h$BC)
SAMN36705973_48h_1 <- SAMN36705973_48h[(SAMN36705973_48h$WT.calls!=0)|(SAMN36705973_48h$MUT.calls!=0),]
SAMN36705973_48h_1 <- SAMN36705973_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_48h_1,"XBP1_GoT/SAMN36705973_48h.csv",quote=F,sep=",")

SAMN36823227_0h <- read.table("SAMN36823227_0h_S9_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_0h) <- sub(";.*$", "", SAMN36823227_0h$BC)
SAMN36823227_0h_1 <- SAMN36823227_0h[(SAMN36823227_0h$WT.calls!=0)|(SAMN36823227_0h$MUT.calls!=0),]
SAMN36823227_0h_1 <- SAMN36823227_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_0h_1,"XBP1_GoT/SAMN36823227_0h.csv",quote=F,sep=",")

SAMN36823227_12h <- read.table("SAMN36823227_12h_S10_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_12h) <- sub(";.*$", "", SAMN36823227_12h$BC)
SAMN36823227_12h_1 <- SAMN36823227_12h[(SAMN36823227_12h$WT.calls!=0)|(SAMN36823227_12h$MUT.calls!=0),]
SAMN36823227_12h_1 <- SAMN36823227_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_12h_1,"XBP1_GoT/SAMN36823227_12h.csv",quote=F,sep=",")

SAMN36823227_48h <- read.table("SAMN36823227_48h_S11_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_48h) <- sub(";.*$", "", SAMN36823227_48h$BC)
SAMN36823227_48h_1 <- SAMN36823227_48h[(SAMN36823227_48h$WT.calls!=0)|(SAMN36823227_48h$MUT.calls!=0),]
SAMN36823227_48h_1 <- SAMN36823227_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_48h_1,"XBP1_GoT/SAMN36823227_48h.csv",quote=F,sep=",")

SAMN37871873_0h <- read.table("SAMN37871873_0h_S12_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_0h) <- sub(";.*$", "", SAMN37871873_0h$BC)
SAMN37871873_0h_1 <- SAMN37871873_0h[(SAMN37871873_0h$WT.calls!=0)|(SAMN37871873_0h$MUT.calls!=0),]
SAMN37871873_0h_1 <- SAMN37871873_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_0h_1,"XBP1_GoT/SAMN37871873_0h.csv",quote=F,sep=",")

SAMN37871873_12h <- read.table("SAMN37871873_12h_S13_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_12h) <- sub(";.*$", "", SAMN37871873_12h$BC)
SAMN37871873_12h_1 <- SAMN37871873_12h[(SAMN37871873_12h$WT.calls!=0)|(SAMN37871873_12h$MUT.calls!=0),]
SAMN37871873_12h_1 <- SAMN37871873_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_12h_1,"XBP1_GoT/SAMN37871873_12h.csv",quote=F,sep=",")

SAMN37871873_48h <- read.table("SAMN37871873_48h_S14_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_48h) <- sub(";.*$", "", SAMN37871873_48h$BC)
SAMN37871873_48h_1 <- SAMN37871873_48h[(SAMN37871873_48h$WT.calls!=0)|(SAMN37871873_48h$MUT.calls!=0),]
SAMN37871873_48h_1 <- SAMN37871873_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_48h_1,"XBP1_GoT/SAMN37871873_48h.csv",quote=F,sep=",")

SAMN39523303_0h <- read.table("SAMN39523303_0h_S15_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_0h) <- sub(";.*$", "", SAMN39523303_0h$BC)
SAMN39523303_0h_1 <- SAMN39523303_0h[(SAMN39523303_0h$WT.calls!=0)|(SAMN39523303_0h$MUT.calls!=0),]
SAMN39523303_0h_1 <- SAMN39523303_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_0h_1,"XBP1_GoT/SAMN39523303_0h.csv",quote=F,sep=",")

SAMN39523303_12h <- read.table("SAMN39523303_12h_S16_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_12h) <- sub(";.*$", "", SAMN39523303_12h$BC)
SAMN39523303_12h_1 <- SAMN39523303_12h[(SAMN39523303_12h$WT.calls!=0)|(SAMN39523303_12h$MUT.calls!=0),]
SAMN39523303_12h_1 <- SAMN39523303_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_12h_1,"XBP1_GoT/SAMN39523303_12h.csv",quote=F,sep=",")

SAMN39523303_48h <- read.table("SAMN39523303_48h_S17_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_48h) <- sub(";.*$", "", SAMN39523303_48h$BC)
SAMN39523303_48h_1 <- SAMN39523303_48h[(SAMN39523303_48h$WT.calls!=0)|(SAMN39523303_48h$MUT.calls!=0),]
SAMN39523303_48h_1 <- SAMN39523303_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_48h_1,"XBP1_GoT/SAMN39523303_48h.csv",quote=F,sep=",")


# exporting metadata from the scRNA-seq data

ERstress.integrated <- readRDS("ERstress.integrated.rds")
meta.data <- ERstress.integrated@meta.data

# prepare the individual GoT results for merging with meta.data
## Rename rownames based on filename to incorporate sample name and write out new CSVs

# Directory containing the input CSVs
dir_path <- "XBP1_GoT"

# Identify all CSV files
files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  
  # Extract base filename without extension
  base <- sub("\\.csv$", "", basename(f))
  
  # Construct prefix: e.g. "SAMN39523303_48h_Tg_"
  prefix <- paste0(base, "_Tg_")
  
  # Read CSV with rownames in first column
  df <- read.csv(
    f,
    header = TRUE,
    row.names = 1  
    )
  
  # Modify rownames: PREFIX + old_rowname + "-1"
  rownames(df) <- paste0(prefix, rownames(df), "-1")
  
  # Output filename
  out_file <- file.path(
    dir_path,
    paste0(base, "_renamed.csv")
  )
  
  # Write CSV
  write.csv(
    df,
    file = out_file,
    quote = FALSE
  )
}

## append all the files together 
dir_path <- "XBP1_GoT"

# Identify all CSV files
files <- list.files(path = dir_path, pattern = "\\_renamed.csv$", full.names = TRUE)
combined <- data.frame()
for (f in files) {
  df <- read.csv(
    f,
    header = TRUE,
    row.names = 1  
  )
  combined <- rbind(combined,df)
}
write.table(combined,file.path(dir_path,"combined.csv"),quote=F,sep=",")

## put everything together
# Left outer join: keeps all rows from df1 and matching rows from df2 (fills NAs for non-matches)
meta.data$cell <- rownames(meta.data)
combined$cell <- rownames(combined)
merged <- merge(meta.data, combined, by = "cell", all.x = TRUE) # by=0 is shorthand for by="row.names"
rownames(merged) = merged$cell

# compute normalized counts for uXBP1 and sXBP1, using the LogNormalize method as in Seurat
## WT call corresponding to uXBP1, MUT call corresponding to sXBP1
merged$norm.WT.calls <- log1p(merged$WT.calls/merged$nCount_RNA * 10000)
merged$norm.MUT.calls <- log1p(merged$MUT.calls/merged$nCount_RNA * 10000)

merged_1 <- na.omit(merged)
write.table(merged,"meta.data.GoT.merged.csv",quote=F,sep=",")

ERstress.integrated[["uXBP1"]] = merged$WT.calls
ERstress.integrated[["sXBP1"]] = merged$MUT.calls
ERstress.integrated[["norm.uXBP1"]] = merged$norm.WT.calls
ERstress.integrated[["norm.sXBP1"]] = merged$norm.MUT.calls

## filter the Seurat data to retain only those with corresponding to GoT data
ERstress.integrated[["CellName"]] = colnames(ERstress.integrated)
ERstress.integrated.sub <- subset(ERstress.integrated, subset = CellName %in% rownames(merged_1) )


# Merge with pseudotime data
beta <- subset(ERstress.integrated.sub,idents="beta")
beta_sce <- readRDS("beta_sub_sce_filtered_MNN.rds")
pseudotime <- data.frame(colnames(beta_sce),beta_sce$slingPseudotime_1)
common_cells <- intersect(beta$CellName, pseudotime$colnames.beta_sce.)
beta.sub <- subset(beta, subset = CellName %in% common_cells)
rownames(pseudotime) = pseudotime$colnames.beta_sce.
pseudotime.sub <- pseudotime[colnames(beta.sub),]
beta.sub$pseudotime <- pseudotime.sub$beta_sce.slingPseudotime_1




# plotting, the plotting method used
library(RColorBrewer)
meta.beta <- beta.sub@meta.data
scales <- brewer.pal(7,"Accent")[1:2]

ggplot(meta.beta, aes(x = pseudotime, y = norm.uXBP1, color=sex)) +
 # geom_point() +
  geom_smooth() +
  scale_color_manual(values=scales)+
  ylim(0, 2) +
  theme_classic() 

ggplot(meta.beta, aes(x = pseudotime, y = norm.sXBP1, color=sex)) +
  #geom_point() +
  geom_smooth() +
  scale_color_manual(values=scales)+
  ylim(0, 2) +
  theme_classic()



