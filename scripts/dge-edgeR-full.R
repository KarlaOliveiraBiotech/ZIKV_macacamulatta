rm(list = ls(all.names = TRUE))


########## Load libraries ##########
library("tidyverse")
library("edgeR")
library("DESeq2")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("ashr")
library("DEGreport")
library("ggfortify")
library("GenomicRanges")
library("glue")
library("ggrepel")
library("EnhancedVolcano")
library("datawizard")
library("msigdbr")
library("clusterProfiler")
library("enrichplot")
library("RColorBrewer")
library("fgsea")
library("biomaRt")
########## Load libraries ##########

########## Import data as raw_counts and phenodata as lib_spec ##########
raw_counts <- readr::read_tsv(file = "workflow/results/feature_counts/raw_count_matrix_symbol.tsv", 
                              col_names = TRUE)

lib_spec <- readr::read_tsv(file = "resource/data/phenodata.tsv", 
                            col_names = TRUE)
########## Import data as raw_counts and phenodata as lib_spec ##########



########## Prep lib_spec ##########
spl <- strsplit(x = lib_spec$Sample, split = "_")
spl <- sapply(spl, "[[", 1)
spl <- as.numeric(spl)

lib_spec <- dplyr::mutate(lib_spec, spl)
lib_spec <- lib_spec[order(spl), ]
lib_spec <- lib_spec %>% dplyr::select(-c(5,6))

old_lib_spec <- lib_spec

lib_spec <- lib_spec %>% 
  tibble::column_to_rownames("Sample")
########## Prep lib_spec ##########


########## Decide factors ##########
lib_spec$Timepoint <- factor(lib_spec$Timepoint, 
                             levels = c("Dm14", "D1", "D3", 
                                        "D5", "D7", "D10", "D14"))
lib_spec$Animal <- factor(lib_spec$Animal, 
                          levels = c("A", "B", "C", "D"))
lib_spec$Class <- factor(lib_spec$Class, 
                         levels = c("baseline", "ZIKV"))
########## Decide factors ##########


########## Prep df ##########
df_unsort <- raw_counts
df <- df_unsort[!grepl("ENSMMU", df_unsort$Symbol), ]
df <- df[!duplicated(df$Symbol), ]
anyDuplicated(df)

df <- data.frame(df, row.names = 1)
cn2 <- base::gsub(pattern = "X(.*)", replacement = "\\1", x = colnames(df))
cn2 <- base::gsub(patter = "(00.)", replacement = "\\00-", x = cn2)
colnames(df) <- cn2

order <- old_lib_spec$Sample
df <- df[, order]


all(rownames(lib_spec) == colnames(df))
old_df <- df
########## Prep df ##########

########## Basic analysis ##########
dirOut <- paste0("workflow/results/edgeR/")
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }

###### 1. CPM of raw data ######
cpm <- cpm(df, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
min(cpm)

png(filename = paste0(dirOut, "/log2cpm_boxplot.png"))
boxplot(cpm, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(cpm),col="blue")
title("Boxplots of log2CPMs")
dev.off()


###### 2. PCA of raw data ######
pca <- prcomp(t(cpm))

pcaPlot <- autoplot(pca, data = lib_spec, colour = "Timepoint", shape = "Animal", size = 3) +
  ggtitle("PCA") + 
  theme(plot.title = element_text(hjust=0.5)) # Without clustering

pcaPlot2 <- autoplot(pca, data = lib_spec, colour = "Timepoint", shape = "Animal", size = 4, frame = TRUE) + 
  ggtitle("PCA - clustering") + 
  theme(plot.title = element_text(hjust=0.5)) # With clustering (frame)



ggsave(filename = paste0(dirOut, "/", "PCA_nocluster.png"), plot = pcaPlot)
ggsave(filename = paste0(dirOut, "/", "PCA_nocluster.pdf"), plot = pcaPlot)

ggsave(filename = paste0(dirOut, "/", "PCA_cluster.png"), plot = pcaPlot2)
ggsave(filename = paste0(dirOut, "/", "PCA_cluster.png"), plot = pcaPlot2)
########## Basic analysis ##########

########## DGE - edgeR ##########

###### 1. Factors Preparing ######
tp <- factor(lib_spec$Timepoint, levels = c("Dm14", "D1", "D3", "D5", 
                                            "D7", "D10", "D14"))
animal <- factor(lib_spec$Animal, levels = c("A", "B", "C", "D"))
design <- model.matrix(~ 0 + tp + animal)
###### Factors Preparing ######

###### 2. Matrix Preparing ######
y <- DGEList(counts = df, group = tp)

bfr_samples <- y$samples
keep <- filterByExpr(y = y, 
                     design = design, 
                     group = tp, 
                     min.count = 1, 
                     min.prop = 0.5) 

y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(object = y, method = "TMM")

y <- estimateDisp(y, design) 
###### 2. Matrix Preparing ######

###### 3. Contrasts ######
colnames(design) <- base::gsub(pattern = "tpD", 
                               replacement = "D", 
                               x = colnames(design))

colnames(design) <- base::gsub(pattern = "animal", 
                               replacement = "", 
                               x = colnames(design))

fit <- glmQLFit(y = y, design = design)


contrasts <- makeContrasts(D1_vs_Dm14 = D1 - Dm14,
                           D3_vs_Dm14 = D3 - Dm14,
                           D5_vs_Dm14 = D5 - Dm14,
                           D7_vs_Dm14 = D7 - Dm14,
                           D10_vs_Dm14 = D10 - Dm14,
                           D14_vs_Dm14 = D14 - Dm14,
                           levels = colnames(design))
###### 3. Contrasts ######

###### 4. DEG ######
for(i in colnames(contrasts)){
  qlf <- glmQLFTest(glmfit = fit, contrast = contrasts[,i])
  assign(paste0("qlf_", i), qlf)
  
  tTags <- as.data.frame(topTags(qlf, n = nrow(y), adjust.method = "BH"))
  
  # Direction is calculated based on: 
  # Pvalue and logFC at first
  tTags$DirectionPvalue <- ifelse(test = tTags$logFC > 1 & 
                                    tTags$PValue < 0.05, yes = "Up", no = 
                                    ifelse(test = tTags$logFC < -1 & 
                                             tTags$PValue < 0.05, 
                                           yes = "Down", no = "None"))
  # FDR and logFC second
  tTags$DirectionPadj <- ifelse(test = tTags$logFC > 1 & 
                                  tTags$FDR < 0.05, yes = "Up", no = 
                                  ifelse(test = tTags$logFC < -1 & 
                                           tTags$FDR < 0.05, 
                                         yes = "Down", no = "None"))
  
  
  
  assign(paste0("tTags_", i), tTags)
  
  tTags_up   <- tTags %>% dplyr::filter(DirectionPadj == "Up")
  tTags_down <- tTags %>% dplyr::filter(DirectionPadj == "Down")
  tTags_none <- tTags %>% dplyr::filter(DirectionPadj == "None")
  
  dirOut <- paste0("workflow/results/edgeR/", i)
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }
  
  # Output files - save    
  readr::write_tsv(x = data.frame(Ensembl_id = rownames(tTags), tTags), 
                   file = paste0(dirOut, ".tsv"))
  xlsx::write.xlsx(x = data.frame(Ensembl_id = rownames(tTags), tTags), 
                   file = paste0(dirOut, ".xlsx"), row.names = FALSE)
  
  # Histogram of Pvalues and FDR
  png(paste0(dirOut, "/Hist_PValue.png"))
  hist(tTags$PValue)
  dev.off()
  
  png(paste0(dirOut, "/Hist_FDR.png"))
  hist(tTags$FDR)
  dev.off()
  
  # Manual MA Plot
  p <- ggplot() + 
    geom_point(data = tTags_none, aes(x = logCPM, y = logFC, color = DirectionPadj), 
               size = 2, alpha = 0.5) + 
    geom_point(data = tTags_up, aes(x = logCPM, y = logFC, color = DirectionPadj),
               size = 2, alpha = 0.5) + 
    geom_point(data = tTags_down, aes(x = logCPM, y = logFC, color = DirectionPadj),
               size = 2, alpha = 0.5) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    theme_linedraw() + 
    theme(text = element_text(size = 9))
  
  if(nrow(tTags_down) == 0 & nrow(tTags_up) == 0) {
    p <- p + scale_color_manual(values = c("gray80"))
  } else if(nrow(tTags_down) > 0 & nrow(tTags_up) == 0) {
    p <- p + scale_color_manual(values = c("limegreen", "gray80"))
  } else if(nrow(tTags_down) == 0 & nrow(tTags_up) > 0) {
    p <- p + scale_color_manual(values = c("gray80", "red2"))
  } else {
    p <- p + 
      scale_color_manual(values = c("limegreen", "gray80", "red2"))
  }
  p <- p + ylim(c(-8, 8)) +
    annotate(geom = "text", x = 12, y = 5, 
             label = paste0("Up: ", nrow(tTags_up), " genes"), 
             color = "red2", size = 3) +
    annotate(geom = "text", x = 12, y = -5, 
             label = paste0("Down: ", nrow(tTags_down), " genes"), 
             color = "limegreen", size = 3)
  p
  
  ggsave(filename = paste0(dirOut, "/", i , "_MAplot.png"), 
         width = 5, height = 3)
  ggsave(filename = paste0(dirOut, "/", i, "_MAplot.pdf"), 
         width = 5, height = 3)
  
  
  # Manual Volcano Plot
  pVolcano <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_point(data = tTags_none, 
               aes(x = logFC, y = -log10(PValue), color = DirectionPadj), size = 2, 
               alpha = 0.5) +
    geom_point(data = tTags_up, 
               aes(x = logFC, y = -log10(PValue), color = DirectionPadj), size = 2, 
               alpha = 0.5) +
    geom_point(data = tTags_down, 
               aes(x = logFC, y = -log10(PValue), color = DirectionPadj), size = 2, 
               alpha = 0.5) +
    theme_linedraw() +
    theme(text = element_text(size = 8), panel.grid = element_blank())
  
  if(nrow(tTags_down) == 0 & nrow(tTags_up) == 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("gray80"))
  } else if(nrow(tTags_down) > 0 & nrow(tTags_up) == 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("limegreen", "gray80"))
  } else if(nrow(tTags_down) == 0 & nrow(tTags_up) > 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("gray80", "red2"))
  } else {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("limegreen", "gray80", "red2"))
  }
  pVolcano <- pVolcano + 
    xlim(c(-8, 8)) +
    annotate(geom = "text", x = 4,  y = 0, 
             label = paste0("Up: ", nrow(tTags_up), " genes"), 
             color = "red2", size = 3) +
    annotate(geom = "text", x = -4, y = 0, 
             label = paste0("Down: ", nrow(tTags_down), " genes"), 
             color = "limegreen", size = 3)
  pVolcano
  ggsave(filename = paste0(dirOut, "/", i, "_VolcanoPlot.png"), 
         width = 7, height = 7)
  ggsave(filename = paste0(dirOut, "/", i, "_VolcanoPlot.pdf"), 
         width = 7, height = 7)
  
  
  # Enhanced Volcano Plot
  pEnhanced <- EnhancedVolcano(tTags, lab = rownames(tTags), 
                               x = "logFC", y = "PValue", pCutoff = 0.05, 
                               FCcutoff = 1, pointSize = 3.0, labSize = 2, 
                               title = "Volcano Plot", 
                               subtitle = gsub(pattern = "(.*)(_vs_)(.*)", 
                                               replacement = "\\1 vs \\3", x = i))
  ggsave(filename = paste0(dirOut, "/", i, "_enhancedVolcano.png"))
  
  remove(tTags, tTags_down, tTags_up, tTags_none, qlf)
}