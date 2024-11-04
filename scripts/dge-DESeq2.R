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
cn2 <- base::gsub(pattern = "(00.)", replacement = "\\00-", x = cn2)
colnames(df) <- cn2

order <- old_lib_spec$Sample
df <- df[, order]


all(rownames(lib_spec) == colnames(df))
old_df <- df
########## Prep df ##########



########## DESeq2 - prep ##########
dds <- DESeqDataSetFromMatrix(countData = df, 
                              colData = lib_spec, 
                              design = ~ Animal + Timepoint) 

minimumGeneExpression <- 1 # Decides the minimum expression to the gene be considered

keep <- rowSums(counts(dds) >= 3) >= minimumGeneExpression 
# Filters based on minimun expression

dds <- dds[keep, ]
dds <- DESeq(dds)

resultsNames(dds)
########## DESeq2 - prep ##########



########## DESeq2 - DE ##########
for(Timepoint in lib_spec %>% 
    dplyr::filter(Timepoint != "Dm14") %>%
    dplyr::select(Timepoint) %>%
    unlist(use.names = F) %>% sort() %>%
    unique()) {
  timepoint <- Timepoint
  control <- "Dm14"
  tp_contrasts <- c("Timepoint", timepoint, control)
  res <- results(dds, contrast = tp_contrasts, alpha = 0.05, 
                 pAdjustMethod = "BH")

  
  resultsNames(dds)
  
  resLFC <- as.data.frame(lfcShrink(dds, 
                                    contrast = c("Timepoint", timepoint, control), 
                                    type = "ashr"))
  
  
  resLFC$DirectionPvalue <- ifelse(test = resLFC$log2FoldChange > 1 &
                                     resLFC$pvalue < 0.05 &
                                     resLFC$baseMean > 0,
                                   yes = "Up",
                                   no = ifelse(test = resLFC$log2FoldChange < -1 &
                                                 resLFC$pvalue < 0.05 &
                                                 resLFC$baseMean > 0,
                                               yes = "Down",
                                               no = "None"))
  
  resLFC$DirectionPadj <- ifelse(test = resLFC$log2FoldChange > 1 &
                                   resLFC$padj < 0.05 &
                                   resLFC$baseMean > 0,
                                 yes = "Up",
                                 no = ifelse(test = resLFC$log2FoldChange < -1 &
                                               resLFC$padj < 0.05 &
                                               resLFC$baseMean > 0,
                                             yes = "Down",
                                             no = "None"))
  
  resLFC$DirectionPadj[is.na(resLFC$DirectionPadj)] <- "None"
  #resLFC$DirectionPadj <- replace(resLFC$Direction, is.na(resLFC$Direction), "None")
  
  assign(paste0("resLFC_", timepoint, "_vs_", control), resLFC)
  
  resLFC_up   <- resLFC %>% dplyr::filter(DirectionPadj == "Up")
  resLFC_down <- resLFC %>% dplyr::filter(DirectionPadj == "Down")
  resLFC_none <- resLFC %>% dplyr::filter(DirectionPadj == "None")
  
  dirOut <- paste0("workflow/results/DESeq2/", timepoint, "_vs_", control)
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
  
  # Output files - save    
  readr::write_tsv(x = data.frame(Ensembl_id = rownames(resLFC), resLFC), 
                   file = paste0(dirOut, ".tsv"))
  xlsx::write.xlsx(x = data.frame(Ensembl_id = rownames(resLFC), resLFC), 
                   file = paste0(dirOut, ".xlsx"), row.names = FALSE)
  
  # Histogram of Pvalues and Padj
  png(paste0(dirOut, "/Hist_Pvalue.png"))
  hist(resLFC$pvalue)
  dev.off()
  
  png(paste0(dirOut, "/Hist_FDR.png"))
  hist(resLFC$padj)
  dev.off()   
  
  # Manual MA Plot
  p <- ggplot() + 
    geom_point(data = resLFC_none, 
               aes(x = log10(baseMean), y = log2FoldChange, color = DirectionPadj),
               size = 2, alpha = 0.5) + 
    geom_point(data = resLFC_up, 
               aes(x = log10(baseMean), y = log2FoldChange, color = DirectionPadj), 
               size = 2, alpha = 0.5) + 
    geom_point(data = resLFC_down, 
               aes(x = log10(baseMean), y = log2FoldChange, color = DirectionPadj), 
               size = 2, alpha = 0.5) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_linedraw() + 
    theme(text = element_text(size = 9))
  
  if(nrow(resLFC_down) == 0 & nrow(resLFC_up) == 0) {
    p <- p + scale_color_manual(values = c("gray80"))
  } else if(nrow(resLFC_down) > 0 & nrow(resLFC_up) == 0) {
    p <- p + scale_color_manual(values = c("limegreen", "gray80"))
  } else if(nrow(resLFC_down) == 0 & nrow(resLFC_up) > 0) {
    p <- p + scale_color_manual(values = c("gray80", "red2"))
  } else {
    p <- p + 
      scale_color_manual(values = c("limegreen", "gray80", "red2"))
  }
  
  p <- p + ylim(c(-8, 8)) +
    annotate(geom = "text", x = 3, y = 5, 
             label = paste0("Up: ", nrow(resLFC_up), " genes"), 
             color = "red2", size = 3) +
    annotate(geom = "text", x = 3, y = -5, 
             label = paste0("Down: ", nrow(resLFC_down), " genes"), 
             color = "limegreen", size = 3)
  p
  
  ggsave(filename = paste0(dirOut, "/", timepoint, "_vs_", control , "_MAplot.png"),
         width = 5, height = 3)
  ggsave(filename = paste0(dirOut, "/", timepoint, "_vs_", control, "_MAplot.pdf"),
         width = 5, height = 3)
  
  # Manual Volcano Plot
  pVolcano <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_point(data = resLFC_none, 
               aes(x = log2FoldChange, y = -log10(pvalue), color = DirectionPadj), 
               size = 2, alpha = 0.5) +
    geom_point(data = resLFC_up,
               aes(x = log2FoldChange, y = -log10(pvalue), color = DirectionPadj), 
               size = 2, alpha = 0.5) +
    geom_point(data = resLFC_down, 
               aes(x = log2FoldChange, y = -log10(pvalue), color = DirectionPadj), 
               size = 2, alpha = 0.5) +
    theme_linedraw() +
    theme(text = element_text(size = 8), panel.grid = element_blank())
  
  if(nrow(resLFC_down) == 0 & nrow(resLFC_up) == 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("gray80"))
  } else if(nrow(resLFC_down) > 0 & nrow(resLFC_up) == 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("limegreen", "gray80"))
  } else if(nrow(resLFC_down) == 0 & nrow(resLFC_up) > 0) {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("gray80", "red2"))
  } else {
    pVolcano <- pVolcano + 
      scale_color_manual(values = c("limegreen", "gray80", "red2"))
  }
  pVolcano <- pVolcano + xlim(c(-8, 8)) +
    annotate(geom = "text", x = 4,  y = 0, 
             label = paste0("Up: ", nrow(resLFC_up), " genes"), 
             color = "red2", size = 3) +
    annotate(geom = "text", x = -4, y = 0, 
             label = paste0("Down: ", nrow(resLFC_down), " genes"), 
             color = "limegreen", size = 3)
  pVolcano
  
  ggsave(filename = paste0(dirOut, "/", timepoint, "_vs_", control, "_VolcanoPlot.png"),
         width = 7, height = 7)
  ggsave(filename = paste0(dirOut, "/", timepoint, "_vs_", control, "_VolcanoPlot.pdf"),
         width = 7, height = 7)
  
  # Enhanced Volcano Plot
  pEnhanced <- EnhancedVolcano(resLFC, lab = rownames(resLFC), 
                               x = "log2FoldChange", y = "pvalue", 
                               pCutoff = 0.05, FCcutoff = 1, 
                               pointSize = 3.0, labSize = 2, 
                               xlim = c(-6,6), 
                               title = "Volcano Plot", 
                               subtitle = paste0(timepoint, " vs ", control))
  ggsave(filename = paste0(dirOut, "/", timepoint, "_vs_", control, "_enhancedVolcano.png"))   
  
  remove(resLFC, resLFC_down, resLFC_up, resLFC_none, res)
}
########## DESeq2 - DE ##########



########## Cluster samples ##########
select <- order(rowMeans(counts(dds, normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("Animal","Timepoint")])

annoCol <- list(Timepoint = c(Dm14 = "#005824", 
                              D1 = "#238B45", 
                              D3 = "#41AE76",  
                              D5 = "#66C2A4", 
                              D7 = "#99D8C9", 
                              D10 = "#CCECE6", 
                              D14 = "#EDF8FB"),
                Animal = c(A = "#FEE391", 
                           B = "#FEC44F", 
                           C = "#FE9929", 
                           D = "#EC7014"))

png(paste0("workflow/results/DESeq2/", "Pheatmap_clustCols.png"))
pheatmap(assay(dds)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=TRUE, 
         annotation_col=df, 
         annotation_colors = annoCol)
dev.off()

png(paste0("workflow/results/DESeq2/", "Pheatmap.png"))
pheatmap(assay(dds)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=df,
         annotation_colors = annoCol)
dev.off()
########## Cluster samples ##########

