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


########## Directory Output  ##########
dirOut <- paste0("workflow/results/DGE-comparison/")
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }
########## Directory Output  ##########


###### 1. Factors Preparing ######
tp <- factor(lib_spec$Timepoint, levels = c("Dm14", "D1", "D3", "D5", 
                                            "D7", "D10", "D14"))
animal <- factor(lib_spec$Animal, levels = c("A", "B", "C", "D"))

design_edger <- model.matrix(~ 0 + tp + animal) 
###### 1. Factors Preparing ######

###### 2. Matrix Preparing ######
y <- DGEList(counts = df, group = tp)
y <- normLibSizes(object = y, method = "TMM")
y <- estimateDisp(y, design_edger)
###### 2. Matrix Preparing ######

###### 3. Contrasts ######
colnames(design_edger) <- base::gsub(pattern = "tpD", 
                                     replacement = "D", 
                                     x = colnames(design_edger))

colnames(design_edger) <- base::gsub(pattern = "animal", 
                                     replacement = "", 
                                     x = colnames(design_edger))

fit <- glmQLFit(y = y, design = design_edger)


contrasts <- makeContrasts(D1_vs_Dm14 = D1 - Dm14,
                           D3_vs_Dm14 = D3 - Dm14,
                           D5_vs_Dm14 = D5 - Dm14,
                           D7_vs_Dm14 = D7 - Dm14,
                           D10_vs_Dm14 = D10 - Dm14,
                           D14_vs_Dm14 = D14 - Dm14,
                           levels = colnames(design_edger))
###### 3. Contrasts ######

###### 4. DEG ######
for(i in colnames(contrasts)){
  qlf <- glmQLFTest(glmfit = fit, contrast = contrasts[,i])
  assign(paste0("qlf_", i), qlf)
  
  tTags <- as.data.frame(topTags(qlf, n = nrow(y), adjust.method = "BH"))
  
  tTags$DirectionPvalue <- ifelse(test = tTags$logFC > 1 & 
                              tTags$PValue < 0.05, yes = "Up", no = 
                              ifelse(test = tTags$PValue < -1 & 
                                       tTags$FDR < 0.05, 
                                     yes = "Down", no = "None"))
  
  tTags$Direction <- ifelse(test = tTags$logFC > 1 & 
                                  tTags$FDR < 0.05, yes = "Up", no = 
                                  ifelse(test = tTags$logFC < -1 & 
                                           tTags$FDR < 0.05, 
                                         yes = "Down", no = "None"))
  
  
  
  assign(paste0("tTags_", i), tTags)
}
###### 4. DEG ######
########## DGE - edgeR ##########


########## DESeq2 - prep ##########
dds <- DESeqDataSetFromMatrix(countData = df, 
                              colData = lib_spec, 
                              design = ~ Animal + Timepoint) 
dds <- DESeq(dds)
resultsNames(dds)
########## DESeq2 - prep ##########


########## DGE - DESeq2 ##########
for(Timepoint in lib_spec %>% 
    dplyr::filter(Timepoint != "Dm14") %>%
    dplyr::select(Timepoint) %>%
    unlist(use.names = F) %>% sort() %>%
    unique()) {
  timepoint <- Timepoint
  control <- "Dm14"
  tp_contrasts <- c("Timepoint", timepoint, control)
  
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
  
  resLFC$Direction <- ifelse(test = resLFC$log2FoldChange > 1 &
                                   resLFC$padj < 0.05 &
                                   resLFC$baseMean > 0,
                                 yes = "Up",
                                 no = ifelse(test = resLFC$log2FoldChange < -1 &
                                               resLFC$padj < 0.05 &
                                               resLFC$baseMean > 0,
                                             yes = "Down",
                                             no = "None"))
  
  resLFC$Direction <- replace(resLFC$Direction, is.na(resLFC$Direction), "None")
  
  assign(paste0("resLFC_", timepoint, "_vs_", control), resLFC)
}  
########## DGE - DESeq2 ##########


## D1
D1_edgeR <- tTags_D1_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D1_edgeR <- rownames_to_column(D1_edgeR)

D1_DESEq2 <- resLFC_D1_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D1_DESEq2 <- rownames_to_column(D1_DESEq2)

D1_inner <- inner_join(x = D1_edgeR, y = D1_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D1_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D1 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))

p 
ggsave(filename = paste0(dirOut, "/", "D1", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D1", "_comparison.pdf"),
       width = 5, height = 5)

## D3
D3_edgeR <- tTags_D3_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D3_edgeR <- rownames_to_column(D3_edgeR)

D3_DESEq2 <- resLFC_D3_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D3_DESEq2 <- rownames_to_column(D3_DESEq2)

D3_inner <- inner_join(x = D3_edgeR, y = D3_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D3_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D3 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))
p

ggsave(filename = paste0(dirOut, "/", "D3", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D3", "_comparison.pdf"),
       width = 5, height = 5)


## D5
D5_edgeR <- tTags_D5_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D5_edgeR <- rownames_to_column(D5_edgeR)

D5_DESEq2 <- resLFC_D5_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D5_DESEq2 <- rownames_to_column(D5_DESEq2)

D5_inner <- inner_join(x = D5_edgeR, y = D5_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D5_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D5 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))
p 
ggsave(filename = paste0(dirOut, "/", "D5", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D5", "_comparison.pdf"),
       width = 5, height = 5)


## D7
D7_edgeR <- tTags_D7_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D7_edgeR <- rownames_to_column(D7_edgeR)

D7_DESEq2 <- resLFC_D7_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D7_DESEq2 <- rownames_to_column(D7_DESEq2)

D7_inner <- inner_join(x = D7_edgeR, y = D7_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D7_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D7 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))
p 
ggsave(filename = paste0(dirOut, "/", "D7", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D7", "_comparison.pdf"),
       width = 5, height = 5)


## D10
D10_edgeR <- tTags_D10_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D10_edgeR <- rownames_to_column(D10_edgeR)

D10_DESEq2 <- resLFC_D10_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D10_DESEq2 <- rownames_to_column(D10_DESEq2)

D10_inner <- inner_join(x = D10_edgeR, y = D10_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D10_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D10 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))
p 
ggsave(filename = paste0(dirOut, "/", "D10", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D10", "_comparison.pdf"),
       width = 5, height = 5)


## D14
D14_edgeR <- tTags_D14_vs_Dm14 %>% dplyr::select(logFC, PValue, FDR, Direction)
D14_edgeR <- rownames_to_column(D14_edgeR)

D14_DESEq2 <- resLFC_D14_vs_Dm14 %>% dplyr::select(log2FoldChange, pvalue, padj, Direction)
D14_DESEq2 <- rownames_to_column(D14_DESEq2)

D14_inner <- inner_join(x = D14_edgeR, y = D14_DESEq2, by = "rowname")

p <- ggplot() + 
  geom_point(data = D14_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("D14 vs Dm14") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-10,10))
p 
ggsave(filename = paste0(dirOut, "/", "D14", "_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "D14", "_comparison.pdf"),
       width = 5, height = 5)





cols <- c("down", "up", "none")

## D1
D1_sum_down_edgeR <- sum(D1_edgeR$Direction == "Down")
D1_sum_up_edgeR <- sum(D1_edgeR$Direction == "Up")
D1_sum_none_edgeR <- sum(D1_edgeR$Direction == "None")

D1_sum_down_DESeq <- sum(D1_DESEq2$Direction == "Down")
D1_sum_up_DESeq <- sum(D1_DESEq2$Direction == "Up")
D1_sum_none_DESeq <- sum(D1_DESEq2$Direction == "None")

edger <- c(D1_sum_down_edgeR, D1_sum_up_edgeR, D1_sum_none_edgeR)
deseq <- c(D1_sum_down_DESeq, D1_sum_up_DESeq, D1_sum_none_DESeq)

D1 <- rbind(edger, deseq)
colnames(D1) <- c("down_D1", "up_D1", "none_D1")



## D3
D3_sum_down_edgeR <- sum(D3_edgeR$Direction == "Down")
D3_sum_up_edgeR <- sum(D3_edgeR$Direction == "Up")
D3_sum_none_edgeR <- sum(D3_edgeR$Direction == "None")

D3_sum_down_DESeq <- sum(D3_DESEq2$Direction == "Down")
D3_sum_up_DESeq <- sum(D3_DESEq2$Direction == "Up")
D3_sum_none_DESeq <- sum(D3_DESEq2$Direction == "None")

edger <- c(D3_sum_down_edgeR, D3_sum_up_edgeR, D3_sum_none_edgeR)
deseq <- c(D3_sum_down_DESeq, D3_sum_up_DESeq, D3_sum_none_DESeq)

D3 <- rbind(edger, deseq)
colnames(D3) <- c("down_D3", "up_D3", "none_D3")



## D5
D5_sum_down_edgeR <- sum(D5_edgeR$Direction == "Down")
D5_sum_up_edgeR <- sum(D5_edgeR$Direction == "Up")
D5_sum_none_edgeR <- sum(D5_edgeR$Direction == "None")

D5_sum_down_DESeq <- sum(D5_DESEq2$Direction == "Down")
D5_sum_up_DESeq <- sum(D5_DESEq2$Direction == "Up")
D5_sum_none_DESeq <- sum(D5_DESEq2$Direction == "None")

edger <- c(D5_sum_down_edgeR, D5_sum_up_edgeR, D5_sum_none_edgeR)
deseq <- c(D5_sum_down_DESeq, D5_sum_up_DESeq, D5_sum_none_DESeq)

D5 <- rbind(edger, deseq)
colnames(D5) <- c("down_D5", "up_D5", "none_D5")



## D7
D7_sum_down_edgeR <- sum(D7_edgeR$Direction == "Down")
D7_sum_up_edgeR <- sum(D7_edgeR$Direction == "Up")
D7_sum_none_edgeR <- sum(D7_edgeR$Direction == "None")

D7_sum_down_DESeq <- sum(D7_DESEq2$Direction == "Down")
D7_sum_up_DESeq <- sum(D7_DESEq2$Direction == "Up")
D7_sum_none_DESeq <- sum(D7_DESEq2$Direction == "None")

edger <- c(D7_sum_down_edgeR, D7_sum_up_edgeR, D7_sum_none_edgeR)
deseq <- c(D7_sum_down_DESeq, D7_sum_up_DESeq, D7_sum_none_DESeq)

D7 <- rbind(edger, deseq)
colnames(D7) <- c("down_D7", "up_D7", "none_D7")



## D10
D10_sum_down_edgeR <- sum(D10_edgeR$Direction == "Down")
D10_sum_up_edgeR <- sum(D10_edgeR$Direction == "Up")
D10_sum_none_edgeR <- sum(D10_edgeR$Direction == "None")

D10_sum_down_DESeq <- sum(D10_DESEq2$Direction == "Down")
D10_sum_up_DESeq <- sum(D10_DESEq2$Direction == "Up")
D10_sum_none_DESeq <- sum(D10_DESEq2$Direction == "None")

edger <- c(D10_sum_down_edgeR, D10_sum_up_edgeR, D10_sum_none_edgeR)
deseq <- c(D10_sum_down_DESeq, D10_sum_up_DESeq, D10_sum_none_DESeq)

D10 <- rbind(edger, deseq)
colnames(D10) <- c("down_D10", "up_D10", "none_D10")



## D14
D14_sum_down_edgeR <- sum(D14_edgeR$Direction == "Down")
D14_sum_up_edgeR <- sum(D14_edgeR$Direction == "Up")
D14_sum_none_edgeR <- sum(D14_edgeR$Direction == "None")

D14_sum_down_DESeq <- sum(D14_DESEq2$Direction == "Down")
D14_sum_up_DESeq <- sum(D14_DESEq2$Direction == "Up")
D14_sum_none_DESeq <- sum(D14_DESEq2$Direction == "None")

edger <- c(D14_sum_down_edgeR, D14_sum_up_edgeR, D14_sum_none_edgeR)
deseq <- c(D14_sum_down_DESeq, D14_sum_up_DESeq, D14_sum_none_DESeq)

D14 <- rbind(edger, deseq)
colnames(D14) <- c("down_D14", "up_D14", "none_D14")


dge_comparison <- cbind(D1, D3, D5, D7, D10, D14)
readr::write_tsv(x = as.data.frame(dge_comparison), 
                 file = "workflow/results/DGE-comparison/dge_comparison.tsv")