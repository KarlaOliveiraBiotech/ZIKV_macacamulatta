rm(list = ls(all.names = TRUE))


########## Load libraries ##########
library("tidyverse")
library("edgeR")
library("dplyr")
library("readr")
library("ggplot2")
library("limma")
library("reshape2")
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
dirOut <- paste0("workflow/results/DESeq2/")
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }

###### 1. CPM of raw data ######
cpm <- cpm(df, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
min(cpm)

png(filename = paste0(dirOut, "/log2cpm_boxplot.png"))
boxplot(cpm, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(cpm),col="blue")
title("Boxplots of log2CPMs")
dev.off()


cpm_normalized <- normalizeQuantiles(cpm)
cpm_melted <- melt(cpm_normalized)


p <- ggplot(cpm_melted, aes(x = Var2, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplot of Normalized CPM Data", x = "Samples", y = "Normalized Expression Level") +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p

ggsave(filename = paste0(dirOut, "/", "log2cpm_normalized.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "log2cpm_normalized.pdf"),
       width = 5, height = 5)

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