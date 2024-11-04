rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
library("CEMiTool")
########## Load libraries ##########


D1 <- readr::read_tsv("workflow/results/DESeq2/D1_vs_Dm14.tsv")
D3 <- readr::read_tsv("workflow/results/DESeq2/D3_vs_Dm14.tsv")
D5 <- readr::read_tsv("workflow/results/DESeq2/D5_vs_Dm14.tsv")
D7 <- readr::read_tsv("workflow/results/DESeq2/D7_vs_Dm14.tsv")
D10 <- readr::read_tsv("workflow/results/DESeq2/D10_vs_Dm14.tsv")
D14 <- readr::read_tsv("workflow/results/DESeq2/D14_vs_Dm14.tsv")

D1_sub <- D1 %>% dplyr::select(Ensembl_id, log2FoldChange)
D3_sub <- D3 %>% dplyr::select(Ensembl_id, log2FoldChange)
D5_sub <- D5 %>% dplyr::select(Ensembl_id, log2FoldChange)
D7_sub <- D7 %>% dplyr::select(Ensembl_id, log2FoldChange)
D10_sub <- D10 %>% dplyr::select(Ensembl_id, log2FoldChange)
D14_sub <- D14 %>% dplyr::select(Ensembl_id, log2FoldChange)


df_list <- list(D1_sub, D3_sub, D5_sub, D7_sub, D10_sub, D14_sub)
df_all <- df_list %>% reduce(inner_join, by='Ensembl_id')

names(df_all) <- c("Symbol", "D1", "D3", "D5", "D7", "D10", "D14")

dirOut <- "workflow/results/DESeq2/subsets/"
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }

data.table::fwrite(x = df_all,
                   file = paste0(dirOut, "all_data_log2FC.tsv"),
                   sep = "\t")

