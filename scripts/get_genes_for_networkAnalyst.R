rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
########## Load libraries ##########

# OBS: After running this script, remove 'V1" from first line and checks if there are quotations. Remove them.


########## Path and list of folders ##########
folder_leadingEdges <- "workflow/results/DESeq2/fGSEA/leading_edges/reactome_level3/"
files_leadinEdges <- list.files(path = folder_leadingEdges, pattern = ".tsv")
########## Path and list of folders ##########


for(file in files_leadinEdges){
  df <- read_tsv(paste0(folder_leadingEdges, file))
  
  filename <- sub(pattern = "(.*)(_reactome_level3_immune.tsv)", 
                        replacement = "\\1", basename(file))

  df_split <- strsplit(as.character(df[,8]), "[[:blank:]]+")

  symbols <- sapply(df_split, unique) 
  symbols <- gsub(pattern = '"', replacement = '', x = symbols)  
  symbols <- gsub(pattern = '""', replacement = '', x = symbols) 
  symbols <- gsub(pattern = ',', replacement = '', x = symbols)  
  symbols <- gsub(pattern = ')', replacement = '', x = symbols)  
  symbols <- gsub(pattern = 'c(', replacement = '', x = symbols, fixed = TRUE)
  
  symbols_df <- as.data.table(symbols)
  
  dirOut <- "workflow/results/DESeq2/fGSEA/leading_edges/reactome_level3/only_genes/"
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
  
  data.table::fwrite(x = symbols_df, 
                     file = paste0(dirOut, filename, "_genes.tsv"), 
                     sep = "\t")

}
