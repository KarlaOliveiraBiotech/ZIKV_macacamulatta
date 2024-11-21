rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
library('gplots')
########## Load libraries ##########

########## Prepares the comparison dataframe ##########
# Reads reactome_level3_with_level1 file
db_reactome_level3_1 <- readr::read_tsv(file = "../../resource/ref/curated/ReactomePathwaysLevel3WithLevel1.tsv")

# Subsets immune system pathways
immune_system <- db_reactome_level3_1 %>% 
  dplyr::filter(Level1 == "Immune_System") %>% 
  dplyr::select(Pathway) %>% 
  base::unlist(use.names = TRUE)
########## File for comparison (immune system) ##########

########## Path and list of folders ##########
folder_fgsea <- "../results/DESeq2/fGSEA/db_reactome_level3/"


files_fgsea <- list.files(path = folder_fgsea, pattern = "_fGSEA.tsv")
########## Path and list of folders ##########

########## Creates files for enriched immune system pathway ##########
for(file_enrich in files_fgsea){
  df_fgsea <- read_tsv(paste0(folder_fgsea, file_enrich))
  
  filename_fgsea <- sub(pattern = "(.*)(_fGSEA).tsv", 
                        replacement = "\\1", basename(file_enrich))
  
  # Selects results from immune system pathways only
  df_immune <- df_fgsea %>% 
    dplyr::filter(pathway %in% immune_system)
  
  # Save results
  dirOut <- "../results/DESeq2/fGSEA/leading_edges/reactome_level3/"
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }

  data.table::fwrite(x = df_immune,
                     file = paste0(dirOut, filename_fgsea, "_reactome_level3_immune.tsv"),
                     sep = "\t")
  
  # Subsets pathway and NES from df
  df_immune_sub <- df_immune %>%
    dplyr::select(pathway, NES)
  
  # Assigns subset to different timeepoints, accordingly
  assign(filename_fgsea, df_immune_sub)
}
########## Creates files for enriched immune system pathway ##########

########## Plot ##########
# Creates a list of dataframes
list_df <- lapply(ls(pattern = "_vs_Dm14"), get)

# Makes a dataframe of NES values and immune system pathway
df_all <- list_df %>% 
  purrr::reduce(inner_join, by='pathway') %>% 
  tibble::column_to_rownames('pathway')

# Names samples
names(df_all) <- c("D1", "D3", "D5", "D7", "D10", "D14")

# Creates ans saves pheatmap plot
png(filename = paste0("../results/DESeq2/fGSEA/leading_edges/reactome_level3/", "Pheatmap_Enriched_ImmuneSystem.png"),
    width = 1400, height = 1300)
pheatmap(mat = df_all, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         main = "NES values from Enriched Immune System Pathway",
         legend = TRUE, 
         fontsize = 18)
dev.off()
########## Plot ##########








