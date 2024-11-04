rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
########## Load libraries ##########


########## Listing pathways ##########
db_btm <- "db_btm"
db_hallmarks <- "db_hallmarks"
db_kegg_no_disease <- "db_kegg_no_disease"
db_reactome_all <- "db_reactome_all"
db_reactome_level3 <- "db_reactome_level3"
db_reactome_level3_1 <- "db_reactome_level3_1"
########## Listing pathways ##########

# Choose pathway (db)
chosen_db <- db_reactome_level3


########## Path and list of folders ##########
folder_deseq <- "workflow/results/DESeq2/"
files_deseq <- list.files(path = folder_deseq, pattern = ".tsv")

folder_fgsea <- paste0("workflow/results/DESeq2/fGSEA/", chosen_db, "/")
files_fgsea <- list.files(path = folder_fgsea, pattern = "_fGSEA.tsv")
########## Path and list of folders ##########


########## Check up-regulated genes in enriched files ##########
for(file_deg in files_deseq){
  df_deseq <- readr::read_tsv(paste0(folder_deseq, file_deg))
  names(df_deseq)[1] <- "Symbol"
  
  filename_deseq <- sub(pattern = "(.*).tsv", 
                        replacement = "\\1", 
                        basename(file_deg))  
  
  for(file_enrich in files_fgsea){
    df_fgsea <- read_tsv(paste0(folder_fgsea, file_enrich))
    
    filename_fgsea <- sub(pattern = "(.*)(_fGSEA).tsv", 
                          replacement = "\\1", basename(file_enrich))
    
    
    if(filename_deseq == filename_fgsea){
      #message("This is DEG: ", filename_deseq)
      #message("This is GSEA: ", filename_fgsea)
      #print("They are equal!")
      
      # Get leading edge genes from 10 pathways 
      fgsea_up10 <- strsplit(as.character(df_fgsea[1:10,8]), "[[:blank:]]+")
      fgsea_up_unique <- sapply(fgsea_up10, unique)
      
      # Get DEG (up and down) for comparison
      deseq_up <- df_deseq %>% 
        filter(DirectionPadj == "Up")
      
      deseq_down <- df_deseq %>% 
        filter(DirectionPadj == "Down")
      
      # Compare
      up_checked <- deseq_up[deseq_up$Symbol %in% fgsea_up_unique, ]
      down_checked <- deseq_down[deseq_down$Symbol %in% fgsea_up_unique, ]
      
      
      dirOut <- paste0("workflow/results/DESeq2/fGSEA/leading_edges/", chosen_db, "/")
      if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
       
      # Saves leading edges output files 
      data.table::fwrite(x = up_checked, 
                          file = paste0(dirOut, filename_deseq, "_upchecked.tsv"), 
                          sep = "\t")
      
      data.table::fwrite(x = down_checked, 
                         file = paste0(dirOut, filename_deseq, "_downchecked.tsv"), 
                         sep = "\t")
    }
  }
}
########## Check up-regulated genes in enriched files ##########