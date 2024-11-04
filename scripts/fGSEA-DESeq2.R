rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("msigdbr")
library("RColorBrewer")
library("fgsea")
library("biomaRt")
library("data.table")
########## Load libraries ##########

########## Path and list of folders ##########
folder <- "workflow/results/DESeq2/"
files <- list.files(path = folder, pattern = ".tsv")
########## Path and list of folders ##########

########## Set functions ##########
### Function: Gets file names
pathway_function <- function(pathway) {
  deparse(substitute(pathway))
}
########## Set functions ##########


########## Possible pathways ##########
db_hallmarks      <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/h.all.v2024.1.Hs.symbols.gmt")
db_go_all         <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/c5.all.v2024.1.Hs.symbols.gmt")
db_reactome_all   <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
db_go_bp          <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/c5.go.bp.v2024.1.Hs.symbols.gmt")        

db_reactome_level3    <- gmtPathways("resource/ref/curated/ReactomePathwaysLevel3.gmt")
db_kegg_no_disease    <- gmtPathways("resource/ref/curated/KEGG_pathways_NO_diseases.gmt")
db_btm                <- gmtPathways("resource/ref/curated/BTM_for_GSEA_20131008.gmt")
########## Possible pathways ##########


########## Choose pathway ##########
chosen_pathway <- db_btm # Choose pathway to be tested according to previous names
pathway_name <- pathway_function(db_btm) # Choose pathway to be tested according to previous names
########## Choose pathway ##########

########## Analysis ##########
for(file in files){
  # Takes file's name based on folder and list above
  filename <- sub(pattern = "(.*).tsv", replacement = "\\1", basename(file)) 
  df <- read_tsv(paste0(folder, file), col_names = TRUE)
  
  res <- df %>% 
    dplyr::select(Ensembl_id, log2FoldChange, pvalue, padj) 
  names(res)[1] <- "SYMBOL"
  
  # Removes duplication of gene symbol in df 
  res_uniq <- res[!duplicated(res[, "SYMBOL"]), ] 
  
  # Creates a rank of pi-pvalue
  ranks <- data.frame(SYMBOL = res_uniq$SYMBOL, 
                      Log2FC = res_uniq$log2FoldChange * -log10(res_uniq$pvalue))
  
  # Orders rank and remove NAs
  ranks <- ranks[order(-ranks$Log2FC), ]
  ranks <- ranks[!is.na(ranks$Log2FC), ]
  
  # Transform rank in list
  ranks2 <- deframe(ranks)
  
  # Performs fgsea analysis
  fgseaRes <- fgsea(pathways = chosen_pathway,
                    stats = ranks2, 
                    minSize = 15, 
                    maxSize = 500)
  
  # Creates folder for general outputs
  dirOut <- paste0("workflow/results/DESeq2/fGSEA/", pathway_name, "/")
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
  
  # Saves fGSEA output files 
  data.table::fwrite(x = fgseaRes[order(-NES),], 
                     file = paste0(dirOut, filename, "_fGSEA.tsv"), 
                     sep = "\t", 
                     sep2 = c("", " ", ""))
  
  # readr::write_tsv(x = data.frame(fgseaRes[order(-NES), ]),
  #                  file = paste0(dirOut, filename, "_fGSEA.tsv"), 
  #                  sep = "\t", sep2 = c("", " ", ""))
  
  xlsx::write.xlsx(x = data.frame(fgseaRes[order(-NES), ][, 1:7]),
                   file = paste0(dirOut, filename, "_fGSEA.xlsx"), 
                   row.names = FALSE)
  

  
  # Does over-representation pathways
  fg <- names(head(ranks2[order(ranks2, decreasing=TRUE)], 200))
  bg <- names(ranks2)
  foraRes <- fora(genes = fg, universe = bg, pathways = chosen_pathway)
  
  # Saves output files from over-representation analysis
  data.table::fwrite(x = foraRes, 
                     file = paste0(dirOut, filename, "_foraRes.tsv"), 
                     sep = "\t", 
                     sep2 = c("", " ", ""))
  
  
  # readr::write_tsv(x = data.frame(foraRes),
  #                  file = paste0(dirOut, filename, "_foraRes.tsv"))
  
  xlsx::write.xlsx(x = data.frame(fgseaRes[, 1:7]),
                   file = paste0(dirOut, filename, "_foraRes.xlsx"), 
                   row.names = FALSE)
  
  # Creates folder for plots
  dirPlot <- paste0("workflow/results/DESeq2/fGSEA/", pathway_name, "/Plots/")
  if(!file.exists(dirPlot)) { dir.create(path = dirPlot, recursive = TRUE) }
  
  # Plots - for 5 first pathways 
  ## Plot 1
  top1 <- plotEnrichment(chosen_pathway[[fgseaRes[order(-NES)]$pathway[1]]], 
                         ranks2) + 
    labs(title = fgseaRes[order(-NES)]$pathway[1])
  ggsave(filename = paste0(dirPlot, filename, "_top1_EnrichPlot.png"), 
         plot = top1)
  
  
  ## Plot 2
  top2 <- plotEnrichment(chosen_pathway[[fgseaRes[order(-NES)]$pathway[2]]], 
                 ranks2) + 
    labs(title = fgseaRes[order(-NES)]$pathway[2])
  ggsave(filename = paste0(dirPlot, filename, "_top2_EnrichPlot.png"), 
         plot = top2)
  
  ## Plot 3
  top3 <- plotEnrichment(chosen_pathway[[fgseaRes[order(-NES)]$pathway[3]]], 
                 ranks2) + 
    labs(title = fgseaRes[order(-NES)]$pathway[3])
  ggsave(filename = paste0(dirPlot, filename, "_top3_EnrichPlot.png"), 
         plot = top3)
  
  ## Plot 4
  top4 <- plotEnrichment(chosen_pathway[[fgseaRes[order(-NES)]$pathway[4]]], 
                 ranks2) + 
    labs(title = fgseaRes[order(-NES)]$pathway[4])
  ggsave(filename = paste0(dirPlot, filename, "_top4_EnrichPlot.png"), 
         plot = top4)
  
  ## Plot 5
  top5 <- plotEnrichment(chosen_pathway[[fgseaRes[order(-NES)]$pathway[5]]],
                 ranks2) +
    labs(title = fgseaRes[order(-NES)]$pathway[5])
  ggsave(filename = paste0(dirPlot, filename, "_top5_EnrichPlot.png"), 
         plot = top5)

  # Plot for multiple pathways (up and down)
  topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=10), pathway]
  topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
   
  plotMultiple <- plotGseaTable(pathways = chosen_pathway[topPathways], 
                                stats = ranks2, 
                                fgseaRes = fgseaRes, 
                                gseaParam=0.5, 
                                pathwayLabelStyle = list(size = 6, color = "black"), 
                                valueStyle = list(size = 6, color = "black"), 
                                headerLabelStyle = list(size = 8, color = "black"))
  ggsave(filename = paste0(dirPlot, filename, "_gseaTable.png"), plot = plotMultiple)

}
########## Analysis ##########

