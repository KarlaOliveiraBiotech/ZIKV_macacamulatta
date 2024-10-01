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
library("ggupset")
library("plyr")
########## Load libraries ##########

dirIn <- "workflow/results/DESeq2/"
files_dge <- list.files(path = dirIn, pattern = ".tsv")

file_gmt <- "resource/ref/msigdb_v2024.1.Hs_GMTs/c5.all.v2024.1.Hs.symbols.gmt"
rds_path <- "resource/ref/msigdb_v2024.1.Hs_GMTs/" 

#files_dge <- list.files(path = dirIn, pattern = ".tsv")  

for(file in files_dge){
  df <- read_tsv(file = paste0(dirIn, file))
  
  names(df)[1] <- "gene_symbol" 

  genes_in_data <- df$gene_symbol

  bg_full <- read.gmt(file_gmt)
  bg_subset <- bg_full[bg_full$gene %in% genes_in_data, ]

  rds <- "bg_go.RDS"
  saveRDS(bg_subset, paste0(rds_path, rds))

  df_diffexp <- df[df$DirectionPadj != "None", ]
  df_diffexp <- df_diffexp %>%
    dplyr::select(gene_symbol, pvalue, padj, log2FoldChange, DirectionPadj)


  deg_results_list <- split(df_diffexp, df_diffexp$DirectionPadj)


  background_genes <- 'GO'
  bg_genes <- readRDS("resource/ref/msigdb_v2024.1.Hs_GMTs/bg_go.RDS")
  padj_cutoff <- 0.05
  genecount_cutoff <- 5

  dirOut <- paste0("workflow/results/DESeq2/clusterProfiler/")
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }

  filename_cut <- base::gsub(pattern = ".tsv", replacement = "", x = file)
  filename <- paste0(dirOut, filename_cut, "_", background_genes)

  res <- lapply(names(deg_results_list),
                function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                     TERM2GENE = bg_genes))

  names(res) <- names(deg_results_list)

  res_df <- ldply(res, data.frame)
  res_df <- res_df %>%
    dplyr::mutate(minuslog10padj = -log10(p.adjust))

  target <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff])
  res_df <- res_df[res_df$ID %in% target, ]


  print('Saving clusterprofiler results')
  write.csv(x = res_df, 
            file = paste0(filename, '_resclusterp.csv'),
            row.names = FALSE)

  enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)

  
  
  #barplot(enrichres, showCategory = 15)
  #dotplot(enrichres, showCategory = 15) + ggtitle("D1 vs Dm14")

  #enrichres2 <- pairwise_termsim(enrichres)
  #treeplot(enrichres2)
}
