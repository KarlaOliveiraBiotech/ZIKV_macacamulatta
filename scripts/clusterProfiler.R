rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("msigdbr")
library("clusterProfiler")
library("enrichplot")
library("plyr")
library("ggrepel")
########## Load libraries ##########


########## Set functions ##########
### Function: Barplot
enrichres_barPlot <- function(x, y){
  tryCatch({
    barPlot <- graphics::barplot(x, showCategory = 15) + ggtitle(paste0(filename_cut, " - ", y, "-regulated genes"))
    ggsave(filename = paste0(dirPlot, filename_cut, "_", y, "_barPlot.png"), plot = barPlot)
  }, error = function(msg){
    return(NA)
  })
  }

### Function: Dotplot
enrichres_dotPlot = function(x, y){
  tryCatch({
    dotPlot <- enrichplot::dotplot(x, showCategory = 15) + ggtitle(paste0(filename_cut, " - ", y, "-regulated genes"))
    ggsave(filename = paste0(dirPlot, filename_cut, "_", y, "_dotPlot.png"), plot = dotPlot)},
    error = function(msg){
      return(NA)
    })
  }

### Function: Cnetplot
enrichres_cnetPlot = function(x, y){
  tryCatch({
    cnetPlot <- enrichplot::cnetplot(x = x, showCategory = 5) +
      ggtitle(paste0(filename_cut, " - ", y, "-regulated genes"))
    ggsave(filename = paste0(dirPlot, filename_cut, "_", y, "_cnetPlot.png"), plot = cnetPlot)},
    error = function(msg){
      return(NA)
    })
  }

### Function: Heatplot
enrichres_heatPlot <- function(x, y){
  tryCatch({
    heatPlot <- enrichplot::heatplot(x, showCategory = 5) + ggtitle(paste0(filename_cut, " - ", y, "-regulated genes"))
    ggsave(filename = paste0(dirPlot, filename_cut, "_", y, "_heatPlot.png"), plot = heatPlot)},
    error = function(msg){
      return(NA)
    })
  }

### Function: Treeplot
enrichres_treePlot = function(x, y){
  tryCatch({
    enrichres2 <- pairwise_termsim(x)
    treePlot <- enrichplot::treeplot(enriches2) + ggtitle(paste0(filename_cut, " - ", y, "-regulated genes"))
    ggsave(filename = paste0(dirPlot, filename_cut, "_", y, "_treePlot.png"), plot = treePlot)}, 
    error = function(msg){
      return(NA)
    })
  }

### Function: Gets file names
gets_name <- function(gmt) {
  deparse(substitute(gmt))
  }
########## Set functions ##########


########## Listing GMTs ##########
db_hallmarks    <- "resource/ref/msigdb_v2024.1.Hs_GMTs/h.all.v2024.1.Hs.symbols.gmt"
db_go_all       <- "resource/ref/msigdb_v2024.1.Hs_GMTs/c5.all.v2024.1.Hs.symbols.gmt"
db_reactome_all <- "resource/ref/msigdb_v2024.1.Hs_GMTs/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
db_go_bp        <- "resource/ref/msigdb_v2024.1.Hs_GMTs/c5.go.bp.v2024.1.Hs.symbols.gmt"        

db_reactome_level3    <- "resource/ref/curated/ReactomePathwaysLevel3.gmt"
db_reactome_level3_1  <- "resource/ref/curated/ReactomePathwaysLevel3WithLevel1.tsv"
db_kegg_no_disease    <- "resource/ref/curated/KEGG_pathways_NO_diseases.gmt"
db_btm                <- "resource/ref/curated/BTM_for_GSEA_20131008.gmt"
########## Possible pathways ##########


########## Choose background genes based on GMT's list##########
file_gmt <- db_reactome_level3 # Choose GMT pathway to be tested according to previous names
gmt_name <- gets_name(db_reactome_level3) # Choose pathway to be tested according to previous names
########## Choose background genes ##########


########## Set folders and path ##########
# For input
dirIn <- "workflow/results/DESeq2/"
files_dge <- list.files(path = dirIn, pattern = ".tsv")

# For saving RDS file
rds_path <- "resource/ref/msigdb_v2024.1.Hs_GMTs/" 
########## Set folders and path ##########


########## Analysis ##########
for(file in files_dge){
  df <- read_tsv(file = paste0(dirIn, file))
  #df <- read_tsv("workflow/results/DESeq2/D14_vs_Dm14.tsv")
  names(df)[1] <- "gene_symbol" 
  
  genes_in_data <- df$gene_symbol
  
  bg_full <- read.gmt(file_gmt) # Read the full GMT file
  bg_subset <- bg_full[bg_full$gene %in% genes_in_data, ] # Subset GMT present in the current sample
  background_genes <- gmt_name
  
  rds <- paste0(background_genes, ".RDS")  # Background genes' file name
  saveRDS(bg_subset, paste0(rds_path, rds))
  
  
  df_diffexp <- df[df$DirectionPadj != "None", ] # Remove non-significant genes
  
  df_diffexp <- df_diffexp %>%
    dplyr::select(gene_symbol, 
                  pvalue, 
                  padj, 
                  log2FoldChange, 
                  DirectionPadj)
  
  
  deg_results_list <- split(df_diffexp, df_diffexp$DirectionPadj)
  
  
  bg_genes <- readRDS(paste0("resource/ref/msigdb_v2024.1.Hs_GMTs/", rds))
  padj_cutoff <- 0.05
  genecount_cutoff <- 5 # Minimum number of genes in the pathway, used to filter out pathways
  
  background_genes <- gmt_name
  
  
  # Set folder to save output
  dirOut <- paste0("workflow/results/DESeq2/clusterProfiler/", background_genes, "/")
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
  
  #filename_cut <- "D14_vs_Dm14"
  filename_cut <- base::gsub(pattern = ".tsv", replacement = "", x = file)
  filename <- paste0(dirOut, filename_cut, "_", background_genes)
  
  tryCatch({
    res <- lapply(names(deg_results_list),
                  function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                       TERM2GENE = bg_genes))

    names(res) <- names(deg_results_list)

    res_df <- ldply(res, data.frame)
    res_df <- res_df %>%
      dplyr::mutate(minuslog10padj = -log10(p.adjust))

    target <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff])
    res_df <- res_df[res_df$ID %in% target, ]

    names(res_df)[1] <- "diffexpressed"

  
    # Save Results
    print(paste0('Saving clusterprofiler results: ', filename_cut))
    write.csv(x = res_df, 
              file = paste0(filename, '_resclusterp.csv'),
              row.names = FALSE)

    # ORA up-regulated genes
    res_df_up <- res_df %>% 
      filter(diffexpressed == 'Up') %>% 
      dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
  
    rownames(res_df_up) <- res_df_up$ID
  
    enrichres_up <- new("enrichResult",
                        readable = FALSE,
                        result = res_df_up,
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

  
    # ORA down-regulated genes
    res_df_down <- res_df %>% 
      filter(diffexpressed == 'Down') %>% 
      dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
  
    rownames(res_df_down) <- res_df_down$ID
  
    enrichres_down <- new("enrichResult",
                          readable = FALSE,
                          result = res_df_down,
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
  
  
    # Set folder to save plots
    dirPlot <- paste0("workflow/results/DESeq2/clusterProfiler/", background_genes, "/Plots/")
    if(!file.exists(dirPlot)) { dir.create(path = dirPlot, recursive = TRUE) }
  
  
    # Calling plots and saving 
    ## Up regulataed genes 
    print(paste0('Saving clusterprofiler results for up-regulated genes: ', filename_cut))
  
    enrichres_barPlot(enrichres_up, "up")
    enrichres_dotPlot(enrichres_up, "up")
    enrichres_cnetPlot(enrichres_up, "up")
    enrichres_heatPlot(enrichres_up, "up")
    enrichres_treePlot(enrichres_up, "up")
  
    ## Down regulataed genes 
    print(paste0('Saving clusterprofiler results for down-regulated genes: ', filename_cut))
  
    enrichres_barPlot(enrichres_down, "down")
    enrichres_dotPlot(enrichres_down, "down")
    enrichres_cnetPlot(enrichres_down, "down")
    enrichres_heatPlot(enrichres_down, "down")
    enrichres_treePlot(enrichres_down, "down")}, 
    
    error = function(msg){
      return(NA)
    })
}
########## Analysis ##########





