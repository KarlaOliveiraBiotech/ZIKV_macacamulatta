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

folder <- "workflow/results/edgeR/"
files <- list.files(path = folder, pattern = ".tsv")

mart <- useDataset("mmulatta_gene_ensembl", 
                   mart = useMart("ensembl"))

bm <- getBM(attributes = c("ensembl_gene_id", 
                           "external_gene_name", 
                           "hsapiens_homolog_associated_gene_name"), 
            mart = mart)

ens2symbol <- bm %>% 
  dplyr::select(ENSEMBL = ensembl_gene_id,SYMBOL = hsapiens_homolog_associated_gene_name) %>% 
  mutate(SYMBOL = if_else(condition = SYMBOL == "", true = ENSEMBL, false = SYMBOL))

ens2symbol <- ens2symbol[!grepl("ENSMMU", ens2symbol$SYMBOL), ]



db_hallmarks  <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/h.all.v2024.1.Hs.symbols.gmt")
db_go         <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/c5.all.v2024.1.Hs.symbols.gmt")
db_reactome   <- gmtPathways("resource/ref/msigdb_v2024.1.Hs_GMTs/c2.cp.reactome.v2024.1.Hs.symbols.gmt")

for(file in files){
  filename <- sub(pattern = "(.*).tsv", replacement = "\\1", basename(file))
  df <- read_tsv(paste0(folder, file))
  
  res <- inner_join(x = df, ens2symbol, by = c("Ensembl_id" =  "ENSEMBL")) %>% 
    dplyr::select(SYMBOL, logFC, PValue, FDR) %>% 
    na.omit() %>% 
    distinct()
  
  res_uniq <- res[!duplicated(res[, "SYMBOL"]), ]
  
  ranks <- data.frame(SYMBOL = res_uniq$SYMBOL, Log2FC = res_uniq$logFC * -log10(res_uniq$PValue))
  
  ranks <- ranks[order(-ranks$Log2FC), ]
  ranks <- ranks[!is.na(ranks$Log2FC), ]
  
  ranks2 <- deframe(ranks)
  
  fgseaRes <- fgsea(pathways = db_go, stats = ranks2)
  
  dirOut <- paste0("workflow/results/edgeR/fGSEA/")
  if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }
  
  
  # Output files - save
  readr::write_tsv(x = data.frame(fgseaRes),
                   file = paste0(dirOut, filename, "_fGSEA.tsv"))
  xlsx::write.xlsx(x = data.frame(fgseaRes[, 1:7]),
                   file = paste0(dirOut, filename, "_fGSEA.xlsx"), row.names = FALSE)

  
  # Tidy results
  fgseaResTidy <- fgseaRes %>% 
    as_tibble() %>% 
    arrange(desc(NES))
  
  # Presenting in a nice table
  fgseaResTidy %>% 
    dplyr::select(-c(leadingEdge, ES, log2err)) %>% 
    arrange(padj) %>% 
    DT::datatable()
  
}
