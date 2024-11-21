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


########## Import data  ##########
raw_counts <- readr::read_tsv(file = "workflow/results/feature_counts/raw_count_matrix_symbol.tsv", 
                              col_names = TRUE)

lib_spec <- readr::read_tsv(file = "resource/data/phenodata.tsv", 
                            col_names = TRUE)
db <- as.data.frame(read_gmt("resource/ref/curated/ReactomePathwaysLevel3.gmt"))


int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
########## Import data  ##########


########## Prep data ##########
spl <- strsplit(x = lib_spec$Sample, split = "_")
spl <- sapply(spl, "[[", 1)
spl <- as.numeric(spl)

lib_spec <- dplyr::mutate(lib_spec, spl)
lib_spec <- lib_spec[order(spl), ]
lib_spec <- lib_spec %>% dplyr::select(-c(5,6))

old_lib_spec <- lib_spec

expr_unsort <- raw_counts
expr <- expr_unsort[!grepl("ENSMMU", expr_unsort$Symbol), ]
expr <- expr[!duplicated(expr$Symbol), ]
anyDuplicated(expr)

# Set gene symbols as rownames and order samples
expr <- data.frame(expr, row.names = 1)
cn2 <- base::gsub(pattern = "X(.*)", replacement = "\\1", x = colnames(expr))
cn2 <- base::gsub(pattern = "(00.)", replacement = "\\00-", x = cn2)
colnames(expr) <- cn2

order <- old_lib_spec$Sample
expr <- expr[, order]

# Creates annotation df
annot <- lib_spec %>% 
  dplyr::select(Sample, Timepoint) %>% 
  as.data.frame()

# Removes duplicated genes
expr2 <- expr[!duplicated(expr), ]

# Normalize expression 
expr_norm <- edgeR::cpm(y = expr2, 
                        log = TRUE, 
                        prior.count = 2, 
                        normalized.lib.sizes = TRUE) %>% 
  as.data.frame()

# Formats samples' name
colnames(expr_norm) <- paste0("S", colnames(expr_norm))
annot$Sample <- paste0("S", annot$Sample)
########## Prep data ##########



# Force to be dataframes
# expr_norm <- as.data.frame(expr_norm)
# annot <- as.data.frame(annot)
# db <- as.data.frame(db)


########## Run CEMiTool ##########
cem <- cemitool(expr = expr_norm, annot = annot, gmt = db, interactions = int_df,
                filter = TRUE, apply_vst = TRUE, network_type = "signed",
                cor_method = "pearson", force_beta = TRUE,
                sample_name_column = "Sample", class_column = "Timepoint",
                gsea_min_size = 15, gsea_max_size = 500,
                merge_similar = FALSE, verbose = TRUE)



# Creates dirOut and saves results accordingly
dirOut <- "workflow/results/CEMiTool/"
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }  

generate_report(cem,   force = T, directory = dirOut)
write_files(cem,       force = T, directory = paste0(dirOut, "/Tables"))
save_plots(cem, "all", force = T, directory = paste0(dirOut, "/Plots"))
diagnostic_report(cem, force = T, directory = paste0(dirOut, "/Diagnostic"))
########## Run CEMiTool ##########

