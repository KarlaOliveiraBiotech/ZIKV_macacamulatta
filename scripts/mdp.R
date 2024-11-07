rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
library("CEMiTool")
library("mdp")
########## Load libraries ##########


########## Import data  ##########
raw_counts <- readr::read_tsv(file = "workflow/results/feature_counts/raw_count_matrix_symbol.tsv", 
                              col_names = TRUE)

lib_spec <- readr::read_tsv(file = "resource/data/phenodata.tsv", 
                            col_names = TRUE)

pathways <- read_gmt("resource/ref/curated/ReactomePathwaysLevel3.gmt")
pathways_list <- as.list(pathways)
########## Import data  ##########


########## Prep data ##########
spl <- strsplit(x = lib_spec$Sample, split = "_")
spl <- sapply(spl, "[[", 1)
spl <- as.numeric(spl)

lib_spec <- dplyr::mutate(lib_spec, spl)
lib_spec <- lib_spec[order(spl), ]
lib_spec <- lib_spec %>% dplyr::select(-c(5,6))

old_lib_spec <- lib_spec

# Creates an annotation object
annot <- lib_spec %>% 
  dplyr::select(Sample, Timepoint)

# Add "S" to samples' name: MDP requirement
annot$Sample <- paste0("S", annot$Sample)
names(annot) <- c("Sample", "Class")

# Coerce annotation to dataframe format
annot <- data.frame(annot)

# Removes non-symbol genes 
df_unsort <- raw_counts
df <- df_unsort[!grepl("ENSMMU", df_unsort$Symbol), ]

# Removes duplicated genes
df <- df[!duplicated(df$Symbol), ]
anyDuplicated(df)

# Transform gene symbols to rownames
df <- data.frame(df, row.names = 1)
cn2 <- base::gsub(pattern = "X(.*)", replacement = "\\1", x = colnames(df))
cn2 <- base::gsub(pattern = "(00.)", replacement = "\\00-", x = cn2)
colnames(df) <- cn2

# Sorts expression dataframe by samples' names
order <- old_lib_spec$Sample
df <- df[, order]

# Add "S" to samples' names: MDP requirement
colnames(df) <- paste0("S", colnames(df))

# Coerce expression table as dataframe format
df <- data.frame(df)
########## Prep data ##########


########## MDP analysis ##########
# Creates output directory when there is pathway list
dirOut_path <- "workflow/results/mdp/path/"
if(!file.exists(dirOut_path)) { dir.create(path = dirOut_path, recursive = TRUE) }

# Runs MDP when there is pathway list
mdp_results <- mdp(data = df, 
                   pdata = annot, 
                   control_lab = "Dm14", 
                   directory = dirOut_path, 
                   pathways = pathways_list, 
                   print = TRUE,
                   measure = "mean", 
                   std = 2, 
                   fraction_genes = 0.25,
                   save_tables = TRUE, 
                   file_name = "")


# Creates output directory when there is no pathway list
dirOut <- "workflow/results/mdp/no_path/"
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }

# Runs MDP when there is no pathway list
mdp_results <- mdp(data = df, 
                   pdata = annot, 
                   control_lab = "Dm14", 
                   directory = dirOut, 
                   print = TRUE,
                   measure = "mean", 
                   std = 2, 
                   fraction_genes = 0.25,
                   save_tables = TRUE, 
                   file_name = "no_pathways_")
########## MDP analysis ##########
