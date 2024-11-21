rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("dplyr")
library("readr")
library("ggplot2")
library("pheatmap")
library("data.table")
library("edgeR")
library("mdp")
library("CEMiTool")
########## Load libraries ##########


########## Import data  ##########
raw_counts <- readr::read_tsv(file = "workflow/results/feature_counts/raw_count_matrix_symbol.tsv")

annot <- readr::read_tsv(file = "resource/data/phenodata.tsv")

pathways <- read_gmt("resource/ref/curated/ReactomePathwaysLevel3.gmt")
pathways_list <- as.list(pathways)

db_reactome_level3_1 <- readr::read_tsv("resource/ref/curated/ReactomePathwaysLevel3WithLevel1.tsv")
########## Import data  ##########


# Subsets immune system pathways
immune_system <- db_reactome_level3_1 %>% 
  dplyr::filter(Level1 == "Immune_System") %>% 
  dplyr::select(Pathway) %>% 
  base::unlist(use.names = TRUE)


db_immune <- pathways %>% 
  dplyr::filter(term %in% immune_system)

immune_list <- as.list(db_immune)
########## File for comparison (immune system) ##########



########## Prep data ##########
spl <- strsplit(x = annot$Sample, split = "_")
spl <- sapply(spl, "[[", 1)
spl <- as.numeric(spl)
annot <- dplyr::mutate(annot, spl)
annot <- annot[order(spl), ]
annot <- annot %>% dplyr::select(-c(5,6))
old_annot <- annot

# Removes non-symbol genes
expr_unsort <- raw_counts
expr <- expr_unsort[!grepl("ENSMMU", expr_unsort$Symbol), ]

# Removes duplicated genes
expr <- expr[!duplicated(expr$Symbol), ]
anyDuplicated(expr)

# Transform "Symbol" into rownames
expr <- expr %>% 
  tibble::column_to_rownames("Symbol") %>% 
  as.data.frame()

# Normalizes counts by CPM
expr_cpm <- cpm(y = expr, 
              log = TRUE, 
              prior.count = 2, 
              normalized.lib.sizes = TRUE)

expr <- expr_cpm %>% 
  as.data.frame()

# Sorts expression dataframe by samples' names
order <- old_annot$Sample
expr <- expr[, order]

# Adds "S" to sample' names (necessary for MDP analysis)
colnames(expr) <- paste0("S", colnames(expr))
annot$Sample <- paste0("S", annot$Sample)

# Adjusts inputs for MDP
annot$Sample <- gsub(pattern = "-", replacement = "\\.", annot$Sample)
colnames(expr) <- gsub(pattern = "-", replacement = "\\.", x = colnames(expr))

names(annot)[2] <- "ClassBckp"

# Adds 01--07 to order timepoints
annot$New <- c("01-Dm14", "01-Dm14", "01-Dm14", "01-Dm14", 
                  "02-D1", "02-D1", "02-D1", "02-D1", 
                  "03-D3", "03-D3", "03-D3", "03-D3", 
                  "04-D5", "04-D5", "04-D5", 
                  "05-D7", "05-D7", "05-D7", "05-D7", 
                  "06-D10", "06-D10", "06-D10", "06-D10", 
                  "07-D14", "07-D14", "07-D14", "07-D14")

names(annot)[5] <- "Class"
########## Prep data ##########


########## MDP analysis ##########
# Creates output directory when there is pathway list
dirOut_path <- "workflow/results/mdp/path/"
if(!file.exists(dirOut_path)) { dir.create(path = dirOut_path, recursive = TRUE) }

mdp_results <- mdp(data = expr, 
                   pdata = annot, 
                   control_lab = "01-Dm14",
                   directory = dirOut_path, 
                   pathways = pathways_list, 
                   print = TRUE,
                   measure = "median", 
                   std = 2, 
                   fraction_genes = 0.25,
                   save_tables = TRUE, 
                   file_name = "pathways_")

# Creates output directory when there is no pathway list
dirOut <- "workflow/results/mdp/no_path/"
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = TRUE) }

# Runs MDP when there is no pathway list
mdp_results <- mdp(data = expr, 
                   pdata = annot, 
                   control_lab = "01-Dm14", 
                   directory = dirOut, 
                   print = TRUE,
                   measure = "median", 
                   std = 2, 
                   fraction_genes = 0.25,
                   save_tables = TRUE, 
                   file_name = "no_pathways_")


# Creates output directory when there is immune genes
dirOut_immune <- "workflow/results/mdp/immune/"
if(!file.exists(dirOut_immune)) { dir.create(path = dirOut_immune, recursive = TRUE) }

mdp_results <- mdp(data = expr, 
                   pdata = annot, 
                   control_lab = "01-Dm14", 
                   directory = dirOut_immune,
                   pathways = immune_list, 
                   print = TRUE,
                   measure = "median", 
                   std = 2, 
                   fraction_genes = 0.25,
                   save_tables = TRUE, 
                   file_name = "immune_")
########## MDP analysis ##########