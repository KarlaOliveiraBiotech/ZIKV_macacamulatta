## Importing libraries
library(tidyverse)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)


## Importing data and specifications
raw_counts <- readr::read_tsv("results/feature_counts/raw_count_matrix_ensemblid.tsv",
                              col_names = TRUE, skip = 1)

sample_spec <- readr::read_tsv(file = "resource/data/phenodata.tsv", col_names = TRUE)


## Preparing data
df <- base::subset(raw_counts, select = -7)


### Smart way to change columns' names
cn <- base::gsub(pattern = "../results/bowtie2/mapped/bam/(.*)_S(.*)_good_clean_2M.bam", 
                 replacement = "\\1", x = colnames(df))
colnames(df) <- cn

### Creating a df of counts that will be read by EdgeR
count_df <- df %>% 
  dplyr::select(-c(2:6))


### Sorting df
count_df_sorted <- subset(count_df, 
                          select = c(1, 16, 27, 8, 12, 14, 
                                     5, 4, 3, 6, 9, 20, 24, 
                                     2, 22, 7, 25, 18, 15, 
                                     11, 23, 26, 13, 19, 10, 
                                     17, 21, 28
                                     )
                          )


### Sorting sample spec
splitted <- strsplit(sample_spec$Sample, "_") 
splitted <- sapply(splitted, "[[",1)
splitted <- as.numeric(splitted)

sample_spec_order <- dplyr::mutate(sample_spec, splitted) 
sample_spec_order <- sample_spec_order[order(splitted),] 

sample_spec_order <- sample_spec_order %>% 
  select(-Filename) %>% 
  select(-spl)
sample_spec_order

group <- factor(paste(sample_spec_order$Class, sample_spec_order$Timepoint, 
                      sep = "."))
sample_spec_order <- cbind(sample_spec_order, group = group)

sample_spec_order


design <- model.matrix(~ 0 + group)
y <- DGEList(counts = count_df_sorted, group = group)
keep <- filterByExpr(y, group = "Dm14") # verficar se está certo
y <- y[keep, , keep.lib.size = FALSE]
y <- normLibSizes(y)



my_contrasts <- makeContrasts(D1_vs_Dm14 = groupZIKV.D1 - groupbaseline.Dm14, 
                              D3_vs_Dm14 = groupZIKV.D3 - groupbaseline.Dm14, 
                              D5_vs_Dm14 = groupZIKV.D5 - groupbaseline.Dm14, 
                              D7_vs_Dm14 = groupZIKV.D7 - groupbaseline.Dm14, 
                              D10_vs_Dm14 = groupZIKV.D10 - groupbaseline.Dm14, 
                              D14_vs_Dm14 = groupZIKV.D14 - groupbaseline.Dm14, 
                              levels = design)



fit <- glmQLFit(y, design)
colnames(fit)

D1_vs_Dm14 <- glmQLFTest(fit, contrast = my_contrasts[, "D1_vs_Dm14"])
topTags(D1_vs_Dm14)

### OR 
### Contrasts will be done according to order gave by colnames(fit). 
### For example

### colnames(fit)
### [1] "groupbaseline.Dm14" "groupZIKV.D1"       "groupZIKV.D10"     
### [4] "groupZIKV.D14"      "groupZIKV.D3"       "groupZIKV.D5"      
### [7] "groupZIKV.D7"

### To compare D1_vs_Dm14 and D3_vs_Dm14 in this way use: 

### D1_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0))
### D3_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 1, 0, 0))


___________
### D1_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0))
### D3_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 1, 0, 0))
### D5_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 0, 1, 0))
### D7_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 0, 0, 0, 1))
### D10_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 1, 0, 0, 0, 0))
### D14_vs_Dm14 <- glmQLFTest(fit, contrast = c(-1, 0, 0, 1, 0, 0, 0))

