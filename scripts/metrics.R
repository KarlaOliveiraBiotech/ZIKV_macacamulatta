# R scripts must finsh in capital R
# install.packages('tidyverse')  # important to comment the "install". It didn't work otherwise
library(tidyverse)


df <- read.table(file = snakemake@input[[1]], header = TRUE)  # same structure
print(df)

df %>%
count(Sample)

df %>%
count(Class)

df %>%
count(Timepoint)