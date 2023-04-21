library(tidyverse)

# data processing

asv_table <- read.csv("Data/ASV_counts_5mincount_withTax_NoNorm_2021.05.21_IE_v2.csv")

metadata <- read.csv("Data/Sample metadata_2021.07.csv")

# pull out taxonomy

taxonomy <- asv_table[,1:2]

taxonomy <- separate(taxonomy, phylogeny, sep = ";", into = c("A","B","C", "D", "E", "F", "G"), remove = FALSE)

# asv table

hash <- asv_table$Hash

samples <- colnames(asv_table[,-c(1:2)])

asv_table <- t(asv_table[,-c(1:2)]) %>% as.data.frame()

colnames(asv_table) <- hash
rownames(asv_table) <- samples

# match metadata

sample_name <- str_split_fixed(rownames(asv_table), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]

fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)

sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names

rownames(asv_table) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(asv_table), metadata$sample_name),]

save(asv_table, taxonomy, metadata, file = "Data/clean_data_2021.Rdata")

