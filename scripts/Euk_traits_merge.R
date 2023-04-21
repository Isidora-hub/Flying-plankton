#### Merging old with new eukaryote table to complete trait table

####load new ASV with taxonomy
asv_new <- read.table("Data/ASV_counts_5mincount_withTax_NoNorm_2021.05.21_IE_v2.csv", sep = ";", header = TRUE)

#install.packages("stringr")
library("stringr")
library("tidyverse")

####Separate
asv_new_tax <- str_split_fixed(asv_new$phylogeny, ";", 7)

asv_tax_table <- cbind(asv_new, Tax = asv_new_tax) 

asv_tax_table$phylogeny <- NULL

levels(as.factor(asv_tax_table$Tax.1))

#######only euk
euk<-subset(asv_tax_table, Tax.1 !="Archaea" & Tax.1 !="Bacteria" & Tax.1 !="Mitochondria")
euk<-subset(euk, Tax.2 !=" Fungi" & Tax.2 != " Metazoa" & Tax.2 != " NA") ###take in consideration space between " and F
euk<-subset(euk, Tax.3 !=" Embryophyceae") ###take in consideration space between " and E
euk<-subset(euk, Tax.6 !=" Ulva" & Tax.6 !=" Chara") ###take in consideration space between " and E


####upload old trait table
Traits_drive<- read.csv("~/Documents/Post-doc/02 Shurin Lab/Flying plankton/Analysis/Data 2019/Euk_traits_mincount5k_ASV50.csv", sep=",") ###this is the table of traits filled on the google drive

#######merge old with new
Traits_new <- merge(Traits_drive, euk, by = "Hash", all.x = FALSE, all.y = TRUE)
write.csv(Traits_new, file = "New_Traits_euk_drive.csv")

save(Traits_new, Traits_drive, file = "Data/Euk_Trait.Rdata")



