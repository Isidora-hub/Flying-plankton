# Analysis figures 1 and 2
##### figure 1e is not plotted here since it is a screenshot from CAAF

library(tidyverse)
library(lubridate)
library(VennDiagram)
library(RColorBrewer)
library(ggplotify)
library(patchwork)
library(ggrepel)
library(ggmap)
library(ggspatial)
library(vegan)
library(multcompView)
library(ragg)
library(tmap)
library(rstudioapi)


setwd("~/Data")

locs <- read.csv("Station_ID_location.csv")
load("clean_data_2021.Rdata")


metadata$location <- vapply(strsplit(metadata$sample_name, "_", fixed = TRUE), "[", "", 1)

metadata$container[which(metadata$type == "environmental_aerosol" & metadata$container == "not_applicable")] <- metadata$location[which(metadata$type == "environmental_aerosol" & metadata$container == "not_applicable")]

metadata$location[which(metadata$type == "experimental" | metadata$type == "environmental_aerosol")] <- "Within Site"

summary <- metadata %>% group_by(type, location, container) %>% 
  summarise(samples = n_distinct(sample_name))

summary$lat <- locs$Latitude[match(summary$location,locs$Station.ID)]
summary$long <- locs$Longitude[match(summary$location,locs$Station.ID)]

# clean text
summary$type[which(summary$type == "experimental")] <- "Tank"
summary$type[which(summary$type == "environmental_aerosol")] <- "Aerosol"
summary$type[which(summary$type == "environmental_water")] <- "Environmental Water"
summary$type[which(summary$type == "environmental_soil")] <- "Environmental Soil"


library(rstudioapi)
register_google(key = "XXXXXXXXXXXX") # need an updated key for this to work



########Figure 1a
s_map3 <- get_googlemap(c(-115.7,33.2),
                        zoom = 10, maptype = "satellite")
sat_map3 <-ggmap(s_map3)
sat_map3


#######Figure 1b
s_map <- get_googlemap(center = c(lon = -115.59, lat = 33.19),
                       zoom = 12, maptype = "satellite", 
                       path = "&style=feature:all|element:labels|visibility:off")

sat_map <- ggmap(s_map) + 
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "black")) 

sat_map




######Figure 1c
map <- get_googlemap(center = c(lon = -115.59, lat = 33.19),
                     zoom = 12, maptype = "terrain", 
                     path = "&style=feature:all|element:labels|visibility:off")

p_map <- ggmap(map) 

p3 <- p_map  + 
  geom_point(data = summary %>% filter(type == "Environmental Water"),
             aes(x = long, y = lat, size = samples, fill = type), pch = 21) +
  geom_text(data = summary %>% filter(type == "Environmental Water"),
             aes(x = long, y = lat, label = samples, fill = type), fontface = "bold") +
  scale_size_continuous(trans = "log10", range = c(2,8), limits = c(1,615)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = c(0.2,0.02),
        legend.direction = "horizontal",
        legend.justification = c(0,0),
        legend.background = element_blank()) +
  labs(fill = "", size = "# of Samples", x = "Longitude", y = "Latitude") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#56B4E9"))
p3


########Figure 1d
p4 <- p_map  + 
  geom_point(data = summary %>% filter(type == "Environmental Soil"),
             aes(x = long, y = lat, size = samples, fill = type), pch = 21) +
  geom_text(data = summary %>% filter(type == "Environmental Soil"),
            aes(x = long, y = lat, label = samples, fill = type), fontface = "bold", color = "white") +
  scale_size_continuous(trans = "log10", range = c(2,8), limits = c(1,615)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = c(0.2,0.02),
        legend.direction = "horizontal",
        legend.justification = c(0,0),
        legend.background = element_blank()) +
  labs(fill = "", size = "# of Samples", x = "Longitude", y = "Latitude") +
  guides(fill = "none") +
  scale_fill_manual(values = c("#8B4513"))
p4




##########Venn diagram and diversity analysis

bact <- taxonomy$Hash[which(taxonomy$A == "Bacteria")]
bact_tab <- asv_table[,match(bact,colnames(asv_table))]

# Clean Euks a bit more for only protists
euk_tax <- taxonomy %>% filter(A %in% c("Eukaryota", "Chloroplast"),
                               B %in% c(" Apicomplexa", " Apusomonadidae", "Centroheliozoa", " Cercozoa", 
                                        " Chlorophyta", " Choanoflagellida", " Ciliophora", " Conosa"," Cryptophyta",
                                        " Dinoflagellata", " Discoba", " Glaucophyta", " Hacrobia_X",
                                        " Haptophyta", " Hilomonadea"," Katablepharidophyta", " Lobosa", " Mesomycetozoa", " Metamonada",
                                        " Ochrophyta", " Perkinsea", " Protalveolata_X", " Stramenopiles_X", 
                                        " Streptophyta", " Telonemia"),
                               !C %in% c(" Charophyceae"," Embryophyceae", " Klebsormidiophyceae"," Phaeophyceae", " Ulvophyceae"))

euks <- euk_tax$Hash
euk_tab <- asv_table[,match(euks,colnames(asv_table))]

shannon_div <- diversity(asv_table)

shannon_bact <- diversity(bact_tab)
shannon_euks <- diversity(euk_tab)

asv_table$sample <- rownames(asv_table)

metadata$shannon <- shannon_div[match(metadata$sample_name, asv_table$sample)]

metadata$shannon_bac <- shannon_bact[match(metadata$sample_name, asv_table$sample)]
metadata$shannon_euk <- shannon_euks[match(metadata$sample_name, asv_table$sample)]


######Figure 2a
# Venn Diagram
piv_table <- asv_table %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

piv_table$type <- metadata$type[match(piv_table$sample,metadata$sample_name)]

piv_table <- piv_table %>% group_by(type,Hash) %>% summarise(sum_reads = sum(reads)) %>% 
  filter(sum_reads != 0)

exper <- piv_table$Hash[which(piv_table$type == "experimental")]
aero <- piv_table$Hash[which(piv_table$type == "environmental_aerosol")]
water <- piv_table$Hash[which(piv_table$type == "environmental_water")]
soil <- piv_table$Hash[which(piv_table$type == "environmental_soil")]

out <- venn.diagram(
  x = list(exper, aero, water,soil),
  category.names = c("Tank colonizers" , "Aerial dispersers" ,
                     "Aquatic source", "Terrestrial source"),
  filename = NULL,
  units = "in",
  height = 10, 
  width = 6, 
  resolution = 400,
  compression = "lzw",
  
  col=c("#009E73", '#E69F00', '#56B4E9', '#8B4513'),
  fill = c(alpha("#009E73",0.6), alpha('#E69F00',0.6),
           alpha('#56B4E9',0.6),alpha('#8B4513',0.6)),
  
  cex = 0.75,
  fontfamily = "sans",
  fontface = "plain",
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 0.75,
  cat.pos = c(-15,15,-20,20),
  cat.dist = c(0.22,0.22,0.15,0.15),
  cat.fontfamily = "sans",
  rotation.degree = 0,
  margin = 0.05
)

wrap_elements(grobTree(out))




metadata$type[which(metadata$type == "experimental")] <- "Tank\ncolonizers"
metadata$type[which(metadata$type == "environmental_aerosol")] <- "Aerial\ndispersers"
metadata$type[which(metadata$type == "environmental_water")] <- "Aquatic\nsource"
metadata$type[which(metadata$type == "environmental_soil")] <- "Terrestrial\nsource"

metadata$type <- factor(metadata$type,
                        levels = c("Terrestrial\nsource","Aquatic\nsource", "Aerial\ndispersers", "Tank\ncolonizers"))

bac_aov <- aov(shannon_bac ~ type, data = metadata)
euk_aov <- aov(shannon_euk ~ type, data = metadata)

# do manually since we want highest mean to be first letter
bac_letters <- multcompLetters(TukeyHSD(bac_aov)$type[,4])$Letters %>% toupper()
bac_letters[1:4] <- c("B","B","C","A")

euk_letters <- multcompLetters(TukeyHSD(euk_aov)$type[,4])$Letters %>% toupper()

metadata$bac_letters <- bac_letters[match(metadata$type, names(bac_letters))]
metadata$euk_letters <- euk_letters[match(metadata$type, names(euk_letters))]


######Figure 2b
p6 <- ggplot(metadata, aes(x = type, y = shannon_bac, fill = type)) +
  geom_violin(draw_quantiles = 0.5) + geom_jitter(alpha = 0.5, width = 0.25, size = 0.5) +
  geom_text(aes(x = type, y = max(shannon_bac + 0.4), label = bac_letters)) +
  scale_fill_manual(values = c(c( '#8B4513','#56B4E9','#E69F00', "#009E73"))) +
  labs(x = "", y = "Bacterial Shannon Diversity", fill = "") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "none")

######Figure 2c
p7 <- ggplot(metadata, aes(x = type, y = shannon_euk, fill = type)) +
  geom_violin(draw_quantiles = 0.5) + geom_jitter(alpha = 0.5, width = 0.25, size = 0.5) +
  geom_text(aes(x = type, y = max(shannon_euk + 0.4), label = euk_letters)) +
  scale_fill_manual(values = c(c( '#8B4513','#56B4E9','#E69F00', "#009E73"))) +
  labs(x = "", y = "Eukaryotic Shannon Diversity", fill = "") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "none")





