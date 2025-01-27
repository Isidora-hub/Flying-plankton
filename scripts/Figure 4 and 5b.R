
######CLEANED VERSION

#####Eukaryotes Taxonomy and Traits Analysis
#####using two main csv files, one with traits associated to each ASV (Traits_Feb2022.csv) 
####and a second one associated to taxonomy information (EukTax_Feb2022.csv)
##### for the traits dataset MBFG is not used, still included in this code

library(ade4)
library(ggplot2)
library(vegan)
library(ape)
library(devtools)
library(tidyverse) 
library(DataExplorer)
library(reshape2)
library(dplyr)
library(rlang)
library(RVAideMemoire)
library(agricolae)
library(grid)
library(tidyr)
library(compositions)


setwd("~/Data")
load("clean_data_2021.Rdata")


######Figure 4a
###########################
#######Taxa radar plots

Euk<- read.csv("~/Data/EukTax_Feb2022.csv", sep=",")
names(Euk)

#####keeping only taxa info and counts
Euk<-Euk[,-c(9:61)]

Euk[,c(2:8)]<-lapply(Euk[,c(2:8)], as.factor)



####Tax_3 (or division)
Euk_div<-Euk[,-c(1:3,5:8)]
Euk_div_agg<-aggregate(.~Tax.3, Euk_div, sum)
Euk_div_agg2<-as.data.frame(t(as.matrix(Euk_div_agg[,c(2:1078)]))) 
names(Euk_div_agg2)<-levels(as.factor(Euk_div_agg$Tax.3))
Euk_div_agg2 <- Euk_div_agg2[which(rowSums(Euk_div_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Euk_div_agg2)) {
  print(i)
  Euk_div_agg2[i,] <- as.numeric(Euk_div_agg2[i,])/sum(as.numeric(Euk_div_agg2[i,]), na.rm = TRUE)*100
  
}

# match metadata
sample_name <- str_split_fixed(rownames(Euk_div_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Euk_div_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Euk_div_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Euk_div_agg2$sample_name<-rownames(Euk_div_agg2)


data_div<-dplyr::full_join(Euk_div_agg2, data2, by="sample_name")
data_div[,c(1:77)]<-lapply(data_div[,c(1:77)], as.numeric)

data_div <- data_div[,-c(78)]

data_div_avg<-aggregate(.~type, data=data_div, mean)

#keep only taxa above 10%
data_div_avg10<-data_div_avg[,c(1,2,7,12,19,21,68,73)]


######radarchart
#######over 10%
data11 <- as.data.frame(t(data_div_avg10)) #changes columns to rows, here we need MBFG on rows and sample types in columns
colnames(data11) <-c("Aerosol", "Soil","Water","Tank" ) ####spaces help to center side labels in radar plot
data11 <-data11[-c(1),]



data11 <- rbind(rep(60,4) , rep(0,4) , data11) # To use the fmsb package, add 2 lines to the dataframe: the max and min of each variable to show on the plot

data11[,c(1:4)]<-lapply(data11[,c(1:4)], as.numeric)#spaces apear after transpose, this fix the problem, first transform to numeric
data11 <-as.data.frame(data11) #then back to dataframe


###radarchart
# Library
library(fmsb)

# Customize with legend

# Color vector
colors_border= c("darkorchid4","gold","chartreuse3","red", "darkblue", "indianred4","chartreuse4") 

pdf("Radar_Euk_Division10_nofungiNA_March23.pdf", width = 14, height = 8, pointsize=12)
radarchart(data11  , axistype=1,
           #custom polygon
           pcol=colors_border, plwd=4 , plty=1,
           #custom the grid
           cglcol="black", cglty=1, axislabcol="black", cglwd=0.8,caxislabels=c("0%","15%","30%","45%","60%"),
           #custom labels
           vlcex=2.4, vlabels=c("Aerial dispersers", "Terrestrial source                    ","Aquatic source","                 Tank colonizers" )) ####spaces help to center side labels in radar plot


# Add a legend
legend(x=1, y=1.4, legend = rownames(data11[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.8, pt.cex=3)
dev.off()



#######Figure 4b
#########Taxa PcoA
######remember to enter row names
rownames(Euk) <- Euk$Hash


Euk<-Euk[,-c(1:8)]
Euk <- Euk[which(rowSums(Euk) > 0), ] #######deletes rows summ 0
Euk <- Euk[,which(colSums(Euk) > 0) ] #######deletes columns summ 0

#######transform to relative abundances

for (i in 1:ncol(Euk)) {
  print(i)
  Euk[,i] <- as.numeric(Euk[,i])/sum(as.numeric(Euk[,i]), na.rm = TRUE)*100
  
}


data<-t(Euk)
data_Euk <-as.data.frame(data)



#### PCoA taxa Euk

library(RVAideMemoire)

######REMEMBER TO CHANGE NUMBER OF VARIABLES in trait_df below

compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  trait_df = trait_df[, c(1:3129)]
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  
  # Add values of covariances inside object
  given_pcoa$U <- U
  
  return(given_pcoa)
}

# match metadata
sample_name <- str_split_fixed(rownames(data_Euk), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(data_Euk) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(data_Euk), metadata$sample_name),]



######distance matrix
data_dm<-vegan::vegdist(data_Euk,method="bray",diag=TRUE,upper=TRUE)
data_distancematrix<-as.matrix(data_dm)

#####PcoA from ape is not working, I try other solution proposed

source("https://raw.githubusercontent.com/emmanuelparadis/ape/master/R/pcoa.R")

result<-pcoa(data_distancematrix)
summary(result)

result$values
pcoa_vectors<-result$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$sample_name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL


head(pcoa_vectors)


dim(metadata)
dim(pcoa_vectors)
data<-dplyr::left_join(pcoa_vectors,metadata,by="sample_name") ####left_join(): includes all rows in x.



length(data$type)
vegan::adonis2(data_dm ~ data$type)


#Call:
#vegan::adonis(formula = data_dm ~ data$type) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data$type    3     48.01 16.0037  40.102 0.10459  0.001 ***
#Residuals 1030    411.04  0.3991         0.89541        
#Total     1033    459.05                 1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Plotting the PCoA
trait_pcoa_arrows = compute_arrows(result, data_Euk)



# /!\ Be sure to change $vectors to a data.frame before putting it in ggplot2 /!\

vectors <- as.data.frame(trait_pcoa_arrows$vectors)
vectors$sample_name<-rownames(vectors)
rownames(vectors)<-NULL


head(vectors)
vectors<-dplyr::left_join(vectors,metadata,by="sample_name")
head(vectors)

####ARROWS!!!!
#arrows_df<-as.data.frame(trait_pcoa_arrows$U/10000) #this was the original, can be customized
#arrows_df<-as.data.frame(trait_pcoa_arrows$U) #####this is the data frame that will help customize and plot the arrows in the PcoA
arrows_df<-as.data.frame(trait_pcoa_arrows$U/100) 

barplot(result$values$Relative_eig[1:10]) # plot the eigenvalues and interpret
PCOAaxes <- result$vectors[,c(1,2)] ## Extract the plot scores from first two PCoA axes (if you need them):
# Now plot a bar plot of relative eigenvalues. This is the percentage variance explained by each axis
barplot((result$values$Eigenvalues)/sum(result$values$Eigenvalues)) 

# Calculate the percent of variance explained by first two axes
sum(as.vector((result$values$Eigenvalues)/sum(result$values$Eigenvalues))[1:2])
#[1] 0.2353058
#or by each axes 
sum(as.vector((result$values$Eigenvalues)/sum(result$values$Eigenvalues))[1])
sum(as.vector((result$values$Eigenvalues)/sum(result$values$Eigenvalues))[2])
#[1] 0.1511735
#[1] 0.08413225


##########Plot

p1_2 <- ggplot(NULL, aes(Axis.1, Axis.2)) +
  geom_point(data = vectors,size=2,aes(color=type)) +
  stat_ellipse(data = vectors, geom = "polygon", aes(color=type), alpha = 0, show.legend = FALSE, level = 0.95, size=1.1) +
  #geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,mapping=aes(xend=Axis.1,yend=Axis.2),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.7,yend=-0.28),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.37,yend=-0.137),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.096,yend=0.646),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.06,yend=0.39),arrow=arrow(length=unit(3,"mm")))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 25),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=17),legend.text=element_text(size=17), plot.title=element_text(size=30,face="bold"))+
  scale_shape_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c(15,16,17,18))+
  scale_color_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  xlab("\nAxis 1 (15.1% of variance)")+
  ylab("Axis 2 (8.4% of variance)\n")+
  annotate("text", x=-0.78, y=0.75,label=paste("list(Adonis:~F['3,1033']==40.1", ", p==0.001", ", R^{2}==0.105)"),parse=TRUE,hjust=0.0,size=9.5)+
  annotate("text", x=-0.6,y=-0.3, label="Cymbellaceae",size=9)+
  annotate("text",x=-0.55,y=-0.15, label="Amphora sp.",size=9)+
  annotate("text",x=-0.15,y=0.65, label="Cymbellaceae",size=9)+
  annotate("text",x=-0.1,y=0.4, label="Amphora coffeaeformis",size =9)


p1_2 +labs(title = "Eukaryotes PCoA (taxonomy)")


pdf("Euk_taxa_PcoA_March23.pdf", width = 12, height = 8, pointsize=12)
p1_2 +labs(title = "Eukaryotes PCoA (taxonomy)")
dev.off()





Traits<- read.csv("~/Data/Traits_Feb2022.csv", sep=",")

Traits[,c(3:17)]<-lapply(Traits[,c(3:17)], as.factor)

#######Trophy

Trophy<-Traits[,-c(1:3,5:17)]
Trophy_agg<-aggregate(.~Trophy, Trophy, sum)
Trophy_agg2<-as.data.frame(t(as.matrix(Trophy_agg[,c(2:1078)]))) 
names(Trophy_agg2)<-levels(as.factor(Trophy_agg$Trophy))
Trophy_agg2 <- Trophy_agg2[which(rowSums(Trophy_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Trophy_agg2)) {
  print(i)
  Trophy_agg2[i,] <- as.numeric(Trophy_agg2[i,])/sum(as.numeric(Trophy_agg2[i,]), na.rm = TRUE)*100
  
}


# match metadata
sample_name <- str_split_fixed(rownames(Trophy_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Trophy_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Trophy_agg2), metadata$sample_name),]



##########
data2 <- metadata[,c(2,4)]
Trophy_agg2$sample_name<-rownames(Trophy_agg2)


data_Trophy<-dplyr::full_join(Trophy_agg2, data2, by="sample_name")
colnames(data_Trophy) <-c("NA","A", "H", "LTR", "M","sample_name","type")


#######Size Median
Size_Md<-Traits[,-c(1:5,7:17)]
Size_Md_agg<-aggregate(.~Md.Size.Classification, Size_Md, sum)
Size_Md_agg2<-as.data.frame(t(as.matrix(Size_Md_agg[,c(2:1078)]))) 
names(Size_Md_agg2)<-levels(as.factor(Size_Md_agg$Md.Size.Classification))

Size_Md_agg2 <- Size_Md_agg2[which(rowSums(Size_Md_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Size_Md_agg2)) {
  print(i)
  Size_Md_agg2[i,] <- as.numeric(Size_Md_agg2[i,])/sum(as.numeric(Size_Md_agg2[i,]), na.rm = TRUE)*100
  
}

# match metadata
sample_name <- str_split_fixed(rownames(Size_Md_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Size_Md_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Size_Md_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Size_Md_agg2$sample_name<-rownames(Size_Md_agg2)


data_Size_Md<-dplyr::full_join(Size_Md_agg2, data2, by="sample_name")
colnames(data_Size_Md) <-c("NA","LTR","macro","micro","nano", "pico","sample_name","type")

#######Shape
Shape<-Traits[,-c(1:6,8:17)]
Shape_agg<-aggregate(.~Shape.Groups, Shape, sum)
Shape_agg2<-as.data.frame(t(as.matrix(Shape_agg[,c(2:1078)]))) 
names(Shape_agg2)<-levels(as.factor(Shape_agg$Shape))
Shape_agg2 <- Shape_agg2[which(rowSums(Shape_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Shape_agg2)) {
  print(i)
  Shape_agg2[i,] <- as.numeric(Shape_agg2[i,])/sum(as.numeric(Shape_agg2[i,]), na.rm = TRUE)*100
  
}
#########


# match metadata
sample_name <- str_split_fixed(rownames(Shape_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Shape_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Shape_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Shape_agg2$sample_name<-rownames(Shape_agg2)


data_Shape<-dplyr::full_join(Shape_agg2, data2, by="sample_name")
colnames(data_Shape) <-c("NA","cone","cylinder", "ellipsoid","LTR", "parallelepiped", "pyriform", "scoop", "sphere", "trapezoid", "variable","sample_name","type")


#######Motility
Motility<-Traits[,-c(1:7,9:17)]
Motility_agg<-aggregate(.~Motility, Motility, sum)
Motility_agg2<-as.data.frame(t(as.matrix(Motility_agg[,c(2:1078)]))) 
names(Motility_agg2)<-levels(as.factor(Motility_agg$Motility))
Motility_agg2 <- Motility_agg2[which(rowSums(Motility_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Motility_agg2)) {
  print(i)
  Motility_agg2[i,] <- as.numeric(Motility_agg2[i,])/sum(as.numeric(Motility_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata

sample_name <- str_split_fixed(rownames(Motility_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Motility_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Motility_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Motility_agg2$sample_name<-rownames(Motility_agg2)


data_Motility<-dplyr::full_join(Motility_agg2, data2, by="sample_name")
colnames(data_Motility) <-c("NA","attached","floater", "gliding", "LTR", "swimmer","sample_name","type")



#######Cover  
Cover<-Traits[,-c(1:8,10:17)]
Cover_agg<-aggregate(.~Cover, Cover, sum) # the first is for the variable, the second is the data
Cover_agg2<-as.data.frame(t(as.matrix(Cover_agg[,c(2:1078)]))) 
names(Cover_agg2)<-levels(as.factor(Cover_agg$Cover))
Cover_agg2 <- Cover_agg2[which(rowSums(Cover_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Cover_agg2)) {
  print(i)
  Cover_agg2[i,] <- as.numeric(Cover_agg2[i,])/sum(as.numeric(Cover_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata

sample_name <- str_split_fixed(rownames(Cover_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Cover_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Cover_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Cover_agg2$sample_name<-rownames(Cover_agg2)


data_Cover<-dplyr::full_join(Cover_agg2, data2, by="sample_name")
colnames(data_Cover) <-c("NA","Calcareous","LTR", "Naked", "Organic", "Siliceous","sample_name","type")


#######Colony  
Colony<-Traits[,-c(1:9,11:17)]
Colony_agg<-aggregate(.~Colony, Colony, sum) # the first is for the variable, the second is the data
Colony_agg2<-as.data.frame(t(as.matrix(Colony_agg[,c(2:1078)]))) 
names(Colony_agg2)<-levels(as.factor(Colony_agg$Colony))
Colony_agg2 <- Colony_agg2[which(rowSums(Colony_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Colony_agg2)) {
  print(i)
  Colony_agg2[i,] <- as.numeric(Colony_agg2[i,])/sum(as.numeric(Colony_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata

sample_name <- str_split_fixed(rownames(Colony_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Colony_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Colony_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Colony_agg2$sample_name<-rownames(Colony_agg2)
colnames(Colony_agg2) <-c("NA","No_Colony","Yes_Colony","LTR", "sample_name")


data_Colony<-dplyr::full_join(Colony_agg2, data2, by="sample_name")


#######Mucilage  
Mucilage<-Traits[,-c(1:10,12:17)]
Mucilage_agg<-aggregate(.~Mucilage, Mucilage, sum)
Mucilage_agg2<-as.data.frame(t(as.matrix(Mucilage_agg[,c(2:1078)]))) 
names(Mucilage_agg2)<-levels(as.factor(Mucilage_agg$Mucilage))
Mucilage_agg2 <- Mucilage_agg2[which(rowSums(Mucilage_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Mucilage_agg2)) {
  print(i)
  Mucilage_agg2[i,] <- as.numeric(Mucilage_agg2[i,])/sum(as.numeric(Mucilage_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata

sample_name <- str_split_fixed(rownames(Mucilage_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Mucilage_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Mucilage_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Mucilage_agg2$sample_name<-rownames(Mucilage_agg2)
colnames(Mucilage_agg2) <-c("NA","No_Mucilage","Yes_Mucilage","LTR", "sample_name")


data_Mucilage<-dplyr::full_join(Mucilage_agg2, data2, by="sample_name")


#######Flagella 
Flagella<-Traits[,-c(1:11,13:17)]
Flagella_agg<-aggregate(.~Flagella, Flagella, sum)
Flagella_agg2<-as.data.frame(t(as.matrix(Flagella_agg[,c(2:1078)]))) 
names(Flagella_agg2)<-levels(as.factor(Flagella_agg$Flagella))
Flagella_agg2 <- Flagella_agg2[which(rowSums(Flagella_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Flagella_agg2)) {
  print(i)
  Flagella_agg2[i,] <- as.numeric(Flagella_agg2[i,])/sum(as.numeric(Flagella_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata
sample_name <- str_split_fixed(rownames(Flagella_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Flagella_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Flagella_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Flagella_agg2$sample_name<-rownames(Flagella_agg2)
colnames(Flagella_agg2) <-c("NA","No_Flagella","Yes_Flagella","LTR","sample_name")


data_Flagella<-dplyr::full_join(Flagella_agg2, data2, by="sample_name")


#######Cilia 
Cilia<-Traits[,-c(1:12,14:17)]
Cilia_agg<-aggregate(.~Cilia, Cilia, sum)
Cilia_agg2<-as.data.frame(t(as.matrix(Cilia_agg[,c(2:1078)]))) 
names(Cilia_agg2)<-levels(as.factor(Cilia_agg$Cilia))
Cilia_agg2 <- Cilia_agg2[which(rowSums(Cilia_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Cilia_agg2)) {
  print(i)
  Cilia_agg2[i,] <- as.numeric(Cilia_agg2[i,])/sum(as.numeric(Cilia_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata
sample_name <- str_split_fixed(rownames(Cilia_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Cilia_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Cilia_agg2), metadata$sample_name),]


data2 <- metadata[,c(2,4)]
Cilia_agg2$sample_name<-rownames(Cilia_agg2)
colnames(Cilia_agg2) <-c("NA","No_Cilia","Yes_Cilia","LTR", "sample_name")


data_Cilia<-dplyr::full_join(Cilia_agg2, data2, by="sample_name")


#######Cyst 
Cyst<-Traits[,-c(1:13,15:17)]
Cyst_agg<-aggregate(.~Cyst, Cyst, sum)
Cyst_agg2<-as.data.frame(t(as.matrix(Cyst_agg[,c(2:1078)]))) 
names(Cyst_agg2)<-levels(as.factor(Cyst_agg$Cyst))
Cyst_agg2 <- Cyst_agg2[which(rowSums(Cyst_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Cyst_agg2)) {
  print(i)
  Cyst_agg2[i,] <- as.numeric(Cyst_agg2[i,])/sum(as.numeric(Cyst_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata
sample_name <- str_split_fixed(rownames(Cyst_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Cyst_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Cyst_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Cyst_agg2$sample_name<-rownames(Cyst_agg2)
colnames(Cyst_agg2) <-c("NA","No_Cyst","Yes_Cyst","LTR", "sample_name")


data_Cyst<-dplyr::full_join(Cyst_agg2, data2, by="sample_name")



#######Spore 
Spore<-Traits[,-c(1:14,16:17)]
Spore_agg<-aggregate(.~Spore, Spore, sum)
Spore_agg2<-as.data.frame(t(as.matrix(Spore_agg[,c(2:1078)]))) 
names(Spore_agg2)<-levels(as.factor(Spore_agg$Spore))
Spore_agg2 <- Spore_agg2[which(rowSums(Spore_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Spore_agg2)) {
  print(i)
  Spore_agg2[i,] <- as.numeric(Spore_agg2[i,])/sum(as.numeric(Spore_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata
sample_name <- str_split_fixed(rownames(Spore_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Spore_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Spore_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Spore_agg2$sample_name<-rownames(Spore_agg2)
colnames(Spore_agg2) <-c("NA","No_Spore","Yes_Spore","LTR", "sample_name")


data_Spore<-dplyr::full_join(Spore_agg2, data2, by="sample_name")


#########Cyst&Spore -> Resting Stage (SR)

Spore_Cyst<-Traits[,-c(1:16)]
Spore_Cyst_agg<-aggregate(.~Spore_Cyst, Spore_Cyst, sum)
Spore_Cyst_agg2<-as.data.frame(t(as.matrix(Spore_Cyst_agg[,c(2:1078)]))) 
names(Spore_Cyst_agg2)<-levels(as.factor(Spore_Cyst_agg$Spore))
Spore_Cyst_agg2 <- Spore_Cyst_agg2[which(rowSums(Spore_Cyst_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances (??? to be discussed)
for (i in 1:nrow(Spore_Cyst_agg2)) {
  print(i)
  Spore_Cyst_agg2[i,] <- as.numeric(Spore_Cyst_agg2[i,])/sum(as.numeric(Spore_Cyst_agg2[i,]), na.rm = TRUE)*100
  
}
#########

# match metadata
sample_name <- str_split_fixed(rownames(Spore_Cyst_agg2), "\\.",3)[,3] %>% substr(.,4,100)

# some samples are named with additional dates in the sequencing run
problem_names <- sample_name[which(is.na(match(sample_name, metadata$sample_name)))]
fixed_names <- str_split_fixed(problem_names, "\\.",3)[,3] %>% substr(.,4,100)
sample_name[which(is.na(match(sample_name, metadata$sample_name)))] <- fixed_names
rownames(Spore_Cyst_agg2) <- sample_name   

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Spore_Cyst_agg2), metadata$sample_name),]
data2 <- metadata[,c(2,4)]
Spore_Cyst_agg2$sample_name<-rownames(Spore_Cyst_agg2)
colnames(Spore_Cyst_agg2) <-c("NA","No_RS","Yes_RS","LTR", "sample_name")


data_Spore_Cyst<-dplyr::full_join(Spore_Cyst_agg2, data2, by="sample_name")




#######Figure 4c

####################
######## for PCoA traits

data_all<-dplyr::full_join(Trophy_agg2,Size_Md_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Motility_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Shape_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Cover_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Colony_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Spore_Cyst_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Mucilage_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Flagella_agg2, by="sample_name")
data_all<-dplyr::full_join(data_all, Cilia_agg2, by="sample_name")

######just removinf NA and LTR
#data_all <-data_all[,-c(1,5,10,13,15,16,21,25,27,31,38,40,44,45,48,49,52,53,56,57,60,61,64,65)]


####removing also the No for binary variables No_Cilia, No_Flagella, etc...
#data_all <-data_all[,-c(1,4,7,8,13,17,19,23,30,32,36,37,39:41,43:45,47:49,51:53,55:57,59)] ##including spore and cyst separated
data_all <-data_all[,-c(1,4,7,8,13,17,19,23,30,32,36,37,39:41,43:45,47:49,51:53,55)] ##including spore and cyst as resting stage


rownames(data_all)<-data_all$sample_name
#data_all <-data_all[,-c(7)]
data_all <-data_all[,-c(4)] ###deleting sample name



#######PCOA arrows
library(RVAideMemoire)

compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  trait_df = trait_df[, c(1:29)]
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  
  # Add values of covariances inside object
  given_pcoa$U <- U
  
  return(given_pcoa)
}


data_all <- data_all[which(rowSums(data_all) > 0), ] #######deletes rows summ 0
data_all <- data_all[,which(colSums(data_all) > 0) ] #######deletes columns summ 0

data_all_bray_dm<-vegan::vegdist(data_all,method="bray",diag=TRUE,upper=TRUE)


#data_all_bray_dm_fixed<-stepacross(data_all_bray_dm, path="shortest",toolong=0.75)
#data_all_bray_dm_fixed<-as.matrix(data_all_bray_dm_fixed)

data_all_bray_dm<-as.matrix(data_all_bray_dm)
data_all_pcoa <- pcoa(data_all_bray_dm)
summary(data_all_pcoa)


data_all_pcoa$values
pcoa_vectors<-data_all_pcoa$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$sample_name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL


head(pcoa_vectors)


dim(metadata)
dim(pcoa_vectors)
data<-dplyr::left_join(pcoa_vectors,metadata,by="sample_name") ####left_join(): includes all rows in x.



length(data$type)
vegan::adonis2(data_all_bray_dm ~ data$type)


######Call:
#vegan::adonis(formula = data_all_bray_dm ~ data$type) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#############Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#data$type    3    13.299  4.4330  39.986 0.10459  0.001 ***
#Residuals 1027   113.858  0.1109         0.89541           
#Total     1030   127.157                 1.00000             
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




## Plotting the PCoA
trait_pcoa_arrows = compute_arrows(data_all_pcoa, data_all)



# /!\ Be sure to change $vectors to a data.frame before putting it in ggplot2 /!\

vectors <- as.data.frame(trait_pcoa_arrows$vectors)
vectors$sample_name<-rownames(vectors)
rownames(vectors)<-NULL


head(vectors)
vectors<-dplyr::left_join(vectors,metadata,by="sample_name")
head(vectors)

####ARROWS!!!!
#arrows_df<-as.data.frame(trait_pcoa_arrows$U/10000) #this was the original, can be customized
#arrows_df<-as.data.frame(trait_pcoa_arrows$U) #####this is the data frame that will help customize and plot the arrows in the PcoA
arrows_df<-as.data.frame(trait_pcoa_arrows$U/250) 

barplot(data_all_pcoa$values$Relative_eig[1:10]) # plot the eigenvalues and interpret
PCOAaxes <- data_all_pcoa$vectors[,c(1,2)] ## Extract the plot scores from first two PCoA axes (if you need them):
# Now plot a bar plot of relative eigenvalues. This is the percentage variance explained by each axis
barplot((data_all_pcoa$values$Eigenvalues)/sum(data_all_pcoa$values$Eigenvalues)) 

# Calculate the percent of variance explained by first two axes
sum(as.vector((data_all_pcoa$values$Eigenvalues)/sum(data_all_pcoa$values$Eigenvalues))[1:2])
#[1] 0.6902
#or by each axes 
sum(as.vector((data_all_pcoa$values$Eigenvalues)/sum(data_all_pcoa$values$Eigenvalues))[1])
sum(as.vector((data_all_pcoa$values$Eigenvalues)/sum(data_all_pcoa$values$Eigenvalues))[2])
#[1] 0.4581398
#[1] 0.2320602



#####figure plot
p1_2 <- ggplot(NULL, aes(Axis.1, Axis.2)) +
  geom_point(data = vectors,size=2,aes(color=type)) +
  stat_ellipse(data = vectors, geom = "polygon", aes(color=type), alpha = 0, show.legend = FALSE, level = 0.95, size=1.1) +
  #geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,mapping=aes(xend=Axis.1,yend=Axis.2),arrow=arrow(length=unit(3,"mm")))
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.559,yend=0.4206),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.372,yend=-0.42),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.417,yend=-0.296),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.349,yend=-0.349),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.259,yend=-0.571),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.271,yend=-0.554),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.19,yend=-0.3),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.262,yend=0.508),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.211,yend=0.332),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.048,yend=0.312),arrow=arrow(length=unit(3,"mm")))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 25),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=17),legend.text=element_text(size=17), plot.title=element_text(size=30, face= "bold"))+
  scale_shape_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c(15,16,17,18))+
  scale_color_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  xlab("\nAxis 1 (45.8% of variance)")+
  ylab("Axis 2 (23.2% of variance)\n")+
  annotate("text", x=-0.78, y=0.8,label=paste("list(Adonis:~F['3,1030']==40", ", p==0.001", ", R^{2}==0.105)"),parse=TRUE,hjust=0.0,size=9.5)+
  annotate("text", x=-0.559,y=0.47, label="Si",size=8)+
  annotate("text",x=-0.372,y=-0.45, label="At",size=8)+
  annotate("text",x=0.57,y=-0.2957, label="Swim",size =8)+
  annotate("text",x=0.36,y=-0.37, label="Fla",size =8)+
  annotate("text",x=0.259,y=-0.62, label="Org",size=8)+
  annotate("text",x=0.35,y=-0.554, label="Nano",size =8)+
  annotate("text",x=0.19,y=-0.35, label="Col",size =8)+
  annotate("text",x=0.262,y=0.55, label="Ht",size =8)+
  annotate("text",x=0.211,y=0.39, label="Nk",size =8)+
  annotate("text",x=-0.048,y=0.37, label="Micro",size =8)




p1_2 +labs(title = "Eukaryotes PcoA (phenotypic traits)")


pdf("Euk_traits_PcoA_March23.pdf", width = 12, height = 8, pointsize=12)
p1_2 +labs(title = "Eukaryotes PCoA (phenotypic traits)")
dev.off()


###samples

data2$type
type <-as.factor(data2$type)
experimental <-subset(data2, type=="experimental") #590
environmental_water <-subset(data2, type=="environmental_water") #298
environmental_soil <-subset(data2, type=="environmental_soil") #75
environmental_aerosol <-subset(data2, type=="environmental_aerosol") #71





#########################
########Boxplot##########

#########FIGURE 5b

#######colony
data_Colony <- data_Colony[,-c(5)]
Colony <-melt(data = data_Colony, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Colony)

# specify the factor levels in the order you want
Colony$type <- factor(Colony$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))

#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

###### anova tukey Colony
res.aov.Colony <- aov(Yes_Colony ~ type, data = data_Colony)
summary(res.aov.Colony)

Tukey.Colony<-TukeyHSD(res.aov.Colony)
Tukey.Colony



Colony_plot <- ggplot(subset(Colony,variable=="Yes_Colony"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 35, label="AB",size =6)+
  annotate("text",x = 1.65, y = 30, label="B",size =6)+
  annotate("text",x = 2.65, y = 45, label="A",size =6)+
  annotate("text",x = 3.65, y = 20, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
ylab("Colony")

Colony_plot


pdf("boxplot_Colony_AXIS.pdf", width = 8, height = 4, pointsize=12)
Colony_plot
dev.off()


########cover organic
data_Cover <- data_Cover[,-c(7)]
Cover <-melt(data = data_Cover, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Cover)

# specify the factor levels in the order you want
Cover$type <- factor(Cover$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

#######anova tukey cover organic
res.aov.Cover <- aov(Organic ~ type, data = data_Cover)
summary(res.aov.Cover)

Tukey.Cover<-TukeyHSD(res.aov.Cover)
Tukey.Cover


Organic_plot <- ggplot(subset(Cover,variable=="Organic"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 50, label="B",size =6)+
  annotate("text",x = 1.65, y = 50, label="AB",size =6)+
  annotate("text",x = 2.65, y = 65, label="A",size =6)+
  annotate("text",x = 3.65, y = 50, label="AB",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Organic Cover")

Organic_plot


pdf("boxplot_Organic_AXIS.pdf", width = 8, height = 4, pointsize=12)
Organic_plot
dev.off()


######Resting stage
data_Spore_Cyst <- data_Spore_Cyst[,-c(5)]
Spore_Cyst <-melt(data = data_Spore_Cyst, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Spore_Cyst)

# specify the factor levels in the order you want
Spore_Cyst$type <- factor(Spore_Cyst$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova tukey Spore_Cyst
res.aov.Spore_Cyst <- aov(Yes_RS ~ type, data = data_Spore_Cyst)
summary(res.aov.Spore_Cyst)

Tukey.Spore_Cyst<-TukeyHSD(res.aov.Spore_Cyst)
Tukey.Spore_Cyst


Spore_Cyst_plot <- ggplot(subset(Spore_Cyst,variable=="Yes_RS"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 100, label="C",size =6)+
  annotate("text",x = 1.65, y = 100, label="B",size =6)+
  annotate("text",x = 2.65, y = 105, label="AB",size =6)+
  annotate("text",x = 3.65, y = 105, label="A",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Resting Stage")

Spore_Cyst_plot


pdf("boxplot_Spore_Cyst_AXIS.pdf", width = 8, height = 4, pointsize=12)
Spore_Cyst_plot
dev.off()

######Autotrophy
data_Trophy <- data_Trophy[,-c(6)]
Trophy <-melt(data = data_Trophy, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Trophy)

# specify the factor levels in the order you want
Trophy$type <- factor(Trophy$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova and tukey Trophy
res.aov.Trophy <- aov(A ~ type, data = data_Trophy)
summary(res.aov.Trophy)

Tukey.Trophy<-TukeyHSD(res.aov.Trophy)
Tukey.Trophy


Autotrophy_plot <- ggplot(subset(Trophy,variable=="A"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 70, label="C",size =6)+
  annotate("text",x = 1.65, y = 100, label="A",size =6)+
  annotate("text",x = 2.65, y = 85, label="B",size =6)+
  annotate("text",x = 3.65, y = 100, label="A",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Autotrophy")

Autotrophy_plot


pdf("boxplot_Autotrophy_AXIS.pdf", width = 8, height = 4, pointsize=12)
Autotrophy_plot
dev.off()








############################
#######additional graphs for supplementary figures####

######Spore
data_Spore <- data_Spore[,-c(5)]
Spore <-melt(data = data_Spore, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Spore)

# specify the factor levels in the order you want
Spore$type <- factor(Spore$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova tukey Spore
res.aov.Spore <- aov(Yes_Spore ~ type, data = data_Spore)
summary(res.aov.Spore)

Tukey.Spore<-TukeyHSD(res.aov.Spore)
Tukey.Spore


Spore_plot <- ggplot(subset(Spore,variable=="Yes_Spore"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 85, label="C",size =6)+
  annotate("text",x = 1.65, y = 90, label="B",size =6)+
  annotate("text",x = 2.65, y = 105, label="A",size =6)+
  annotate("text",x = 3.65, y = 105, label="A",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Spore")

Spore_plot


pdf("boxplot_Spore_AXIS.pdf", width = 8, height = 4, pointsize=12)
Spore_plot
dev.off()


######Cilia
data_Cilia <- data_Cilia[,-c(5)]
Cilia <-melt(data = data_Cilia, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cilia or data_Size_Md
str(Cilia)

# specify the factor levels in the order you want
Cilia$type <- factor(Cilia$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova and tukey Cilia
res.aov.Cilia <- aov(Yes_Cilia ~ type, data = data_Cilia)
summary(res.aov.Cilia)

Tukey.Cilia<-TukeyHSD(res.aov.Cilia)
Tukey.Cilia


Cilia_plot <- ggplot(subset(Cilia,variable=="Yes_Cilia"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 80, label="A",size =6)+
  annotate("text",x = 1.65, y = 15, label="B",size =6)+
  annotate("text",x = 2.65, y = 15, label="B",size =6)+
  annotate("text",x = 3.65, y = 10, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
ylab("Cilia")

Cilia_plot


pdf("boxplot_Cilia_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Cilia_plot
dev.off()


#######Cyst
data_Cyst <- data_Cyst[,-c(5)]
Cyst <-melt(data = data_Cyst, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Cyst or data_Size_Md
str(Cyst)

# specify the factor levels in the order you want
Cyst$type <- factor(Cyst$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova and tukey Cyst
res.aov.Cyst <- aov(Yes_Cyst ~ type, data = data_Cyst)
summary(res.aov.Cyst)

Tukey.Cyst<-TukeyHSD(res.aov.Cyst)
Tukey.Cyst


Cyst_plot <- ggplot(subset(Cyst,variable=="Yes_Cyst"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 55, label="A",size =6)+
  annotate("text",x = 1.65, y = 30, label="BC",size =6)+
  annotate("text",x = 2.65, y = 25, label="C",size =6)+
  annotate("text",x = 3.65, y = 35, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
ylab("Cyst")

Cyst_plot


pdf("boxplot_Cyst_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cyst or barplot_Size_Md
Cyst_plot
dev.off()


#########mucilage
data_Mucilage <- data_Mucilage[,-c(5)]
Mucilage <-melt(data = data_Mucilage, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Mucilage or data_Size_Md
str(Mucilage)

# specify the factor levels in the order you want
Mucilage$type <- factor(Mucilage$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank


###### ANova tukey Mucilage
res.aov.Mucilage <- aov(Yes_Mucilage ~ type, data = data_Mucilage)
summary(res.aov.Mucilage)

Tukey.Mucilage<-TukeyHSD(res.aov.Mucilage)
Tukey.Mucilage


Mucilage_plot <- ggplot(subset(Mucilage,variable=="Yes_Mucilage"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 10, label="C",size =6)+
  annotate("text",x = 1.65, y = 35, label="B",size =6)+
  annotate("text",x = 2.65, y = 30, label="C",size =6)+
  annotate("text",x = 3.65, y = 50, label="A",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"), values=c("chocolate4","#56B4E9","#E69F00","#009E73"))+
ylab("Mucilage")

Mucilage_plot


pdf("boxplot_Mucilage_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Mucilage or barplot_Size_Md
Mucilage_plot
dev.off()

######flagella
data_Flagella <- data_Flagella[,-c(5)]
Flagella <-melt(data = data_Flagella, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Flagella or data_Size_Md
str(Flagella)

# specify the factor levels in the order you want
Flagella$type <- factor(Flagella$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Flagella
res.aov.Flagella <- aov(Yes_Flagella ~ type, data = data_Flagella)
summary(res.aov.Flagella)

#######no significant statistical differences for flagella

Flagella_plot <- ggplot(subset(Flagella,variable=="Yes_Flagella"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
 ylab("Flagella")

Flagella_plot


pdf("boxplot_Flagella_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Flagella or barplot_Size_Md
Flagella_plot
dev.off()


########cover others

#######cover calcareous anova tukey
res.aov.Cover <- aov(Calcareous ~ type, data = data_Cover)
summary(res.aov.Cover)

Tukey.Cover<-TukeyHSD(res.aov.Cover)
Tukey.Cover

Calcareous_plot <- ggplot(subset(Cover,variable=="Calcareous"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 5, label="AB",size =6)+
  annotate("text",x = 1.65, y = 5, label="A",size =6)+
  annotate("text",x = 2.65, y = 5, label="AB",size =6)+
  annotate("text",x = 3.65, y = 5, label="B",size =6)+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Calcareous Cover")

Calcareous_plot


pdf("boxplot_Calcareous_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Calcareous_plot
dev.off()



#######naked anova tukey
res.aov.Cover <- aov(Naked ~ type, data = data_Cover)
summary(res.aov.Cover)

Tukey.Cover<-TukeyHSD(res.aov.Cover)
Tukey.Cover

Naked_plot <- ggplot(subset(Cover,variable=="Naked"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 85, label="A",size =6)+
  annotate("text",x = 1.65, y = 25, label="B",size =6)+
  annotate("text",x = 2.65, y = 25, label="B",size =6)+
  annotate("text",x = 3.65, y = 25, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Naked")

Naked_plot


pdf("boxplot_Naked_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Naked_plot
dev.off()

###########others cover
######Anova tukey Siliceous
res.aov.Cover <- aov(Siliceous ~ type, data = data_Cover)
summary(res.aov.Cover)

Tukey.Cover<-TukeyHSD(res.aov.Cover)
Tukey.Cover

Siliceous_plot <- ggplot(subset(Cover,variable=="Siliceous"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 10, label="D",size =6)+
  annotate("text",x = 1.65, y = 75, label="B",size =6)+
  annotate("text",x = 2.65, y = 50, label="C",size =6)+
  annotate("text",x = 3.65, y = 95, label="A",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Siliceous Cover")

Siliceous_plot


pdf("boxplot_Siliceous_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Siliceous_plot
dev.off()

#####Motility others
data_Motility <- data_Motility[,-c(7)]

Motility <-melt(data = data_Motility, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Motility or data_Size_Md
str(Motility)

# specify the factor levels in the order you want
Motility$type <- factor(Motility$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

#######Anova tukey swimmer
res.aov.Motility <- aov(swimmer ~ type, data = data_Motility)
summary(res.aov.Motility)

Tukey.Motility<-TukeyHSD(res.aov.Motility)
Tukey.Motility


Swimmer_plot <- ggplot(subset(Motility,variable=="swimmer"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 90, label="A",size =6)+
  annotate("text",x = 1.65, y = 65, label="B",size =6)+
  annotate("text",x = 2.65, y = 65, label="BC",size =6)+
  annotate("text",x = 3.65, y = 65, label="C",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Swimmer")

Swimmer_plot

pdf("boxplot_Swimmer_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Motility or barplot_Size_Md
Swimmer_plot
dev.off()

##########Anova tukey gliding
res.aov.Motility <- aov(gliding ~ type, data = data_Motility)
summary(res.aov.Motility)

Tukey.Motility<-TukeyHSD(res.aov.Motility)
Tukey.Motility


Gliding_plot <- ggplot(subset(Motility,variable=="gliding"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 20, label="A",size =6)+
  annotate("text",x = 1.65, y = 5, label="C",size =6)+
  annotate("text",x = 2.65, y = 5, label="BC",size =6)+
  annotate("text",x = 3.65, y = 10, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Gliding")

Gliding_plot


pdf("boxplot_Gliding_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Motility or barplot_Size_Md
Gliding_plot
dev.off()

#######anova tukey floater
res.aov.Motility <- aov(floater ~ type, data = data_Motility)
summary(res.aov.Motility)

Tukey.Motility<-TukeyHSD(res.aov.Motility)
Tukey.Motility

Floater_plot <- ggplot(subset(Motility,variable=="floater"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 10, label="C",size =6)+
  annotate("text",x = 1.65, y = 20, label="B",size =6)+
  annotate("text",x = 2.65, y = 10, label="C",size =6)+
  annotate("text",x = 3.65, y = 35, label="A",size =6)+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Floater")

Floater_plot


pdf("boxplot_Floater_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Motility or barplot_Size_Md
Floater_plot
dev.off()


######Anova tukey attached
res.aov.Motility <- aov(attached ~ type, data = data_Motility)
summary(res.aov.Motility)

Tukey.Motility<-TukeyHSD(res.aov.Motility)
Tukey.Motility

Attached_plot <- ggplot(subset(Motility,variable=="attached"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 5, label="B",size =6)+
  annotate("text",x = 1.65, y = 5, label="B",size =6)+
  annotate("text",x = 2.65, y = 5, label="A",size =6)+
  annotate("text",x = 3.65, y = 5, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Attached")

Attached_plot


pdf("boxplot_Attached_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Motility or barplot_Size_Md
Attached_plot
dev.off()

#########Shape
data_Shape <- data_Shape[,-c(12)]
Shape <-melt(data = data_Shape, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Shape or data_Size_Md
str(Shape)

# specify the factor levels in the order you want
Shape$type <- factor(Shape$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

####Better way to include letters
#####Anova tukey cone
res.aov.Shape <- aov(cone ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#####Anova tukey cylinder
res.aov.Shape <- aov(cylinder ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#####Anova tukey ellipsoid
res.aov.Shape <- aov(ellipsoid ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#######Anova tukey parallelepiped
res.aov.Shape <- aov(parallelepiped ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#######Anova tukey pyriform
res.aov.Shape <- aov(pyriform ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#######Anova tukey scoop
res.aov.Shape <- aov(scoop ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#######Anova tukey sphere
res.aov.Shape <- aov(sphere ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

#######Anova tukey trapezoid
res.aov.Shape <- aov(trapezoid ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape


######Anova tukey variable
res.aov.Shape <- aov(variable ~ type, data = data_Shape)
summary(res.aov.Shape)

Tukey.Shape<-TukeyHSD(res.aov.Shape)
Tukey.Shape

######no significant differences in "variable' shape



ann_dat_text<-data.frame(variable=c("cone","cylinder","ellipsoid","parallelepiped","pyriform","scoop", "sphere","trapezoid", "variable"),
                         x=c(2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5),
                         y=c(-3,-3,-3,-3,-3,-3,-3,-3,-3),
                         label=c("B A B B","C B C A","A CB C B","A A B B","A C B C","B A B B","B AB B A","B A B B",""))


Shape_plot <- ggplot(subset(Shape,variable!="NA" & variable!="LTR"), aes(x = type, y= value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape=1, colour="grey27")+
  facet_grid(cols = vars(variable))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 10, angle=45, hjust=1),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=1, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial","environmental_aerosol"="Aerial","environmental_water"="Aquatic","experimental"="Tank"))+
  #scale_fill_manual(name="Environment",labels=c("environmental_soil"="Soil","environmental_aerosol"="Aerosol","environmental_water"="Water","experimental"="Tank"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  ylab("Shape")+ 
  geom_text(
    data = ann_dat_text,
    mapping= aes(x=x, y=y, label=label), size=5)

Shape_plot


pdf("boxplot_Shape_AXIS.pdf", width = 9, height = 5, pointsize=12) #######change name of graph, for example barplot_Shape or barplot_Size_Md
Shape_plot
dev.off()

########Size_Md (Median)
data_Size_Md <- data_Size_Md[,-c(7)]

#######nano no statistical differences
Size_Md <-melt(data = data_Size_Md, id.vars = "type")  #######This will reorder the table, also remember to change name of data, for example data_Size_Md or data_Size_Md
str(Size_Md)

# specify the factor levels in the order you want
Size_Md$type <- factor(Size_Md$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank

######Anova Tukey pico
res.aov.Size_Md <- aov(pico ~ type, data = data_Size_Md)
summary(res.aov.Size_Md)
#######pico size, almost all values zero and no statistical differences

Pico_plot <- ggplot(subset(Size_Md,variable=="pico"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Pico Size")

Pico_plot

png("boxplot_Pico_Nov22.png", width = 8, height = 4, units = 'in', pointsize=12,res = 150) #######change name of graph, for example barplot_Size_Md or barplot_Size_Md
Pico_plot
dev.off()

######Anova Tukey nano
res.aov.Size_Md <- aov(nano ~ type, data = data_Size_Md)
summary(res.aov.Size_Md)
#######nano size no statistical differences

Nano_plot <- ggplot(subset(Size_Md,variable=="nano"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Nano Size")

Nano_plot

pdf("boxplot_Nano_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Size_Md or barplot_Size_Md
Nano_plot
dev.off()

######Anova Tukey micro
res.aov.Size_Md <- aov(micro ~ type, data = data_Size_Md)
summary(res.aov.Size_Md)

Tukey.Size_Md<-TukeyHSD(res.aov.Size_Md)
Tukey.Size_Md



Micro_plot <- ggplot(subset(Size_Md,variable=="micro"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 75, label="A",size =6)+
  annotate("text",x = 1.65, y = 60, label="A",size =6)+
  annotate("text",x = 2.65, y = 50, label="B",size =6)+
  annotate("text",x = 3.65, y = 50, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Micro Size")

Micro_plot


pdf("boxplot_Micro_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Size_Md or barplot_Size_Md
Micro_plot
dev.off()

########Anova tukey macro
res.aov.Size_Md <- aov(macro ~ type, data = data_Size_Md)
summary(res.aov.Size_Md)

Tukey.Size_Md<-TukeyHSD(res.aov.Size_Md)
Tukey.Size_Md

Macro_plot <- ggplot(subset(Size_Md,variable=="macro"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 5, label="B",size =6)+
  annotate("text",x = 1.65, y = 5, label="A",size =6)+
  annotate("text",x = 2.65, y = 5, label="B",size =6)+
  annotate("text",x = 3.65, y = 5, label="B",size =6)+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Macro Size")

Macro_plot


pdf("boxplot_Macro_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Size_Md or barplot_Size_Md
Macro_plot
dev.off()


#######Trophy others

######anova tukey heterotrophy
res.aov.Trophy <- aov(H ~ type, data = data_Trophy)
summary(res.aov.Trophy)

Tukey.Trophy<-TukeyHSD(res.aov.Trophy)
Tukey.Trophy


Heterotrophy_plot <- ggplot(subset(Trophy,variable=="H"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y = 105, label="A",size =6)+
  annotate("text",x = 1.65, y = 25, label="D",size =6)+
  annotate("text",x = 2.65, y = 80, label="B",size =6)+
  annotate("text",x = 3.65, y = 35, label="C",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Heterotrophy")

Heterotrophy_plot


pdf("boxplot_Heterotrophy_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Heterotrophy_plot
dev.off()

######anova tukey Mixotrophy
res.aov.Trophy <- aov(M ~ type, data = data_Trophy)
summary(res.aov.Trophy)

Tukey.Trophy<-TukeyHSD(res.aov.Trophy)
Tukey.Trophy



Mixotrophy_plot <- ggplot(subset(Trophy,variable=="M"), aes(x = type, y= value, fill=type))+
  geom_boxplot(outlier.shape = NA)+ 
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA), legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  annotate("text",x = 0.65, y =2, label="B",size =6)+
  annotate("text",x = 1.65, y =5, label="A",size =6)+
  annotate("text",x = 2.65, y =2, label="B",size =6)+
  annotate("text",x = 3.65, y =2, label="B",size =6)+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))+
  ylab("Mixotrophy")

Mixotrophy_plot


pdf("boxplot_Mixotrophy_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Mixotrophy_plot
dev.off()



#######GRID
# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load all the packages you are missing!
pacman::p_load('knitr', 'lme4', 'lmerTest', 'tidyverse', 'magrittr', 'effects', 'plyr', 'dplyr', 'plotrix', 'car',"gridExtra", "cowplot", "tools", "mgcv", "gratia", "MASS", "stats", "tidymv", "sjstats", "coin", "emmeans")


extract.legend <- get_legend(
  # create some space to the left of the legend
  Mixotrophy_plot + theme(legend.box.margin = margin(0, 0, 0, 10)))


EukTrait_sup <-plot_grid(
  Nano_plot + theme(legend.position = "none"), 
  Micro_plot + theme(legend.position = "none"), 
  Macro_plot+ theme(legend.position = "none"), 
  Calcareous_plot + theme(legend.position = "none"), 
  Naked_plot + theme(legend.position = "none"), 
  Siliceous_plot + theme(legend.position = "none"), 
  Mucilage_plot + theme(legend.position = "none"), 
  Cilia_plot + theme(legend.position = "none"), 
  Flagella_plot + theme(legend.position = "none"), 
  Swimmer_plot + theme(legend.position = "none"), 
  Gliding_plot + theme(legend.position = "none"), 
  Floater_plot + theme(legend.position = "none"), 
  Heterotrophy_plot + theme(legend.position = "none"), 
  Mixotrophy_plot + theme(legend.position = "none"), 
  Attached_plot + theme(legend.position = "none"), 
  Spore_plot + theme(legend.position = "none"), 
  Cyst_plot + theme(legend.position = "none"), extract.legend,
  rel_widths = c(8,8,8), ncol=3, nrow = 6)
ggsave("EukTrait_sup_AXIS.pdf", height=18, width=15)

EukTrait_sup





