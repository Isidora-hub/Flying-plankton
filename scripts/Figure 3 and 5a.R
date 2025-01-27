
######CLEANED VERSION

#####Bacteria Taxonomy and Traits Analysis

library(phyloseq)
library(agricolae)
library(ggplot2)
library(grid)
library(ade4)
library(pheatmap)
library(tidyr)
library(rlang)
library(fmsb)
library(vegan)
library(ape)
library(devtools)
library(RVAideMemoire)
library(dplyr)
library(compositions)

# Loading data
setwd("~/Data")
load("clean_data_2021.Rdata")
x<-taxonomy$Hash
y<-colnames(asv_table)

# Building phyloseq object
mapping_phyloseq<-sample_data(metadata)
mapping_phyloseq$Samples<-mapping_phyloseq$sample_name
mapping_phyloseq$sample_name<-NULL
rownames(mapping_phyloseq)<-mapping_phyloseq$Samples
asv_table<-t(as.matrix(asv_table))
OTU_Table_phyloseq<-otu_table(asv_table, taxa_are_rows = TRUE)

rownames(taxonomy)<-taxonomy$Hash
taxonomy$Hash<-NULL
taxonomy$phylogeny<-NULL
Tax_Table_phyloseq<-as.matrix(taxonomy)
Tax_Table_phyloseq<-tax_table(Tax_Table_phyloseq)

my_phyloseq<-merge_phyloseq(OTU_Table_phyloseq,Tax_Table_phyloseq,mapping_phyloseq)

#Removing non-bacteria, bacteria unclassified, and Chloroplast

OTUs_to_Keep<-read.csv("otu_list_bacteria_only_no_chloroplast_no_bacteria_unclassified.csv",head=F)
my_list<-OTUs_to_Keep$V1
my_list<-my_list[-1]
test <- subset(otu_table(my_phyloseq), rownames(otu_table(my_phyloseq)) %in% my_list)
my_bacterial_phyloseq <- merge_phyloseq(test, tax_table(my_phyloseq), sample_data(my_phyloseq))

#Removing doubletons
bacterial_phyloseq_doubltons_removed <- prune_taxa(taxa_sums(my_bacterial_phyloseq) > 2, my_bacterial_phyloseq)

#Remove samples with fewer than 5000 reads
bacterial_phyloseq_doubltons_removed = prune_samples(sample_sums(bacterial_phyloseq_doubltons_removed)>=5000, bacterial_phyloseq_doubltons_removed)
#Rarefying samples to an even depth of 5032
rare_bacterial_phyloseq<-rarefy_even_depth(bacterial_phyloseq_doubltons_removed,sample.size = min(sample_sums(bacterial_phyloseq_doubltons_removed)),rngseed = 771,replace=FALSE,trimOTUs = TRUE,verbose=TRUE)


output<-as.data.frame(otu_table(rare_bacterial_phyloseq))
#write.csv(output,"otu_table_from_5krarefied_bacterial_phyloseq.csv")

output2<-as.data.frame(sample_data(rare_bacterial_phyloseq))
#write.csv(output2,"sample_table_from_5krarefied_bacterial_phyloseq.csv")

output3<-as.data.frame(tax_table(rare_bacterial_phyloseq))
#write.csv(output3,"taxonomy_table_from_5krarefied_bacterial_phyloseq.csv")

#######here starts the pipeline for RADARCHART
#####create a column with hash names
output$Hash<-rownames(output)
output3$Hash<-rownames(output3)

Prok<-dplyr::full_join(output, output3, by="Hash")

Prok[,c(992:998)]<-lapply(Prok[,c(992:998)], as.factor)

######remember to enter row names
rownames(Prok) <- Prok$Hash

####Tax_2 (or Supergroup, B)
Prok_super<-Prok[,-c(991:992,994:998)]
Prok_super_agg<-aggregate(.~B, Prok_super, sum)
Prok_super_agg2<-as.data.frame(t(as.matrix(Prok_super_agg[,c(2:991)]))) 
names(Prok_super_agg2)<-levels(as.factor(Prok_super_agg$B))
Prok_super_agg2 <- Prok_super_agg2[which(rowSums(Prok_super_agg2) > 0), ] #######deletes rows summ 0

#######transform to relative abundances
for (i in 1:nrow(Prok_super_agg2)) {
  print(i)
  Prok_super_agg2[i,] <- as.numeric(Prok_super_agg2[i,])/sum(as.numeric(Prok_super_agg2[i,]), na.rm = TRUE)*100
  
}
#########

#####create a column with sample_name names
Prok_super_agg2$sample_name<-rownames(Prok_super_agg2)

#reorder metadata table to align with asv table
metadata <- metadata[match(rownames(Prok_super_agg2), metadata$sample_name),]

#####keep only sample_name and type on metadata
metadata2 <- metadata[,c(2,4)]


data_super<-dplyr::full_join(Prok_super_agg2, metadata2, by="sample_name")
data_super[,c(1:34)]<-lapply(data_super[,c(1:34)], as.numeric)

data_super <- data_super[,-c(35)]

data_super_avg<-aggregate(.~type, data=data_super, mean)

#keep only taxa above 10%

data_super_avg10<-data_super_avg[,c(1,3,6,15,20,30)]


data10 <- as.data.frame(t(data_super_avg10)) #changes columns to rows, here we need MBFG on rows and sample types in columns
colnames(data10) <-c("Aerosol", "Soil","Water","Experimental" )
data10 <-data10[-c(1),]


######radarchart

# Figure 3a
# Library
library(fmsb)

data10 <- rbind(rep(60,4) , rep(0,4) , data10) # To use the fmsb package, add 2 lines to the dataframe: the max and min of each variable to show on the plot

data10[,c(1:4)]<-lapply(data10[,c(1:4)], as.numeric)#spaces apear after transpose, this fix the problem, first transform to numeric
data10 <-as.data.frame(data10) #then back to dataframe

radarchart(data10)

# Customize with legend
# Color vector
colors <- c("violet", "black", "#E69F00", "#009E73", "#0072B2")


# PDF plot with default options:
pdf("Radar_Bact_March23.pdf", width = 14, height = 8, pointsize=12)
radarchart(data10  , axistype=1 , 
           #custom polygon
           pcol=colors , plwd=4 , plty=1,
           #custom the grid
           cglcol="black", cglty=1, axislabcol="black", cglwd=0.8,caxislabels=c("0%","15%","30%","45%","60%"),
           #custom labels
           vlcex=2.4, vlabels=c("Aerial dispersers", "Terrestrial source                    ","Aquatic source","                 Tank colonizers" )) ####spaces help to center side labels in radar plot

# Add a legend
legend(x=0.55, y=1.2, legend = rownames(data10[-c(1,2),]), bty = "n", pch=20 , col=colors , text.col = "black", cex=1.8, pt.cex=3)
dev.off()


########### Bray-Curtis Taxonomic Ordinations

data<-as.data.frame(t(otu_table(rare_bacterial_phyloseq)))
data_dm<-vegan::vegdist(data,method="bray",diag=TRUE,upper=TRUE)
data_distancematrix<-as.matrix(data_dm)
data_distancematrix_fixed<-stepacross(data_distancematrix, path="shortest",toolong=0.95)
result<-pcoa(data_distancematrix_fixed)
compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  trait_df = trait_df[, c(1:18990)]#18990 ASVs
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
taxa_pcoa_arrows = compute_arrows(result, data)

pcoa_vectors<-result$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$Name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL

metadata<-as(sample_data(rare_bacterial_phyloseq), "matrix")
metadata<-as.data.frame(metadata)
metadata$Name<-rownames(metadata)
rownames(metadata)<-NULL
data<-dplyr::full_join(pcoa_vectors,metadata,by="Name")
data_dm<-as.matrix(data_dm)
vegan::adonis2(data_dm ~ data$type)




vectors <- as.data.frame(taxa_pcoa_arrows$vectors)
vectors$Name<-rownames(vectors)
rownames(vectors)<-NULL
vectors<-dplyr::full_join(vectors,metadata,by="Name")
arrows_df<-as.data.frame(taxa_pcoa_arrows$U/100)
arrows_df<-arrows_df[,1:2]
arrows_df$variable<-rownames(arrows_df)
arrows_df$sum<-abs(arrows_df$Axis.1)+abs(arrows_df$Axis.2)
rownames(arrows_df[order(arrows_df$sum, decreasing = T),])[1:10]
arrows_df<-(arrows_df[order(arrows_df$sum, decreasing = T),])
arrows_df<-arrows_df[1:10,]
#write.csv(arrows_df,"pcoa taxa arrows phylogengy test.csv")
labels<-read.csv("pcoa taxa arrows phylogeny.csv")
arrows_df$Hash<-rownames(arrows_df)
rownames(arrows_df)<-NULL
arrows_df<-dplyr::full_join(arrows_df,labels,by="Hash")

########get values of variance explained by each of the axis
sum(as.vector((taxa_pcoa_arrows$values$Eigenvalues)/sum(taxa_pcoa_arrows$values$Eigenvalues))[1])
sum(as.vector((taxa_pcoa_arrows$values$Eigenvalues)/sum(taxa_pcoa_arrows$values$Eigenvalues))[2])


# Figure 3b
Prok_tax_plot <- ggplot(NULL, aes(Axis.1, Axis.2)) + 
  geom_point(data = vectors,size=2,aes(color=type)) +
  stat_ellipse(data = vectors, geom = "polygon", aes(color=type), alpha = 0, show.legend = FALSE, level = 0.95, size=1.1) +
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=Axis.1,yend=Axis.2),arrow=arrow(length=unit(3,"mm")))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 25),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=17),legend.text=element_text(size=17), plot.title=element_text(size=30, face="bold"))+
  scale_shape_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c(15,16,17,18))+
  scale_color_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  xlab("\nAxis 1 (42.4% of variance)")+
  ylab("Axis 2 (25.2% of variance)\n")+
  geom_label(data = arrows_df, aes(label = Letter),size=10)+
  annotate("text", x=-4.0, y=2.65,label=paste("list(Adonis:~F['3,989']==64.9", ", p==0.001", ", R^{2}==0.165)"),parse=TRUE,hjust=0.0,size=9.5)

Prok_tax_plot +labs(title = "Bacterial PCoA (taxonomy)")

pdf("Prok_tax_PCoA_March23.pdf", width = 12, height = 8, pointsize=12)
Prok_tax_plot +labs(title = "Bacterial PCoA (taxonomy)")
dev.off()




######### Madin PCoA
####Figure 3c


compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  trait_df = trait_df[, c(1:18)] ####18 traits from Madin
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

data2<-read.csv("Madin_Pipeline_final_results_simple.csv",header=T,row.names=1)

data2$ID<-NULL
data2<-t(data2)
data<-cbind(data2)

data_dm<-vegan::vegdist(data,method="bray",diag=TRUE,upper=TRUE)
data_distancematrix<-as.matrix(data_dm)
data_distancematrix<-stepacross(data_distancematrix, path="shortest", toolong=0.75)
result<-pcoa(data_distancematrix)
trait_pcoa_arrows = compute_arrows(result, data)
trait_pcoa_arrows$U[,1]

pcoa_vectors<-result$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$Name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL

metadata<-as(sample_data(rare_bacterial_phyloseq), "matrix")
metadata<-as.data.frame(metadata)
metadata$Name<-rownames(metadata)
data<-dplyr::full_join(pcoa_vectors,metadata,by="Name")
data<-data[!is.na(data$Axis.1),]
data_distancematrix<-as.matrix(data_distancematrix)
vegan::adonis2(data_distancematrix ~ data$type)

# /!\ Be sure to change $vectors to a data.frame before putting it in ggplot2 /!\
vectors <- as.data.frame(trait_pcoa_arrows$vectors)
vectors$Name<-rownames(vectors)
rownames(vectors)<-NULL
vectors<-dplyr::full_join(vectors,metadata,by="Name")
arrows_df<-as.data.frame(trait_pcoa_arrows$U/1000)
arrows_df<-arrows_df[,1:2]
arrows_df$variable<-rownames(arrows_df)
#write.csv(arrows_df,"arrows_df_madin_revised.csv")

########get values of variance explained by each of the axis
sum(as.vector((trait_pcoa_arrows$values$Eigenvalues)/sum(trait_pcoa_arrows$values$Eigenvalues))[1])
sum(as.vector((trait_pcoa_arrows$values$Eigenvalues)/sum(trait_pcoa_arrows$values$Eigenvalues))[2])


#NEED TO PICK YOUR ARROWS AGAIN!!!
arrows_df<-read.csv("arrows_df_madin_revised_top6.csv",check.names=F)

rownames(arrows_df)<-arrows_df[,1]
arrows_df[,1]<-NULL

madin_plot<-ggplot(NULL, aes(Axis.1, Axis.2)) +
  geom_point(data = vectors,size=2,aes(color=type)) +
  stat_ellipse(data = vectors, geom = "polygon", aes(color=type), alpha = 0, show.legend = FALSE, level = 0.95, size=1.1) +
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=Axis.1,yend=Axis.2),arrow=arrow(length=unit(3,"mm")))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 25),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=17),legend.text=element_text(size=17), plot.title=element_text(size=30, face="bold"))+
  scale_shape_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c(15,16,17,18))+
  scale_color_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  xlab("\nAxis 1 (70.3% of variance)")+
  ylab("Axis 2 (8.0% of variance)\n")+
  annotate("text", x=-1.6, y=0.4,label=paste("list(Adonis:~F['3,989']==50.7", ", p==0.001", ", R^{2}==0.134)"),parse=TRUE,hjust=0.0,size=9.5)+
  annotate("text",label="Cell length",x=0.1,y=0.3,size=9)+
  #annotate("text",label= "Cell length (L)",x=-0.15012631,y=0.08570310,size=6)+
  #annotate("text",label="Doubling time",x=-0.12580879,y=-0.1078773,size=6)+
  annotate("text",label="Rod shape",x=0.376622,y=-0.22,size=9)+
  annotate("text",label="Cell dia",x=0.46,y=-0.11,size=9)+
  annotate("text",label="Mesoph",x=0,y=-0.18,size=9)+
  #annotate("text",label="Cell Diam (L)",x=0.42639730,y=-0.03645539,size=6)+
  annotate("text",label="Salt tol",x=0.400435420,y=-0.18,size=9)+
  annotate("text",label="Spor",x=0,y=0.2,size=9)
#annotate("text",label="Gram Positive",x=0.32508096,y=-0.0592753,size=6)

madin_plot

pdf("Madin_PcoA_March23.pdf", width = 12, height = 8, pointsize=12)
madin_plot +labs(title = "Bacterial PCoA (Madin et al. 2020)")
dev.off()




#### PCoA functional faprotax

#######Figure 3d

compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  trait_df = trait_df[, c(1:91)] #91 traits from faprotax
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

#removed the functions with zero occurrence
data<-read.csv("FAPROTAX_func_table.csv",header=T,row.names = 1)
data$ID<-NULL
data<-t(data)
data<-cbind(data)

data_dm<-vegan::vegdist(data,method="bray",diag=TRUE,upper=TRUE)
data_distancematrix<-as.matrix(data_dm)
data_distancematrix<-stepacross(data_distancematrix, path="shortest", toolong=0.75)
result<-pcoa(data_distancematrix)
trait_pcoa_arrows = compute_arrows(result, data)
trait_pcoa_arrows$U[,1]

pcoa_vectors<-result$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$Name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL

metadata<-as(sample_data(rare_bacterial_phyloseq), "matrix")
metadata<-as.data.frame(metadata)
metadata$Name<-rownames(metadata)
data<-dplyr::full_join(pcoa_vectors,metadata,by="Name")
data<-data[!is.na(data$Axis.1),]
data_distancematrix<-as.matrix(data_distancematrix)
vegan::adonis2(data_distancematrix ~ data$type)

# /!\ Be sure to change $vectors to a data.frame before putting it in ggplot2 /!\
vectors <- as.data.frame(trait_pcoa_arrows$vectors)
vectors$Name<-rownames(vectors)
rownames(vectors)<-NULL
vectors<-dplyr::full_join(vectors,metadata,by="Name")
arrows_df<-as.data.frame(trait_pcoa_arrows$U/1000)
arrows_df<-arrows_df[,1:2]
arrows_df$variable<-rownames(arrows_df)
#write.csv(arrows_df,"arrows_df_faprotax_with_cyanos.csv")

#NEED TO PICK YOUR ARROWS AGAIN!!!
#divided aerobic_chemoheterotrophy by 10 so that can see the other arrows well on the same graph... be sure to note this in the legend
arrows_df<-read.csv("arrows_df_faprotax_with_cyanos_revised.csv",check.names = F)
rownames(arrows_df)<-arrows_df[,1]
arrows_df[,1]<-NULL

########get values of variance explained by each of the axis
sum(as.vector((trait_pcoa_arrows$values$Eigenvalues)/sum(trait_pcoa_arrows$values$Eigenvalues))[1])
sum(as.vector((trait_pcoa_arrows$values$Eigenvalues)/sum(trait_pcoa_arrows$values$Eigenvalues))[2])



faprotax_plot<-ggplot(NULL, aes(Axis.1, Axis.2)) +
  geom_point(data = vectors,size=2,aes(color=type)) +
  stat_ellipse(data = vectors, geom = "polygon", aes(color=type), alpha = 0, show.legend = FALSE, level = 0.95, size=1.1) +
  #geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=Axis.1,yend=Axis.2),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.470732244,yend=-0.360190741),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.009095556,yend=-0.331843426),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=-0.009544398,yend=-0.344202783),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.282381173,yend=-0.281107576),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.424292337,yend=0.328038585),arrow=arrow(length=unit(3,"mm")))+
  geom_segment(data = arrows_df, x=0,y=0,alpha=0.7,size=1,color="red",mapping=aes(xend=0.241334528,yend=-0.314326565),arrow=arrow(length=unit(3,"mm")))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 25),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=17),legend.text=element_text(size=17), plot.title=element_text(size=30, face="bold"))+
  scale_shape_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c(15,16,17,18))+
  scale_color_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("#E69F00","chocolate4","#56B4E9", "#009E73"))+
  xlab("\nAxis 1 (46.3% of variance)")+
  ylab("Axis 2 (20.5% of variance)\n")+
  annotate("text", x=-1.0, y=1.25,label=paste("list(Adonis:~F['3,989']==80.1", ", p==0.001", ", R^{2}==0.196)"),parse=TRUE,hjust=0.0,size=9.5)+
  annotate("text",label="cyano/chloro",x=-0.6  ,y=-0.4, size=9)+
  annotate("text",label="PT",x=0  ,y=-0.4, size=9)+
  annotate("text",label="PHt",x=0  ,y=-0.5, size=9)+
  annotate("text",label="Aer chemoHt",x=0.35 ,y=-0.35, size=9)+
  annotate("text",label="Ferment",x=0.4,y=0.4, size=9)+
  #annotate("text",label="Dark Sulfide Oxidation",x=0.202351949 ,y=-0.217865809, size=3)+
  #annotate("text",label="Dark Oxidation of Sulfur Compounds",x=0.205984682 ,y=-0.237815982, size=3)+
  #annotate("text",label="Dark Sulfur Oxidation",x=0.202351949 ,y=-0.217865809, size=3)+
  annotate("text",label="HC degrad",x=0.37 ,y=-0.45, size=9)
#annotate("text",label="Dark Thiosulfate Oxidation",x=0.087238963 ,y=-0.128352774, size=3)+
#annotate("text",label="Dark Hydrogen Oxidation",x=0.067511943  ,y=0.066513302, size=3)+
#annotate("text",label="Nitrate Reduction",x=0.165053410 ,y=-0.007019182, size=3)+
#annotate("text",label="Methanol Oxidation",x=0.050733196  ,y=0.094051299, size=3)+
#annotate("text",label="Nitrate Respiration",x=0.060395055  ,y=0.039101595, size=3)+
#annotate("text",label="Nitrogen Respiration",x=0.060398099  ,y=0.039104329, size=3)+
#annotate("text",label="Methylotrophy",x=0.049861743  ,y=0.101550511, size=3)

faprotax_plot +labs(title = "Bacterial PCoA (FAPROTAX)")


pdf("Faprotax_PcoA_Mar23.pdf", width = 12, height = 8, pointsize=12)
faprotax_plot +labs(title = "Bacterial PCoA (FAPROTAX)")
dev.off()


#######boxplots 

#######Figure 5a

data2<-read.csv("Madin_Pipeline_final_results_simple.csv",header=T,row.names=1)
data2<-t(data2)

data2<-as.data.frame(data2)
data2$Sample<-rownames(data2)
rownames(data2)<-NULL


metadata<-as(sample_data(rare_bacterial_phyloseq), "matrix")
metadata<-as.data.frame(metadata)
metadata$Sample<-rownames(metadata)
rownames(metadata)<-NULL

data<-dplyr::full_join(data2,metadata,by="Sample")
#data<-data[!is.na(data$Genome_size),]

# specify the factor levels in the order you want
data$type <- factor(data$type, levels = c("environmental_soil", "environmental_water", "environmental_aerosol", "experimental"))
#####ORDER CHANGED  SOIL -> Water -> Aerosol -> Tank


#######Cell_length

result<-aov(data$Cell_length_upper~data$type,data=data)
summary(result)
TukeyHSD(result)


Cell_length_plot<-ggplot(data, aes(x=type,y = Cell_length_upper, fill=type)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  #theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5,size=14), axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),legend.text.align = 1,legend.key=element_rect(fill=NA),axis.line = element_line(colour = "black"),axis.text=element_text(size = 18),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=18),legend.text=element_text(size=18))
  ylab("Cell Length")+
  xlab("")+
  annotate("text",x = 0.65, y = 15000, label="A",size =6)+
  annotate("text",x = 1.65,y = 11000, label="C",size =6 )+
  annotate("text",x = 2.65, y = 13000, label="B",size =6)+
  annotate("text",x = 3.65, y = 13000, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))

Cell_length_plot


pdf("Cell_length_plot_AXIS.pdf", width = 8, height = 4, pointsize=12) #######change name of graph, for example barplot_Cilia or barplot_Size_Md
Cell_length_plot
dev.off()

#####Sporulation

result<-aov(data$Sporulation~data$type,data=data)
summary(result)
TukeyHSD(result)



sporulation_plot<-ggplot(data, aes(x=data$type,y = Sporulation, fill=data$type)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  ylab("Sporulation")+
  xlab("")+
  annotate("text",x = 0.65, y = 1000, label="A",size =6)+
  annotate("text",x = 1.65, y = 250, label="C",size =6)+
  annotate("text",x = 2.65, y = 1000, label="A",size =6)+
  annotate("text",x = 3.65, y = 350, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))
sporulation_plot

pdf("sporulation_plot_AXIS.pdf", width = 8, height = 4, pointsize=12) 
sporulation_plot
dev.off()


####Cell diameter

result<-aov(data$Cell_diameter_upper~data$type,data=data)
summary(result)
TukeyHSD(result)

Cell_diameter_plot<-ggplot(data, aes(x=type,y = Cell_diameter_upper, fill=type)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  ylab("Cell Diameter")+
  xlab("")+
  annotate("text",x = 0.65, y = 4500, label="B",size =6)+
  annotate("text",x = 1.65, y = 1900, label="D",size =6)+
  annotate("text",x = 2.65, y = 4900, label="A",size =6)+
  annotate("text",x = 3.65, y = 2200, label="C",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))
Cell_diameter_plot


pdf("Cell_diameter_plot_AXIS.pdf", width = 8, height = 4, pointsize=12)
Cell_diameter_plot
dev.off()

####Doubling times
result<-aov(data$Doubling_Time_hours~data$type,data=data)
summary(result)
TukeyHSD(result)

Cell_doubling_time_plot<-ggplot(data, aes(x=type,y = Doubling_Time_hours, fill=type)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.2), colour="grey27")+
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red") +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 15),axis.title=element_text(size=25,face="plain"))+
  ylab("Doubling Time")+
  xlab("")+
  annotate("text",x = 0.65, y = 6000, label="A",size =6)+
  annotate("text",x = 1.65, y = 6000, label="A",size =6)+
  annotate("text",x = 2.65, y = 4000, label="B",size =6)+
  annotate("text",x = 3.65, y = 4000, label="B",size =6)+
  scale_x_discrete(name=NULL,labels=c("environmental_soil"="Terrestrial \n source","environmental_aerosol"="Aerial \n dispersers","environmental_water"="Aquatic \n source","experimental"="Tank \n colonizers"))+
  scale_fill_manual(name="Environment",labels=c("environmental_soil"="Terrestrial source","environmental_aerosol"="Aerial dispersers","environmental_water"="Aquatic source","experimental"="Tank colonizers"),values=c("chocolate4","#56B4E9","#E69F00", "#009E73"))
Cell_doubling_time_plot


pdf("Cell_doubling_time_plot_AXIS.pdf", width = 8, height = 4, pointsize=12) 
Cell_doubling_time_plot
dev.off()
