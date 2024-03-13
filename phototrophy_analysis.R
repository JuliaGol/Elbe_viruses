#script for phototrophy analysis in the Elbe
##setting
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("~/IGB_phd/phototrophy")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog"
##Loading files
#gene_cluster_table
clstr <- read.csv(paste0(path,"/GROS22.clstr"), sep="\t", skip=0, header = F)
colnames(clstr) <- c("gene_cluster", "mem", "gene", "length", "pos")
#annotation
annot <- read.csv(paste0(path,"/GROS22.kegg-annotations.tsv"), sep="\t", skip=0, header = TRUE, quote="")
#abundance
geneabund <- read.csv(paste0(path,"/BICEST_cellabund.tsv"), sep="\t", skip=1, header = TRUE, row.names = 1)
#prepare geneabund for merging  
geneabund <- geneabund %>% 
  rownames_to_column("gene_cluster")
#prepare metadata fro merging
metadata <- read.csv(paste0(path,"/SAMEAID_SampleID_simplified.csv"), header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAT.genecount.profile")) %>%
  mutate(data_type="METAT")

metadata2 <- read.csv(paste0(path, "/SAMEAID_SampleID_simplified.csv"), header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAG.genecount.profile")) %>%
  mutate(data_type="METAG")

metadata <- rbind(metadata, metadata2)


#let's find rhodopsins using KO
#here the ids from KEGG database
#K00909 GRK1_7; rhodopsin kinase [EC:2.7.11.14]
#K04641 bop; bacteriorhodopsin
#K04642 hop; halorhodopsin
#K04643 sop; sensory rhodopsin
rhd_clust <- annot %>%  
  filter(KO %in% c("K00909","K04641","K04642","K04643")) %>% 
  left_join(geneabund,  join_by(QUERY == gene_cluster))
#not enough to do meaningful analysis 
#interestigly all rhodopsins were found in viral scaffolds ?

#prepare the table with matadata for plot
meta_rhd_clust <- rhd_clust %>% 
  pivot_longer(!c("QUERY", "KO", "BRITE", "CAZY", "COG", "DISEASE", "DRUG", "ENZYME", "GO", "MODULE", "NETWORK", "PATHWAY", "PUBMED", "RCLASS", "REACTION", "TC", "VP", "DESCRIPTION", "length"), names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata)
meta_rhd_clust$Stromkilometer <- sub(",", ".", meta_rhd_clust$Stromkilometer)
meta_rhd_clust$Stromkilometer <- as.double(meta_rhd_clust$Stromkilometer)

#plot it!
meta_rhd_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot(aes(x=Stromkilometer, y=counts)) + geom_point() +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Rhd counts across different seasons and river sites")
#not enough to do meaningful analysis 
#interestigly all rhodopsins were found in viral scaffolds ?
#usually the hiest variance in the last station

#check different fractions
#plot it!
meta_rhd_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot(aes(x=Stromkilometer, y=counts)) + geom_point() +
  geom_tile()  + facet_wrap(~ Sample_type) 
#nothing interesting

#explore correlations for METAT counts
meta_rhd_clust_matrix <-  meta_rhd_clust  %>% filter(data_type == "METAT") %>%
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>% 
  select( c("counts","Stromkilometer", "SPM_mgperL", "DOC_mg.L", "TN_mg.L", "DIC_mg.L", "PTH_mgperL", "POC_mgperL", "PTC_mgperL", "PTN_mgperL")) 
 
meta_rhd_clust_matrix[ meta_rhd_clust_matrix == ""] <- NA 

for (i in c(3:length(meta_rhd_clust_matrix[1,]))){
  meta_rhd_clust_matrix[,i] <-  sub(",", ".", meta_rhd_clust_matrix[,i]) #replace , to . to get a string which can be converted to double 
}

head(meta_rhd_clust_matrix)
#to numeric 

meta_rhd_clust_matrix <- matrix(apply(meta_rhd_clust_matrix,2, as.numeric), ncol=10)
meta_rhd_clust_matrix[ is.na(meta_rhd_clust_matrix)] <- 0 
res <- cor(meta_rhd_clust_matrix)
round(res, 2)
#veeery strong correlation between "DOC_mg.L", "TN_mg.L", 
#quite strong between  "TN_mg.L", "DIC_mg.L"
var_meta_rhd_clust_station <- meta_rhd_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>% 
  group_by(Station) %>% 
  summarise(var = var(counts))

#variance across stations (similarly the variance is lower n the brakish water)
var_meta_rhd_clust_station  %>% arrange(factor(Station, levels = c("Muehlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuettel","Meedem Grund")))

#repeat similar analysis for phototrophy psbA - K02703

psbA_clust <- annot %>%  
  filter(KO %in% c("K02703")) %>% 
  left_join(geneabund,  join_by(QUERY == gene_cluster))
#not enough to do meaningful analysis 
#interestigly all rhodopsins were found in viral scaffolds ?

#prepare the table with matadata for plot
meta_psbA_clust <- psbA_clust %>% 
  pivot_longer(!c("QUERY", "KO", "BRITE", "CAZY", "COG", "DISEASE", "DRUG", "ENZYME", "GO", "MODULE", "NETWORK", "PATHWAY", "PUBMED", "RCLASS", "REACTION", "TC", "VP", "DESCRIPTION", "length"), names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata)
meta_psbA_clust$Stromkilometer <- sub(",", ".", meta_psbA_clust$Stromkilometer)
meta_psbA_clust$Stromkilometer <- as.double(meta_psbA_clust$Stromkilometer)

#plot it!
meta_psbA_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot(aes(x=Stromkilometer, y=counts)) + ggtitle("psbA counts across different seasons and river sites") + geom_point() + 
  geom_tile() + facet_wrap(~data_type *Sample_date) + geom_smooth() 
#usually the highest variance in the beginning and at the end of the estuary 
#but very similar across  the samples

#check different fractions
#plot it!
meta_psbA_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot(aes(x=Stromkilometer, y=counts)) + geom_point() +
  geom_tile() + facet_wrap(~ Sample_type) + geom_smooth()

#more in particle associated 

#explore correlations
meta_psbA_clust_matrix <-  meta_psbA_clust  %>%  
  filter(data_type == "METAT") %>%
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>% 
  select( c("counts","Stromkilometer", "SPM_mgperL", "DOC_mg.L", "TN_mg.L", "DIC_mg.L", "PTH_mgperL", "POC_mgperL", "PTC_mgperL", "PTN_mgperL")) 

meta_psbA_clust_matrix[ meta_psbA_clust_matrix == ""] <- NA 

for (i in c(3:length(meta_psbA_clust_matrix[1,]))){
  meta_psbA_clust_matrix[,i] <-  sub(",", ".", meta_psbA_clust_matrix[,i]) #replace , to . to get a string which can be converted to double 
}

head(meta_psbA_clust_matrix)
#to numeric 

meta_psbA_clust_matrix <- matrix(apply(meta_psbA_clust_matrix,2, as.numeric), ncol=10)
meta_psbA_clust_matrix[ is.na(meta_psbA_clust_matrix)] <- 0 
res <- cor(meta_psbA_clust_matrix)
round(res, 2)

#variance analysis 
var_meta_psbA_clust <- meta_psbA_clust  %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>% 
  group_by(Station, Sample_date) %>% 
  summarise(var = var(counts))

var_meta_psbA_clust_station <- meta_psbA_clust %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>% 
  group_by(Station) %>% 
  summarise(var = var(counts))

#variance across stations (similarly the variance is lower n the brakish water)
var_meta_psbA_clust_station  %>% arrange(factor(Station, levels = c("Muehlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuettel","Meedem Grund")))

