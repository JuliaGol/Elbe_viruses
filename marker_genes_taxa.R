#setting
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("~/IGB_phd/virus/viral_abundances")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog"
##Loading files
#gene_cluster_table
clstr <- read.csv(paste0(path,"/GROS22.clstr"), sep="\t", skip=0, header = F)
colnames(clstr) <- c("gene_cluster", "mem", "gene", "length", "pos")
#viral genes annotated
vibrant_annot <- read.csv("BICEST.phages.annot",  sep = "\t", header = T)
colnames(vibrant_annot) <- c("genome", "range", "KO", "DESCRIPTION")
vibrant_annot[,"genome"] = sub(">", "", vibrant_annot[,"genome"])
#BICEST annotation of gene clusters
annot <- read.csv(paste0(path,"/GROS22.kegg-annotations.tsv"), sep="\t", skip=0, header = TRUE, quote="")

#taxa assignment for MAGs
taxa <- read.csv(paste0(path,"/GROS22-2.gtdb"), sep="\t", skip=0, header = TRUE)
taxa <- taxa %>% 
  separate(GTDBTK_TAXONOMY, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";")
#abundance
geneabund <- read.csv(paste0(path,"/BICEST_cellabund.tsv"), sep="\t", skip=1, header = TRUE, row.names = 1)
##Load Metadata for later

metadata <- read.csv("SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAT.genecount.profile")) %>%
  mutate(data_type="METAT")


metadata2 <- read.csv("SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAG.genecount.profile"))%>%
  mutate(data_type="METAG")

metadata <- rbind(metadata, metadata2)
metadata <-metadata %>% 
  mutate(sample=str_remove(sampleid, ".genecount.profile")) %>%
  mutate(sample=str_replace(sample, "[.]", "-"))

#drep clusters for finding gene markers 
drep_clusters <- read.csv("Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])

#


#merge the vibrant annotation with gene clusters
vir_genecluster_annot <- vibrant_annot %>% 
  left_join(clstr,  join_by(genome == gene))
#modify the genome column for merging with drep clusters 
vir_genecluster_annot[,"genome"] = sub("_[1-9]+$", "", vir_genecluster_annot[,"genome"])
vibrant_genecluster_drep = merge(drep_clusters, vir_genecluster_annot, by="genome")  %>%
  select(genome, gene_cluster, KO, primary_cluster)

#find primary clusters using all KO
vibrant_annot_drep_marker <- vibrant_genecluster_drep %>%
  filter(KO != "None") %>% 
  group_by(gene_cluster) %>% 
  mutate(count_primary_cluster = n_distinct(primary_cluster)) %>%
  mutate(count_genome = n_distinct(genome)) %>%
  filter(count_primary_cluster == 1)  

vibrant_annot_drep_marker_gene_cluster <- vibrant_annot_drep_marker[,"gene_cluster"]  
#how many primary clusters has marker genes 
n_distinct(vibrant_annot_drep_marker[, "primary_cluster"])/max(vibrant_genecluster_drep[, "primary_cluster"]) ##99,54% of clusters covered!!!
n_distinct(vibrant_annot_drep_marker[,"genome"])/n_distinct(vibrant_genecluster_drep[,"genome"]) #99,58% genomes!!! 

#let make a list of gene_clusters which marker genes belong to 
gene_cluster_marker <- vibrant_annot_drep_marker$gene_cluster

#prepare geneabund for merging  
geneabund <- geneabund %>% 
  rownames_to_column("gene_cluster")

#redo host analysis similarly to the one based on MCP
#merging the current vibrant_annot_drep_marker dataframe with data about gene clusters abundance 
geneabund_marker_gene_cluster <- geneabund[which(geneabund$gene_cluster %in% gene_cluster_marker),]
geneabund_marker_gene_cluster_MAGs <- left_join(geneabund_marker_gene_cluster, clstr) %>%  
  filter(grepl("*_MAG_*", gene)) %>%
  mutate(genome=str_remove(gene, "-scaff.*")) 

geneabund_drep_marker_taxa <- left_join(geneabund_marker_gene_cluster_MAGs, taxa, by=c("genome"="GENOME"), relationship = "many-to-many") %>%
  select(-c("length", "mem", "pos","GTDBTK_VERSION", "GTDBTK_DATABASE"))

geneabund_drep_marker_taxa_plots <- geneabund_drep_marker_taxa %>% 
  select(-c("class", "order", "family", "genus", "species", "gene", "genome"))  %>%
  pivot_longer(!c("domain", "phylum",  "gene_cluster"),  names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata[, c("Station", "Sample_date", "sampleid", "data_type")])
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")

##abundance of the viral hosts, based on the veeery general assumption that possessing cluster marker means being host for viruses bossessing the same gene cluster    
geneabund_drep_marker_taxa_plots %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=factor(Station, levels= c("Muehlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuettel","Meedem Grund")), 
                y=phylum, fill=log(counts))) +
  facet_wrap(~ data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  labs(fill = "Copies/Transcripts per genome", y="Phylum")
