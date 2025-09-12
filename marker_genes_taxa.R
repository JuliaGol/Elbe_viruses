#setting
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("~/IGB_phd/BICEST/virus/viral_abundances")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/workshop_gene_catalog"
##Loading files
#save in the environment
load("~/IGB_phd/BICEST/virus/viral_abundances/marker_genes.RData")
#gene_cluster_table
clstr <- read.csv(paste0(path,"/GROS22.clstr"), sep="\t", skip=0, header = F)
colnames(clstr) <- c("gene_cluster", "mem", "gene", "length", "pos")
#viral genes annotated
# vibrant_annot <- read.csv("BICEST.phages.annot",  sep = "\t", header = T)
# colnames(vibrant_annot) <- c("genome", "range", "KO", "DESCRIPTION")
# vibrant_annot[,"genome"] = sub(">", "", vibrant_annot[,"genome"])
# #BICEST annotation of gene clusters
# annot <- read.csv(paste0(path,"/GROS22.kegg-annotations.tsv"), sep="\t", skip=0, header = TRUE, quote="")

#taxa assignment for MAGs
# taxa <- read.csv(paste0(path,"/GROS22-2.gtdb"), sep="\t", skip=0, header = TRUE)
# taxa <- taxa %>% 
#   separate(GTDBTK_TAXONOMY, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";")
#abundance
# geneabund <- read.csv(paste0(path,"/BICEST_cellabund.tsv"), sep="\t", skip=1, header = TRUE, row.names = 1)

path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) 
#remove rows with NA in Acession number - no molecular data
metadata <- metadata[which(!is.na(metadata$AccessionNumber_TBDSven)),]
#drep clusters for finding gene markers 
drep_clusters <- read.csv("Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])
#read combine quality file 
#check how many of each quality 
#
quality <- read.csv("Combined_quality_summary.tsv", header = TRUE, sep = "\t")
drep_clusters <-merge(drep_clusters, quality[, c("contig_id", "checkv_quality")], by.x = "genome", by.y = "contig_id")
#merge the vibrant annotation with gene clusters
vir_genecluster_annot <- vibrant_annot %>% 
  left_join(clstr,  join_by(genome == gene)) #st wrong here - many NAs in final table
#modify the genome column for merging with drep clusters 
vir_genecluster_annot[,"genome"] = sub("_[1-9]+$", "", vir_genecluster_annot[,"genome"])
vibrant_genecluster_drep = merge(drep_clusters, vir_genecluster_annot, by="genome")  %>%
  select(genome, gene_cluster, KO, primary_cluster, secondary_cluster, DESCRIPTION)
#add abundance data
geneabund <- geneabund %>% 
  rownames_to_column("gene_cluster")
#but first filter out MCPs to reduce calculations
# MCP_list <- c("VOG00035", "VOG00461", "VOG01150", "VOG01000", "VOG01164", "VOG01920", "VOG02155", "VOG02220", "VOG02554", "VOG03097", "VOG03553", "VOG03780", "VOG03813", "VOG33139")
# vibrant_genecluster_drep <- vibrant_genecluster_drep %>% filter(KO %in% MCP_list)
vibrant_genecluster_drep_abund <- left_join(vibrant_genecluster_drep, geneabund, by = "gene_cluster") 
#check how many drep_clusters per one MCP gene cluster
vibrant_genecluster_drep_abund_primary_counts <- vibrant_genecluster_drep_abund %>% group_by (primary_cluster,gene_cluster) %>% summarise(n = n())
#can be more then one 
#how many more then one?
vibrant_genecluster_drep_abund_primary_counts[which(vibrant_genecluster_drep_abund_primary_counts$n >1), ] %>% nrow()
#mean(vibrant_genecluster_drep_abund_primary_counts$n) 3.587838 <0 averege how many clusters
#median(vibrant_genecluster_drep_abund_primary_counts$n) 2
#how many exacly one
vibrant_genecluster_drep_abund_primary_counts[which(vibrant_genecluster_drep_abund_primary_counts$n == 1),] %>% nrow()
#the same for secondary clusters 
vibrant_genecluster_drep_abund_secondary_counts <- vibrant_genecluster_drep_abund %>% group_by(secondary_cluster,gene_cluster) %>% summarise(n = n())
vibrant_genecluster_drep_abund_secondary_counts[which(vibrant_genecluster_drep_abund_secondary_counts$n >1), ] %>% nrow()
mean(vibrant_genecluster_drep_abund_secondary_counts$n) #3.587838 <0 averege how many clusters
median(vibrant_genecluster_drep_abund_secondary_counts$n) #2
#how many exacly one
vibrant_genecluster_drep_abund_secondary_counts[which(vibrant_genecluster_drep_abund_secondary_counts$n <= 1),] %>% nrow()
#how many gene clusters
length(unique(vibrant_genecluster_drep_abund_secondary_counts$gene_cluster)) #2894
#the same resutls for primary and secondary - check it out if is one secondary per primary  
vibrant_genecluster_drep_abund_clusters <- vibrant_genecluster_drep_abund %>% group_by(gene_cluster,secondary_cluster, primary_cluster)
#looks like 
#ok - so MCP not good enougth to have a good resolution
#do we have one per genome?
length(unique(vir_genecluster_annot_abund$genome)) == length(vir_genecluster_annot_abund$genome)
#12737 - yes!
#try heatmap to check if there are expetionally high signals 
#gene clust as rownames to indentify exceptionary high counts
vir_genecluster_annot_abund_heatmap <- vir_genecluster_annot_abund %>% select(matches("GROS|cluster")) %>% distinct()
rownames(vir_genecluster_annot_abund_heatmap) <- vir_genecluster_annot_abund_heatmap$gene_cluster
heatmap(as.matrix(vir_genecluster_annot_abund_heatmap[,2:length(vir_genecluster_annot_abund_heatmap[1,])]), scale = "none")
#check
#save original with all KO
vibrant_genecluster_drep = merge(drep_clusters, vir_genecluster_annot, by="genome")  %>%
  select(genome, gene_cluster, KO, primary_cluster, secondary_cluster, DESCRIPTION)
saveRDS(vibrant_genecluster_drep, "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/vibrant_genecluster_drep.RDS")

#find primary clusters using all KO
vibrant_annot_drep_marker <- vibrant_genecluster_drep %>%
  filter(KO != "None") %>% 
  group_by(gene_cluster) %>% 
  mutate(count_primary_cluster = n_distinct(primary_cluster)) %>%
  filter(count_primary_cluster == 1)  

vibrant_annot_drep_marker_gene_cluster <- vibrant_annot_drep_marker[,"gene_cluster"]  
#how many primary clusters has marker genes 
n_distinct(vibrant_annot_drep_marker[, "primary_cluster"])/max(vibrant_annot_drep_marker[, "primary_cluster"]) #max(vibrant_genecluster_drep[, "primary_cluster"]) ##99,54% of clusters covered!!! (sth wrong 2372 clusters)
n_distinct(vibrant_annot_drep_marker[,"genome"])/length(unique(drep_clusters$genome)) #not correct because it can be that 
#gene wasnt detected at some genomes but were asisgned to cluster
#[1] 10861 0 also strange number ? #n_distinct(vibrant_genecluster_drep[,"genome"]) #99,58% genomes!!! 
#not all genomes in vibrant_annot_drep_marker - keep only valid columns for clarity 
#save only first row from each cluster 
vibrant_annot_drep_marker <- vibrant_annot_drep_marker %>%
  group_by(primary_cluster) %>%
  mutate(count_genome = n_distinct(genome)) %>% #count genomes per cluster (not per gene) for statistics
  arrange(desc(KO)) %>% #arrange KO in descending order so we will favor VOG 
  filter(row_number()==1) %>%
  select(!c("genome")) #remove genome column as this is table sumarising glusters and their marker genes 
#statistics 
n_distinct(vibrant_annot_drep_marker[, "primary_cluster"])/max(vibrant_annot_drep_marker[, "primary_cluster"]) #max(vibrant_genecluster_drep[, "primary_cluster"]) ##99,54% of clusters covered!!! 
sum(vibrant_annot_drep_marker[,"count_genome"])/length(unique(drep_clusters$genome)) #not correct because it can be that some genomes were filter out together 

1-(length(which(!drep_clusters$primary_cluster %in% vibrant_annot_drep_marker$primary_cluster) )/length(unique(drep_clusters$genome)))
#relative frequency of genomes ith assigned marker 0.997514
#drep_clusters$primary_cluster[which(!(drep_clusters$primary_cluster %in% vibrant_annot_drep_marker$primary_cluster))]
saveRDS(vibrant_annot_drep_marker, file = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/vibrant_annot_drep_marker.RDS")
#let make a list of gene_clusters which marker genes belong to 
gene_cluster_marker <- vibrant_annot_drep_marker$gene_cluster
#TODO blast gene_cluster marker sequences 
#prepare geneabund for merging  
#summarise new gene counts list
vibrant_annot_drep_marker %>% group_by(DESCRIPTION) %>% summarise(n = n())  %>% arrange(desc(n))

#redo host analysis similarly to the one based on MCP
#merging the current vibrant_annot_drep_marker dataframe with data about gene clusters abundance 
geneabund_marker_gene_cluster <- geneabund[which(geneabund$gene_cluster %in% gene_cluster_marker),]
vibrant_annot_drep_marker_abund <- vibrant_annot_drep_marker %>% left_join(geneabund, by = "gene_cluster")
saveRDS(vibrant_annot_drep_marker_abund, file ="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/vibrant_annot_drep_marker_abund.RDS")


###interfere with Iphop data
host_prediction_genome = read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/Host_prediction_to_genome_m90.csv")
host_prediction_genus = read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/Host_prediction_to_genus_m90.csv")
#host_prediction_genus 5 times less observations
length(host_prediction_genome$Virus)
# 5424
length(unique(host_prediction_genome$Virus)) 
# 2043 #more different hosts - more specific taxonomy
length(host_prediction_genus$Virus)
# 1119
length(unique(host_prediction_genus$Virus)) 
# 1036
#check - there is a difference in viruses lists - longer in host_prediction_genome
#is host_prediction_genus$Virus  a subset of host_prediction_genome$Virus
host_prediction_genus$Virus %in% host_prediction_genome$Virus #no 
#add cluster number to each virus
host_prediction_genus <- merge(host_prediction_genus, drep_clusters,  by.x = "Virus", by.y = "genome")
#the same with genome file - hosts obtained differently 
host_prediction_genome <- merge(host_prediction_genome, drep_clusters,  by.x = "Virus", by.y = "genome")
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/marker_genes_taxa.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")
#prepare lineage to genus only host_prediction_genome and division to different taxa in the column
host_prediction_genome <- host_prediction_genome %>% 
  separate(Host.taxonomy, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";")
host_prediction_genome$lineage_genus <- paste(host_prediction_genome$domain, host_prediction_genome$phylum, host_prediction_genome$class, host_prediction_genome$order, host_prediction_genome$family,  host_prediction_genome$genus, host_prediction_genome$species, sep=";")

#how many host carry each cluster

host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(lineage_genus)) %>%
  filter(count_hosts == 1)  %>% nrow
#427 
host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(species)) %>%
  filter(count_hosts == 1)  %>% nrow
#490 
host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(lineage_genus)) %>%
  filter(count_hosts > 1)  %>% nrow
#257


host_prediction_genus %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts == 1)  %>% nrow
#373
host_prediction_genus %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts > 1)  %>% nrow
#52  
host_prediction_genus %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts == 0)  %>% nrow
#0
host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(lineage_genus)) %>%
  filter(count_hosts == 0)  %>% nrow
#0
#in host_prediction_genome more observation, thats why we have higher number 

#most of the clusters has one host genus!!!!!
#use only clusters with marker genes and reevaluate 
#merge dataframes to get the virus genome name (Vibrant) and counts of marker gene in the same table 
#vibrant_annot_drep_marker_counts <- left_join(vibrant_annot_drep_marker, geneabund_marker_gene_cluster, by = "gene_cluster")
#using the viral genomes names we can interefere the results with iphop 
host_prediction_genus_counts <- left_join(host_prediction_genus, vibrant_annot_drep_marker_abund[,!colnames(vibrant_annot_drep_marker_abund) %in% c("secondary_cluster", "genome", "count_primary_cluster")], by="primary_cluster")
#the same for genome file 
host_prediction_genome_counts <- left_join(host_prediction_genome, vibrant_annot_drep_marker_abund[,!colnames(vibrant_annot_drep_marker_abund) %in% c("secondary_cluster", "genome", "count_primary_cluster")], by="primary_cluster")

#check how many host in primary cluster now (some are filtered out as they did not have marker gene)
host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>% 
  filter(count_hosts > 1) %>%
  nrow
#52

host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts == 1) %>%
  nrow
#373
host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts < 1) %>%
  nrow
n_distinct(host_prediction_genus_counts[, "primary_cluster"]) 
#425 the same as above


# how to deal in the plot with different hosts coming from the same virus
# focus on virus, so use counts even they are redundant, as we want to know how many viruses have particular host, not how many viruses are in total we want to see dynamic
host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genus_counts)[1:16], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 
#consider phyla not genera
host_prediction_genus_counts_metadata <- host_prediction_genus_counts_metadata %>%
  separate(Host.genus, c("domain", "phylum", "class", "order", "family",  "genus"), ";")

# the same for genome counts 
host_prediction_genome_counts_metadata <- host_prediction_genome_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genome_counts)[1:24], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 
host_prediction_genus_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=Host.genus)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts")  +
  theme(legend.position = "none")
host_prediction_genome_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=lineage_genus)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts - genome output iphop")  +
  theme(legend.position = "none")
#the same plot as above but for phyla
host_prediction_genus_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=phylum)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts")  +
  theme(legend.position = "none")

host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genus_counts)[1:16], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 
#the same for host_predction-genome ...
host_prediction_genus_genome_metadata <- host_prediction_genome_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genome_counts)[1:23], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 

host_prediction_genus_counts_metadata$Station <- factor(host_prediction_genus_counts_metadata$Station, levels = c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund"))
host_prediction_genome_counts_metadata$Station <- factor(host_prediction_genome_counts_metadata$Station, levels = c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund"))

host_prediction_genus_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=Host.genus)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts")

host_prediction_genus_counts_metadata$Sample_date <- factor(host_prediction_genus_counts_metadata$Sample_date, levels= c("May-21", "Jul-21", "Nov-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
host_prediction_genome_counts_metadata$Sample_date <- factor(host_prediction_genome_counts_metadata$Sample_date, levels= c("May-21", "Jul-21", "Nov-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
host_prediction_genus_counts_metadata$Stromkilometer <- factor(round(host_prediction_genus_counts_metadata$Stromkilometer,0))
#check how many clusters have different phyla as hosts (multiple hosts phyla)
host_prediction_genus_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 1) %>%
  nrow                   
#399
host_prediction_genus_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum > 1) %>%
  nrow  
#26
host_prediction_genus_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 0) %>%
  nrow
#keep only distinct phyla per primary cluster 

host_prediction_genus_counts_metadata <- host_prediction_genus_counts_metadata  %>% 
  select(c("phylum", "primary_cluster", "counts", "Station", "Sample_date", "data_type")) %>% 
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov-21") %>%
  #filter(data_type != "METAT") %>%
  mutate(counts=counts/1000) %>%
  mutate(phylum=str_remove(phylum, "p__")) %>% 
  distinct() 
host_prediction_genus_counts_metadata$Sample_date <- factor(host_prediction_genus_counts_metadata$Sample_date, levels=c("May-21", "Jul-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
tiff("geneabund_drep_marker_iphop_genus_plots.tiff", unit="px", width = 700, height = 1000) 
host_prediction_genus_counts_metadata  %>%
  ggplot() + 
  geom_tile(aes(x=Station, y=phylum, fill=log(counts)), na.rm = T) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Log Copies\n/Transcripts per genome", y="phylum") +
  ggtitle("Hosts of viruses based on genus iPHoP prediction") + 
  scale_x_discrete(labels = c("633", "651", "665", "692", "712"), name="Elbe kilometer")
dev.off()

host_prediction_genome_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 1) %>%
  nrow                   
#558
host_prediction_genome_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum > 1) %>%
  nrow  
#126
host_prediction_genome_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 0) %>%
  nrow
#0
host_prediction_genome_counts_metadata <- host_prediction_genome_counts_metadata  %>% 
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov-21") %>%
  # filter(data_type != "METAT") %>%
  mutate(counts=counts/1000) %>%
  mutate(phylum=str_remove(phylum, "p__")) %>% distinct()
  
host_prediction_genome_counts_metadata$Sample_date <- factor(host_prediction_genome_counts_metadata$Sample_date, levels=c("May-21", "Jul-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))

tiff("geneabund_drep_marker_iphop_genome_plots.tiff", unit="px", width = 700, height = 1300) 
host_prediction_genome_counts_metadata  %>% 
  ggplot() + 
  geom_tile(aes(x=Station,y=phylum, fill=log(counts))) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_text("Elbe kilometer"), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Log Copies\n/Transcripts per genome", y="phylum") +
  ggtitle("Hosts of viruses based on genome iPHoP prediction")  + 
  scale_x_discrete(labels = c("633", "651", "665", "692", "712"), name="Elbe kilometer")
dev.off()

# another plot to see how many generalist vs specialist
  host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "Host.genus")) %>% 
  summarise(count_hosts = n_distinct(Host.genus)) %>% arrange(count_hosts) %>% tail()
    
    
