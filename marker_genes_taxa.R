#setting
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("~/IGB_phd/BICEST/virus/viral_abundances")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/workshop_gene_catalog"
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

# metadata <- read.csv("SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
#   mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAT.genecount.profile")) %>%
#   mutate(data_type="METAT")
# 
# 
# metadata2 <- read.csv("SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
#   mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAG.genecount.profile"))%>%
#   mutate(data_type="METAG")
# 
# metadata <- rbind(metadata, metadata2)
# metadata <-metadata %>% 
#   mutate(sample=str_remove(sampleid, ".genecount.profile")) %>%
#   mutate(sample=str_replace(sample, "[.]", "-"))
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) 
#remove rows with NA in Acession number - no molecular data
metadata <- metadata[which(!is.na(metadata$AccessionNumber_TBDSven)),]
#drep clusters for finding gene markers 
drep_clusters <- read.csv("Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])

#


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
MCP_list <- c("VOG00035", "VOG00461", "VOG01150", "VOG01000", "VOG01164", "VOG01920", "VOG02155", "VOG02220", "VOG02554", "VOG03097", "VOG03553", "VOG03780", "VOG03813", "VOG33139")
vibrant_genecluster_drep <- vibrant_genecluster_drep %>% filter(KO %in% MCP_list)
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
saveRDS(vibrant_genecluster_drep, "vibrant_genecluster_drep.RDS")

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
#save only first row from each cluster 
vibrant_annot_drep_marker <- vibrant_annot_drep_marker %>%
  group_by(primary_cluster) %>%
  arrange(desc(KO)) %>% #arrange KO in descending order so we will favour VOG 
  filter(row_number()==1)
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
# geneabund_marker_gene_cluster_MAGs <- left_join(geneabund_marker_gene_cluster, clstr) %>%  
#   filter(grepl("*_MAG_*", gene)) %>%
#   mutate(genome=str_remove(gene, "-scaff.*")) 
geneabund_marker_gene_cluster_MAGs <- left_join(vibrant_annot_drep_marker_abund, clstr["gene", "gene_cluster"], by="gene_cluster") %>%  
  filter(grepl("*_MAG_*", gene)) %>%
  mutate(genome=str_remove(gene, "-scaff.*")) 

geneabund_drep_marker_taxa <- left_join(geneabund_marker_gene_cluster_MAGs, taxa, by=c("genome"="GENOME"), relationship = "many-to-many") %>% select(-c("mem", "pos","GTDBTK_VERSION", "GTDBTK_DATABASE", "length"))
saveRDS(geneabund_drep_marker_taxa, "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/geneabund_drep_marker_taxa.RDS")
#keep only one length column
geneabund_drep_marker_taxa$length.x <- NULL
colnames(geneabund_drep_marker_taxa)[which(colnames(geneabund_drep_marker_taxa) == "length.y")] <- "length"
geneabund_drep_marker_taxa_plots <- geneabund_drep_marker_taxa %>% distinct() %>% 
  select(-c("class", "order", "family", "genus", "species", "gene", "genome"))  %>%
  pivot_longer(col=!c("domain", "phylum",  "gene_cluster", "KO", "primary_cluster", "secondary_cluster", "count_primary_cluster", "count_genome", "length", "DESCRIPTION"),  names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  distinct() %>% merge(metadata[, c("Station", "Sample_date", "sampleid", "data_type", "Sample_type")])
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/marker_genes_taxa.RData")

##abundance of the viral hosts, based on the veeery general assumption that possessing cluster marker means being host for viruses possessing the same gene cluster    
tiff("geneabund_drep_marker_taxa_plots_sample_date_gene_origin_based_taxonomy.tiff", unit="px", width = 700, height = 1000) 
geneabund_drep_marker_taxa_plots %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=factor(Station, levels= c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund")), 
                y=phylum, fill=log(counts))) +
  facet_wrap(~ data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1, size=12), axis.text.y = element_text( size=12)) +
  labs(fill = "Copies/Transcripts per genome", y="Phylum")
dev.off()

tiff("geneabund_drep_marker_taxa_plots_sample_type_gene_origin_based_taxonomy.tiff", unit="px", width = 700, height = 1000) 
geneabund_drep_marker_taxa_plots %>%  
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov 21") %>%
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=factor(Station, levels= c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund")), 
                y=phylum, fill=log(counts))) +
  facet_wrap(~ data_type *Sample_type) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1, size=12), axis.text.y = element_text( size=12)) +
  labs(fill = "Copies/Transcripts per genome", y="Phylum")
dev.off()

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
host_prediction_genus %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts > 1)  %>% nrow
#52  

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
host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(lineage_genus)) %>%
  filter(count_hosts == 1)  %>% nrow
#427 
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
host_prediction_genus_counts <- left_join(host_prediction_genus, vibrant_annot_drep_marker_abund, by=c("Virus" = "genome"))
#the same for genome file 
host_prediction_genome_counts <- left_join(host_prediction_genome, vibrant_annot_drep_marker_abund, by=c("Virus" = "genome"))

#check how many host in primary cluster now (some are filtered out as they did not have marker gene)
host_prediction_genus_counts  %>%
  group_by(primary_cluster.x) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>% 
  filter(count_hosts > 1) %>%
  nrow
#52

host_prediction_genus_counts  %>%
  group_by(primary_cluster.x) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts == 1) %>%
  nrow
#373
host_prediction_genus_counts  %>%
  group_by(primary_cluster.x) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>%
  filter(count_hosts < 1) %>%
  nrow
n_distinct(host_prediction_genus_counts[, "primary_cluster.x"]) 
#425 the same as above

#Still most have exacly one, check how many viruses have more then one host assigned - is it possible in iphop?

host_prediction_genus_counts  %>%
  group_by(Virus) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>% 
  filter(count_hosts > 1) %>%
  nrow
#67 

host_prediction_genus_counts  %>%
  group_by(Virus) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>% 
  filter(count_hosts > 1) %>% arrange(desc(count_hosts))
#the most - 8 hosts assigned 

host_prediction_genus_counts  %>%
  group_by(Virus) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus"))  %>%
  summarise(count_hosts = n_distinct(Host.genus)) %>% 
  filter(count_hosts == 1) %>%
  nrow
# 969
# most specialist?

# how to deal in the plot with different hosts coming from the same virus
# focus on virus, so use counts even they are redundant, as we want to know how many viruses have particular host, not how many viruses are in total we want to see dynamic
host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genus_counts)[1:16], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 
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
#consider phyla not genera
host_prediction_genus_counts_metadata <- host_prediction_genus_counts_metadata %>%
  separate(Host.genus, c("domain", "phylum", "class", "order", "family",  "genus"), ";")
#the same plot as above but for phyla
host_prediction_genus_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=phylum)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts")  +
  theme(legend.position = "none")

host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genus_counts)[1:16], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 
#the same for host_predction-genome ...
host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>% 
  pivot_longer(col = !colnames(host_prediction_genome_counts)[1:23], names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  merge(metadata, by = "sampleid") 

host_prediction_genus_counts_metadata$Station <- factor(host_prediction_genus_counts_metadata$Station, levels = c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund"))
host_prediction_genome_counts_metadata$Station <- factor(host_prediction_genome_counts_metadata$Station, levels = c("Muhlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuttel","Meedem Grund"))

host_prediction_genus_counts_metadata  %>% 
  ggplot(aes(x=Station, y=counts, fill=Host.genus)) + geom_bar( stat = "identity") +
  geom_tile() + facet_wrap(~data_type *Sample_date) + ggtitle("Virus count of particular hosts")

host_prediction_genus_counts_metadata$Sample_date <- factor(host_prediction_genus_counts_metadata$Sample_date, levels= c("May-21", "Jul-21", "Nov-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
host_prediction_genome_counts_metadata$Sample_date <- factor(host_prediction_genome_counts_metadata$Sample_date, levels= c("May-21", "Jul-21", "Nov-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))

tiff("geneabund_drep_marker_iphop_genus_plots.tiff", unit="px", width = 700, height = 1000) 
host_prediction_genus_counts_metadata  %>% 
filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov-21") %>%
  filter(data_type != "METAT") %>%
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=Station,y=phylum, fill=log(counts))) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Copies\n/Transcripts per genome", y="phylum") +
  ggtitle("Hosts of viruses based on genus iPHoP prediction")
dev.off()

tiff("geneabund_drep_marker_iphop_genome_plots.tiff", unit="px", width = 700, height = 1000) 
host_prediction_genome_counts_metadata  %>% 
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov-21") %>%
  filter(data_type != "METAT") %>%
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=Station,y=phylum, fill=log(counts))) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Copies\n/Transcripts per genome", y="phylum") +
  ggtitle("Hosts of viruses based on genome iPHoP prediction")
dev.off()

# another plot to see how many generalist vs specialist
  host_prediction_genus_counts  %>%
  group_by(Virus) %>% 
  select(c("Virus", "primary_cluster.x", "Host.genus")) %>% 
  summarise(count_hosts = n_distinct(Host.genus))  %>% tail()
    
    
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/marker_genes_taxa.RData")
host_prediction_genus_counts_metadata <- host_prediction_genus_counts_metadata %>%
  mutate(sampleid=str_replace(sampleid, pattern="GROS22.[0-9]_", "")) %>%
  mutate(sampleid=str_replace(sampleid, ".genecount.profile", "")) %>%
  filter(Station != 'BunthausSpitze') %>% #remove 21 Nov and Buntspitze as we dont have this station for all
  filter(Sample_date != "Nov 21")

host_prediction_genome_counts_metadata <- host_prediction_genome_counts_metadata %>%
  mutate(sampleid=str_replace(sampleid, pattern="GROS22.[0-9]_", "")) %>%
  mutate(sampleid=str_replace(sampleid, ".genecount.profile", "")) %>%
  filter(Station != 'BunthausSpitze') %>% #remove 21 Nov and Buntspitze as we dont have this station for all
  filter(Sample_date != "Nov 21")
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/marker_genes_taxa.RData")
#lets save it now and continue in the Motus script - as it needs too much computational power
#we will filter out not neccessery phylla there which are subset of those from motus and iphop file
saveRDS(host_prediction_genome_counts_metadata, file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/host_prediction_genome_counts_metadata.RDS") 
saveRDS(host_prediction_genus_counts_metadata, file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/host_prediction_genus_counts_metadata.RDS") 

host_prediction_genus_counts_metadata_mantel <- host_prediction_genus_counts_metadata  %>%
  select(phylum, sampleid, counts) %>%
  group_by(phylum, sampleid) %>% 
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% 
  pivot_wider(names_from="sampleid", values_from="counts", values_fill=0) 

saveRDS(host_prediction_genus_counts_metadata_mantel, file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/host_prediction_genus_counts_metadata_mantel.RDS")
host_prediction_genome_counts_metadata_mantel = readRDS(file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/host_prediction_genome_counts_metadata_mantel.RDS")

#prepare matrix for mantel test 
host_prediction_genus_counts_metadata_mantel <- host_prediction_genus_counts_metadata[, c("phylum", "Virus", "counts", "sampleid")] %>%
  pivot_wider(names_from = sampleid, values_from = counts) 
#host_prediction_genus_counts_metadata_mantel <- as.dataframe(host_prediction_genus_counts_metadata_mantel)
#rownames(host_prediction_genus_counts_metadata_mantel) <- host_prediction_genus_counts_metadata_mantel[,1] 
host_prediction_genome_counts_metadata_mantel <- host_prediction_genome_counts_metadata[, c("phylum", "Virus", "counts", "sampleid")] %>%
  pivot_wider(names_from = sampleid, values_from = counts)

 host_prediction_genus_counts_mantel <- host_prediction_genus_counts %>%
  separate(Host.genus, c("domain", "phylum", "class", "order", "family",  "genus"), ";") %>%
  select(-c("domain", "class", "order", "family",  "genus"))  %>%  filter(!is.na(phylum)) %>% 
  group_by(phylum) %>% head()
  mutate_each(funs(sum)) %>% distinct()
host_prediction_genus_counts_mantel <- as.data.frame(host_prediction_genus_counts_mantel)
rownames(host_prediction_genus_counts_mantel) <- host_prediction_genus_counts_mantel$phylum

geneabund_drep_marker_taxa_mantel <- geneabund_drep_marker_taxa %>% filter(!is.na(phylum)) %>% group_by(phylum) %>%
  select(-c("genome", "gene_cluster", "KO","primary_cluster",  "secondary_cluster","DESCRIPTION", "count_primary_cluster", "count_genome", "length.x", "length.y", "domain", "class", "order", "family", "genus", "species", "gene"))  %>%  
  mutate_each(funs(sum)) %>% distinct()

geneabund_drep_marker_taxa_mantel <- as.data.frame(geneabund_drep_marker_taxa_mantel)
rownames(geneabund_drep_marker_taxa_mantel) <- geneabund_drep_marker_taxa_mantel$phylum
#merge prepared data frames
mantel_merge <- merge(host_prediction_genus_counts_mantel, geneabund_drep_marker_taxa_mantel, 
                           by = 'row.names', all = TRUE) 
#after merging take each part again and transform into matrix to have the same rows
matrix_iphop <- as.numeric(as.matrix(mantel_merge[, 2:length(colnames(host_prediction_genus_counts_mantel))+1]))
matrix_markers <-  as.numeric(as.matrix(mantel_merge[,(length(colnames(host_prediction_genus_counts_mantel))+2):(length(mantel_merge[1,])-1)]))
#mantel test between two differently obtained hosts tables 
library("vegan")
iphop.dist <- vegdist(matrix_iphop, "euclidean", na.rm = T)
markers.dist <- vegdist(matrix_markers, "euclidean", na.rm = T)

mantel(iphop.dist, markers.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
####
#to do mantel test - between viruses classified by hosts and hosts
#
library(vegan)

iphop.dist <- vegdist(host_prediction_genus_counts[, 17:length(host_prediction_genus_counts[1,])], "euclidean", na.rm = T)
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/marker_genes_taxa.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/marker_genes_taxa.RData")
#distance matrix viruses
#distance matrixes hosts - use 16s as a marker gene? do we have it? or also find individual genes
#mantel test for for matrixes - (use the same columns order)
#conclusion between driving diversity of viruses and hosts    
#for differently obtained hosts matrixes we can try - 1st distance between samples 
#then mantel test or just mantel test but  it have to be for non quatratic matrixes
#check if the same taxonomy host_prediction_genus and host_prediction_genome
#save for motus code analysis
saveRDS(host_prediction_genome, file = "host_prediction_genome.RDS")
host_prediction_genome$lineage_genus == host_prediction_genus$Host.genus 
length(host_prediction_genome$lineage_genus)
length(host_prediction_genus$Host.genus)
length(unique(host_prediction_genus$Host.genus))
length(unique(host_prediction_genome$lineage_genus))### I they are different genus use only hosts-based tools and genome use host-based and phage-based tools
#lets stick with host_prediction_genus as this was used before  
vibrant_annot_drep_marker %>% group_by(DESCRIPTION) %>%
  summarise(DESCRIPTION, counts = n()) %>% distinct() %>% 
  arrange(desc(counts)) %>% print(n = 25)

vir_genecluster_annot %>% filter(KO!="NA" & KO!="None") %>% filter(gene_cluster=="GROS22_PROK_MAG_VIR_0000717556") %>% head(n=50)
#vibrant_annot_drep_marker <- check this for vibrant annotation 
###vir_genecluster_annot %>% filter(KO!="NA" & KO!="None") %>% group_by(gene_cluster) %>% summarise(gene_cluster, counts=n_distinct(genome))  %>% filter(counts > 1) %>% arrange(desc(counts))
#vir_genecluster_annot %>% filter(KO!="NA" & KO!="None") %>% filter(gene_cluster=="GROS22_PROK_MAG_VIR_0000239754")
#check MCPs across all sammples - if therse is biases due to conserved sequences
MCP_list <- c("VOG00035", "VOG00461", "VOG01150", "VOG01000", "VOG01164", "VOG01920", "VOG02155", "VOG02220", "VOG02554", "VOG03097", "VOG03553", "VOG03780", "VOG03813", "VOG33139")
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