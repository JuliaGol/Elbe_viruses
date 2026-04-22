#setting
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("your/path/to/Elbe_viruses")
##Loading files
#load if you saved
#load("your/path/to/marker_genes.RData")
#gene_cluster_table
clstr <- read.csv("data/GROS22.clstr", sep="\t", skip=0, header = F)
colnames(clstr) <- c("gene_cluster", "mem", "gene", "length", "pos")
#viral genes annotated
vibrant_annot <- read.csv("data/BICEST.phages.annot",  sep = "\t", header = T)
colnames(vibrant_annot) <- c("genome", "range", "KO", "DESCRIPTION")
vibrant_annot[,"genome"] = sub(">", "", vibrant_annot[,"genome"])
#taxa assignment for MAGs
taxa <- read.csv("data/GROS22-2.gtdb", sep="\t", skip=0, header = TRUE)
taxa <- taxa %>%
  separate(GTDBTK_TAXONOMY, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";")
#abundance
geneabund <- read.csv("data/BICEST_cellabund.tsv", sep="\t", skip=1, header = TRUE, row.names = 1)

metadata <- read.csv("data/PhysicochemicalParameters_mod3.csv", header=TRUE, sep=",", row.names=1) 
#remove rows with NA in Acession number - no molecular data
metadata <- metadata[which(!is.na(metadata$AccessionNumber_TBDSven)),]
#unify stromkilometer
metadata$Stromkilometer <- round(metadata$Stromkilometer)
metadata$Stromkilometer[which(metadata$Stromkilometer==715)] <- 712 
metadata$Stromkilometer[which(metadata$Stromkilometer==694)] <- 692 
metadata$Stromkilometer[which(metadata$Stromkilometer==666)] <- 665 
metadata$Stromkilometer[which(metadata$Stromkilometer==652)] <- 651
metadata$Stromkilometer <- as.factor(metadata$Stromkilometer)  
#drep clusters for finding gene markers 
drep_clusters <- read.csv("data/Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])
#read combine quality file 
#check how many of each quality 
quality <- read.csv("data/VIBRANT/Combined_quality_summary.tsv", header = TRUE, sep = "\t")
drep_clusters <-merge(drep_clusters, quality[, c("contig_id", "checkv_quality")], by.x = "genome", by.y = "contig_id")
#merge the vibrant annotation with gene clusters
vir_genecluster_annot <- vibrant_annot %>% 
  left_join(clstr,  join_by(genome == gene)) #st wrong here - many NAs in final table
#modify the genome column for merging with drep clusters 
vir_genecluster_annot[,"genome"] = sub("_[1-9]+$", "", vir_genecluster_annot[,"genome"])
vibrant_genecluster_drep = merge(drep_clusters, vir_genecluster_annot, by="genome")  %>%
  select(genome, gene_cluster, KO, primary_cluster, secondary_cluster, DESCRIPTION)
vibrant_genecluster_drep_abund <- left_join(vibrant_genecluster_drep, geneabund, by = c("gene_cluster" = rownames(geneabund))) 
#check how many drep_clusters per one gene cluster
vibrant_genecluster_drep_abund_primary_counts <- vibrant_genecluster_drep_abund %>%
  group_by(primary_cluster) %>% summarise(n = n_distinct(gene_cluster))
#how many more then one?
vibrant_genecluster_drep_abund_primary_counts[which(vibrant_genecluster_drep_abund_primary_counts$n >1), ] %>% nrow()
#how many exactly one
vibrant_genecluster_drep_abund_primary_counts[which(vibrant_genecluster_drep_abund_primary_counts$n == 1),] %>% nrow()
#the same for secondary clusters 
vibrant_genecluster_drep_abund_secondary_counts <- vibrant_genecluster_drep_abund %>%
  group_by(secondary_cluster) %>% summarise(n = n_distinct(gene_cluster))
vibrant_genecluster_drep_abund_secondary_counts[which(vibrant_genecluster_drep_abund_secondary_counts$n >1), ] %>% nrow()
#how many exactly one
vibrant_genecluster_drep_abund_secondary_counts[which(vibrant_genecluster_drep_abund_secondary_counts$n <= 1),] %>% nrow()

#find primary clusters using all KO
vibrant_annot_drep_marker <- vibrant_genecluster_drep %>%
  filter(KO != "None") %>% 
  group_by(gene_cluster) %>% 
  mutate(count_primary_cluster = n_distinct(primary_cluster)) %>%
  filter(count_primary_cluster == 1)  

vibrant_annot_drep_marker_gene_cluster <- vibrant_annot_drep_marker[,"gene_cluster"]  
#how many primary clusters has marker genes 
n_distinct(vibrant_annot_drep_marker[, "primary_cluster"])/max(vibrant_annot_drep_marker[, "primary_cluster"])  ##99,54% of clusters covered!!! 
#max(vibrant_genecluster_drep[, "primary_cluster"])
n_distinct(vibrant_annot_drep_marker[,"genome"])/n_distinct(vibrant_genecluster_drep[,"genome"]) #not correct because it can be that 
#99,58% genomes!!! -gene wasn't detected at some genomes but were assigned to cluster
#not all genomes in vibrant_annot_drep_marker - keep only valid columns for clarity 
#save only first row from each cluster 
vibrant_annot_drep_marker <- vibrant_annot_drep_marker %>%
  group_by(primary_cluster) %>%
  mutate(count_genome = n_distinct(genome)) %>% #count genomes per cluster (not per gene) for statistics
  arrange(desc(KO)) %>% #arrange KO in descending order so we will favor VOG 
  filter(row_number()==1) %>%
  select(!c("genome")) #remove genome column as this is table summarising clusters and their marker genes 
#statistics 
n_distinct(vibrant_annot_drep_marker[, "primary_cluster"])/max(vibrant_annot_drep_marker[, "primary_cluster"]) #max(vibrant_genecluster_drep[, "primary_cluster"]) ##99,54% of clusters covered!!! 
sum(vibrant_annot_drep_marker[,"count_genome"])/length(unique(drep_clusters$genome)) #not correct because it can be that some genomes were filter out together 
saveRDS(vibrant_annot_drep_marker, file = "vibrant_annot_drep_marker.RDS")
#let make a list of gene_clusters which marker genes belong to 
gene_cluster_marker <- vibrant_annot_drep_marker$gene_cluster
#summarise new gene counts list
vibrant_annot_drep_marker %>% group_by(DESCRIPTION) %>% summarise(n = n())  %>% arrange(desc(n))

#redo host analysis similarly to the one based on MCP
#merging the current vibrant_annot_drep_marker dataframe with data about gene clusters abundance 
geneabund_marker_gene_cluster <- geneabund[which(geneabund$gene_cluster %in% gene_cluster_marker),]
vibrant_annot_drep_marker_abund <- vibrant_annot_drep_marker %>% left_join(geneabund, by = "gene_cluster")
saveRDS(vibrant_annot_drep_marker_abund, file ="vibrant_annot_drep_marker_abund.RDS")


###interfere with Iphop data
host_prediction_genome = read.csv("data/Host_prediction_to_genome_m90.csv")
host_prediction_genus = read.csv("data/Host_prediction_to_genus_m90.csv")
#load eukaryotic host prediction
eukryotic_host_prediction_genome = read.csv("data/eukvirus_hits_with_tax.tsv", sep = '\t')
##add cluster number to each virus
host_prediction_genus <- merge(host_prediction_genus, drep_clusters,  by.x = "Virus", by.y = "genome")
#the same with genome file - hosts obtained differently 
host_prediction_genome <- merge(host_prediction_genome, drep_clusters,  by.x = "Virus", by.y = "genome")
#the same for eukaryotic
eukryotic_host_prediction_genome <- merge(eukryotic_host_prediction_genome, drep_clusters,  by.x = "query_full", by.y = "genome")
####select phages
#load  taxonomy information to filter out eukaryotic viruses - 
geneabund_drep_marker_taxonomy_long <- readRDS(file="data/geneabund_drep_marker_taxonomy_long.RDS")
#add cluster number to each virus
host_prediction_genus <- merge(host_prediction_genus, geneabund_drep_marker_taxonomy_long[,c("level_3", "primary_cluster")],  by = "primary_cluster")
#the same with genome file - hosts obtained differently 
host_prediction_genome <- merge(host_prediction_genome, geneabund_drep_marker_taxonomy_long[,c("level_3", "primary_cluster")],  by = "primary_cluster")

#prepare lineage to genus only host_prediction_genome and division to different taxa in the column
host_prediction_genus_phages <- host_prediction_genus %>% 
  separate(Host.genus, c("domain", "phylum", "class", "order", "family",  "genus"), ";") %>%
  filter(host_prediction_genus$level_3 == "Heunggongvirae")
host_prediction_genome_phages <- host_prediction_genome %>% 
  separate(Host.taxonomy, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";") %>%
  filter(host_prediction_genome$level_3== "Heunggongvirae")
#the same for eukaryotic virus
eukryotic_host_prediction_genome <- eukryotic_host_prediction_genome %>% 
  separate(host_lineage, c("host.level_1", "host.level_2", "host.level_3", "host.level_4", "host.level_5",  "host.level_6", "host.level_7", "host.level_8", "host.level_9", "host.level_10"), "; ", remove=F) 
#keep euk hosts
eukryotic_host_prediction_genome <- eukryotic_host_prediction_genome %>% filter(host.level_2=="Eukaryota")
#compare virus prediction
eukryotic_host_prediction_genome <- eukryotic_host_prediction_genome %>% 
  separate(virus_lineage, c("virus.domain", "virus.phylum", "virus.class", "virus.order", "virus.family",  "virus.genus", "virus.species"), "; ", remove=F)

#paste lineage to genus level in order to compare
host_prediction_genome_phages$lineage_genus <- paste(host_prediction_genome_phages$domain, host_prediction_genome_phages$phylum, host_prediction_genome_phages$class, host_prediction_genome_phages$order, host_prediction_genome_phages$family,  host_prediction_genome_phages$genus, sep=";")

#how many host carry each cluster
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts == 1)  %>% nrow
#496
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(Host.taxonomy)) %>%
  filter(count_hosts == 1)  %>% nrow
#407
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  mutate(species=case_when(
    species=="s__" ~ NA_character_,
    TRUE  ~ species)) %>% 
  summarise(count_hosts = n_distinct(species, na.rm=T)) %>%
  filter(count_hosts == 1)  %>% nrow
#155
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(phylum)) %>%
  filter(count_hosts == 1)  %>% nrow
#533
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(phylum)) %>%
  filter(count_hosts > 1)  %>% nrow
#122
host_prediction_genus_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts == 1)  %>% nrow
#374
host_prediction_genus_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts > 1)  %>% nrow
#40  
host_prediction_genus_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(phylum)) %>%
  filter(count_hosts == 1)  %>% nrow
#388
host_prediction_genome_phages %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts == 0)  %>% nrow
#0
#for eukaryotic 
#how many host carry each cluster

eukryotic_host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(host_lineage)) %>%
  filter(count_hosts == 1)  %>% nrow
#9 
eukryotic_host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(host.level_10)) %>%
  filter(count_hosts == 1)  %>% nrow
#11
eukryotic_host_prediction_genome %>% 
  group_by(primary_cluster) %>% 
  summarise(count_hosts = n_distinct(host_lineage)) %>%
  filter(count_hosts > 1)  %>% nrow
#14

#check how many host in primary cluster now (some are filtered out as they did not have marker gene)
host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "genus"))  %>%
  summarise(count_hosts = n_distinct(genus)) %>% 
  filter(count_hosts > 1) %>%
  nrow
#40

host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "genus"))  %>%
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts == 1) %>%
  nrow
#374
host_prediction_genus_counts  %>%
  group_by(primary_cluster) %>% 
  select(c("Virus", "primary_cluster", "genus"))  %>%
  summarise(count_hosts = n_distinct(genus)) %>%
  filter(count_hosts < 1) %>%
  nrow
n_distinct(host_prediction_genus_counts[, "primary_cluster"]) 
#414 the same as above

#merge dataframes to get the virus genome name (Vibrant) and counts of marker gene in the same table
vibrant_annot_drep_marker_abund <- left_join(vibrant_annot_drep_marker, geneabund_marker_gene_cluster, by = "gene_cluster")

host_prediction_genus_counts <- left_join(host_prediction_genus_phages, vibrant_annot_drep_marker_abund[,!colnames(vibrant_annot_drep_marker_abund) %in% c("secondary_cluster", "genome", "count_primary_cluster")], by="primary_cluster")
#the same for genome file
host_prediction_genome_counts <- left_join(host_prediction_genome_phages, vibrant_annot_drep_marker_abund[,!colnames(vibrant_annot_drep_marker_abund) %in% c("secondary_cluster", "genome", "count_primary_cluster")], by="primary_cluster")
#and eukaryotic viruses
eukryotic_host_prediction_genome_counts <- left_join(eukryotic_host_prediction_genome, vibrant_annot_drep_marker_abund[,!colnames(vibrant_annot_drep_marker_abund) %in% c("secondary_cluster", "genome", "count_primary_cluster")], by="primary_cluster")

# how to deal in the plot with different hosts coming from the same virus
# focus on virus, so use counts even they are redundant, as we want to know how many viruses have particular host, not how many viruses are in total we want to see dynamic
host_prediction_genus_counts_metadata <- host_prediction_genus_counts %>%
  pivot_longer(col = starts_with("GROS"), names_to = "sampleid", values_to = "counts", values_drop_na = T) %>%
  merge(metadata, by = "sampleid")
# the same for genome counts
host_prediction_genome_counts_metadata <- host_prediction_genome_counts %>% 
  select(-c("class", "order", "family", "genus", "lineage_genus")) %>% distinct() %>% 
  pivot_longer(col = starts_with("GROS"), names_to = "sampleid", values_to = "counts", values_drop_na = T) %>%
  merge(metadata, by = "sampleid")
#and eukaryotes
eukryotic_host_prediction_genome_counts_metadata <- eukryotic_host_prediction_genome_counts %>%
  pivot_longer(col = starts_with("GROS"), names_to = "sampleid", values_to = "counts", values_drop_na = T) %>%
  merge(metadata, by = "sampleid")

#check how many clusters have different phyla as hosts (multiple hosts phyla)

host_prediction_genus_counts %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 1) %>%
  nrow                   
#388
host_prediction_genus_counts %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum > 1) %>%
  nrow  
#26
host_prediction_genus_counts %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 0) %>%
  nrow
#keep only distinct phyla per primary cluster 

host_prediction_genus_counts_metadata <- host_prediction_genus_counts_metadata  %>% 
  select(c("phylum", "primary_cluster", "counts", "Station", "Stromkilometer", "Sample_date", "data_type")) %>% 
  filter(Station != 'BunthausSpitze') %>% 
  filter(Sample_date != "Nov-21") %>%
  droplevels() %>%
  #filter(data_type != "METAT") %>%
  mutate(counts=counts/1000) %>%
  mutate(phylum=str_remove(phylum, "p__")) %>% 
  distinct() 
host_prediction_genus_counts_metadata$Sample_date <- factor(host_prediction_genus_counts_metadata$Sample_date, levels=c("May-21", "Jul-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
tiff("geneabund_drep_marker_iphop_genus_plots.tiff", unit="cm", width = 20, height = 30, res=300) 
host_prediction_genus_counts_metadata  %>%
  ggplot() + 
  geom_tile(aes(x=as.factor(Stromkilometer), y=phylum, fill=log(counts)), na.rm = T) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_text(), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Log Copies\n/Transcripts per genome", y="phylum", x="Elbe kilometer") +
  ggtitle("Hosts of viruses based on genus iPHoP prediction") #+ 
 # scale_x_discrete(labels = c("633", "651", "665", "692", "712"), name="Elbe kilometer")
dev.off()
host_prediction_genome_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum == 1) %>%
  nrow                   
#533
host_prediction_genome_counts_metadata %>%
  group_by(primary_cluster) %>% 
  select(c("primary_cluster", "phylum"))  %>%
  summarise(count_phylum = n_distinct(phylum)) %>%
  filter(count_phylum > 1) %>%
  nrow  
#122
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
  droplevels() %>%
  mutate(counts=counts/1000) %>%
  mutate(phylum=str_remove(phylum, "p__")) %>% distinct()
  
host_prediction_genome_counts_metadata$Sample_date <- factor(host_prediction_genome_counts_metadata$Sample_date, levels=c("May-21", "Jul-21", "Feb-22", "May-22", "Jun-22", "Nov-22"))
host_prediction_genome_counts_metadata$phylum[which(host_prediction_genome_counts_metadata$phylum=="")] <- "No Prediction"
ordered_phylum <-sort(unique(host_prediction_genome_counts_metadata$phylum[which(host_prediction_genome_counts_metadata$phylum!="No Prediction")]))
host_prediction_genome_counts_metadata$phylum <- factor(host_prediction_genome_counts_metadata$phylum, levels = c(ordered_phylum, "No Prediction"))
tiff("geneabund_drep_marker_iphop_genome_plots_phylum.tiff", unit="cm", width = 20, height = 30, res=300) 
host_prediction_genome_counts_metadata  %>%
  ggplot() +
  geom_tile(aes(x=as.factor(Stromkilometer),y=phylum, fill=log(counts))) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_text("Elbe kilometer"), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Log Copies\n/Transcripts per genome", y="phylum", x="Elbe kilometer") +
  ggtitle("Hosts of viruses based on genome iPHoP prediction")  #+
 # scale_x_discrete(labels = c("633", "651", "665", "692", "712"), name="Elbe kilometer")
dev.off()

# another plot to see how many generalist vs specialist

eukryotic_host_prediction_genome_counts_metadata %>% 
    group_by(primary_cluster) %>% 
    summarise(count_hosts = n_distinct(host_lineage)) %>%
    filter(count_hosts == 1)  %>% nrow
  #9 
eukryotic_host_prediction_genome_counts_metadata %>% 
    group_by(primary_cluster) %>% 
    summarise(count_hosts = n_distinct(host.level_10)) %>%
    filter(count_hosts == 1)  %>% nrow
  #11
eukryotic_host_prediction_genome_counts_metadata %>% 
    group_by(primary_cluster) %>% 
    summarise(count_hosts = n_distinct(host_lineage)) %>%
    filter(count_hosts > 1)  %>% nrow
  #14
  
#eukaryotes
tiff("geneabund_drep_marker_iphop_euk_plots_taxon8.tiff", unit="cm", width = 20, height = 30, res=300)
eukryotic_host_prediction_genome_counts_metadata  %>% filter(!is.na(host.level_8)) %>%
  ggplot() +
  geom_tile(aes(x=Stromkilometer,y=host.level_8, fill=log(counts))) +
  facet_wrap(~data_type *Sample_date) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_text("Elbe kilometer"), axis.text.x = element_text(angle=45, hjust=1), axis.text.y = element_text(size=12)) +
  labs(fill = "Log Copies\n/Transcripts per genome", y="taxon", x="Elbe kilometer") +
  ggtitle("Hosts of eukaryotic viruses")
dev.off()