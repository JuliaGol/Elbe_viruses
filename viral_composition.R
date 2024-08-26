library("dplyr")
library("statip")
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viruses_composition")
virus_bicest_blast_taxonomy <- read.csv("virus_bicest_blast_taxonomy.tsv", sep="\t", header = F)
colnames(virus_bicest_blast_taxonomy) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "staxids", "sscinames", "scomnames", "sblastnames", "sskingdoms")
length(unique(virus_bicest_blast_taxonomy$qseqid)) # 62 why so little

#filter out e-value more then 0.001 
virus_bicest_blast_taxonomy <- virus_bicest_blast_taxonomy %>% filter(evalue<0.001)
virus_bicest_blast_taxonomy_grouped <- virus_bicest_blast_taxonomy %>% group_by(qseqid) %>%  summarise(mfv_kingdoms = mfv(sskingdoms)) 
#the most common 
virus_bicest_blast_taxonomy_grouped_sci <- virus_bicest_blast_taxonomy %>% group_by(qseqid) %>%  summarise(mfv_sscinames = mfv(sscinames)) 
#or use the one with the lowest e-value 
virus_bicest_blast_taxonomy_grouped_low_e_value <- virus_bicest_blast_taxonomy %>% group_by(qseqid) %>%  summarise(min_evalue = min(evalue), sscinames = sscinames) %>%
  group_by(qseqid) %>%  summarise(mfv_sscinames = mfv(sscinames)) #if multiple find the most common
#use the last one and the counts to plot a taxonomy
metadata <- read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog/SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAG.genecount.profile"))%>%
  mutate(data_type="METAG")
geneabund_drep_marker_taxa <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/geneabund_drep_marker_taxa.RDS")


sum(virus_bicest_blast_taxonomy_grouped[,"mfv_kingdoms"] == "Viruses")
sum(virus_bicest_blast_taxonomy_grouped[,"mfv_kingdoms"] == "Bacteria")
length(unique(virus_bicest_blast_taxonomy$qseqid)) # 43 after filtering out low evalue why so little
#very little was classified 
#grep ">" /scratch/jgolebiowska/viruses_BICEST/VIRUSES/drep_virus_all.fna  | wc -l
#2416
mfv(virus_bicest_blast_taxonomy$sskingdoms) #Bacteria