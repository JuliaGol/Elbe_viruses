library(tidyverse)

setwd("~/IGB_phd/virus/viral_abundances")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog"
##Load Kegg annot file and clstr  (skip if already loaded in get_data)
annot <- read.csv(paste0(path,"/GROS22.kegg-annotations.tsv"), sep="\t", skip=0, header = TRUE, row.names = 1, quote="")
clstr <- read.csv(paste0(path,"/GROS22.clstr"), sep="\t", skip=0, header = FALSE)
colnames(clstr) <- c("gene_cluster", "mem", "gene", "length", "pos")
clstr <- clstr %>% 
  select(gene_cluster, gene)


##Load Vibrant annotations, note annotation at gene not gene_cluster level so we need to resolve this.
vibrant_annot <- read.csv("BICEST.phages.annot", header = FALSE, sep = "\t")
colnames(vibrant_annot) <- c("gene", "range", "KO", "DESCRIPTION")

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


##We have around 8Mio genes total but some of those are unannotated so we dont use them here. 
##Fix gene name and remove genes without annotation
vibrant_clean <- vibrant_annot %>% 
  select(gene, KO, DESCRIPTION) %>% 
  mutate(gene = str_remove(gene, ">")) %>% 
  filter(KO != "None")
#Cleanup
rm(vibrant_annot) 
##We are down to around 2.7Mio annotation genes.
##Ok lets assign genes to gene clusters and get annotations for the gene clusters
head(vibrant_clean)
vir_gene_genecluster_annot <- vibrant_clean %>% 
  left_join(clstr)

##We actually have some genes annotated as KEGG terms so we can use this to add some info from the annot file (assuming its there, but it also adds the extra colums to make it easier)

annot <- annot %>% rownames_to_column("gene_cluster") #my edit

vir_kegg_annots <- vir_gene_genecluster_annot %>% 
  select(KO) %>%
  filter(grepl("K", KO)) %>%
  distinct() %>% 
 left_join(select(annot, !gene_cluster)) %>%
 distinct()
  
vir_gene_genecluster_annot_new <- vir_gene_genecluster_annot %>%
  left_join(vir_kegg_annots, by=c("KO", "DESCRIPTION"), keep=NULL) 
  
  head(vir_gene_genecluster_annot_new)
##We need to check to see if there is a single annotation or multiple different annotatations for each gene cluster. 
  
  vir_gene_genecluster_annot_new %>% 
  select(!gene) %>% 
  group_by(gene_cluster) %>% 
  summarise(count = n_distinct(KO)) %>% 
  filter(count > 1)

##We see that for 604K annotated gene_clusters we have around 26K clusters that have more than one distinct annotation 
 
##We can use the mutate and paste to concatenate the non distinct annotations
vir_genecluster_annot <- vir_gene_genecluster_annot_new %>% 
  select(!gene) %>%  
  distinct() %>%
  group_by(gene_cluster) %>% 
  mutate(across(everything(), ~ paste0(., collapse = ";")))

###my edit
vir_gene_genecluster_annot_concat <- vir_gene_genecluster_annot_new %>% 
  distinct() %>%
  group_by(gene_cluster) %>% 
  mutate(across(everything(), ~ paste0(., collapse = ";")))

saveRDS(vir_gene_genecluster_annot_concat, "vir_gene_genecluster_annot_concat.rds")

##Lets check again to make sure we have a single annotation for each cluster

vir_genecluster_annot %>% 
  group_by(gene_cluster) %>% 
  summarise(count = n_distinct(KO)) %>% 
  filter(count > 1)

##Great the result is empty so we can move on! 
##Now we need to format the vibrant annotations so that they match the kegg annotations
##Lets check what columns we have in the annot file

colnames(annot)

colnames(vir_genecluster_annot)
## ok thats a few, but gene cluster (at least for me is rownames, so I will fix)

#annot <- annot %>% rownames_to_column("gene_cluster") # I run it earlier

##Lets reorder the columns so that they match
vir_genecluster_annot <- vir_genecluster_annot %>%
  select(gene_cluster, everything())%>%
  relocate(DESCRIPTION, .after=everything())

saveRDS(vir_genecluster_annot, "vir_genecluster_annot.rds")

colnames(annot)
colnames(vir_genecluster_annot)

##Ok now we can see that the two match

annot_mags_vibrant <- rbind(annot, vir_genecluster_annot)
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")
##Lets check again for duplications because there were some annotations for viral genes already existing

annot_mags_vibrant %>% 
  group_by(gene_cluster) %>% 
  summarise(count = n_distinct(KO)) %>% 
  filter(count > 1)

##Yeah so there is an overlap of around 8k genes, from 1 Mio, lets fix this, 

annot_mags_vibrant <- annot_mags_vibrant %>% 
  group_by(gene_cluster) %>% 
  mutate(across(everything(), ~ paste0(., collapse = ";"))) #what does it? collapsing results?
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")
annot_mags_vibrant %>% 
  group_by(gene_cluster) %>% 
  summarise(count = n_distinct(KO)) %>% 
  filter(count > 1)

head(annot_mags_vibrant)
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")

##Great now there is a good chance that this isnt a perfect solution, but it should do what we want. 

##Lets check now to see how many gene clusters are annotated as capsid proteins, before this was 1

##VOGDB list of putative Major Capsid Protein Orthologous Groups. Some are specific others not. 
#MCP_list <- c("VOG00035", "VOG00461", "VOG01150", "VOG01000", "VOG01164", "VOG01920", "VOG02155", "VOG02220", "VOG02554", "VOG03097", "VOG03553", "VOG03780", "VOG03813", "VOG33139")
MCP_list <- c("VOG00035", "VOG01150")


annot_mags_vibrant %>% 
  ungroup() %>%
  filter(grepl(paste(MCP_list, collapse = "|"), KO)) %>% 
  group_by(KO) %>%
  summarise(count = n_distinct(gene_cluster)) %>%
  summarise(sum(count))

##So we can identify 2894 Unique MCP proteins, lets see if we can filter also with KEGG some false positives - Do you know how many dreplicated viral genomes?

MCP_summary <- annot_mags_vibrant %>% 
  ungroup() %>%
  filter(grepl(paste(MCP_list, collapse = "|"), KO)) %>% 
  group_by(KO, DESCRIPTION) %>%
  summarise(count = n_distinct(gene_cluster)) 

##Possible contaminations:::VOG00654;K03546 - DNA repair, K06881;VOG01164 - Phosphatase, K07313;VOG01150 - Phosphatase, K07491;VOG02220 - transposase, VOG00461;K03546 - DNA repair
##Putative hits:: K07097;VOG01164 - unknown protein, 

##Im not sure what to remove and its not clear why there are mismatches between VOGDB and KEGG with some proteins annotatied as exonuc, phosphatase and capsid protein. Check mVIR paper. 

##Lets continue assuming we have a list of gene_clusters, the rest of this pipeline essentially follows the pipeline that we used before for amoA 

BICEST_MCP_taxa_sample <- annot_mags_vibrant %>% 
  ungroup() %>%
  filter(grepl(paste(MCP_list, collapse = "|"), KO)) %>% 
  select(gene_cluster, KO) %>%
  left_join(clstr, relationship = "many-to-many") %>%
  select(KO, gene_cluster, gene) %>%
  mutate(genome=if_else(grepl("MAG", gene), str_remove(gene, "-scaff.*"), str_remove(gene, "_length.*"))) %>%
  mutate(sample=if_else(grepl("MAG", gene), str_remove(genome, "_MAG.*"), str_remove(genome, "_MET.*"))) %>%
  left_join(taxa, by=c("genome"="GENOME"), relationship = "many-to-many") %>%
  left_join(metadata, relationship = "many-to-many") %>%
  mutate(mag=if_else(grepl("Bacteria", domain), 1, 0)) %>% 
  distinct()

BICEST_MCP_taxa_sample_summary <- BICEST_MCP_taxa_sample %>%
  filter(KO == "VOG00035" | KO == "VOG01150") %>%
  group_by(gene_cluster, KO) %>%
  summarise(sample_count = n_distinct(sample), scaf_count = n_distinct(genome), MAG_hits = sum(mag)) %>%
  arrange(scaf_count)

sum(BICEST_MCP_taxa_sample_summary$scaf_count) - sum(BICEST_MCP_taxa_sample_summary$MAG_hits)
save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")

##The previous output gives us an overview for each gene_cluster (MCP), which proportion are found in the viral contigs and MAGs 
##There are a bunch of examples, so sometimes we find the capsid protein in a lot of MAGs but only on one scaffold (not assembled), 
##We find a lot of capsid proteins on scaffolds that we never find in MAGs
## There are also some cases where we find a bunch of scaffolds and they appear 1-2 times in a MAG. 
##Just be careful with not overintepreting these results on their own but use them also in combi with phist/iphop. 
##For Nitrospira and Methylopumilus where there seems to be quite a few MAGs containing capsid protein hits then it might be worthwhile to check these MAGs 
##and see if we can find other viral genes/scaffolds. (https://github.com/SushiLab/mVIRs) or by repeating this pipeline with other genes. 

geneabund <- geneabund %>% 
  rownames_to_column("gene_cluster")

BICEST_MCP_abund <- BICEST_MCP_taxa_sample_summary %>% 
  select(gene_cluster, KO) %>% 
  left_join(geneabund, by=c("gene_cluster")) %>%
  ungroup() %>%
  select(., KO, contains("GROS")) %>% 
  pivot_longer(!KO, names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  left_join(metadata) %>%
  filter(Sample_date != "Nov 21") %>%
  filter(Station != 'BunthausSpitze') %>% 
  filter(!grepl(";", KO)) ##Removes terms containing multiple VOG/KEGG terms
  
BICEST_MCP_abund %>% 
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=factor(Station, levels= c("Muehlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuettel","Meedem Grund")), 
                y=KO, fill=log(counts))) +
  facet_wrap(~ data_type *Sample_date ) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  labs(fill = "Copies/Transcripts per genome", y="Major Capsid Protein VOGs")
  
saveRDS(BICEST_MCP_abund, "BICEST_MCP_abund.rds")
#######By host
MCP_taxa <- BICEST_MCP_taxa_sample %>% 
  filter(grepl("MAG", gene)) %>%
  select(gene_cluster, KO, phylum)

BICEST_MCP_taxa_abund <- BICEST_MCP_taxa_sample_summary %>% 
  select(gene_cluster, KO) %>% 
  left_join(geneabund, by=c("gene_cluster")) %>%
  ungroup() %>%
  select(., gene_cluster, contains("GROS")) %>% 
  pivot_longer(!gene_cluster, names_to = "sampleid", values_to = "counts", values_drop_na = 0) %>%
  left_join(metadata) %>%
  left_join(MCP_taxa, by="gene_cluster", relationship = "many-to-many") %>%
  filter(Sample_date != "Nov 21") %>%
  filter(Station != 'BunthausSpitze') %>% 
  filter(!grepl(";", KO)) %>%##Removes terms containing multiple VOG/KEGG terms
  group_by(phylum) %>%
  mutate(total = sum(counts)) %>%
  filter(total > 550) %>%
  ungroup()
  

BICEST_MCP_taxa_abund %>% 
  mutate(counts=counts/1000) %>%
  ggplot() + 
  geom_tile(aes(x=factor(Station, levels= c("Muehlenberger Loch","Twielenfleth", "Schwarztonnensand", "Brunsbuettel","Meedem Grund")), 
                y=phylum, fill=log(counts))) +
  facet_wrap(~ data_type *Sample_date ) +
  scale_fill_gradient(low="darkblue", na.value = 'darkblue', high="red") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  labs(fill = "Copies/Transcripts per genome", y="Phylum")

save.image("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")
load("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/viral_abundance.RData")
