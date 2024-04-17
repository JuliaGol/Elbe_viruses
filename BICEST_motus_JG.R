rm(list=ls())

require(RCurl)
require(tidyverse)
require(R.utils)
require(phyloseq)
require(ggplot2)
require(plotly)

path <- getwd()

h <- curl::new_handle()
curl::handle_setopt(
  handle = h,
  httpauth = 1,
  userpwd = "BICEST:vJGDYLs7"
)

##Download Data##
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/MOTUs")
url_new_motus <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/motus/mag_2_new_motus" 
url_ref_motus <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/motus/mag_2_existing_motus" 
url_tax_bac <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/gtdb/gtdbtk.bac120.summary.R214.tsv"
url_tax_arc <-"https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/gtdb/gtdbtk.ar53.summary.R214.tsv"
url_tax_motus <- "https://zenodo.org/records/10275750/files/mOTUs3.1.0.genome_metadata.tsv.gz?download=1"
url_checkm <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/metadata/GROS22-2.prok.genomes_2_quality.tsv.gz"
url_checkm2 <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/metadata/GROS22-1.prok.genomes_2_quality.tsv.gz"
motus_pkg_url <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/motus/motus/GROS22_mOTUs3.1.grosextended.motus"
drep_url <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/drep/drep-099/data_tables/Cdb.csv.gz"
url_16S_1 <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/metadata/GROS22-1.prok.barrnap.0.9.tsv.gz"
url_16S_2 <- "https://sunagawalab.ethz.ch/share/BICEST/GROS22/MAGs/metadata/GROS22-2.prok.barrnap.0.9.tsv.gz"

curl::curl_download(url_new_motus,paste(path,"/BICEST_mag_2_new_motus", sep=""), handle =  h)
curl::curl_download(url_ref_motus,paste(path,"/BICEST_mag_2_existing_motus", sep=""), handle =  h)
curl::curl_download(url_tax_bac,paste(path,"/BICEST_gtdbtk.bac120.summary.R214.tsv", sep=""), handle =  h)
curl::curl_download(url_tax_arc,paste(path,"/BICEST_gtdbtk.ar53.summary.R214.tsv", sep=""), handle =  h)
curl::curl_download(url_tax_motus,paste(path,"/mOTUs3.1.0.genome_metadata.tsv.gz", sep=""), handle =  h)
curl::curl_download(url_checkm,paste(path,"/GROS22-2.prok.genomes_2_quality.tsv.gz", sep=""), handle =  h)
curl::curl_download(url_checkm2,paste(path,"/GROS22-1.prok.genomes_2_quality.tsv.gz", sep=""), handle =  h)
curl::curl_download(motus_pkg_url,paste(path,"/GROS22_mOTUs3.1.grosextended.motus", sep=""), handle =  h)
curl::curl_download(drep_url,paste(path,"/Cdb.csv.gz", sep=""), handle =  h)
curl::curl_download(url_16S_1,paste(path,"/GROS22-1.prok.barrnap.0.9.tsv.gz", sep=""), handle =  h)
curl::curl_download(url_16S_2,paste(path,"/GROS22-2.prok.barrnap.0.9.tsv.gz", sep=""), handle =  h)

gunzip(paste0(path,"/mOTUs3.1.0.genome_metadata.tsv.gz"))
gunzip(paste0(path,"/Cdb.csv.gz"))
gunzip(paste0(path,"/GROS22-2.prok.genomes_2_quality.tsv.gz"))
gunzip(paste0(path,"/GROS22-1.prok.genomes_2_quality.tsv.gz"))
gunzip(paste0(path,"/GROS22-1.prok.barrnap.0.9.tsv.gz"))
gunzip(paste0(path,"/GROS22-2.prok.barrnap.0.9.tsv.gz"))


##Read mOTUs##

motus <- read.csv(paste0(path,"/GROS22_mOTUs3.1.grosextended.motus"), sep="\t", skip=2, header = TRUE, row.names = 1) ##Theres some junk in the first two lines
colnames(motus) <- gsub(".*_S", "S", colnames(motus) ) ##here im just removing the projectid for whatever reason. If you want to add this to gene_cat or remove here
rownames(motus) <- gsub(".* ", "", gsub("\\[|\\]", "", rownames(motus))) ##mOTUs contain a taxonomy followed by [mOTU_id] i just want mOTU_id

##Now we want to create  a list which matches each genome to its mOTU 
membership.new = read.csv(paste0(path, "/BICEST_mag_2_new_motus"), sep = "\t", header = FALSE) %>%
  rename(mOTU = V1, reps=V2) %>%
  separate_wider_delim(reps, ";", names_sep = "split", too_few=c("align_start")) %>% 
  pivot_longer(!mOTU, values_to = "user_genome",values_drop_na = TRUE) %>% 
  select(user_genome, mOTU)
membership.ref = read.csv(paste0(path, "/BICEST_mag_2_existing_motus"), sep = "\t", header = FALSE, quote="") %>% 
  select(V1, V2)%>%
  rename(., user_genome=V1, mOTU=V2)
membership <- rbind(membership.ref, membership.new)
nrow(membership)## 13765
##Read mOTUs Taxonomic Data

bac_gtdb = read.csv(paste0(path, "/BICEST_gtdbtk.bac120.summary.R214.tsv"), sep = "\t", header = TRUE) 
arc_gtdb = read.csv(paste0(path, "/BICEST_gtdbtk.ar53.summary.R214.tsv"), sep = "\t", header = TRUE)
gtdb = rbind(bac_gtdb, arc_gtdb) %>% 
  select(user_genome, classification )##Read in both archaea and bacteria
nrow(gtdb)## 13765
###Read in dRep data 
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances")
drep <- read.csv("Cdb.csv") %>% 
  select(genome, secondary_cluster)%>%
  rename(user_genome=genome, drep_cluster = secondary_cluster) %>% 
  mutate(user_genome = gsub(".fa", "", user_genome))
nrow(drep) ## 13765
##Read in qa data 

qa <- read.csv(paste0(path,"/GROS22-2.prok.genomes_2_quality.tsv"), sep="\t", quote="") %>% 
  rbind(., read.csv(paste0(path,"/GROS22-1.prok.genomes_2_quality.tsv"), sep="\t", quote="")) %>% 
  rename(user_genome = GENOME, completeness=COMPLETENESS, contamination=CONTAMINATION, genome_size=GENOME_SIZE) %>% 
  select(user_genome, completeness, contamination, N50, genome_size)
nrow(qa) ## 13765


##Read in barrnap (16S data)

barrnap <- read.csv(paste0(path, "/GROS22-1.prok.barrnap.0.9.tsv"), sep="\t", quote="") %>% 
  rbind(., read.csv(paste0(path, "/GROS22-2.prok.barrnap.0.9.tsv"), sep="\t", quote="")) %>% 
  rename(user_genome = genome, full_16S_rRNA = X16S_rRNA, partial_16S_rRNA = X16S_rRNA_partial) %>%
  select(user_genome, full_16S_rRNA, partial_16S_rRNA)

###Combine all together
genome_stats <- left_join(membership, gtdb) %>% 
  left_join(drep) %>% 
  left_join(qa) %>%
  left_join(barrnap)
nrow(genome_stats) ## 13765

######

###Stats####
motu_stats <- genome_stats %>%
  mutate(type = if_else(grepl("NotEnoughMGs", mOTU), "Unassigned", if_else(grepl("ELB", mOTU), "BICEST", "mOTUs.3.1.0"))) 

###If you want to you can use this opportunity to reassign some of your MAGs from Unassigned using dRep clusters###
motu_stats %>% 
  group_by(drep_cluster, mOTU) %>%
  summarise(n_mags = n_distinct(user_genome))

##You can already see ere that for instance 1000_1 has a MAG annotated as NotEnoughMGs, but this probably is ELBEEST_794
##However 1002_1  has 6 MAgs that are NotEnoughMGs, but then 2 that are ELBEEST_87 and ELBEEST_88. It would be nice to know which.
##This can also work in reverse, if we write this to output and search for ext_mOTU_v31_19264, we find around 10 drep clusters, 
##one of these has 141 and the rest 1 MAG each. Here I would tend towards assuming that each of these should be associated to this single mOTU
##Like i said before bad quality MAGs might not cluster with dREp but we find the MGs.
##Just be careful, if you want to modify the easiest way would be
motu_stats %>% 
  group_by(mOTU) %>% 
  summarise(n_drep = n_distinct(drep_cluster)) %>%
  filter(n_drep >1) %>% 
  arrange(desc(n_drep))

##This gives you a list of all the mOTUs with lots of drep clusters, 

##Double check the data 
motu_stats[which(motu_stats$mOTU == "ELBEEST_35"),] %>% 
  select(mOTU, drep_cluster, user_genome) %>% 
  group_by(drep_cluster) %>% 
  summarise(n_mags = n_distinct(user_genome)) %>% 
  arrange(desc(n_mags))

##Now a quick one liner will rename all the drep_clusters to match the major mOTU cluster
motu_stats[which(motu_stats$mOTU == "ELBEEST_35"), "drep_cluster"] <- "568_8" ##This will remove the additional drep clusters
##Now just rerun the code block above to see that it worked and move on to the next one. 

##Now thats sorted lets see about those NotEnoughMG annots, here we are only interested in those that have 2 unique mOTU, ie un-/classified

motu_stats %>% 
  group_by(drep_cluster, mOTU) %>%
  summarise(n_mags = n_distinct(user_genome)) %>% 
  ungroup() %>%
  group_by(drep_cluster) %>% 
  summarise(n_mOTU = n_distinct(mOTU), n_mags=sum(n_mags)) %>%
  filter(n_mOTU ==2) %>% 
  arrange(desc(n_mags))

##This gives you a list of all the drep_clusters with lots of mOTUs, 

##Double check the data 
motu_stats[which(motu_stats$drep_cluster == "616_1"),] %>% 
  select(mOTU, drep_cluster, user_genome) %>% 
  group_by(mOTU) %>% 
  summarise(n_mags = n_distinct(user_genome)) %>% 
  arrange(desc(n_mags))

##So we can see here that there are 134 assigned to ELBEEST_33 and 3 unassigned, lets go ahead and assign these. 
##Again a quick one liner
motu_stats[which(motu_stats$drep_cluster == "616_1"), "mOTU"] <- "ELBEEST_33" ##This will remove the additional drep clusters
##rerun the chunk above and we can see that now all are assigned to ELBEEST_33. 

###Lets plot the data to see how its looking. 
###We dont want to plot all groups so lets see which ones to focus on 
motu_stats %>% 
  separate(classification, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";") %>%
  group_by(phylum) %>%
  summarise(n_mags = n_distinct(user_genome)) %>% 
  arrange(desc(n_mags))

##Lets select the top 7 to plot and assign the rest to the "Other" bin

motu_stats %>% 
  separate(classification, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";") %>% 
  mutate(Phylum_plot = ifelse(grepl("p__Pseudomonadota|p__Actinomycetota|p__Bacteroidota|p__Verrucomicrobiota|p__Patescibacteria|p__Planctomycetota|p__Chloroflexota", phylum), phylum, "Other")) %>%
  mutate(Phylum_plot = gsub("p__", "", Phylum_plot)) %>%
  group_by(type, Phylum_plot) %>%
  summarise(n_mags = n_distinct(user_genome), PRESENT = sum(partial_16S_rRNA, full_16S_rRNA), ABSENT = n_mags-PRESENT) %>% 
  pivot_longer(cols = c("PRESENT", "ABSENT"), names_to = "rRNA", values_to ="counts") %>% 
  ggplot(aes(y=factor(type, levels=c("Unassigned", "mOTUs.3.1.0", "BICEST")), 
             x = n_mags, 
             alpha = factor(rRNA, levels=c("PRESENT","ABSENT")),
             fill=factor(type, levels=c("Unassigned", "mOTUs.3.1.0", "BICEST")))) + 
  geom_col() + scale_alpha_manual(values = c(1,0.5)) +
  theme_bw() +
  facet_grid(rows=vars(factor(Phylum_plot, levels=c("Actinomycetota", "Bacteroidota", "Chloroflexota", "Patescibacteria", "Planctomycetota",  "Pseudomonadota", "Verrucomicrobiota", "Other"))), 
             drop = TRUE, 
             scales="free", 
             switch = "y", 
             as.table = TRUE
  ) + xlab(label = "MAGs with and without partial 16S rRNA") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.text = element_text(angle = 90),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0,'lines'),
        panel.border = element_blank()) + 
  labs (fill= "Database Sournce", alpha = "Contains full/partial 16S rRNA") +
  guides(fill = guide_legend(order=2, ncol=3, title.position = "top", 
                             label.hjust = 0, title.vjust = 1,
                             keywidth = unit(0.5, 'cm'), override.aes=list(size = 1)),
         alpha = guide_legend(order=2, ncol=3, title.position = "top", 
                              label.hjust = 0, title.vjust = 1,
                              keywidth = unit(0.5, 'cm')))

#####Ok so we have all that information together. Use the above examples also if you want to go ahead and plot for instance genome completeness, size etc

##Now lets set ourselves up for some ecology!!! woohoo

##First we need to create a single taxonomy for each of the new ELBEEST/BICEST mOTUs

genome_stats %>% 
  select(mOTU) %>% 
  unique() %>% 
  nrow() ## So I get 1172 unique mOTUs, including 874 ELBEEST mOTUs. Depending on your cleaning steps (drep etc) the following might change a bit
##particularly when we look at unique classifications. Lets keep that in mind and go forward

genome_stats %>% 
  select(mOTU, classification) %>% 
  unique() %>% 
  nrow() ## 1946, so now we have added classification and so we have a few hundred examples where there is more than one classification for a mOTU

genome_stats %>% 
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  arrange(desc(count))

##Ok that was a false alarm as we have 732 classifications for unclassified mOTUs

genome_stats %>% 
  filter(mOTU != "NotEnoughMGs") %>%
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  arrange(desc(count)) %>% 
  nrow() ##ok so we need to resolve 40, this might differ when you add unclassified mOTUs into mOTU clusters above. 

##Ok so lets go through and see if we can do this relatively quickly 

motu_to_check <- genome_stats %>% 
  filter(mOTU != "NotEnoughMGs") %>%
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  select(mOTU) %>% 
  pull()

##Ok so we have our list of 40 mOTUs that are conflicting

genome_stats %>% 
  filter(mOTU %in% motu_to_check) %>% 
  select(mOTU, classification, user_genome) %>% 
  group_by(mOTU, classification) %>%
  summarise(count = n_distinct(user_genome)) %>%
  arrange(mOTU) %>% 
  print(n=100)
  
##Here we see at a first look that usually there is a more complete and less complete taxonomy. 
##Wracking my brain to automate this but I think we need to do a manual curate so lets set that up. 

##Just a quick note, how do we make sense of different level tax assignments. Lets say there are 10 MAGs belonging to a mOTU and 9 are assigned to species and 1 not. 
##This is an easy fix, ok lets just take the majority and go to species level. 
##Now what about if we have 9 not assigned and 1 is. Well this might be because only 1 was HQ and able to be assigned to species level. Shouldn't we trust the info from our HQ mag. I say yes!
##Of course with 10 it would be better if a few MAGs were assigned to species not just 1. 
##Generally most of your analysis will occur at the mOTU level or higher taxonomy (family or even phylum ) so dont stress about the need to assign each mOTU to species level. 

##You will also get mOTUs that are assigned to the same species but are two distinct mOTUs, so its not all that simple.

##code hidden here
genome_stats[which(genome_stats$mOTU == "ELBEEST_124" & genome_stats$classification == "d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__Nitrospiraceae;g__Nitrospira_F;s__"), "classification"] <- "d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__Nitrospiraceae;g__Nitrospira_F;s__Nitrospira_F sp919902665"
#The next example has two genus annotations for one mOTU, but the major is assigned to a single so we take that one 
genome_stats[which(genome_stats$mOTU == "ELBEEST_125" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__SZUA-36;g__;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__SZUA-36;g__JABDPF01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_125" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__SZUA-36;g__JAHEFT01;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales;f__SZUA-36;g__JABDPF01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_174" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__Halioglobus;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_193" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__UBA3439;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__UBA3439;s__UBA3439 sp016778825"
genome_stats[which(genome_stats$mOTU == "ELBEEST_194" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__UBA3439;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__UBA3439;s__UBA3439 sp016778825"
genome_stats[which(genome_stats$mOTU == "ELBEEST_200" & genome_stats$classification == "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__"), "classification"] <- "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sulfonica"
genome_stats[which(genome_stats$mOTU == "ELBEEST_221" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Porticoccaceae;g__HTCC2207;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Porticoccaceae;g__HTCC2207;s__HTCC2207 sp002685195"
##Again here we have a genus conflict with 12/14 mOTUs assigned to a single genus which we will take
genome_stats[which(genome_stats$mOTU == "ELBEEST_272" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__SG8-39;g__CAILKO01;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__SG8-39;g__SG8-39;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_310" & genome_stats$classification == "d__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__Acidimicrobiales;f__Ilumatobacteraceae;g__UBA3006;s__"), "classification"] <- "d__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__Acidimicrobiales;f__Ilumatobacteraceae;g__UBA3006;s__UBA3006 sp903850275"
genome_stats[which(genome_stats$mOTU == "ELBEEST_362" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_A;g__JAHDXY01;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_A;g__JAHDXY01;s__JAHDXY01 sp023257915"
genome_stats[which(genome_stats$mOTU == "ELBEEST_370" & genome_stats$classification == "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Pedosphaerales;f__AAA164-E04;g__;s__"), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Pedosphaerales;f__AAA164-E04;g__AAA164-E04;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_376" & genome_stats$classification == "d__Bacteria;p__Desulfobacterota_B;c__Binatia;o__HRBIN30;f__;g__;s__"), "classification"] <- "d__Bacteria;p__Desulfobacterota_B;c__Binatia;o__HRBIN30;f__JAGDMS01;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_403" & genome_stats$classification == "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__UBA3478;s__"), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__UBA3478;s__UBA3478 sp016780435"
genome_stats[which(genome_stats$mOTU == "ELBEEST_412" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_A;g__UBA2463;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_A;g__UBA2463;s__UBA2463 sp945901825"
genome_stats[which(genome_stats$mOTU == "ELBEEST_427" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Micavibrionales;f__UBA2020;g__;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Micavibrionales;f__UBA2020;g__UBA2020;s__"
##here the most are unclassifed at the species level with 3 and 1 assigned to two genera. Lets be rational and take up to family 
genome_stats[which(genome_stats$mOTU == "ELBEEST_439" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__PHCI01;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_439" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__JAKFVA01;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__;s__"
##Again two genera 5/1 split
genome_stats[which(genome_stats$mOTU == "ELBEEST_455" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__Rhodoferax_A;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_B;g__CAIKVZ01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_458" & genome_stats$classification == "d__Bacteria;p__Patescibacteria;c__Paceibacteria;o__UBA9983_A;f__JAACPR01;g__;s__"), "classification"] <- "d__Bacteria;p__Patescibacteria;c__Paceibacteria;o__UBA9983_A;f__JAACPR01;g__JAGOTN01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_492" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudohongiellaceae;g__;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudohongiellaceae;g__CAILUG01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_504" & genome_stats$classification == "d__Bacteria;p__Myxococcota;c__UBA727;o__UBA727;f__JABDBI01;g__;s__"), "classification"] <- "d__Bacteria;p__Myxococcota;c__UBA727;o__UBA727;f__JABDBI01;g__JAHFEG01;s__"
##Here we have a family split 10/4. im tempted to take to order level but for now will select dominant family
genome_stats[which(genome_stats$mOTU == "ELBEEST_514" & genome_stats$classification == "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__DEV007;g__;s__"), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__SKLO01;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_517" & genome_stats$classification == "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__;g__;s__"), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__SKLO01;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_519" & genome_stats$classification == "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__;s__"), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__Lutibacter;s__"
##genus 6/1 split
genome_stats[which(genome_stats$mOTU == "ELBEEST_542" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Methylophilaceae;g__BACL14;s__BACL14 sp905181685"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Methylophilaceae;g__BACL14;s__BACL14 sp019823025"
genome_stats[which(genome_stats$mOTU == "ELBEEST_58" & genome_stats$classification == "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__"), "classification"] <- "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp014190015"
genome_stats[which(genome_stats$mOTU == "ELBEEST_585" & genome_stats$classification == "d__Bacteria;p__Acidobacteriota;c__Vicinamibacteria;o__Vicinamibacterales;f__UBA2999;g__;s__"), "classification"] <- "d__Bacteria;p__Acidobacteriota;c__Vicinamibacteria;o__Vicinamibacterales;f__UBA2999;g__CADEFD01;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_59" & genome_stats$classification == "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__"), "classification"] <- "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp014190015"
##Family 3/1 split 
genome_stats[which(genome_stats$mOTU == "ELBEEST_601" & genome_stats$classification == "d__Bacteria;p__Patescibacteria;c__Paceibacteria;o__UBA9983_A;f__W02-35-19;g__;s__"), "classification"] <- "d__Bacteria;p__Patescibacteria;c__Paceibacteria;o__UBA9983_A;f__JAGLPS01;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_623" & genome_stats$classification == "d__Bacteria;p__Bdellovibrionota;c__UBA1018;o__UBA1018;f__UBA1018;g__;s__"), "classification"] <- "d__Bacteria;p__Bdellovibrionota;c__UBA1018;o__UBA1018;f__UBA1018;g__CAINWE01;s__"
##Genus 2/1 split, definately going up to family level
genome_stats[which(genome_stats$mOTU == "ELBEEST_700" & genome_stats$classification == "d__Bacteria;p__Cyanobacteriota;c__Vampirovibrionia;o__LMEP-6097;f__LMEP-6097;g__JAJTHA01;s__"), "classification"] <- "d__Bacteria;p__Cyanobacteriota;c__Vampirovibrionia;o__LMEP-6097;f__LMEP-6097;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_700" & genome_stats$classification == "d__Bacteria;p__Cyanobacteriota;c__Vampirovibrionia;o__LMEP-6097;f__LMEP-6097;g__CAIYXB01;s__"), "classification"] <- "d__Bacteria;p__Cyanobacteriota;c__Vampirovibrionia;o__LMEP-6097;f__LMEP-6097;g__;s__"
genome_stats[which(genome_stats$mOTU == "ELBEEST_779" & genome_stats$classification == "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rickettsiales;f__UBA1997;g__;s__"), "classification"] <- "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rickettsiales;f__UBA1997;g__UBA2645;s__"


##So that is all of the BICEST mOTUs, we will deal with the mOTUs3.1 seperately because they are built from MAGs that might have a better taxonomic assignment than ours. 
##

genome_stats %>% 
  filter(mOTU != "NotEnoughMGs") %>%
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  select(mOTU) %>% 
  pull()

##Ok so now we only have a few conflicts relating to the ext_mOTUs, first lets extract the classificaiton for the BICEST mOTUs

elbe_taxa <- genome_stats %>% 
  filter(grepl("ELBEEST", mOTU)) %>% 
  select(mOTU, classification) %>% 
  group_by(mOTU) %>%
  unique()
##Quick sanity check

length(unique(elbe_taxa$mOTU)) ##1046
nrow(elbe_taxa) ##1046

##Ok so now we have a single taxonomy for each BICEST mOTU!!

##Lets deal with the reference mOTUs. 

motus_taxa <- read.table(paste0(path, "/mOTUs3.1.0.genome_metadata.tsv"), sep="\t", header=TRUE, quote="") %>% 
  select(GENOME..MOTU, GTDB.R207) %>% 
  rename(mOTU = GENOME..MOTU, classification = GTDB.R207)
nrow(motus_taxa) ## so we have taxonomy for 700K MAGs/genomes!!
length(unique(motus_taxa$mOTU)) ## This is around 34195 mOTUs

##Lets see how many are in our profile

motus %>% 
  rownames_to_column("mOTU") %>% 
  select(mOTU) %>%
  filter(!grepl("ELBE", mOTU)) %>% 
  nrow() ## so we have 34341 mOTUs 

##This is a very similar number to the total number of mOTUs so maybe there are quite some zero values

refmOTUs_in_profile <- motus %>% 
  rownames_to_column("mOTU") %>% 
  pivot_longer(!mOTU, names_to = "samples", values_to = "counts") %>% 
  group_by(mOTU) %>%
  mutate(total = sum(counts)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>% 
  select(mOTU)  %>% 
  unique() %>%
  filter(!grepl("ELBE", mOTU)) %>% 
  pull() ##%>%
  #length() ## Ok so now we have only 2072 ref mOTUs, also for sanity we can remove ! from filter and see that all 1046 BICEST mOTUs are also detected in our profiles.  

##Ok lets check some taxonomy!!  

##First lets remove all of the information for mOTUs not in our profiles      
  
motus_taxa <- motus_taxa %>% 
    rownames_to_column("user_genome") %>% 
    filter(mOTU %in% refmOTUs_in_profile) %>% 
    select(mOTU, classification, user_genome)

nrow(motus_taxa) ##Ok so we are down from 700K genomes to 126k. 
length(unique(motus_taxa$mOTU)) ##Now we lost 11 mOTUs. This is meta_mOTUs which dont have taxonomy associated with them, they are built on MG profiles not MAGs (check paper for details) - we will fix later

##Lets add also our genome taxonomy information incase this helps 

motus_taxa <- genome_stats %>% 
  filter(mOTU %in% motu_to_check) %>% 
  select(mOTU, classification, user_genome) %>% 
  rbind(motus_taxa) 

##So we added around 180 MAGs to the 126K (obviously not much) mOTUs3.1.0 metadata file from our own dataset

##Now we make again a list of taxa with conflict
motu_to_check <- motus_taxa %>%
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  select(mOTU) %>% 
  pull()

##Great so now we have a list of around 200 mOTUs to check. Lets do this!!

motus_taxa %>% 
  filter(mOTU %in% motu_to_check) %>% 
  select(mOTU, classification, user_genome) %>% 
  group_by(mOTU, classification) %>%
  summarise(count = n_distinct(user_genome)) %>%
  arrange(mOTU) %>% 
  print(n=800)

##Now I just realised that this could be a lot more efficient than what we did above.

##code hidden here
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_14660"), "classification"] <- "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosarchaeum;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_15383"), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Atopobiaceae;g__NM07-P-09;s__NM07-P-09 sp004793665"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_15977"), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__S36-B12;g__S36-B12;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_16111"), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Rhodoferax;s__Rhodoferax sp014190675"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_16115"), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Crocinitomicaceae;g__UBA952;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_16131"), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Polynucleobacter;s__Polynucleobacter sp009927445"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_17177"), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae;g__Pelagibacter;s__Pelagibacter sp902565935"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_17878"), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19170" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Crocinitomicaceae;g__M0103;s__M0103 sp903930085"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19171" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp009702965"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19172" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Nanopelagicus;s__Nanopelagicus sp001437855"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19173" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Nanopelagicus;s__Nanopelagicus sp001437855"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19207" ), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__V1-33;g__CAJBME01;s__CAJBME01 sp903954065"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19221" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Saprospiraceae;g__BJGN01;s__BJGN01 sp014190515"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19232" ), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Opitutales;f__UBA953;g__UBA953;s__UBA953 sp004293385"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19239" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Limnohabitans_A;s__Limnohabitans_A sp001517545"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19252" ), "classification"] <- "dd__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp014190015"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19259" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Nanopelagicus;s__Nanopelagicus sp003569185"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19260" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Methylophilaceae;g__Methylotenera;s__Methylotenera sp903951385"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19262" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila vernalis"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19273" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__Sediminibacterium;s__Sediminibacterium sp017987615"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19287" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__Flavobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19290" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Crocinitomicaceae;g__Fluviicola;s__Fluviicola sp017983675"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19295" ), "classification"] <- "d__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp014190015"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19298" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Hyphomonadaceae;g__Hyphomonas;s__Hyphomonas sp016124495"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19383" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Schleiferiaceae;g__TMED14;s__TMED14 sp005786915"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19415" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Cytophagales;f__Cyclobacteriaceae;g__ELB16-189;s__ELB16-189 sp016787665"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19437" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Thermoleophilia;o__Gaiellales;f__F1-60-MAGs149;g__F1-60-MAGs149;s__F1-60-MAGs149 sp903960645"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19455" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__Planktophila sp903839585"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19461" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Crocinitomicaceae;g__Fluviicola;s__Fluviicola sp017983675"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19467" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingorhabdus_B;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_19767" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__MED-G116;s__MED-G116 sp004212735"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_21813" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__S36-B12;g__S36-B12;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_22053" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Nitrincolaceae;g__ASP10-02a;s__ASP10-02a sp002335115"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_21767" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Puniceispirillales;f__Puniceispirillaceae;g__MED-G116;s__MED-G116 sp004212735"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_28513" ), "classification"] <- "d__Bacteria;p__Chloroflexota;c__Limnocylindria;o__Limnocylindrales;f__Limnocylindraceae;g__Limnocylindrus;s__Limnocylindrus sp903937225"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_28707" ), "classification"] <- "d__Bacteria;p__Planctomycetota;c__Planctomycetia;o__Gemmatales;f__Gemmataceae;g__UBA969;s__"
motus_taxa[which(motus_taxa$mOTU == "ext_mOTU_v31_28709" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__OLB11;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12332" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri_A"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12366" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter faecis"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12572" ), "classification"] <- "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__Nitrosopumilus sp002690535"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12608" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__HIMB59;f__HIMB59;g__HIMB59;s__HIMB59 sp902529795"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12714" ), "classification"] <- "d__;p__;c__;o__;f__;g__;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12781" ), "classification"] <- "d__Archaea;p__Thermoplasmatota;c__Poseidoniia;o__Poseidoniales;f__Poseidoniaceae;g__MGIIa-L2;s__MGIIa-L2 sp002171315"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12833" ), "classification"] <- "d__Bacteria;p__;c__;o__;f__;g__;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_12945" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Rhodothermia;o__Balneolales;f__Balneolaceae;g__UBA1275;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13018" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Puniceispirillales;f__AAA536-G10;g__AAA536-G10;s__AAA536-G10 sp016777705"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13056" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__Luminiphilus;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13185" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__MED-G52;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13209" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__Luminiphilus;s__Luminiphilus sp905182485"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13257" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__GCA-2697345;s__GCA-2697345 sp002697345"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13622" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__GCA-002733185;s__GCA-002733185 sp002713705"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_13782" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__UBA10066;g__MED-G20;s__MED-G20 sp016780365"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_14148" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__;f__;g__;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_14237" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Halieaceae;g__Luminiphilus;s__"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_14284" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__SAR86;f__SAR86;g__GCA-2707915;s__GCA-2707915 sp016777005"
motus_taxa[which(motus_taxa$mOTU == "meta_mOTU_v31_14430" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__;g__;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00049" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00050" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00051" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00052" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00056" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00085" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00086" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00095" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00122" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00129" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00131" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00133" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00140" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00150" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00164" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00167" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00168" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00169" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00170" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00188" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00201" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00206" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00215" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00227" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00259" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter baumannii"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00267" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00283" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00296" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00312" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00328" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00344" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus hominis"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00424" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Janthinobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00444" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Curvibacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00453" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Acidovorax_A;s_"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00459" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Comamonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00474" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Variovorax;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00483" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Variovorax;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00485" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Variovorax;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00487" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Variovorax;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00520" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00521" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Duganella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00547" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Achromobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00576" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Photobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00728" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Micromonosporaceae;g__Micromonospora;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00738" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Micromonosporaceae;g__Micromonospora;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00772" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__Mesorhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00800" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium;s__Cutibacterium acnes"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00802" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Corynebacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00853" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00855" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00870" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00871" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00873" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00881" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00882" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Comamonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00964" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Micrococcaceae;g__Micrococcus;s__Micrococcus luteus"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00996" ), "classification"] <- "d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter_B;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_00999" ), "classification"] <- "d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01010" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;f__Rickettsiaceae;g__Rickettsia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01011" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;f__Rickettsiaceae;g__Rickettsia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01035" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Oleiphilaceae;g__Marinobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01043" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__Rhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01071" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Brevundimonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01089" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Yersinia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01096" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Pseudoduganella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01141" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Beijerinckiaceae;g__Methylobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01170" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Macellibacteroides;s__Macellibacteroides fermentans"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01185" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01191" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01211" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01213" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus;s_"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01234" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Paraburkholderia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01241" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Cupriavidus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01242" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Cupriavidus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01243" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Cupriavidus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01274" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01337" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01338" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01361" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Chromatiaceae;g__Marichromatium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01404" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01513" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Cereibacter_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01518" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Azospirillales;f__Azospirillaceae;g__Azospirillum;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01547" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Sphingobacteriales;f__Sphingobacteriaceae;g__Pedobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01556" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01692" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01693" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01786" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01789" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingomonas;s__Sphingomonas aquatilis"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01861" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pseudoalteromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01862" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pseudoalteromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01866" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pseudoalteromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01874" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Microbacteriaceae;g__Clavibacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01920" ), "classification"] <- "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__PCC-6307;f__Cyanobiaceae;g__Synechococcus_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01940" ), "classification"] <- "d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01941" ), "classification"] <- "d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_01953" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02019" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02037" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Weeksellaceae;g__Elizabethkingia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02129" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__Aurantimonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02154" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02348" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Stappiaceae;g__Pannonibacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02367" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02398" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02485" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02580" ), "classification"] <- "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcystaceae;g__Microcystis;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02588" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_K;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02642" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Polaromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02650" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02696" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Brevundimonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02700" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02783" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_D;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02795" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02801" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02804" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Rhodanobacteraceae;g__Rhodanobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_02810" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03205" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03341" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Dysosmobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03342" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03567" ), "classification"] <- "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcoleaceae;g__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03572" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Tissierellales;f__Peptoniphilaceae;g__Finegoldia;s_"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03577" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03591" ), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03640" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03657" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03701" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03761" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Aeromonadaceae;g__Aeromonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03850" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_03854" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Sphaerotilus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04025" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Oleiphilaceae;g__Marinobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04078" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04204" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Ascidiaceihabitans;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04231" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Arenimonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04265" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae;g__Pelagibacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04339" ), "classification"] <- "d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04469" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04470" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04516" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04596" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Rhizobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04616" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Microbacteriaceae;g__Plantibacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04737" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04792" ), "classification"] <- "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Sphingobacteriales;f__Sphingobacteriaceae;g__Sphingobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04802" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Corynebacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_04827" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Halomonadaceae;g__Halomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05125" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Pseudonocardiaceae;g__Pseudonocardia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05354" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05380" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05437" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Phenylobacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05522" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__;g__;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05841" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Marinomonadaceae;g__Marinomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_05975" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06002" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Moraxella_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06109" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06146" ), "classification"] <- "d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Aliarcobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06160" ), "classification"] <- "d__Bacteria;p__Myxococcota;c__Myxococcia;o__Myxococcales;f__Anaeromyxobacteraceae;g__Anaeromyxobacter;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06277" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06494" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae;g__IMCC9063;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_06514" ), "classification"] <- "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Pedosphaerales;f__AAA164-E04;g__AAA164-E04;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_07361" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__SAR86;f__D2472;g__D2472;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_09389" ), "classification"] <- "d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Romboutsia;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_10657" ), "classification"] <- "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae;g__Pelagibacter_A;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_12136" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Nanopelagicus;s__"
motus_taxa[which(motus_taxa$mOTU == "ref_mOTU_v31_12137" ), "classification"] <- "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Nanopelagicales;f__Nanopelagicaceae;g__Planktophila;s__"


motus_taxa %>%
  select(mOTU, classification) %>% 
  group_by(mOTU) %>% 
  summarise(count = n_distinct(classification)) %>%
  filter(count > 1) %>% 
  select(mOTU) %>% 
  pull()

##Hurrah so now we have a single taxonomy for each of the mOTUs3.1.0 db mOTUs. 

motu_final_taxa <- motus_taxa %>% 
  select(mOTU, classification) %>% 
  group_by(mOTU) %>%
  unique()

##Quick sanity check

length(unique(motu_final_taxa$mOTU)) ##2063
nrow(motu_final_taxa) ##2063

##If this doesnt look right reload mOTUs db and try again, i was having some issues.

elbe_motu_taxa <- rbind(elbe_taxa, motu_final_taxa) %>% 
  separate(classification, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";")

##Lets also do what we didnt before and create a simplified mOTU table

motus_final <- motus %>% 
  rownames_to_column("mOTU") %>% 
  pivot_longer(!mOTU, names_to = "samples", values_to = "counts") %>% 
  group_by(mOTU) %>%
  mutate(total = sum(counts)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>% 
  pivot_wider(names_from = "samples", values_from = "counts")

nrow(motus_final) #3118 which is longer than our 3109 taxa, this includes some meta_motus that dont have taxonomy, as well as unclassified seqs, lets fix this.

elbe_motu_taxa <- motus_final %>%
  select(mOTU) %>% 
  filter(!mOTU %in% elbe_motu_taxa$mOTU) %>% 
  mutate(classification = "d__;p__;c__;o__;f__;g__;s__") %>%
  separate(classification, c("domain", "phylum", "class", "order", "family",  "genus", "species"), ";") %>%
  rbind(.,elbe_motu_taxa) 
  
##Great so now we have taxonomy for the 3118 mOTUs (inc unassigned) that occur in our profile


#All we are missing is the metadata, i steal code chunk from the functional profile##
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog/") #edited by JG
path= "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog/"
metadata <- read.csv("SAMEAID_SampleID.csv", header=TRUE, sep=";", dec=",") %>% 
  mutate(sampleid=paste0(BioSample,"_METAT")) %>%
  mutate(data_type="METAT")

metadata2 <- read.csv("SAMEAID_SampleID.csv", header=TRUE, sep=";", dec=",") %>% 
  mutate(sampleid=paste0(BioSample,"_METAG"))%>%
  mutate(data_type="METAG")

metadata <- rbind(metadata, metadata2) 

##There are some missing values for replicates and I want to simplify headers using good practice (lowercases, no special char, english)
##JG  change order more selected some hashed
metadata <- metadata %>% 
  select(sampleid, station=Station, station_km = Stromkilometer, sample_type = Sample_type, date=Sample_date, o2_sat=Sat_O2_TBDHereon, wtemp=Temperature_TBDHereon, 
         salinity=Salinity_TBDHereon, turbidty=Turbidity_TBDHereon, ph=pH_TBDHereon, o2_conc=O2_TBDHereon, spm_mg_l=SPM_mgperL, doc_mg_l=DOC_mg.L, tn_mg_l=TN_mg.L,
         dic_mg_l=DIC_mg.L, si_mg_l=Silicate_mg.L, nh4_mg_l=Ammonium_mg.L, no2_mg_l=Nitrite_mg.L, no3_mg_l=Nitrate_mg.L, BioSample, ERANumber, data_type)#, din_um=Total_DIN_ÂµM, srp_ml_l=SRP_mgperL,
         #tdp_mg_l=TotalDissolvedPhosphate_mg.L, nh4_um = Ammonium_ÂµM, no2_um = Nitrite_ÂµM, no3_um = Nitrate_ÂµM, srp_um = SRP_ÂµM, tdp_um = Phosphate_ÂµM,
         #si_um = Silicate_ÂµM, doc_um = DOC_uM.L, dic_um = DIC_uM.L, respiration_o2_ug_l_h=RespirationRate_O2ug.L.h, poc_mg_l = POC_mgperL, ptc_mg_l = PTC_mgperL, 
         #ptn_mg_l = PTN_mgperL, pth_mg_l = PTH_mgperL

##Sometimes the replicates miss values, to correct we will group by station, date, sample_type(fraction), then calculate the mean value ignoring na. 
##If there isnt a measurement it returns a NaN value

metadata <- metadata %>% 
  group_by(station, date, sample_type) %>%
  mutate(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>% 
  ungroup()


##Lets create a dataframe to hold all the data (yes this can be done with phyloseq but phyloseq has its own problems and we usually end up exporting tables back out for vegan anyway so we can work from here )

motus_final_taxa_metadata <- motus_final %>%
  pivot_longer(!mOTU, names_to = "sampleid", values_to = "counts") %>% 
  left_join(., elbe_motu_taxa) %>% 
  left_join(., metadata)

##Lets save what we have produced and cleanup some of the rubbish
##Its a good idea to hash this out so we dont accidentally overwrite anything
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/MOTUs")
saveRDS(motus_final_taxa_metadata, file="motus_final_taxa_metadata.RDS")
saveRDS(motu_stats, file="motu_stats.RDS")
saveRDS(metadata, file="metadata.RDS")

rm(list=ls())

motus_final_taxa_metadata <- readRDS(file="motus_final_taxa_metadata.RDS")
metadata <- readRDS(file="metadata.RDS")
##Lets check out experimental scheme

motus_final_taxa_metadata %>% 
  select(sampleid, data_type, date, station, sample_type) %>% 
  group_by(data_type, sample_type, date) %>% 
  summarise(count=n_distinct(sampleid)) %>% 
  print(n=30) 

##So when we consider fraction and date we see that Nov 21 has few samples, we should consider removing

motus_final_taxa_metadata %>% 
  select(sampleid, data_type, date, station, sample_type) %>% 
  group_by(data_type, sample_type, station) %>% 
  summarise(count=n_distinct(sampleid)) %>% 
  print(n=30) 

##By station we see that Kollmar, SeemanshÃ¶ft and BunthausSpitze also were sampled infrequently, lets also remove these. 

motus_final_taxa_metadata <- motus_final_taxa_metadata %>% 
  filter(!grepl("Bunthaus|Kollmar|Seemans", station), !grepl("Nov 21", date))


####Lets see how our counts look across all the samples, im going to create sub_files for plotting rather than piping directly to ggplot because it takes time
##when we want to adjust visuals

sample_specific_counts <- motus_final_taxa_metadata %>% 
  group_by(sampleid) %>%
  mutate(sample_sums=sum(counts)) %>% 
  select(sampleid, data_type, sample_type, sample_sums)

sample_specific_counts %>%
  ggplot() +
  geom_point(aes(x=sampleid, y=sample_sums)) 

##Ok there are clearly some samples with far fewer mOTU counts. Lets see how this looks when we compare metaT and metaG, also for the different size fractions

sample_specific_counts %>%
  ggplot() +
  geom_point(aes(x=sampleid, y=sample_sums)) +
  facet_wrap(~data_type*sample_type, scales="free_y", nrow=2)

###So take homes. We shouldnt compare metaG an metaT because the counts are an order of magnitude. You can look at metaT later. For now i will ignore. 
###Should we compare FL and PA (maybe but also maybe not), 
###There is an argument for comparing FL and PA but generally you want to explore what are the drivers of composition independent of fraction anyway. 
##Lets split the files, i will keep heavy and light together, splitting more reduces the number of pairwise comparisons which also helps.

colnames(motus_final_taxa_metadata)

motus_metadata_metag_fl <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG" & sample_type == "Free_living")

motus_metadata_metag_pa <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG" & sample_type != "Free_living")

##Lets start with free-living by first checking the sample_sums again

motus_metadata_metag_fl %>% 
  group_by(sampleid) %>%
  mutate(sample_sums=sum(counts)) %>%
  select(sampleid, station, date, sample_sums) %>%
  unique() %>% 
  arrange(sample_sums) %>% 
  print(n=61)

##So we can see here that counts from Feb 22 are also quite low from all free-living samples. 

motus_metadata_metag_pa %>% 
  group_by(sampleid) %>%
  mutate(sample_sums=sum(counts)) %>%
  select(sampleid, station, date, sample_sums) %>%
  unique() %>% 
  arrange(sample_sums) %>% 
  print(n=61)

##So we can see here that counts from Feb 22 are also quite low from all particle-associated samples. These are in theory the only winter samples. 
##These samples generally had lower read counts. Anyway lets compare this also with the functional profiles to make an assessment. 
##We dont want to throw away a whole month of samples if the data are biologically interesting but also we shouldnt force a perspective based on sampling depth

##For now I will remove them

#motus_final_taxa_metadata <- motus_final_taxa_metadata %>%  ### Do no remove for mentel test
#  filter(!grepl("Feb 22", date)) #sample_date to date

motus_metadata_metag_fl <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG" & sample_type == "Free_living")

motus_metadata_metag_pa <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG" & sample_type != "Free_living")

##Lets check again

motus_metadata_metag_fl %>% 
  group_by(sampleid) %>%
  mutate(sample_sums=sum(counts)) %>%
  select(sampleid, station, date, sample_sums) %>%
  unique() %>% 
  arrange(sample_sums) %>% 
  print(n=61)

##Theres two samples with low counts (lower than 14k), lets remove these also

motus_metadata_metag_fl <- motus_metadata_metag_fl %>% 
  filter(!grepl("SAMEA112714820_METAG|SAMEA110290254_METAG", sampleid))

##and the particle associated

motus_metadata_metag_pa %>% 
  group_by(sampleid) %>%
  mutate(sample_sums=sum(counts)) %>%
  select(sampleid, station, date, sample_sums) %>%
  unique() %>% 
  arrange(sample_sums) %>% 
  print(n=61)

#Here the sums are genearlly lower with a plateau around 8k reads. For now we will leave everything and come back if these cause issues. 

###Ok lets calculate some distances, 

##First we want to covert our data back to an OTU table 

motus_fl_df <- motus_metadata_metag_fl %>%
  select(mOTU, sampleid, counts) %>% #problem with pivot more then one value per group of count and sample try to leave distinct lines - works
  distinct() %>% 
  pivot_wider(names_from="mOTU", values_from="counts", values_fill=0) %>% 
  column_to_rownames("sampleid") %>%
  as.data.frame()

##Lets check this works first with vegdist which doesnt resample
library(vegan)
library(glue)
motus_fl_norare_dist_matrix <- vegdist(motus_fl_df, method="bray")

##great that worked, we can ordinate the data 

pcoa <- cmdscale(motus_fl_norare_dist_matrix, eig=TRUE, add=TRUE)
positions_fl <- pcoa$points
colnames(positions_fl) <- c("pcoa1", "pcoa2")
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
            glue("PCo Axis 2 ({pretty_pe[2]}%)"))

positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color = date, shape=station)) +
  geom_point(size=4) +
  labs(x=labels[1], y=labels[2]) + 
  theme_classic()+ theme(legend.position = "right") +
  guides(colour = guide_legend(title.position="top", title.hjust = 0, title = "Date "),
         shape = guide_legend(title.position="top", title.hjust = 0, title = "Station"))

##Noice, so here we see a strong station effect (PC1, marine left and brackish right) and a strong season effect PC2 (May bottom, Jun-Nov top)

##Lets run this with avgdist which does rarefaction and subsampling https://www.biorxiv.org/content/biorxiv/early/2023/06/26/2023.06.23.546312.full.pdf

motus_fl_rare_dist_matrix <- avgdist(motus_fl_df, dmethod="bray", sample=2400)
pcoa <- cmdscale(motus_fl_rare_dist_matrix, eig=TRUE, add=TRUE)
positions_fl <- pcoa$points
colnames(positions_fl) <- c("pcoa1", "pcoa2")
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
            glue("PCo Axis 2 ({pretty_pe[2]}%)"))

positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color = date, shape=station)) +
  geom_point(size=4) +
  labs(x=labels[1], y=labels[2]) + 
  theme_classic()+ theme(legend.position = "right") +
  guides(colour = guide_legend(title.position="top", title.hjust = 0, title = "Date "),
         shape = guide_legend(title.position="top", title.hjust = 0, title = "Station"))

##So when we apply rarefaction to the data the overall trend remains the same but we are explaining much more variation in the data.
##In particular we see that the replicates are much closer to one another suggesting nearly all variation is due to sampling depth. 

##Go back and recreate the motus_fl_df, this time do not remove the samples from Feb 22. We know they have less sequences so lets check what impact rarefaction has. 
##Make sure to change sample depth from 14000, to 2400

##So when we run vegdist with the Feb22 samples we see that they all cluster separately from the others, independent of station. This was true for gene catalog

##When we run avgdist with 2400 reads we maintain the same structure as before without Feb22 so thats not a problem
##The FEb22 samples no longer stand out, they actually group woth the Mai samples to some extent. 
##We can also see now the gradient in salinity although it is less pronounced. This resembles Jun22 and Mai 21 where BB is less resolved.

##We can test our environmental factors against this

##First lets check what we actually have complete data for, we dont want to remove too many samples. 

nrow(positions_fl)


env_params <- positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  select(-pcoa1, -pcoa2) %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>%
#   filter(!is.na(wtemp)) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(~ ! any(is.character(.))) %>% 
  select(-station_km) %>% 
  colnames()## Ok thats not the worst list, lets export this so we can read it into our mantel function

##Lets also do a quick correlation matrix with these parameters 
library(rstatix)
env_params_corr <- positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>% 
#  filter(!is.na(wtemp)) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(~ ! any(is.character(.))) %>% 
  cor_mat(method = "spearman") %>%
  as.matrix() %>%
  as_tibble() %>%
  pivot_longer(-rowname)

env_params_corr_p <- positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>% 
#  filter(!is.na(wtemp)) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(~ ! any(is.character(.))) %>% 
  cor_pmat(method = "spearman") %>%
  as.matrix() %>%
  as_tibble() %>%
  pivot_longer(-rowname) %>%
  select(rowname, name, p_value = value) %>% 
  cbind(., corr=env_params_corr$value) %>% 
  filter(rowname < name) %>%
  mutate(p_value = as.numeric(p_value), corr=as.numeric(corr))

env_params_corr_p %>% 
  ggplot(aes(x=rowname,y=name)) +
  geom_point(aes(size=corr, colour=corr, ), shape=15) +
  geom_point(aes(), shape=0, size = 8.5) + 
  scale_colour_gradient2(low = "red", mid = "white", high = "darkgreen", midpoint=0) + 
  guides(colour = guide_colorbar(title = "Mantel's r", frame.colour = "black", barwidth = 1.5)) + #title.theme = element_text(angle=90, vjust=1))) +
  scale_size_continuous(guide="none") +
  scale_x_discrete()+ theme_bw() +
  scale_y_discrete(position = "right") +
  coord_fixed(ratio = 1) +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), line =element_blank(),
        strip.background = element_blank()) +
  theme(rect = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1.1, vjust=1.1), 
        legend.position = "right") +
  geom_text(vjust=0.76,hjust=0.49, aes(label = paste0(c("","*")[(abs(p_value) <= .05)+1], c("","*")[(abs(p_value) <= .001)+1], c("","*")[(abs(p_value) <= .001)+1])))

##Plot salinity v pc1

positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>% 
  filter(!is.na(wtemp)) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(~ ! any(is.character(.))) %>% 
  ggplot(aes(x=pcoa1, y=o2_sat)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x)

##Obviously some critical data, salinity temp etc are missing. 

##We will create a function to test our env_params. Note that it subsamples the dist_matrix so we can test for example wtemp where data is missing

mantel_function = function(parameter){
motus_fl_rare_dist_matrix_mantel <- motus_fl_rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>% 
  left_join(., metadata, by = c("sample" = "sampleid"))%>% 
  left_join(., metadata, by = c("name" = "sampleid")) %>% 
  filter(paste0(parameter,".x") != "NA", paste0(parameter,".y") != "NA") %>%  ##  This allows us to remove missing data, but ideally you shouldnt have any
  select(sample, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-sample) %>% 
  as.dist()  

motus_fl_env_dist_matrix_mantel <- motus_fl_rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>% 
  left_join(., metadata, by = c("sample" = "sampleid"))%>% 
  left_join(., metadata, by = c("name" = "sampleid")) %>% 
  filter(paste0(parameter,".x") != "NA", paste0(parameter,".y") != "NA") %>%  ##  This allows us to remove missing data, but ideally you shouldnt have any
  select(sample, paste0(parameter,".x")) %>%
  unique() %>%
  column_to_rownames("sample") %>%
  dist(method="euclidean")

mantel_out = mantel(motus_fl_rare_dist_matrix_mantel, motus_fl_env_dist_matrix_mantel, method = "spearman", permutations = 9999, na.rm = TRUE)
cbind(env=parameter, mantel=mantel_out$statistic, sig=mantel_out$signif)
}

##Now we can write a little loop to calculate mantel test statistics for each of our parameters

out = data.frame(env=NA, mantel=NA, sig=NA)
for(env in env_params){
out <- rbind(out, mantel_function(parameter = env))
}

out %>% 
  mutate(sig=as.numeric(sig)) %>%
  filter(sig < 0.05) %>%
  arrange(desc(mantel))

##So we can see now that DOC is the most well correlated, followed by TN, DIC and nitrate. 

##Lets quickly run this again with temperature 

positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  select(-pcoa1, -pcoa2) %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>% 
  select(station, date, wtemp) %>% 
  filter(is.na(wtemp)) ##%>%
#  nrow() ## so we are missing six samples for wtemp, removing Nov22 #JG commented out
mantel_function(parameter = "wtemp")

## We see that wtemp is similar correlated as is DOC

##Lets quickly run this again with temperature 

positions_fl %>% 
  as_tibble(rownames = "sampleid") %>%
  select(-pcoa1, -pcoa2) %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>% 
  select(station, date, salinity) %>% 
  filter(is.na(salinity)) ##%>%
#  nrow() ## so we are missing six samples for wtemp, removing Nov22 
mantel_function(parameter = "salinity")

##Salinity however is massive, 0.7 compared to 0.3 for DOC and wtemp. Again we are missing six samples here so I suggest you go and get these sorted. 


###Revisiting the size fractions, 

##Considering we observed how sub-sampling at 2400 was sufficient to reproduce the same patterns we saw when sampling at 14000 we can prob compare all sample_types

motus_df <- motus_final_taxa_metadata %>%
  filter(data_type == "METAG") %>%
  distinct() %>% #once again distinct to make it works
  select(mOTU, sampleid, counts) %>%
  pivot_wider(names_from="mOTU", values_from="counts", values_fill=0) %>% 
  column_to_rownames("sampleid") %>%
  as.data.frame()

motus_rare_dist_matrix <- avgdist(motus_df, dmethod="bray", sample=2400)
pcoa <- cmdscale(motus_rare_dist_matrix, eig=TRUE, add=TRUE)
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
            glue("PCo Axis 2 ({pretty_pe[2]}%)"))

positions %>% 
  as_tibble(rownames = "sampleid") %>%
  left_join(metadata, by = c("sampleid" = "sampleid")) %>%
  mutate(station_type = paste(station, sample_type)) %>%
  ggplot(aes(x=pcoa1, y=pcoa2, colour = date, shape=station, alpha=sample_type)) +
  geom_point(size=4) +
  labs(x=labels[1], y=labels[2]) + 
  theme_classic()+ theme(legend.position = "right") +
  guides(colour = guide_legend(title.position="top", title.hjust = 0, title = "Date "),
         shape = guide_legend(title.position="top", title.hjust = 0, title = "Station"))

##Ok so there actually does not seem to be a huge variation between the different fractions

library(ggtern)
library(microbiomeutilities)

#so to do a ternary plot we actually need to create a physeq object really quick 
motus_table <- motus_final_taxa_metadata %>%
  filter(data_type == "METAG") %>%
  select(mOTU, sampleid, counts) %>%
  group_by(sampleid) %>% 
  mutate(total=sum(counts)) %>% 
  filter(total > 2400) %>%  ##Here we just remove those with low counts
  select(-total) %>% 
  ungroup() %>% 
  distinct() %>%
  pivot_wider(names_from = sampleid, values_from = counts, values_fill = 0) %>% 
  column_to_rownames("mOTU") %>% 
  phyloseq::otu_table(taxa_are_rows = TRUE)
motus_taxa <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG") %>%
  select(mOTU, phylum, class, order, family, genus, species) %>% 
  unique() %>% 
  column_to_rownames("mOTU") %>% 
  as.matrix() %>%
  phyloseq::tax_table()
motus_sampledata <- motus_final_taxa_metadata %>% 
  filter(data_type == "METAG") %>%
  select(-phylum, -class, -order, -family, -genus, -species, -counts, -mOTU, -domain) %>% 
  unique() %>%  
  column_to_rownames("sampleid") %>% 
  sample_data()
bicest_motu_ps <- phyloseq(motus_table, motus_taxa, motus_sampledata)
bicest_motu_ps

##Great now lets make a ternary plot showing the differences between the fractions, 
tern_table <- prep_ternary(
  bicest_motu_ps,
  abund.thres = 1e-04, ##We can adjust this if we feel there are too few/many taxa being plotted
  prev.thres = 0.1, ##We can adjust this if we feel there are too few/many taxa being plotted
  group = "sample_type",
  level = "lowest"
)
tern_table %>%
  mutate(phylum = ifelse(grepl("Pseudomonadota|Actinomycetota|Bacteroidota|Verrucomicrobiota|Patescibacteria|Planctomycetota|Chloroflexota|Nitrospirota", phylum), phylum, "Other")) %>%
  ggtern(., aes(Free_living,Heavy_fraction,Light_fraction)) + 
  facet_wrap(~phylum, nrow = 3, strip.position = "bottom") +
  geom_point(size=4, aes()) +theme_rgbw() +
  labs(title = "Sample_type") +
  theme(legend.position = "none", tern.axis.arrow.text.T=element_text(vjust = -0.5), tern.axis.arrow.text.L=element_text(vjust = -0.5), 
        tern.axis.arrow.text.R=element_text(vjust = 0.5), strip.background = element_blank()) + theme_hidetitles()

##So we see that there are few mOTUs that are unique to the free-living. We see quite some that are enriched for the fractions (right side) 
##Nitrospirota for example appear to be completely missing from the free-living which would point to particles as being critical for nitrogen cycling 


sample_data(bicest_motu_ps) <- sample_data(bicest_motu_ps) %>% 
  as_tibble(rownames="sampleid") %>% 
  mutate(salinity_fac = cut(salinity, breaks = c(0,5,18,40), labels = c("oligohaline", "mesohaline", "polyhaline"), include.lowest = T, right = F)) %>% 
  filter(!is.na(salinity_fac)) %>% 
  column_to_rownames("sampleid") %>% 
  sample_data()

##Great now lets make a ternary plot showing the differences between salinity levels, 
##For this im going to define some categorical levels for salinity 
tern_table <- prep_ternary(
  bicest_motu_ps,
  abund.thres = 1e-04, ##We can adjust this if we feel there are too few/many taxa being plotted
  prev.thres = 0.1, ##We can adjust this if we feel there are too few/many taxa being plotted
  group = "salinity_fac",
  level = "lowest"
)
tern_table %>%
  mutate(phylum = ifelse(grepl("Pseudomonadota|Actinomycetota|Bacteroidota|Verrucomicrobiota|Patescibacteria|Planctomycetota|Chloroflexota|Nitrospirota", phylum), phylum, "Other")) %>%
  ggtern(., aes(oligohaline,polyhaline,mesohaline)) + 
  facet_wrap(~phylum, nrow = 3, strip.position = "bottom") +
  geom_point(size=4, aes()) +theme_rgbw() +
  labs(title = "Sample_type") +
  theme(legend.position = "none", tern.axis.arrow.text.T=element_text(vjust = -0.5), tern.axis.arrow.text.L=element_text(vjust = -0.5), 
        tern.axis.arrow.text.R=element_text(vjust = 0.5), strip.background = element_blank()) + theme_hidetitles()

##Now this is quite different. We have few examples where there are mOTUs common to all sites. Nitrospirota is strictly oligohaline(freshwater)
#JG edit:
#What you will need to do is to link 3-4 things
#You need to link a virus with a mag (iphop) - done, but taxonomy use genus file to be consistent with previous analysis!!!!
host_prediction_genome = read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/Host_prediction_to_genome_m90.csv")#You need to link a mag with a motu 
host_prediction_genus = read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/Host_prediction_to_genus_m90.csv")#You need to link a mag with a motu 
#made lineage for search 
motus_final_taxa_metadata$lineage <- paste(motus_final_taxa_metadata$domain, motus_final_taxa_metadata$phylum, motus_final_taxa_metadata$class, motus_final_taxa_metadata$order, motus_final_taxa_metadata$family,  motus_final_taxa_metadata$genus, motus_final_taxa_metadata$species, sep=";")
motus_hosts_taxa_metadata = motus_final_taxa_metadata[motus_final_taxa_metadata$lineage %in% host_prediction_genome$Host.taxonomy,]

motus_final_taxa_metadata$lineage <- paste(motus_final_taxa_metadata$domain, motus_final_taxa_metadata$phylum, motus_final_taxa_metadata$class, motus_final_taxa_metadata$order, motus_final_taxa_metadata$family,  motus_final_taxa_metadata$genus, sep=";")
motus_hosts_genus_taxa_metadata = motus_final_taxa_metadata[motus_final_taxa_metadata$lineage %in% host_prediction_genus$Host.genus,]

#And then you can take the counts of the motu
#group by phylum
motus_hosts_taxa_metadata_wide <- motus_hosts_taxa_metadata %>%
  filter(station != 'BunthausSpitze') %>% 
  filter(date != "Nov 21") %>%
  select(phylum, sampleid, counts) %>%
  group_by(phylum, sampleid) %>% 
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% 
  pivot_wider(names_from="sampleid", values_from="counts", values_fill=0) 
#read viruses data counts
host_prediction_genus_counts_metadata_mantel <- readRDS(file="C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/host_prediction_genus_counts_metadata_mantel.RDS")
samples_diff <- colnames(motus_hosts_taxa_metadata_wide)[!(colnames(motus_hosts_taxa_metadata_wide) %in% colnames(host_prediction_genus_counts_metadata_mantel))]
samples_diff #zero check if the same
colnames(motus_hosts_taxa_metadata_wide) == colnames(host_prediction_genus_counts_metadata_mantel) #there are equal
motus_hosts_taxa_metadata_wide$phylum == host_prediction_genus_counts_metadata_mantel$phylum 
#there are different lets take subset 
motus_hosts_taxa_metadata_wide$phylum[!(motus_hosts_taxa_metadata_wide$phylum %in% host_prediction_genus_counts_metadata_mantel$phylum)] #motus phyla should be subset - check
#different check taxonomy from two files host_prediction_genome, host_prediction_genus
#now take subset 
subsetphylum = motus_hosts_taxa_metadata_wide$phylum[(motus_hosts_taxa_metadata_wide$phylum %in% host_prediction_genus_counts_metadata_mantel$phylum)] #motus phyla should be subset - check
#take only phyla rom subset now
motus_hosts_taxa_metadata_wide <- motus_hosts_taxa_metadata_wide[motus_hosts_taxa_metadata_wide$phylum %in% subsetphylum,]
host_prediction_genus_counts_metadata_mantel <- host_prediction_genus_counts_metadata_mantel[host_prediction_genus_counts_metadata_mantel$phylum %in% subsetphylum,]
#distance for each 
library("vegan")
saveRDS(motus_hosts_taxa_metadata_wide, file="motus_hosts_taxa_metadata_wide.RDS")
saveRDS(host_prediction_genus_counts_metadata_mantel, file="host_prediction_genus_counts_metadata_mantel.RDS")
motus_hosts_taxa_metadata_wide <- readRDS(file="C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/MOTUs/motus_hosts_taxa_metadata_wide.RDS")
host_prediction_genus_counts_metadata_mantel <- readRDS(, file="C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/MOTUs/host_prediction_genus_counts_metadata_mantel.RDS")

motus_hosts_taxa_metadata_wide$phylum == host_prediction_genus_counts_metadata_mantel$phylum #no there are the same - needed taxonomy from  genus file  
motus.dist <- vegdist(motus_hosts_taxa_metadata_wide[, 2:length(colnames(motus_hosts_taxa_metadata_wide))], "euclidean")
virus.dist <- vegdist(host_prediction_genus_counts_metadata_mantel[, 2:length(colnames(motus_hosts_taxa_metadata_wide))], "euclidean")
#mantel test!!!!
mantel(motus.dist, virus.dist, method = "spearman", permutations = 9999, na.rm = TRUE)#non significant
#check out other phyla 
motus_final_taxa_metadata$lineage <- paste(motus_final_taxa_metadata$domain, motus_final_taxa_metadata$phylum, motus_final_taxa_metadata$class, motus_final_taxa_metadata$order, motus_final_taxa_metadata$family,  motus_final_taxa_metadata$genus, motus_final_taxa_metadata$species, sep=";")

#check other taxonomy assignment from genome files - produce hosts assignment based on "genome" output from iphop - check out why different lines number?   