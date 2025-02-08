library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/Elbe_data_overview")
tiara_stats <- read.csv("tiara_out_all_clean.txt", header = F, sep ="\t")
colnames(tiara_stats) <- c("sequence_id", "class_fst_stage", "class_snd_stage", "org", "bac", "arc", "euk", "unk1",  "pla", "unk2",  "mit")
tiara_stats[,"sample_id"] = substr(tiara_stats$sequence_id, 10, 23)
tiara_stats[,"sample_type"] = substr(tiara_stats$sequence_id, 25, 29)
tiara_stats %>% group_by(sample_id, class_fst_stage) %>% summarise(count_fst_classes = n())
summary_class_fst <- tiara_stats %>% count(sample_id, class_fst_stage)
summary_class_snd <- tiara_stats %>% count(sample_id, class_snd_stage)
colnames(summary_class_fst) <- c("sample_id", "class", "counts")
colnames(summary_class_snd) <- c("sample_id", "class", "counts")
summary_class_fst <- summary_class_fst[which(summary_class_fst$class != "organelle"),] 
summary_class_snd[which(summary_class_snd$class == "unknown"),"class"] <- "unknown_org"
summary_class_snd <- summary_class_snd[which(summary_class_snd$class != "n/a"),] 
summary_class <- rbind(summary_class_fst, summary_class_snd)
#abundance
summary_class %>% ggplot(aes(x=sample_id, y=counts, fill = class)) + geom_bar(position="stack", stat="identity")  
#rel abundance
summary_class <- summary_class %>% group_by(sample_id) %>% 
  mutate(percent=round(counts/sum(counts)*100, 2)) %>% 
  ungroup
summary_class %>% ggplot(aes(x=sample_id, y=percent, fill = class)) + geom_bar(position="stack", stat="identity")  
#some samples have enriched eukaryots
#check if this is seasonal or there is other factor
#read metadata 

metadata <- read.csv("~/IGB_phd/virus/viral_abundances/SAMEAID_SampleID_simplified.csv", header=TRUE, sep=";") %>% 
  mutate(sampleid=paste0(ProjectID,"_",BioSample,"_METAG.genecount.profile"))%>%
  mutate(data_type="METAG")

metadata <-metadata %>% 
  mutate(sample=str_remove(sampleid, ".genecount.profile")) %>%
  mutate(sample=str_replace(sample, "[.]", "-")) %>% filter(Station != "Seemanshöft") %>% filter(Station != "Kollmar")

metadata[,"Station"] <- factor(metadata$Station, levels=c( "BunthausSpitze", "Muehlenberger Loch" ,"Schwarztonnensand", "Twielenfleth", "Brunsbuettel", "Meedem Grund"))
summary_class <- summary_class %>% left_join(metadata , by = c("sample_id" = "BioSample"))

summary_class %>% ggplot(aes(x=sample_id, y=percent, fill = class)) + geom_bar(position="stack", stat="identity") + facet_grid(~ Sample_date)
summary_class %>% ggplot(aes(x=sample_id, y=percent, fill = class)) + geom_bar(position="stack", stat="identity") + facet_grid(~ Sample_type)
summary_class %>% ggplot(aes(x=sample_id, y=percent, fill = class)) + geom_bar(position="stack", stat="identity") + facet_grid(~ Station)
summary_class %>% ggplot(aes(x=sample_id, y=percent, fill = class)) + geom_bar(position="stack", stat="identity") + facet_grid(~ StationNumber)
tiff("tiara_bicest_2021.tiff")
p1 <- summary_class %>% ggplot(aes(x=Station, y=percent, fill = class)) + geom_bar(position="stack", stat="identity") +
  ylab("no. contings") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(~ Sample_date)
dev.off()
#more in May and PA fractions Bunthaus and Meedemgrund - beginning and the end of the estuary