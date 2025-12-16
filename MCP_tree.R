library("ggplot2")
library("ggtree")
library("ape")
library("dplyr")
library("stringr")
library("RColorBrewer")
library("cowplot")
library("phytools")
library("treeio")
library("ggrepel")
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/MPC_tree")
MCP_annot <- read.csv(file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/MPC_tree/headers_MCP.txt", sep="\t", header=F)
MCP_annot$V1 = gsub("_[1-9]$", "", MCP_annot$V1)
MCP_annot$V1 = gsub(">", "", MCP_annot$V1)
MCP_annot$V1 <- gsub("=", "_", MCP_annot$V1)
# branch lengths as edge labels
#metadata
drep_clusters_fin_taxonomy_unified <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/drep_clusters_fin_taxonomy_unified.RDS")
geneabund_drep_marker_taxonomy_long <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/geneabund_drep_marker_taxonomy_long.RDS")

drep_clusters <- read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])
drep_clusters$genome <- gsub("=", "_", drep_clusters$genome)
drep_clusters_fin_taxonomy_unified_genome <- merge(drep_clusters[c("genome", "primary_cluster")], drep_clusters_fin_taxonomy_unified, by="primary_cluster")
##### MCP_E
#lets plot different trees separately and colour by salinity niche - check again the order of tip coloring 
#and  check if there are any tendencies according to salinity 
tree_MCP_E <- read.iqtree(file="Combined_phagegenes_E_MCP.aln.contree")
#REFSEQ major capsid protein
dup_MCP <- read.csv(file="duplicated.detail_E_MCP.txt", sep='\t', header=F)
dup_MCP$V2 <-  gsub("=", "_", dup_MCP$V2)
#separate into vectors
new_dup_tip <-strsplit(dup_MCP$V2, split= ", ", fixed=TRUE)
#remove without first element (leave one)
for (i in 1:length(new_dup_tip)){
  vec=new_dup_tip[[i]]
  tree_MCP_E <- drop.tip(tree_MCP_E, vec[2:length(vec)])}
# read.iqtree(file)
tree_MCP_E@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_E@phylo$tip.label)
metadata <- as.data.frame(tree_MCP_E@phylo$tip.label)
colnames(metadata) <- c("genome")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata_params <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) %>% select(!c("data_type", "sampleid")) %>% distinct()
metadata <- metadata %>% mutate(AccessionNumber_TBDSven = str_extract(genome, "SAMEA[0-9]+"))
metadata <- left_join(metadata, metadata_params, by="AccessionNumber_TBDSven") 
metadata <- left_join(metadata, MCP_annot, by=c("genome"="V1")) 
#add taxonomic information 
#unify with genome format in data from with taxonomy information 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
tip_na <- metadata[which(is.na(metadata$Salinity_level)),"genome"]
#drop from tree
tree_MCP_E <- drop.tip(tree_MCP_E, tip_na)
#drop from metadata
metadata <- metadata[which(!metadata$genome%in%tip_na),]

colours <- brewer.pal(n = 4, name = "Blues") 
tip_fill <- colours[as.factor(metadata$Salinity_level)]
tip_colors <-tip_fill
tip_colors[which(tip_colors=="#EFF3FF")] <- "black"
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 

# Build base tree
plot_tree_MCP_E <- ggtree(tree_MCP_E,
                                ladderize = T,
                                size = 0.3,
                                #branch.length = "none"
                          ) +
  geom_tippoint(fill = tip_fill, color = tip_colors, shape = 21, size = 2) +
  ggtitle("Major Capsid Protein E")  # + geom_nodelab(aes(label = node, label.size = 1), color="red") + geom_tiplab2(align = TRUE, aes(label=NA), linesize = 0.4) 

plot_tree_MCP_E
# Create df with named rownames
df <- data.frame(Salinity_level = factor(metadata$Salinity_level,
                                         levels = c("freshwater", "oligohaline", "mesohaline", "polyhaline")))
rownames(df) <- tree_MCP_E@phylo$tip.label
colours_named <- setNames(colours, levels(df$Salinity_level))


# Define your clade data
cladeda <- data.frame(id = c(601, 852, 989), type =c("lower salinity", "higher salinity", "lower salinity") )
p3 <- plot_tree_MCP_E +
  geom_cladelab(
    data = cladeda,
    mapping = aes(node = id, label = type),
    textcolour = "black",
    barsize = 1,
    color = "black",
    extend = 0.2,
    clip="off"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 0),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  ) + hexpand(ratio = 0.05)
#hilights behind the tree so they are visible 
plot_tree_MCP_E_annot <- p3 +
   geom_hilight(data=cladeda, aes(node=id, fill=type),type = "roundrect", fill="#C8BFBF") + geom_tree(linewidth=0.5)  + 
  geom_tippoint(fill = tip_fill, color = tip_colors, shape = 21, size = 2)  + coord_cartesian(clip="off")

tiff("tree_MCP_E_annot.tiff" )
plot_tree_MCP_E_annot
dev.off()
#REFSEQ major capsid protein
dup_MCP <- read.csv(file="duplicated.detail_REFSEQ_MCP.txt", sep='\t', header=F)
dup_MCP$V2 <-  gsub("=", "_", dup_MCP$V2)

#separate into vectors
new_dup_tip <-strsplit(dup_MCP$V2, split= ", ", fixed=TRUE)
tree_MCP_REFSEQ <- read.iqtree(file="Combined_phagegenes_REFSEQ_MCP.aln.contree")

#remove without first element (leave one)
for (i in 1:length(new_dup_tip)){
  vec=new_dup_tip[[i]]
  tree_MCP_REFSEQ <- drop.tip(tree_MCP_REFSEQ, vec[2:length(vec)])}

tree_MCP_REFSEQ@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_REFSEQ@phylo$tip.label)

metadata <- as.data.frame(tree_MCP_REFSEQ@phylo$tip.label)
colnames(metadata) <- c("genome")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata_params <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) 
#parameters for the same accession number should be the same - remove data type column and keep unique rows only 
#so that we dont end up with many-to-many relationship
metadata_params <- unique(metadata_params[, colnames(metadata_params)[which(colnames(metadata_params)!="data_type" & colnames(metadata_params)!="sampleid")]])
metadata <- metadata %>% mutate(AccessionNumber_TBDSven = str_extract(genome, "SAMEA[0-9]+"))
metadata <- left_join(metadata, metadata_params, by="AccessionNumber_TBDSven") 
metadata <- left_join(metadata, MCP_annot, by=c("genome"="V1")) 
#add taxonomic information 
#unify with genome format in data from with taxonomy information 
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
colours <- brewer.pal(n = 4, name = "Blues") 
tip_colors <- colours[as.factor(metadata$Salinity_level)]

library(ggnewscale)

#"sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17"
#read duplicated sequences
dup_MCP <- read.csv(file="duplicated.detail_gp17_DNA_MCP.txt", sep='\t', header=F)
dup_MCP$V2 <-  gsub("=", "_", dup_MCP$V2)

#separate into vectors
new_dup_tip <-strsplit(dup_MCP$V2, split= ", ", fixed=TRUE)


tree_MCP_gp17 <- read.iqtree(file="Combined_phagegenes_gp17_MCP.aln.contree")
#remove without first element (leave one)
for (i in 1:length(new_dup_tip)){
  vec=new_dup_tip[[i]]
  tree_MCP_gp17 <- drop.tip(tree_MCP_gp17, vec[2:length(vec)])}

tree_MCP_gp17@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_gp17@phylo$tip.label)

metadata <- as.data.frame(tree_MCP_gp17@phylo$tip.label)
colnames(metadata) <- c("genome")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata_params <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) 
#parameters for the same accession number should be the same - remove data type column and keep unique rows only 
#so that we dont end up with many-to-many relationship
metadata_params <- unique(metadata_params[, colnames(metadata_params)[which(colnames(metadata_params)!="data_type" & colnames(metadata_params)!="sampleid")]])
metadata <- metadata %>% mutate(AccessionNumber_TBDSven = str_extract(genome, "SAMEA[0-9]+"))
metadata <- left_join(metadata, metadata_params, by="AccessionNumber_TBDSven") 
metadata <- left_join(metadata, MCP_annot, by=c("genome"="V1")) 
#add taxonomic information 
#unify with genome format in data from with taxonomy information 
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
colours <- brewer.pal(n = 4, name = "Blues") 
tip_colors <- colours[as.factor(metadata$Salinity_level)]

plot_tree_MCP_gp17 <- ggtree(tree_MCP_gp17, ladderize = T, size = 0.3 , branch.length="none") + 
  geom_tippoint(fill=tip_colors, color=tip_colors, shape=21, size=3) + 
  geom_tiplab2(align = TRUE, aes(label=NA), linesize = 0.4) +
  geom_nodelab(aes(label = node), size = 2)

plot_tree_MCP_gp17
#hard to find big clusters of high /low salinity levels 


############################
#"Large eukaryotic DNA virus major capsid protein" 
dup_MCP <- read.csv(file="duplicated.detail_euk_DNA_MCP.txt", sep='\t', header=F)
dup_MCP$V2 <-  gsub("=", "_", dup_MCP$V2)

#separate into vectors
new_dup_tip <-strsplit(dup_MCP$V2, split= ", ", fixed=TRUE)


tree_MCP_euk_DNA <- read.iqtree(file="Combined_phagegenes_euk_DNA_MCP.aln.contree")
#remove without first element (leave one)
for (i in 1:length(new_dup_tip)){
  vec=new_dup_tip[[i]]
  print(length(vec))
  tree_MCP_euk_DNA <- drop.tip(tree_MCP_euk_DNA, vec[2:length(vec)])
}
tree_MCP_euk_DNA@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_euk_DNA@phylo$tip.label)

metadata <- as.data.frame(tree_MCP_euk_DNA@phylo$tip.label)
colnames(metadata) <- c("genome")
metadata <- metadata %>% mutate(AccessionNumber_TBDSven = str_extract(genome, "SAMEA[0-9]+"))
metadata <- left_join(metadata, metadata_params, by="AccessionNumber_TBDSven") 
metadata <- left_join(metadata, MCP_annot, by=c("genome"="V1")) 
#add taxonomic information 
#unify with genome format in data fram with taxonomy information 
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
#remove proteins without salinity data 
tip_na <- metadata[ which(is.na(metadata$Salinity_level)),"genome"]
#drop from tree
tree_MCP_euk_DNA <- drop.tip(tree_MCP_euk_DNA, tip_na)
#drop from metadata
metadata <- metadata[ which(metadata$genome!=tip_na),]

colours <- brewer.pal(n = 4, name = "Blues") 
tip_fill <- colours[as.factor(metadata$Salinity_level)]
tip_colors <-tip_fill
tip_colors[which(tip_colors=="#EFF3FF")] <- "black"
tree_MCP_euk_DNA@phylo$tip.label <- c(1:length(tree_MCP_euk_DNA@phylo$tip.label))
rownames(metadata) <- tree_MCP_euk_DNA@phylo$tip.label
#find tip without salinity information
plot_tree_MCP_euk_DNA <- ggtree(tree_MCP_euk_DNA,
                                ladderize = T,
                                size = 0.3,
                                #branch.length = "none"
) +
  geom_tippoint(fill = tip_fill, color = tip_colors, shape = 21, size = 2) +
  ggtitle("Large eukaryotic DNA virus major capsid protein") #+ geom_nodelab(aes(label = node, label.size = 1), color="red") + geom_tiplab2(align = TRUE, aes(label=NA), linesize = 0.4) 

plot_tree_MCP_euk_DNA
# Create df with named rownames
df <- data.frame(Salinity_level = factor(metadata$Salinity_level,
                                         levels = c("freshwater", "oligohaline", "mesohaline", "polyhaline")))
rownames(df) <- tree_MCP_euk_DNA@phylo$tip.label
colours_named <- setNames(colours, levels(df$Salinity_level))

# Define your clade data
cladeda <- data.frame(id = c(123,   212, 207, 170, 165, 140), type =c("lower salinity", "higher salinity", "lower salinity", "higher salinity","higher salinity","lower salinity"))
# 
p3 <- plot_tree_MCP_euk_DNA +
  geom_cladelab(
    data = cladeda,
    mapping = aes(node = id, label = type),
    textcolour = "black",
    barsize = 1,
    color = "black",
    extend = 0.2,
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 0),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )
#hilights behind the tree so they are visible 
plot_tree_MCP_euk_DNA_annot <- p3 +
  geom_hilight(data=cladeda, aes(node=id, fill=type),type = "roundrect", fill="#C8BFBF") + geom_tree(linewidth=0.5)  + 
  geom_tippoint(fill = tip_fill, color = tip_colors, shape = 21, size = 2) 
tiff("tree_MCP_euk_DNA_annot.tiff" )
plot_tree_MCP_euk_DNA_annot
dev.off()

#combine together with IP 
IP_hist_all <- readRDS(file="~/IGB_phd/BICEST/virus/viral_abundances/IP_hist_all.rds")
IP_hist_MCP <- readRDS(file="~/IGB_phd/BICEST/virus/viral_abundances/IP_hist_MCP.rds")
first_row <- plot_grid(IP_hist_MCP, IP_hist_all,rel_widths = c(1,2), labels = c("A", "B"))
second_row <-plot_grid(plot_tree_MCP_E_annot, plot_tree_MCP_euk_DNA_annot, labels = c("C", "D"))
plot_MCP_all<-plot_grid( first_row, second_row, label_fontfamily ="sans", rel_heights = c(1,3), rows=2)
tiff("plot_MCP_all.tiff",unit="cm", width = 32, height=44, res=300)
plot_MCP_all
dev.off()
svg("plot_MCP_all.svg", width = 16, height=22)
plot_MCP_all
dev.off()

pdf(file="plot_MCP_all.pdf", width = 16, height=22)
plot_MCP_all
dev.off()


