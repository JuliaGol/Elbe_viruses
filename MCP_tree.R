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
#detach("package:phyloseq", unload=TRUE) #conflict with ape
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

#lets plot different trees separately and colour by salinity niche - check again the order of tip coloring 
#and  check if there are any tendencies according to salinity 
tree_MCP_E <- read.iqtree(file="Combined_phagegenes_E_MCP.aln.contree")
# read.iqtree(file)
tree_MCP_E@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_E@phylo$tip.label)

metadata <- as.data.frame(tree_MCP_E@phylo$tip.label)
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
#unify with genome format in data fram with taxonomy information 
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
colours <- brewer.pal(n = 4, name = "Blues") 
tip_colors <- colours[as.factor(metadata$Salinity_level)]

colors <- tip_colors

plot_tree_MCP_E <- ggtree(tree_MCP_E, ladderize = FALSE ) + geom_tippoint(fill=tip_colors, color=colors, shape=21, size=2) # + geom_tiplab(fill = tip_colors,
                                                                                #  color = NA, geom = "label",
                                                                                #  label.padding = unit(0.15, "lines"), # amount of padding around the labels
                                                                                #  label.size = 0, align = F)  # size of label border
tiff("tree_MCP_E_sep.tiff" )
tree_MCP_E_sep <- plot_tree_MCP_E + coord_cartesian(clip="off") + ggtitle("Phage major capsid protein E") + theme(plot.title = element_text(hjust = 0.5)) + geom_nodelab(aes(subset = node %in% c(1054,1462,1677,1463, 1461, 1055, 1045)), size = 5) + theme(panel.background = element_rect(fill="lightgrey", size =1.5))
dev.off()

# rest mostly oligohaline 
#do the same for rest 
#REFSEQ major capsid protein
tree_MCP_REFSEQ <- read.iqtree(file="Combined_phagegenes_REFSEQ_MCP.aln.contree")
# read.iqtree(file)
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


plot_tree_MCP_REFSEQ <- ggtree(tree_MCP_REFSEQ, ladderize = FALSE ) + geom_tippoint(fill=tip_colors, color=tip_colors, shape=21, size=2) #+ geom_tiplab(fill = tip_colors,
tree_MCP_REFSEQ_sep_nodes <- plot_tree_MCP_REFSEQ + coord_cartesian(clip="off") + ggtitle("REFSEQ major capsid protein") + theme(plot.title = element_text(hjust = 0.5)) + geom_nodelab(aes(label = node, subset = as.numeric(node)>300 ), size = 1) + theme(panel.background = element_rect(fill="lightgrey", size =1.5))
tiff("tree_MCP_REFSEQ_sep.tiff" )
tree_MCP_REFSEQ_sep <- plot_tree_MCP_REFSEQ + coord_cartesian(clip="off") + ggtitle("REFSEQ major capsid protein") + theme(plot.title = element_text(hjust = 0.5)) + geom_nodelab(aes(subset = node %in% c()), size = 5) + theme(panel.background = element_rect(fill="lightgrey", size =1.5))
dev.off()


#"sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17"
tree_MCP_gp17 <- read.iqtree(file="Combined_phagegenes_gp17_MCP.aln.contree")
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

plot_tree_MCP_gp17 <- ggtree(tree_MCP_gp17, ladderize = FALSE ) + geom_tippoint(fill=tip_colors, color=tip_colors, shape=21, size=2) # + geom_tiplab(fill = tip_colors,
tiff("tree_MCP_gp17_sep.tiff" )
tree_MCP_gp17_sep <- plot_tree_MCP_gp17 + coord_cartesian(clip="off") + ggtitle("Probable major capsid protein gp17") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.background = element_rect(fill="lightgrey", size =1.5))
dev.off()
#"Large eukaryotic DNA virus major capsid protein" 
tree_MCP_euk_DNA <- read.iqtree(file="Combined_phagegenes_euk_DNA_MCP.aln.contree")
tree_MCP_euk_DNA@phylo$tip.label <-  gsub("_fragment.*$", "", tree_MCP_euk_DNA@phylo$tip.label)

metadata <- as.data.frame(tree_MCP_euk_DNA@phylo$tip.label)
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
#unify with genome format in data fram with taxonomy information 
metadata$genome = gsub("_[1-9]$", "", metadata$genome)
metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome[,c("genome", "level_5")], by="genome") 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
colours <- brewer.pal(n = 4, name = "Blues") 
tip_colors <- colours[as.factor(metadata$Salinity_level)]

plot_tree_MCP_euk_DNA <- ggtree(tree_MCP_euk_DNA, ladderize = FALSE ) + geom_tippoint(fill=tip_colors, color=tip_colors, shape=21, size=3)  #+ geom_nodelab(aes(label = node, label.size = 1)) + geom_tiplab(fill = tip_colors,

tiff("tree_MCP_euk_DNA_sep.tiff" )

tree_MCP_euk_DNA_sep <- plot_tree_MCP_euk_DNA + coord_cartesian(clip="off") + ggtitle("Large eukaryotic DNA virus major capsid protein") + theme(plot.title = element_text(hjust = 0.5)) + geom_nodelab(aes(subset = node %in% c(306, 274, 225, 236, 241), label.size = 3)) + theme(panel.background = element_rect(fill="lightgrey", size =1.5))
dev.off()
#clear  oligohaline, mesohaline clusters
#put all together 
all_MPC_trees <- plot_grid(tree_MCP_E_sep, tree_MCP_euk_DNA_sep, tree_MCP_REFSEQ_sep, tree_MCP_gp17_sep,  ncol=2, align="none", rel_heights = c(1,2))
tiff("MPC_phylogeny_grid_tippoint.tiff", width = 1000, height=1000)
all_MPC_trees
dev.off()
