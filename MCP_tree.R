library("ggplot2")
library("ggtree")
library("ape")
library("dplyr")
library(stringr)
library(RColorBrewer)

setwd("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/")
nwk <- ("Combined_phagegenes_GT_MCP.aln.treefile")
MCP_annot <- read.csv(file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/MPC_tree/headers_MCP.txt", sep="\t", header=F)
MCP_annot$V1 = gsub("_[1-9]$", "", MCP_annot$V1)
MCP_annot$V1 = gsub(">", "", MCP_annot$V1)
MCP_annot$V1 <- gsub("=", "_", MCP_annot$V1)
tree <- read.tree(nwk)
#metadata
drep_clusters_fin_taxonomy_unified <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/drep_clusters_fin_taxonomy_unified.RDS")
geneabund_drep_marker_taxonomy_long <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/geneabund_drep_marker_taxonomy_long.RDS")
# tiff("MCPtree.tiff",    width = 1080, height = 10080, units = "px", )
# plot(tree, tip.color = "black")
# dev.off()
# #
drep_clusters <- read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])
# tree.tip_target = gsub("_[1-9]$", "", tree$tip.label)
drep_clusters_fin_taxonomy_unified_genome <- merge(drep_clusters[c("genome", "primary_cluster")], drep_clusters_fin_taxonomy_unified, by="primary_cluster")
#tip label change with sub _nr
# drep_clusters_fin_taxonomy_unified_genome <- drep_clusters_fin_taxonomy_unified_genome %>% arrange(factor(genome, levels = tree.tip_target))
#unify with names in the tree
# drep_clusters_fin_taxonomy_unified_genome$genome <- gsub("=", "_", drep_clusters_fin_taxonomy_unified_genome$genome)

#modyfy at tree 
#tree$tip.label <- tree.tip_target
metadata <- as.data.frame(tree$tip.label)
colnames(metadata) <- c("genome")
# metadata <- left_join(metadata, drep_clusters_fin_taxonomy_unified_genome, by="genome")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/metadata/"
metadata_params <- read.csv(paste0(path,"/PhysicochemicalParameters_mod3.csv"), header=TRUE, sep=",", row.names=1) 
#parameters for the same accession number should be the same - remove data type column and keep unique rows only 
#so that we dontend up with many-to-many relationship
metadata_params <- unique(metadata_params[, colnames(metadata_params)[which(colnames(metadata_params)!="data_type" & colnames(metadata_params)!="sampleid")]])
metadata <- metadata %>% mutate(AccessionNumber_TBDSven = str_extract(genome, "SAMEA[0-9]+"))

metadata <- left_join(metadata, metadata_params, by="AccessionNumber_TBDSven") 
metadata <- left_join(metadata, MCP_annot, by=c("genome"="V1")) 
metadata <- metadata %>% mutate(Salinity_level = case_when(Salinity_PSU <= .5  ~ 'freshwater',  Salinity_PSU > .5 & Salinity_PSU <= 5.5  ~ 'oligohaline', Salinity_PSU > 5.5 & Salinity_PSU <= 18  ~ 'mesohaline', Salinity_PSU > 18  ~ 'polyhaline'))
MCP_annot$meta
metadata$Salinity_level <-factor(metadata$Salinity_level, levels=c( "freshwater", "oligohaline", "mesohaline",  "polyhaline" ))
colours <- brewer.pal(n = 4, name = "Blues") 

p1 <- ggtree(tree, aes(color=metadata$Salinity_level)) + theme_tree2() 
tiff("MCPtree.tiff",    width = 1080, height = 10080, units = "px", )
plot(tree, tip.color = colours[as.factor(metadata$Salinity_level)])
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()

coloursSet1 <-  brewer.pal(n=6, name = "Set1")
tiff("MCPtree_clades.tiff",    width = 1080, height = 10080, units = "px", )
plot(tree, tip.color = coloursSet1[as.factor(metadata$V5)])
nodelabels(cex=0.5)
dev.off()

findMRCA(tree, tree$tip.label[c(1,4)])
#lets plot different trees separatly and colour by salinity niche - check again the order of tip coloring 
#and  check if there are any tendencies according to salinity 
rows_MCP_E <- which(metadata$V5=="Phage major capsid protein E" & !is.na(metadata$Salinity_level))
tree_MCP_E <- keep.tip(tree, tip=metadata$genome[rows_MCP_E])
#sanity check
plot(tree_MCP_E,  tip.color = colours[as.factor(metadata$V5[which(metadata$V5=="Phage major capsid protein E"  & !is.na(metadata$Salinity_level))])])
legend("bottomleft",
       legend = levels(metadata$V5), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()    

tip_colours <- colours[as.factor(metadata$Salinity_level[rows_MCP_E])]                                         
tiff("tree_MCP_E.tiff")
plot(tree_MCP_E, tip.color=tip_colours)
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
# rest mostly oligohaline 
#do the same for rest 
#REFSEQ major capsid protein
rows_MCP_REFSEQ <- which(metadata$V5=="REFSEQ major capsid protein" & !is.na(metadata$Salinity_level))
tree_MCP_REFSEQ <- keep.tip(tree, tip=metadata$genome[rows_MCP_REFSEQ])
tip_colours <- colours[as.factor(metadata$Salinity_level[rows_MCP_REFSEQ])]  
tiff("tree_MCP_REFSEQ.tiff")
plot(tree_MCP_REFSEQ, tip.color=tip_colours )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#"sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17"
rows_MCP_REFSEQ <- which(metadata$V5=="sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17" & !is.na(metadata$Salinity_level))
tree_MCP_gp17 <- keep.tip(tree, tip=metadata$genome[rows_MCP_REFSEQ])
tip_colours <- colours[as.factor(metadata$Salinity_level[rows_MCP_REFSEQ])]     
tiff("tree_MCP_gp17.tiff")
plot(tree_MCP_gp17, tip.color=tip_colours )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#"Large eukaryotic DNA virus major capsid protein" 
rows_MCP_large_DNA  <- which(metadata$V5=="Large eukaryotic DNA virus major capsid protein" & !is.na(metadata$Salinity_level))
tree_MCP_large_DNA <- keep.tip(tree, tip=metadata$genome[rows_MCP_large_DNA])
tip_colours <- colours[as.factor(metadata$Salinity_level[rows_MCP_large_DNA])]                                         
tiff("tree_large_DNA.tiff")
plot(tree_MCP_large_DNA, tip.color=tip_colours )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#clear  oligohaline, mesohaline clusters

#"REFSEQ putative major capsid protein"  - too little nodes
rows_MCP_REFSEQ_putative  <- which(metadata$V5=="REFSEQ putative major capsid protein" & !is.na(metadata$Salinity_level))
tree_MCP_REFSEQ_putative <- keep.tip(tree, tip=metadata$genome[rows_MCP_REFSEQ_putative])
tip_colours <- colours[as.factor(metadata$Salinity_level[rows_MCP_REFSEQ_putative])]                                         
tiff("tree_MCP_REFSEQ_putative.tiff")
plot(tree_MCP_REFSEQ_putative, tip.color=tip_colours )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#freshwater only!
#check order 
metadata$genome[which(metadata$V5=="REFSEQ putative major capsid protein")]
#tiff("MCPtree.tiff",    width = 1080, height = 10080, units = "px", )
small_metadata <- metadata[which(!is.na(metadata$Salinity_niche)),]
#take only one per cluster so we dont get artificial clades
small_metadata <- metadata %>% group_by(primary_cluster) %>% slice_head() %>% ungroup()
small_metadata$Fraction_niche <-as.factor(small_metadata$Fraction_niche)
small_tree <- keep.tip(tree, tip=small_metadata$genome[which(!is.na(small_metadata$Salinity_niche))])
metadata <- as.data.frame(tree.tip_target)
as.factor(small_metadata$Salinity_niche)
colours=c(brewer.pal(n = length(unique(drep_clusters_fin_taxonomy_unified$Salinity_niche))-2, name = "Blues"), "yellow", "orange")

tiff("salinity_tree.tiff")
plot(small_tree, tip.color = colours[as.factor(small_metadata$Salinity_niche)])
legend("bottomleft",
       legend = levels(small_metadata$Salinity_niche), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
dev.off()


tiff("fraction_tree.tiff")
colours=c(brewer.pal(n = length(unique(small_metadata$Fraction_niche)), name = "Set2"))
plot(small_tree, tip.color = colours[as.factor(small_metadata$Fraction_niche)])
legend("bottomleft",
       legend = levels(small_metadata$Fraction_niche), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
dev.off()

colours=c(brewer.pal(n = length(unique(drep_clusters_fin_taxonomy_unified$Fraction_niche)), name = "Set1"))

plot(small_tree, tip.color = colours[as.factor(small_metadata$level_5)])
