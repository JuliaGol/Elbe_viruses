library("ggplot2")
library("ggtree")
library("ape")
library("dplyr")
library("stringr")
library("RColorBrewer")
library("cowplot")
#detach("package:phyloseq", unload=TRUE) #conflict with ape
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/MPC_tree")
nwk <- ("Combined_phagegenes_GT_MCP.aln.treefile")
MCP_annot <- read.csv(file="C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes/MPC_tree/headers_MCP.txt", sep="\t", header=F)
MCP_annot$V1 = gsub("_[1-9]$", "", MCP_annot$V1)
MCP_annot$V1 = gsub(">", "", MCP_annot$V1)
MCP_annot$V1 <- gsub("=", "_", MCP_annot$V1)
tree <- ape::read.tree(nwk)
#metadata
drep_clusters_fin_taxonomy_unified <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/drep_clusters_fin_taxonomy_unified.RDS")
geneabund_drep_marker_taxonomy_long <- readRDS("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/CCA_marker_gene/geneabund_drep_marker_taxonomy_long.RDS")

drep_clusters <- read.csv("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_abundances/Cdb.csv", header = TRUE, sep = ",")
drep_clusters[,"genome"] = sub(".fna", "", drep_clusters[,"genome"])
drep_clusters$genome <- gsub("=", "_", drep_clusters$genome)
drep_clusters_fin_taxonomy_unified_genome <- merge(drep_clusters[c("genome", "primary_cluster")], drep_clusters_fin_taxonomy_unified, by="primary_cluster")

#modyfy at tree 
metadata <- as.data.frame(tree$tip.label)
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
#visualise clades of different MCPs
coloursSet1 <-  brewer.pal(n=6, name = "Set1")
tiff("MCPtree_clades.tiff",    width = 1080, height = 10080, units = "px", )
plot(tree, tip.color = coloursSet1[as.factor(metadata$V5)])
nodelabels(cex=0.5)
dev.off()

#lets plot different trees separately and colour by salinity niche - check again the order of tip coloring 
#and  check if there are any tendencies according to salinity 
tree_MCP_E <- keep.tip(tree, tip=metadata$genome[rows_MCP_E])
rows_MCP_E <- which(metadata$V5=="Phage major capsid protein E" & !is.na(metadata$Salinity_level))
tip_colors <- colours[as.factor(metadata$Salinity_level[rows_MCP_E])]
plot_tree_MCP_E <- ggtree(tree_MCP_E, ladderize = FALSE ) + geom_tiplab(fill = tip_colors,,
                                                                                  color = "black", geom = "label", 
                                                                                  label.padding = unit(0.15, "lines"), # amount of padding around the labels
                                                                                  label.size = 0, align = TRUE)  # size of label border
plot_tree_MCP_E + ggtitle("Phage major capsid protein E") + theme(plot.title = element_text(hjust = 0.5))
  


tiff("tree_MCP_E.tiff")
plot(tree_MCP_E, tip.color=tip_colors)
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
tip_colors <- colours[as.factor(metadata$Salinity_level[rows_MCP_REFSEQ])] 

plot_tree_MCP_REFSEQ <- ggtree(tree_MCP_REFSEQ, ladderize = FALSE ) + geom_tiplab(fill = tip_colors,,
                                                                                  color = "black", geom = "label", 
                                                                                  label.padding = unit(0.15, "lines"), # amount of padding around the labels
                                                                                  label.size = 0, align = TRUE)  + ggtitle("REFSEQ major capsid protein") + theme(plot.title = element_text(hjust = 0.5))

tiff("tree_MCP_REFSEQ.tiff")
plot(tree_MCP_REFSEQ, tip.color=tip_colors )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#"sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17"
rows_MCP_gp17 <- which(metadata$V5=="sp|O64210|CAPSD_BPMD2 Probable major capsid protein gp17" & !is.na(metadata$Salinity_level))
tree_MCP_gp17 <- keep.tip(tree, tip=metadata$genome[rows_MCP_gp17])
tip_colors <- colours[as.factor(metadata$Salinity_level[rows_MCP_gp17])]     

plot_tree_MCP_gp17 <- ggtree(tree_MCP_gp17, ladderize = FALSE ) + geom_tiplab(fill = tip_colors,
                                                                        color = "black", geom = "label", 
                                                                        label.padding = unit(0.15, "lines"), # amount of padding around the labels
                                                                        label.size = 0, size=0.5, align = TRUE)  + ggtitle("Probable major capsid protein gp17") + theme(plot.title = element_text(hjust = 0.5))# size of label border

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
tip_colors <- colours[as.factor(metadata$Salinity_level[rows_MCP_large_DNA])]    
plot_tree_MCP_large_DNA <- ggtree(tree_MCP_large_DNA, ladderize = FALSE ) + geom_tiplab(fill = tip_colors,
                                                                              color = "black", geom = "label", 
                                                                              label.padding = unit(0.15, "lines"), # amount of padding around the labels
                                                                              label.size = 0, align = TRUE)  # size of label border
plot_tree_MCP_large_DNA <-  plot_tree_MCP_large_DNA + ggtitle("Large eukaryotic DNA virus major capsid protein") + theme(plot.title = element_text(hjust = 0.5))

tiff("tree_large_DNA.tiff")
plot(tree_MCP_large_DNA, tip.color=tip_colors )
legend("bottomleft",
       legend = levels(metadata$Salinity_level), #Name of groups
       fill = colours,       # Color of the squares
       border = "black", # Color of the border of the squares
       cex =.6, #sets legend size
       xpd= TRUE) #places outside plot area
nodelabels()
dev.off()
#clear  oligohaline, mesohaline clusters
#put all together 
all_MPC_trees <- plot_grid(plot_tree_MCP_E, plot_tree_MCP_REFSEQ, plot_tree_MCP_gp17, plot_tree_MCP_large_DNA, ncol=2, aligh = "v")
all_MPC_trees
#visualise with heatmap 
#REFSEQ 
df <- as.data.frame(metadata[rows_MCP_REFSEQ, c("Salinity_level")])
colnames(df) <- c("Salinity_level")
rownames(df) <- tree_MCP_REFSEQ$tip.label
df2 <- as.data.frame(metadata[rows_MCP_REFSEQ, c("Sample_type")])
colnames(df2) <- c("Sample_type")
rownames(df2) <- tree_MCP_REFSEQ$tip.label 
plot_tree_MCP_REFSEQ <- ggtree(tree_MCP_REFSEQ, 
                               ladderize = FALSE)
colours <- c("freshwater" = "#EFF3FF", "oligohaline" = "#BDD7E7", "mesohaline" = "#6BAED6", "polyhaline"= "#2171B5") 
p1 <- gheatmap(plot_tree_MCP_REFSEQ,df, colnames = FALSE, color=NA) + scale_fill_manual(values=colours)  + ggtitle("REFSEQ major capsid protein") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
library(ggnewscale)
p12 <- p1 + new_scale_fill()
gheatmap(p12, df2, colnames = FALSE, offset=5, color=NA) + scale_fill_brewer(palette = "Set1") 

#E
df <- as.data.frame(metadata[rows_MCP_E, "Salinity_level"])
colnames(df) <- c("Salinity_level")
rownames(df) <- tree_MCP_E$tip.label

plot_tree_MCP_E <- ggtree(tree_MCP_E, 
                               ladderize = FALSE)
p2 <- gheatmap(plot_tree_MCP_E, df,  colnames = FALSE, color=NA) + scale_fill_manual(values=colours) + ggtitle("Phage major capsid protein E") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
df2 <- as.data.frame(metadata[rows_MCP_E, c("Sample_type")])
colnames(df2) <- c("Sample_type")
rownames(df2) <- tree_MCP_E$tip.label 
p22 <- p2 + new_scale_fill()
gheatmap(p22, df2, colnames = FALSE, offset=5, color=NA) + scale_fill_brewer(palette = "Set1") 

#_large_DNA
df <- as.data.frame(metadata[rows_MCP_large_DNA, "Salinity_level"])
colnames(df) <- c("Salinity_level")
rownames(df) <- tree_MCP_large_DNA$tip.label

plot_tree_MCP_large_DNA <- ggtree(tree_MCP_large_DNA, 
                               ladderize = FALSE)
colours <- c("freshwater" = "#EFF3FF", "oligohaline" = "#BDD7E7", "mesohaline" = "#6BAED6", "polyhaline"= "#2171B5") 
df2 <- as.data.frame(metadata[rows_MCP_large_DNA, c("Sample_type")])
colnames(df2) <- c("Sample_type")
rownames(df2) <- tree_MCP_large_DNA$tip.label 

p3 <- gheatmap(plot_tree_MCP_large_DNA, df, colnames = FALSE, color=NA) + scale_fill_manual(values=colours) + ggtitle("Large eukaryotic DNA virus major capsid protein") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p32 <- p3 + new_scale_fill()
gheatmap(p32, df2, colnames = FALSE, offset=5, color=NA) + scale_fill_brewer(palette = "Set1") 

#gp_17
df <- as.data.frame(metadata[rows_MCP_gp17, "Salinity_level"])
#df$Salinity_level <- factor(df$Salinity_level , levels=c("freshwater", "oligohaline", "mesohaline", "polyhaline"))
colnames(df) <- c("Salinity_level")
rownames(df) <- tree_MCP_gp17$tip.label

plot_tree_MCP_gp17 <- ggtree(tree_MCP_gp17, 
                                  ladderize = FALSE)
df2 <- as.data.frame(metadata[rows_MCP_gp17, c("Sample_type")])
colnames(df2) <- c("Sample_type")
rownames(df2) <- tree_MCP_gp17$tip.label 

colours <- c("freshwater" = "#EFF3FF", "oligohaline" = "#BDD7E7", "mesohaline" = "#6BAED6", "polyhaline"= "#2171B5") 
p4 <- gheatmap(plot_tree_MCP_gp17, df, offset=0, colnames = FALSE, color=NA) + scale_fill_manual(values=colours) + ggtitle("Probable major capsid protein gp17") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 
p42 <- p4 + new_scale_fill()
gheatmap(p42, df2, colnames = FALSE, offset=5, color=NA) + scale_fill_brewer(palette = "Set1") 

tiff("MPC_phylogeny_grid.tiff", width = 1000, height=1000)
plot_grid(p1,p2,p3,p4, ncol=2, align ="none")
dev.off()

#"REFSEQ putative major capsid protein"  - too little nodes
rows_MCP_REFSEQ_putative  <- which(metadata$V5=="REFSEQ putative major capsid protein" & !is.na(metadata$Salinity_level))
tree_MCP_REFSEQ_putative <- keep.tip(tree, tip=metadata$genome[rows_MCP_REFSEQ_putative])
#only two nodes - too little to plot tree                                      
#freshwater only!
