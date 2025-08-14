library("phyloseq")
library("ggplot2") 
library(cowplot)

library("dplyr")
library(vegan)
library("tidyr")
library("ape")

setwd("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/mesocosm_warming")
count_tab <- read.table("feature-table-rare.tsv", header=T,
                        check.names=F, sep="\t", row.names = 1)
sample <- read.table("metadata.txt", header=T,
                     check.names=F, sep="\t", row.names = 1)
#change the data to get better legend 
# sample$`temperature elevation`[which(sample$`temperature elevation`=="0")] <- "tank delta 0"
# sample$`temperature elevation`[which(sample$`temperature elevation`=="2")] <- "tank delta 2"
# sample$`temperature elevation`[which(sample$`temperature elevation`=="4")] <- "tank delta 4"
sample$source <- sample$`temperature elevation`
sample$source[which(sample$`temperature elevation`=="-")] <- "field"
sample$source[which(sample$`temperature elevation`=="0")] <- "tank delta 0"
sample$source[which(sample$`temperature elevation`=="2")] <- "tank delta 2"
sample$source[which(sample$`temperature elevation`=="4")] <- "tank delta 4"
sample$`temperature elevation`[which(sample$`temperature elevation`=="-")] <- 0
sample$`temperature elevation` <- as.numeric(sample$`temperature elevation`)
colnames(sample) <- c("sample_nr", "sampling_date", "sampling_site", "fraction", "tank_or_field", 
                      "replicate", "type", "DNA_or_RNA", "DNA_RNA_label","`temperature elevation`", "source")
#export distance matrix from qime2
weighted_unifrac_distance_matrix <- read.table("weighted_unifrac_distance_matrix.tsv", header=T,
                                               check.names=F, sep="\t", row.names = 1)

# dist.unifrac(abund, tree, weighted = FALSE, normalise = TRUE, num_threads = 1L)

weighted_unifrac_dist <- as.dist(weighted_unifrac_distance_matrix) 
taxa <- row.names(count_tab)
otu <- count_tab#[,5:length(count_tab[1,])]
treefile <- read.tree("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/tree.nwk")

#divide into two experiments
#DNA
# otu_DNA <- otu %>% select(contains("D"))
otu_RNA <- otu %>% select(contains("R"))
samples_first <- row.names(sample)[which(sample$`sampling date` == "31/07/2024")]
#filter out abundance zero
sum_vec <- as.vector(apply(otu_RNA, 1, sum))

#sample_infected <- sample[1:28,1:4]

sample_RNA <- sample[which(sample$DNA_RNA_label == "R"),]

#pick only distances for RNA sample so PCA is not confounded
RNA_samples_names <- row.names(sample_RNA)
otu_RNA <- otu[, which(colnames(otu) %in% RNA_samples_names)]
RNAcol <- which(colnames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
RNArow <- which(rownames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
otu_RNA.nonzero <- otu_RNA[which(sum_vec!=0),]
otu_RNA.nonzero <- as.matrix(otu_RNA.nonzero)
taxa_RNA.nonzero <- taxa[which(sum_vec!=0)]
taxa_RNA.nonzero <- as.matrix(taxa_RNA.nonzero)
row.names(taxa_RNA.nonzero) <- taxa_RNA.nonzero[,1]



ps_RNA <- phyloseq(otu_table(otu_RNA.nonzero, taxa_are_rows=TRUE), 
                        sample_data(sample_RNA), 
                        tax_table(taxa_RNA.nonzero), treefile)

# ord.pcoa.bray <- ordinate(ps_RNA, method="PCoA", distance="bray", na.rm = T)
# plot_ordination(ps_RNA, ord.pcoa.bray, color="temperature_elevation")
# plot_ordination(ps_RNA, ord.pcoa.bray, color="fraction")
# plot_ordination(ps_RNA, ord.pcoa.bray, color="tank_or_field")
# ord.nmds.bray <- ordinate(ps_RNA, method="NMDS", distance="bray")
# plot_ordination(ps_RNA, ord.nmds.bray, color="temperature_elevation")
# plot_ordination(ps_RNA, ord.nmds.bray, color="fraction")
# plot_ordination(ps_RNA, ord.nmds.bray, color="tank_or_field")
# ord.cca.bray <- ordinate(ps_RNA, method="CCA", distance="bray")
# plot_ordination(ps_RNA, ord.cca.bray, color="temperature_elevation")
# plot_ordination(ps_RNA, ord.cca.bray, color="fraction")
# plot_ordination(ps_RNA, ord.cca.bray, color="tank_or_field")
#try unifrac
RNA_row_dist <- which(rownames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
RNA_col_dist <- which(colnames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
weighted_unifrac_dist_RNA <- as.dist(weighted_unifrac_distance_matrix[RNA_row_dist, RNA_col_dist]) 
# ord.pcoa.weighted_unifrac_dist <- ordinate(ps_RNA, method="PCoA", distance=weighted_unifrac_dist_RNA, na.rm = T)
# dist_matrix <- distance(ps_RNA, method = "unifrac", weighted = TRUE)
##### works!!!
#repeat for all combination RNA RNA FL PA
ps.ord.weighted_unifrac.phyoseq<- ordinate(ps_RNA, method="PCoA", distance="unifrac", weighted=TRUE)
plot_ordination(ps_RNA, ps.ord.weighted_unifrac.phyoseq, color="temperature_elevation")
plot_ordination(ps_RNA, ps.ord.weighted_unifrac.phyoseq, color="fraction")



ps_RNA_data<-merge(data.frame(ps.ord.weighted_unifrac.phyoseq[["vectors"]]), data.frame(sample_RNA), by=0)
ggplot(data=ps_RNA_data, aes(x=Axis.1, y=Axis.2, color=fraction, shape = temperature_elevation)) + geom_point(size= 3) + labs(x="PC1", y="PC2")
ps.ord.weighted_unifrac.phyoseq<- ordinate(ps_RNA, method="NMDS", distance="unifrac", weighted=TRUE)


tiff("RNA_temperature_elevation.tiff")
plot_ordination(ps_RNA, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
dev.off()
tiff("RNA_fraction.tiff")
plot_ordination(ps_RNA, ord.pcoa.weighted_unifrac_dist, color="fraction")
dev.off()
tiff("tank_or_field.tiff")
plot_ordination(ps_RNA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
dev.off()
require("ggrepel")
ggplot(data=p$data, aes(x=Axis.1, y=Axis.2, label=sample_nr, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + geom_text(hjust=0, vjust=0) + stat_ellipse(type="norm", level=0.97) + labs(x="PC1", y= "PC2")
#110 is clustering badly - remove as an outlier
tiff("RNA_fraction_temp_elevation.tiff")
RNA_samples_names <- RNA_samples_names[which(RNA_samples_names!="JG_110R")]
RNA_row_dist <- which(rownames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
RNA_col_dist <- which(colnames(weighted_unifrac_distance_matrix) %in% RNA_samples_names)
weighted_unifrac_dist_RNA <- as.dist(weighted_unifrac_distance_matrix[RNA_row_dist, RNA_col_dist]) 
p<-plot_ordination(ps_RNA, weighted_unifrac_dist_RNA, color="tank_or_field")
ggplot(data=p$data, aes(x=JG_100R, y=JG_101R, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y= "PC2")
dev.off()
#DNA
taxa_DNA.nonzero <- taxa[which(sum_vec!=0)]
taxa_DNA.nonzero <- as.matrix(taxa_DNA.nonzero)
sample_DNA <- sample[which(sample$DNA_RNA_label == "D"),]

DNA_samples_names <- row.names(sample_DNA)
otu_DNA <- otu[, which(colnames(otu) %in% DNA_samples_names)]
sum_vec <- as.vector(apply(otu_DNA, 1, sum))
otu_DNA.nonzero <- otu_DNA[which(sum_vec!=0),]
otu_DNA.nonzero <- as.matrix(otu_DNA.nonzero)
DNAcol <- which(colnames(weighted_unifrac_distance_matrix) %in% DNA_samples_names)
DNArow <- which(rownames(weighted_unifrac_distance_matrix) %in% DNA_samples_names)

weighted_unifrac_dist_DNA <- as.dist(weighted_unifrac_distance_matrix[DNArow, DNAcol]) 

row.names(taxa_DNA.nonzero) <- taxa_DNA.nonzero[,1]
ps_DNA <- phyloseq(otu_table(otu_DNA.nonzero, taxa_are_rows=TRUE), 
                   sample_data(sample_DNA), 
                   tax_table(taxa_DNA.nonzero))

#try unifrac
ord.pcoa.weighted_unifrac_dist <- ordinate(ps_DNA, method="PCoA", distance=weighted_unifrac_dist_DNA, na.rm = T)
tiff("DNA_temperature_elevation.tiff")
plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
dev.off()
tiff("DNA_fraction.tiff")
plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="fraction")
dev.off()
tiff("DNA_tank_or_field.tiff")
plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
dev.off()
plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="DNA_or_RNA")
p<-plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
ggplot(data=p$data, aes(x=JG_100D, y=JG_101D, label=sample_nr, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + geom_text(hjust=0, vjust=0) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y= "PC2")
#remove 110 - outlier 

p<-plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
ggplot(data=p$data, aes(x=JG_100D, y=JG_101D, label=sample_nr, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + geom_text(hjust=0, vjust=0) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y= "PC2")

tiff("DNA_fraction_temp_elevation.tiff")
DNA_samples_names <- DNA_samples_names[which(DNA_samples_names!="JG_110D")]
DNA_row_dist <- which(rownames(weighted_unifrac_distance_matrix) %in% DNA_samples_names)
DNA_col_dist <- which(colnames(weighted_unifrac_distance_matrix) %in% DNA_samples_names)
weighted_unifrac_dist_DNA <- as.dist(weighted_unifrac_distance_matrix[DNA_row_dist, DNA_col_dist]) 
p<-plot_ordination(ps_DNA, weighted_unifrac_dist_DNA, color="tank_or_field")
ggplot(data=p$data, aes(x=JG_100D, y=JG_101D, color=temperature_elevation, shape = fraction)) + geom_point(size= 3)+ stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y= "PC2")
dev.off()

tiff("DNA_fraction_temp_elevation.tiff")
p<-plot_ordination(ps_DNA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
ggplot(data=p$data, aes(x=Axis.1, y=Axis.2, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + stat_ellipse(type="norm", level=0.97) + labs(x="PC1", y= "PC2")
dev.off()
ggplot(data=p$data, aes(x=Axis.1, y=Axis.2, label=sample_nr, color=temperature_elevation, shape = fraction)) + geom_point(size= 3) + geom_text(hjust=0, vjust=0) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y= "PC2")

#check first phase of experiment
#31/07/2025
sample_first_day <- sample[which(sample$sampling_date == "31/07/2024"),]
sample_first_day_names <- row.names(sample_first_day)
otu_first_day <- otu[,which(colnames(otu) %in% sample_first_day_names)]
sum_vec <- as.vector(apply(otu_first_day, 1, sum))
otu_first_day.nonzero <- otu_first_day[which(sum_vec!=0),]
otu_first_day.nonzero <- as.matrix(otu_first_day.nonzero)
taxa_first_day.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day.nonzero <- as.matrix(taxa_first_day.nonzero)

first_day_samples_names <- row.names(sample_first_day)
first_day_samples_names <- sample_first_day[which(first_day_samples_names!="JG_110D" & first_day_samples_names!="JG_110R")]
first_daycol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_samples_names)
first_dayrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_samples_names)

weighted_unifrac_dist_first_day <- as.dist(weighted_unifrac_distance_matrix[first_dayrow, first_daycol]) 

row.names(taxa_first_day.nonzero) <- taxa_first_day.nonzero[,1]
ps_first_day <- phyloseq(otu_table(otu_first_day.nonzero, taxa_are_rows=TRUE), 
                   sample_data(sample_first_day), 
                   tax_table(taxa_first_day.nonzero))

#try unifrac
ord.pcoa.weighted_unifrac_dist <- ordinate(ps_first_day, method="PCoA", distance=weighted_unifrac_dist_first_day, na.rm = T)
plot_ordination(ps_first_day, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
plot_ordination(ps_first_day, ord.pcoa.weighted_unifrac_dist, color="fraction")
plot_ordination(ps_first_day, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
plot_ordination(ps_first_day, ord.pcoa.weighted_unifrac_dist, color="DNA_or_RNA")

#separated fractions
#FL
#first turn
sample_first_day_FL <- sample[which((sample$sampling_date == "29/07/2024" | sample$sampling_date == "31/07/2024") & sample$fraction == "FL"),]
sample_first_day_FL_names <- row.names(sample_first_day_FL)

#remove 110 paired samples 
sample_first_day_FL_names <- row.names(sample_first_day_FL)
sample_first_day_FL_names <- sample_first_day_FL_names[which(sample_first_day_FL_names!="JG_110D" & sample_first_day_FL_names!="JG_110R")]
otu_first_day_FL <- otu[,which(colnames(otu) %in% sample_first_day_FL_names)]
sum_vec <- as.vector(apply(otu_first_day_FL, 1, sum))
otu_first_day_FL.nonzero <- otu_first_day_FL[which(sum_vec!=0),]
otu_first_day_FL.nonzero <- as.matrix(otu_first_day_FL.nonzero)
taxa_first_day_FL.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_FL.nonzero <- as.matrix(taxa_first_day_FL.nonzero)

first_day_FL_samples_names <- row.names(sample_first_day_FL)
first_day_FLcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_FL_samples_names)
first_day_FLrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_FL_samples_names)

weighted_unifrac_dist_first_day_FL <- as.dist(weighted_unifrac_distance_matrix[first_day_FLrow, first_day_FLcol]) 

row.names(taxa_first_day_FL.nonzero) <- taxa_first_day_FL.nonzero[,1]
ps_first_day_FL <- phyloseq(otu_table(otu_first_day_FL.nonzero, taxa_are_rows=TRUE), 
                         sample_data(sample_first_day_FL), 
                         tax_table(taxa_first_day_FL.nonzero))

#try unifrac
ord.pcoa.weighted_unifrac_dist <- ordinate(ps_first_day_FL, method="PCoA", distance=weighted_unifrac_dist_first_day_FL, na.rm = T)
plot_ordination(ps_first_day_FL, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
plot_ordination(ps_first_day_FL, ord.pcoa.weighted_unifrac_dist, color="fraction")
plot_ordination(ps_first_day_FL, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
plot_ordination(ps_first_day_FL, ord.pcoa.weighted_unifrac_dist, color="DNA_or_RNA")
tiff("DNA_RNA_FL_fraction_temp_elevation.tiff")
p<-plot_ordination(ps_first_day_FL, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
ggplot(data=p$data, aes(x=Axis.1, y=Axis.2, color=temperature_elevation, shape = DNA_or_RNA)) + geom_point(size= 3) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y="PC2")
dev.off()
manovaFL <- manova(data=data.frame(p$data), formula = cbind(Axis.1, Axis.2) ~ temperature_elevation*DNA_or_RNA)
summary(manovaFL)        # Pillai's trace, Wilks' Lambda, etc.
summary.aov(manovaFL)    #
aov_pc1 <- aov(Axis.1 ~ temperature_elevation*DNA_or_RNA, data = p$data)
TukeyHSD(aov_pc1)
aov_pc2 <- aov(Axis.2 ~ temperature_elevation*DNA_or_RNA, data = p$data)
TukeyHSD(aov_pc2)
#PA
#first turn
#DNA
sample_first_day_PA_DNA_first_run <- sample[which((sample$sampling_date == "29/07/2024" | sample$sampling_date == "31/07/2024") & sample$fraction == "PA" & sample$DNA_or_RNA == "DNA"),]
sample_first_day_PA_DNA_first_run_names <- row.names(sample_first_day_PA_DNA_first_run)
otu_first_day_PA_DNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_PA_DNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_PA_DNA_first_run, 1, sum))
otu_first_day_PA_DNA_first_run.nonzero <- otu_first_day_PA_DNA_first_run[which(sum_vec!=0),]
otu_first_day_PA_DNA_first_run.nonzero <- as.matrix(otu_first_day_PA_DNA_first_run.nonzero)
taxa_first_day_PA_DNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_PA_DNA_first_run.nonzero <- as.matrix(taxa_first_day_PA_DNA_first_run.nonzero)

first_day_PA_DNA_first_run_samples_names <- row.names(sample_first_day_PA_DNA_first_run)
#remove 110
first_day_PA_DNA_first_run_samples_names <- first_day_PA_DNA_first_run_samples_names[which(first_day_PA_DNA_first_run_samples_names!="JG_110D" & first_day_PA_DNA_first_run_samples_names!="JG_110R")]
first_day_PA_DNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_PA_DNA_first_run_samples_names)
first_day_PA_DNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_PA_DNA_first_run_samples_names)

weighted_unifrac_dist_first_day_PA_DNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_PA_DNA_first_runrow, first_day_PA_DNA_first_runcol]) 

row.names(taxa_first_day_PA_DNA_first_run.nonzero) <- taxa_first_day_PA_DNA_first_run.nonzero[,1]
ps_first_day_PA_DNA_first_run <- phyloseq(otu_table(otu_first_day_PA_DNA_first_run.nonzero, taxa_are_rows=TRUE), 
                            sample_data(sample_first_day_PA_DNA_first_run), 
                            tax_table(taxa_first_day_PA_DNA_first_run.nonzero), treefile)

ord.pcoa.unweighted_unifrac_dist.first_day_PA_DNA_first_run <- ordinate(ps_first_day_PA_DNA_first_run, method="PCoA", distance="unifrac", weighted=FALSE)
DNA_PA_plot <- plot_ordination(ps_first_day_PA_DNA_first_run, ord.pcoa.unweighted_unifrac_dist.first_day_PA_DNA_first_run, color="source") + theme(legend.position ="none") +  stat_ellipse(level = 0.95)
DNA_PA_plot_leg <- plot_ordination(ps_first_day_PA_DNA_first_run, ord.pcoa.unweighted_unifrac_dist.first_day_PA_DNA_first_run, color="source")
ord.pcoa.weighted_unifrac_dist.first_day_PA_DNA_first_run <- ordinate(ps_first_day_PA_DNA_first_run, method="PCoA", distance=weighted_unifrac_dist_first_day_PA_DNA_first_run)
DNA_PA_plot_weighted <- plot_ordination(ps_first_day_PA_DNA_first_run, ord.pcoa.weighted_unifrac_dist.first_day_PA_DNA_first_run, color="source")  + theme(legend.position ="none") + stat_ellipse(level = 0.95)


#FL
#first turn
#DNA
sample_first_day_FL_DNA_first_run <- sample[which((sample$sampling_date == "29/07/2024" | sample$sampling_date == "31/07/2024") & sample$fraction == "FL" & sample$DNA_or_RNA == "DNA"),]
sample_first_day_FL_DNA_first_run_names <- row.names(sample_first_day_FL_DNA_first_run)
otu_first_day_FL_DNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_FL_DNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_FL_DNA_first_run, 1, sum))
otu_first_day_FL_DNA_first_run.nonzero <- otu_first_day_FL_DNA_first_run[which(sum_vec!=0),]
otu_first_day_FL_DNA_first_run.nonzero <- as.matrix(otu_first_day_FL_DNA_first_run.nonzero)
taxa_first_day_FL_DNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_FL_DNA_first_run.nonzero <- as.matrix(taxa_first_day_FL_DNA_first_run.nonzero)

first_day_FL_DNA_first_run_samples_names <- row.names(sample_first_day_FL_DNA_first_run)
#remove 110
first_day_FL_DNA_first_run_samples_names <- first_day_FL_DNA_first_run_samples_names[which(first_day_FL_DNA_first_run_samples_names!="JG_110D" & first_day_FL_DNA_first_run_samples_names!="JG_110R")]
first_day_FL_DNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_FL_DNA_first_run_samples_names)
first_day_FL_DNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_FL_DNA_first_run_samples_names)

weighted_unifrac_dist_first_day_FL_DNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_FL_DNA_first_runrow, first_day_FL_DNA_first_runcol]) 

row.names(taxa_first_day_FL_DNA_first_run.nonzero) <- taxa_first_day_FL_DNA_first_run.nonzero[,1]
ps_first_day_FL_DNA_first_run <- phyloseq(otu_table(otu_first_day_FL_DNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_FL_DNA_first_run), 
                                          tax_table(taxa_first_day_FL_DNA_first_run.nonzero), treefile)

ord.pcoa.unweighted_unifrac_dist.first_day_FL_DNA_first_run <- ordinate(ps_first_day_FL_DNA_first_run, method="PCoA", distance="unifrac", weighted=FALSE)
DNA_FL_plot <- plot_ordination(ps_first_day_FL_DNA_first_run, ord.pcoa.unweighted_unifrac_dist.first_day_FL_DNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)

ord.pcoa.weighted_unifrac_dist.first_day_FL_DNA_first_run <- ordinate(ps_first_day_FL_DNA_first_run, method="PCoA", distance=weighted_unifrac_dist_first_day_FL_DNA_first_run)
DNA_FL_plot_weighted <- plot_ordination(ps_first_day_FL_DNA_first_run, ord.pcoa.weighted_unifrac_dist.first_day_FL_DNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)


#PA
#first turn
#RNA
sample_first_day_PA_RNA_first_run <- sample[which((sample$sampling_date == "29/07/2024" | sample$sampling_date == "31/07/2024") & sample$fraction == "PA" & sample$DNA_or_RNA == "RNA"),]
sample_first_day_PA_RNA_first_run_names <- row.names(sample_first_day_PA_RNA_first_run)
otu_first_day_PA_RNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_PA_RNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_PA_RNA_first_run, 1, sum))
otu_first_day_PA_RNA_first_run.nonzero <- otu_first_day_PA_RNA_first_run[which(sum_vec!=0),]
otu_first_day_PA_RNA_first_run.nonzero <- as.matrix(otu_first_day_PA_RNA_first_run.nonzero)
taxa_first_day_PA_RNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_PA_RNA_first_run.nonzero <- as.matrix(taxa_first_day_PA_RNA_first_run.nonzero)

first_day_PA_RNA_first_run_samples_names <- row.names(sample_first_day_PA_RNA_first_run)
#remove 110
first_day_PA_RNA_first_run_samples_names <- first_day_PA_RNA_first_run_samples_names[which(first_day_PA_RNA_first_run_samples_names!="JG_110D" & first_day_PA_RNA_first_run_samples_names!="JG_110R")]
first_day_PA_RNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_PA_RNA_first_run_samples_names)
first_day_PA_RNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_PA_RNA_first_run_samples_names)

weighted_unifrac_dist_first_day_PA_RNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_PA_RNA_first_runrow, first_day_PA_RNA_first_runcol]) 

row.names(taxa_first_day_PA_RNA_first_run.nonzero) <- taxa_first_day_PA_RNA_first_run.nonzero[,1]
ps_first_day_PA_RNA_first_run <- phyloseq(otu_table(otu_first_day_PA_RNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_PA_RNA_first_run), 
                                          tax_table(taxa_first_day_PA_RNA_first_run.nonzero), treefile)

ord.pcoa.unweighted_unifrac_dist.first_day_PA_RNA_first_run <- ordinate(ps_first_day_PA_RNA_first_run, method="PCoA", distance="unifrac", weighted=FALSE)
RNA_PA_plot <- plot_ordination(ps_first_day_PA_RNA_first_run, ord.pcoa.unweighted_unifrac_dist.first_day_PA_RNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)
ord.pcoa.weighted_unifrac_dist.first_day_PA_RNA_first_run <- ordinate(ps_first_day_PA_RNA_first_run, method="PCoA", distance=weighted_unifrac_dist_first_day_PA_RNA_first_run)
RNA_PA_plot_weighted <- plot_ordination(ps_first_day_PA_RNA_first_run, ord.pcoa.weighted_unifrac_dist.first_day_PA_RNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)


#FL
#first turn
#RNA
sample_first_day_FL_RNA_first_run <- sample[which((sample$sampling_date == "29/07/2024" | sample$sampling_date == "31/07/2024") & sample$fraction == "FL" & sample$DNA_or_RNA == "RNA"),]
sample_first_day_FL_RNA_first_run_names <- row.names(sample_first_day_FL_RNA_first_run)
otu_first_day_FL_RNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_FL_RNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_FL_RNA_first_run, 1, sum))
otu_first_day_FL_RNA_first_run.nonzero <- otu_first_day_FL_RNA_first_run[which(sum_vec!=0),]
otu_first_day_FL_RNA_first_run.nonzero <- as.matrix(otu_first_day_FL_RNA_first_run.nonzero)
taxa_first_day_FL_RNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_FL_RNA_first_run.nonzero <- as.matrix(taxa_first_day_FL_RNA_first_run.nonzero)

first_day_FL_RNA_first_run_samples_names <- row.names(sample_first_day_FL_RNA_first_run)
#remove 110
first_day_FL_RNA_first_run_samples_names <- first_day_FL_RNA_first_run_samples_names[which(first_day_FL_RNA_first_run_samples_names!="JG_110D" & first_day_FL_RNA_first_run_samples_names!="JG_110R")]
first_day_FL_RNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_FL_RNA_first_run_samples_names)
first_day_FL_RNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_FL_RNA_first_run_samples_names)

weighted_unifrac_dist_first_day_FL_RNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_FL_RNA_first_runrow, first_day_FL_RNA_first_runcol]) 

row.names(taxa_first_day_FL_RNA_first_run.nonzero) <- taxa_first_day_FL_RNA_first_run.nonzero[,1]
ps_first_day_FL_RNA_first_run <- phyloseq(otu_table(otu_first_day_FL_RNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_FL_RNA_first_run), 
                                          tax_table(taxa_first_day_FL_RNA_first_run.nonzero), treefile)

ord.pcoa.unweighted_unifrac_dist.first_day_FL_RNA_first_run <- ordinate(ps_first_day_FL_RNA_first_run, method="PCoA", distance="unifrac", weighted=FALSE)
RNA_FL_plot <- plot_ordination(ps_first_day_FL_RNA_first_run, ord.pcoa.unweighted_unifrac_dist.first_day_FL_RNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)
ord.pcoa.weighted_unifrac_dist.first_day_FL_RNA_first_run <- ordinate(ps_first_day_FL_RNA_first_run, method="PCoA", distance=weighted_unifrac_dist_first_day_FL_RNA_first_run) 
RNA_FL_plot_weighted <- plot_ordination(ps_first_day_FL_RNA_first_run, ord.pcoa.weighted_unifrac_dist.first_day_FL_RNA_first_run, color="source") + theme(legend.position ="none") + stat_ellipse(level = 0.95)

#plot all together 
title_gg_unweighted <- ggdraw() + draw_label("Unweighted UniFrac", fontface='bold', hjust=0.5)
title_gg_weighted <- ggdraw() + draw_label("Weighted UniFrac", fontface='bold', hjust=0.5)
#unweighted
unweigted_no_legend <- plot_grid(RNA_PA_plot + theme(plot.margin = unit(c(2,2,2,2), "lines")), DNA_PA_plot + theme(plot.margin = unit(c(2,2,2,2), "lines"))  + theme(legend.position ="none"), RNA_FL_plot + theme(plot.margin = unit(c(2,2,2,2), "lines")), DNA_FL_plot + theme(plot.margin = unit(c(2,2,2,2), "lines")), labels = c('A) RNA PA', 'B) DNA PA', 'C) RNA FL', 'D) DNA FL'), label_size = 12) + ggtitle("Unweigted unifrac")
legend <- get_legend(
  # create some space to the left of the legend
  DNA_PA_plot_leg
)

tiff("unweighted_DNA_RNA.tiff")
unweighted_legend <- plot_grid(unweigted_no_legend, legend, rel_widths = c(4, 1), ncol=2)
unweighted_legend <- plot_grid(title_gg_unweighted, unweighted_legend, ncol=1, rel_heights=c(0.1, 1), align = "h")
unweighted_legend
dev.off()
#weighted
tiff("weighted_DNA_RNA.tiff")
weigthed_no_legend <- plot_grid(RNA_PA_plot_weighted + theme(plot.margin = unit(c(2,2,2,2), "lines")), DNA_PA_plot_weighted + theme(plot.margin = unit(c(2,2,2,2), "lines")), RNA_FL_plot_weighted + theme(plot.margin = unit(c(2,2,2,2), "lines")), DNA_FL_plot_weighted + theme(plot.margin = unit(c(2,2,2,2), "lines")), labels = c('A) RNA PA', 'B) DNA PA', 'C) RNA FL', 'D) DNA FL'), label_size = 12)
weigthed_legend <- plot_grid(weigthed_no_legend, legend, ncol=2, rel_widths = c(4, 1))
weigthed_legend <- plot_grid(title_gg_weighted, weigthed_legend, ncol=1, rel_heights=c(0.1, 1), aligh = "h")
weigthed_legend 
dev.off()

#mantel test for weigthed, unweighthed RNA, DNA, FL, PA, but with out field samples - prepare all matrixes
#DNA PA
sample_first_day_PA_DNA_first_run <- sample[which( sample$sampling_date == "31/07/2024" & sample$fraction == "PA" & sample$DNA_or_RNA == "DNA"),]
sample_first_day_PA_DNA_first_run_names <- row.names(sample_first_day_PA_DNA_first_run)
otu_first_day_PA_DNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_PA_DNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_PA_DNA_first_run, 1, sum))
otu_first_day_PA_DNA_first_run.nonzero <- otu_first_day_PA_DNA_first_run[which(sum_vec!=0),]
otu_first_day_PA_DNA_first_run.nonzero <- as.matrix(otu_first_day_PA_DNA_first_run.nonzero)
taxa_first_day_PA_DNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_PA_DNA_first_run.nonzero <- as.matrix(taxa_first_day_PA_DNA_first_run.nonzero)

first_day_PA_DNA_first_run_samples_names <- row.names(sample_first_day_PA_DNA_first_run)
#remove 110
first_day_PA_DNA_first_run_samples_names <- first_day_PA_DNA_first_run_samples_names[which(first_day_PA_DNA_first_run_samples_names!="JG_110D" & first_day_PA_DNA_first_run_samples_names!="JG_110R")]
first_day_PA_DNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_PA_DNA_first_run_samples_names)
first_day_PA_DNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_PA_DNA_first_run_samples_names)

weighted_unifrac_dist_first_day_PA_DNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_PA_DNA_first_runrow, first_day_PA_DNA_first_runcol]) 

row.names(taxa_first_day_PA_DNA_first_run.nonzero) <- taxa_first_day_PA_DNA_first_run.nonzero[,1]
ps_first_day_PA_DNA_first_run <- phyloseq(otu_table(otu_first_day_PA_DNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_PA_DNA_first_run), 
                                          tax_table(taxa_first_day_PA_DNA_first_run.nonzero), treefile)
unweighted_unifrac_dist_first_day_PA_DNA_first_run <-distance(ps_first_day_PA_DNA_first_run, "unifrac", type = "samples")
#DNA FL
sample_first_day_FL_DNA_first_run <- sample[which( sample$sampling_date == "31/07/2024" & sample$fraction == "FL" & sample$DNA_or_RNA == "DNA" & row.names(sample) != "JG_110D"),] #check again
sample_first_day_FL_DNA_first_run_names <- row.names(sample_first_day_FL_DNA_first_run)
otu_first_day_FL_DNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_FL_DNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_FL_DNA_first_run, 1, sum))
otu_first_day_FL_DNA_first_run.nonzero <- otu_first_day_FL_DNA_first_run[which(sum_vec!=0),]
otu_first_day_FL_DNA_first_run.nonzero <- as.matrix(otu_first_day_FL_DNA_first_run.nonzero)
taxa_first_day_FL_DNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_FL_DNA_first_run.nonzero <- as.matrix(taxa_first_day_FL_DNA_first_run.nonzero)

first_day_FL_DNA_first_run_samples_names <- row.names(sample_first_day_FL_DNA_first_run)
#remove 110
first_day_FL_DNA_first_run_samples_names <- first_day_FL_DNA_first_run_samples_names[which(first_day_FL_DNA_first_run_samples_names!="JG_110D" & first_day_FL_DNA_first_run_samples_names!="JG_110R")]
first_day_FL_DNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_FL_DNA_first_run_samples_names)
first_day_FL_DNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_FL_DNA_first_run_samples_names)

weighted_unifrac_dist_first_day_FL_DNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_FL_DNA_first_runrow, first_day_FL_DNA_first_runcol]) 

row.names(taxa_first_day_FL_DNA_first_run.nonzero) <- taxa_first_day_FL_DNA_first_run.nonzero[,1]
ps_first_day_FL_DNA_first_run <- phyloseq(otu_table(otu_first_day_FL_DNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_FL_DNA_first_run), 
                                          tax_table(taxa_first_day_FL_DNA_first_run.nonzero), treefile)
unweighted_unifrac_dist_first_day_FL_DNA_first_run <-distance(ps_first_day_FL_DNA_first_run, "unifrac", type = "samples")

#RNA PA
sample_first_day_PA_RNA_first_run <- sample[which( sample$sampling_date == "31/07/2024" & sample$fraction == "PA" & sample$DNA_or_RNA == "RNA"),]
sample_first_day_PA_RNA_first_run_names <- row.names(sample_first_day_PA_RNA_first_run)
otu_first_day_PA_RNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_PA_RNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_PA_RNA_first_run, 1, sum))
otu_first_day_PA_RNA_first_run.nonzero <- otu_first_day_PA_RNA_first_run[which(sum_vec!=0),]
otu_first_day_PA_RNA_first_run.nonzero <- as.matrix(otu_first_day_PA_RNA_first_run.nonzero)
taxa_first_day_PA_RNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_PA_RNA_first_run.nonzero <- as.matrix(taxa_first_day_PA_RNA_first_run.nonzero)

first_day_PA_RNA_first_run_samples_names <- row.names(sample_first_day_PA_RNA_first_run)
#remove 110
first_day_PA_RNA_first_run_samples_names <- first_day_PA_RNA_first_run_samples_names[which(first_day_PA_RNA_first_run_samples_names!="JG_110D" & first_day_PA_RNA_first_run_samples_names!="JG_110R")]
first_day_PA_RNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_PA_RNA_first_run_samples_names)
first_day_PA_RNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_PA_RNA_first_run_samples_names)

weighted_unifrac_dist_first_day_PA_RNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_PA_RNA_first_runrow, first_day_PA_RNA_first_runcol]) 

row.names(taxa_first_day_PA_RNA_first_run.nonzero) <- taxa_first_day_PA_RNA_first_run.nonzero[,1]
ps_first_day_PA_RNA_first_run <- phyloseq(otu_table(otu_first_day_PA_RNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_PA_RNA_first_run), 
                                          tax_table(taxa_first_day_PA_RNA_first_run.nonzero), treefile)
unweighted_unifrac_dist_first_day_PA_RNA_first_run <-distance(ps_first_day_PA_RNA_first_run, "unifrac", type = "samples")
#RNA FL
#remove 110
first_day_FL_RNA_first_run_samples_names <- first_day_FL_RNA_first_run_samples_names[which(first_day_FL_RNA_first_run_samples_names!="JG_110D" & first_day_FL_RNA_first_run_samples_names!="JG_110R")]
first_day_FL_RNA_first_runcol <- which(colnames(weighted_unifrac_distance_matrix) %in% first_day_FL_RNA_first_run_samples_names)
first_day_FL_RNA_first_runrow <- which(rownames(weighted_unifrac_distance_matrix) %in% first_day_FL_RNA_first_run_samples_names)

weighted_unifrac_dist_first_day_FL_RNA_first_run <- as.dist(weighted_unifrac_distance_matrix[first_day_FL_RNA_first_runrow, first_day_FL_RNA_first_runcol]) 
sample_first_day_FL_RNA_first_run <- sample[which( sample$sampling_date == "31/07/2024" & sample$fraction == "FL" & sample$DNA_or_RNA == "RNA" & row.names(sample) != "JG_110R"),]
sample_first_day_FL_RNA_first_run_names <- row.names(sample_first_day_FL_RNA_first_run)
otu_first_day_FL_RNA_first_run <- otu[,which(colnames(otu) %in% sample_first_day_FL_RNA_first_run_names)]
sum_vec <- as.vector(apply(otu_first_day_FL_RNA_first_run, 1, sum))
otu_first_day_FL_RNA_first_run.nonzero <- otu_first_day_FL_RNA_first_run[which(sum_vec!=0),]
otu_first_day_FL_RNA_first_run.nonzero <- as.matrix(otu_first_day_FL_RNA_first_run.nonzero)
taxa_first_day_FL_RNA_first_run.nonzero <- taxa[which(sum_vec!=0)]
taxa_first_day_FL_RNA_first_run.nonzero <- as.matrix(taxa_first_day_FL_RNA_first_run.nonzero)

first_day_FL_RNA_first_run_samples_names <- row.names(sample_first_day_FL_RNA_first_run)

row.names(taxa_first_day_FL_RNA_first_run.nonzero) <- taxa_first_day_FL_RNA_first_run.nonzero[,1]
ps_first_day_FL_RNA_first_run <- phyloseq(otu_table(otu_first_day_FL_RNA_first_run.nonzero, taxa_are_rows=TRUE), 
                                          sample_data(sample_first_day_FL_RNA_first_run), 
                                          tax_table(taxa_first_day_FL_RNA_first_run.nonzero), treefile)
unweighted_unifrac_dist_first_day_FL_RNA_first_run <-distance(ps_first_day_FL_RNA_first_run, "unifrac", type = "samples")
#weighted dist
weighted_unifrac_dist_first_day_PA_DNA_first_run
weighted_unifrac_dist_first_day_FL_DNA_first_run
weighted_unifrac_dist_first_day_PA_RNA_first_run
weighted_unifrac_dist_first_day_FL_RNA_first_run
#unweighted dist 
unweighted_unifrac_dist_first_day_PA_DNA_first_run
unweighted_unifrac_dist_first_day_FL_DNA_first_run
unweighted_unifrac_dist_first_day_PA_RNA_first_run
unweighted_unifrac_dist_first_day_FL_RNA_first_run
# dist of temp 
#dist of temp 
dist_temp_DNA_PA <- vegdist(sample[which(row.names(sample) %in% colnames(as.matrix(unweighted_unifrac_dist_first_day_PA_DNA_first_run))), 10]+20,"bray", na.rm = T)
dist_temp_DNA_FL <- vegdist(sample[which(row.names(sample) %in% colnames(as.matrix(unweighted_unifrac_dist_first_day_FL_DNA_first_run))), 10]+20,"bray", na.rm = T)
dist_temp_RNA_PA <- vegdist(sample[which(row.names(sample) %in% colnames(as.matrix(unweighted_unifrac_dist_first_day_PA_RNA_first_run))), 10]+20,"bray", na.rm = T)
dist_temp_RNA_FL <- vegdist(sample[which(row.names(sample) %in% colnames(as.matrix(unweighted_unifrac_dist_first_day_FL_RNA_first_run))), 10]+20,"bray", na.rm = T)
# rna and dna are coming from the same tanks so matrixes are the same  (checked) there can be differences between FL and PA samples (e.g some didnt work in sequencing  - outlier removal)
dist_temp_DNA_PA == dist_temp_DNA_FL
dist_temp_DNA_PA == dist_temp_RNA_PA
dist_temp_DNA_FL == dist_temp_RNA_FL
list.PA.dists <- list(weighted_unifrac_dist_first_day_PA_DNA_first_run, weighted_unifrac_dist_first_day_PA_RNA_first_run, 
                        unweighted_unifrac_dist_first_day_PA_DNA_first_run, unweighted_unifrac_dist_first_day_PA_RNA_first_run)

list.FL.dists <- list(weighted_unifrac_dist_first_day_FL_DNA_first_run, weighted_unifrac_dist_first_day_FL_RNA_first_run, 
                        unweighted_unifrac_dist_first_day_FL_DNA_first_run, unweighted_unifrac_dist_first_day_FL_RNA_first_run)
res_mantel_sig <- matrix(nrow=2, ncol=4, dimnames = list(c("PA", "FL"), c("weighted Unifrac DNA", "weighted Unifrac RNA" , "unweighted Unifrac DNA", "unweighted Unifrac RNA" )))
res_mantel_R <- matrix(nrow=2, ncol=4, dimnames = list(c("PA", "FL"), c("weighted Unifrac DNA", "weighted Unifrac RNA" , "unweighted Unifrac DNA", "unweighted Unifrac RNA" )))

for (i in c(1:length(list.PA.dists))){
dist.matrix <- list.PA.dists[[i]]
res_PA <- mantel(dist.matrix, dist_temp_DNA_PA, method = "spearman", permutations = 9999, na.rm = TRUE)
res_mantel_sig[1,i] <- res_PA$signif
res_mantel_R[1,i] <- res_PA$statistic
dist.matrix <- list.FL.dists[[i]]
res_FL <- mantel(dist.matrix, dist_temp_DNA_PA, method = "spearman", permutations = 9999, na.rm = TRUE)
res_mantel_sig[2,i] <- res_FL$signif
res_mantel_R[2,i] <- res_FL$statistic
}#sth wrong with dimentions unweighted FL ne more samples 


#try unifrac
# ord.pcoa.weighted_unifrac_dist <- ordinate(ps_first_day_PA, method="PCoA", distance=weighted_unifrac_dist_first_day_PA, na.rm = T)
# plot_ordination(ps_first_day_PA, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
# plot_ordination(ps_first_day_PA, ord.pcoa.weighted_unifrac_dist, color="fraction")
# plot_ordination(ps_first_day_PA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
# plot_ordination(ps_first_day_PA, ord.pcoa.weighted_unifrac_dist, color="DNA_or_RNA")
# tiff("DNA_RNA_PA_fraction_temp_elevation.tiff")
# p<-plot_ordination(ps_first_day_PA, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
# ggplot(data=p$data, aes(x=Axis.1, y=Axis.2, color=temperature_elevation, shape = DNA_or_RNA)) + geom_point(size= 3) + stat_ellipse(type="norm", level=0.95) + labs(x="PC1", y="PC2")
# dev.off()
#PA

manovaPA <- manova(data=data.frame(p$data), formula = cbind(Axis.1, Axis.2) ~ temperature_elevation*DNA_or_RNA)
summary(manovaPA)        # Pillai's trace, Wilks' Lambda, etc.
summary.aov(manovaPA)    #
aov_pc1 <- aov(Axis.1 ~ temperature_elevation*DNA_or_RNA, data = p$data)
TukeyHSD(aov_pc1)
aov_pc2 <- aov(Axis.2 ~ temperature_elevation*DNA_or_RNA, data = p$data)
TukeyHSD(aov_pc2)
#sediment with a second batch
sample_sed <- sample[which(sample$sampling_date == "21/08/2024" | sample$sampling_date == "23/08/2024"),]
sample_sed_names <- row.names(sample_sed)
otu_sed <- otu[,which(colnames(otu) %in% sample_sed_names)]
sum_vec <- as.vector(apply(otu_sed, 1, sum))
otu_sed.nonzero <- otu_sed[which(sum_vec!=0),]
otu_sed.nonzero <- as.matrix(otu_sed.nonzero)
taxa_sed.nonzero <- taxa[which(sum_vec!=0)]
taxa_sed.nonzero <- as.matrix(taxa_sed.nonzero)

sed_samples_names <- row.names(sample_sed)
sedcol <- which(colnames(weighted_unifrac_distance_matrix) %in% sed_samples_names)
sedrow <- which(rownames(weighted_unifrac_distance_matrix) %in% sed_samples_names)

weighted_unifrac_dist_sed <- as.dist(weighted_unifrac_distance_matrix[sedrow, sedcol]) 

row.names(taxa_sed.nonzero) <- taxa_sed.nonzero[,1]
ps_sed <- phyloseq(otu_table(otu_sed.nonzero, taxa_are_rows=TRUE), 
                            sample_data(sample_sed), 
                            tax_table(taxa_sed.nonzero))

#try unifrac
ord.pcoa.weighted_unifrac_dist <- ordinate(ps_sed, method="PCoA", distance=weighted_unifrac_dist_sed, na.rm = T)
plot_ordination(ps_sed, ord.pcoa.weighted_unifrac_dist, color="temperature_elevation")
plot_ordination(ps_sed, ord.pcoa.weighted_unifrac_dist, color="fraction")
plot_ordination(ps_sed, ord.pcoa.weighted_unifrac_dist, color="tank_or_field")
plot_ordination(ps_sed, ord.pcoa.weighted_unifrac_dist, color="DNA_or_RNA")
plot_ordination(ps_sed, ord.pcoa.weighted_unifrac_dist, color="replicate")

library(RColorBrewer)
colours <- rep(brewer.pal(n = 8, name = "Dark2"),10 )
taxa_level7 <- read.table("level-7.csv", header=T,
                        check.names=F, sep=",", row.names = 1)
taxa_level7_long <- taxa_level7 %>%
  pivot_longer(!colnames(taxa_level7)[43:52], names_to = "taxa", values_to = "count") %>% filter(count>0)
taxa_level7_long_DNA <- taxa_level7_long[which(taxa_level7_long$`DNA or RNA`=="DNA" & taxa_level7_long$type=="Water"), ]
taxa_level7_long_DNA <-  taxa_level7_long_DNA %>% mutate(alias= paste( fraction, `temperature elevation`, `sample nr`,sep = '_'))
DNA_taxa<- ggplot(data=taxa_level7_long_DNA, aes(x=alias, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#repeat for level 3
taxa_level3 <- read.table("level-3.csv", header=T,
                          check.names=F, sep=",", row.names = 1)
taxa_level3_long <- taxa_level3 %>%
  pivot_longer(!colnames(taxa_level3)[17:26], names_to = "taxa", values_to = "count") %>% filter(count>0)
taxa_level3_long_DNA <- taxa_level3_long[which(taxa_level3_long$`DNA or RNA`=="DNA" & taxa_level3_long$type=="Water"), ]
taxa_level3_long_DNA <-  taxa_level3_long_DNA %>% mutate(alias= paste( fraction, `temperature elevation`, `sample nr`, `sampling date`), sep = '_')
DNA_taxa<- ggplot(data=taxa_level3_long_DNA, aes(x=alias, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#comapre with RNA
taxa_level3_long_RNA <- taxa_level3_long[which(taxa_level3_long$`DNA or RNA`=="RNA" & taxa_level3_long$type=="Water"), ]
taxa_level3_long_RNA <-  taxa_level3_long_RNA %>% mutate(alias= paste( fraction, `temperature elevation`, `sample nr`, `sampling date`), sep = '_')
RNA_taxa<- ggplot(data=taxa_level3_long_RNA, aes(x=alias, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#facet RNA and dna
#comapre with RNA
taxa_level3_long <- taxa_level3_long[which(taxa_level3_long$type=="Water" & (taxa_level3_long$`sampling date`=="29/07/2024" | taxa_level3_long$`sampling date`=="31/07/2024" )),]
taxa_level3_long <-  taxa_level3_long %>% mutate(alias= paste( fraction, `temperature elevation`, `sample nr`, `tank or field`), sep = '_')
taxa <- ggplot(data=taxa_level3_long, aes(x=alias, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(rows = vars(`DNA or RNA`))
tiff("taxa.tiff")
taxa
dev.off()

taxa_level3_long <- taxa_level3 %>%
  pivot_longer(!colnames(taxa_level3)[17:26], names_to = "taxa", values_to = "count") %>% filter(count>0)
taxa_level3_long_sed <- taxa_level3_long[which(taxa_level3_long$`type`=="Sediment"),]
taxa_level3_long_sed <-  taxa_level3_long_sed %>% mutate(alias= paste( fraction, `temperature elevation`, `sample nr`), sep = '_')
taxa_sed <- ggplot(data=taxa_level3_long_sed, aes(x=alias, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(rows = vars(`DNA or RNA`))
tiff("taxa_sed.tiff")
taxa_sed
dev.off()

taxa_level3_long_sed_vs_water <- taxa_level3_long[which(taxa_level3_long$`sampling date`=="23/08/2024" | taxa_level3_long$`sampling date`=="21/08/2024" ),]
taxa_level3_long_sed_vs_water <-  taxa_level3_long_sed_vs_water %>% mutate(alias= paste(fraction, `temperature elevation`, `sample nr`), sep = '_')
#repair mistake in metadata
taxa_level3_long_sed_vs_water$sampling_site[which(taxa_level3_long_sed_vs_water$sampling_site=="1B")] <-"1C"
taxa_level3_long_sed_vs_water$sampling_site[which(taxa_level3_long_sed_vs_water$sampling_site=="3J")] <-"3I"

taxa_sed_vs_water <- ggplot(data=taxa_level3_long_sed_vs_water, aes(x=sampling_site, y=count, fill=taxa)) +
  geom_bar(position="fill", stat="identity") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = colours) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(rows = vars(`type`))
tiff("taxa_sed_vs_water.tiff")
taxa_sed_vs_water
dev.off()


