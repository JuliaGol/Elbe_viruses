#setup
library(ggplot2)
library(FactoMineR)
library(factoextra)
#library(missMDA)
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/PCA_outputs")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog"
metadata <- read.csv(paste0(path,"/metadata_all.csv"))
#read PCA results calculated on the server 
#PCA_expression <- read.csv("PCA_expression.csv")
#PCA_vir_abund <- read.csv("PCA_vir_abund.csv")
#PCA_vir_transcription <- read.csv("PCA_vir_transcription.csv")
#set first column as sample ids name
#colnames(PCA_expression)[1] <- "sampleid"
#colnames(PCA_vir_transcription)[1] <- "sampleid"
#colnames(PCA_vir_abund)[1] <- "sampleid"




#PCA
PCA_log_expression <- PCA(t_log_paired_vir_metat_paired, graph=F)
PCA_log_transcript <- PCA(t_log_paired_vir_metat_paired, graph=F)
PCA_log_abund <- PCA(t_log_paired_vir_metag_paired, graph=F)

saveRDS(PCA_log_expression, file="PCA_log_expression.RDS")
saveRDS(PCA_log_transcript, file="PCA_log_transcript.RDS")
saveRDS(PCA_log_abund, file="PCA_log_abund.RDS")

#find what are two outliers 

#take only numeric data 
metadata_num <-metadata[,c("Sat_O2_TBDHereon", "Temperature_TBDHereon", "Salinity_TBDHereon", "Turbidity_TBDHereon", "pH_TBDHereon", "O2_TBDHereon", "SPM_mgperL", "DOC_mg.L", "TN_mg.L",
                           "DIC_mg.L", "Silicate_mg.L", "Ammonium_mg.L", "Nitrate_mg.L", "Total_DIN_µM", "TotalDissolvedPhosphate_mg.L", "Ammonium_µM",
                           "Nitrite_µM", "Nitrate_µM", "SRP_µM", "Phosphate_µM", "RespirationRate_O2ug.L.h", "POC_mgperL", "PTC_mgperL", "PTN_mgperL",
                           "PTH_mgperL", "dCH4_nM", "Chlorophyll", "dCO2_uM")]
#and ensure they have point instead of comma 
#and there are set as numeric
metadata_num <- apply(apply(metadata_num, 2, gsub, patt=",", replace="."), 2, as.numeric)
metadata_num <- cbind(metadata_num, metadata$sampleid) 
colnames(metadata_num)[length(colnames(metadata_num))] <- "sampleid"   
t_norm_expression_paired_vir <- merge(t_norm_expression_paired_vir, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames

PCA_log_transcript <- readRDS(file="PCA_log_transcript.RDS")
PCA_log_abund <- readRDS(file="PCA_log_abund.RDS")
t_log_paired_vir_metag_paired <- merge(t_log_paired_vir_metag_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames
t_log_paired_vir_metat_paired <- merge(t_log_paired_vir_metat_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames
#visualise log abundance 
png(filename="PCA_log_transcript=Sample_type_Station.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_type, shape=Station)) + geom_point() + ggtitle("Viral log-transcripon")
dev.off()
png(filename="PCA_log_transcript_Sample_date_Station.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-transcripon")
dev.off()
png(filename="PCA_log_transcript_Station_Sample_date.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_date)) + geom_point() + ggtitle("Viral log-transcripon")
dev.off()
png(filename="PCA_log_transcript_Station_Sample_type.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-transcripon")
dev.off()
png(filename="PCA_log_abund_Sample_date_Station.png")
ggplot(PCA_log_abund_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()
png(filename="PCA_log_abund_Station_Sample_type.png")
ggplot(PCA_log_abund_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()

#remove sampleid from colnames
#pca_result <- PCA(t_norm_expression_paired_vir, scale.unit = TRUE, graph = FALSE)
#try PCA for missing values 
#based on https://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html
# estimate number of components
nb <- estim_ncpPCA(t_norm_expression_paired_vir,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)
nb$ncp #2
res.comp <- imputePCA(t_norm_expression_paired_vir, ncp = nb$ncp)
#merge with metadata 
imp <- merge(res.comp$completeObs, metadata[, c("Sample_type", "Station", "Sample_date", "data_type")], by.x = "row.names", by.y = "sampleid") #merge by rownames
res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = nb$ncp, graph=FALSE)

plot(res.pca, hab=12, lab="quali")
#cannot run with Nas values
#there Nas values in each column  
#apply(apply(metadata_num,2, is.na),2,sum)
# Create a biplot using FactoMineR and factoextra
#try to overcome this with some fancy math
fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "blue", # Variables color
                col.ind = iris$Station, # Individuals color by groups
                palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = TRUE, ellipse.level = 0.95)



#add metadata to PCAs 
PCA_expression <- merge(PCA_expression, metadata)
PCA_vir_transcription <- merge(PCA_vir_transcription, metadata)
PCA_vir_abund <- merge(PCA_vir_abund, metadata)

#PCA expression
png(filename="PCA_expression_Sample_type_Station.png")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Sample_type, shape=Station)) + geom_point() + ggtitle("Viral expression")
dev.off()
png(filename="PCA_expression_Station_Sample_type.png")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral expression")
dev.off()
#November much separated form Mai and Jun - can also batch effect
png(filename="PCA_expression_Sample_date_Station.png")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral expression")
dev.off()

#sanity check
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=data_type, shape=Station)) + geom_point() + ggtitle("Viral expression")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Stromkilometer)) + geom_point() + ggtitle("Viral expression")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Sample_type)) + geom_point() + ggtitle("Viral expression")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Sat_O2_TBDHereon)) + geom_point() + ggtitle("Viral expression")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=Temperature_TBDHereon)) + geom_point() + ggtitle("Viral expression")
PCA_expression$Salinity_TBDHereon <- as.numeric(sub(",", ".", PCA_expression$Salinity_TBDHereon)) 
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=DOC_mg.L)) + geom_point() + ggtitle("Viral expression")
ggplot(PCA_expression, aes(x=PC1, y=PC2, colour=DIC_mg.L)) + geom_point() + ggtitle("Viral expression")

#PCA abundance
png(filename="PCA_abund_Sample_type.png")
ggplot(PCA_vir_abund, aes(x=PC1, y=PC2, colour=Sample_type)) + geom_point() + ggtitle("Viral gene counts")
dev.off()
png(filename="PCA_abund_Sample_type.png")
ggplot(PCA_vir_abund, aes(x=PC1, y=PC2, colour=Station)) + geom_point() + ggtitle("Viral gene counts")
dev.off()
png(filename="PCA_abund_data_type.png")
ggplot(PCA_vir_abund, aes(x=PC1, y=PC2, colour=data_type)) + geom_point() + ggtitle("Viral gene counts")
dev.off()

#PCA transcription
png(filename="PCA_transcription_sample_type.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=Sample_type)) + geom_point() + ggtitle("Viral transcription")
dev.off()
#outliers 
#PC1=-1E07 PC2=-3E6
#check what it is
PCA_vir_transcription[which(PCA_vir_transcription$PC1<(-1*10^7)),]
PCA_vir_transcription[which(PCA_vir_transcription$PC2<(-3*10^6)),]
outliers <- PCA_vir_transcription[which(PCA_vir_transcription$PC2<(-3*10^6)  | PCA_vir_transcription$PC1<(-1*10^7)),]$sampleid

#lets excluded it those - PCA is very sensitive to outliers 
png(filename="PCA_abund_Station.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=Station)) + geom_point() + ggtitle("Viral transcription")
dev.off()
png(filename="PCA_abund_Sample_date.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=Sample_date)) + geom_point() + ggtitle("Viral transcription")
dev.off()
png(filename="PCA_abund_data_typee.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=data_type)) + geom_point() + ggtitle("Viral transcription")
dev.off()
save.image("viruses_PCA.RData")


####Split into fractions
#read objects
log_paired_vir_metag_paired <- readRDS("log_paired_vir_metag_paired.RDS")
log_paired_vir_metat_paired <- readRDS("log_paired_vir_metat_paired.RDS")
norm_expression_paired_vir <- readRDS("norm_expression_paired_vir.RDS")
#transpose
t_norm_expression_paired_vir <- t(norm_expression_paired_vir)
t_log_paired_vir_metag_paired <- t(log_paired_vir_metag_paired)
t_log_paired_vir_metat_paired <- t(log_paired_vir_metat_paired)
#remove outliers 
t_norm_expression_paired_vir <- t_norm_expression_paired_vir[which(!rownames(t_norm_expression_paired_vir) %in% outliers),]
t_log_paired_vir_metag_paired <- t_log_paired_vir_metag_paired[which(!rownames(t_log_paired_vir_metag_paired) %in% outliers),]
t_log_paired_vir_metat_paired <- t_log_paired_vir_metat_paired[which(!rownames(t_log_paired_vir_metat_paired) %in% outliers),]


#add metadata to split data into fractions
#t_log_paired_vir_metag_paired <- merge(t_log_paired_vir_metag_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames
#t_log_paired_vir_metat_paired <- merge(t_log_paired_vir_metat_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames

t_norm_expression_paired_vir_meta <- merge(t_norm_expression_paired_vir, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_norm_expression_paired_vir_meta) <- t_norm_expression_paired_vir_meta$Row.names
t_norm_expression_paired_vir_meta[,c("Row.names")] <- NULL
t_log_paired_vir_metag_paired_meta <- merge(t_log_paired_vir_metag_paired, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_log_paired_vir_metag_paired_meta) <- t_log_paired_vir_metag_paired_meta$Row.names
t_log_paired_vir_metag_paired_meta[,c("Row.names")] <- NULL
t_log_paired_vir_metat_paired_meta <- merge(t_log_paired_vir_metat_paired, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_log_paired_vir_metat_paired_meta) <- t_log_paired_vir_metat_paired_meta$Row.names
t_log_paired_vir_metat_paired_meta[,c("Row.names")] <- NULL
#expression
PCA_log_expression <- PCA(t_norm_expression_paired_vir, graph=F)
PCA_log_expression.eig.val <- get_eigenvalue(PCA_log_expression)
#14% an 13% for PC1 ad PC2 respectively 
#expression FL
log_expression_FL <- t_norm_expression_paired_vir_meta[which(t_norm_expression_paired_vir_meta$Sample_type == "Free_living"),!names(t_norm_expression_paired_vir_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_expression_FL <- PCA(log_expression_FL, graph=F)
PCA_log_expression_FL.eig.val <- get_eigenvalue(PCA_log_expression_FL)
#higher variance explained
PCA_log_expression_FL_meta <-  merge(PCA_log_expression_FL$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_expression_FL_meta, file="PCA_log_expression_FL.RDS")
png(filename="PCA_log_expression_FL_sample_date_station.png")
ggplot(PCA_log_expression_FL_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-expression Free-living") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_expression_FL.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_expression_FL.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()

#expression Light_fraction
log_expression_LF <- t_norm_expression_paired_vir_meta[which(t_norm_expression_paired_vir_meta$Sample_type == "Light_fraction"),!names(t_norm_expression_paired_vir_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_expression_LF <- PCA(log_expression_LF, graph=F)
PCA_log_expression_LF.eig.val <- get_eigenvalue(PCA_log_expression_LF)
PCA_log_expression_LF_meta <-  merge(PCA_log_expression_LF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_expression_LF_meta, file="PCA_log_expression_LF.RDS")
png(filename="PCA_log_expression_LF_sample_date_station.png")
ggplot(PCA_log_expression_LF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-expression Light Fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_expression_LF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_expression_LF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#expression Heavy_fraction
log_expression_HF <- t_norm_expression_paired_vir_meta[which(t_norm_expression_paired_vir_meta$Sample_type == "Heavy_fraction"),!names(t_norm_expression_paired_vir_meta) %in% c("Row.names","Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_expression_HF <- PCA(log_expression_HF, graph=F)
PCA_log_expression_HF.eig.val <- get_eigenvalue(PCA_log_expression_HF)
PCA_log_expression_HF_meta <-  merge(PCA_log_expression_HF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_expression_HF_meta, file="PCA_log_expression_HF.RDS")
png(filename="PCA_log_expression_HF_sample_date_station.png")
ggplot(PCA_log_expression_HF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-expression Heavy Fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_expression_HF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_expression_HF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#log gene counts
PCA_log_gene_abund <- PCA(t_log_paired_vir_metag_paired, graph=F)
PCA_log_gene_abund.eig.val <- get_eigenvalue(PCA_log_gene_abund)
#gene abund  FL
log_gene_abund_FL <- t_log_paired_vir_metag_paired_meta[which(t_log_paired_vir_metag_paired_meta$Sample_type == "Free_living"),!names(t_log_paired_vir_metag_paired_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_FL <- PCA(log_gene_abund_FL, graph=F)
PCA_log_gene_abund_FL.eig.val <- get_eigenvalue(PCA_log_gene_abund_FL)
#higher variance explained
PCA_log_gene_abund_FL_meta <-  merge(PCA_log_gene_abund_FL$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_FL_meta, file="PCA_log_gene_abund_FL.RDS")
png(filename="PCA_log_expression_FL_sample_date_station.png")
ggplot(PCA_log_gene_abund_FL_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell Free-living") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_expression_FL.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_expression_FL.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()

#log transcription