#setup
library(ggplot2)
library(FactoMineR)
library(factoextra)
#library(missMDA)
#set up import and upload files/objects
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

log_paired_vir_metag_paired <- readRDS("log_paired_vir_metag_paired.RDS")
log_paired_vir_metat_paired <- readRDS("log_paired_vir_metat_paired.RDS")
log_vir_metag <- readRDS("log_vir_metag.RDS")
norm_expression_paired_vir <- readRDS("norm_expression_paired_vir.RDS")
#transpose
t_norm_expression_paired_vir <- t(norm_expression_paired_vir)
t_log_paired_vir_metag_paired <- t(log_paired_vir_metag_paired)
t_log_paired_vir_metat_paired <- t(log_paired_vir_metat_paired)
t_log_vir_metag <- t(log_vir_metag)

#PCA
PCA_log_expression <- PCA(t_norm_expression_paired_vir, graph=F)
PCA_log_transcript <- PCA(t_log_paired_vir_metat_paired, graph=F)
PCA_log_abund <- PCA(t_log_paired_vir_metag_paired, graph=F)
PCA_log_abund_all <- PCA(t_log_vir_metag, graph=F)

saveRDS(PCA_log_expression, file="PCA_log_expression.RDS")
saveRDS(PCA_log_transcript, file="PCA_log_transcript.RDS")
saveRDS(PCA_log_abund, file="PCA_log_abund.RDS")


#take only numeric data 
#metadata_num <-metadata[,c("Sat_O2_TBDHereon", "Temperature_TBDHereon", "Salinity_TBDHereon", "Turbidity_TBDHereon", "pH_TBDHereon", "O2_TBDHereon", "SPM_mgperL", "DOC_mg.L", "TN_mg.L",
                           "DIC_mg.L", "Silicate_mg.L", "Ammonium_mg.L", "Nitrate_mg.L", "Total_DIN_µM", "TotalDissolvedPhosphate_mg.L", "Ammonium_µM",
                           "Nitrite_µM", "Nitrate_µM", "SRP_µM", "Phosphate_µM", "RespirationRate_O2ug.L.h", "POC_mgperL", "PTC_mgperL", "PTN_mgperL",
                           "PTH_mgperL", "dCH4_nM", "Chlorophyll", "dCO2_uM")]
#and ensure they have point instead of comma 
#and there are set as numeric
# metadata_num <- apply(apply(metadata_num, 2, gsub, patt=",", replace="."), 2, as.numeric)
# metadata_num <- cbind(metadata_num, metadata$sampleid) 
# colnames(metadata_num)[length(colnames(metadata_num))] <- "sampleid"   
# 
# #numerical data merged 
# t_norm_expression_paired_vir <- merge(t_norm_expression_paired_vir, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames
# t_log_paired_vir_metag_paired <- merge(t_log_paired_vir_metag_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames
# t_log_paired_vir_metat_paired <- merge(t_log_paired_vir_metat_paired, metadata_num, by.x = "row.names", by.y = "sampleid") #merge by rownames

#PCA with metadata for plotting
PCA_log_transcript_meta <- merge(PCA_log_transcript$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
PCA_log_abund_meta <- merge(PCA_log_abund$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
PCA_log_abund_all_meta <- merge(PCA_log_abund_all$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
PCA_log_abund_all_meta <- PCA_log_abund_all_meta[which(PCA_log_abund_all_meta$Station != "Seemanshöft" & PCA_log_abund_all_meta$Station != "Kollmar"),] 
PCA_log_expression_meta <- merge(PCA_log_expression$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames 

#visualise log abundance 
png(filename="PCA_plots/PCA_log_transcript_Sample_type_Station.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_type, shape=Station)) + geom_point() + ggtitle("Viral log-transcription")
dev.off()
png(filename="PCA_plots/PCA_log_transcript_Sample_date_Station.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-transcription")
dev.off()
png(filename="PCA_plots/PCA_log_transcript_Station_Sample_date.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_date)) + geom_point() + ggtitle("Viral log-transcription")
dev.off()
png(filename="PCA_plots/PCA_log_transcript_Station_Sample_type.png")
ggplot(PCA_log_transcript_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-transcription")
dev.off()
png(filename="PCA_plots/PCA_log_abund_Sample_date_Station.png")
ggplot(PCA_log_abund_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()
png(filename="PCA_plots/PCA_log_abund_Station_Sample_type.png")
ggplot(PCA_log_abund_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()
png(filename="PCA_plots/PCA_log_abund_all_Sample_date_Station.png")
ggplot(PCA_log_abund_all_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()
png(filename="PCA_plots/PCA_log_abund_all_Station_Sample_type.png")
ggplot(PCA_log_abund_all_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-gene counts")
dev.off()
png(filename="PCA_plots/PCA_log_expression_Sample_date_Station.png")
ggplot(PCA_log_expression_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) + geom_point() + ggtitle("Viral log-expression counts")
dev.off()
png(filename="PCA_plots/PCA_log_expression_Station_Sample_type.png")
ggplot(PCA_log_expression_meta, aes(x=Dim.1, y=Dim.2, colour=Station, shape=Sample_type)) + geom_point() + ggtitle("Viral log-expression counts")
dev.off()
#remove some objects to empty some space 
rm(PCA_log_transcript)
rm(PCA_log_transcript_meta)
rm(PCA_log_expression)
rm(PCA_log_expression_meta)
rm(PCA_log_abund)
rm(PCA_log_abund_meta)
rm(PCA_log_abund_all)
rm(PCA_log_abund_all_meta)
####Split into fractions
t_norm_expression_paired_vir_meta <- merge(t_norm_expression_paired_vir, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_norm_expression_paired_vir_meta) <- t_norm_expression_paired_vir_meta$Row.names
t_norm_expression_paired_vir_meta[,c("Row.names")] <- NULL
t_log_paired_vir_metag_paired_meta <- merge(t_log_paired_vir_metag_paired, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_log_paired_vir_metag_paired_meta) <- t_log_paired_vir_metag_paired_meta$Row.names
t_log_paired_vir_metag_paired_meta[,c("Row.names")] <- NULL
t_log_paired_vir_metat_paired_meta <- merge(t_log_paired_vir_metat_paired, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_log_paired_vir_metat_paired_meta) <- t_log_paired_vir_metat_paired_meta$Row.names
t_log_paired_vir_metat_paired_meta[,c("Row.names")] <- NULL
t_log_vir_metag_meta <- merge(t_log_vir_metag, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid") #merge by rownames
rownames(t_log_vir_metag_meta) <- t_log_vir_metag_meta$Row.names
t_log_paired_vir_metag_paired_meta[,c("Row.names")] <- NULL
#expression
PCA_log_expression.eig.val <- get_eigenvalue(PCA_log_expression)
#expression FL
log_expression_FL <- t_norm_expression_paired_vir_meta[which(t_norm_expression_paired_vir_meta$Sample_type == "Free_living"),!names(t_norm_expression_paired_vir_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_expression_FL <- PCA(log_expression_FL, graph=F)
PCA_log_expression_FL.eig.val <- get_eigenvalue(PCA_log_expression_FL)
#higher variance explained
PCA_log_expression_FL_meta <-  merge(PCA_log_expression_FL$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_expression_FL_meta, file="PCA_log_expression_FL.RDS")
png(filename="PCA_plots/PCA_log_expression_FL_sample_date_station.png")
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
png(filename="PCA_plots/PCA_log_expression_LF_sample_date_station.png")
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
png(filename="PCA_plots/PCA_log_expression_HF_sample_date_station.png")
ggplot(PCA_log_expression_HF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-expression Heavy Fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_expression_HF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_expression_HF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#paired abundance
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
png(filename="PCA_plots/PCA_log_gene_abund_FL_sample_date_station.png")
ggplot(PCA_log_gene_abund_FL_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell Free-living") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_FL.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_FL.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#gene abund  LF
log_gene_abund_LF <- t_log_paired_vir_metag_paired_meta[which(t_log_paired_vir_metag_paired_meta$Sample_type == "Light_fraction"),!names(t_log_paired_vir_metag_paired_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_LF <- PCA(log_gene_abund_LF, graph=F)
PCA_log_gene_abund_LF.eig.val <- get_eigenvalue(PCA_log_gene_abund_LF)
#higher variance explained
PCA_log_gene_abund_LF_meta <-  merge(PCA_log_gene_abund_LF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_LF_meta, file="PCA_log_gene_abund_FL.RDS")
png(filename="PCA_plots/PCA_log_gene_abund_LF_sample_date_station.png")
ggplot(PCA_log_gene_abund_LF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell light fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_LF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_LF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#gene abund  HF
log_gene_abund_HF <- t_log_paired_vir_metag_paired_meta[which(t_log_paired_vir_metag_paired_meta$Sample_type == "Heavy_fraction"),!names(t_log_paired_vir_metag_paired_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_HF <- PCA(log_gene_abund_HF, graph=F)
PCA_log_gene_abund_HF.eig.val <- get_eigenvalue(PCA_log_gene_abund_HF)
#higher variance explained
PCA_log_gene_abund_HF_meta <-  merge(PCA_log_gene_abund_HF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_HF_meta, file="PCA_log_gene_abund_FL.RDS")
png(filename="PCA_plots/PCA_log_gene_abund_HF_sample_date_station.png")
ggplot(PCA_log_gene_abund_HF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell heavy fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_HF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_HF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#all abundance

#gene abund  FL
log_gene_abund_all_FL <- t_log_vir_metag_meta[which(t_log_vir_metag_meta$Sample_type == "Free_living"),!names(t_log_vir_metag_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_all_FL <- PCA(log_gene_abund_all_FL, graph=F)
PCA_log_gene_abund_all_FL.eig.val <- get_eigenvalue(PCA_log_gene_abund_all_FL)
#higher variance explained
PCA_log_gene_abund_all_FL_meta <-  merge(PCA_log_gene_abund_all_FL$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_all_FL_meta, file="PCA_log_gene_abund_all_FL.RDS")
png(filename="PCA_plots/PCA_log_gene_abund_all_FL_sample_date_station.png")
ggplot(PCA_log_gene_abund_all_FL_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell Free-living") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_all_FL.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_all_FL.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#gene abund  LF
log_gene_abund_all_LF <- t_log_vir_metag_meta[which(t_log_vir_metag_meta$Sample_type == "Light_fraction"),!names(t_log_vir_metag_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_all_LF <- PCA(log_gene_abund_all_LF, graph=F)
PCA_log_gene_abund_all_LF.eig.val <- get_eigenvalue(PCA_log_gene_abund_all_LF)
#higher variance explained
PCA_log_gene_abund_all_LF_meta <-  merge(PCA_log_gene_abund_all_LF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_all_LF_meta, file="PCA_log_gene_abund_all_FL.RDS")
png(filename="PCA_plots/PCA_log_gene_abund_all_LF_sample_date_station.png")
ggplot(PCA_log_gene_abund_all_LF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell light fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_all_LF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_all_LF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
#gene abund  HF
log_gene_abund_all_HF <- t_log_vir_metag_meta[which(t_log_vir_metag_meta$Sample_type == "Heavy_fraction"),!names(t_log_vir_metag_meta) %in% c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")]
PCA_log_gene_abund_all_HF <- PCA(log_gene_abund_all_HF, graph=F)
PCA_log_gene_abund_all_HF.eig.val <- get_eigenvalue(PCA_log_gene_abund_all_HF)
#higher variance explained
PCA_log_gene_abund_all_HF_meta <-  merge(PCA_log_gene_abund_all_HF$ind$coord, metadata[, c("Sample_type", "Station", "Sample_date", "data_type", "sampleid")], by.x = "row.names", by.y = "sampleid")
saveRDS(PCA_log_gene_abund_all_HF_meta, file="PCA_log_gene_abund_all_FL.RDS")
png(filename="PCA_plots/PCA_log_gene_abund_all_HF_sample_date_station.png")
ggplot(PCA_log_gene_abund_all_HF_meta, aes(x=Dim.1, y=Dim.2, colour=Sample_date, shape=Station)) +
  geom_point() +
  ggtitle("Viral log-gene counts per cell heavy fraction") + 
  xlab(paste0(paste0("PC1 ", as.character(round(PCA_log_gene_abund_all_HF.eig.val[1,c("variance.percent")],2))), "%")) +
  ylab(paste0(paste0("PC2 ", as.character(round(PCA_log_gene_abund_all_HF.eig.val[2,c("variance.percent")],2))), "%"))
dev.off()
