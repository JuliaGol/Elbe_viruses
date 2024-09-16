#Check PCA
library(ggplot2)
setwd("C:/Users/jgolebiowska/Documents/IGB_phd/virus/viral_abundances/PCA_outputs")
norm_expression_paired_vir <- readRDS("norm_expression_paired_vir.RDS")
path = "C:/Users/jgolebiowska/Documents/IGB_phd/workshop_gene_catalog"
metadata <- read.csv(paste0(path,"/metadata_all.csv"))
#read PCA results calculated on the server 
PCA_expression <- read.csv("PCA_expression.csv")
PCA_vir_abund <- read.csv("PCA_vir_abund.csv")
PCA_vir_transcription <- read.csv("PCA_vir_transcription.csv")
#set first column as sample ids name
colnames(PCA_expression)[1] <- "sampleid"
colnames(PCA_vir_transcription)[1] <- "sampleid"
colnames(PCA_vir_abund)[1] <- "sampleid"
install.packages("factoextra")
library(FactoMineR)
library(factoextra)
library(missMDA)

# Perform PCA using FactoMineR for normalised expression data - similar to  PCA_expression 
# but environmental data included
t_norm_expression_paired_vir <- t(norm_expression_paired_vir)
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
png(filename="PCA_abund_Station.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=Station)) + geom_point() + ggtitle("Viral transcription")
dev.off()
png(filename="PCA_abund_Sample_date.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=Sample_date)) + geom_point() + ggtitle("Viral transcription")
dev.off()
png(filename="PCA_abund_data_typee.png")
ggplot(PCA_vir_transcription, aes(x=PC1, y=PC2, colour=data_type)) + geom_point() + ggtitle("Viral transcription")
dev.off()
