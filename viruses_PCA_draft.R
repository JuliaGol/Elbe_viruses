#run PCA on the server for a huge object
library(ggplot2)
setwd("../PCA_vir")
#PCA on expression
norm_expression_paired_vir <- readRDS(file="norm_expression_paired_vir.RDS")
t_norm_expression_paired_vir<- t(norm_expression_paired_vir)
head(t_norm_expression_paired_vir)
B <- t_norm_expression_paired_vir
B_svd <- svd(B)
B_V <- B_svd$v
B_Sigma <- B_svd$d
B_U <- B_svd$u
SG <- B %*% B_V
SG <- data.frame(SG)
colnames(SG) <- paste0('PC', 1:length(SG[1,])) 
#save as csv file
write.csv(SG, file = "PCA_expression.csv")
rownames(SG) <- rownames(t_norm_expression_paired_vir)
png(filename="var_explain_pc.png")
plot(svd1$d^2/sum(svd1$d^2), pch=19, xlab="Singluar vector", ylab="Variance explained")
dev.off()
#transcription
vir_metat <- readRDS(file="vir_metat.RDS")
t_vir_metat<- t(vir_metat)
head(t_vir_metat)
B <- t_vir_metat
B_svd <- svd(B)
B_V <- B_svd$v
B_Sigma <- B_svd$d
B_U <- B_svd$u
SG <- B %*% B_V
SG <- data.frame(SG)
colnames(SG) <- paste0('PC', 1:length(SG[1,])) 
#save as csv file
write.csv(SG, file = "PCA_vir_transcription.csv")
rownames(SG) <- rownames(t_vir_metat)
#copies per cell/abundance
vir_metag <- readRDS(file="vir_metag.RDS")
t_vir_metag<- t(vir_metag)
head(t_vir_metag)
B <- t_vir_metag
B_svd <- svd(B)
B_V <- B_svd$v
B_Sigma <- B_svd$d
B_U <- B_svd$u
SG <- B %*% B_V
SG <- data.frame(SG)
colnames(SG) <- paste0('PC', 1:length(SG[1,])) 
#save as csv file
write.csv(SG, file = "PCA_vir_abund.csv")
rownames(SG) <- rownames(t_vir_metag)
#MCP - major capsid protein - abundance 
#MCP - transcription
#MCP - expression
