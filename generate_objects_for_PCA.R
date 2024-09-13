#group by gene_cluster to optimise calculations
#chose only columns which we need
#select columns with metat 
metat_cols <- colnames(vir_genecluster_annot)[which(grepl("METAT", colnames(vir_genecluster_annot)))]
#find paired metag
paired_METAG <- sort(gsub("METAT", "METAG", metat_cols)) 
#check and save those which are existing
paired_METAG <- sort(paired_METAG[which(paired_METAG %in% colnames(vir_genecluster_annot))])
#save respected METAT columns 
paired_METAT <- sort(gsub("METAG", "METAT", paired_METAG))
#check if the same order
sum(gsub("METAG", "METAT", colnames(paired_METAG)) == colnames(paired_METAT)) == length(colnames(paired_METAT))
#select only needed columns and summarise by gene cluster 
vir_genecluster_annot_cols <- colnames(vir_genecluster_annot)[which(grepl("GROS|gene_cluster", colnames(vir_genecluster_annot)))]

vir_genecluster_annot_sum <- vir_genecluster_annot[, vir_genecluster_annot_cols] %>% 
  group_by(gene_cluster) %>% summarise_at(vir_genecluster_annot_cols[-1] , sum)
#
vir_metag  <- vir_genecluster_annot_sum %>% select(contains("METAG"))
vir_metat  <- vir_genecluster_annot_sum %>% select(contains("METAT"))
saveRDS(vir_metag,  file="vir_metag.RDS")
saveRDS(vir_metat,  file="vir_metat.RDS")

vir_genecluster_annot_sum_paired <- vir_genecluster_annot_sum[, c("gene_cluster", paired_METAG, paired_METAT)]
vir_metag_paired <- vir_genecluster_annot_sum_paired %>% select(contains("METAG"))
vir_metat_paired <- vir_genecluster_annot_sum_paired %>% select(contains("METAT"))
#
saveRDS(vir_metag_paired,  file="vir_metag_paired.RDS")
saveRDS(vir_metat_paired,  file="vir_metat_paired.RDS")
#find respective metagenomes for metatranscriptomes vir_metat_paired
paired_METAG <- sort(gsub("METAT", "METAG", colnames(vir_metat_paired))) ####why????
colnames_vir_metag <- sort(colnames(vir_metag_paired))
paired_vir_metag <- vir_metag_paired %>% select(contains(paired_METAG)) #strange some METAT do not have 
#METAG
#then select metat which has paired METAG
paired_METAT <- sort(gsub("METAG", "METAT", colnames(paired_vir_metag_paired)))
paired_vir_metat_paired <- vir_metat_paired[,paired_METAT]
#check the order 
#sort it
#== sub("METAT", "METAG", colnames(vir_metat_paired))
#log transform both
#add small pseudocount avoid problems wth NaN values coming from zeros
log_paired_vir_metag_paired <- log(paired_vir_metag_paired + 0.000001)
log_paired_vir_metat_paired <- log(paired_vir_metat_paired + 0.000001)
#normalization interpreted as gene expression
norm_expression_paired_vir <- log_paired_vir_metat_paired - log_paired_vir_metag_paired
saveRDS(norm_expression_paired_vir, file="norm_expression_paired_vir.RDS")