#install.packages("seqinr")
library(seqinr)
install.packages(c("foreach", "doParallel"))
library(doParallel)
library(foreach)
# specifies that you will be manually setting the cores
options(future.availablecores.methods = "mc.cores")
# SLURM_CPUS_PER_TASK is the amount of cores specified in the job environment
options(mc.cores = Sys.getenv("SLURM_CPUS_PER_TASK"))
cluster <- makeCluster(mc.cores)
registerDoParallel(cluster)
#setwd("C:/Users/jgolebiowska/Documents/IGB_phd/BICEST/virus/viral_genes")
myProts <- read.fasta(file = 'cd /scratch/jgolebiowska/viruses_BICEST/Combined_phagegenes_GT.faa', seqtype = "AA")
headers <- c()
IP <- c()
# Use foreach and %dopar% to run the loop in parallel
results <- foreach(prot=iter(myProts), .combine='rbind') %dopar% {
  # Store the results
  IP <- seqinr::computePI(prot)
  header <- attr(prot,"name")
  x <- c(header, IP)
}
# Don't fotget to stop the cluster
stopCluster(cl = cluster)

df <- data.frame(results)
colnames(df) <- c("header", "IP")
rownames(df) <- NULL
write.csv(df, file = "/scratch/jgolebiowska/viruses_BICEST/Isolectric_proteom.csv")
