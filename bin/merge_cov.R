args <- commandArgs(TRUE)

cov <- args[1]
pc <- args[2]

covar <- read.table(cov,header=T)
pcs <- read.table(pc)

covar_pcs <- merge(covar,pcs,by=c(1,2))

write.table(covar_pcs, file="cov_pcs.txt",row.names=FALSE, quote=FALSE)




