args <- commandArgs(TRUE)
root <- args[1]
pheno <- args[2]
PCAEVEC<-read.table(paste(root,".pca.evec",sep=""),head=T)
PCAEVEC[,1] <- PCAEVEC[,2]
colnames(PCAEVEC)<-c("FID","IID","PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40", "PC41", "PC42", "PC43", "PC44", "PC45", "PC46", "PC47", "PC48", "PC49", "PC50", "PC51", "PC52", "PC53", "PC54", "PC55", "PC56", "PC57", "PC58", "PC59", "PC60", "PC61", "PC62", "PC63", "PC64", "PC65", "PC66", "PC67", "PC68", "PC69", "PC70", "PC71", "PC72", "PC73", "PC74", "PC75", "PC76", "PC77", "PC78", "PC79", "PC80", "PC81", "PC82", "PC83", "PC84", "PC85", "PC86", "PC87", "PC88", "PC89", "PC90", "PC91", "PC92", "PC93", "PC94", "PC95", "PC96", "PC97", "PC98", "PC99", "PC100", "Pheno")
PHENOTYPE<-read.table(pheno,head=T)
PCAPHENO<-merge(PCAEVEC,PHENOTYPE, by = "FID")
colnames(PCAPHENO)[c(2,105)]<- c("IID","Pheno")
PCAPHENO[,c(104,106)]<- NULL

sink(paste(root,".PC_Output_Associations_FULL.txt",sep=""))
for (i in 1:100) {
DATA<-as.data.frame(PCAPHENO[,c(3:(i+2))])
print(summary(lm(PCAPHENO[,104] ~ ., data=DATA)))
}
sink()
