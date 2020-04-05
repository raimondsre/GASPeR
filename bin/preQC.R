#Running example:
#If Phenotype other than third column of pheno file, EDIT THE CODE. Else:

#R --file=/home/raimondsre/GWAS_pipeline/bin/preQC.R --args RigaGWAS_2017.fam FULL_pheno_T2D.txt pheno_file #this command will take in RigaGWAS_2017.fam file as well as file names RigaGWAS_2017.pheno.txt. If second argument is specified, programm will look for column name called exactly as such and concatenate this column with from pheno file with .fam file. Third argument is phenotype file path. If this is not specified, programm assumes file is equal to first argument + .pheno.txt.

#Additional code is left to enhance customisability



args <- commandArgs(TRUE)
root <- args[1]
pheno_name <- args[3]
pheno <- args[2]


if (is.na(pheno)) {pheno <- read.delim(paste(root, ".pheno.txt", sep=""), header = TRUE, sep=" ")} else {pheno <- read.delim(pheno, header = TRUE, sep=" ")}
if (is.na(pheno_name)) {pheno_name <- colnames(pheno)[3]}
fam <- read.table(root, sep = "\t")
head(fam)
#Arrange rows of pheno so that order matches .fam file one. If pheno file lacks person, "NA" value will be assigned.
newPheno <- pheno[match(fam[,2], pheno[,2]),]
#Column-bind .fam and pheno files. .fam file take all columns till "sex", while pheno file take in all columns after phenotype column.
fam <- cbind(fam[,c(1:5)],newPheno[,which(colnames(pheno)==pheno_name)])
fam[,c(3:4)] <- 0
#Copy IID to FID to replace just the order that was there before. It will help in combining datasets.
fam[,1] <- fam[,2]
write.table(fam,file = root, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)




