#Phenotype manipulation
rm(list=ls())
pheno <- read.table("GWAS_all_pheno20.02.2019.txt", sep = "\t")

#Change first row to be header.
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

pheno <- header.true(pheno)
plink <- pheno[,c(1,1,11,25, 26, 13,14:23)]
names(plink)[c(1,2)] <- c("FID", "IID")

write.table(plink, file="pheno.txt", sep = "\t", col.names = TRUE, quote=FALSE)

pcs <- read.table("RigaGWAS_2017.common.filtered.hw_dropped.LD_IBD.pop_strat.PC_for_covariates.txt", header=TRUE)
X <- as.matrix(plink[,c(3:length(names(plink)))])
Xscaled <- scale(as.numeric(X))
a <- cor(t(Xscaled))
heatmap(a)

X <- X[, apply(X, 2, function(i) var(i) > 0)]


pcs[,3]
plot(density(pcs[,3], bw = 0.007))


#Formatted table
install.packages("rlang")
install.packages("sjPlot")
install.packages("stargazer")
install.packages("effects")
install.packages("openxlsx")
install.packages("tangram")
install.packages("table1")
install.packages("MatchIt")
install.packages("data.table")
library(data.table) #to combine ethnicities into one column
library(MatchIt) #for p-value calc in table1
library(boot) #data frame for table1
library(ggplot2)
library(effects)
library(sjPlot) #For scientific tables, ex
library(stargazer) #for scientific tables
library(openxlsx)
library(tangram) #for as.categorical function
library(table1)

gwas <- read.table("PHENO_FULL.txt", header=T)
gwas[gwas== -1] <- NA
gwas[!is.na(gwas[,4]) & gwas[,4]==0 ,4] <- NA
gwas[!is.na(gwas[,5]) & gwas[,5]==0 ,5] <- NA

lv384 <- read.table("Latvia384.fam", sep="\t")
R2017 <- read.table("RigaGWAS_2017.fam", sep = "\t")
expres475 <- read.table("LV475_express.fam")
eurhd475 <- read.table("LV475_eurhd.fam")
df <- list()

for (i in 1:6) df[[i]] <- read.xlsx(sheet = i,"ML_GWAS_controls_atsifreets.xlsx")
df[[1]] <- df[[1]][,c(1:8,ncol(df[[1]]),9:ncol(df[[1]]))]
colnames(df[[1]])[c(9,18,19,20,21)] <- c("U_ALIQUOT_NUMBER","X18","X19","X20","X21")
df[[6]][22] <- NULL
SAMPLE_ID <- Reduce(rbind, df)
write.table(SAMPLE_ID, file="LAB_SAMPLE_ID.txt", sep="\t")
bare_sample_ID <- SAMPLE_ID[,c("SAMPLE_ID", "Gender", "Age", "X18")]
#People who doesn't have a translation from LV**KONLUKD****** code to numbered one in file "ML_GWAS_controls_atsifreets.xlsx"
write.table(eurhd475[!eurhd475[,2] %in% bare_sample_ID$X18,], file="eurhd475_removed.noConversionForSAMPLE_ID.txt", quote=FALSE)
#Replacing LC**KONLUKD**** code to numbered
expres475[,2] <- bare_sample_ID[match(expres475[,2], bare_sample_ID[,"X18"]),"SAMPLE_ID"]
eurhd475[,2] <- as.character(as.factor(eurhd475[,2]))
eurhd475[,2] <- unlist(sapply(eurhd475[,2] , function(i) {if (i %in% bare_sample_ID[,"X18"]) {eurhd475[eurhd475[,2]==i,2] <- subset(bare_sample_ID, X18 == i)[1]} else {i}}))

#Creating .fam file
write.table(expres475,file = "LV475_express.fam", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#Eurhd needs a special treatment to account for NAs
eurhd475[,6] <- gwas[match(eurhd475[,2], gwas$IID),"diagnosis"]
eurhd475[,1] <- eurhd475[,2]
write.table(eurhd475,file = "LV475_eurhd.fam", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(eurhd475[],file = "LV475_eurhd.PresentInConversionTable.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(eurhd475[,c(2,2)], file = "LV475_eurhd.noConversionRemoved.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#plink --bfile LV475_eurhd --keep LV475_eurhd.PresentInConversionTable.txt --make-bed --out plink 
#mv plink.{bim,bed,fam} LV475_eurhd.{bim,bed,fam}

#phenotype table for each dataset
pheno <- list()
pheno[[1]] <- gwas[gwas$IID %in% R2017[,2],]
pheno[[2]] <- gwas[gwas$IID %in% lv384[,2],]
pheno[[3]] <- gwas[gwas$IID %in% expres475[,2],]
pheno[[4]] <- gwas[gwas$IID %in% eurhd475[,2],]
pheno[[5]] <- Reduce(rbind,pheno)


formula <- list()
formula <- sapply(c(1:5),function(i) paste(colnames(pheno[[i]])[18],paste(colnames(pheno[[i]])[c(3:17,19:ncol(pheno[[i]]))], collapse=" + "), sep= " ~ -1 + "))
pheno[[5]][,7] <- as.logical(pheno[[5]][,7])

glm <- lapply(c(5), function(i) glm(formula[[i]], data=pheno[[i]]))

stargazer(pheno[[5]][,c(-1,-2,-17)],
          column.labels = c("RigaGWAS-2017","Latvia384","All"),
          type="text",
          initial.zero = FALSE,
          ord.intercepts = FALSE,
          column.sep.width = "15pt",
          decomal.mark = " ",
          dep.var.caption = "Type II Diabetes mellitus",
          dep.var.labels =  "",
          digits = 2,
          digits.extra = 2,
          keep.stat = c("n"),
          notes = c("\"All\" include Latvia384, LV475-eurhd, LV475-express"),
          omit.table.layout = c( "#"),
          zero.component = FALSE,
          summary.logical = TRUE,
          summary.stat = c("mean","sd","min","max","n"),
          nobs = TRUE,
          mean.sd = TRUE
          #out = "table_full.html"
          )

stargazer(lapplypheno[[i]],
          type="text")










#table1 scintific table for gwas covariates
phenoF <- pheno[[5]]

#Making ethnicity column
#column for all ethnicities
d2 <- phenoF[,c(7:16)]
setDT(d2)
cols <- names(d2)
#update values to colunumbers
d2[, (cols) := as.data.table( ifelse( d2 == '0', NA, col(d2) ) )]
#final output
d2[, .(ethnicity = fcoalesce( d2 ) ) ][]
Ethnicity <- d2[, .(ethnicity = fcoalesce( d2 ) ) ][]
#combine with the phenotype file
phenoEthn <- cbind(phenoF, Ethnicity)

#Converting columns to labeled factors
phenoEthn$diagnosis <- factor(phenoEthn$diagnosis, levels=c(1,2),labels=c("Control", "Case"))
phenoEthn$ethnicity <- factor(phenoEthn$ethnicity, levels=c(1,2,3,4,5,6,7,8,9,10), labels=c("Latv", "Ru","By","Ukr","Polis","Ebrejs", "Liet", "Roma", "Est", "Ger"))
phenoEthn$smoking_status <- factor(phenoEthn$smoking_status)
phenoEthn$sex <- factor(phenoEthn$sex)
str(phenoEthn)

label(phenoEthn$diagnosis) <- "Type II Diabetes mellitus"
label(phenoEthn$ethnicity) <- "Ethnicity"
label(phenoEthn$age) <- "Age"
label(phenoEthn$height) <- "Height"
label(phenoEthn$weight) <- "Weight"
label(phenoEthn$smoking_status) <- "Smoking status"

#phenoEthn[,c(1,2,7:16)]
rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- phenoEthn[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ phenoEthn$diagnosis)$p.value
        } else {
            p <- chisq.test(table(y, droplevels(phenoEthn$diagnosis)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}

rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
}


names(phenoF)
table1(~ sex + age + height +weight + smoking_status + ethnicity
       | diagnosis,
       data=phenoEthn,
       overall= F,
       topclass = "Rtable1-zebra Rtable1-times",
       render = rndr,
       render.strat = rndr.strat
       )

str(phenoE)

