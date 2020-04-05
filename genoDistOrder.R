library(snpStats)
args <- commandArgs(TRUE)
root <- args[1]



#converting plink to R acceptable format
pathC <- paste(root, c(".bed", ".bim", ".fam" ), sep = "")
RigaGWAS <- read.plink(pathC[1], pathC[2], pathC[3])
colnames(RigaGWAS$fam) <- c("pedigree", "member", "father", "mother", "sex", "affected")
colnames(RigaGWAS$map) <- c("chr", "SNP", "location", "position", "A1", "A2")
rownames(RigaGWAS$genotypes) <- RigaGWAS$fam$member
#producing dataset for further analysis
SNP <- RigaGWAS$genotypes
MAP <- RigaGWAS$map
FAM <- RigaGWAS$fam
genData <- list(SNP = SNP, MAP = MAP, FAM = FAM)

#skripts, kas izveido sarakstu or SNP ar nperecizu genoma distancu secibu
 b <- data.frame()
      for (i in 1:length(unique(genData$MAP$chr))) {
        a <- genData$MAP[genData$MAP$chr == i,]
        print(i)
        for (j in 1:(length(a$location)-1)) {
	  tryCatch({
	  if (is.na(a$location[j])) {
            b <- rbind(b,a[j+1,])}
          else if (a$location[j] > a$location[j+1]) {
            k <- j+1
            while (length(a$location) >= k & a$location[j] > a$location[k]){
              b <- rbind(b,a[k,])
              k <- k+1
             }
            }
	}, error=function(e){})
        }
      }
      c <- b$SNP[!duplicated(b$SNP)]
      genData$MAP <- genData$MAP[!genData$MAP$SNP %in% c,]

#atkartojam veilreiz, lai parbauditu pari palikušos
      b <- data.frame()
      for (i in 1:length(unique(genData$MAP$chr))) {
        a <- genData$MAP[genData$MAP$chr == i,]
        print(i)
        for (j in 1:(length(a$location)-1)) {
	tryCatch({ 
         if (is.na(a$location[j])) {
            b <- rbind(b,a[j+1,])}
          else if (a$location[j] > a$location[j+1]) {
            k <- j+1
            while (length(a$location) >= k & a$location[j] > a$location[k]){
              b <- rbind(b,a[k,])
              k <- k+1
            }
          }
	}, error=function(e){})
        }
      }
      c <- rbind(c, b$SNP[!duplicated(b$SNP)])
      genData$MAP <- genData$MAP[!genData$MAP$SNP %in% c,]
      
      write.table(c, paste(root,".snpoutoforder.txt", sep=""), quote=F)
      
      
