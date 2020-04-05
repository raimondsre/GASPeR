args <- commandArgs(TRUE)
root <- args[1]
PCAEVEC<-read.table(paste(root,".LD_pop_strat.pca.evec_RENAMED",sep=""), head=T)
colnames(PCAEVEC) <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Pop")

a <- length(unique(PCAEVEC$Pop))

library(ggplot2)
color <- c("#4ec3ff",
"#6a76ff",
"#8ff011",
"#ff74f1",
"#01dd77",
"#cd006b",
"#dce200",
"#001b46",
"#cf7800",
"#c9ff54",
"#7f0010",
"#00f3f3",
"#e40025",
"#5ec400",
"#00453a",
"#734900")
color <- color[1:a]
pdf(paste("ALL_DATASETS.LD_pop_strat_PCA.pdf",sep=""))
 
with(PCAEVEC, qplot(PC1,PC2,colour=Pop,ylim=c(-0.05,0.05),xlim=c(-0.05,0.07))+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC1,PC2,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC1,PC4,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC1,PC5,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC2,PC3,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC2,PC4,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC2,PC5,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC3,PC4,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC3,PC5,colour=Pop)+scale_color_manual(values=color))
with(PCAEVEC, qplot(PC4,PC5,colour=Pop)+scale_color_manual(values=color))
dev.off()
