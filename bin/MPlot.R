# List of packages for session
.packages = c("ggplot2","ggrepel","dplyr")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#if (!require('qqman')) install.packages('qqman')


args <- commandArgs(TRUE)
root <- args[1]
pheno <- args[2]
assoc <- read.table(paste(root, ".whole.post_imp.ChrPosA1A2.post_imputation_conc_analysis.assoc.logistic", sep = ""), header=T)


assoc <- assoc[!is.na(assoc$P),]
gwasResults <- assoc 
gwasResults <- gwasResults[-log(gwasResults$P) > 4,]
highlights <- gwasResults$SNP
gwasResults <- gwasResults[order(gwasResults$P),]
gwasResults$SNP[duplicated(gwasResults$SNP)] <- NA
# Prepare the dataset
{don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(as.numeric(BP))) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( Chromosomes=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% highlights, "yes", "no")) %>%
    mutate( is_annotate=ifelse(-log10(P)>4.5 , "yes", "no")) %>%
    mutate( OR=OR)
  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(Chromosomes) + min(Chromosomes) ) / 2 )
  # Make the plot
  manhattan_plot <- ggplot(don, aes(x=Chromosomes, y=-log10(P))) +
    
    # Show all points
    geom_hline(yintercept=c(-log10(1e-05),-log10(5e-08)), linetype = c(4,2),size= 0.5, color=c("grey", "black"), alpha=(0.8) ) +
    
    geom_point( aes(color=as.factor(CHR)),size=don$OR,  alpha=0.7) +
    
    scale_color_manual(values = rep(c("grey58", "grey16"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0.2)) +     # remove space between plot area and x axis
    
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes", size=don$OR), alpha=0.5) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2, alpha=0.6) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.key=element_rect(fill = NA)
    ) +
    guides(color=FALSE)
}


#QQ PLOT
#Grid for association analysis
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, col ="#E69F00") +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}


inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- (median(chisq) / qchisq(0.5, 1))
  lambda
}

qq_plot <- gg_qqplot(assoc$P) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("λ = %.2f", inflation(assoc$P)),
    size = 3
  ) +
  theme(
    axis.ticks = element_line(size = 0.2),
    panel.grid = element_blank(),
    axis.text=element_text(size=10),
    axis.title=element_text(size=10)
    # panel.grid = element_line(size = 0.5, color = "grey80")
  ) 






pdf(paste(root, "_manhattan_plot.pdf", sep=""))
manhattan_plot
dev.off()

pdf(paste(root, "_qq.pdf", sep=""))
qq_plot
dev.off()
