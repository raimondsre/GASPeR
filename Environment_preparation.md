# Pipeline purpose:
 - inputs ```plink``` format ```.bim .bed .fam``` file
 - performs:
   - quality control, imputation, assocation test
   - polygenic risk score (PRS) modeling by randomly dividing input in two sets for association test & model creation
 - outputs:
   - summary statistics
   - individual PRS scores
   - manhattan plot & qq plot
   - PRS graphs (e.g quantile plot)




 
# Environment preparation

Following steps have to be completed sequentially :heavy_exclamation_mark: 


1. Prepare required data files:
    - ```plink.bim```
      - no duplicated SNPs
      - only autosomes (1 to 22 chromosome)
    - ```plink.fam```
      - no duplicated Individuals
      - first and second columns are equal as no relatives allowed
    - ```plink.bed```
    - ```phenotype.txt```
      - contains three colums names "FID", "IID", "Phenotype"
      - binary phenotype (1 - control, 2 - case)
    - ```covariates.txt```
      - first two column names are "FID", "IID"
      - all covariates present will be included
      - Five principal components will be automatically calculated and added to covariates after ```imputation```

2. Required software:
    - Pipeline is executed in ```linux```, no root (admin privileges) are necessary
    - Should be run on cluster. Default parameters expect ```pbs``` system ```Torque 6.1.1.1/ Moab 9.1.1```. Change it in ```nextflow.configure``` files ```executor``` parameter
    - ```nextflow```19.10+ [get](https://www.nextflow.io/docs/latest/getstarted.html)
    - ```plink``` 1.9 [get](https://www.cog-genomics.org/plink/)
    - ```R``` 3.4+
      - packages ```ggplot2```, ```dplyr```
    - ```impute2```v2 [get](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
    - ```shapeit```v2 [get](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
    - Eigensoft 6.1.4 [get](https://data.broadinstitute.org/alkesgroup/EIGENSOFT/)
      - ```convertf``` ```smartpca.perl``` tools are used
      - These may not work after installation, here is a possible [fix](https://www.biostars.org/p/173436/)
    - ```python```
    - ```PRSice``` 2.2+ [get](https://www.prsice.info/#executable-downloads)
    
    
    
      
    
    
  - No duplicated individuals

a <- read.table()
  

