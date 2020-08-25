# Pipeline purpose:
 - inputs ```plink``` format ```.bim .bed .fam``` file
 - performs:
   - quality control, imputation, assocation test
   - polygenic risk score (PRS) modeling by randomly dividing input in two sets for association test & model creation
 - outputs:
   - summary statistics
   - individual PRS scores
   - manhattan plot & qq plot
   - PRS graphs (e.g. quantile plot)




 
# Building environment

Following steps have to be completed sequentially :heavy_exclamation_mark: 


1. Prepare required data files:
    - ```plink.bim```
      - no duplicated SNPs
      - only autosomes (1 to 22 chromosome)
    - ```plink.fam```
      - no duplicated Individuals
      - first and second columns are equal as no relatives allowed
      - no missing sex allowed
    - ```plink.bed```
    - ```phenotype.txt```
      - contains three colums names "FID", "IID", "Phenotype"
      - binary phenotype (1 - control, 2 - case)
      - no rownames, space separated
    - ```covariates.txt```
      - first two column names are "FID", "IID"
      - no rownames, space separated
      - all covariates present will be included
      - Five principal components will be automatically calculated and added to covariates after ```imputation```

2. Required platform and specs:
    - Pipeline is executed in ```linux```, no root (admin privileges) are necessary
    - Optimally run on cluster. Default parameters expect ```pbs``` system ```Torque 6.1.1.1/ Moab 9.1.1```. Change it in ```nextflow.configure``` files ```executor``` parameter
    - One analysis of 1000 individual ```plink``` files may take up to one Terabite of memory. After analysis, generated files in ```work``` directory must be deleted manually.
    
3. Required software:
    - ```nextflow```19.10+ [get](https://www.nextflow.io/docs/latest/getstarted.html)
    - ```plink``` 1.9 [get](https://www.cog-genomics.org/plink/)
    - ```R``` 3.4+
      - packages ```ggplot2```, ```dplyr```
    - ```impute2```v2 [get](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
    - ```shapeit```v2 [get](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
    - ```Eigensoft``` 6.1.4 [get](https://data.broadinstitute.org/alkesgroup/EIGENSOFT/)
      - ```convertf``` ```smartpca.perl``` tools are used
      - These may not work after installation, here is a possible [fix](https://www.biostars.org/p/173436/)
    - ```python``` 3.7.3
    - ```PRSice``` 2.2+ [get](https://www.prsice.info/#executable-downloads)
      - Source code files ```PRSice_linux```, ```PRSice.R``` must be located in ```working directory``` as pipeline runs these directly

4. Ensure each tool functions properly
5. Choose ```working directory```
6. ```working directory```must contain:
    - ```input_data``` directory
      - which must contain ```plink.bed```, ```plink.bim```, ```plink.fam```, ```covariates.txt```, ```phenotype.txt``` files
    - put ```PRSice_linux``` and ```PRSice.R``` files in ```working directory```
    - put ```nextflow.nf``` pipeline script in ```working directory```
    - 1000 genome reference folder ```ALL.integrated_phase1_SHAPEIT_16-06-14.nomono```
      - first, use ```wget``` with [this](https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz) link
      - use ```tar -vxzf``` to unzip it
    - clone /bin folder here
      - ```chmod +x``` every script file located there

7. Run the pipeline in ```working directory```: 
    - **```$ nextflow pipeline.nf --input plink```**
    
**done**
      
    
    

  

