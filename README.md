# FULL GWAS & PRS pipeline
Nextflow pipeline for full GWAS quality control, imputation, association analysis and PRS model creation&amp;visualisation

To perform full analysis, prepare environment and simply type in the following command line:

```$ nextflow pipeline.nf --input plink```

where ```plink``` is a prefix of .bim .bed .map files.



![Manhattan_plot](https://github.com/raimondsre/GWAS-PRS-Piepeline/blob/master/github_example.png?raw=true)

On of the outputs, namely, manhattan plot of association summary statistics. Circle size signify odds ratio.
