# Full GWAS & PRS pipeline
Nextflow pipeline for full GWAS quality control, imputation, association analysis and polygenic risk score model creation&amp;visualisation

To perform full analysis, follow [environment preparation](https://github.com/raimondsre/GWAS-PRS-Piepeline/blob/master/Environment_preparation.md) steps and simply type in the following command line:

```$ nextflow pipeline.nf --input plink```

where ```plink``` is a prefix of .bim .bed .map files.



![Manhattan_plot](https://github.com/raimondsre/GWAS-PRS-Piepeline/blob/master/github_example.png?raw=true)

_One of the outputs, namely, manhattan plot of association summary statistics. Circle size signifies odds ratio._
