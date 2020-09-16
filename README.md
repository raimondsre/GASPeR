# GASPeR: Genome wide Association Studie and Polygenic Risk analysis pipeline software with cluster support
Nextflow pipeline for full GWAS quality control, imputation, association analysis and polygenic risk score model creation&amp;visualisation

To perform full analysis, follow [environment preparation](https://github.com/raimondsre/GWAS-PRS-Piepeline/blob/master/Environment_preparation.md) steps and simply type in the following command line:

```$ nextflow pipeline.nf --input plink```

where ```plink``` is a prefix of .bim .bed .map files.



![Manhattan_plot](https://github.com/raimondsre/GWAS-PRS-Piepeline/blob/master/github_example.png?raw=true)

_One of the outputs, namely, manhattan plot of association summary statistics. Circle size signifies odds ratio._

### Acknowledgement

Parameters as well as scripts supporting pipeline were selected according to Coleman et al. (2016) and supporting repository [gwas_scripts](https://github.com/JoniColeman/gwas_scripts) by [Joni Coleman](https://github.com/JoniColeman). Imputation part of main pipeline contains significant code snippets from [InSilicoDB](https://github.com/InSilicoDB/snp-imputation-nf) by [GenePlaza](https://github.com/InSilicoDB).
