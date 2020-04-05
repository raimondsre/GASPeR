root=$1
root2=$2
bin=$3
pheno=$4


convertf -p <(printf "genotypename: $root.bed
snpname: $root.bim
indivname: $root.fam
outputformat: EIGENSTRAT
genotypeoutname: $root.pop_strat.eigenstratgeno
snpoutname: $root.pop_strat.snp
indivoutname: $root.pop_strat.ind")

smartpca.perl \
-i $root.pop_strat.eigenstratgeno \
-a $root.pop_strat.snp \
-b $root.pop_strat.ind \
-o $root.pop_strat.pca \
-p $root.pop_strat.plot \
-e $root.pop_strat.eval \
-l $root.pop_strat_smartpca.log \
-m 0 \
-t 100 \
-k 100 \
-s 6

sed -i -e 's/^[ \t]*//' -e 's/:/ /g' $root.pop_strat.pca.evec

R --file="${bin}PC-VS-OUTCOME_IN_R_SHORT.R" --args $root.pop_strat "${pheno}"


R --file="${bin}PC-VS-OUTCOME_IN_R_FULL.R" --args $root.pop_strat "${pheno}"


smartpca.perl \
-i $root.pop_strat.eigenstratgeno  \
-a $root.pop_strat.snp        \
-b $root.pop_strat.ind         \
-o $root.pop_strat_outliers.pca \
-p $root.pop_strat_outliers.plot \
-e $root.pop_strat_outliers.eval  \
-l $root.pop_strat_outliers_smartpca.log \
-m 5 \
-t 3 \
-k 100 \
-s 6

sed -i -e 's/^[ \t]*//' -e 's/:/ /g' $root.pop_strat_outliers.pca.evec
R --file="${bin}PlotPCs.R" --args $root.pop_strat 1 2
R --file="${bin}PlotPCs.R" --args $root.pop_strat_outliers 1 2

awk '/REMOVED/ {print $3}' $root.pop_strat_outliers_smartpca.log | sed 's/:/ /g' > $root.pop_strat_outliers.outliers

plink \
--bfile $root \
--remove $root.pop_strat_outliers.outliers \
--make-bed \
--out $root.LD_pop_strat

plink \
--bfile $root2 \
--remove $root.pop_strat_outliers.outliers \
--make-bed \
--out $root2.pop_strat

convertf -p <(printf "genotypename: $root.LD_pop_strat.bed
snpname: $root.LD_pop_strat.bim
indivname: $root.LD_pop_strat.fam
outputformat: EIGENSTRAT
genotypeoutname: $root.PCS_for_covariates.eigenstratgeno
snpoutname: $root.PCS_for_covariates.snp
indivoutname: $root.PCS_for_covariates.ind")

smartpca.perl \
-i $root.PCS_for_covariates.eigenstratgeno \
-a $root.PCS_for_covariates.snp \
-b $root.PCS_for_covariates.ind \
-o $root.PCS_for_covariates.pca \
-p $root.PCS_for_covariates.plot \
-e $root.PCS_for_covariates.eval \
-l $root.PCS_for_covariates_smartpca.log \
-m 0 \
-t 100 \
-k 100 \
-s 6 \

#R --file="${bin}PC-VS-OUTCOME_IN_R_SHORT.R" --args "${root}.PCS_for_covariates"

