
//Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

params.workingDir = System.getProperty("user.dir")
params.input_dir = "${params.workingDir}/input_data"
params.input = "plink"
params.pheno_file = "phenotype.txt"
params.pheno = "${params.input_dir}/${params.pheno_file}"
inpat = "${params.input_dir}/${params.input}"
params.output_dir = "${params.workingDir}/output_data/${params.input}"
params.output = "${params.workingDir}/output_data"

raw_ch       = Channel.create()
bim_ch       = Channel.create()


bin = "${params.workingDir}/bin/"
binChannel = Channel.value("${bin} ")

process executePermissionForScripts {
        script:
	"""
	
        chmod +x "${bin}Iterative_Missingness.sh" "${bin}highLDregions4bim_b37.awk" "${bin}IndividualIBD.R" "${bin}PCA.sh" "${bin}PC-VS-OUTCOME_IN_R_SHORT.R" "${bin}PC-VS-OUTCOME_IN_R_FULL.R" "${bin}PlotPCs.R" "${bin}IdHets.R" "${bin}makeChunks.sh"
	mkdir -p ${params.output}
	mkdir -p ${params.output_dir}
        """
}

Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }\
      .separate(raw_ch, bim_ch) { a -> [a,a[1]] }

// Filter for common SNPs
process maf {
	publishDir params.output_dir

	input:
	set file(bed), file(bim), file(fam) from raw_ch
	
	output:
	tuple file("${maf}.bed"), file("${maf}.bim"),file("${maf}.fam")\
	into (qc1)
	file("${base}_removed_snp.bim")	

	script:
	base = bed.baseName
	maf = "${base}.common"
	
	"""
	plink --bfile $base --pheno ${params.pheno} --maf 0.01 --make-bed --out $maf
	plink --bfile $base --exclude "${maf}.bim" --make-bed  --out ${base}_removed_snp
	"""
}

// Iterative Missingness, removing SNPs with call rate <90%
process iteratCallRate {
	publishDir params.output_dir
	input:
	set file(bed), file(bim), file(fam) from qc1
	output:
	tuple file("${im}.bed"), file("${im}.bim"), file("${im}.fam")\
	into (qc2, qc2_A) 
 	file "${base}_removed_snp_im.bim"
	script:
	base = bed.baseName
	im = "${base}.filtered"
	"""
	${bin}Iterative_Missingness.sh 90 99 1  $base 
	plink --bfile $base --exclude "${base}.filtered.bim" --make-bed --out ${base}_removed_snp_im
	"""	
}

//checking if all missing SNPs and individuals have been dropped, first check missing rates per Individual & per genome
process checkIfMissingDroped {
        publishDir params.output_dir
	input:
	set file(bed), file(bim), file(fam) from qc2_A
	
	output:
	file("${base}_missing.lmiss") 
	file("${base}_missingness.head")
	script:
	base = bed.baseName
	
	"""
	plink --bfile $base --missing --out ${base}_missing
	sort -k 5 -gr "${base}_missing.lmiss" | head > ${base}_missingness.head
	"""	
}

//Assess SNPs for deviation from Hardy-Weinberg equilibrium, removing individuals with deviation significance >0.00001
process HWEfilt{
        publishDir params.output_dir
        input:
        set file(bed), file(bim), file(fam) from qc2

        output:
        tuple file("${base}.hw_dropped.bed"),file("${base}.hw_dropped.bim"),file("${base}.hw_dropped.fam") into qc3
        file("${base}_removed_hwe.bim")
        script:
        base = bed.baseName
        """
        plink --bfile $base --hardy --out ${base}.hw_p_values
        plink --bfile $base --hwe 0.00001 --make-bed --out ${base}.hw_dropped
        plink --bfile $base --exclude "${base}.hw_dropped.bim" --make-bed --out ${base}_removed_hwe
        """

}

//pruning LD
process pruningLD {
        publishDir params.output_dir
	input:
        set file(bed), file(bim), file(fam) from qc3

        output:
        tuple file("${base}.IBD_cleaned.bed"), file("${base}.IBD_cleaned.bim"), file("${base}.IBD_cleaned.fam") into qc4
	tuple file("${base}.LD_IBD.bed"), file("${base}.LD_IBD.bim"), file("${base}.LD_IBD.fam") into qc4_PCA	
	file("${base}.IBD_outliers.txt")	
	shell:
	base = bed.baseName
        '''
	plink --bfile !{base} --indep-pairwise 1500 150 0.2 --out !{base}.LD_one
	plink --bfile !{base} --extract !{base}.LD_one.prune.in --make-bed --out !{base}.LD_two
	awk -f "!{bin}highLDregions4bim_b37.awk" !{base}.LD_two.bim > !{base}.highLDexcludes
	awk '($1<1)|| ($1>22) {print $2}' !{base}.LD_two.bim > !{base}.autosomeexcludes
	cat !{base}.highLDexcludes !{base}.autosomeexcludes > !{base}.highLD_and_autosomal_excludes
	plink --bfile !{base}.LD_two --exclude !{base}.highLD_and_autosomal_excludes --make-bed --out !{base}.LD_three
	plink --bfile !{base}.LD_three --genome --make-bed --out !{base}.IBD
	awk '$10 >= 0.1875 {print $1, $2}' !{base}.IBD.genome > !{base}.IBD_outliers.txt
	plink --bfile !{base}.IBD --remove !{base}.IBD_outliers.txt --make-bed --out !{base}.no_close_relatives
	R --file="!{bin}IndividualIBD.R" --args "!{base}" 1
	plink --bfile !{base}.LD_three -remove !{base}.IBD_INDIV_outliers.txt --make-bed --out !{base}.LD_IBD
	plink --bfile !{base} --remove !{base}.IBD_outliers.txt --make-bed --out !{base}.IBD_cleaned
	'''     

}

process PCA {
        publishDir params.output_dir
	input:
        set file(bed), file(bim), file(fam) from qc4
        file(plinks) from qc4_PCA

        output:
        tuple file("${LD_IBD}.LD_pop_strat.bed"), file("${LD_IBD}.LD_pop_strat.bim"), file("${LD_IBD}.LD_pop_strat.fam") into (qc5_PCA)
        tuple file("${IBD_cleaned}.pop_strat.bed"), file("${IBD_cleaned}.pop_strat.bim"), file("${IBD_cleaned}.pop_strat.fam") into qc5
        file("${LD_IBD}.pop_strat.PC_Output_Associations_FULL.txt")
        file("${LD_IBD}.pop_strat.PC_Output_Associations_SHORT.txt")
	file("${LD_IBD}.pop_strat_PC1_PC2.pdf")
        file("${LD_IBD}.pop_strat_outliers_PC1_PC2.pdf")
        file("${LD_IBD}.pop_strat_outliers.outliers")
	file("${LD_IBD}.pop_strat_outliers.pca.evec") into outliers_pca_evec
	file("${LD_IBD}.pop_strat.PC_for_covariates.txt") into PC_for_cov
	script:
	IBD_cleaned = bed.baseName
        LD_IBD = plinks[0].baseName
        """
	R --file=${bin}genoDistOrder.R --args ${LD_IBD}
	plink --bfile ${LD_IBD} --exclude ${LD_IBD}.snpoutoforder.txt --make-bed --out ${LD_IBD}


	${bin}PCA.sh ${LD_IBD} ${IBD_cleaned} "${bin}" "${params.pheno}"



        """
}

process heteroziogistyTes {
	publishDir params.output_dir
	input:
        set file(bed), file(bim), file(fam) from qc5
        file(plinks) from qc5_PCA
        output:
        tuple file("${base}.het_cleaned.bed"), file("${base}.het_cleaned.bim"), file("${base}.het_cleaned.fam") into qc6
	tuple file("${LD}.LD_het_cleaned.bed"), file("${LD}.LD_het_cleaned.bim"), file("${LD}.LD_het_cleaned.fam") into qc6_LD
        file("${LD}.for_impute.sample") into sample
	script:
        base = bed.baseName
        LD = plinks[0].baseName
        outliers = "${LD}.het.LD_het_outliers_sample_exclude"
        """
	plink --bfile ${LD} --ibc --out ${LD}.het
        R --file="${bin}IdHets.R" --args ${LD}.het
        plink --bfile ${LD} --remove $outliers --make-bed --out ${LD}.LD_het_cleaned

        plink --bfile ${base} --remove $outliers --make-bed --out ${base}.het_cleaned
	plink --bfile ${LD}.LD_het_cleaned --recode oxford --out ${LD}.for_impute
	
	"""
}

///////////////IMPUTATION///////////////////

params.genticMapDir = "ALL.integrated_phase1_SHAPEIT_16-06-14.nomono"
params.chromosomeSizesFile ="${bin}b37.chrom.sizes"
params.referenceHapsFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz"
params.referenceLegendFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"
params.publishDirPath = "${params.output_dir}/imputation-results"


chromosomesList = 1..22

db_path = file(params.genticMapDir)



//Imputation Starts Here//

process plink {


  input:
  tuple file(bed), file(bim), file(fam) from qc6
  each chromosome from chromosomesList

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into plinkOutChan
  script:
  base = bed.baseName
  """
  plink --bfile ${base} --chr $chromosome --make-bed --out chr${chromosome}
  plink -bfile chr${chromosome} --list-duplicate-vars ids-only suppress-first
  [[ -e "plink.dupvar" ]] && plink --bfile chr${chromosome} --exclude plink.dupvar --make-bed --out chr${chromosome}
  """

}

process shapeitCheck {
  validExitStatus 0,1,2
  errorStrategy 'ignore'

  input:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from plinkOutChan
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into shapitCheckChan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )

  """
  shapeit -check --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam --input-ref $hapFile $legendFile $sampleFile --output-log chr${chromosome}.alignments
  """

}

process shapeit {
  errorStrategy 'ignore'
  memory '10 GB'
  input:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from shapitCheckChan
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample") into shapeitChan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )
  excludeFile = "chr${chromosome}.alignments.snp.strand.exclude"
  mapFile = file( db_path.name + "/" + sprintf(params.referenceGeneticMapPattern, chromosome) )

  """
  shapeit --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam --input-ref $hapFile $legendFile $sampleFile --exclude-snp $excludeFile --input-map $mapFile -O chr${chromosome}.phased
  """

}


imputeChromChunckChannel = shapeitChan.flatMap { chromosome, hapsFile, sampleFile ->
   def results = []

   def chunks = getChromosomeChunkPairs(getChromosomeSize(file(params.chromosomeSizesFile), chromosome))
  
   chunks.each { chunkStart, chunkEnd ->
     results.push( [ chromosome, hapsFile, sampleFile, chunkStart, chunkEnd] )
   }

   return results
}

process impute2 {

  errorStrategy 'ignore'
  memory '10 GB'
  input:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample"), val(chunkStart), val(chunkEnd) from imputeChromChunckChannel
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed") into impute2Chan
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed_info") into impute2Chan2
  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )
  mapFile = file( db_path.name + "/" + sprintf(params.referenceGeneticMapPattern, chromosome) )

  """
  impute2 -filt_rules_l 'eur.maf==0' 'type!=SNP' -use_prephased_g -known_haps_g chr${chromosome}.phased.haps -h $hapFile -l $legendFile -m $mapFile -int $chunkStart $chunkEnd -Ne 20000 -o chr${chromosome}-${chunkStart}-${chunkEnd}.imputed

  #sometimes there are no SNP's in a region
  if [ ! -f "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed" ]; then
    touch "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed";
  fi

  """
}



impute2List = impute2Chan.toSortedList() //gives a dataFlow instance, nee to get the val property of it
impute2List = impute2List.val

impute2Map = [:]

impute2List.each { chrom, file ->

  if ( !impute2Map.containsKey(chrom) ) {
   impute2Map.put(chrom, [])
  }
  impute2Map.get(chrom).add(file)

}

impute2MapChannel = Channel.create()

impute2Map.each { chrom, fileList ->
  impute2MapChannel.bind([chrom, fileList])
}

impute2MapChannel.close()


impute2List2 = impute2Chan2.toSortedList() //gives a dataFlow instance, nee to get the val property of it
impute2List2 = impute2List2.val

impute2Map2 = [:]

impute2List2.each { chrom, file ->

  if ( !impute2Map2.containsKey(chrom) ) {
   impute2Map2.put(chrom, [])
  }
  impute2Map2.get(chrom).add(file)

}

impute2MapChannel2 = Channel.create()

impute2Map2.each { chrom, fileList ->
  impute2MapChannel2.bind([chrom, fileList])
}

impute2MapChannel2.close()

input=params.input

process impute2Concat {

  scratch true
  errorStrategy 'ignore'
  input:
  set val(chromosome), file(imputedFiles) from impute2MapChannel
  output:
  file("${input}.temp${chromosome}.impute") into addedChrNum
  shell:
  '''
  cat !{imputedFiles} > chr!{chromosome}.imputed
  awk -v "chr=!{chromosome}" '{ gsub("---",chr,$1); print $0 }' chr!{chromosome}.imputed > "!{input}.temp!{chromosome}.impute"
  '''
}

process impute2Concat2 {

  scratch true
  errorStrategy 'ignore'
  input:
  set val(chromosome), file(imputedFiles_info) from impute2MapChannel2
  output:
  file("${input}.temp${chromosome}.impute_info") into addedChrNum_info
  shell:
  '''
  cat !{imputedFiles_info} > chr!{chromosome}.imputed_info
  awk -v "chr=!{chromosome}" '{ gsub("---",chr,$1); print $0 }' chr!{chromosome}.imputed_info > "!{input}.temp!{chromosome}.impute_info"
  '''
}



////////POST IMPUTATION//////////



process mergeChromosomes{
  publishDir params.output_dir
  input:
  file("${input}.temp*.impute") from addedChrNum.toList()
  file("${input}.temp*.impute_info") from addedChrNum_info.toList()
  output:
  file "${input}.whole.genome" into genome
  file "${input}.whole.genome_info" into genome_info
  script:
  """
  cat ${input}.temp*.impute > "${input}.whole.genome"
  cat ${input}.temp*.impute_info > "${input}.whole.genome_info"

  """
}

process filterByInfo {
  publishDir params.output_dir
  input:
  file "${input}.whole.genome_info" from genome_info
  file "${input}.whole.genome" from genome
  file("*sample") from sample
  output:
  tuple file("${input}.whole.post_imputation.bim"), file("${input}.whole.post_imputation.bed"), file("${input}.whole.post_imputation.fam") into plink_post_imp
  shell:
  '''
  awk '$7 >= 0.8' "!{input}.whole.genome_info" > "!{input}.whole.genome_info.filtered"
  awk 'FNR==NR { a[$2]; next } $2 in a' "!{input}.whole.genome_info.filtered" "!{input}.whole.genome" > "!{input}.whole.genome.filteredByInfo"
  plink --gen "!{input}.whole.genome.filteredByInfo" --sample *sample --hard-call-threshold 0.1 --make-bed --out "!{input}.whole.post_imputation"

  '''
}


process postImpQC {
  publishDir params.output_dir
  input:
  tuple file("${input}.whole.post_imputation.bim"), file("${input}.whole.post_imputation.bed"), file("${input}.whole.post_imputation.fam") from plink_post_imp
  output:
  tuple file("${input}.whole.post_imp.ChrPosA1A2.bim"), file("${input}.whole.post_imp.ChrPosA1A2.bed"),file("${input}.whole.post_imp.ChrPosA1A2.fam") into (chrposa1a2, chrposa1a2_forPRS)
  shell:
  '''
  plink --bfile "!{input}.whole.post_imputation" --maf 0.01 --make-bed --out "!{input}.whole.post_imp.common"
  plink --bfile "!{input}.whole.post_imp.common" --geno 0.99 --make-bed --out "!{input}.whole.post_imp.common.filtered"
  awk 'BEGIN {OFS="\t"} {$2=$1 ":" $4 ":" $5 ":" $6; print $0}' "!{input}.whole.post_imp.common.filtered.bim" > "!{input}.whole.post_imp.ChrPosA1A2.bim"
  cp "!{input}.whole.post_imp.common.filtered.bed" "!{input}.whole.post_imp.ChrPosA1A2.bed"
  cp "!{input}.whole.post_imp.common.filtered.fam" "!{input}.whole.post_imp.ChrPosA1A2.fam"

  '''
}


process associationTest {
  publishDir "${params.output_dir}/imputation-results"
  input:
  tuple file("${input}.whole.post_imp.ChrPosA1A2.bim"), file("${input}.whole.post_imp.ChrPosA1A2.bed"),file("${input}.whole.post_imp.ChrPosA1A2.fam") from chrposa1a2
  output:
  file ("*assoc*") into GWAS_results
  
  script:
  base = "${input}.whole.post_imp.ChrPosA1A2"
  covar = "${params.workingDir}/input_data/covariates.txt"
  """
  plink --bfile $base --pca 5 --out pca 
  R --file=${bin}merge_cov.R --args $covar pca.eigenvec
  

  plink \
--bfile $base \
--logistic  \
--covar cov_pcs.txt \
--pheno ${params.pheno}\
--hide-covar \
--out "${base}.post_imputation_conc_analysis"

   """
}


process gwasGraphs {
  storeDir params.output_dir
  input:
  file GWAS_results
  output:
  file("*manhattan_plot*")
  file("*qq_plot*")
  script:
  """
  ##First argument will be paired with ".whole.post_imp.ChrPosA1A2.post_imputation_conc_analysis.assoc.logistic" and read in. Second argument can be specified for title of the table.
  R --file=${bin}MPlot.R --args ${params.input}

  """
}

process PRS {
  storeDir params.output_dir
  input:
  tuple file("${input}.whole.post_imp.ChrPosA1A2.bim"), file("${input}.whole.post_imp.ChrPosA1A2.bed"),file("${input}.whole.post_imp.ChrPosA1A2.fam") from chrposa1a2_forPRS
  
  output:
  file("*BARPLOT*")
  file("*HIGH-RES*")
  file("*QUANTILES*")
  file("*best")
  script:
  base = "${input}.whole.post_imp.ChrPosA1A2"
  covar = "${params.workingDir}/input_data/covariates.txt"
  """
  gawk 'BEGIN {srand()} {f = FILENAME (rand() <= 0.5 ? ".base" : ".target"); print > f}' ${base}.fam
  plink --bfile ${base} --keep ${base}.fam.base --make-bed --out ${base}.base
  plink --bfile ${base} --keep ${base}.fam.target --make-bed --out ${base}.target
  plink --bfile ${base}.target --pca 5 --out target
  plink --bfile ${base}.base --pca 5 --out base

  R --file=${bin}merge_cov.R --args $covar base.eigenvec
plink \
--bfile ${base}.base \
--logistic  \
--covar cov_pcs.txt \
--pheno ${params.pheno}\
--hide-covar \
--out ${base}

  Rscript ${params.workingDir}/PRSice.R --prsice ${params.workingDir}/PRSice_linux --base ${base}.assoc.logistic --cov target.eigenvec --target ${base}.target --pheno ${params.pheno} --quantile 10 --out ${base}


  """
}



// ----==== utility methods ====----


def getChromosomeSize( chromosomeSizesFile, chromosome ) {

    def result = 0

    chromosomeSizesFile.splitEachLine("\t") {fields ->

        def genomeId
        def path
        if ( fields[0].trim() == "${chromosome}" ) {
          //println "in if"
          result = fields[1].trim().toInteger()
          return
        }

    }

    result
}


def getChromosomeChunkPairs ( chromosomeSize, chunkSize=5000000 ) {

  def result = []
  def numberOfChunks = chromosomeSize / chunkSize
  def remainder = chromosomeSize % chunkSize

  1.upto(numberOfChunks) {
    endPosition = it * chunkSize
    startPosition = (endPosition - chunkSize) + 1
    result = result + [[startPosition, endPosition ]]
  }

  if ( remainder > 0 ) {
    result = result + [[endPosition + 1 , endPosition + remainder ]]
  }

  result
}








