#!/bin/bash

####aliases
alias bwa="/Users/bencekover/bioinformatics/bwa/bwa"
alias samtools="/Users/bencekover/bioinformatics/samtools-1.16.1/samtools"
alias gatk="/Users/bencekover/bioinformatics/gatk-4.3.0.0/gatk" 
alias picard="java -jar /Users/bencekover/bioinformatics/picard/picard.jar"
alias bcftools="/Users/bencekover/bioinformatics/bcftools/bcftools.h"

input1=/Users/bencekover/bioinformatics/projects/segregants/R4_raw_snps.vcf
input2=/Users/bencekover/bioinformatics/projects/segregants/R34_raw_snps.vcf
input3=/Users/bencekover/bioinformatics/projects/segregants/R44_raw_snps.vcf
input4=/Users/bencekover/bioinformatics/projects/segregants/R45_raw_snps.vcf
input5=/Users/bencekover/bioinformatics/projects/segregants/R48_raw_snps.vcf
input6=/Users/bencekover/bioinformatics/projects/segregants/R49_raw_snps.vcf
input7=/Users/bencekover/bioinformatics/projects/segregants/R50_raw_snps.vcf
input8=/Users/bencekover/bioinformatics/projects/segregants/R51_raw_snps.vcf
input9=/Users/bencekover/bioinformatics/projects/segregants/R52_raw_snps.vcf
input10=/Users/bencekover/bioinformatics/projects/segregants/R53_raw_snps.vcf
input11=/Users/bencekover/bioinformatics/projects/segregants/R54_raw_snps.vcf
input12=/Users/bencekover/bioinformatics/projects/snp_jb50/results_raw.vcf
input13=/Users/bencekover/bioinformatics/projects/snp_jb759/results_raw.vcf


ref=/Users/bencekover/bioinformatics/projects/snp_jb759/Schizosaccharomyces_pombe_all_chromosomes.fa





 gatk CombineGVCFs \
   -R $ref\
   --variant $input1 \
   --variant $input2 \
    --variant $input3 \
    --variant $input4 \
    --variant $input5 \
    --variant $input6 \
    --variant $input7 \
    --variant $input8 \
    --variant $input9 \
    --variant $input10 \
    --variant $input11 \
    --variant $input12 \
    --variant $input13 \
   -O cohort.g.vcf.gz


  gatk GenotypeGVCFs \
    -R $ref \
    -V cohort.g.vcf.gz \
    -O cohort.vcf.gz  \
    -ploidy 1