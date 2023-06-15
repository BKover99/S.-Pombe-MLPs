#!/bin/bash

####aliases
alias bwa="/Users/bencekover/bioinformatics/bwa/bwa"
alias samtools="/Users/bencekover/bioinformatics/samtools-1.16.1/samtools"
alias gatk="/Users/bencekover/bioinformatics/gatk-4.3.0.0/gatk" 
alias picard="java -jar /Users/bencekover/bioinformatics/picard/picard.jar"
alias bcftools="/Users/bencekover/bioinformatics/bcftools/bcftools.h"


for sample in R4 R34 R44 R45 R48 R49 R50 R51 R52 R53 R54
do
  echo "Processing sample: $sample"
  
  cd /Users/bencekover/bioinformatics/projects/segregants/40-798766616/00_fastq/
  input1=${sample}_R1_001.fastq.gz
  input2=${sample}_R2_001.fastq.gz
  output1=${sample}_R1_001_trimmed.fastq.gz
  output2=${sample}_R2_001_trimmed.fastq.gz

  cutadapt -a AGATCGGAAGAGC\
   -A AGATCGGAAGAGC\
   -o ${output1} -p ${output2}\
   ${input1} ${input2}
  
  #### setting working directory and assigning file paths to variables
  
  ref="/Users/bencekover/bioinformatics/projects/snp_jb759/Schizosaccharomyces_pombe_all_chromosomes.fa"
  seq1=${output1}
  seq2=${output2}
  variation="/Users/bencekover/bioinformatics/projects/snp_jb759/Schizosaccharomycespombesnps.20121212.filt3c.all.snpsonly.anno4.vcf"

  
  #####building reference index
  bwa index ${ref}

  #####creating alignment
  bwa mem ${ref} ${seq1} ${seq2} > alignment.sam

  rm ${seq1}
  rm ${seq2}

  #####sorting
  samtools sort alignment.sam -o sorted_alignment.bam
  rm alignment.sam

  #####assigning to one read group
  picard AddOrReplaceReadGroups \
         I=sorted_alignment.bam \
         O=sorted_alignment_single_rg.bam \
         RGID=4 \
         RGLB=lib1 \
         RGPL=ILLUMINA \
         RGPU=unit1 \
         RGSM=20
  rm sorted_alignment.bam

  #####GATK mark duplicates
  gatk MarkDuplicates \
        I=sorted_alignment_single_rg.bam \
        O=marked_sorted_alignment_single_rg.bam \
        M=marked_dup_metrics.txt

  rm sorted_alignment_single_rg.bam

  #####GATK recalibrating base quality scores
  Gate IndexFeatureFile -I ${variation}

#####GATK recalibrating base quality scores (cont'd)
gatk BaseRecalibrator \
-R ${ref} \
-I marked_sorted_alignment_single_rg.bam \
--known-sites ${variation} \
-O recal_data.table

#####Apply base quality score recalibration
gatk ApplyBQSR \
-R ${ref} \
-I marked_sorted_alignment_single_rg.bam \
--bqsr-recal-file recal_data.table \
-O recal_reads.bam 

rm marked_sorted_alignment_single_rg.bam

#####index the bam file
samtools index recal_reads.bam

#####call SNPs
gatk HaplotypeCaller \
-R ${ref} \
-I recal_reads.bam \
-O ${sample}_raw_snps.vcf \
-ERC GVCF \
--max-reads-per-alignment-start 0 \
--native-pair-hmm-threads 4 


rm recal_reads.bam

done

echo "Pipeline complete"



