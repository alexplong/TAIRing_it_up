#!/usr/bin/bash

# This example is utilizing samtools to perform the SNP calling
# Our script here involves only looking at Chromosome 4
SAMPLEFILE = Source.info
REFERENCE = TAIR10

# Make this an array job, or just do the first file by default
LINE=$PBS_ARRAYID
if [ ! $LINE ]; then
LINE=$1
fi

if [ ! $LINE ]; then
LINE=1
fi

# Deal with the multiple files
ROW=`head -n $LINE $SAMPLEFILE | tail -n 1`
SAMPLE=`echo "$ROW" | awk '{print $1}'`
URL=`echo "$ROW" | awk '{print $2}'`

# If chromosome 4 is not yet isolated from the reference
if [ ! -f $REFERENCE.Chr.fa ]; then
  module load samtools
  samtools faidx $REFERENCE.fa Chr4 > $REFERENCE.Chr4.fa # pull only Chr4 sequence from Fasta file
fi

# Download BAM files if they haven't been downloaded yet
# Download files from the Source file if not already available
# Followed by isolation of chromosome 4, then sorting and indexing
if [ ! -f $SAMPLE.bam ]; then
  module load samtools
  module load bcftools
  curl $URL > $SAMPLE.bam
  samtools view $SAMPLE.bam Chr4 > $SAMPLE.Chr4.bam # Isolate Chr4 for analysis
  samtools sort $SAMPLE.Chr4.bam $SAMPLE.Chr4.sorted.bam #sort BAM file
  samtools index $SAMPLE.Chr4.sorted.bam $SAMPLE.Chr4.sorted.bai # index bam
  samtools mpileup -u -f $REFERENCEfa $SAMPLE.sorted.bam > $SAMPLE.Chr4.bcf # Call varients
  bcftools view $SAMPLE.Chr4.bcf > $SAMPLE.Chr4.vcf # Compile a readable VCF, including reference files
  bcftools view -v snps $SAMPLE.Chr4.bcf > $SAMPLE.Chr4.vcf #Compile a readable VCF, displaying only SNPs
fi
