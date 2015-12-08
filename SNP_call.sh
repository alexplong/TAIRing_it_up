#!/usr/bin/bash

########
# This is the initial segment of code to be run for our project.
# This code assumes that you want only to analyze chromosome 4 (picked because it had the smallest file size) #ommited
# IMPORTANT: the reference genome file that was downloaded from TAIR had to be modified
# Specifically the chromosomes had to be named "Chr1" instead of "1" to match the annotation in the BAM files
# Be aware that we also had to modify the order of the chloroplast and mitochondrial genomes, to match the order in the BAM files
# We did this manually using a texteditor before the files were read into the program here.
# As a final note, the last line with HaplotypeCaller specifies an amount of memory for java, this can't be called from the
#	environment, so we hard coded it as 32gb

# Sample file provides 2-column tab delimited file containing:
# Sample name and URL
SAMPLEFILE=Source.info

# Reference Genome
# Modify to match annotation in BAM files

REFERENCE=TAIR10

GATK=/opt/GATK/3.4.0/GenomeAnalysisTK.jar
PICARD=/opt/picard/1.81


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

# Index the reference
if [ ! -f $REFERENCE.dict ]; then
  module load samtools
  java -jar $PICARD/CreateSequenceDictionary.jar R=$REFERENCE.fa O=$REFERENCE.dict
  samtools faidx $REFERENCE.fa
fi

# Download files from the Source file if not already available
if [ ! -f $SAMPLE.bam ]; then
  curl $URL > $SAMPLE.bam

# Index the BAM file

  if [ ! -f $SAMPLE.bai ]; then
    java -jar $PICARD/BuildBamIndex.jar I=$SAMPLE.bam
  fi
fi

# Call SNP variants

if [ ! -f $SAMPLE.g.vcf ]; then
  java -Xmx32g -jar $GATK -T HaplotypeCaller -R $REFERENCE.fa -I $SAMPLE.bam -ERC GVCF -o $SAMPLE.g.vcf -nct $PBS_NUM_PPN -L
fi
