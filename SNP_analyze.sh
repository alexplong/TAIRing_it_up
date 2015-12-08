#!/usr/bin/bash

# This is the second segment of code necessary to run our project's analysis
# This script should not be run until the first code file has been completely run as it depends on output from that script

# Combine  output files and perform a joint genotyping using GATK

REFERENCE=TAIR10

GATK=/opt/GATK/3.4.0/GenomeAnalysisTK.jar
snpeff=/opt/snpeff/3.4/snpEff.jar
SNPdat=


# Get the g.vcf files and add “--variant” between them
N=`echo ls *.g.vcf | sed 's/ / --variant /g' | sed 's/ls //'`

#Combine all of the gvcf files into one
if [ ! -f All.g.vcf ]; then
  java -Xmx32g -jar $GATK -T CombineGVCFs -R $REFERENCE.fa $N -o All.g.vcf
fi

#do joint genotyping of that file
if [ ! -f All.vcf ]; then
  java -Xmx32g -jar $GATK -T GenotypeGVCFs -R $REFERENCE.fa --variant All.g.vcf -o All.vcf -nt $PBS_NUM_PPN
fi

# Create a more readable .tab file and fasta alignment of the variants
# Run script to identify unique and shared SNPs
# Note that private SNPs script will require heavy modification if used for anything else
# or if more samples are added
if [ ! -f All.tab ]; then
  module load vcftools
  vcf-to-tab  < All.vcf > All.tab
  perl vcf_tab_to_fasta_alignmentv1.pl -i All.tab
  python Private_SNPs.py > private_SNPs.tab
fi

# Get TAIR gff files
if [ ! -f TAIR10_GFF3genes.gff ]; then
	curl ftp://ftp.arabidopsis.org/home/tair//Maps/gbrowse_data/TAIR10/TAIR10_GFF3_genes.gff
fi

# Convert GFF to GTF. Necesarry for SNPdat
if [ ! -f TAIR10.gtf ]; then
  python GFF2GTF.py TAIR10_GFF3genes.gff > TAIR10.gtf
fi

#I cant figure out how to use the snpEff that's on the cluster, so I put the .jar, .config, and reference genomes files in the SNPCalling folder
# snpEff
java -Xmx32g -jar -snpEff.jar athalianaTair10 -v  ALL.vcf  -t $PBS_NUM_PPN  > ALL.snpEff.vcf


#Use snpsift to extract the missense, nonsense and silent variants
cat -Xmx32g ALL.snpEff.vcf | ./vcfEffOnePerLine.pl |java -jar SnpSift.jar  extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "GEN[Bur_0].GT" "GEN[Can_0].GT" "GEN[Ler_0].GT" "GEN[Tsu_0].GT" | grep -v "stream" |grep -v "intergenic" |uniq > ALL.snpSift.txt

#An example of how to manually count the number of missense and silent variants in the Bur accession
burmv=$(cat ALL.snpSift.txt | cut -f1,2,5,6 | awk '$4!="0/0"' | awk '$3=="missense_variant"' | uniq | wc -l)
bursv=$(cat ALL.snpSift.txt | cut -f1,2,5,6 | awk '$4!="0/0"' | awk '$3=="synonymous_variant"' | uniq | wc -l)

# Run SNPdat to identify (non)synonymous mutations
# Requires REF & GTF file
# Please not that SNPdat is extremely slow
# It's recommended that the vcf output file be filtered first to regions of interest before running SNPdat
if [ ! -f ALL.vcf.output ]; then
  perl SNPdat_v1.0.5.pl -i ALL.vcf -f $REFERENCE.fa -g TAIR10.GTF
fi

# blah test git
