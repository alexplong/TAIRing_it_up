# TAIRing_it_up
GEN220 SNP Group Project

Team members: Alex Rajewski, Alex Plong, Joseph Carrillo, Chrissy Dodge

Our project revolves around 2 main scripts and should be run in order.
  1. SNP_call should be run first, followed by
  2. SNP_analyze

The header for SNP_call contains these lines  
```
SAMPLEFILE=Source.info
REFERENCE=TAIR10
GATK=/opt/GATK/3.4.0/GenomeAnalysisTK.jar
PICARD=/opt/picard/1.81
```

The Source.info file is a 2-column tab-delimited file, with the sample name and URL for its associated BAM file.
The reference file should look like the example below. Ensure that the names of the chromosomes associated in the BAM file also match your reference. This check is only necessary if you did not align the sequences yourself and downloaded already aligned BAM files from a study.

For example, for our run it looks like this:
```
## Source.info
bur_0	ftp://ftp.sra.ebi.ac.uk/vol1/ERA023/ERA023479/bam/Bur_0_bur_PII.bam
Can_0	ftp://ftp.sra.ebi.ac.uk/vol1/ERA023/ERA023479/bam/Can_0_can_PII.bam
Ler_0	ftp://ftp.sra.ebi.ac.uk/vol1/ERA023/ERA023479/bam/Ler_0_ler_PII.bam
Tsu_0	ftp://ftp.sra.ebi.ac.uk/vol1/ERA023/ERA023479/bam/Tsu_0_tsu_PII.bam

## TAIR10.fa
>Chr1
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
GAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT
TGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTA
CATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTA
```
Some additional notes
When running GATK's HaplotypeCaller, we hard coded the amount of memory the program should utilize (for our case, 32gb). Please adjust accordingly
```
java -Xmx32g -jar $GATK -T HaplotypeCaller -R $REFERENCE.fa -I $SAMPLE.bam -ERC GVCF -o $SAMPLE.g.vcf -nct $PBS_NUM_PPN -L
```
