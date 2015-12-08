#!usr/bin/perl -w
use strict;
use bytes;

my $version = "SNPdat v1.0.5";

print "\n$version\n\n";

print "start time:\n";
system("date");

#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------
#
#	Section 1 - Specify data and check for inconsistencies with input/define global amino acid table
#
#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------

#Print flash page of all arguments that can be supplied for SNPdat if no arguments are supplied
if(@ARGV == 0){
	print "\n\nError no arguments supplied.\n\n";
	
	&help(0);
	
	exit;
}

#check options for the help menu
for(my $i = 0; $i < @ARGV; $i++){
	if($ARGV[$i] eq '-h' or $ARGV[$i] eq '--help'){
		
		&help(4);
		exit;
		
	}
	
	if($ARGV[$i] eq '-v' or $ARGV[$i] eq '--version'){
		print "\nThis version of SNPdat is\n$version\n\n";
		exit;
	}
}

#to check if flags that SNPdat will not recognise
my %elements = (
	'-i' => 'input',
	'-g' => 'GTF',
	'-f' => 'FASTA',
	'-d' => 'dbSNP',
	'-o' => 'output',
	'-s' => 'summary',
	'-x' => 'invoke feature bounday crossing',
	'-h' => 'HELP',		#these are not necessary here as if they are provided the program will exit before here anyway
	'--help' => 'HELP',

);



#check for flags that SNPdat doesnt recognise
for(my $i = 0; $i < @ARGV; $i++){
	if(exists $elements{$ARGV[$i]}){
		$i++;
		next;
	}else{
		print "\nWarning $ARGV[$i] not an appropriate flag!\n\n";
		
		&help(0);
	}
}

#	Variable for storing information for the Summary output file
#
my $summary_gtf_feat = 0;			# Total number of features in the GTF - the number of lines
my %summary_features;				# The available features in the GTF
my %summary_chrs_in_GTF;			# Used to get total number of chromosomes in the GTF file

my $summary_user_snps = 0;			# Total number of queried SNPs - the number of input lines
my %summary_user_chrs;				# The queried chromosomes
my $summary_annotated = 0;			# Number of SNPs annotated
my $summary_snp_skipped = 0;		# The number of skipped SNPs - no GTF or incorrect input
my $summary_in_feature_snps = 0;	# Number of annotated SNPs that were in a feature
my $summary_nonsyn_annotations = 0;	# Number of non-synonymous annotations
my $summary_noGTF_snp = 0;			# Number of SNPs with no GTF information
my $summary_nofasta_snps = 0;		# Number of SNPs with no FASTA information



#globally defined hash table of the amino acid Single Letter Codes (SLC) and their corresponding DNA codes
my %SLC = (
	'TCA' => 'S', # Serine
	'TCC' => 'S', # Serine
	'TCG' => 'S', # Serine
	'TCT' => 'S', # Serine
	'TTC' => 'F', # Phenylalanine
	'TTT' => 'F', # Phenylalanine
	'TTA' => 'L', # Leucine
	'TTG' => 'L', # Leucine
	'TAC' => 'Y', # Tyrosine
	'TAT' => 'Y', # Tyrosine
	'TAA' => '-', # Stop
	'TAG' => '-', # Stop
	'TGC' => 'C', # Cysteine
	'TGT' => 'C', # Cysteine
	'TGA' => '-', # Stop
	'TGG' => 'W', # Tryptophan
	'CTA' => 'L', # Leucine
	'CTC' => 'L', # Leucine
	'CTG' => 'L', # Leucine
	'CTT' => 'L', # Leucine
	'CCA' => 'P', # Proline
	'CCC' => 'P', # Proline
	'CCG' => 'P', # Proline
	'CCT' => 'P', # Proline
	'CAC' => 'H', # Histidine
	'CAT' => 'H', # Histidine
	'CAA' => 'Q', # Glutamine
	'CAG' => 'Q', # Glutamine
	'CGA' => 'R', # Arginine
	'CGC' => 'R', # Arginine
	'CGG' => 'R', # Arginine
	'CGT' => 'R', # Arginine
	'ATA' => 'I', # Isoleucine
	'ATC' => 'I', # Isoleucine
	'ATT' => 'I', # Isoleucine
	'ATG' => 'M', # Methionine
	'ACA' => 'T', # Threonine
	'ACC' => 'T', # Threonine
	'ACG' => 'T', # Threonine
	'ACT' => 'T', # Threonine
	'AAC' => 'N', # Asparagine
	'AAT' => 'N', # Asparagine
	'AAA' => 'K', # Lysine
	'AAG' => 'K', # Lysine
	'AGC' => 'S', # Serine
	'AGT' => 'S', # Serine
	'AGA' => 'R', # Arginine
	'AGG' => 'R', # Arginine
	'GTA' => 'V', # Valine
	'GTC' => 'V', # Valine
	'GTG' => 'V', # Valine
	'GTT' => 'V', # Valine
	'GCA' => 'A', # Alanine
	'GCC' => 'A', # Alanine
	'GCG' => 'A', # Alanine
	'GCT' => 'A', # Alanine
	'GAC' => 'D', # Aspartic Acid
	'GAT' => 'D', # Aspartic Acid
	'GAA' => 'E', # Glutamic Acid
	'GAG' => 'E', # Glutamic Acid
	'GGA' => 'G', # Glycine
	'GGC' => 'G', # Glycine
	'GGG' => 'G', # Glycine
	'GGT' => 'G', # Glycine
	'A??' => 'X', # Unknown
	'T??' => 'X', # Unknown
	'G??' => 'X', # Unknown
	'C??' => 'X', # Unknown
	'?A?' => 'X', # Unknown
	'?T?' => 'X', # Unknown
	'?G?' => 'X', # Unknown
	'?C?' => 'X', # Unknown
	'??A' => 'X', # Unknown
	'??T' => 'X', # Unknown
	'??G' => 'X', # Unknown
	'??C' => 'X', # Unknown
	'AA?' => 'X', # Unknown
	'AT?' => 'X', # Unknown
	'AC?' => 'X', # Unknown
	'AG?' => 'X', # Unknown
	'TA?' => 'X', # Unknown
	'TT?' => 'X', # Unknown
	'TC?' => 'X', # Unknown
	'TG?' => 'X', # Unknown
	'CA?' => 'X', # Unknown
	'CT?' => 'X', # Unknown
	'CC?' => 'X', # Unknown
	'CG?' => 'X', # Unknown
	'GA?' => 'X', # Unknown
	'GT?' => 'X', # Unknown
	'GC?' => 'X', # Unknown
	'GG?' => 'X', # Unknown
	'A?A' => 'X', # Unknown
	'A?T' => 'X', # Unknown
	'A?C' => 'X', # Unknown
	'A?G' => 'X', # Unknown
	'T?A' => 'X', # Unknown
	'T?T' => 'X', # Unknown
	'T?C' => 'X', # Unknown
	'T?G' => 'X', # Unknown
	'C?A' => 'X', # Unknown
	'C?T' => 'X', # Unknown
	'C?C' => 'X', # Unknown
	'C?G' => 'X', # Unknown
	'G?A' => 'X', # Unknown
	'G?T' => 'X', # Unknown
	'G?C' => 'X', # Unknown
	'G?G' => 'X', # Unknown
	'?AA' => 'X', # Unknown
	'?AT' => 'X', # Unknown
	'?AC' => 'X', # Unknown
	'?AG' => 'X', # Unknown
	'?TA' => 'X', # Unknown
	'?TT' => 'X', # Unknown
	'?TC' => 'X', # Unknown
	'?TG' => 'X', # Unknown
	'?CA' => 'X', # Unknown
	'?CT' => 'X', # Unknown
	'?CC' => 'X', # Unknown
	'?CG' => 'X', # Unknown
	'?GA' => 'X', # Unknown
	'?GT' => 'X', # Unknown
	'?GC' => 'X', # Unknown
	'?GG' => 'X', # Unknown
	'???' => 'X', # Unknown
	
);


#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------
#
#	Section 2 - Check for correct number of arguments and open output file
#
#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------

#
#	This will be sample code to check to see if flags have been given for the input files
#	currently still need error handling for command line arguments not contained in snpdat script
#	also need to error handle when multiple switches are given with no arguments
#

#file input check variables - 0 for not set; 1 for set; 2+ error for assigned to many
my $user_input = 0;		# -i
my $gtf_file = 0;		# -g
my $fasta_file = 0;		# -f
my $output_file = 0;	# -o
my $db_snp_file = 0;	# -d
my $summary_output = 0;	# -s
my $exon_only = 0;		# -x


#hash containing the command line switches and arguments (files)
my %snpdat_files;
my @x_options;		#these are for the -x option only. These is the option given with -x
my %x_features;		#These are the features specified by the user with the -x option

#traverse the command line arguments and assign
for(my $i = 0; $i <@ARGV; $i++){	#$i+=2 should be better than ++ because this way we can always ensure that a file follows the switch - it is still possible to confuse the script
	
	
	#switch statement to assign command line switches to files
	SWITCH:{
		$ARGV[$i] eq '-i' && do { $user_input++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-g' && do { $gtf_file++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-f' && do { $fasta_file++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-o' && do { $output_file++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-d' && do { $db_snp_file++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-s' && do { $summary_output++; $snpdat_files{$ARGV[$i]} = $ARGV[$i+1]; $i++ };
		$ARGV[$i] eq '-x' && do { $exon_only++; @x_options = split(",", $ARGV[$i+1]); $i++};
		
	}
	
	#maybe should use a set of if-elsif-else statements - this will also allow an error handling for ambigous switches

}


for(my $i = 0; $i < @x_options; $i++){
	$x_features{$x_options[$i]} = $x_options[$i];
}

#switch statement to test that user command line correct or read correctly
#need to test this code a bit more to ensure all of the switches are tested
SWITCH:{
	
	$user_input != 1 && do { die &help(1); };
	$gtf_file != 1 && do { die &help(2); };
	$fasta_file != 1 && do { die &help(3); };
	$output_file != 1 && do { $snpdat_files{-o} = $snpdat_files{-i}.".output"; print "NOTE - output file not specified - OUTPUT file is now \'$snpdat_files{-o}\'\n"; };
	$db_snp_file != 1 && do {print "NOTE - No dbSNP processed FLAT file given\n";};
	$summary_output != 1 && do { $snpdat_files{-s} = $snpdat_files{-i}.".summary"; print "NOTE - summary file not specified - summary file is now \'$snpdat_files{-s}\'\n"; };
	$exon_only > 0 && do {print "NOTE - You have chosen to invoke the feature boundary retrieval option.\nThis is recommended for advanced users only. See Manual/website for more information\n";};

}

#open the Summary output file
open(SUMMARY, ">$snpdat_files{-s}") or die "Error cannot open the summary file for writing\n";
#print this version to the summary file
print SUMMARY "The is $version\n";
#print the version of SNPdat that is in use

print SUMMARY <<END;
This is a summary result file generated by SNPdat

The following options were supplied to SNPdat:

END

print SUMMARY "Input file (-i):      \t\"$snpdat_files{-i}\"\n";
print SUMMARY "Output file (-o):     \t\"$snpdat_files{-o}\"\n";
print SUMMARY "GTF file (-g):        \t\"$snpdat_files{-g}\"\n";
print SUMMARY "FASTA file (-f):      \t\"$snpdat_files{-f}\"\n";
print SUMMARY "This summary file(-s):\t\"$snpdat_files{-s}\"\n";

if($exon_only > 0){
	print SUMMARY "NOTE - You have chosen to invoke the feature boundary retrieval option.\nThis is recommended for advanced users only. See Manual/website for more information\n";
	print SUMMARY "Codons crossing the feature boundary for these features only will have information retrieved from the up/downstream feature of the same type:\n";
	
	for(my $i = 0; $i < @x_options; $i++){
		print SUMMARY "$x_options[$i]\n";
	}
	
}

print SUMMARY "\nThis is a brief summary of your results:";

#Summary statistics are printed at the end to the SUMMARY file

#open program output file
open(OUTPUT, ">$snpdat_files{-o}") or die "Error! Could not open output file\n";

#print output file headers
print OUTPUT "Chromosome Number\t";
print OUTPUT "SNP Position\t";
print OUTPUT "Within a Feature\t";
print OUTPUT "Region\t";
print OUTPUT "Distance to nearest feature\t";
print OUTPUT "Feature\t";
print OUTPUT "Number of different features for this SNP\t";
print OUTPUT "Number of annotations for this feature\t";
print OUTPUT "Start of current feature\t";
print OUTPUT "End of current feature\t";
print OUTPUT "gene ID containing the current feature\t";
print OUTPUT "gene name containing the current feature\t";
print OUTPUT "transcript ID containing the current feature\t";
print OUTPUT "transcript name containing the current feature\t";
print OUTPUT "Annotated Exon [Rank/Total]\t";
print OUTPUT "Strand sense\t";
print OUTPUT "Annotated reading Frame\t";
print OUTPUT "Estimated Reading Frame\t";
print OUTPUT "Estimated Number of Stop Codons\t";
print OUTPUT "Codon\t";
print OUTPUT "Amino Acid\t";
print OUTPUT "synonymous\t";
print OUTPUT "Protein ID containing the current feature\t";

if($db_snp_file > 0){
	print OUTPUT "RSidentifier\t";
}

print OUTPUT "Additional Notes\n";

#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------
#
#	Section 3 - Parsing necessary files
#
#--------------------------------------------------
# + + + + + + + + + + + + + + + + + + + + + + + + + 
#--------------------------------------------------

#--------------------------------------------------
#
#	Open GTF File and pass into memory
#
#--------------------------------------------------

#Open the annotation file - GTF format
open(GTF, "$snpdat_files{-g}") or die "Error! No GTF annotation file specified!\n";

#==================================================
#
#	GTF contains 9 fields each seperated by a single tab
#	Field 9 is the attributes field
#	This is a field with descriptive additional information on the feature such as transcript_id, exon_number (not mandatory) 
#	attributes are seperated by a semi-colon and a single space
#
#==================================================

my %seqname;				#must be a chromosome or scaffold (usually a chromosome)
my %source;					#usually the name of the program that generated this feature
my %feature;				#the name of this type of feature eg "CDS", "start_codon","stop_codon", "exon"
my %start;					#Starting position of the feature - bases start at 1
my %end;					#Ending position of the feature (inclusive)
my %score;					#Score between 0-1000 to determine the level of grey in which the feature is displayed. If no score enter "."
my %strand;					#+/- or "." for missing (dont know)
my %frame;					#0-1 to represent the reading frame of the first base (for coding exons only, if not coding then should be set to ".")
my %gene_id;				#gene_id value (mandatory first attribute to field 9)
my %transcript_id;			#transcipt_id value (mandatory second attribute to field 9)
my %exon_number;			#not necessary for script to run but will check to see if that information is contained in the annotation file (Should be in field 9)
my %gene_name;				#not necessary for script to run but will check to see if that information is contained in the annotation file (Should be in field 9)
my %transcript_name;		#not necessary for script to run but will check to see if that information is contained in the annotation file (Should be in field 9)
my %protein_id;				#not necessary for script to run but will check to see if that information is contained in the annotation file (Should be in field 9)

#create a list of unique chromosome names from the GTF
#this will be used to check user input as it is read
my %chr_list;

#create a list with the maximum value of an exon number at a gene_id
my %gene_list;

#keep track of the transcript max and min positions
#used to identify SNPs that are intergenic or genic
my %transcript_max;
my %transcript_min;

##additional data for the start and end points for all transcripts
##This can later be used to get the number of exons for each transcript and the the exon number of each feature
##for the current transcript. Will just need to make sure the feature used to buiild this structure is exon
my %transcript_data;
my %transcript_start;
my %transcript_end;

#These will track transcript specific exon data
my %transcript_exon_min;
my %transcript_exon_max;
my %transcript_exon_start;
my %transcript_exon_end;
my %transcript_exon_data;

#Read into memory the GTF file (annotation file)
while(<GTF>){
	
	$summary_gtf_feat++;

	$_ =~ s/\r|\n//g;					#the same as chomp but also dos carriage returns

	my @temp = split("\t", $_);			#Split the GTF fields
	
	#remove any whitspace
	foreach(@temp){
		$_ =~ s/\s*//;
		
	}


	$summary_features{$temp[2]} = 0;
	$summary_chrs_in_GTF{$temp[0]} = 0;


	$temp[0] =~ s/^.*?\"//;
	$temp[0] =~ s/\".*$//;
	
	$temp[1] =~ s/^.*?\"//;
	$temp[1] =~ s/\".*$//;
	
	$temp[2] =~ s/^.*?\"//;
	$temp[2] =~ s/\".*$//;
	
	$temp[3] =~ s/^.*?\"//;
	$temp[3] =~ s/\".*$//;
	
	$temp[4] =~ s/^.*?\"//;
	$temp[4] =~ s/\".*$//;
	
	$temp[5] =~ s/^.*?\"//;
	$temp[5] =~ s/\".*$//;
	
	$temp[6] =~ s/^.*?\"//;
	$temp[6] =~ s/\".*$//;
	
	$temp[7] =~ s/^.*?\"//;
	$temp[7] =~ s/\".*$//;

	$temp[0] = "\U$temp[0]";	#convert seqname (chromosome name) to upper case - ensure matching with Sequence IDs

	#push first 8 fields of information onto respective arrays
	push(@{$seqname{$temp[0]}}, $temp[0]);
	push(@{$source{$temp[0]}}, $temp[1]);
	push(@{$feature{$temp[0]}}, $temp[2]);
	push(@{$start{$temp[0]}}, $temp[3]);
	push(@{$end{$temp[0]}}, $temp[4]);
	push(@{$score{$temp[0]}}, $temp[5]);
	push(@{$strand{$temp[0]}}, $temp[6]);
	push(@{$frame{$temp[0]}}, $temp[7]);
	
	#split the 9th field (attributes field)
	my @attributes = split("; ", $temp[8]);

	$attributes[0] =~ s/^.*?\"//;
	$attributes[0] =~ s/\".*$//;
	
	$attributes[1] =~ s/^.*?\"//;
	$attributes[1] =~ s/\".*$//;

	#the first 2 attributes in the 9th field are mandatory so push them onto respective arrays
	push(@{$gene_id{$temp[0]}}, $attributes[0]);
	push(@{$transcript_id{$temp[0]}}, $attributes[1]);
	
	my $size = @attributes;

	#these are the only non-mandatory attributes that the script will look for - assumes there is a set name for each type
	#the script will only print these to out file - will not use these for further calculations
	#suggestions for other fields are welcome
	my $exon_test = 0;
	my $transcript_test = 0;
	my $gene_test = 0;
	my $protein_test = 0;

	#figure out if the current line has these non-mandatory fields - if so push onto respective array
	for(my $i = 2; $i < $size; $i++){
		
		if($attributes[$i] =~ m/exon_number/ and $exon_test == 0){
			$attributes[$i] =~ s/^.*?\"//;
			$attributes[$i] =~ s/\".*$//;
			push(@{$exon_number{$temp[0]}}, $attributes[$i]);
			$exon_test = 1;
			
			#keep track of the max exon number for a gene id
			#still need to code a what to do if exon number field is not included
			if(exists $gene_list{$attributes[0]}){
				if($attributes[$i] > $gene_list{$attributes[0]}){
					$gene_list{$attributes[0]} = $attributes[$i];
				}
			}else{
				$gene_list{$attributes[0]} = $attributes[$i];
			}
			
		}
		
		if($attributes[$i] =~ m/gene_name/ and $gene_test == 0){
			$attributes[$i] =~ s/^.*?\"//;
			$attributes[$i] =~ s/\".*$//;
			push(@{$gene_name{$temp[0]}}, $attributes[$i]);
			$gene_test = 1;
			
		}
		
		if($attributes[$i] =~ m/transcript_name/ and $transcript_test == 0){
			$attributes[$i] =~ s/^.*?\"//;
			$attributes[$i] =~ s/\".*$//;
			push(@{$transcript_name{$temp[0]}}, $attributes[$i]);
			$transcript_test = 1;
			
		}
		
		if($attributes[$i] =~ m/protein_id/ and $protein_test == 0){
			$attributes[$i] =~ s/^.*?\"//;
			$attributes[$i] =~ s/\".*$//;
			push(@{$protein_id{$temp[0]}}, $attributes[$i]);
			$protein_test = 1;
			
		}
	}

	#*If the current line does not have any of the non-mandatory fields the script will push 'NA' onto the respective array
	if($exon_test == 0){
		push(@{$exon_number{$temp[0]}}, "NA");
		
	}
	
	if($gene_test == 0){
		push(@{$gene_name{$temp[0]}}, "NA");
		
	}
	
	if($transcript_test == 0){
		push(@{$transcript_name{$temp[0]}}, "NA");
		
	}
	
	if($protein_test == 0){
		push(@{$protein_id{$temp[0]}}, "NA");
		
	}
	
	#create list of unique chromosomes from GTF
	$chr_list{$temp[0]} = 1;	
	
	#keep track of transcript min and max position
	if(exists $transcript_max{$temp[0]}{$attributes[1]}){
		if($temp[4] > $transcript_max{$temp[0]}{$attributes[1]}){
			$transcript_max{$temp[0]}{$attributes[1]} = $temp[4];
		}
		
		if($temp[3] < $transcript_min{$temp[0]}{$attributes[1]}){
			$transcript_min{$temp[0]}{$attributes[1]} = $temp[3];
		}
	}else{
		$transcript_max{$temp[0]}{$attributes[1]} = $temp[4];
		$transcript_min{$temp[0]}{$attributes[1]} = $temp[3];
	}
	
	
	push(@{$transcript_data{$temp[0]}{$attributes[1]}}, $attributes[1]);
	push(@{$transcript_start{$temp[0]}{$attributes[1]}}, $temp[3]);
	push(@{$transcript_end{$temp[0]}{$attributes[1]}}, $temp[4]);
	
	#keep track of transcript information for only the exons!
	#keep track of min and max exon position in a gene transcript
	
	if(exists $x_features{$temp[2]}){
		
		if(exists $transcript_exon_max{$temp[0]}{$attributes[1]}{$temp[2]}){
			
			if($temp[4] > $transcript_exon_max{$temp[0]}{$attributes[1]}{$temp[2]}){
				$transcript_exon_max{$temp[0]}{$attributes[1]}{$temp[2]} = $temp[4];	#max exon value for each transcript
			}
			
			if($temp[3] < $transcript_exon_min{$temp[0]}{$attributes[1]}{$temp[2]}){
				$transcript_exon_min{$temp[0]}{$attributes[1]}{$temp[2]} = $temp[3];	#min exon value for each transcript
			}
			
		}else{
			$transcript_exon_max{$temp[0]}{$attributes[1]}{$temp[2]} = $temp[4];
			$transcript_exon_min{$temp[0]}{$attributes[1]}{$temp[2]} = $temp[3];
		}
		
		push(@{$transcript_exon_start{$temp[0]}{$attributes[1]}{$temp[2]}}, $temp[3]);	#exon starts for each transcript
		push(@{$transcript_exon_end{$temp[0]}{$attributes[1]}{$temp[2]}}, $temp[4]);	#exon ends for each transcript
		push(@{$transcript_exon_data{$temp[0]}{$attributes[1]}{$temp[2]}}, $attributes[1]); #transcript IDs for each exon
		
	}
	
}

close(GTF);		#close GTF file

print "GTF parsed:\n";
system("date");


#--------------------------------------------------
#
#	Open User Input File and pass into memory
#
#--------------------------------------------------


#==================================================
#
#	The user input should be a tab-delimited 3 column file
#	The first column should be the chromosome - make sure this matches the first field in the annotation file and the identifier in the sequence file
#	The second column is the location in bp on the chromosome (positions begin at 1)
#	The third column is the actual mutation at that position - one of a,c,t,g
#
#==================================================

#open user input file
open(INPUT, "$snpdat_files{-i}") or die "Error! No user input file specified\n";

my @uchrnum;
my @usnppos;
my @usnpchange;
my @umirror;
my @uindel;
my %chrlist;

#keep track of positions for the dbSNP ref file
my %db_snp_pos;

my @nchr;					#new list of CHR that do not have sequence data
my @npos;					#new list of SNp positions that do not have sequence data
my @nmut;					#new list of user mutations that do not have sequence data
my @nmirror;				#array to store a value of 1 or 0 to tell me if the current line is within a feature
my @nindel;

my $vcf = 0;
my $tab = 0;

#read into memory the user input file
while(<INPUT>){
	$summary_user_snps++;
	
	$_ =~ s/\r|\n//g;							#the same as chomp but also dos carriage returns
	
	if($. == 1){
		if($_ =~ m/^##fileformat\=vcf/i){
			$vcf++;
		}else{
			$tab++;
		}
	}
	
	#skip lines with a leading '#'
	next if($_ =~ /^\#/);
	
	######################################
	#
	#	Will need to change the below code to check for vcf format
	#	might want to place the below (until sequence file sending to arrays) in a subroutine and just call the sub
	#
	#
	#	or otherwise (a quick and dirty solution that should work fine) is to nest all of the below in a for loop
	#	place the current line input into three array variable and loop through the array
	#
	######################################
	
	#split the current line of the input file into iarray (input array - this is only a temp array)
	my @iarray = split("\t", $_);
	my $chr;
	my $pos;
	my $mut;
	
	my $isize = @iarray;
	
	if($tab !=0){
		#then we are using the simple tab-delimited file structure
		# if in here then implement the code already in SNPdat
		
		$chr = $iarray[0];
		$pos = $iarray[1];
		$mut = $iarray[2];
		
	}elsif($vcf != 0){
		#then we are using the VCF file format for SNPs
		# if here check the ref and alt variant size difference. 
		#  if the alt has length > 1 OR ref has length > 1 then they can be pushed onto the no seq info array
		#   if the alternative has a comma, split on the comma and create a new set of variables, check their length
		#    if length of new variables >1 then pass to no seq info array with chr and pos (all to be treated as seperate queries)
		#     
		
		$chr = $iarray[0];
		$pos = $iarray[1];
		$mut = $iarray[4];
		
	}
	
	$summary_user_chrs{$chr} = 0;
	
	$chr = "\U$chr";							#convert user chromosome input to upper case to match sequence and annotation file
	$mut = "\U$mut";							#convert user mutation input to upper case to match sequence and codon data
	
	$db_snp_pos{$chr}{$pos} = 'NA';
	
	my $chr_test = 0;
	my $pos_test = 0;
	my $mut_test = 0;

	#split mut based on comma - this is to handle VCF input for multiple bases
	my @varray = split(",", $mut);

	for(my $i = 0; $i < @varray; $i++){
		
		#check input data for errors
		
		#check chromosome
		if(exists $chr_list{$chr}){
			
			$chr_test = 1;
			
		}else{
			print OUTPUT "$chr\t";
			print OUTPUT "$pos\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			print OUTPUT "NA\t";
			
			if($db_snp_file > 0){
				print OUTPUT "NA\t";
			}
			
			print OUTPUT "User chromosome ($chr) not found in GTF";
			$summary_snp_skipped++;
			$summary_noGTF_snp++;
			
		}
		
		#check position
		if($chr_test == 1){
			if($pos =~ m/[a-z]/i){
				print OUTPUT "$chr\t";
				print OUTPUT "$pos\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
				print OUTPUT "NA\t";
			
			if($db_snp_file > 0){
				print OUTPUT "2NA\t";
			}
			
				print OUTPUT "User mutation position ($pos) not valid - Not a number";
				
				$summary_snp_skipped++;
				
			}else{
				if($pos <= 0){
					print OUTPUT "$chr\t";
					print OUTPUT "$pos\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
					print OUTPUT "NA\t";
				
				if($db_snp_file > 0){
					print OUTPUT "NA\t";
				}
				
					print OUTPUT "User mutation position ($pos) not valid - negative or zero";
					
					$summary_snp_skipped++;
					
				}else{
					$pos_test = 1;
					
				}
			}
			
		}elsif($chr_test == 0){
			if($pos =~ m/[a-z]/i){
				
				print OUTPUT ";User mutation position ($pos) not valid - Not a number";
				
			}else{
				if($pos <= 0){
					
					print OUTPUT ";User mutation position ($pos) not valid - negative or zero";
					
				}else{
					$pos_test = 1;
				
				}
			}
			
		}
		
		######################################
		#
		#	Will need a change on the below code
		#	
		#	add a elsif statement to check the mutation input for the user
		#
		######################################

#		my @varray = split(",", $mut);
#		for(my $i = 0; $i < @varray; $i++){
			
		#This will handle base deletions given in the VCF that are one base long. They are not SNPs and thus sequence annotation will not be done
		if($vcf != 0){
			if(length($varray[$i]) == length($iarray[3])){
				if($pos_test == 0 or $chr_test == 0){
					if($varray[$i] =~ m/a|t|c|g/i and length($varray[$i]) == 1){
						
						$mut_test = 1;
						print OUTPUT "\n";
						
					}else{
						
						print OUTPUT ";1User mutation ($varray[$i]) does not represent a valid base\n";
						
					}
				}else{	# - need to put an elsif length($mut) != 1 or eq '-' {push onto the no sequence information arrays at line 792}
					#push onto respective arrays if the chromosome and position are ok - does not matter about the mutation at this point
					push(@uchrnum, $chr);
					push(@usnppos, $pos);
					push(@usnpchange, $varray[$i]);
					push(@umirror, 0);
					push(@uindel, 0);
					
					#define existence of Hash key
					$chrlist{$chr} = 0;
				}
			}else{
				if($pos_test == 0 or $chr_test == 0){
					if($varray[$i] =~ m/a|t|c|g/i and length($varray[$i]) == 1){
						
						$mut_test = 1;
						print OUTPUT ";User mutation ($varray[$i]) is not a SNP (May be an InDel)\n";
						
					}else{
						
						print OUTPUT ";1User mutation ($varray[$i]) does not represent a valid base\n";
						
					}
				}else{
					
					push(@nchr, $chr);
					push(@nmut, $pos);
					push(@npos, $varray[$i]);
					push(@nmirror, 0);
					push(@nindel, 1);
					
					#define existence of Hash key
					$chrlist{$chr} = 0;
				}
			}
		}else{
			########
			
			
			if($pos_test == 0 or $chr_test == 0){
				if($varray[$i] =~ m/a|t|c|g/i and length($varray[$i]) == 1){
					
					$mut_test = 1;
					print OUTPUT "\n";
					
				}else{
					
					print OUTPUT ";1User mutation ($varray[$i]) does not represent a valid base\n";
					
				}
			}else{	
			
				push(@uchrnum, $chr);
				push(@usnppos, $pos);
				push(@usnpchange, $varray[$i]);
				push(@umirror, 0);
				push(@uindel, 0);
				
				#define existence of Hash key
				$chrlist{$chr} = 0;
			}
			
			
			
			
			########
		}
	}
}

my $usize = @uchrnum;
$summary_annotated = @uchrnum;

if($db_snp_file > 0){
	open(DBSNP, "$snpdat_files{-d}") or die "Error opening the processed dbSNP file\n";
	
	while(<DBSNP>){
		$_ =~ s/\r|\n//g;
		$_ =~ s/\s*$//;
		my ($dbchr, $dbpos, $dbrsid) = split("\t", $_);
		
		$dbchr = "\U$dbchr";
		
		if(exists $db_snp_pos{$dbchr}{$dbpos}){
			$db_snp_pos{$dbchr}{$dbpos} = $dbrsid;
		}
	}
	
	close(DBSNP);
}

#close user input file
close(INPUT);

print "User input parsed:\n";
system("date");

#-------------------------------------------
#
#	Slurp into memory the sequence (FASTA) file by chromosome
#
#-------------------------------------------

#===========================================
#
#	Fasta format begins with a single-line description followed by lines of sequence data
#	The description line is distinguished from the sequence data by a ">" at the start of the line
#	The word following the ">" symbol is the sequence identifier and the rest is the description (Both optional)
#	There should be no space between the ">" and first letter of the identifier
#
#===========================================


#store the current line seperator
my $linesep = $/;

#change the default line seperator
$/ = ">";

open(SEQUENCE, "$snpdat_files{-f}") or die "Error opening file\n";

my $chr_count = 0;
my $chr_numb = scalar keys %chr_list;	#this is also the number of chromosomes that will be processed


while(<SEQUENCE>){
	
	next if ($. == 1);
	
	#skip the first line as this is the information before the first > (This value should be null)
	if($. != 1){
		my ($temp1, $temp2) = split("\n", $_, 2);
		
		#will neeed to change this to split on whitespace and take only the first value
		#this is because standard FASTA files have the header and information for the feature type
		#this information is seperated typically by a space. This will ensure tthat the correct
		#part of the header is utilised
		
		#	*Chris recommended no need to do this
		
		$temp1 = "\U$temp1";
		
#		my @temp_val = split()
		
		
		
		#
		#	Include code to remove any leading chr 
		#	ie $temp1 =~ s/^chr//i'
		#
		#	Also consider removing everthing after the whitespace
		#	ie $temp1 =~ s/\s+.*$//g;
		
		#	*Also Chris recommended not to do this
		
		
		#only pass a sequence to memory if the sequence will be used by user input data
		
		#	code in an escape here if the number of processed chromosomes is equal to the number of chromosomes in the input file
		#	Hash %chrlist contains all the user input chromosomes
		
		if(exists $chrlist{$temp1}){
			chomp $temp2;
			chomp $temp1;
			$temp2 =~ s/>[.]*$//g;
			$temp2 =~ s/\s*//g;
			
			$temp1 = "\U$temp1";
			$temp2 = "\U$temp2";
			
			my $ref = \$temp2;
			
			&main($temp1, $ref);
			$chr_count++;
		}
	}
	
	last if($chr_count == $chr_numb);
}

#reset the default line seperator
$/ = $linesep;

close(SEQUENCE);	#close FASTA file

print "Finished analysing all SNPs with sequence information\n:";
system("date");

#------------------------------
#
#Section 3 - Now look at SNPs that have no sequence data available
#
#------------------------------

#arrays to store the user snps (that have no sequence info) that dont occure within a feature
my @nnonfeatchr;
my @nnonfeatpos;
my @nnonfeatmut;
my @nnonfeatindel;

#arrays to store the user SNPs (that have no sequence info) that occure within a feature
my @nfeatchr;
my @nfeatpos;
my @nfeatmut;
my @nfeatindel;

#pass information that doesnt have sequence data to new arrays
#for definitions of the @n[arrays] jump to line before reading user input 
for(my $i = 0; $i < $usize; $i++){

	if($umirror[$i] == 0){
		push(@nchr, $uchrnum[$i]);
		push(@nmut, $usnpchange[$i]);
		push(@npos, $usnppos[$i]);
		push(@nmirror, 0);
		push(@nindel, $uindel[$i]);	#zero here is for SNPs, 1 is for indels (including indels of size 1)
		
	}
}

$summary_nofasta_snps = @nchr;

#find out which the SNPs without sequence information occur within a feature
my $nsize = @nchr;
for(my $i = 0; $i < $nsize; $i++){

	my $gtfsize = @{$seqname{$nchr[$i]}};

	for(my $j = 0; $j < $gtfsize; $j++){
		
		if($npos[$i] > $start{$nchr[$i]}[$j] and $npos[$i] < $end{$nchr[$i]}[$j]){
			$nmirror[$i] = 1;			#1 - indicates that the SNP occurs in a feature
			$summary_in_feature_snps++;
		}
		
	}
	
}

#push SNPs with no seqeunce data onto relevant arrays
for(my $i = 0; $i < $nsize; $i++){
	
	if($nmirror[$i] == 0){
		push(@nnonfeatchr, $nchr[$i]);
		push(@nnonfeatpos, $npos[$i]);
		push(@nnonfeatmut, $nmut[$i]);
		push(@nnonfeatindel, $nindel[$i]);
		
	}elsif($nmirror[$i] == 1){
		push(@nfeatchr, $nchr[$i]);
		push(@nfeatpos, $npos[$i]);
		push(@nfeatmut, $nmut[$i]);
		push(@nfeatindel, $nindel[$i]);
		
	}
}

#now traverse each of the above arrays pulling out the relevant information from the gtf

#traversing the feature array
my $nfeatsize = @nfeatchr;

for(my $i = 0; $i < $nfeatsize; $i++){

	#count the number of annotations at this snp position
	#use sub routine feature check
	
	my @featurecheck = &featurecheck($nfeatchr[$i], $nfeatpos[$i]);
	
	my %features;							#hash to keep a list of the features that it is contained within - this will always be unique
	my $fsize = @featurecheck;
	
	#create the unique list of features
	for(my $j = 1; $j < $fsize; $j++){
		$features{$feature{$nfeatchr[$i]}[$featurecheck[$j]]} = 0;
		
	}
	
	#count the number of times each feature appears for the current SNP - will provide a total for the feature
	for(my $j = 1; $j < $fsize; $j++){
		
		if(exists $features{$feature{$nfeatchr[$i]}[$featurecheck[$j]]}){
			$features{$feature{$nfeatchr[$i]}[$featurecheck[$j]]}++;
			
		}
	}
	
	#total number of different features that SNP is in
	my @numfeature = keys %features;
	my $numsize = @numfeature;
	
	#traverse the user input to find and print out the current number of the annotation then increment
	foreach my $key (keys %features){
		
		my $counter = 0;
		my $fstrand;								#this is the strand sense for the current feature
		
		for(my $j = 1; $j < $fsize; $j++){
			
			if($feature{$nfeatchr[$i]}[$featurecheck[$j]] eq $key){
				
				$counter++;
				
				#print everything exon from the annotations
				print OUTPUT "$nfeatchr[$i]\t";																		#user chr - note probably should be the listchr array
				print OUTPUT "$nfeatpos[$i]\t";																		#user SNP position
				print OUTPUT "Y\t";																					#whether or not in a feature
				print OUTPUT "Exonic\t";
				print OUTPUT "NA\t";																				#distance to nearest feature
				print OUTPUT "$feature{$nfeatchr[$i]}[$featurecheck[$j]]\t";														#feature of SNP
				print OUTPUT "$numsize\t";																			#number of different features the SNP is in
				print OUTPUT "[$counter/$features{$key}]\t";														#number of annotations for that feature
				print OUTPUT "$start{$nfeatchr[$i]}[$featurecheck[$j]]\t";															#start of current feature
				print OUTPUT "$end{$nfeatchr[$i]}[$featurecheck[$j]]\t";															#end of current feature
				print OUTPUT "$gene_id{$nfeatchr[$i]}[$featurecheck[$j]]\t";														#gene_id for the current feature
				print OUTPUT "$gene_name{$nfeatchr[$i]}[$featurecheck[$j]]\t";														#gene_name
				print OUTPUT "$transcript_id{$nfeatchr[$i]}[$featurecheck[$j]]\t";													#transcript_id
				print OUTPUT "$transcript_name{$nfeatchr[$i]}[$featurecheck[$j]]\t";												#transcript_name
				print OUTPUT "[$exon_number{$nfeatchr[$i]}[$featurecheck[$j]]"."/$gene_list{$gene_id{$nfeatchr[$i]}[$featurecheck[$j]]}]\t";		#annotated exon_number
				
				if($strand{$nfeatchr[$i]}[$featurecheck[$j]] ne "\."){
					print OUTPUT "$strand{$nfeatchr[$i]}[$featurecheck[$j]]\t";
					
				}else{
					print OUTPUT "$strand{$nfeatchr[$i]}[$featurecheck[$j]]**\t";
					
				}
				
				if($frame{$nfeatchr[$i]}[$featurecheck[$j]] eq "\."){
					print OUTPUT "NA\t";
					
				}elsif($frame{$nfeatchr[$i]}[$featurecheck[$j]] == 0 or $frame{$nfeatchr[$i]}[$featurecheck[$j]] == 1 or $frame{$nfeatchr[$i]}[$featurecheck[$j]] == 2){
					
					if($strand{$nfeatchr[$i]}[$featurecheck[$j]] eq "+"){						#positive strand
						my $rf_change = ($frame{$nfeatchr[$i]}[$featurecheck[$j]]+1);			#increment the GTF reading frame value so we never have a -0 reading frame
						
						print OUTPUT "$rf_change\t";				#annotated reading frame
						
					}else{														#negative strand
						my $rf_change = ($frame{$nfeatchr[$i]}[$featurecheck[$j]]+1);
						
						print OUTPUT "-"."$rf_change\t";	
						
					}
				}
				
				print OUTPUT "NA\t";								#cant estimate reading frame
				print OUTPUT "NA\t";								#cant estimate number of stop codons
				print OUTPUT "NA\t";								#cant estimate codon
				print OUTPUT "NA\t";								#cant estimate amino acid
				print OUTPUT "NA\t";								#cant estimate synonymous
				print OUTPUT "$protein_id{$nfeatchr[$i]}[$featurecheck[$j]]\t";
		
		if($db_snp_file > 0){
			print OUTPUT "$db_snp_pos{$nfeatchr[$i]}{$nfeatpos[$i]}\t";
		}
		
		if($nfeatindel[$i] == 0){
				print OUTPUT "There was no sequence information available for SNPs on this chromosome";
		}#else{
		#	print OUTPUT "Warning the user is an indel - no sequence annotation performed";
		#}
				if($nfeatmut[$i] =~ m/a|t|c|g/i and length($nfeatmut[$i]) == 1 and $nfeatindel[$i] == 0){
					print OUTPUT "\n";
				}elsif($nfeatmut[$i] =~ m/a|t|c|g/i and length($nfeatmut[$i]) == 1 and $nfeatindel[$i] == 1){
					print OUTPUT "Warning user mutation ($nfeatmut[$i]) is not a valid base (Its an InDel)\n";
				}else{
					print OUTPUT "; Warning user mutation ($nfeatmut[$i]) is not a valid base\n";
				}
				
			}
		}
	}
}

#------------------------------
#
#	traversing the non feature array - not within a feature
#	Need to find distances to nearest features first before I can call subroutine 'gtfnf'
#
#-------------------------------

my $nnonfeatchrsize = @nnonfeatchr;
for(my $i = 0; $i < $nnonfeatchrsize; $i++){
	my $distst = 999999999999;						#set a max dist (always larger than a possible dist for genome)
	my $distend = 999999999999;
	my $asize = @{$seqname{$nnonfeatchr[$i]}};							#number of lines in GTF array
	
	
	my $trans_count = 0;		##found out if the current SNP is intronic or intergenic (>0 = intronic)
	
	#
	#	get the transcripts that the current SNP is a member of
	#
	
	my @trans_member;	#keep track of the transcripts that the SNP is a member of
	
	##put counter here to keep track of whether the SNP is within a transcript or not
	
	##this is a check to find all transcripts that each SNP is found within the min and max boundary of
	foreach my $key (keys %{$transcript_max{$nnonfeatchr[$i]}}){
		
		if($nnonfeatpos[$i] > $transcript_min{$nnonfeatchr[$i]}{$key} and $nnonfeatpos[$i] < $transcript_max{$nnonfeatchr[$i]}{$key}){
			push(@trans_member, $key);
			$trans_count++;
		}
		
	}
	
	#
	#
	#
	#
	#
	
	##trans_count > 0 - intronic else intergenic
	
	
	##Now get the distance to the feature within the transcripts for @trans_member
	##the transcript_data, transcript_start, transcript_end hashes will be used
	#
	#
	#If statement here to check whether or not the counter has been incremented
	#	if it has then use the transcript search otherwise use the already existing code
	#
	##it must be intronic to make it in here
	if($trans_count > 0){
		##need size of array for features within transcript
		for(my $a = 0; $a <@trans_member; $a++){
			my $trans_size = @{$transcript_data{$nnonfeatchr[$i]}{$trans_member[$a]}};
			
			#dist to nearest start point of a feature within the transcript
			for(my $b = 0; $b < $trans_size; $b++){
				if(($transcript_start{$nnonfeatchr[$i]}{$trans_member[$a]}[$b] - $nnonfeatpos[$i]) < $distst and ($transcript_start{$nnonfeatchr[$i]}{$trans_member[$a]}[$b] - $nnonfeatpos[$i]) > 0){
					$distst = ($transcript_start{$nnonfeatchr[$i]}{$trans_member[$a]}[$b] - $nnonfeatpos[$i]);
				}
			}
			
			#dist to nearest end point of feature within transcript
			for(my $b = 0; $b < $trans_size; $b++){
				if(($nnonfeatpos[$i] - $transcript_end{$nnonfeatchr[$i]}{$trans_member[$a]}[$b]) < $distend and ($nnonfeatpos[$i] - $transcript_end{$nnonfeatchr[$i]}{$trans_member[$a]}[$b]) > 0){
					$distend = ($nnonfeatpos[$i] - $transcript_end{$nnonfeatchr[$i]}{$trans_member[$a]}[$b]);
				}
			}
			
			##Need to call gtfnf after finding the smallest distance
			##add a parameter to the subroutine to tell it if the passed value is intonic or intergenic
			
			my $chr = $nnonfeatchr[$i];
			my $pos = $nnonfeatpos[$i];
			my $mut = $nnonfeatmut[$i];
			
			my $indel = $nnonfeatindel[$i];	#need to track whether or not the mutation is an indel
			
			if($distend == $distst){
				&gtfnf(0, 1, $distend, $distst, $chr, $pos, $mut, $indel, 0);
				
			}elsif($distend > $distst){
				&gtfnf(1, 1, $distst, $chr, $pos, $mut, $indel, 0);
				
			}elsif($distst > $distend){
				&gtfnf(2, 1, $distend, $chr, $pos, $mut, $indel, 0);
				
			}
			
		}
		
		
	}else{
	
	
		#get dist to nearest start point
		for(my $a = 0; $a < $asize; $a++){
			
			if(($start{$nnonfeatchr[$i]}[$a] - $nnonfeatpos[$i]) < $distst and ($start{$nnonfeatchr[$i]}[$a] - $nnonfeatpos[$i]) > 0){
				$distst = ($start{$nnonfeatchr[$i]}[$a] - $nnonfeatpos[$i]);
				
			}
		}
		
		#get dist to nearest end point
		for(my $a = 0; $a < $asize; $a++){
			
			if(($nnonfeatpos[$i] - $end{$nnonfeatchr[$i]}[$a]) < $distend and ($nnonfeatpos[$i] - $end{$nnonfeatchr[$i]}[$a]) > 0){
				$distend = ($nnonfeatpos[$i] - $end{$nnonfeatchr[$i]}[$a]);
				
			}
		}
		
		my $chr = $nnonfeatchr[$i];
		my $pos = $nnonfeatpos[$i];
		my $mut = $nnonfeatmut[$i];
		my $indel = $nnonfeatindel[$i];
		
		#in gtf sub need to print out cols from dist to 
		if($distend == $distst){
			&gtfnf(0, 0, $distend, $distst, $chr, $pos, $mut, $indel, 1);
			
		}elsif($distend > $distst){
			&gtfnf(1, 0, $distst, $chr, $pos, $mut, $indel, 1);
			
		}elsif($distst > $distend){
			&gtfnf(2, 0, $distend, $chr, $pos, $mut, $indel, 1);
			
		}
	}
}

print "Finished analysing all SNPs with no sequence information\n:";
system("date");

#-----------------------------
#
#	print explanation of symbols here
#
#-----------------------------

#	*** - The position at this location is four fold degenerate
#	**** - Codon could not be confirmed as synonymous or nonsynonymous due to missing sequence information
#	*^ - the SNP was eqidistance from 2 features
#	** - no annotated strand sense so strand sense was assumed to be '+'


#only data with a sequence is passed to here
sub main{

	#pass the current chromosome and sequence from the sequence file to the analysis code
	my ($val, $ref2) = @_;

	#these will be list of u input SNPs that match the current sequence
	my @ulistchr;
	my @ulistpos;
	my @ulistmut;
	my @ulistindel;

	my $size = @uchrnum;
	for(my $i = 0; $i < $size; $i++){
		
		if($uchrnum[$i] eq $val){
			$umirror[$i] = 1;
			#change the mirror to 1 to show that sequence information was available for this SNP
			
			push(@ulistchr, $uchrnum[$i]);
			push(@ulistpos, $usnppos[$i]);
			push(@ulistmut, $usnpchange[$i]);
			push(@ulistindel, $uindel[$i]);
			
		}
	}

	my $ulist = @ulistchr;

	for(my $i = 0; $i < $ulist; $i++){
		#check if the user input current line is within a feature or not
		#if $featurecheck[0] == 0 than the SNP is not within a feature
		
		#could always create a mirror array to do the search just once!
		my @featurecheck = &featurecheck($ulistchr[$i], $ulistpos[$i]);	
		#featurecheck contains the index positions for each relevant feature
		
		if($featurecheck[0] == 0){
			#SNP position is not within a feature
			
			my $distst = 999999999999;		#set a max dist (always larger than a possible dist for genome)
			my $distend = 999999999999;
			my $asize = @{$seqname{$ulistchr[$i]}};		#number of lines in GTF array for the current chromosome
			my $trans_count = 0;		##found out if the current SNP is intronic or intergenic (>0 = intronic)
			
			#
			#	get the transcripts that the current SNP is a member of
			#
			
			my @trans_member;	#keep track of the transcripts that the SNP is a member of
			
			##put counter here to keep track of whether the SNP is within a transcript or not
			
			##this is a check to find all transcripts that each SNP is found within the min and max boundary of
			foreach my $key (keys %{$transcript_max{$ulistchr[$i]}}){
				
				if($ulistpos[$i] > $transcript_min{$ulistchr[$i]}{$key} and $ulistpos[$i] < $transcript_max{$ulistchr[$i]}{$key}){
					push(@trans_member, $key);
					$trans_count++;
				}
				
			}
			
			##trans_count > 0 - intronic else intergenic
			
			
			##Now get the distance to the feature within the transcripts for @trans_member
			##the transcript_data, transcript_start, transcript_end hashes will be used
			#
			#
			#If statement here to check whether or not the counter has been incremented
			#	if it has then use the transcript search otherwise use the already existing code
			#
			##it must be intronic to make it in here
			if($trans_count > 0){
				##need size of array for features within transcript
				for(my $a = 0; $a <@trans_member; $a++){
					my $trans_size = @{$transcript_data{$ulistchr[$i]}{$trans_member[$a]}};
					
					#dist to nearest start point of a feature within the transcript
					for(my $b = 0; $b < $trans_size; $b++){
						if(($transcript_start{$ulistchr[$i]}{$trans_member[$a]}[$b] - $ulistpos[$i]) < $distst and ($transcript_start{$ulistchr[$i]}{$trans_member[$a]}[$b] - $ulistpos[$i]) > 0){
							$distst = ($transcript_start{$ulistchr[$i]}{$trans_member[$a]}[$b] - $ulistpos[$i]);
						}
					}
					
					#dist to nearest end point of feature within transcript
					for(my $b = 0; $b < $trans_size; $b++){
						if(($ulistpos[$i] - $transcript_end{$ulistchr[$i]}{$trans_member[$a]}[$b]) < $distend and ($ulistpos[$i] - $transcript_end{$ulistchr[$i]}{$trans_member[$a]}[$b]) > 0){
							$distend = ($ulistpos[$i] - $transcript_end{$ulistchr[$i]}{$trans_member[$a]}[$b]);
						}
					}
					
					##Need to call gtfnf after finding the smallest distance
					##add a parameter to the subroutine to tell it if the passed value is intonic or intergenic
					
					my $chr = $ulistchr[$i];
					my $pos = $ulistpos[$i];
					my $mut = $ulistmut[$i];
					
					my $indel = $ulistindel[$i];	#need to track whether or not the mutation is an indel
					
					if($distend == $distst){
						&gtfnf(0, 1, $distend, $distst, $chr, $pos, $mut, $indel, 0);
						
					}elsif($distend > $distst){
						&gtfnf(1, 1, $distst, $chr, $pos, $mut, $indel, 0);
						
					}elsif($distst > $distend){
						&gtfnf(2, 1, $distend, $chr, $pos, $mut, $indel, 0);
						
					}
					
				}
				
				
			}else{
				
				##It must be intergenic to make it in here
				
				#get dist to nearest start point
				for(my $a = 0; $a < $asize; $a++){
					
					if(($start{$ulistchr[$i]}[$a] - $ulistpos[$i]) < $distst and ($start{$ulistchr[$i]}[$a] - $ulistpos[$i]) > 0){
						$distst = ($start{$ulistchr[$i]}[$a] - $ulistpos[$i]);
					}
				}
				
				#get dist to nearest end point
				for(my $a = 0; $a < $asize; $a++){
					
					if(($ulistpos[$i] - $end{$ulistchr[$i]}[$a]) < $distend and ($ulistpos[$i] - $end{$ulistchr[$i]}[$a]) > 0){
						$distend = ($ulistpos[$i] - $end{$ulistchr[$i]}[$a]);
						
					}
				}
				
				my $chr = $ulistchr[$i];
				my $pos = $ulistpos[$i];
				my $mut = $ulistmut[$i];
				
				
				my $indel = $ulistindel[$i];	#need to track whether or not the mutation is an indel
				
				
				#in gtf sub need to print out cols from dist to 
				if($distend == $distst){
					&gtfnf(0, 1, $distend, $distst, $chr, $pos, $mut, $indel, 1);
					
				}elsif($distend > $distst){
					&gtfnf(1, 1, $distst, $chr, $pos, $mut, $indel, 1);
					
				}elsif($distst > $distend){
					&gtfnf(2, 1, $distend, $chr, $pos, $mut, $indel, 1);
					
				}
			}
			
		}else{
			#SNP position is within a feature
			my %features;	#hash to keep a list of the features that it is contained within - this will always be unique
			my $fsize = @featurecheck;
			$summary_in_feature_snps++;
			
			#create the unique list of features
			for(my $j = 1; $j < $fsize; $j++){
				$features{$feature{$ulistchr[$i]}[$featurecheck[$j]]} = 0;
				
			}
			
			#count the number of times each feature appears for the current SNP - will provide a total for the feature
			for(my $j = 1; $j < $fsize; $j++){
				
				if(exists $features{$feature{$ulistchr[$i]}[$featurecheck[$j]]}){
					$features{$feature{$ulistchr[$i]}[$featurecheck[$j]]}++;
					
				}
			}
			
			#total number of different features that SNP is in
			my @numfeature = keys %features;
			my $numsize = @numfeature;
			
			#traverse the user input to find and print out the current number of the annotation then increment
			foreach my $key (keys %features){
				
				my $counter = 0;
				my $fstrand;	#this is the strand sense for the current feature
				
				for(my $j = 1; $j < $fsize; $j++){
					
					if($feature{$ulistchr[$i]}[$featurecheck[$j]] eq $key){	#I think this should be $feature[$featurecheck[$j]]
						$counter++;
						
						#print everything exon from the annotations
						print OUTPUT "$ulistchr[$i]\t";																		#user chr - note probably should be the listchr array
						print OUTPUT "$ulistpos[$i]\t";																		#user SNP position
						print OUTPUT "Y\t";																					#Whether or not in a feature
						print OUTPUT "Exonic\t";
						print OUTPUT "NA\t";																				#distance to nearest feature
						print OUTPUT "$feature{$ulistchr[$i]}[$featurecheck[$j]]\t";														#feature of SNP
						print OUTPUT "$numsize\t";																			#number of different features the SNP is in
						print OUTPUT "[$counter/$features{$key}]\t";														#number of annotations for that feature
						print OUTPUT "$start{$ulistchr[$i]}[$featurecheck[$j]]\t";															#start of current feature
						print OUTPUT "$end{$ulistchr[$i]}[$featurecheck[$j]]\t";															#end of current feature
						print OUTPUT "$gene_id{$ulistchr[$i]}[$featurecheck[$j]]\t";														#gene_id for the current feature
						print OUTPUT "$gene_name{$ulistchr[$i]}[$featurecheck[$j]]\t";														#gene_name
						print OUTPUT "$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]\t";													#transcript_id
						print OUTPUT "$transcript_name{$ulistchr[$i]}[$featurecheck[$j]]\t";												#transcript_name
						print OUTPUT "[$exon_number{$ulistchr[$i]}[$featurecheck[$j]]"."/$gene_list{$gene_id{$ulistchr[$i]}[$featurecheck[$j]]}]\t";		#annotated exon_number
						
						#If the annotated strand sense is not given the script will assume it is '+'
						if($strand{$ulistchr[$i]}[$featurecheck[$j]] eq "+"){
							print OUTPUT "$strand{$ulistchr[$i]}[$featurecheck[$j]]\t";			#annotated strand sense
							$fstrand = 1;				#1 means the feature is on the positive strand
							
						}elsif($strand{$ulistchr[$i]}[$featurecheck[$j]] eq "-" ){
							print OUTPUT "$strand{$ulistchr[$i]}[$featurecheck[$j]]\t";
							$fstrand = 0;				#0 means the feature is on the negative strand
							
						}else{
							print OUTPUT "+**\t";
							$fstrand = 1;				#1 means the feature is on the positive strand
							
						}
						
						#check to see if there is an annotated reading frame
						my $rfcheck = 0;
						
						if($frame{$ulistchr[$i]}[$featurecheck[$j]] eq "\."){
							print OUTPUT "NA\t";
							$rfcheck = 1;				#$rfcheck == 1 then there was no annotated reading frame
							
						}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 0 or $frame{$ulistchr[$i]}[$featurecheck[$j]] == 1 or $frame{$ulistchr[$i]}[$featurecheck[$j]] == 2){
							
							if($fstrand == 1){				#positive strand
								my $rf_change = ($frame{$ulistchr[$i]}[$featurecheck[$j]]+1);
								print OUTPUT "$rf_change\t";				#annotated reading frame
								
							}else{							#negative strand
								my $rf_change = ($frame{$ulistchr[$i]}[$featurecheck[$j]]+1);
								print OUTPUT "-"."$rf_change\t";	
								
							}
							
							print OUTPUT "NA\t";			#print NA for estimated reading frame - no need to estimate as there is an annotated reading frame
							
						}
						
						if($ulistmut[$i] =~ m/a|t|c|g/i and length($ulistmut[$i]) == 1){
							
							#
							#	Extract the sequence for the current feature
							#
							
							#
							#	Need to get the mod of the start position of the feature modulo 3
							#	
							
							#	The feature start - start(mod3) = substring offset
							
							
							#annotated start
							my $start = $start{$ulistchr[$i]}[$featurecheck[$j]];
							
							
							my $rel_start_snp = $ulistpos[$i] - $start;
							$rel_start_snp++;
							
							my $rel_end_snp = $end{$ulistchr[$i]}[$featurecheck[$j]] - $ulistpos[$i];
							$rel_end_snp++;
							
							#annotated end
							my $end = $end{$ulistchr[$i]}[$featurecheck[$j]];
							
							my $modulo = $start%3;
							my $new_start;
							my $start_grab = 0;
							#estimate start o sequence to extract
							#I think that this information should be taken from the previous exon it will now be done with the -x option
							if($modulo == 2){
								$new_start = ($start-1);	#should be last base of previous exon otherwise ?
								$start_grab = 1;
								$rel_start_snp += 1;
							}elsif($modulo == 0){
								$new_start = ($start-2);	#should be last 2 bases of previous exon otherwise ??						
								$start_grab = 2;
								$rel_start_snp += 2;
							}else{
								$new_start = $start;		#no change
							}
							
							
							
							my $temp_dist = ($end{$ulistchr[$i]}[$featurecheck[$j]] - $new_start)+1;
							my $val =  $temp_dist%3;
							
							my $new_end;
							my $end_grab = 0;
							#I think that this information should be taken from the next exon
							if($val == 1){
								$new_end = $end{$ulistchr[$i]}[$featurecheck[$j]]+2;	#should be first two bases of next exon otherwise ?
								$end_grab = 2;
								$rel_end_snp += 2;
							}elsif($val == 2){
								$new_end = $end{$ulistchr[$i]}[$featurecheck[$j]]+1;	#should be first base of next exon otherwise ?						
								$end_grab = 1;
								$rel_end_snp += 1;
							}else{
								$new_end = $end{$ulistchr[$i]}[$featurecheck[$j]]	#no change
							}
							
							my $new_dist_final = ($new_end - $new_start)+1;
							my $sub_start = $new_start - 1;
							my $sub_end = $new_end - 1;
							
							#print "Start grab is $start_grab\tEnd grab is $end_grab\tlength mod3 is $val\tsubtraction coordinates: $sub_start\t$sub_end\n";
							
							my $new_snp_dist = ($ulistpos[$i] - $new_start)+1;
							
							my $val2 = $$ref2;
							
							#perl substring begins at zero
							my $sequence2 = substr($val2, $sub_start, $new_dist_final);		#extract the relevant piece of the sequence
							#print "$sequence2\n";
							
							#now depending on the value of start_grab and end_grab
							#substitute start and end pieces with ? or ??
							#substitute start and end pieces from the nex/previous exon nucleotides
							#	only with the -x option in force
							
							my $start_e_seq;
							my $pre_pos;
							my $end_e_seq;
							my $pre_pos_end;
							my $stemp1 = "??";
							my $stemp2 = "??";
							
							#maybe an if statement for $val and/or $modulo
							
							#below: The program will see if the -x option is in use
							#if it is, it will grab the nucleotides either side of the exon
							#that is requires to ensure that the sequence is in the rf1
							#than readjust the SNP coordinate relative to the start of the sequence
							
							
							#need a counter for special cases where the amino acid crosses the exon boundary
							if($exon_only == 0){
								
								#not correct = need
								
								if($start_grab == 1){
									
									$sequence2 =~ s/^./?/;
									
								}elsif($start_grab == 2){
									
									$sequence2 =~ s/^../??/;
								}
								
								
								if($end_grab == 1){
									
									$sequence2 =~ s/.$/?/;
									
								}elsif($end_grab == 2){
									$sequence2 =~ s/..$/??/;
								}
								
								$stemp1 = "??";
								$stemp2 = "??";
								
							}elsif(exists $x_features{$feature{$ulistchr[$i]}[$featurecheck[$j]]}){
								
								if($start_grab == 1){
									
									#bases before the first exon will be ?
									if($start - $transcript_exon_min{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} < 2){
										$sequence2 =~ s/^./?/;
										$stemp1 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										
										for(my $x = 0; $x < @{$transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											if($nonzero == 0){
												
												$diff = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												if($diff_t > 0){
													$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
												
											}
										}
										
										$start_e_seq = substr($val2, $pre_pos-1, 1);
										$stemp1 = substr($val2, $pre_pos-1-1,1);
										
										$sequence2 =~ s/^./$start_e_seq/;
									}
									
								}elsif($start_grab == 2){
									#bases before the first exon will be ?
									if($start - $transcript_exon_min{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} < 2){
										$sequence2 =~ s/^../??/;
										$stemp1 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										my $pre_pos;
										
										for(my $x = 0; $x < @{$transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											
											
											if($nonzero == 0){
												
												$diff = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff_t > 0){
													$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
											}
										}
										
										
										$start_e_seq = substr($val2, $pre_pos-1-1, 2);
										$sequence2 =~ s/^../$start_e_seq/;
									}
								}elsif($start_grab == 0){
									if($start - $transcript_exon_min{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} < 2){
										
										$stemp1 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										my $pre_pos;
										
										for(my $x = 0; $x < @{$transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											
											
											if($nonzero == 0){
												
												$diff = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $ulistpos[$i] - $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff_t > 0){
													$pre_pos = $transcript_exon_end{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
											}
										}
										
										
										$start_e_seq = substr($val2, $pre_pos-1-1, 2);
										$stemp1 = $start_e_seq;
									}
								}
								
								
								if($end_grab == 1){
									#bases before the first exon will be ?
									if($transcript_exon_max{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} - $end < 2){
										$sequence2 =~ s/.$/?/;
										$stemp2 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										
										for(my $x = 0; $x < @{$transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											if($nonzero == 0){
												
												$diff = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												if($diff_t > 0){
													$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
												
											}
										}
										
										$start_e_seq = substr($val2, $pre_pos-1, 1);
										
										$sequence2 =~ s/.$/$start_e_seq/;
										
										$stemp2 = substr($val2, $pre_pos, 1);
									}
									
								}elsif($end_grab == 2){
									#bases before the first exon will be ?
									if($transcript_exon_max{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} - $end < 2){
										$sequence2 =~ s/..$/??/;
										$stemp2 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										my $pre_pos;
										
										for(my $x = 0; $x < @{$transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											
											
											if($nonzero == 0){
												
												$diff = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												
												if($diff_t > 0){
													$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
											}
										}
										
										
										$start_e_seq = substr($val2, $pre_pos-1, 2);
										$sequence2 =~ s/..$/$start_e_seq/;
									}
									
								}elsif($end_grab == 0){
									#bases before the first exon will be ?
									if($transcript_exon_max{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]} - $end < 2){
										$stemp2 = "??";
									}else{
										my $diff;
										my $nonzero = 0;
										my $pre_pos;
										
										for(my $x = 0; $x < @{$transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}}; $x++){
											
											
											if($nonzero == 0){
												
												$diff = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x];
												
												if($diff > 0){
													$nonzero = 1;
												}else{
													undef $diff;
													undef $pre_pos;
												}
												
												
											}else{
												my $diff_t = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] - $ulistpos[$i];
												
												if($diff_t > 0){
													$pre_pos = $transcript_exon_start{$ulistchr[$i]}{$transcript_id{$ulistchr[$i]}[$featurecheck[$j]]}{$feature{$ulistchr[$i]}[$featurecheck[$j]]}[$x] if ($diff_t < $diff);
													$diff = $diff_t if ($diff_t < $diff);
													
												}
												
											}
										}
										
										
										$start_e_seq = substr($val2, $pre_pos-1, 2);
										$stemp2 = $start_e_seq;
									}
								}
							}else{
								
								if($start_grab == 1){
									
									$sequence2 =~ s/^./?/;
									
								}elsif($start_grab == 2){
									
									$sequence2 =~ s/^../??/;
								}
								
								
								if($end_grab == 1){
									
									$sequence2 =~ s/.$/?/;
									
								}elsif($end_grab == 2){
									$sequence2 =~ s/..$/??/;
								}
								
								$stemp1 = "??";
								$stemp2 = "??";
							}
							
							
							#change of $sequence2 should occur here based on the previous comments
							#bases from next/previous exons should be used
							
							#cat ?? to the end of the sequence as these are definitely unknown (Off the chromosome)
							my $val3 = "??"."$val2"."??";
							
							#want to extract the 2 nucleotides either side of the feature sequence - might need to cat these to end of sequence
							
							#my $stemp1 = substr($val3, $sub_start, 2);	#this has been deprecated 22-04-13 as it is taking the sequence from the intron/intergenic region
							#my $stemp1 = "??";
							
							
							
							
							
							
							#NEED to set stemp1 and stemp2 to the bases at the next exon if the -x option is in use
							
							
							
							
							#stemp1 and stemp2 should become the bases from the last and next exon if start_grab and end_grab are 0
							
							
							
							
							
							
							
							
							
							
							
							
							
							# no need to add plus 2 to the position as the start of the feature will be plus two from the original start
							# and because we want to cat just two positions if we dont change the start position we should be fine
							# basically the original start position is now the start of the 2 nucleotides before the feature start for 
							# the sequence $val3
							
							#my $stemp2 = substr($val3, ($sub_end+1)+2, 2);#this has been deprecated 22-04-13 as it is taking the sequence from the intron/intergenic region
							#my $stemp2 = "??";
							
							#this is at the end - +1 so we start at the nucleotide beside the last in the sequence
							
							#+2 to both because ?? has been concatonated to the start of the sequence
							
							$val3 = undef;
							
							#$sequence2 is the extracted sequence
							
							my $rsequence = &reversecom($sequence2);
							
							#the nucleotide attachments reverse complimented
							my $rtemp1 = &reversecom($stemp1);	#this will go at the end
							my $rtemp2 = &reversecom($stemp2);	#this will go at the start
							
							#check for annotated reading frame
							my @codoninfo1;
							my @codoninfo2;
							my @codoninfo3;
							my @codoninfo4;
							my @codoninfo5;
							my @codoninfo6;
							my $fframe;
							my $seqcodon;		#the sequence that wil be used to get the codon
							
							if($rfcheck == 1){	#there was an no annotated reading frame
								
								#the below array positions at [0] will have the total number of stop codons in the sequence (excluding the last codon if it was a stop)
								#the below array positions at [1] will have a value to tell you if the last codon was a stop; 1 for yes, 0 for no
								@codoninfo1 = &codoncounter(0, $sequence2, 1);
								@codoninfo2 = &codoncounter(1, $sequence2, 1);
								@codoninfo3 = &codoncounter(2, $sequence2, 1);
								
								#count stop codons on the reverse compliment
								@codoninfo4 = &codoncounter(0, $rsequence, 1);
								@codoninfo5 = &codoncounter(1, $rsequence, 1);
								@codoninfo6 = &codoncounter(2, $rsequence, 1);
								
								#find the reading frame with the least number of stop codons
								#array1 - array of stop codon numbers for each reading frame
								#array2 - keep track of which reading frame goes where
								#array3 - keep track of whether or not this reading frames last codon was a stop
								my @array1 = ($codoninfo1[0], $codoninfo2[0], $codoninfo3[0], $codoninfo4[0], $codoninfo5[0], $codoninfo6[0]);
								my @array2 = (1, 2, 3, 4, 5, 6);
								my @array3 = ($codoninfo1[1], $codoninfo2[1], $codoninfo3[1], $codoninfo4[1], $codoninfo5[1], $codoninfo6[1]);
								
								#simple sort to find what the min number of codons was
								for(my $p = 0; $p < 5; $p++){
									
									for(my $q = 0; $q < 5; $q++){
										
										if($array1[$q] > $array1[$q+1]){
											my $temp = $array1[$q+1];
											my $temp2 = $array2[$q+1];
											my $temp3 = $array3[$q+1];
											$array1[$q+1] = $array1[$q];
											$array2[$q+1] = $array2[$q];
											$array3[$q+1] = $array3[$q];
											$array1[$q] = $temp;
											$array2[$q] = $temp2;
											$array3[$q] = $temp3;
											
										}
									}
								}
								
								my $mincodonnum = $array1[0];
								
								if($array1[0] == $array1[1]){			#in the event of a tie between reading frames
								
								#if there is a tie between reading frames take the reading frame closest to the strand sense
									print OUTPUT "*";
									
									#create a hash with only the reading frames that match the min number
									my %thash;
									
									for(my $r = 0; $r < 6; $r++){
										
										if($array1[$r] == $mincodonnum){
											#create a hash to check for order preference but also set its value to whether or not the last codon was a stop in that reading frame
											$thash{$array2[$r]} = $array3[$r];
										}
									}
									
									my $laststopcheck;
									#if $fstrand == 1 - feature is on the positive strand ; $fstrand == 0 - feature is on the negative strand
									#print in this order of preference
									if($fstrand == 1){
										
										if(exists $thash{1}){
											$laststopcheck = $thash{1};
											$fframe = 1;
											print OUTPUT "1\t";
											$seqcodon = $sequence2;
											
										}elsif(exists $thash{2}){
											$laststopcheck = $thash{2};
											$fframe = 2;
											print OUTPUT "2\t";
											$seqcodon = $sequence2;
											
										}elsif(exists $thash{3}){
											$laststopcheck = $thash{3};
											$fframe = 3;
											print OUTPUT "3\t";
											$seqcodon = $sequence2;
											
										}elsif(exists $thash{4}){
											$laststopcheck = $thash{4};
											$fframe = 4;
											print OUTPUT "-"."1\t";
											$seqcodon = $rsequence;
											
										}elsif(exists $thash{5}){
											$laststopcheck = $thash{5};
											$fframe = 5;
											print OUTPUT "-"."2\t";
											$seqcodon = $rsequence;
											
										}elsif(exists $thash{6}){
											$laststopcheck = $thash{6};
											$fframe = 6;
											print OUTPUT "-"."3\t";
											$seqcodon = $rsequence;
											
										}
										
										
									#feature is on the negative strand so print in this order of preference
									}elsif($fstrand == 0){
										
										if(exists $thash{4}){
											$laststopcheck = $thash{4};
											$fframe = 4;
											print OUTPUT "-"."1\t";
											$seqcodon = $rsequence;
											
										}elsif(exists $thash{5}){
											$laststopcheck = $thash{5};
											$fframe = 5;
											print OUTPUT "-"."2\t";
											$seqcodon = $rsequence;
											
										}elsif(exists $thash{6}){
											$laststopcheck = $thash{6};
											$fframe = 6;
											print OUTPUT "-"."3\t";
											$seqcodon = $rsequence;
											
										}elsif(exists $thash{1}){
											$laststopcheck = $thash{1};
											$fframe = 1;
											print OUTPUT "1\t";
											$seqcodon = $sequence2;
											
										}elsif(exists $thash{2}){
											$laststopcheck = $thash{2};
											$fframe = 2;
											print OUTPUT "2\t";
											$seqcodon = $sequence2;
											
										}elsif(exists $thash{3}){
											$laststopcheck = $thash{3};
											$fframe = 3;
											print OUTPUT "3\t";
											$seqcodon = $sequence2;
											
										}
									}
									
									#print number of stop codons
									if($laststopcheck == 1){
										print OUTPUT "^"."$mincodonnum\t";
										
									}else{
										print OUTPUT "$mincodonnum\t";
										
									}
									
								#no tie between reading frames
								}else{
									
									#print reading frame
									if($array2[0] == 1 or $array2[0] == 2 or $array2[0] == 3){
										print OUTPUT "$array2[0]\t";
										$seqcodon = $sequence2;
										
									}else{
										my $arframe;
										
										if($array2[0] == 4){
											$arframe = 1;
											
										}elsif($array2[0] == 5){
											$arframe = 2;
										
										}elsif($array2[0] == 6){
											$arframe = 3;
										
										}
										
										print OUTPUT "-"."$arframe\t";
										$seqcodon = $rsequence;
										
									}
									
									#print number of stop codons
									if($array3[0] == 1){
										print OUTPUT "^"."$array1[0]\t";
										
									}else{
										print OUTPUT "$array1[0]\t";
										
									}
									$fframe = $array2[0];
									
									#check the reading frame and see if the positive or reverse compliment of the sequence is used
								}
								
								
								
							}else{
								#there was an annotated reading frame
								#count number of stop codons in that reading frame
								#fstrand == 1 -> positive strand
								if($frame{$ulistchr[$i]}[$featurecheck[$j]] == 0 and $fstrand == 1){
									$fframe = 1;				#reading frame 0
									@codoninfo1 = &codoncounter(0, $sequence2, 0);
									$seqcodon = $sequence2;
									
									#print codon number for the reading frame
									if($codoninfo1[1] == 1){
										print OUTPUT "^"."$codoninfo1[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo1[0]\t";			#prints just the number of stop codons
										
									}
									
								}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 1 and $fstrand == 1){
									$fframe = 2;				#reading frame 1
									@codoninfo2 = &codoncounter(1, $sequence2, 0);
									$seqcodon = $sequence2;
									
									#print codon number for the reading frame
									if($codoninfo2[1] == 1){
										print OUTPUT "^"."$codoninfo2[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo2[0]\t";			#prints just the number of stop codons
										
									}
									
								}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 2 and $fstrand == 1){
									$fframe = 3;				#reading frame 2
									@codoninfo3 = &codoncounter(2, $sequence2, 0);
									$seqcodon = $sequence2;
									
									#print codon number for the reading frame
									if($codoninfo3[1] == 1){
										print OUTPUT "^"."$codoninfo3[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo3[0]\t";			#prints just the number of stop codons
										
									}
									
								}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 0 and $fstrand == 0){
									$fframe = 4;				#reading frame -1
									@codoninfo4 = &codoncounter(0, $rsequence, 0);
									$seqcodon = $rsequence;
									
									#print codon number for the reading frame
									if($codoninfo4[1] == 1){
										print OUTPUT "^"."$codoninfo4[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo4[0]\t";			#prints just the number of stop codons
										
									}
									
								}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 1 and $fstrand == 0){
									$fframe = 5;				#reading frame -2
									@codoninfo5 = &codoncounter(1, $rsequence, 0);
									$seqcodon = $rsequence;
									
									#print codon number for the reading frame
									if($codoninfo5[1] == 1){
										print OUTPUT "^"."$codoninfo5[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo5[0]\t";			#prints just the number of stop codons
										
									}
									
								}elsif($frame{$ulistchr[$i]}[$featurecheck[$j]] == 2 and $fstrand == 0){
									$fframe = 6;				#reading frame -3
									@codoninfo6 = &codoncounter(2, $rsequence, 0);
									$seqcodon = $rsequence;
									
									#print codon number for the reading frame
									if($codoninfo6[1] == 1){
										print OUTPUT "^"."$codoninfo6[0]\t";		#^ - denotes that the last codon was a stop and ignored for counter purposes
										
									}else{
										print OUTPUT "$codoninfo6[0]\t";			#prints just the number of stop codons
										
									}
									
								}
								
							}
							
							#
							#	Once correct sequence has been identified will pull out the correct codon and find the position in the codon
							#
							
							#check to see if feature has a negative reading frame (ie on the negative strand)
							#Will also give us the exact position of the SNP in the current sequence
							#the user SNP is always given as been on the positive strand
							
							#	$usnpuse has been deprecated since 30-04-13
							my $usnpuse;			#the position of the snp in the current sequence
							my $ncodon;
							
							#$usnpuse is the relative position of the snp to the start of the sequence
							#using the gtf annotation start positions
							my $act_len = 0;
							
							if($fframe == 4 or $fframe == 5 or $fframe == 6){
								
								#get rid of the +1 as the dist includes the offset position
								#change 29-04-13
								$usnpuse = (length($seqcodon)-$new_snp_dist);
								#new_snp_dist is the relative position of the SNP to the start of the sequence
								$ncodon = $rel_end_snp - 1;
								$act_len = $ncodon;
								
							}elsif($fframe == 1 or $fframe == 2 or $fframe == 3){
								
								#need to increment as substr offset == 0
								#change 29-04-13
								$usnpuse = $new_snp_dist + 1;
								$ncodon = $rel_start_snp - 1;
								$act_len = $ncodon;
							}
							
							my $cpos = codonpos($fframe, $usnpuse);			#will contain the position in the codon that the SNP appears
							
							#$seqcodon is the sequence to use
							#$fframe is the frame to use
							
							#extract the correct codon based on position in codon
							
							if($fstrand == 1){
								$seqcodon = "$stemp1"."$seqcodon"."$stemp2";
								$act_len += length($stemp1);
								
							}elsif($fstrand == 0){
								$seqcodon = "$rtemp2"."$seqcodon"."$rtemp1";
								$act_len += length($rtemp1);
							}
							
							
							#$usnpuse has been deprecated from this point on as of 30-04-13
							$usnpuse = (($usnpuse + 2)-1);
							#+2 as the stemp1 and stemp2 have been added; -1 as the off set is zero
							
							# 	$codon has been deprecated since 30-04-13
							#my $ucodon;								#this is the actual codon
							
							my $new_extract;
							
							if($cpos == 1){
								#$ucodon = substr($seqcodon, $usnpuse, 3);
								$new_extract = substr($seqcodon, $act_len, 3);
								
							}elsif($cpos == 2){
								#$ucodon = substr($seqcodon, ($usnpuse-1), 3);
								$act_len = $act_len - 1;
								$new_extract = substr($seqcodon, $act_len, 3);
								
								
							}elsif($cpos == 3){
								#$ucodon = substr($seqcodon, ($usnpuse-2), 3);
								$act_len = $act_len - 2;
								$new_extract = substr($seqcodon, $act_len, 3);
								
							}
							
							#print "usnp codon is $ucodon\tnew snp codon is $new_extract\t";
							
							#if reading frame is negative will need the reverse compliment of the user input mutation
							my $umut;								#the actual mutation that will be used
							
							#I think usnpchange should be ulistmut
							if($fframe == 4 or $fframe == 5 or $fframe == 6){
								$umut = reversecom($ulistmut[$i]);
								
							}else{
								$umut = $ulistmut[$i];
								
							}
							
							#	$ucodon has been deprecated since 30-04-13
							
							#	deprecated since 30-04-13
							#
							#my $sub1 = substr($ucodon, 0, 1);
							#my $sub2 = substr($ucodon, 1, 1);
							#my $sub3 = substr($ucodon, 2, 1);
							
							#break the codon into seperate nucleotides
							my $nsub1 = substr($new_extract, 0, 1);
							my $nsub2 = substr($new_extract, 1, 1);
							my $nsub3 = substr($new_extract, 2, 1);
							
							#	deprecated since 30-04-13
							#
							#my $codonold = $ucodon;
							
							my $codonold = $new_extract;
							my $codonnew;
							
							#print to output the relevant codon including the mutation in the correct position
							if($cpos == 1){
								#	deprecated since 30-04-13
								#print OUTPUT "[$sub1/"."$umut]"."$sub2"."$sub3\t";
								#$codonnew = "$umut"."$sub2"."$sub3";
								#print "new_amino: [$nsub1"."/$umut]"."$nsub2"."$nsub3\t";
								#print "old_amino: [$sub1"."/$umut]"."$sub2"."$sub3\n";
								
								print OUTPUT "[$nsub1"."/$umut]"."$nsub2"."$nsub3\t";
								$codonnew = "$umut"."$nsub2"."$nsub3";
								
							}elsif($cpos == 2){
								#	deprecated since 30-04-13
								#print OUTPUT "$sub1"."[$sub2/"."$umut]"."$sub3\t";
								#$codonnew = "$sub1"."$umut"."$sub3";
								#print "new_amino: $nsub1"."[$nsub2/"."$umut]"."$nsub3\t";
								#print "old_amino: $sub1"."[$sub2"."/$umut]"."$sub3\n";
								
								print OUTPUT "$nsub1"."[$nsub2/"."$umut]"."$nsub3\t";
								$codonnew = "$nsub1"."$umut"."$nsub3";
								
							}elsif($cpos == 3){
								#	deprecated since 30-04-13
								#print OUTPUT "$sub1"."$sub2"."[$sub3/"."$umut]\t";
								#$codonnew = "$sub1"."$sub2"."$umut";
								#print "new_amino: $nsub1"."$nsub2"."[$nsub3"."/$umut]\t";
								#print "old_amino: $sub1"."$sub2"."[$sub3"."/$umut]\n";
								
								print OUTPUT "$nsub1"."$nsub2"."[$nsub3"."/$umut]\t";
								$codonnew = "$nsub1"."$nsub2"."$umut";
								
							}
							
							#print to output the old and new amino acids
							print OUTPUT "["."$SLC{$codonold}"."/"."$SLC{$codonnew}"."]\t";
							
							
							if($SLC{$codonold} =~ m/X/g or $SLC{$codonnew} =~ m/X/g){
								print OUTPUT "NA\t";
								
							}elsif($SLC{$codonold} eq $SLC{$codonnew}){
								print OUTPUT "Y\t";
								
							}else{
								print OUTPUT "N\t";
								$summary_nonsyn_annotations++;
								
							}
							
							print OUTPUT "$protein_id{$ulistchr[$i]}[$featurecheck[$j]]\t";		#protein_id
							
							if($db_snp_file > 0){
								print OUTPUT "$db_snp_pos{$ulistchr[$i]}{$ulistpos[$i]}\t";
							}
							
						#	print "$db_snp_pos{$ulistchr[$i]}{$ulistpos[$i]}\n";
							
							print OUTPUT "\n";										#additional notes - there is none here
							
						}else{
							print OUTPUT "NA\t";
							print OUTPUT "NA\t";
							print OUTPUT "NA\t";
							print OUTPUT "NA\t";
							print OUTPUT "NA\t";
							print OUTPUT "NA\t";
							print OUTPUT "$protein_id{$ulistchr[$i]}[$j]\t";
							
							if($db_snp_file > 0){
								print OUTPUT "$db_snp_pos{$ulistchr[$i]}{$ulistpos[$i]}\t";
							}
							
							print OUTPUT "Warning user mutation ($ulistmut[$i]) is not a valid base\n";
						}
					}
				}
			}
		}
	}
}



#want to put any notes here!

#my $summary_gtf_feat = 0;			# Total number of features in the GTF - the number of lines
#my %summary_features;				# The available features in the GTF
#my %summary_chrs_in_GTF;			# Used to get total number of chromosomes in the GTF file
#my $summary_user_snps = 0;			# Total number of queried SNPs - the number of input lines
#my %summary_user_chrs;				# The queried chromosomes
#my $summary_annotated = 0;			# Number of SNPs annotated
#my $summary_snp_skipped = 0;		# The number of skipped SNPs - no GTF or incorrect input
#my $summary_in_feature_snps = 0;	# Number of annotated SNPs that were in a feature
#my $summary_nonsyn_annotations = 0;	# Number of non-synonymous annotations
#my $summary_noGTF_snp = 0;			# Number of SNPs with no GTF information
#my $summary_nofasta_snps = 0;		# Number of SNPs with no FASTA information

my @summary_user_chrs_list = keys %summary_user_chrs;
my $summary_user_chr_size = @summary_user_chrs_list;

print SUMMARY "Total number of queried SNPs: $summary_user_snps\n";
print SUMMARY "Total number of SNPs annotated: $summary_annotated\n";
print SUMMARY "Total number of queried chromosomes: $summary_user_chr_size\n";
print SUMMARY "Total number of skipped SNPs: $summary_snp_skipped\n";
print SUMMARY "Total number of annotated SNPs in a feature: $summary_in_feature_snps\n";
print SUMMARY "Total number of non-synonymous annotations (NOTE - a single SNP can have more than one non-synonymous annotation): $summary_nonsyn_annotations\n";
print SUMMARY "Total number of SNPs with no GTF information: $summary_noGTF_snp\n";
print SUMMARY "Total number of SNPs with no FASTA information: $summary_nofasta_snps\n\n";

# GTF information
#

my @summary_gtf_chrs_list = keys %summary_chrs_in_GTF;
my $summary_gtf_chr = @summary_gtf_chrs_list;
print SUMMARY "Chromosomes in the GTF file: $summary_gtf_chr\n";
print SUMMARY "Number of features in the GTF file $summary_gtf_feat\n";
print SUMMARY "Available features to annotate SNPs to:\n";

foreach my $key (keys %summary_features){
	print SUMMARY "\t$key\n";
}
print "\n";

close(SUMMARY);

sub gtfnf{
	my $size = @_;
	my $check;
	my $check2;
	my $dist1;
	my $dist2;
	my $chr;
	my $pos;
	my $mut;
	my $indel;
	my $intron;
	
	if($_[0] == 0){
		$check = $_[0];
		$check2 = $_[1];
		$dist1 = $_[2];
		$dist2 = $_[3];
		$chr = $_[4];
		$pos = $_[5];
		$mut = $_[6];
		$indel = $_[7];
		$intron = $_[8];
		
	}elsif($_[0] == 1){
		$check = $_[0];
		$check2 = $_[1];
		$dist1 = $_[2];
		$chr = $_[3];
		$pos = $_[4];
		$mut = $_[5];
		$indel = $_[6];
		$intron = $_[7];
		
	}elsif($_[0] == 2){
		$check = $_[0];
		$check2 = $_[1];
		$dist1 = $_[2];
		$chr = $_[3];
		$pos = $_[4];
		$mut = $_[5];
		$indel = $_[6];
		$intron = $_[7];
	}
	
	
	if($_[0] == 0){	#smallest dist to start and end are the same @_ to sub routine is 2,1,dist,chr,pos
		my $temp = ($pos+$dist1);
		my $tsize = @{$seqname{$chr}};
		for(my $y = 0; $y < $tsize; $y++){
			if($start{$chr}[$y] == $temp){
				$check++;
				print OUTPUT "$seqname{$chr}[$y]\t";											#user chr - note probably should be the listchr array
				print OUTPUT "$pos\t";													#user SNP position
				print OUTPUT "N\t";														#not within a feature
				
				if($intron == 0){
					print OUTPUT "Intronic\t"
				}else{
					print OUTPUT "Intergenic\t";
				}
				
				print OUTPUT "$dist1*^\t";												#distance to nearest feature
				print OUTPUT "$feature{$chr}[$y]\t";											#feature of SNP
				print OUTPUT "NA\t";													#number of different features the SNP is in
				print OUTPUT "NA\t";													#number of annotations for that feature
				print OUTPUT "$start{$chr}[$y]\t";											#start of current feature
				print OUTPUT "$end{$chr}[$y]\t";												#end of current feature
				print OUTPUT "$gene_id{$chr}[$y]\t";											#gene_id for the current feature
				print OUTPUT "$gene_name{$chr}[$y]\t";										#gene_name
				print OUTPUT "$transcript_id{$chr}[$y]\t";									#transcript_id
				print OUTPUT "$transcript_name{$chr}[$y]\t";									#transcript_name
				print OUTPUT "[$exon_number{$chr}[$y]"."/$gene_list{$gene_id{$chr}[$y]}]\t";		#annotated exon_number
				
				if($strand{$chr}[$y] ne "\."){
					print OUTPUT "$strand{$chr}[$y]\t";
				}else{
					print OUTPUT "$strand{$chr}[$y]**\t";
				}
				
				if($frame{$chr}[$y] eq "\."){
					print OUTPUT "NA\t";
					
				}elsif($frame{$chr}[$y] == 0 or $frame{$chr}[$y] == 1 or $frame{$chr}[$y] == 2){
				
					if($strand{$chr}[$y] eq "+"){							#positive strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "$rf_change\t";				#annotated reading frame
						
					}else{											#negative strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "-"."$rf_change\t";	
						
					}
				}
				
				print OUTPUT "NA\t";				#cant estimate reading frame
				print OUTPUT "NA\t";				#cant estimate number of stop codons
				print OUTPUT "NA\t";				#cant estimate codon
				print OUTPUT "NA\t";				#cant estimate amino acid
				print OUTPUT "NA\t";				#cant estimate synonymous
				print OUTPUT "$protein_id{$chr}[$y]\t";
		
		if($db_snp_file > 0){
			print OUTPUT "$db_snp_pos{$chr}{$pos}\t";
		}
		
				
				if($check2 == 0){
					print OUTPUT "There was no sequence information available for SNPs on this chromosome";				#no additional notes
					
				}
				
				if($mut =~ m/a|t|c|g/i and length($mut) == 1){
					if($indel == 1){
						print OUTPUT "; Warning the queried mutation is an indel\n";
					}else{
						print OUTPUT "\n";
					}
				}else{
					if($check2 == 0){
						print OUTPUT "; Warning user mutation ($mut) is not a valid base\n";
					}else{
						print OUTPUT "Warning user mutation ($mut) is not a valid base\n";
					}
				}
			}
			
		my $temp2 = ($pos-$dist2);
			if($end{$chr}[$y] == $temp2){
			
				print OUTPUT "$seqname{$chr}[$y]\t";											#user chr - note probably should be the listchr array
				print OUTPUT "$pos\t";													#user SNP position
				print OUTPUT "N\t";														#not within a feature
				
#				if($pos < $transcript_max{$transcript_id{$chr}[$y]} and $pos > $transcript_min{$transcript_id{$chr}[$y]}){
#					print OUTPUT "Intronic\t";
#				}else{
#					print OUTPUT "Intergenic\t";
#				}
				
				
				if($intron == 0){
					print OUTPUT "Intronic\t"
				}else{
					print OUTPUT "Intergenic\t";
				}
				
				
				print OUTPUT "$dist2*^\t";												#distance to nearest feature
				print OUTPUT "$feature{$chr}[$y]\t";											#feature of SNP
				print OUTPUT "NA\t";													#number of different features the SNP is in
				print OUTPUT "NA\t";													#number of annotations for that feature
				print OUTPUT "$start{$chr}[$y]\t";											#start of current feature
				print OUTPUT "$end{$chr}[$y]\t";												#end of current feature
				print OUTPUT "$gene_id{$chr}[$y]\t";											#gene_id for the current feature
				print OUTPUT "$gene_name{$chr}[$y]\t";										#gene_name
				print OUTPUT "$transcript_id{$chr}[$y]\t";									#transcript_id
				print OUTPUT "$transcript_name{$chr}[$y]\t";									#transcript_name
				print OUTPUT "[$exon_number{$chr}[$y]"."/$gene_list{$gene_id{$chr}[$y]}]\t";		#annotated exon_number
				
				if($strand{$chr}[$y] eq "\."){
					print OUTPUT "$strand{$chr}[$y]**\t";
					
				}else{
					print OUTPUT "$strand{$chr}[$y]\t";
					
				}
				
				if($frame{$chr}[$y] eq "\."){
					print OUTPUT "NA\t";
					
				}elsif($frame{$chr}[$y] == 0 or $frame{$chr}[$y] == 1 or $frame{$chr}[$y] == 2){
					
					if($strand{$chr}[$y] eq "+"){							#positive strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "$rf_change\t";				#annotated reading frame
						
					}else{											#negative strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "-"."$rf_change\t";	
						
					}
				}
				
				print OUTPUT "NA\t";				#cant estimate reading frame
				print OUTPUT "NA\t";				#cant estimate number of stop codons
				print OUTPUT "NA\t";				#cant estimate codon
				print OUTPUT "NA\t";				#cant estimate amino acid
				print OUTPUT "NA\t";				#cant estimate synonymous
				print OUTPUT "$protein_id{$chr}[$y]\t";
		
		if($db_snp_file > 0){
			print OUTPUT "$db_snp_pos{$chr}{$pos}\t";
		}
		
				
				if($check2 == 0){
					print OUTPUT "There was no sequence information available for SNPs on this chromosome";				#no additional notes
					
				}
				
				if($mut =~ m/a|t|c|g/i and length($mut) == 1){
					if($indel == 1){
						print OUTPUT "; Warning the queried mutation is an indel\n";
					}else{
						print OUTPUT "\n";
					}
				}else{
					if($check2 == 0){
						print OUTPUT "; Warning user mutation ($mut) is not a valid base\n";
					}else{
						print OUTPUT "Warning user mutation ($mut) is not a valid base\n";
					}
				}
				
			}
		}
	}elsif($_[0] == 1){	#dist to start is smallest @_ to sub routine is 2,1,dist,chr,pos
		my $tsize = @{$seqname{$chr}};
		for(my $y = 0; $y < $tsize; $y++){
			
			my $temp = ($pos+$dist1);
			
			if($start{$chr}[$y] == $temp){
				$check++;
				print OUTPUT "$seqname{$chr}[$y]\t";											#user chr - note probably should be the listchr array
				print OUTPUT "$pos\t";													#user SNP position
				print OUTPUT "N\t";														#not within a feature
				
#				if($pos < $transcript_max{$transcript_id{$chr}[$y]} and $pos > $transcript_min{$transcript_id{$chr}[$y]}){
#					print OUTPUT "Intronic\t";
#				}else{
#					print OUTPUT "Intergenic\t";
#				}
				
				
				if($intron == 0){
					print OUTPUT "Intronic\t"
				}else{
					print OUTPUT "Intergenic\t";
				}
				
				
				print OUTPUT "$dist1\t";												#distance to nearest feature
				print OUTPUT "$feature{$chr}[$y]\t";											#feature of SNP
				print OUTPUT "NA\t";													#number of different features the SNP is in
				print OUTPUT "NA\t";													#number of annotations for that feature
				print OUTPUT "$start{$chr}[$y]\t";											#start of current feature
				print OUTPUT "$end{$chr}[$y]\t";												#end of current feature
				print OUTPUT "$gene_id{$chr}[$y]\t";											#gene_id for the current feature
				print OUTPUT "$gene_name{$chr}[$y]\t";										#gene_name
				print OUTPUT "$transcript_id{$chr}[$y]\t";									#transcript_id
				print OUTPUT "$transcript_name{$chr}[$y]\t";									#transcript_name
				print OUTPUT "[$exon_number{$chr}[$y]"."/$gene_list{$gene_id{$chr}[$y]}]\t";		#annotated exon_number
				
				if($strand{$chr}[$y] ne "\."){
					print OUTPUT "$strand{$chr}[$y]\t";
				}else{
					print OUTPUT "$strand{$chr}[$y]**\t";
				}
				
				if($frame{$chr}[$y] eq "\."){
					print OUTPUT "NA\t";
					
				}elsif($frame{$chr}[$y] == 0 or $frame{$chr}[$y] == 1 or $frame{$chr}[$y] == 2){
					if($strand{$chr}[$y] eq "+"){							#positive strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "$rf_change\t";				#annotated reading frame
						
					}else{											#negative strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "-"."$rf_change\t";	
						
					}
				}
				
				print OUTPUT "NA\t";			#cant estimate reading frame
				print OUTPUT "NA\t";			#cant estimate number of stop codons
				print OUTPUT "NA\t";			#cant estimate codon
				print OUTPUT "NA\t";			#cant estimate amino acid
				print OUTPUT "NA\t";			#cant estimate synonymous
				print OUTPUT "$protein_id{$chr}[$y]\t";
		
		if($db_snp_file > 0){
			print OUTPUT "$db_snp_pos{$chr}{$pos}\t";
		}
		
				
				if($check2 == 0){
					print OUTPUT "There was no sequence information available for SNPs on this chromosome";				#no additional notes
				}
				
				if($mut =~ m/a|t|c|g/i and length($mut) == 1){
					if($indel == 1){
						print OUTPUT "; Warning the queried mutation is an indel\n";
					}else{
						print OUTPUT "\n";
					}
				}else{
					if($check2 == 0){
						print OUTPUT "; Warning user mutation ($mut) is not a valid base\n";
					}else{
						print OUTPUT "Warning user mutation ($mut) is not a valid base\n";
					}
				}				
			}
		}
		
	}elsif($_[0] == 2){	#dist to end is smallest @_ to sub routine is 2,1,dist,chr,pos
		
		my $temp = ($pos-$dist1);
		my $tsize = @{$seqname{$chr}};
		
		for(my $y = 0; $y < $tsize; $y++){
			
			if($end{$chr}[$y] == $temp){
				
				$check++;
				print OUTPUT "$seqname{$chr}[$y]\t";											#user chr - note probably should be the listchr array
				print OUTPUT "$pos\t";													#user SNP position
				print OUTPUT "N\t";														#not within a feature
				
#				if($pos < $transcript_max{$transcript_id{$chr}[$y]} and $pos > $transcript_min{$transcript_id{$chr}[$y]}){
#					print OUTPUT "Intronic\t";
#				}else{
#					print OUTPUT "Intergenic\t";
#				}
				
				
				if($intron == 0){
					print OUTPUT "Intronic\t"
				}else{
					print OUTPUT "Intergenic\t";
				}
				
				
				print OUTPUT "$dist1\t";												#distance to nearest feature
				print OUTPUT "$feature{$chr}[$y]\t";											#feature of SNP
				print OUTPUT "NA\t";													#number of different features the SNP is in
				print OUTPUT "NA\t";													#number of annotations for that feature
				print OUTPUT "$start{$chr}[$y]\t";											#start of current feature
				print OUTPUT "$end{$chr}[$y]\t";												#end of current feature
				print OUTPUT "$gene_id{$chr}[$y]\t";											#gene_id for the current feature
				print OUTPUT "$gene_name{$chr}[$y]\t";										#gene_name
				print OUTPUT "$transcript_id{$chr}[$y]\t";									#transcript_id
				print OUTPUT "$transcript_name{$chr}[$y]\t";									#transcript_name
				print OUTPUT "[$exon_number{$chr}[$y]"."/$gene_list{$gene_id{$chr}[$y]}]\t";		#annotated exon_number
				
				#check to see if strand has a value
				if($strand{$chr}[$y] ne "\."){
					print OUTPUT "$strand{$chr}[$y]\t";
					
				}else{
					print OUTPUT "$strand{$chr}[$y]**\t";
					
				}
				
				if($frame{$chr}[$y] eq "\."){
					print OUTPUT "NA\t";
					
				}elsif($frame{$chr}[$y] == 0 or $frame{$chr}[$y] == 1 or $frame{$chr}[$y] == 2){
					
					if($strand{$chr}[$y] eq "+"){							#positive strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "$rf_change\t";				#annotated reading frame
						
					}else{											#negative strand
						my $rf_change = ($frame{$chr}[$y]+1);
						print OUTPUT "-"."$rf_change\t";	
						
					}
					
				}
				
				print OUTPUT "NA\t";				#cant estimate reading frame
				print OUTPUT "NA\t";				#cant estimate number of stop codons
				print OUTPUT "NA\t";				#cant estimate codon
				print OUTPUT "NA\t";				#cant estimate amino acid
				print OUTPUT "NA\t";				#cant estimate synonymous
				print OUTPUT "$protein_id{$chr}[$y]\t";
		
		if($db_snp_file > 0){
			print OUTPUT "$db_snp_pos{$chr}{$pos}\t";
		}
		
				
				if($check2 == 0){
					print OUTPUT "There was no sequence information available for SNPs on this chromosome";				#no additional notes
					
				}
				
				if($mut =~ m/a|t|c|g/i and length($mut) == 1){
					if($indel == 1){
						print OUTPUT "; Warning the queried mutation is an indel\n";
					}else{
						print OUTPUT "\n";
					}
				}else{
					if($check2 == 0){
						print OUTPUT "; Warning user mutation ($mut) is not a valid base\n";
					}else{
						print OUTPUT "Warning user mutation ($mut) is not a valid base\n";
					}
				}				
			}
		}
	}

}

#sub to count the number of codons in the sequence - returns a 1 or 0 as last value to tell you if the last codon was a stop. 1 - yes, 0 - no
sub codoncounter{
	my ($rf, $seq, $test) = @_;
	
	#get length of current reading frame
	my $length = (length($seq)-$rf);
	
	#check divisibility of sequence length by 3
	my $rem = $length%3;
	
	#if there is a remainder concatonate ?'s to end to ensure divisibility
	if($rem == 1){
		$seq = $seq."??";
		
	}elsif($rem == 2){
		$seq = $seq."?";
		
	}
	
	#get length of current reading frame with necessary ?'s attached
	my $dist = (length($seq)-$rf);
	my $stopcount = 0;								#counter to count the number of stop codons
	my $laststopcount = 0;							#set to zero to say that last codon is not a stop
	my $lastcodon;									#following next for loop will store the last codon in the sequence
	
	#iterate through sequence in 3's counting stop codons in the sequence
	for(my $i = $rf; $i < ($dist+$rf); $i+=3){
		my $codon = substr($seq, $i, 3);
		
		#count number of stop codons - start reading from rf value and to end of sequence (including ?'s when necessary)
		if($codon eq "TAA" or $codon eq "TGA" or $codon eq "TAG"){
			$stopcount++;
		}
		
		#store the value of the codon at each iteration (overwriting each time)
		#this will then keep just the last codon in memory (outside the for loop)
		$lastcodon = $codon;
	}
	
	if($lastcodon eq "TAA" or $lastcodon eq "TGA" or $lastcodon eq "TAG"){
		$laststopcount++;		#tell the script that the last codon in the sequence was a stop
		$stopcount--;			#decrement the total number of stops in the sequence by 1 because the last codon is a stop (we want to ignore this stop)
	}

	return($stopcount, $laststopcount);

}

#sub to get the reverse compliment of a given sequence
sub reversecom{
	my ($seqvar1) = @_;

	my $seqvar = reverse($seqvar1);
	
	$seqvar =~ s/a/1/ig;
	$seqvar =~ s/t/2/ig;
	$seqvar =~ s/c/3/ig;
	$seqvar =~ s/g/4/ig;
	
	$seqvar =~ s/1/T/ig;
	$seqvar =~ s/2/A/ig;
	$seqvar =~ s/3/G/ig;
	$seqvar =~ s/4/C/ig;

	return $seqvar;
}

#sub to find the locations of features that contain the SNP 
#returns a value to say if that SNP was found in a feature followed the feature locations that it was found in
sub featurecheck{
	my ($chr, $snppos) = @_;
	my $check = 0;				#check to ensure SNPs are in a feature
	my @featureindex;
	my $ssize = @{$seqname{$chr}};

	#find the index positions that satisfy the search (that mean the SNP is within a feature)
	for(my $i = 0; $i < $ssize; $i++){
		if($snppos >= $start{$chr}[$i] && $snppos <= $end{$chr}[$i]){
			$check++;
			push(@featureindex, $i);
		}
	}
	
	#return check and a list of gtf index positions that indicate the SNP is within a feature (at that index position)
	#check will be zero if the SNP is not within any feature
	return($check, @featureindex);

}

#sub to get the position in the codon that the SNP occurs in
sub codonpos{
	my ($rf, $position) = @_;
	
	my $codonrem = ($position % 3);
	
	my $snpcodonpos;
	
	if($rf == 1 or $rf == 4){
		if($codonrem == 0){
			$snpcodonpos = 3;
			
		}elsif($codonrem == 1){
			$snpcodonpos = 1;
		
		}elsif($codonrem == 2){
			$snpcodonpos = 2;
		
		}
		
	}elsif($rf == 2 or $rf == 5){
		if($codonrem == 0){
			$snpcodonpos = 2;
		
		}elsif($codonrem == 1){
			$snpcodonpos = 3;
		
		}elsif($codonrem == 2){
			$snpcodonpos = 1;
		
		}
		
	}elsif($rf == 3 or $rf == 6){
		if($codonrem == 0){
			$snpcodonpos = 1;
		
		}elsif($codonrem == 1){
			$snpcodonpos = 2;
		
		}elsif($codonrem == 2){
			$snpcodonpos = 3;
		
		}
		
	}

	return $snpcodonpos;

}


sub help{

print "\n$version\n\n";

if($_[0] == 1){
	print "\nError input file name ambiguous\nEither not set or declared twice!\n\n";
	print "----------------------------------------------------------\n";
}elsif($_[0] == 2){
	print "\nError GTF file name ambiguous\nEither not set or declared twice!\n\n";
	print "----------------------------------------------------------\n";
}elsif($_[0] == 3){
	print "\nError FASTA file name ambiguous\nEither not set or declared twice!\n\n";
	print "----------------------------------------------------------\n";
}

print STDOUT <<END;

SNPdat is a high throughput analysis tool that can provide a comprehensive annotation of both novel and known single nucleotide polymorphisms (SNPs).

SNPdat requires that each file is specified when running the program. There are 3 mandatory file definitions.

Usage:
perl SNPdat -i Input_file -f Fasta_file -g Gene_Transfer_File

Required:

-i	Input file
-g	Gene transfer file (GTF)
-f	FASTA formated sequence file

Optional:

-d	a dbSNP ASN_FLAT file processed using SNPdat_parse_dbsnp.pl (optional)
-s	a file containing a summary of the queried SNPs (optional)
	NOTE:If no output file is specified, results will be printed to 'Input_file.summary'
	
-o	output_file specified by the user (optional)
	NOTE:If no output file is specified, results will be printed to 'Input_file.output'
	
Advanced:
	
-x	retrieve sequence information from the next/previous feature should a codon cross that boundary.
	
	User can specify a comma separated list of features from the GTF. This is case-sensitive.
	This is only recommended for advanced users who understand what it does.
	
	By default this is not set. See website/manual for more information.
	
	USAGE:
	-x feature1,feature2
	
	e.g.
	-x exon
	-x CDS
	-x exon,CDS
	
Info:

-h	This wonderful help page
-v	This version of SNPdat

For more instuctions see the SNPdat webage:
http://code.google.com/p/snpdat/

END

}


print "SNPdat finished analysing all SNPs:\n";
system("date");
print "view file '$snpdat_files{-o}' for results\n";

