PrimerMapper
Author: Damien O'Halloran, The George Washington University, 2015

PrimerMapper

Provides the following: 
1. a GUI to facilitate primer design and visualization
2. Traditional primer design as well as SNP and allele specific primer design
3. Visualization of primer distribution across each sequence and entire input file
4. Returns primers specific to the entire input file using specificty and mis-match options 
5. Remote sequence access from GenBank and dbSNP
6. Primer BLAST facility against multiple NCBI databases 
7. Generates primer dimer scores for all primers generated to facilitate multiplex PCR expts
8. Provides browser visualization of primer maps and permits the user to draw new primers

To run, execute as follows: 
>perl PrimerMapper_driver.pl 
## --> must have web_blast.pl script in PATH in order to perform remote BLAST at NCBI
## --> must have cols_to_rows.pl script in PATH in order to visualize map in browser


To install:

1. Download and extract the PrimerMapper.zip file
2. >tar -xzvf PrimerMapper.zip 
3. The extracted dir will be called PrimerMapper  
4. cd PrimerMapper
5. to install:
   perl Makefile.PL
   make 
   make test
   make install  
   
## --> PC platform may require 'dmake' or 'nmake' 
   
   
PrimerMapper Defaults

default settings are as follows:

####Primer Design Defaults
my $five_prime_end  = "150";
my $three_prime_end = "150";
my $kmer_max        = "25";
my $kmer_min        = "18";
my $clamp           = "N";
my $higher_gc       = "60";
my $lower_gc        = "40";
my $upper_tm        = "68";
my $lower_tm        = "55";
my $spec            = "N";
my $mis             = "2";
my $salt            = "50";
my $DNA_conc        = "100";
my $selfie_cuttof   = "9";
my $outputfile   = "primers.tsv";

####SNP Primer Design Defaults
my $five_prime_SNP    = "120";
my $three_prime_SNP   = "120";
my $kmer_max_SNP      = "25";
my $kmer_min_SNP      = "18";
my $clamp_SNP         = "N";
my $higher_gc_SNP     = "60";
my $lower_gc_SNP      = "40";
my $upper_tm_SNP      = "68";
my $lower_tm_SNP      = "55";
my $spec_SNP          = "N";
my $mis_SNP           = "2";
my $salt_SNP          = "50";
my $DNA_conc_SNP      = "100";
my $selfie_cuttof_SNP = "9";
my $outputfile_SNP = "SNP_primers.tsv";
my $snp_distance      = "30";

####Other Defaults
my $min_size_seq    = "300";
my $format          = ".png";

the defaults will be overwritten if a parameter is changed
   


