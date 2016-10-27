#
# Main module for PrimerMapper
#
# Please direct questions and support issues to <https://github.com/dohalloran/PrimerMapper/issues>
#
# Author: Damien O'Halloran, The George Washington University, 2015
#
# GNU GENERAL PUBLIC LICENSE
#
# POD documentation before the code

=head1 NAME
PrimerMapper - runs the primer design side of things

=head2 SYNOPSIS
  use Tkx;
  use LWP::Simple;
  use Bio::SeqIO;
  use PrimerMapperUtilities;
  use PrimerMapperGraphics;
  
  for Usage, run by creating a PrimerMapper object within PrimerMapper_driver.pl
  my $tmp = PrimerMapper->new();
  
Primer Design Defaults
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
    my $outputfile      = "primers.tsv";
    my $min_size_seq    = "300";
    my $format          = ".png";
    my $repetition      = "N";
SNP Primer Design Defaults
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
    my $outputfile_SNP    = "SNP_primers.tsv";
    my $snp_distance      = "30";
Other Defaults
    my $min_size_seq    = "300";
    my $format          = ".png";

the defaults will be overwritten if a parameter is changed

=head3 DESCRIPTION
This object produces the following:

1. a GUI to facilitate primer design and visualization
2. Provides browser visualization of primer maps and permits the user to draw new primers
3. Traditional primer design as well as SNP and allele specific primer design
4. Visualization of primer distribution across each sequence and entire input file
5. Returns primers specific to the entire input file using specificty and mis-match options 
6. Remote sequence access from GenBank and dbSNP
7. Primer BLAST facility against multiple NCBI databases 
8. Generates primer dimer scores for all primers generated to facilitate multiplex PCR expts

=head4 REFERENCES
The design of the code and it's execution relies heaviliy upon the BioPerl toolkit:
Stajich, J. E. et al. The Bioperl toolkit: Perl modules for the life sciences. Genome
Res. 12, 1611-1618 (2002)

=head5 Support 
All contributions are welcome

=head6 Reporting Bugs
Report bugs to the PrimerMapper bug tracking system to help keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:
  https://github.com/dohalloran/PrimerMapper/issues
  
=head7 APPENDIX
The rest of the documentation details each of the subs and code

=cut


package PrimerMapper;

#############################################################
#
# PrimerMapper
# Written by Damien O'Halloran
# The George Washington University
# Summer 2015
#
#############################################################

use warnings;
use strict;

#### Load Modules
use Tkx;
use LWP::Simple;
use Cwd;
use Bio::SeqIO;
use PrimerMapperUtilities;
use PrimerMapperGraphics;
####

#########################
#########################
#########################

our $VERSION = '2.0';

#########################
#########################
sub new {
    my $class = shift;
    return bless {}, $class;
}
#########################

####Global Primer Design Defaults
my $fasta;
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
my $outputfile      = "primers.tsv";
my $min_size_seq    = "300";
my $format          = ".png";
my $repetition      = "N";
###########################

####Global Strings, Arrays, and Hashes
my $outfile;
my $outputfile_html;
my $out_single = "single_viewer.txt";
my $specificty;
my $seq_in = "NM_001307521.1,GI: 228485043";
my @array_length;
my @array_name;
my @gene_length;
my @increment;
my @primers;
my $foo;
my $Tm2;
my $Tm;
my $kmer;
my $id;
my $out_image;
my $kmer_diff = $kmer_max - $kmer_min;
my $blaster;
###########################

####Global SNP Primer Design Defaults
my $fasta_SNP;
my $fasta_SNP_NEW;
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
my $outputfile_SNP    = "SNP_primers.tsv";
my $snp_distance      = "30";
###########################

####Global SNP Strings and Hashes
my $kmer_diff_SNP = $kmer_max_SNP - $kmer_min_SNP;
my $kmer_SNP;
my $snp;
my $outfile_SNP;
my $out_single_SNP = "single_viewer_SNP.txt";
my $specificty_SNP;
my $seq_in_SNP = "rs744373,rs11136000,rs3764650";
my @array_length_SNP;
my @array_name_SNP;
my @gene_length_SNP;
my @increment_SNP;
my @SNP_position;
my @SNP_primers;
my $foo_SNP;
my $Tm2_SNP;
my $Tm_SNP;
my $id_SNP;
my $out_image_SNP;
###########################

####allele specific primer variables
my $Tm3_SNP;
my $Tm4_SNP;
my $outputfile_SNP_AS = "Allele_specific_primers.tsv";
my $degen;
###########################

####################################
####Interface Begins####
####################################
=head1 Load the Tkx interface
 Title   :  Tkx::MainLoop();
 Usage   :  my $main_window = Tkx::widget->new(".");
 Function:  Loads the Tkx interface 
 Returns :  GUI 
=cut

my $main_window = Tkx::widget->new(".");
$main_window->g_wm_title("PrimerMapper");
my $outer_frame = $main_window->new_ttk__frame( -padding => "3 3 12 12" );
$outer_frame->g_grid( -column => 0, -row => 0, -sticky => "nwes" );

my $input_frame = $outer_frame->new_ttk__labelframe( -text => "INPUT" );
$input_frame->new_ttk__frame( -padding => "3 3 12 12" );
$input_frame->g_grid(
    -column     => 0,
    -columnspan => 5,
    -rowspan    => 1,
    -row        => 0,
    -sticky     => "nwes"
);

my $input_ = $input_frame->new_ttk__entry( -width => 40, -textvariable => \$fasta );
$input_->g_grid( -column => 1, -row => 0, -sticky => "we" );
my $load_file_ = $input_frame->new_ttk__button( -text => "Load File", -command => \&loadfile );
$load_file_->g_grid( -column => 0, -row => 0, -sticky => "w" );
my $load_seq_ = $input_frame->new_ttk__entry( -width => 40, -textvariable => \$seq_in );
$load_seq_->g_grid( -column => 4, -row => 0, -sticky => "w" );
my $get_seqs_ = $input_frame->new_ttk__button(
    -text    => "Get Sequences",
    -command => sub { get_sequences(); }
);
$get_seqs_->g_grid( -column => 2, -row => 0, -sticky => "w" );

########## SNP TOP
my $outer_frame_snp = $outer_frame->new_ttk__labelframe( -text => "SNP INPUT" );
$outer_frame_snp->new_ttk__frame( -padding => "3 3 12 12" );
$outer_frame_snp->g_grid(
    -column     => 0,
    -columnspan => 5,
    -rowspan    => 1,
    -row        => 14,
    -sticky     => "nwes",
    -pady       => 10
);

my $input_snp_ = $outer_frame_snp->new_ttk__entry( -width => 40, -textvariable => \$fasta_SNP );
$input_snp_->g_grid( -column => 1, -row => 14, -sticky => "w" );
my $load_snp_ =
  $outer_frame_snp->new_ttk__button( -text => "Load SNP File", -command => \&loadfile_SNP );
$load_snp_->g_grid( -column => 0, -row => 14, -sticky => "w" );
my $load_seq_snp_ = $outer_frame_snp->new_ttk__entry( -width => 40, -textvariable => \$seq_in_SNP );
$load_seq_snp_->g_grid( -column => 4, -row => 14, -sticky => "w" );
my $get_snp_seqs_ = $outer_frame_snp->new_ttk__button(
    -text    => "Get SNP Sequences",
    -command => sub { get_sequences_SNP(); }
);
$get_snp_seqs_->g_grid( -column => 2, -row => 14, -sticky => "w" );

my $design_frame_ = $outer_frame->new_ttk__labelframe( -text => "PRIMER DESIGN" );
$design_frame_->new_ttk__frame( -padding => "3 3 12 12" );
$design_frame_->g_grid(
    -column     => 0,
    -columnspan => 2,
    -rowspan    => 1,
    -row        => 16,
    -sticky     => "nwes"
);

my $five_prime_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$five_prime_SNP
);
$five_prime_snp_->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $three_prime_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$three_prime_SNP
);
$three_prime_snp_->g_grid( -column => 5, -row => 2, -sticky => "we" );

my $kmer_max_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$kmer_max_SNP );
$kmer_max_snp_->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $kmer_min_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$kmer_min_SNP );
$kmer_min_snp_->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $upper_gc_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$higher_gc_SNP );
$upper_gc_snp_->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $lower_gc_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$lower_gc_SNP );
$lower_gc_snp_->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $upper_tm_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$upper_tm_SNP );
$upper_tm_snp_->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $lower_tm_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$lower_tm_SNP );
$lower_tm_snp_->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $gc_clamp_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$clamp_SNP );
$gc_clamp_snp_->g_grid( -column => 2, -row => 6, -sticky => "e" );

my $salt_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$salt_SNP );
$salt_snp_->g_grid( -column => 2, -row => 7, -sticky => "we" );

my $specific_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$spec_SNP );
$specific_snp_->g_grid( -column => 5, -row => 6, -sticky => "w" );

my $mis_match_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$mis_SNP );
$mis_match_snp_->g_grid( -column => 5, -row => 7, -sticky => "w" );

my $dna_concentration_snp_ = $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$DNA_conc_SNP );
$dna_concentration_snp_->g_grid( -column => 7, -row => 2, -sticky => "we" );

my $outfile_snp_ =
  $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$outputfile_SNP );
$outfile_snp_->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $self_comp_snp_ =
  $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$selfie_cuttof_SNP );
$self_comp_snp_->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $format_snp_ = $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$format );
$format_snp_->g_grid( -column => 7, -row => 6, -sticky => "we" );

my $snp_distance_ =
  $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$snp_distance );
$snp_distance_->g_grid( -column => 7, -row => 5, -sticky => "we" );

my $min_size_snp_ =
  $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$min_size_seq );
$min_size_snp_->g_grid( -column => 7, -row => 7, -sticky => "we" );

my $distance_from_snp_ = $design_frame_->new_ttk__label( -text => "5' up from SNP" );
$distance_from_snp_->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $distance_down_from_snp_ = $design_frame_->new_ttk__label( -text => "3' down from SNP" );
$distance_down_from_snp_->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $max_primer_len_snp_ = $design_frame_->new_ttk__label( -text => "Primer length max" );
$max_primer_len_snp_->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $min_primer_len_snp_ = $design_frame_->new_ttk__label( -text => "Primer length min" );
$min_primer_len_snp_->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $max_gc_snp_ = $design_frame_->new_ttk__label( -text => "Upper GC%" );
$max_gc_snp_->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $min_gc_snp_ = $design_frame_->new_ttk__label( -text => "Lower GC%" );
$min_gc_snp_->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $max_tm_snp_ = $design_frame_->new_ttk__label( -text => "Upper Tm" );
$max_tm_snp_->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $min_tm_snp_ = $design_frame_->new_ttk__label( -text => "Lower Tm" );
$min_tm_snp_->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $clamp_snp_ = $design_frame_->new_ttk__label( -text => "GC clamp (Y/N)" );
$clamp_snp_->g_grid( -column => 1, -row => 6, -sticky => "w" );

my $sequence_specific_snp_ = $design_frame_->new_ttk__label( -text => "input specificity (Y/N)" );
$sequence_specific_snp_->g_grid( -column => 4, -row => 6, -sticky => "w" );

my $sequence_mismatch_snp_ = $design_frame_->new_ttk__label( -text => "mis-matches" );
$sequence_mismatch_snp_->g_grid( -column => 4, -row => 7, -sticky => "w" );

my $nacl_concentration_snp_ = $design_frame_->new_ttk__label( -text => "salt concentration (mM)" );
$nacl_concentration_snp_->g_grid( -column => 1, -row => 7, -sticky => "w" );

my $dna_snp_ = $design_frame_->new_ttk__label( -text => "DNA concentration (nM)" );
$dna_snp_->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $snp_outfile_ = $design_frame_->new_ttk__label( -text => "SNP Primer file" );
$snp_outfile_->g_grid( -column => 6, -row => 3, -sticky => "w" );

my $self_snp_ = $design_frame_->new_ttk__label( -text => "Self-complementarity" );
$self_snp_->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $min_basepairs_from_snp_ = $design_frame_->new_ttk__label( -text => "Minimum distance from SNP" );
$min_basepairs_from_snp_->g_grid( -column => 6, -row => 5, -sticky => "w" );

my $outfile_format_snp_ = $design_frame_->new_ttk__label( -text => "Graphic format (.png or .gif)" );
$outfile_format_snp_->g_grid( -column => 6, -row => 6, -sticky => "w" );

my $min_input_size_snp_ = $design_frame_->new_ttk__label( -text => "Minimum size sequence" );
$min_input_size_snp_->g_grid( -column => 6, -row => 7, -sticky => "w" );

foreach ( Tkx::SplitList( $design_frame_->g_winfo_children ) ) {
    Tkx::grid_configure( $_, -padx => 5, -pady => 5 );
}

my $run_snp_frame_ =
  $outer_frame->new_ttk__labelframe( -text => "RUN", -width => 10, -height => 60 );
$run_snp_frame_->new_ttk__frame( -padding => "3 3 12 12" );
$run_snp_frame_->g_grid( -column => 0, -row => 28, -sticky => "nwes" );


my $design_snp_frame_ = $run_snp_frame_->new_ttk__button(
    -text    => "1: Design SNP Primers and Graphics",
    -command => sub {
        calculate_SNP();
        graphics_all_primers_SNP( \@array_length_SNP, \@array_name_SNP,
            \@gene_length_SNP, $format, $snp );
        graphics_single_view_SNP( \@gene_length_SNP, $format, $out_single_SNP );
    }
);
$design_snp_frame_->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $multiple_snp_ = $run_snp_frame_->new_ttk__button(
    -text    => "2: Multiplex PCR dimer scores",
    -command => sub { primer_dimer_SNP(); }
);
$multiple_snp_->g_grid( -column => 4, -row => 9, -sticky => "s" );

my $clean_up_snp_ = $run_snp_frame_->new_ttk__button(
    -text    => "3: Clean-up (run last)",
    -command => sub { clean_up(); }
);
$clean_up_snp_->g_grid( -column => 5, -row => 9, -sticky => "s" );

########## SNP BTM

my $out_frame_primers = $outer_frame->new_ttk__labelframe( -text => "PRIMER DESIGN" );
$out_frame_primers->new_ttk__frame( -padding => "3 3 12 12" );
$out_frame_primers->g_grid(
    -column     => 0,
    -columnspan => 2,
    -rowspan    => 1,
    -row        => 2,
    -sticky     => "nwes"
);

my $five_primer_ending_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$five_prime_end
);
$five_primer_ending_->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $three_primer_ending_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$three_prime_end
);
$three_primer_ending_->g_grid( -column => 5, -row => 2, -sticky => "we" );

my $max_primer_size_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$kmer_max );
$max_primer_size_allowed_->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $min_primer_size_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$kmer_min );
$min_primer_size_allowed_->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $max_gc_content_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$higher_gc );
$max_gc_content_allowed_->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $min_gc_content_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$lower_gc );
$min_gc_content_allowed_->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $max_primer_tm_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$upper_tm );
$max_primer_tm_allowed_->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $min_primer_tm_allowed_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$lower_tm );
$min_primer_tm_allowed_->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $add_gc_clamp_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$clamp );
$add_gc_clamp_->g_grid( -column => 2, -row => 6, -sticky => "e" );

my $sodium_chloride_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$salt );
$sodium_chloride_->g_grid( -column => 2, -row => 7, -sticky => "we" );

my $input_specific__ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$spec );
$input_specific__->g_grid( -column => 5, -row => 6, -sticky => "w" );

my $permit_mis_matches_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$mis );
$permit_mis_matches_->g_grid( -column => 5, -row => 7, -sticky => "w" );

my $dna_concentr_ = $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$DNA_conc );
$dna_concentr_->g_grid( -column => 7, -row => 2, -sticky => "we" );

my $outfile__ = $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$outputfile );
$outfile__->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $self_complement_score_ =
  $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$selfie_cuttof );
$self_complement_score_->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $format__ = $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$format );
$format__->g_grid( -column => 7, -row => 5, -sticky => "we" );

my $minimum_sequence_size_allowed_ =
  $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$min_size_seq );
$minimum_sequence_size_allowed_->g_grid( -column => 7, -row => 6, -sticky => "we" );

my $repetitive_seq_allowed_ =
  $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$repetition );
$repetitive_seq_allowed_->g_grid( -column => 7, -row => 7, -sticky => "we" );

my $run_frame__ =
  $outer_frame->new_ttk__labelframe( -text => "RUN", -width => 10, -height => 40 );
$run_frame__->new_ttk__frame( -padding => "3 3 12 12" );
$run_frame__->g_grid( -column => 0, -row => 8, -sticky => "nwes" );

my $blast_frame_ = $main_window->new_ttk__labelframe(
    -text => "BLAST PRIMERS - insert primers here in fastA format" );
$blast_frame_->new_ttk__frame( -padding => "3 3 12 12" );
$blast_frame_->g_grid( -column => 0, -row => 10, -sticky => "nwes" );

my $text = $blast_frame_->new_tk__text( -width => 80, -height => 6, -wrap => "none" );
$text->g_grid;
my $thetext;

my $exit_ = $main_window->new_ttk__button(
    -text    => "EXIT",
    -command => sub { exit_program(); }
);
$exit_->g_grid( -column => 0, -row => 14, -sticky => "s" );

my $which_db_ = $main_window->new_ttk__labelframe( -text => "SELECT DATABASE:" );
$which_db_->g_grid( -column => 0, -row => 12, -sticky => "nwes" );

my $nt_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "Nucleotide collection (nt)",
    -variable => \$blaster,
    -value    => "nt"
);
my $human_gDNA_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "Genomic Human",
    -variable => \$blaster,
    -value    => "human_genomic"
);
my $other_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "Genomic others",
    -variable => \$blaster,
    -value    => "other_genomic"
);
my $est_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "EST others",
    -variable => \$blaster,
    -value    => "est_others"
);
my $est_human_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "EST Human",
    -variable => \$blaster,
    -value    => "est_human"
);
my $est_mouse_db_ = $which_db_->new_ttk__radiobutton(
    -text     => "EST Mouse",
    -variable => \$blaster,
    -value    => "est_mouse"
);

$nt_db_->g_grid( -column => 0, -row => 13, -sticky => "w" );
$human_gDNA_db_->g_grid( -column => 1, -row => 13, -sticky => "w" );
$other_db_->g_grid( -column => 2, -row => 13, -sticky => "w" );
$est_db_->g_grid( -column => 3, -row => 13, -sticky => "w" );
$est_human_db_->g_grid( -column => 4, -row => 13, -sticky => "w" );
$est_mouse_db_->g_grid( -column => 5, -row => 13, -sticky => "w" );

my $blast_run_ = $which_db_->new_ttk__button(
    -text    => "BLAST",
    -command => sub { web_blast_ncbi(); }
);
$blast_run_->g_grid( -column => 0, -row => 14, -sticky => "w" );

my $get_the_primers_ = $run_frame__->new_ttk__button(
    -text    => "1: Design Primers and Graphics (run this first)",
    -command => sub {
        calculate();
        graphics_all_primers( \@array_length, \@array_name, \@gene_length,
            $format );
        graphics_single_view( \@gene_length, $format, $out_single );
    }
);
$get_the_primers_->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $multiplex__ = $run_frame__->new_ttk__button(
    -text    => "2: Multiplex PCR dimer scores",
    -command => sub { primer_dimer(); }
);
$multiplex__->g_grid( -column => 4, -row => 9, -sticky => "s" );

my $cleaner__ = $run_frame__->new_ttk__button(
    -text    => "3: Clean-up (run last)",
    -command => sub { clean_up(); }
);
$cleaner__->g_grid( -column => 5, -row => 9, -sticky => "s" );

my $abc = $out_frame_primers->new_ttk__label( -text => "5' search area" );
$abc->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $bca = $out_frame_primers->new_ttk__label( -text => "3' search area" );
$bca->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $cba = $out_frame_primers->new_ttk__label( -text => "Primer length max" );
$cba->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $cab = $out_frame_primers->new_ttk__label( -text => "Primer length min" );
$cab->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $def = $out_frame_primers->new_ttk__label( -text => "Upper GC%" );
$def->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $input_g = $out_frame_primers->new_ttk__label( -text => "Lower GC%" );
$input_g->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $hij = $out_frame_primers->new_ttk__label( -text => "Upper Tm" );
$hij->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $ijk = $out_frame_primers->new_ttk__label( -text => "Lower Tm" );
$ijk->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $jkl = $out_frame_primers->new_ttk__label( -text => "GC clamp (Y/N)" );
$jkl->g_grid( -column => 1, -row => 6, -sticky => "w" );

my $klm = $out_frame_primers->new_ttk__label( -text => "input specificity (Y/N)" );
$klm->g_grid( -column => 4, -row => 6, -sticky => "w" );

my $lmn = $out_frame_primers->new_ttk__label( -text => "mis-matches" );
$lmn->g_grid( -column => 4, -row => 7, -sticky => "w" );

my $lzn = $out_frame_primers->new_ttk__label( -text => "salt concentration (mM)" );
$lzn->g_grid( -column => 1, -row => 7, -sticky => "w" );

my $lnf = $out_frame_primers->new_ttk__label( -text => "DNA concentration (nM)" );
$lnf->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $ltf = $out_frame_primers->new_ttk__label( -text => "Primer text file" );
$ltf->g_grid( -column => 6, -row => 3, -sticky => "w" );

my $lqf = $out_frame_primers->new_ttk__label( -text => "Self-complementarity" );
$lqf->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $alqf = $out_frame_primers->new_ttk__label( -text => "Graphic format (.png or .gif)" );
$alqf->g_grid( -column => 6, -row => 5, -sticky => "w" );

my $galqf = $out_frame_primers->new_ttk__label( -text => "Minimum size sequence" );
$galqf->g_grid( -column => 6, -row => 6, -sticky => "w" );

my $htmlgalqf = $out_frame_primers->new_ttk__label( -text => "Repetitive seq (Y/N)" );
$htmlgalqf->g_grid( -column => 6, -row => 7, -sticky => "w" );

foreach ( Tkx::SplitList( $out_frame_primers->g_winfo_children ) ) {
    Tkx::grid_configure( $_, -padx => 5, -pady => 5 );
}

sub exit_program {
    Tkx::tk___messageBox( -message => "Exiting PrimerMapper" );
    $main_window->g_destroy;
}
####################################
####Interface Ends####
####################################
####################################

####Collect Seqs remotely
=head1 get_sequences
 Title   :  get_sequences
 Usage   :  -command => sub { get_sequences(); }
 Function:  uses NIH e-utilities to retrieve remotely sequences in from GenBank
 Returns :  fastA formatted sequences 
=cut

sub get_sequences {
    $seq_in =~ s/\s//g;

    my $url =
      'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
      . $seq_in
      . '&rettype=fasta&retmode=text';

    my $fileList;
    $fileList = get($url);

    my $new_seqs = 'sequences.fasta';
    open SEQS, ">$new_seqs" or die;

    print SEQS "$fileList";

    close SEQS;
    print "\nSequences collected.\n";

}

####Collect SNP Seqs remotely
=head1 get_sequences_SNP
 Title   :  get_sequences_SNP
 Usage   :  -command => sub { get_sequences_SNP(); }
 Function:  uses NIH e-utilities to retrieve remotely SNP sequences in dbSNP
 Returns :  SNP formatted sequences 
=cut


sub get_sequences_SNP {
    $seq_in_SNP =~ s/\s//g;

    my $url =
        'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id='
      . $seq_in_SNP
      . '&rettype=fasta&retmode=text';

    my $fileList;
    $fileList = get($url);

    my $new_seqs = 'SNP_sequences.fasta';
    open SEQS, ">$new_seqs" or die;

    print SEQS "$fileList";

    close SEQS;
    print "\nSequences collected.\n";
}

####Load Seqs into PrimerMapper
=head1 loadfile
 Title   :  loadfile
 Usage   :  $fasta = Tkx::tk___getOpenFile();
 Function:  loads sequences into PrimerMapper program 
=cut


sub loadfile {
    $fasta = Tkx::tk___getOpenFile();

    $outfile = 'specificty.txt';
    open OUT, ">>$outfile";

    open IN, "<$fasta";
    while (<IN>) {
        print OUT if ( $_ !~ /(>.+)/ );
    }
    close IN;
    close OUT;

    #stringify the specificity file
    open OUT, $outfile or die "Couldn't open file: $!";
    $specificty = do { local $/; <OUT> };
    while (<OUT>) {
        $specificty .= $_;
    }
    close OUT;
    $specificty =~ s/\s//g;
    $specificty =~ s/[^agctn]//ig;

    print "\nSequences Loaded.\n";

}

####Load SNP Seqs into PrimerMapper
=head1 loadfile_SNP
 Title   :  loadfile_SNP
 Usage   :  $fasta_SNP = Tkx::tk___getOpenFile();
 Function:  loads SNP sequences into PrimerMapper program 
=cut


sub loadfile_SNP {
    $fasta_SNP = Tkx::tk___getOpenFile();

    $outfile_SNP = 'specificty_SNP.txt';
    open OUT_SNP, ">>$outfile_SNP ";

    open IN_SNP, "<$fasta_SNP";
    while (<IN_SNP>) {
        print OUT_SNP if ( $_ !~ />.+/ );
    }
    close IN_SNP;
    close OUT_SNP;

    #stringify the specificity file
    open OUT_SNP, $outfile_SNP or die "Couldn't open file: $!";
    $specificty_SNP = do { local $/; <OUT_SNP> };
    while (<OUT_SNP>) {
        $specificty_SNP .= $_;
    }
    close OUT_SNP;
    $specificty_SNP =~ s/\s//g;
    $specificty_SNP =~ s/[^agctn]//ig;

    print "\nSequences Loaded.\n";

}

####Parse fasta sequence using BioSeqIO
sub calculate {
    my $seqio = Bio::SeqIO->new(
        -file   => $fasta,
        -format => "fasta",
    );

    print "\nDesigning and validating primers....\n";

    while ( my $seqobj = $seqio->next_seq() ) {
        my $sequence = $seqobj->seq();
        $sequence =~ s/[^agctn]//ig;
        $sequence =~ s/\s//g;
        $sequence =~ tr/a-z/A-Z/;
        $id = $seqobj->id();
        $id =~ tr/a-zA-Z0-9//cd;
        my $len_seq  = length $sequence;
        my $id_uniq  = $len_seq . "'" . $id;
        my $len_uniq = $id . "'" . $len_seq;
        push @increment, $len_seq;
        $foo = eval join '+', @increment;

        #declare output file
        $out_image       = "GRAPHIC_$id.txt";
        $outputfile_html = "canvas_$id.txt";

    #skip to next if selected 5' or 3' length is longer than the sequence length
        if ( $five_prime_end > length $sequence ) {
            next;
        }
        if ( $three_prime_end > length $sequence ) {
            next;
        }

        #skip to next if sequence length is <100bps
        if ( length $sequence < $min_size_seq ) {
            next;
        }

        open FUSIONFILE, ">>$outputfile";
        print FUSIONFILE "Header\tStart\tTm(degC)\tSequence\tSelf-comp\tGC%\n";
        close FUSIONFILE;

        open CANVASFILE, ">>$outputfile_html";
        print CANVASFILE
          "FS\tLF\tRS\tLR\t>\t#\tSEQ\n\t\t\t\t$id\t$len_seq\t$sequence\n";
        close CANVASFILE;

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## FORWARD PRIMER
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>

######## FORWARD PRIMER
        #start counting
        my $start = 1;
        for ( my $i = $start - 1 ; $i < $five_prime_end - 1 ; $i += 1 ) {
            for ( my $a = $kmer_max ; $a >= $kmer_min ; $a-- ) {

                #$kmer = int( rand($kmer_diff) ) + $kmer_min;
                $_ = substr( $sequence, $i, $a );

                #get self complementarity score
                my $revF = reverse($_);
                $revF =~ tr/ATGCatgc/TACGtacg/;
                my $selfie_score = selfie( $_, $revF );

                #Count Gs and Cs
                my $countGC = tr/GCgc//;

                #Calculate percent GC
                my $percentGC = 100 * $countGC / $a;
                my $percentGCrounded = sprintf( "%0.1f", $percentGC );

                #calculate Tm
                if ( $a <= 36 ) {
                    $Tm = calcTm( $_, $DNA_conc, $salt );
                }
                else {
                    $Tm = calclongTm( $_, $DNA_conc, $salt, $percentGCrounded );
                }
                my $Tmrounded = sprintf( "%0.1f", $Tm );

                my $hairpin = calcdG($_);

                my $primer_end = $i + $a;

                #capture matches
                my $number_matches = 0;
                if ( $spec eq "Y" ) {
                    my $mis_mismatch = mismatch_pattern( $_, $mis );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty );
                    $number_matches = @approximate_matches;
                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched

                ############HTML

                if (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $spec eq "N" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $spec eq "N" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $number_matches < 2
                    && $spec eq "Y" )
                {
                    print FUSIONFILE
"$id\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$_";
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS);
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$i\t$foo\t$i\t$primer_end\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                }
            }
        }

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## REVERSE PRIMER
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>

######## REVERSE PRIMER

        #start counting for reverse primer
        for (
            my $j = length($sequence) - $three_prime_end ;
            $j < length($sequence) ;
            $j += 1
          )
        {
            for ( my $y = $kmer_max_SNP ; $y >= $kmer_min_SNP ; $y-- ) {

                #$kmer = int( rand($kmer_diff) ) + $kmer_min;
                $_ = substr( $sequence, $j, $y );

                #rev comp
                my $revR = reverse($_);
                $revR =~ tr/ATGCatgc/TACGtacg/;

                #get self complementarity score
                my $selfie_scoreR = selfie( $_, $revR );

                #Count Gs and Cs
                my $count_GC = tr/GCgc//;

                #Calculate percent GC
                my $percent_GC = 100 * $count_GC / $y;
                my $percentGC_rounded = sprintf( "%0.1f", $percent_GC );

                #calculate Tm
                if ( $y <= 36 ) {
                    $Tm2 = calcTm( $_, $DNA_conc, $salt );
                }
                else {
                    $Tm2 =
                      calclongTm( $_, $DNA_conc, $salt, $percentGC_rounded );
                }
                my $Tm_rounded = sprintf( "%0.1f", $Tm2 );

                my $hairpin_r = calcdG($_);

                my $primer_start_R = $j + $y;

                #capture matches
                my $number_matches_R = 0;
                if ( $spec eq "Y" ) {
                    my $mis_mismatch = mismatch_pattern( $_, $mis );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty );
                    $number_matches_R = @approximate_matches;
                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched

                ############html
                if (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $spec eq "N" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $spec eq "N" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }

                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $number_matches_R < 2
                    && $spec eq "Y" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" )
                {
                    open( OLIGOS, ">>$out_image" ) or die;
                    print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS);
                    push @array_length, $len_uniq;
                    push @array_name,   $id_uniq;
                    push @primers,      "$id:$revR";
                    open( OUTS, ">>$out_single" ) or die;
                    print OUTS "$j\t$foo\t$primer_start_R\t$j\n";
                    close(OUTS);
                    open( CANVASFILE, ">>$outputfile_html" )
                      or die;
                    print CANVASFILE "\t\t$j,\t$y,\n";
                    close(CANVASFILE);

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
            }
        }
    }

    foreach my $fp ( glob("canvas_*.txt") ) {
        open my $fh, "<", $fp or die;
        my $commander = "perl cols_to_rows.pl $fp";
        system($commander);
        close $fh;
    }
    print "\nPrimer Design completed.\n";
}

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## START SNP PRIMER DESIGN
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
sub calculate_SNP {

    $fasta_SNP_NEW = "processed_SNP.fasta";
    open INNER, '<', $fasta_SNP     or die "Can't read old file: $!";
    open OUTER, '>', $fasta_SNP_NEW or die "Can't write new file: $!";

    while (<INNER>) {
        s/ //g;
        print OUTER $_;
    }

    close OUTER;
    close INNER;

    unlink $fasta_SNP;

    my $seqio = Bio::SeqIO->new(
        -file   => $fasta_SNP_NEW,
        -format => "fasta",
    );

    print
"\nPreprocessing SNP file....\n....Designing and validating SNP primers....\n";

    while ( my $seqobj = $seqio->next_seq() ) {
        my $sequence = $seqobj->seq();
        $id_SNP = $seqobj->id();

        if ( $id_SNP =~ m/pos=(.+).len/i ) {
            $snp = $1;
        }

        ###################
        if ( $sequence =~ m/R/ ) {
            $degen = "type_1";
        }
        elsif ( $sequence =~ m/Y/ ) {
            $degen = "type_2";
        }
        elsif ( $sequence =~ m/S/ ) {
            $degen = "type_3";
        }
        elsif ( $sequence =~ m/W/ ) {
            $degen = "type_4";
        }
        elsif ( $sequence =~ m/K/ ) {
            $degen = "type_5";
        }
        elsif ( $sequence =~ m/M/ ) {
            $degen = "type_6";
        }
        elsif ( $sequence =~ m/B/ ) {
            $degen = "type_7";
        }
        elsif ( $sequence =~ m/D/ ) {
            $degen = "type_8";
        }
        elsif ( $sequence =~ m/H/ ) {
            $degen = "type_9";
        }
        elsif ( $sequence =~ m/V/ ) {
            $degen = "type_10";
        }

        # Degenerate SNP code:
##-->R	A or G
##-->Y	C or T
##-->S	G or C
##-->W	A or T
##-->K	G or T
##-->M	A or C
##-->B	C or G or T
##-->D	A or G or T
##-->H	A or C or T
##-->V	A or C or G

        ###################

        $sequence =~ s/[^agctn]//ig;
        $sequence =~ s/\s//g;
        $sequence =~ tr/a-z/A-Z/;

        $id_SNP =~ tr/a-zA-Z0-9//cd;
        my $len_seq      = length $sequence;
        my $id_uniq_SNP  = $len_seq . "'" . $id_SNP;
        my $len_uniq_SNP = $id_SNP . "'" . $len_seq;
        push @increment_SNP, $len_seq;
        $foo_SNP = eval join '+', @increment_SNP;

        #declare output file
        $out_image_SNP = "GRAPHIC_$id_SNP.txt";

    #skip to next if selected 5' or 3' length is longer than the sequence length
        if ( $five_prime_SNP > $snp ) {
            next;
        }
        if ( $three_prime_SNP > ( $len_seq - $snp ) ) {
            next;
        }

        #skip to next if sequence length is <100bps
        if ( length $sequence < $min_size_seq ) {
            next;
        }

        open FUSIONFILE_SNP, ">>$outputfile_SNP ";
        print FUSIONFILE_SNP
          "Header\tStart\tTm(degC)\tSequence\tSelf-comp\tGC%\n";
        close FUSIONFILE_SNP;

        open AS_SNP, ">>$outputfile_SNP_AS ";
        print AS_SNP
"Header\tTm(degC)\tAllele Specific Primer Sequences\t\tSelf-comp\tHairpin\tGC%\n";
        close AS_SNP;

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## FORWARD SNP PRIMER
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>

        #start counting
        my $upstream = $snp - $five_prime_SNP;

        for ( my $i = $upstream ; $i < $snp - $snp_distance ; $i += 1 ) {
            for ( my $q = $kmer_max_SNP ; $q >= $kmer_min_SNP ; $q-- ) {

                #$kmer_SNP = int( rand($kmer_diff_SNP) ) + $kmer_min_SNP;
                $_ = substr( $sequence, $i, $q );

                #get self complementarity score
                my $revF = reverse($_);
                $revF =~ tr/ATGCatgc/TACGtacg/;

                my $selfie_score = selfie( $_, $revF );

                #Count Gs and Cs
                my $countGC = tr/GCgc//;

                #Calculate percent GC
                my $percentGC = 100 * $countGC / $q;
                my $percentGCrounded = sprintf( "%0.1f", $percentGC );

                #calculate Tm
                if ( $q <= 36 ) {
                    $Tm_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
                }
                else {
                    $Tm_SNP = calclongTm( $_, $DNA_conc_SNP,
                        $salt_SNP, $percentGCrounded );
                }
                my $Tmrounded = sprintf( "%0.1f", $Tm_SNP );

                my $hairpin = calcdG($_);

                my $primer_end = $i + $q;

                #capture matches
                my $number_matches = 0;
                if ( $spec_SNP eq "Y" ) {
                    my $mis_mismatch = mismatch_pattern( $_, $mis_SNP );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty_SNP );
                    $number_matches = @approximate_matches;
                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                if (   open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image" )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $clamp_SNP eq "N"
                    && $spec_SNP eq "N" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                if (   open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp_SNP eq "Y"
                    && $number_matches < 2
                    && $spec_SNP eq "Y" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp_SNP eq "Y"
                    && $number_matches < 2
                    && $spec_SNP eq "Y" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && calcRepeat( $_, $repetition ) == 1
                    && $hairpin > "-9"
                    && $clamp_SNP eq "N"
                    && $number_matches < 2
                    && $spec_SNP eq "Y" )
                {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP " )
                      or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                }
            }
        }

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## REVERSE SNP PRIMER
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>

        #start counting for reverse primer
        my $downstream = $snp + $three_prime_SNP;

        for ( my $j = $snp + $snp_distance ; $j < $downstream ; $j += 1 ) {
            for ( my $c = $kmer_max_SNP ; $c >= $kmer_min_SNP ; $c-- ) {

                #$kmer_SNP = int( rand($kmer_diff_SNP) ) + $kmer_min_SNP;
                $_ = substr( $sequence, $j, $c );

                #rev comp
                my $revR = reverse($_);
                $revR =~ tr/ATGCatgc/TACGtacg/;

                #get self complementarity score
                my $selfie_scoreR = selfie( $_, $revR );

                #Count Gs and Cs
                my $count_GC = tr/GCgc//;

                #Calculate percent GC
                my $percent_GC = 100 * $count_GC / $c;
                my $percentGC_rounded = sprintf( "%0.1f", $percent_GC );

                #calculate Tm
                if ( $c <= 36 ) {
                    $Tm2_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
                }
                else {
                    $Tm2_SNP =
                      calclongTm( $_, $DNA_conc_SNP,
                        $salt_SNP, $percentGC_rounded );
                }
                my $Tm_rounded = sprintf( "%0.1f", $Tm2_SNP );

                my $hairpin_r = calcdG($_);

                my $primer_start_R = $j + $c;

                #capture matches
                my $number_matches_R = 0;
                if ( $spec_SNP eq "Y" ) {
                    my $mis_mismatch = mismatch_pattern( $_, $mis_SNP );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty_SNP );
                    $number_matches_R = @approximate_matches;
                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                if (   open( FUSIONFILE_SNP, ">>$outputfile_SNP" )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $clamp_SNP eq "N"
                    && $spec_SNP eq "N" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP" )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $clamp_SNP eq "N"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP" )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp_SNP eq "Y"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tm_rounded ge $lower_tm_SNP
                    && $Tm_rounded le $upper_tm_SNP
                    && $percent_GC ge $lower_gc_SNP
                    && $percent_GC le $higher_gc_SNP
                    && $selfie_scoreR < $selfie_cuttof_SNP
                    && calcRepeat( $revR, $repetition ) == 1
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp_SNP eq "Y"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                {
                    open( OLIGOS_SNP, ">>$out_image_SNP" )
                      or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" )
                      or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                }
            }
        }
########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## ALLELE SPECIFIC PRIMER
########################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>
################################################################################################################################################################>>>>>>>>>>>>>>>>>>>>>>>>

        ####################______Allele-Specific_BEGIN

        for ( my $p = $kmer_min_SNP ; $p <= $kmer_max_SNP ; $p++ ) {
            $_ = substr( $sequence, $snp - 1, $p );

            #get self complementarity score
            my $revFA = reverse($_);
            $revFA =~ tr/ATGCatgc/TACGtacg/;

            my $selfie_scoreA = selfie( $_, $revFA );

            #Count Gs and Cs
            my $countGCA = tr/GCgc//;

            #Calculate percent GC
            my $percentGCA = 100 * $countGCA / $p;
            my $percentGCroundedA = sprintf( "%0.1f", $percentGCA );

            #calculate Tm
            if ( $p <= 36 ) {
                $Tm3_SNP = calcTm( $revFA, $DNA_conc_SNP, $salt_SNP );
            }
            else {
                $Tm3_SNP =
                  calclongTm( $revFA, $DNA_conc_SNP,
                    $salt_SNP, $percentGCroundedA );
            }
            my $TmroundedA = sprintf( "%0.1f", $Tm3_SNP );

            my $hairpinA = calcdG($revFA);

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.G\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.C\t$revFA.G\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.G\t$revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_10" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.C\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

            }

        }

###################

        my $string_to_snp = $snp - 2;
        my $new_seq = substr( $sequence, 1, $string_to_snp );

        for ( my $u = $kmer_max_SNP ; $u >= $kmer_min_SNP ; $u-- ) {
            $_ = substr( $new_seq, -$u );

            #get self complementarity score
            my $revFAB = reverse($_);
            $revFAB =~ tr/ATGCatgc/TACGtacg/;

            my $selfie_scoreAB = selfie( $_, $revFAB );

            #Count Gs and Cs
            my $countGCAB = tr/GCgc//;

            #Calculate percent GC
            my $percentGCAB = 100 * $countGCAB / $u;
            my $percentGCroundedAB = sprintf( "%0.1f", $percentGCAB );

            #calculate Tm
            if ( $u <= 36 ) {
                $Tm4_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
            }
            else {
                $Tm4_SNP =
                  calclongTm( $_, $DNA_conc_SNP, $salt_SNP,
                    $percentGCroundedAB );
            }
            my $TmroundedAB = sprintf( "%0.1f", $Tm4_SNP );

            my $hairpinAB = calcdG($_);

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.G\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.C\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.G\t$_.C\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.C\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.C\t$_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tReverse:$_.A\t$_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.C\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_10" )
            {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.C\t$_.G\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

            }

        }

        ################

    }

    ################### ----> Allele-Specific_END
    print "\nSNP Primer Design completed.\n";
}

Tkx::MainLoop();

####################################
####################################
####################################
####################################
####################################
=head1 primer_dimer
 Title   :  primer_dimer
 Usage   :  -command => sub { primer_dimer(\@primers); }
 Function:  implements a combinations without replacements algorithm (n choose k) for all primers
            (both forward and reverse) to calculate cross-complementarity primer-dimer scores
 Returns :  primer_dimer_scores.tsv
=cut

sub primer_dimer {
    print
      "\nCalculating global primer dimer score file for multiplex PCR....\n";

    my $primer_dimer_input = "primers_combos.txt";
    my $primer_dimer;

    open IN, ">$primer_dimer_input" or die "Couldn't open file: $!";

    my $strings = \@primers;

    sub combine;

    print IN "@$_\n" for combine $strings, 2;
    close IN;

#--GET_Primer_Dimer_Scores--
    my $primer_dimer_score = "primer_dimer_scores.tsv";
    open IN,    "<", "$primer_dimer_input" or die "Couldn't open file: $!";
    open OUTER, '>', $primer_dimer_score   or die "Can't write new file: $!";
    print OUTER
"BASED ON IN SILICO AND PCR TESTING, PRIMER DIMER SCORES SHOULD BE LESS THAN OR EQUAL TO 12 \n\n";
    while (<IN>) {

        if ( $_ =~ m/(.+):(.+)\s(.+):(.+)\n/i ) {
            $primer_dimer = selfie( $2, $4 );
            print OUTER "$1:$2\t$3:$4\tDimer_score: $primer_dimer\n";
        }
    }
    close OUTER;
    close IN;

    unlink $primer_dimer_input;

    print "\nGlobal Primer Dimer values collected.\n";
}

####################################
=head1 primer_dimer_SNP
 Title   :  primer_dimer_SNP
 Usage   :  -command => sub { primer_dimer_SNP(\@SNP_primers); }
 Function:  implements a combinations without replacements algorithm (n choose k) for all SNP primers
            (both forward and reverse) to calculate cross-complementarity SNP primer-dimer scores
 Returns :  primer_dimer_scores_SNP.tsv
=cut

sub primer_dimer_SNP {

    print
"\nCalculating global SNP primer dimer score file for multiplex PCR....\n";

    my $primer_dimer_input = "primers_combos.txt";
    my $primer_dimer;

    open IN, ">$primer_dimer_input" or die "Couldn't open file: $!";

    my $strings = \@SNP_primers;

    sub combine;

    print IN "@$_\n" for combine $strings, 2;
    close IN;

#--GET_Primer_Dimer_Scores--
    my $primer_dimer_score = "primer_dimer_scores_SNP.tsv";
    open IN,    "<", "$primer_dimer_input" or die "Couldn't open file: $!";
    open OUTER, '>', $primer_dimer_score   or die "Can't write new file: $!";
    print OUTER
"BASED ON IN SILICO AND PCR TESTING, PRIMER DIMER SCORES SHOULD BE LESS THAN OR EQUAL TO 12 \n\n";
    while (<IN>) {

        if ( $_ =~ m/(.+):(.+)\s(.+):(.+)\n/i ) {
            push @SNP_primers, "F:$_";
            $primer_dimer = selfie( $2, $4 );
            print OUTER "$1:$2\t$3:$4\tDimer_score: $primer_dimer\n";
        }
    }
    close OUTER;
    close IN;

    unlink $primer_dimer_input;

    print "\nGlobal Primer Dimer values collected.\n";
}


sub combine {

    my ( $list, $n ) = @_;
    die "Insufficient list members for primer dimer analysis"
      if $n > @$list;

    return map [$_], @$list if $n <= 1;

    my @comb;

    for ( my $i = 0 ; $i + $n <= @$list ; ++$i ) {
        my $val  = $list->[$i];
        my @rest = @$list[ $i + 1 .. $#$list ];
        push @comb, [ $val, @$_ ] for combine \@rest, $n - 1;
    }

    return @comb;
}

####################################
####################################

=pod 
=WEB-BLAST
Note: must have web_blast.pl script in PATH for BLAST feature 
=cut

####################################
####################################
####BLAST SUBROUTINE####
####################################
####################################
####################################
####################################
####################################
=head1 web_blast_ncbi
 Title   :  web_blast_ncbi
 Usage   :  -command => sub { web_blast_ncbi(); }
 Function:  runs remote BLAST from PrimerMapper using the script web_blast.pl
 Returns :  BLAST result in single text file 
=cut


sub web_blast_ncbi {

    print "\nBLAST Analysis of $blaster is running....\n";

    $thetext = $text->get( "1.0", "end" );
    my $blastfile = "blast_primers.fasta";
    open( FILE, ">>$blastfile" ) or print "unable to open file";
    print FILE $thetext;

    my $command = "perl web_blast.pl blastn $blaster blast_primers.fasta";
    system($command);

    close FILE;

    print "\nBLAST Analysis completed.\n";
}

####################################
####################################
1;

