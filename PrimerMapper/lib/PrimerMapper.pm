#
# Main module for PrimerMapper
#
#       uses the software primer3
#       https://sourceforge.net/projects/primer3/
#       by Steve Rozen et al.
#       and the Bio::Graphics modules within Bioperl
#       by Lincoln Stein
#       and by Ewan Birney
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
    my $min_prod_size  = "150";
    my $max_prod_size = "150";
    my $kmer_max        = "25";
    my $kmer_min        = "18";
    my $clamp           = "0";
    my $higher_gc       = "60";
    my $lower_gc        = "40";
    my $upper_tm        = "68";
    my $lower_tm        = "55";
    my $salt            = "50";
    my $outputfile      = "primers.tsv";
    my $format          = ".png";
SNP Primer Design Defaults
    my $five_prime_SNP    = "120";
    my $three_prime_SNP   = "120";
    my $kmer_max_SNP      = "25";
    my $kmer_min_SNP      = "18";
    my $clamp_SNP         = "0";
    my $higher_gc_SNP     = "60";
    my $lower_gc_SNP      = "40";
    my $upper_tm_SNP      = "68";
    my $lower_tm_SNP      = "55";
    my $salt_SNP          = "50";
    my $snp_distance      = "30";

the defaults will be overwritten if a parameter is changed

=head3 DESCRIPTION
This object produces the following:

1. a GUI to facilitate primer design and visualization
2. Traditional primer design as well as SNP and allele specific primer design
3. Visualization of primer distribution across each sequence and entire input file
4. Remote sequence access from GenBank and dbSNP
5. Primer BLAST facility against multiple NCBI databases 

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

our $VERSION = '3.0';

#########################
#########################
sub new {
    my $class = shift;
    return bless {}, $class;
}
#########################

####Global Primer Design Defaults
my $fasta;
my $min_prod_size = "100";
my $max_prod_size = "800";
my $kmer_max      = "25";
my $kmer_min      = "18";
my $clamp         = "0";
my $higher_gc     = "60";
my $lower_gc      = "40";
my $upper_tm      = "68";
my $lower_tm      = "55";
my $salt          = "50";
my $format        = ".png";
my $outputfile    = "primers.tsv";
###########################

####Global Strings, Arrays, and Hashesou
my $outfile;
my $out_single = "single_viewer.txt";
my $seq_in     = "NM_001307521.1,FJ455617.1";
my @array_length;
my @array_name;
my @gene_length;
my @increment;
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
my $five_prime_SNP  = "180";
my $three_prime_SNP = "180";
my $kmer_max_SNP    = "25";
my $kmer_min_SNP    = "18";
my $clamp_SNP       = "0";
my $higher_gc_SNP   = "60";
my $lower_gc_SNP    = "40";
my $upper_tm_SNP    = "68";
my $lower_tm_SNP    = "55";
my $salt_SNP        = "50";
my $outputfile_SNP  = "SNP_primers.tsv";
my $snp_distance    = "30";
###########################

####Global SNP Strings and Hashes
my $kmer_diff_SNP = $kmer_max_SNP - $kmer_min_SNP;
my $kmer_SNP;
my $snp;
my $outfile_SNP;
my $out_single_SNP = "single_viewer_SNP.txt";
my $seq_in_SNP     = "rs744373,rs11136000,rs3764650";
my @array_length_SNP;
my @array_name_SNP;
my @gene_length_SNP;
my @increment_SNP;
my @SNP_position;
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

my $input_ =
  $input_frame->new_ttk__entry( -width => 40, -textvariable => \$fasta );
$input_->g_grid( -column => 1, -row => 0, -sticky => "we" );
my $load_file_ =
  $input_frame->new_ttk__button( -text => "Load File", -command => \&loadfile );
$load_file_->g_grid( -column => 0, -row => 0, -sticky => "w" );
my $load_seq_ =
  $input_frame->new_ttk__entry( -width => 40, -textvariable => \$seq_in );
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

my $input_snp_ = $outer_frame_snp->new_ttk__entry(
    -width        => 40,
    -textvariable => \$fasta_SNP
);
$input_snp_->g_grid( -column => 1, -row => 14, -sticky => "w" );
my $load_snp_ = $outer_frame_snp->new_ttk__button(
    -text    => "Load SNP File",
    -command => \&loadfile_SNP
);
$load_snp_->g_grid( -column => 0, -row => 14, -sticky => "w" );
my $load_seq_snp_ = $outer_frame_snp->new_ttk__entry(
    -width        => 40,
    -textvariable => \$seq_in_SNP
);
$load_seq_snp_->g_grid( -column => 4, -row => 14, -sticky => "w" );
my $get_snp_seqs_ = $outer_frame_snp->new_ttk__button(
    -text    => "Get SNP Sequences",
    -command => sub { get_sequences_SNP(); }
);
$get_snp_seqs_->g_grid( -column => 2, -row => 14, -sticky => "w" );

my $design_frame_ =
  $outer_frame->new_ttk__labelframe( -text => "PRIMER DESIGN" );
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

my $kmer_max_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$kmer_max_SNP
);
$kmer_max_snp_->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $kmer_min_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$kmer_min_SNP
);
$kmer_min_snp_->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $upper_gc_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$higher_gc_SNP
);
$upper_gc_snp_->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $lower_gc_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$lower_gc_SNP
);
$lower_gc_snp_->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $upper_tm_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$upper_tm_SNP
);
$upper_tm_snp_->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $lower_tm_snp_ = $design_frame_->new_ttk__entry(
    -width        => 3,
    -textvariable => \$lower_tm_SNP
);
$lower_tm_snp_->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $gc_clamp_snp_ =
  $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$clamp_SNP );
$gc_clamp_snp_->g_grid( -column => 7, -row => 5, -sticky => "e" );

my $salt_snp_ =
  $design_frame_->new_ttk__entry( -width => 3, -textvariable => \$salt_SNP );
$salt_snp_->g_grid( -column => 7, -row => 2, -sticky => "we" );

my $format_snp_ =
  $design_frame_->new_ttk__entry( -width => 10, -textvariable => \$format );
$format_snp_->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $snp_distance_ = $design_frame_->new_ttk__entry(
    -width        => 10,
    -textvariable => \$snp_distance
);
$snp_distance_->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $distance_from_snp_ =
  $design_frame_->new_ttk__label( -text => "5' up from SNP" );
$distance_from_snp_->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $distance_down_from_snp_ =
  $design_frame_->new_ttk__label( -text => "3' down from SNP" );
$distance_down_from_snp_->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $max_primer_len_snp_ =
  $design_frame_->new_ttk__label( -text => "Primer length max" );
$max_primer_len_snp_->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $min_primer_len_snp_ =
  $design_frame_->new_ttk__label( -text => "Primer length min" );
$min_primer_len_snp_->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $max_gc_snp_ = $design_frame_->new_ttk__label( -text => "Upper GC%" );
$max_gc_snp_->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $min_gc_snp_ = $design_frame_->new_ttk__label( -text => "Lower GC%" );
$min_gc_snp_->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $max_tm_snp_ = $design_frame_->new_ttk__label( -text => "Upper Tm" );
$max_tm_snp_->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $min_tm_snp_ = $design_frame_->new_ttk__label( -text => "Lower Tm" );
$min_tm_snp_->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $clamp_snp_ = $design_frame_->new_ttk__label( -text => "GC clamp (1/0)" );
$clamp_snp_->g_grid( -column => 6, -row => 5, -sticky => "w" );

my $nacl_concentration_snp_ =
  $design_frame_->new_ttk__label( -text => "salt concentration (mM)" );
$nacl_concentration_snp_->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $min_basepairs_from_snp_ =
  $design_frame_->new_ttk__label( -text => "Min. distance from SNP" );
$min_basepairs_from_snp_->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $outfile_format_snp_ =
  $design_frame_->new_ttk__label( -text => "png or gif" );
$outfile_format_snp_->g_grid( -column => 6, -row => 3, -sticky => "w" );

foreach ( Tkx::SplitList( $design_frame_->g_winfo_children ) ) {
    Tkx::grid_configure( $_, -padx => 5, -pady => 5 );
}

my $run_snp_frame_ = $outer_frame->new_ttk__labelframe(
    -text   => "RUN",
    -width  => 10,
    -height => 60
);
$run_snp_frame_->new_ttk__frame( -padding => "3 3 12 12" );
$run_snp_frame_->g_grid( -column => 0, -row => 28, -sticky => "nwes" );

my $design_snp_frame_ = $run_snp_frame_->new_ttk__button(
    -text    => "Design SNP Primers and Graphics",
    -command => sub {
        calculate_SNP();
        graphics_all_primers_SNP( \@array_length_SNP, \@array_name_SNP,
            \@gene_length_SNP, $format, $snp );
        graphics_single_view_SNP( \@gene_length_SNP, $format, $out_single_SNP );
    }
);
$design_snp_frame_->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $clean_up_snp_ = $run_snp_frame_->new_ttk__button(
    -text    => "Clean-up",
    -command => sub { clean_up(); }
);
$clean_up_snp_->g_grid( -column => 5, -row => 9, -sticky => "s" );

########## SNP BTM

my $out_frame_primers =
  $outer_frame->new_ttk__labelframe( -text => "PRIMER DESIGN" );
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
    -textvariable => \$min_prod_size
);
$five_primer_ending_->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $three_primer_ending_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$max_prod_size
);
$three_primer_ending_->g_grid( -column => 5, -row => 2, -sticky => "we" );

my $max_primer_size_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$kmer_max
);
$max_primer_size_allowed_->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $min_primer_size_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$kmer_min
);
$min_primer_size_allowed_->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $max_gc_content_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$higher_gc
);
$max_gc_content_allowed_->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $min_gc_content_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$lower_gc
);
$min_gc_content_allowed_->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $max_primer_tm_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$upper_tm
);
$max_primer_tm_allowed_->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $min_primer_tm_allowed_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$lower_tm
);
$min_primer_tm_allowed_->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $add_gc_clamp_ =
  $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$clamp );
$add_gc_clamp_->g_grid( -column => 7, -row => 2, -sticky => "e" );

my $sodium_chloride_ =
  $out_frame_primers->new_ttk__entry( -width => 3, -textvariable => \$salt );
$sodium_chloride_->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $outfile__ = $out_frame_primers->new_ttk__entry(
    -width        => 10,
    -textvariable => \$outputfile
);
$outfile__->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $format__ =
  $out_frame_primers->new_ttk__entry( -width => 10, -textvariable => \$format );
$format__->g_grid( -column => 7, -row => 5, -sticky => "we" );

my $run_frame__ = $outer_frame->new_ttk__labelframe(
    -text   => "RUN",
    -width  => 10,
    -height => 40
);
$run_frame__->new_ttk__frame( -padding => "3 3 12 12" );
$run_frame__->g_grid( -column => 0, -row => 8, -sticky => "nwes" );

my $blast_frame_ = $main_window->new_ttk__labelframe(
    -text => "BLAST PRIMERS - insert primers here in fastA format" );
$blast_frame_->new_ttk__frame( -padding => "3 3 12 12" );
$blast_frame_->g_grid( -column => 0, -row => 10, -sticky => "nwes" );

my $text =
  $blast_frame_->new_tk__text( -width => 80, -height => 3, -wrap => "none" );
$text->g_grid;
my $thetext;

my $exit_ = $main_window->new_ttk__button(
    -text    => "EXIT",
    -command => sub { exit_program(); }
);
$exit_->g_grid( -column => 0, -row => 14, -sticky => "s" );

my $which_db_ =
  $main_window->new_ttk__labelframe( -text => "SELECT DATABASE:" );
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
    -text    => "Design Primers and Graphics",
    -command => sub {
        calculate();
        graphics_all_primers( \@array_length, \@array_name, \@gene_length,
            $format );
        graphics_single_view( \@gene_length, $format, $out_single );
    }
);
$get_the_primers_->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $cleaner__ = $run_frame__->new_ttk__button(
    -text    => "Clean-up",
    -command => sub { clean_up(); }
);
$cleaner__->g_grid( -column => 5, -row => 9, -sticky => "s" );

my $search_area_five__prime_ =
  $out_frame_primers->new_ttk__label( -text => "Max. product size" );
$search_area_five__prime_->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $search_area_three__prime_ =
  $out_frame_primers->new_ttk__label( -text => "Min. product size" );
$search_area_three__prime_->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $max_permitted_oligo_len__ =
  $out_frame_primers->new_ttk__label( -text => "Primer length max" );
$max_permitted_oligo_len__->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $min_permitted_oligo_len__ =
  $out_frame_primers->new_ttk__label( -text => "Primer length min" );
$min_permitted_oligo_len__->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $max_permitted_oligo_gc__ =
  $out_frame_primers->new_ttk__label( -text => "Upper GC%" );
$max_permitted_oligo_gc__->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $min_permitted_oligo_gc__ =
  $out_frame_primers->new_ttk__label( -text => "Lower GC%" );
$min_permitted_oligo_gc__->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $max_permitted_oligo_tm__ =
  $out_frame_primers->new_ttk__label( -text => "Upper Tm" );
$max_permitted_oligo_tm__->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $min_permitted_oligo_tm__ =
  $out_frame_primers->new_ttk__label( -text => "Lower Tm" );
$min_permitted_oligo_tm__->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $clamp_permitted_oligo__ =
  $out_frame_primers->new_ttk__label( -text => "GC clamp (1/0)" );
$clamp_permitted_oligo__->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $nacl_permitted_oligo__ =
  $out_frame_primers->new_ttk__label( -text => "salt conc. (mM)" );
$nacl_permitted_oligo__->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $regular_oligo_text_outfile__ =
  $out_frame_primers->new_ttk__label( -text => "Primer text file" );
$regular_oligo_text_outfile__->g_grid(
    -column => 6,
    -row    => 3,
    -sticky => "w"
);

my $regular_oligo_graphic_format_outfile__ =
  $out_frame_primers->new_ttk__label( -text => "png or gif" );
$regular_oligo_graphic_format_outfile__->g_grid(
    -column => 6,
    -row    => 5,
    -sticky => "w"
);

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
        $out_image = "GRAPHIC_$id.txt";

        my $p3 = $id."primer3_output";
        open my $fh, ">>", $p3 or die;
        close $fh;

        my $new_sequence =
            "SEQUENCE_TEMPLATE="
          . $sequence . "\n"
          . "PRIMER_TASK=generic" . "\n"
          . "PRIMER_PICK_LEFT_PRIMER=1" . "\n"
          . "PRIMER_PICK_RIGHT_PRIMER=1" . "\n"
          . "PRIMER_PRODUCT_SIZE_RANGE="
          . $min_prod_size . "-"
          . $max_prod_size . "\n"
          . "PRIMER_SALT_MONOVALENT="
          . $salt . "\n"
          . "PRIMER_GC_CLAMP="
          . $clamp . "\n"
          . "PRIMER_MAX_TM="
          . $upper_tm . "\n"
          . "PRIMER_MIN_TM="
          . $lower_tm . "\n"
          . "PRIMER_MIN_SIZE="
          . $kmer_min . "\n"
          . "PRIMER_MAX_SIZE="
          . $kmer_max . "\n"
          . "PRIMER_MAX_GC="
          . $higher_gc . "\n"
          . "PRIMER_MIN_GC="
          . $lower_gc . "\n"
          . "PRIMER_MAX_NS_ACCEPTED=1" . "\n"
          . "PRIMER_EXPLAIN_FLAG=1" . "\n" . "=";

        my $temp_seq = "primer3_temp.txt";

        open my $p3_fh, ">", $temp_seq or die;
        print $p3_fh($new_sequence);

        if ( $^O eq 'MSWin32' ) {
            my $command = "primer3_core.exe -format_output < $temp_seq >> $p3";
            system($command);
        }
        else {
            my $command = "primer3_core -format_output < $temp_seq >> $p3";
            system($command);
        }

        open $fh, "<", $p3 or die;
        my $count = 0;
        while ( my $row = <$fh> ) {

            #print $row;
            if ( $row =~
/LEFT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                my $add = $1 + $2;
                open my $fh_oligo, ">>", $out_image or die;
                print $fh_oligo(
                    $1 . "\t" . $2 . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_oligo;
                open my $fh_outs, ">>", $out_single or die;
                print $fh_outs(
                    $1 . "\t" . $foo . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_outs;
            }
            elsif ( $row =~
/RIGHT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                my $minus = $1 - $2;
                open my $fh_oligo, ">>", $out_image or die;
                print $fh_oligo(
                    $1 . "\t" . $2 . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_oligo;
                open my $fh_outs, ">>", $out_single or die;
                print $fh_outs(
                    $1 . "\t" . $foo . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_outs;
            }
        }
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

        open AS_SNP, ">>$outputfile_SNP_AS ";
        print AS_SNP
          "Header\tTm(degC)\tAllele Specific Primer Sequences\t\t\tGC%\n";
        close AS_SNP;

########################################>>>>>>>>>>>>>>>>>>>>>>>>
######## SNP PRIMERS

        my $short_id = substr $id_SNP, 0, 20;
        my $p3_SNP = $short_id . "primer3_output_SNP";
        open my $fh_SNP, ">>", $p3_SNP or die;
        print $fh_SNP "Sequence Name: " . $short_id . "\n";
        close $fh_SNP;
        my $anchor      = $snp - 150;
        my $dist_before = $five_prime_SNP - $snp_distance;
        my $dist_after  = $snp + $snp_distance;

        my $new_sequence_SNP =
            "SEQUENCE_TEMPLATE="
          . $sequence . "\n"
          . "SEQUENCE_TARGET="
          . $snp . ",1" . "\n"
          . "PRIMER_TASK=generic" . "\n"
          . "PRIMER_PICK_LEFT_PRIMER=1" . "\n"
          . "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST="
          . $anchor . ","
          . $dist_before . ","
          . $dist_after . ","
          . $three_prime_SNP . "\n"
          . "PRIMER_PRODUCT_SIZE_RANGE=100-1500" . "\n"
          . "PRIMER_MAX_TM="
          . $upper_tm_SNP . "\n"
          . "PRIMER_MIN_TM="
          . $lower_tm_SNP . "\n"
          . "PRIMER_GC_CLAMP="
          . $clamp_SNP . "\n"
          . "PRIMER_SALT_MONOVALENT="
          . $salt_SNP . "\n"
          . "PRIMER_MAX_GC="
          . $higher_gc_SNP . "\n"
          . "PRIMER_MIN_GC="
          . $lower_gc_SNP . "\n"
          . "PRIMER_PICK_RIGHT_PRIMER=1" . "\n"
          . "PRIMER_MIN_SIZE="
          . $kmer_min_SNP . "\n"
          . "PRIMER_MAX_SIZE="
          . $kmer_max_SNP . "\n"
          . "PRIMER_MAX_NS_ACCEPTED=1" . "\n"
          . "PRIMER_EXPLAIN_FLAG=1" . "\n" . "=";

        my $temp_seq_SNP = "primer3_temp.txt";

        open my $p3_fh_SNP, ">", $temp_seq_SNP or die;
        print $p3_fh_SNP($new_sequence_SNP);

        if ( $^O eq 'MSWin32' ) {
            my $command_SNP =
              "primer3_core.exe -format_output < $temp_seq_SNP >> $p3_SNP";
            system($command_SNP);
        }
        else {
            my $command_SNP =
              "primer3_core -format_output < $temp_seq_SNP >> $p3_SNP";
            system($command_SNP);
        }

        open $fh_SNP, "<", $p3_SNP or die;
        my $count = 0;
        while ( my $row = <$fh_SNP> ) {

            #print $row;
            if ( $row =~
/LEFT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length_SNP, $len_uniq_SNP;
                push @array_name_SNP,   $id_uniq_SNP;
                my $add = $1 + $2;
                open my $fh_oligo_SNP, ">>", $out_image_SNP or die;
                print $fh_oligo_SNP(
                    $snp . "\t" . $2 . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_oligo_SNP;
                open my $fh_outs_SNP, ">>", $out_single_SNP or die;
                print $fh_outs_SNP(
                    $1 . "\t" . $foo_SNP . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_outs_SNP;
            }
            elsif ( $row =~
/RIGHT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length_SNP, $len_uniq_SNP;
                push @array_name_SNP,   $id_uniq_SNP;
                my $minus = $1 - $2;
                open my $fh_oligo_SNP, ">>", $out_image_SNP or die;
                print $fh_oligo_SNP(
                    $snp . "\t" . $2 . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_oligo_SNP;
                open my $fh_outs_SNP, ">>", $out_single_SNP or die;
                print $fh_outs_SNP(
                    $1 . "\t" . $foo_SNP . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_outs_SNP;
            }
        }

        print "\nSNP Primer Design completed.\n";

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

            #Count Gs and Cs
            my $countGCA = tr/GCgc//;

            #Calculate percent GC
            my $percentGCA = 100 * $countGCA / $p;
            my $percentGCroundedA = sprintf( "%0.1f", $percentGCA );

            #calculate Tm
            if ( $p <= 36 ) {
                $Tm3_SNP = calcTm( $revFA, 100, $salt_SNP );
            }
            else {
                $Tm3_SNP =
                  calclongTm( $revFA, 100, $salt_SNP, $percentGCroundedA );
            }
            my $TmroundedA = sprintf( "%0.1f", $Tm3_SNP );

            if ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.G\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.C\t$revFA.G\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.C\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.G\t$revFA.C\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.A\t$percentGCroundedA\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_10" )
            {
                print AS_SNP
"$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.C\t$percentGCroundedA\n";

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

            #Count Gs and Cs
            my $countGCAB = tr/GCgc//;

            #Calculate percent GC
            my $percentGCAB = 100 * $countGCAB / $u;
            my $percentGCroundedAB = sprintf( "%0.1f", $percentGCAB );

            #calculate Tm
            if ( $u <= 36 ) {
                $Tm4_SNP = calcTm( $_, 100, $salt_SNP );
            }
            else {
                $Tm4_SNP =
                  calclongTm( $_, 100, $salt_SNP, $percentGCroundedAB );
            }
            my $TmroundedAB = sprintf( "%0.1f", $Tm4_SNP );

            if ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.A\t$_.G\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.C\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.G\t$_.C\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.A\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.G\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.A\t$_.C\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.C\t$_.G\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tReverse:$_.A\t$_.G\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.A\t$_.C\t$_.T\t$percentGCroundedAB\n";

            }

            elsif ( open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_10" )
            {
                print AS_SNP
"$short_id\t$TmroundedAB \tForward: $_.A\t$_.C\t$_.G\t$percentGCroundedAB\n";

            }

        }

        ################
    }
    ################### ----> Allele-Specific_END
    print "\nSNP Primer Design completed.\n";
}

Tkx::MainLoop();

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

