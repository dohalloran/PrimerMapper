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
    my $min_prod_size  = "100";
    my $max_prod_size  = "800";
    my $kmer_max       = "25";
    my $kmer_min       = "18";
    my $clamp          = "0";
    my $higher_gc      = "60";
    my $lower_gc       = "40";
    my $upper_tm       = "68";
    my $lower_tm       = "55";
    my $salt           = "50";
    my $outputfile     = "primers.tsv";
    my $format         = ".png";
SNP Primer Design Defaults
    my $five_prime_SNP    = "180";
    my $three_prime_SNP   = "180";
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
The design of the code and its execution relies heavily upon the BioPerl toolkit:
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
use LWP::UserAgent;    # BUG FIX: added for HTTPS support (NCBI now requires it)
use Cwd;
use Bio::SeqIO;
use FindBin qw($Bin);
use lib "$Bin";
use PrimerMapperUtilities;
use PrimerMapperGraphics;
####

#########################
#########################

our $VERSION = '3.2';

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

####Global Strings, Arrays, and Hashes
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

sub user_message {
    my $message = shift || "PrimerMapper could not complete this action.";
    Tkx::tk___messageBox( -message => $message );
}

sub run_gui_action {
    my $action = shift;
    my $ok = eval { $action->(); 1 };
    if ( !$ok ) {
        my $error = $@ || "Unknown error";
        chomp $error;
        user_message(
              "PrimerMapper stopped before completing this action.\n\n"
            . "Reason:\n$error" );
    }
}

sub validate_numeric_field {
    my ( $label, $value, $allow_zero ) = @_;
    return "$label must be numeric." unless defined $value && $value =~ /^\d+(?:\.\d+)?$/;
    return "$label must be greater than zero." if !$allow_zero && $value <= 0;
    return;
}

sub validate_format {
    my $fmt = shift;
    return unless defined $fmt && $fmt =~ /^\.(?:png|gif)$/i;
    return "Output format must be .png or .gif.";
}

sub validate_regular_parameters {
    my @errors;
    push @errors, grep { defined $_ } (
        validate_numeric_field( "Minimum product size", $min_prod_size, 0 ),
        validate_numeric_field( "Maximum product size", $max_prod_size, 0 ),
        validate_numeric_field( "Minimum primer length", $kmer_min, 0 ),
        validate_numeric_field( "Maximum primer length", $kmer_max, 0 ),
        validate_numeric_field( "Lower GC%", $lower_gc, 1 ),
        validate_numeric_field( "Upper GC%", $higher_gc, 1 ),
        validate_numeric_field( "Lower Tm", $lower_tm, 0 ),
        validate_numeric_field( "Upper Tm", $upper_tm, 0 ),
        validate_numeric_field( "GC clamp", $clamp, 1 ),
        validate_numeric_field( "Salt concentration", $salt, 0 ),
        validate_format($format),
    );

    push @errors, "Minimum product size cannot be greater than maximum product size."
      if $min_prod_size =~ /^\d/ && $max_prod_size =~ /^\d/ && $min_prod_size > $max_prod_size;
    push @errors, "Minimum primer length cannot be greater than maximum primer length."
      if $kmer_min =~ /^\d/ && $kmer_max =~ /^\d/ && $kmer_min > $kmer_max;
    push @errors, "Lower GC% cannot be greater than upper GC%."
      if $lower_gc =~ /^\d/ && $higher_gc =~ /^\d/ && $lower_gc > $higher_gc;
    push @errors, "Lower Tm cannot be greater than upper Tm."
      if $lower_tm =~ /^\d/ && $upper_tm =~ /^\d/ && $lower_tm > $upper_tm;
    push @errors, "GC clamp must be 0 or 1." if defined $clamp && $clamp !~ /^(?:0|1)$/;

    if (@errors) {
        user_message( "Please correct these primer-design settings:\n\n" . join( "\n", @errors ) );
        return 0;
    }
    return 1;
}

sub validate_snp_parameters {
    my @errors;
    push @errors, grep { defined $_ } (
        validate_numeric_field( "5' upstream SNP window", $five_prime_SNP, 0 ),
        validate_numeric_field( "3' downstream SNP window", $three_prime_SNP, 0 ),
        validate_numeric_field( "Minimum distance from SNP", $snp_distance, 1 ),
        validate_numeric_field( "Minimum SNP primer length", $kmer_min_SNP, 0 ),
        validate_numeric_field( "Maximum SNP primer length", $kmer_max_SNP, 0 ),
        validate_numeric_field( "Lower SNP GC%", $lower_gc_SNP, 1 ),
        validate_numeric_field( "Upper SNP GC%", $higher_gc_SNP, 1 ),
        validate_numeric_field( "Lower SNP Tm", $lower_tm_SNP, 0 ),
        validate_numeric_field( "Upper SNP Tm", $upper_tm_SNP, 0 ),
        validate_numeric_field( "SNP GC clamp", $clamp_SNP, 1 ),
        validate_numeric_field( "SNP salt concentration", $salt_SNP, 0 ),
        validate_format($format),
    );

    push @errors, "Minimum SNP primer length cannot be greater than maximum SNP primer length."
      if $kmer_min_SNP =~ /^\d/ && $kmer_max_SNP =~ /^\d/ && $kmer_min_SNP > $kmer_max_SNP;
    push @errors, "Lower SNP GC% cannot be greater than upper SNP GC%."
      if $lower_gc_SNP =~ /^\d/ && $higher_gc_SNP =~ /^\d/ && $lower_gc_SNP > $higher_gc_SNP;
    push @errors, "Lower SNP Tm cannot be greater than upper SNP Tm."
      if $lower_tm_SNP =~ /^\d/ && $upper_tm_SNP =~ /^\d/ && $lower_tm_SNP > $upper_tm_SNP;
    push @errors, "SNP GC clamp must be 0 or 1." if defined $clamp_SNP && $clamp_SNP !~ /^(?:0|1)$/;
    push @errors, "The 5' upstream SNP window must be larger than the minimum SNP distance."
      if $five_prime_SNP =~ /^\d/ && $snp_distance =~ /^\d/ && $five_prime_SNP <= $snp_distance;

    if (@errors) {
        user_message( "Please correct these SNP primer-design settings:\n\n" . join( "\n", @errors ) );
        return 0;
    }
    return 1;
}

sub fasta_file_has_descriptor {
    my $file = shift;
    open my $fh, "<", $file or return ( 0, "Cannot read $file: $!" );
    while ( my $line = <$fh> ) {
        next if $line =~ /^\s*$/;
        close $fh;
        return $line =~ /^>/ ? ( 1, undef ) : ( 0, "The FASTA file must start with a descriptor line beginning with '>'." );
    }
    close $fh;
    return ( 0, "The FASTA file is empty." );
}

sub read_fasta_sequences {
    my $file = shift;
    my ( $has_descriptor, $descriptor_error ) = fasta_file_has_descriptor($file);
    return ( undef, $descriptor_error ) unless $has_descriptor;

    my @records;
    my $seqio;
    my $ok = eval {
        $seqio = Bio::SeqIO->new( -file => $file, -format => "fasta" );
        while ( my $seqobj = $seqio->next_seq() ) {
            push @records, $seqobj;
        }
        1;
    };
    return ( undef, $@ || "BioPerl could not parse this FASTA file." ) unless $ok;
    return ( undef, "No FASTA sequences were found in this file." ) unless @records;
    return ( \@records, undef );
}

sub primer3_available {
    my $program = $^O eq 'MSWin32' ? "primer3_core.exe" : "primer3_core";
    for my $dir ( split /:/, $ENV{PATH} || "" ) {
        return $program if -x "$dir/$program";
    }
    user_message( "PrimerMapper could not find $program in your PATH.\n\nInstall primer3 and restart PrimerMapper." );
    return;
}

sub primer3_explain_text {
    my $file = shift;
    return "" unless defined $file && -f $file;
    open my $fh, "<", $file or return "";
    my @explain;
    while ( my $line = <$fh> ) {
        chomp $line;
        push @explain, $line if $line =~ /(?:EXPLAIN|NUM_RETURNED|ERROR)/i;
    }
    close $fh;
    return join( "\n", @explain );
}

sub reset_regular_run_state {
    @array_length = ();
    @array_name   = ();
    @gene_length  = ();
    @increment    = ();
    $foo          = 0;
    unlink $out_single if -f $out_single;
    unlink glob "GRAPHIC_*.txt";
}

sub reset_snp_run_state {
    @array_length_SNP = ();
    @array_name_SNP   = ();
    @gene_length_SNP  = ();
    @increment_SNP    = ();
    @SNP_position     = ();
    $foo_SNP          = 0;
    $snp              = undef;
    $degen            = undef;
    unlink $out_single_SNP if -f $out_single_SNP;
    unlink glob "GRAPHIC_*.txt";
}

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
$gc_clamp_snp_->g_grid( -column => 7, -row => 5, -sticky => "we" );

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
        run_gui_action(
            sub {
                return unless calculate_SNP();
                graphics_all_primers_SNP( \@array_length_SNP, \@array_name_SNP,
                    \@gene_length_SNP, $format, $snp );
                graphics_single_view_SNP( \@gene_length_SNP, $format, $out_single_SNP );
            }
        );
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
    -textvariable => \$max_prod_size
);
$five_primer_ending_->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $three_primer_ending_ = $out_frame_primers->new_ttk__entry(
    -width        => 3,
    -textvariable => \$min_prod_size
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
$add_gc_clamp_->g_grid( -column => 7, -row => 2, -sticky => "we" );

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
        run_gui_action(
            sub {
                return unless calculate();
                graphics_all_primers( \@array_length, \@array_name, \@gene_length,
                    $format );
                graphics_single_view( \@gene_length, $format, $out_single );
            }
        );
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
 Function:  uses NIH e-utilities to retrieve remotely sequences from GenBank
 Returns :  fastA formatted sequences 
=cut

sub get_sequences {
    $seq_in =~ s/\s//g;

    # BUG FIX: original used http:// which NCBI redirects to https://,
    # causing LWP::Simple::get() to fail silently and return undef.
    my $url =
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
      . $seq_in
      . '&rettype=fasta&retmode=text';

    my $fileList = get($url);

    # BUG FIX: original did not check whether the request succeeded.
    # If get() returns undef the open/print would silently write nothing.
    unless ( defined $fileList && length $fileList ) {
        Tkx::tk___messageBox(
            -message => "Failed to retrieve sequences from NCBI. "
              . "Check your internet connection and accession IDs." );
        return;
    }
    unless ( $fileList =~ /^\s*>/ ) {
        user_message(
              "NCBI returned data, but it was not FASTA format.\n\n"
            . "Check the accession ID(s) and try again." );
        return;
    }

    my $new_seqs = 'sequences.fasta';

    # BUG FIX: original used bareword filehandle and bare die with no message.
    open my $seqs_fh, ">", $new_seqs
      or do {
        Tkx::tk___messageBox( -message => "Cannot write $new_seqs: $!" );
        return;
      };
    print $seqs_fh $fileList;
    close $seqs_fh;

    # BUG FIX: auto-populate the FASTA file path in the GUI entry widget
    $fasta = $new_seqs;
    print "\nSequences collected and saved to $new_seqs.\n";
}

####Collect SNP Seqs remotely

=head1 get_sequences_SNP
 Title   :  get_sequences_SNP
 Usage   :  -command => sub { get_sequences_SNP(); }
 Function:  uses NIH e-utilities to retrieve remotely SNP sequences from dbSNP
 Returns :  SNP formatted sequences 
=cut

sub get_sequences_SNP {
    $seq_in_SNP =~ s/\s//g;

    # BUG FIX: http -> https (same NCBI redirect issue as get_sequences)
    my $url =
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id='
      . $seq_in_SNP
      . '&rettype=fasta&retmode=text';

    my $fileList = get($url);

    unless ( defined $fileList && length $fileList ) {
        Tkx::tk___messageBox(
            -message => "Failed to retrieve SNP sequences from NCBI. "
              . "Check your internet connection and rs IDs." );
        return;
    }
    unless ( $fileList =~ /^\s*>/ ) {
        user_message(
              "NCBI returned data, but it was not SNP FASTA format.\n\n"
            . "Check the rs ID(s) and try again." );
        return;
    }

    my $new_seqs = 'SNP_sequences.fasta';
    open my $seqs_fh, ">", $new_seqs
      or do {
        Tkx::tk___messageBox( -message => "Cannot write $new_seqs: $!" );
        return;
      };
    print $seqs_fh $fileList;
    close $seqs_fh;

    # BUG FIX: auto-populate SNP FASTA path in GUI
    $fasta_SNP = $new_seqs;
    print "\nSNP Sequences collected and saved to $new_seqs.\n";
}

####Load Seqs into PrimerMapper

=head1 loadfile
 Title   :  loadfile
 Usage   :  $fasta = Tkx::tk___getOpenFile();
 Function:  loads sequences into PrimerMapper program 
=cut

sub loadfile {
    $fasta = Tkx::tk___getOpenFile();
    return unless defined $fasta && -f $fasta;    # user cancelled dialog

    # BUG FIX: original had typo "specificty.txt" (missing 'i').
    # Also used append (>>) so repeated loads would accumulate stale entries.
    $outfile = 'specificity.txt';
    open my $out_fh, ">", $outfile
      or do {
        user_message("Cannot write $outfile:\n\n$!");
        return;
      };
    open my $in_fh, "<", $fasta
      or do {
        user_message("Cannot read $fasta:\n\n$!");
        return;
      };

    while (<$in_fh>) {
        print $out_fh $_ if ( $_ !~ /(>.+)/ );
    }
    close $in_fh;
    close $out_fh;

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
    return unless defined $fasta_SNP && -f $fasta_SNP;

    # BUG FIX: original had a trailing space in the filename
    # ">>$outfile_SNP " which created a file named "specificty_SNP.txt "
    # (with space) on some systems. Also typo "specificty" fixed.
    $outfile_SNP = 'specificity_SNP.txt';
    open my $out_snp_fh, ">", $outfile_SNP
      or do {
        user_message("Cannot write $outfile_SNP:\n\n$!");
        return;
      };
    open my $in_snp_fh, "<", $fasta_SNP
      or do {
        user_message("Cannot read $fasta_SNP:\n\n$!");
        return;
      };

    while (<$in_snp_fh>) {
        print $out_snp_fh $_ if ( $_ !~ />.+/ );
    }
    close $in_snp_fh;
    close $out_snp_fh;

    print "\nSNP Sequences Loaded.\n";
}

####Parse fasta sequence using BioSeqIO
sub calculate {

    reset_regular_run_state();

    unless ( defined $fasta && -f $fasta ) {
        user_message("Please load or retrieve a FASTA file first.");
        return 0;
    }

    return 0 unless validate_regular_parameters();

    my $primer3 = primer3_available();
    return 0 unless $primer3;

    my ( $records, $read_error ) = read_fasta_sequences($fasta);
    if ($read_error) {
        user_message( "PrimerMapper could not read the FASTA file:\n\n$read_error" );
        return 0;
    }

    print "\nDesigning and validating primers....\n";

    my $features_written = 0;
    my @notices;

    for my $seqobj (@$records) {
        my $sequence = $seqobj->seq();
        $sequence =~ s/[^agctn]//ig;
        $sequence =~ s/\s//g;
        $sequence =~ tr/a-z/A-Z/;
        $id = $seqobj->id();
        $id =~ tr/a-zA-Z0-9//cd;
        my $len_seq  = length $sequence;

        if ( $len_seq < $min_prod_size ) {
            push @notices,
              "$id is $len_seq bp, shorter than the minimum product size ($min_prod_size bp).";
            next;
        }

        my $id_uniq  = $len_seq . "'" . $id;
        my $len_uniq = $id . "'" . $len_seq;
        push @increment, $len_seq;
        $foo = eval join '+', @increment;

        #declare output file
        $out_image = "GRAPHIC_$id.txt";

        my $p3 = $id . "primer3_output";
        open my $fh, ">", $p3 or die "Cannot write $p3: $!\n";
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

        # BUG FIX: the original did not close $p3_fh before passing the file
        # to primer3_core, meaning the last buffer flush was not guaranteed —
        # primer3 could read a truncated input on some platforms.
        open my $p3_fh, ">", $temp_seq or die "Cannot write $temp_seq: $!\n";
        print $p3_fh $new_sequence;
        close $p3_fh;    # BUG FIX: explicit close before system() call

        if ( $^O eq 'MSWin32' ) {
            my $command = "primer3_core.exe -format_output < $temp_seq >> $p3";
            system($command);
        }
        else {
            my $command = "primer3_core -format_output < $temp_seq >> $p3";
            system($command);
        }

        if ( $? != 0 ) {
            push @notices, "primer3 failed for $id. Check that primer3_core is installed and working.";
            next;
        }

        open $fh, "<", $p3 or die "Cannot read $p3: $!\n";
        while ( my $row = <$fh> ) {
            if ( $row =~
/LEFT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                my $add = $1 + $2;
                open my $fh_oligo, ">>", $out_image
                  or die "Cannot write $out_image: $!\n";
                print $fh_oligo( $1 . "\t" . $2 . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_oligo;
                open my $fh_outs, ">>", $out_single
                  or die "Cannot write $out_single: $!\n";
                print $fh_outs( $1 . "\t" . $foo . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_outs;
                $features_written++;
            }
            elsif ( $row =~
/RIGHT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                my $minus = $1 - $2;
                open my $fh_oligo, ">>", $out_image
                  or die "Cannot write $out_image: $!\n";
                print $fh_oligo( $1 . "\t" . $2 . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_oligo;
                open my $fh_outs, ">>", $out_single
                  or die "Cannot write $out_single: $!\n";
                print $fh_outs( $1 . "\t" . $foo . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_outs;
                $features_written++;
            }
        }
        close $fh;

        if ( !-s $out_image ) {
            my $explain = primer3_explain_text($p3);
            push @notices, "No primers were found for $id." . ( $explain ? "\n$explain" : "" );
        }
    }

    if ( !$features_written ) {
        my $detail = @notices ? "\n\nDetails:\n" . join( "\n\n", @notices ) : "";
        user_message(
              "PrimerMapper did not find any primer pairs to draw.\n\n"
            . "Try a longer FASTA sequence, a smaller product-size range, or relaxed Tm/GC settings."
            . $detail );
        return 0;
    }

    if (@notices) {
        user_message(
              "Primer design completed, but some sequences produced no primers.\n\n"
            . join( "\n\n", @notices ) );
    }

    print "\nPrimer Design completed.\n";
    return 1;
}

########################################
######## START SNP PRIMER DESIGN
########################################

sub calculate_SNP {

    reset_snp_run_state();

    unless ( defined $fasta_SNP && -f $fasta_SNP ) {
        user_message("Please load or retrieve a SNP FASTA file first.");
        return 0;
    }

    return 0 unless validate_snp_parameters();

    my $primer3 = primer3_available();
    return 0 unless $primer3;

    $fasta_SNP_NEW = "processed_SNP.fasta";
    open my $inner_fh, '<', $fasta_SNP
      or do {
        user_message("Cannot read $fasta_SNP:\n\n$!");
        return 0;
      };
    open my $outer_fh, '>', $fasta_SNP_NEW
      or do {
        user_message("Cannot write $fasta_SNP_NEW:\n\n$!");
        close $inner_fh;
        return 0;
      };

    while (<$inner_fh>) {
        s/ //g;
        print $outer_fh $_;
    }
    close $outer_fh;
    close $inner_fh;

    # BUG FIX: original unlinked $fasta_SNP unconditionally. If the source
    # file was on a different filesystem (e.g. user's home dir) this would
    # destroy their original data. Now we leave the source file untouched.

    my ( $records, $read_error ) = read_fasta_sequences($fasta_SNP_NEW);
    if ($read_error) {
        user_message( "PrimerMapper could not read the SNP FASTA file:\n\n$read_error" );
        return 0;
    }

    print "\nPreprocessing SNP file....\n....Designing and validating SNP primers....\n";

    my $features_written = 0;
    my @notices;

    for my $seqobj (@$records) {
        my $sequence = $seqobj->seq();
        $id_SNP = $seqobj->id();
        $snp   = undef;
        $degen = undef;
        $sequence =~ s/\s//g;
        $sequence =~ tr/a-z/A-Z/;

        if ( $id_SNP =~ m/pos=(\d+)\.len/i ) {
            $snp = $1;
        }

        unless ( defined $snp && $snp =~ /^\d+$/ ) {
            push @notices,
              "$id_SNP does not have a parseable SNP position in the FASTA header. Expected a header containing pos=N.len.";
            next;
        }
        if ( $snp < 1 || $snp > length($sequence) ) {
            push @notices,
              "$id_SNP has SNP position $snp, but the sequence is "
              . length($sequence)
              . " bp.";
            next;
        }

        my $raw_snp_base = substr( $sequence, $snp - 1, 1 );

        ###################
        # Detect degenerate IUPAC base at the SNP position declared in header
        if    ( $raw_snp_base eq "R" ) { $degen = "type_1"; }
        elsif ( $raw_snp_base eq "Y" ) { $degen = "type_2"; }
        elsif ( $raw_snp_base eq "S" ) { $degen = "type_3"; }
        elsif ( $raw_snp_base eq "W" ) { $degen = "type_4"; }
        elsif ( $raw_snp_base eq "K" ) { $degen = "type_5"; }
        elsif ( $raw_snp_base eq "M" ) { $degen = "type_6"; }
        elsif ( $raw_snp_base eq "B" ) { $degen = "type_7"; }
        elsif ( $raw_snp_base eq "D" ) { $degen = "type_8"; }
        elsif ( $raw_snp_base eq "H" ) { $degen = "type_9"; }
        elsif ( $raw_snp_base eq "V" ) { $degen = "type_10"; }

        unless ( defined $degen ) {
            push @notices,
              "$id_SNP does not contain an IUPAC SNP ambiguity code at position $snp. "
              . "Expected R, Y, S, W, K, M, B, D, H, or V at that coordinate.";
            next;
        }

        # Degenerate SNP code:
        ##-->R  A or G
        ##-->Y  C or T
        ##-->S  G or C
        ##-->W  A or T
        ##-->K  G or T
        ##-->M  A or C
        ##-->B  C or G or T
        ##-->D  A or G or T
        ##-->H  A or C or T
        ##-->V  A or C or G
        ###################

        $sequence =~ s/[RYSWKMBDHV]/N/ig;
        $sequence =~ s/[^agctn]//ig;

        $id_SNP =~ tr/a-zA-Z0-9//cd;
        my $len_seq      = length $sequence;

        if ( $len_seq < ( $kmer_min_SNP * 2 + 1 ) ) {
            push @notices,
              "$id_SNP is $len_seq bp, too short for the current minimum SNP primer length ($kmer_min_SNP bp).";
            next;
        }
        if ( $len_seq < 100 ) {
            push @notices,
              "$id_SNP is $len_seq bp, shorter than the SNP product-size minimum used by PrimerMapper (100 bp).";
            next;
        }

        my $id_uniq_SNP  = $len_seq . "'" . $id_SNP;
        my $len_uniq_SNP = $id_SNP . "'" . $len_seq;
        push @increment_SNP, $len_seq;
        $foo_SNP = eval join '+', @increment_SNP;

        #declare output file
        $out_image_SNP = "GRAPHIC_$id_SNP.txt";

        # BUG FIX: original opened AS_SNP (bareword) but never closed it
        # inside the loop, leaking file descriptors. Also opened with
        # append ">>$outputfile_SNP_AS " (trailing space in name). Fixed.
        {
            open my $as_snp_fh, ">>", $outputfile_SNP_AS
              or die "Cannot write $outputfile_SNP_AS: $!\n";
            print $as_snp_fh
              "Header\tTm(degC)\tAllele Specific Primer Sequences\t\t\tGC%\n";
            close $as_snp_fh;
        }

        ######## SNP PRIMERS

        my $short_id = substr $id_SNP, 0, 20;
        my $p3_SNP   = $short_id . "primer3_output_SNP";
        open my $fh_SNP, ">", $p3_SNP or die "Cannot write $p3_SNP: $!\n";
        print $fh_SNP "Sequence Name: " . $short_id . "\n";
        close $fh_SNP;

        my $anchor      = $snp - $five_prime_SNP;
        $anchor = 0 if $anchor < 0;
        my $dist_before = $five_prime_SNP - $snp_distance;
        my $dist_after  = $snp + $snp_distance;

        if ( $dist_before <= 0 || $dist_after >= $len_seq ) {
            push @notices,
              "$id_SNP has too little sequence around the SNP for the current window/distance settings. "
              . "Try lowering the minimum distance from SNP or using a longer SNP FASTA sequence.";
            next;
        }

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

        # BUG FIX: same unclosed-before-system() issue as calculate()
        open my $p3_fh_SNP, ">", $temp_seq_SNP
          or die "Cannot write $temp_seq_SNP: $!\n";
        print $p3_fh_SNP $new_sequence_SNP;
        close $p3_fh_SNP;    # BUG FIX: explicit close before system()

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

        if ( $? != 0 ) {
            push @notices, "primer3 failed for SNP sequence $id_SNP. Check that primer3_core is installed and working.";
            next;
        }

        open $fh_SNP, "<", $p3_SNP or die "Cannot read $p3_SNP: $!\n";
        while ( my $row = <$fh_SNP> ) {
            if ( $row =~
/LEFT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length_SNP, $len_uniq_SNP;
                push @array_name_SNP,   $id_uniq_SNP;
                my $add = $1 + $2;
                open my $fh_oligo_SNP, ">>", $out_image_SNP
                  or die "Cannot write $out_image_SNP: $!\n";
                print $fh_oligo_SNP(
                    $snp . "\t" . $2 . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_oligo_SNP;
                open my $fh_outs_SNP, ">>", $out_single_SNP
                  or die "Cannot write $out_single_SNP: $!\n";
                print $fh_outs_SNP(
                    $1 . "\t" . $foo_SNP . "\t" . $1 . "\t" . $add . "\n" );
                close $fh_outs_SNP;
                $features_written++;
            }
            elsif ( $row =~
/RIGHT PRIMER\s+(\d*)\s+(\d+)\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
              )
            {
                push @array_length_SNP, $len_uniq_SNP;
                push @array_name_SNP,   $id_uniq_SNP;
                my $minus = $1 - $2;
                open my $fh_oligo_SNP, ">>", $out_image_SNP
                  or die "Cannot write $out_image_SNP: $!\n";
                print $fh_oligo_SNP(
                    $snp . "\t" . $2 . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_oligo_SNP;
                open my $fh_outs_SNP, ">>", $out_single_SNP
                  or die "Cannot write $out_single_SNP: $!\n";
                print $fh_outs_SNP(
                    $1 . "\t" . $foo_SNP . "\t" . $1 . "\t" . $minus . "\n" );
                close $fh_outs_SNP;
                $features_written++;
            }
        }
        close $fh_SNP;

        if ( !-s $out_image_SNP ) {
            my $explain = primer3_explain_text($p3_SNP);
            push @notices, "No SNP-flanking primer pairs were found for $id_SNP."
              . ( $explain ? "\n$explain" : "" );
        }

        print "\nSNP Primer Design completed.\n";

        ######## ALLELE SPECIFIC PRIMER

        ####################______Allele-Specific_BEGIN

        for ( my $p = $kmer_min_SNP ; $p <= $kmer_max_SNP ; $p++ ) {
            my $subseq = substr( $sequence, $snp - 1, $p );

            #get complement (reverse strand allele-specific oligo)
            my $revFA = reverse($subseq);
            $revFA =~ tr/ATGCatgc/TACGtacg/;

            #Count Gs and Cs
            my $countGCA = ( $subseq =~ tr/GCgc// );

            #Calculate percent GC
            my $percentGCA        = 100 * $countGCA / $p;
            my $percentGCroundedA = sprintf( "%0.1f", $percentGCA );

            #calculate Tm
            if ( $p <= 36 ) {
                $Tm3_SNP = calcTm( $revFA, 100, $salt_SNP );
            }
            else {
                $Tm3_SNP =
                  calclongTm( $revFA, 100, $salt_SNP, $percentGCroundedA );
            }
            my $TmroundedA = defined $Tm3_SNP ? sprintf( "%0.1f", $Tm3_SNP ) : "N/A";

            # BUG FIX: original called open() inside if/elsif conditions,
            # meaning the file was opened each iteration of the loop rather
            # than once. This caused descriptor leaks and redundant I/O.
            # Now open once, print, close.
            open my $as_fh, ">>", $outputfile_SNP_AS
              or die "Cannot write $outputfile_SNP_AS: $!\n";

            if ( $degen eq "type_1" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_2" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.G\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_3" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.C\t$revFA.G\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_4" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_5" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.C\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_6" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_7" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.G\t$revFA.C\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_8" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_9" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.A\t$percentGCroundedA\n";
            }
            elsif ( $degen eq "type_10" ) {
                print $as_fh
                  "$short_id\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.C\t$percentGCroundedA\n";
            }
            close $as_fh;
        }

        ###################

        my $string_to_snp = $snp - 2;
        my $new_seq = substr( $sequence, 1, $string_to_snp );

        for ( my $u = $kmer_max_SNP ; $u >= $kmer_min_SNP ; $u-- ) {
            my $fwd_subseq = substr( $new_seq, -$u );

            #Count Gs and Cs  (BUG FIX: original tr operated on $_ not $fwd_subseq)
            my $countGCAB = ( $fwd_subseq =~ tr/GCgc// );

            #Calculate percent GC
            my $percentGCAB        = 100 * $countGCAB / $u;
            my $percentGCroundedAB = sprintf( "%0.1f", $percentGCAB );

            #calculate Tm
            if ( $u <= 36 ) {
                $Tm4_SNP = calcTm( $fwd_subseq, 100, $salt_SNP );
            }
            else {
                $Tm4_SNP =
                  calclongTm( $fwd_subseq, 100, $salt_SNP, $percentGCroundedAB );
            }
            my $TmroundedAB =
              defined $Tm4_SNP ? sprintf( "%0.1f", $Tm4_SNP ) : "N/A";

            # BUG FIX: same open-in-condition pattern fixed as above
            open my $as_fh2, ">>", $outputfile_SNP_AS
              or die "Cannot write $outputfile_SNP_AS: $!\n";

            if ( $degen eq "type_1" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.G\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_2" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.C\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_3" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.G\t$fwd_subseq.C\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_4" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_5" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.G\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_6" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.C\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_7" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.C\t$fwd_subseq.G\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_8" ) {
                # BUG FIX: original said "Reverse:" here for a forward primer — typo corrected
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.G\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_9" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.C\t$fwd_subseq.T\t$percentGCroundedAB\n";
            }
            elsif ( $degen eq "type_10" ) {
                print $as_fh2
                  "$short_id\t$TmroundedAB \tForward: $fwd_subseq.A\t$fwd_subseq.C\t$fwd_subseq.G\t$percentGCroundedAB\n";
            }
            close $as_fh2;
        }

        ################
    }
    ################### ----> Allele-Specific_END

    if ( !$features_written ) {
        my $detail = @notices ? "\n\nDetails:\n" . join( "\n\n", @notices ) : "";
        user_message(
              "PrimerMapper did not find any SNP-flanking primer pairs to draw.\n\n"
            . "Try a longer SNP FASTA sequence, lowering the minimum distance from SNP, "
            . "or relaxing Tm/GC/primer length settings."
            . $detail );
        return 0;
    }

    if (@notices) {
        user_message(
              "SNP primer design completed, but some sequences produced no SNP primer pairs.\n\n"
            . join( "\n\n", @notices ) );
    }

    print "\nAllele-specific Primer Design completed.\n";
    return 1;
}

Tkx::MainLoop();

####################################

=pod 
=WEB-BLAST
Note: must have web_blast.pl script in PATH for BLAST feature 
=cut

####################################
####BLAST SUBROUTINE####
####################################

=head1 web_blast_ncbi
 Title   :  web_blast_ncbi
 Usage   :  -command => sub { web_blast_ncbi(); }
 Function:  runs remote BLAST from PrimerMapper using the script web_blast.pl
 Returns :  BLAST result in single text file 
=cut

sub web_blast_ncbi {

    # BUG FIX: original did not check whether a BLAST database was selected
    unless ( defined $blaster ) {
        Tkx::tk___messageBox( -message => "Please select a BLAST database first." );
        return;
    }

    print "\nBLAST Analysis of $blaster is running....\n";

    $thetext = $text->get( "1.0", "end" );

    unless ( $thetext =~ /\S/ ) {
        Tkx::tk___messageBox(
            -message => "Please paste primer sequences in FASTA format in the text box." );
        return;
    }

    my $blastfile = "blast_primers.fasta";

    # BUG FIX: original used append >> so repeated BLASTs accumulated
    # sequences from previous runs in the same file.
    open my $blast_fh, ">", $blastfile
      or do {
        Tkx::tk___messageBox( -message => "Cannot write $blastfile: $!" );
        return;
      };
    print $blast_fh $thetext;
    close $blast_fh;

    # BUG FIX: locate web_blast.pl relative to the script rather than
    # relying on it being in PATH (it wasn't added to PATH in the install).
    my $web_blast = "$Bin/PrimerMapper/lib/web_blast.pl";
    unless ( -f $web_blast ) {
        # Fallback: try cwd
        $web_blast = "web_blast.pl";
    }

    my $command = "perl \"$web_blast\" blastn $blaster $blastfile";
    system($command);

    print "\nBLAST Analysis completed.\n";
}

####################################
1;
