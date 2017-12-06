#
# Utility module for PrimerMapper
#
# Please direct questions and support issues to <https://github.com/dohalloran/PrimerMapper/issues>
#
# Author: Damien O'Halloran, The George Washington University, 2015
#
# GNU GENERAL PUBLIC LICENSE
#
# POD documentation before the code

=head1 NAME
PrimerMapperUtilities - runs various subroutines mainly for primer design

=head2 SYNOPSIS
  use base 'Exporter';
  use Cwd;
  
  for Usage, run by creating a PrimerMapper object within PrimerMapper_driver.pl
  my $tmp = PrimerMapper->new();
  
=head3 DESCRIPTION
The ouput from this package produces graphical primer::sequence
maps of primer distribution across each sequence and the entire input file

=head4 REFERENCES
The design of the code is based around the following:
Rychlik, W. & Rhoads, R. E. A computer program for choosing optimal
oligonucleotides for filter hybridization, sequencing and in vitro amplification of DNA.
Nucleic Acids Res. 17, 8543-8551 (1989)
Li, K. et al. Novel computational methods for increasing PCR primer design
effectiveness in directed sequencing. BMC Bioinformatics 9, 191-2105-9-191 (2008)

=head5 Support 
All contributions are welcome

=head6 Reporting Bugs
Report bugs to the PrimerMapper bug tracking system to help keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:
  https://github.com/dohalloran/PrimerMapper/issues
  
=head7 APPENDIX
The rest of the documentation details each of the object methods.


=cut

package PrimerMapperUtilities;

use strict;
use warnings;
use base 'Exporter';
use Cwd;

our @EXPORT =
  qw/ calcTm calclongTm clean_up /;

####Global Hashes
my %nn_s = (
    "AA" => 240,
    "AC" => 173,
    "AG" => 208,
    "AT" => 239,
    "AN" => 215,
    "CA" => 129,
    "CC" => 266,
    "CG" => 278,
    "CT" => 208,
    "CN" => 220,
    "GA" => 135,
    "GC" => 267,
    "GG" => 266,
    "GT" => 173,
    "GN" => 210,
    "TA" => 169,
    "TC" => 135,
    "TG" => 129,
    "TT" => 240,
    "TN" => 168,
    "NA" => 168,
    "NC" => 210,
    "NG" => 220,
    "NT" => 215,
    "NN" => 203,
    "aa" => 240,
    "ac" => 173,
    "ag" => 208,
    "at" => 239,
    "an" => 215,
    "ca" => 129,
    "cc" => 266,
    "cg" => 278,
    "ct" => 208,
    "cn" => 220,
    "ga" => 135,
    "gc" => 267,
    "gg" => 266,
    "gt" => 173,
    "gn" => 210,
    "ta" => 169,
    "tc" => 135,
    "tg" => 129,
    "tt" => 240,
    "tn" => 168,
    "na" => 168,
    "nc" => 210,
    "ng" => 220,
    "nt" => 215,
    "nn" => 203
);
my %nn_h = (
    "AA" => 91,
    "AC" => 65,
    "AG" => 78,
    "AT" => 86,
    "AN" => 80,
    "CA" => 58,
    "CC" => 110,
    "CG" => 119,
    "CT" => 78,
    "CN" => 91,
    "GA" => 56,
    "GC" => 111,
    "GG" => 110,
    "GT" => 65,
    "GN" => 85,
    "TA" => 60,
    "TC" => 56,
    "TG" => 58,
    "TT" => 91,
    "TN" => 66,
    "NA" => 66,
    "NC" => 85,
    "NG" => 91,
    "NT" => 80,
    "NN" => 80,
    "aa" => 91,
    "ac" => 65,
    "ag" => 78,
    "at" => 86,
    "an" => 80,
    "ca" => 58,
    "cc" => 110,
    "cg" => 119,
    "ct" => 78,
    "cn" => 91,
    "ga" => 56,
    "gc" => 111,
    "gg" => 110,
    "gt" => 65,
    "gn" => 85,
    "ta" => 60,
    "tc" => 56,
    "tg" => 58,
    "tt" => 91,
    "tn" => 66,
    "na" => 66,
    "nc" => 85,
    "ng" => 91,
    "nt" => 80,
    "nn" => 80
);

####################################
####################################
=head1 calcTm 
 Title   :  calcTm 
 Usage   :  $Tm = calcTm( $_, $DNA_conc, $salt );
 Function:  calculates Tm if primer is <=36bps 
 Returns :  primer Tm
=cut

sub calcTm {
    my $sequence = shift;
    my $DNA_nM   = shift;
    my $K_mM     = shift;
    my $dH       = 0;
    my $dS       = 108;
    my $seq_len  = length($sequence);
    my $i;

    # Compute dH and dS
    for ( $i = 0 ; $i < ( $seq_len - 1 ) ; $i++ ) {
        my $pair = substr( $_, $i, 2 );
        $dH += $nn_h{$pair};
        $dS += $nn_s{$pair};
    }

    $dH *= -100.0;
    $dS *= -0.1;

    return $dH / ( $dS + 1.987 * log( $DNA_nM / 4000000000.0 ) ) - 273.15 +
      16.6 * ( log( $K_mM / 1000.0 ) / log(10) );
}

####################################
=head1 calclongTm 
 Title   :  calclongTm 
 Usage   :  $Tm = calclongTm( $_, $DNA_conc, $salt, $percentGCrounded );
 Function:  calculates Tm if primer is >36bps 
 Returns :  primer Tm
=cut

sub calclongTm {
    my $sequence           = shift;
    my $DNA_nM             = shift;
    my $K_mM               = shift;
    my $imported_gcPercent = shift;
    my $new_GC             = $imported_gcPercent / 100;
    return 81.5 + ( 16.6 * ( log( $K_mM / 1000.0 ) / log(10) ) ) +
      ( 41.0 * $new_GC ) - ( 600.0 / length($sequence) );
}

####################################
=head1 clean_up 
 Title   :  clean_up 
 Usage   :  -command => sub { clean_up(); }
 Function:  removes temporary files from current working directory
 Returns :  deletes .txt files generated during program execution
=cut

sub clean_up {
    my $dir = cwd();
    unlink glob "$dir/*.txt";

    print "\nDirectory has been cleaned.\n";
}

####################################

1;
