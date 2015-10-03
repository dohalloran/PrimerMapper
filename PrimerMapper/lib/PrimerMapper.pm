package PrimerMapper;

#############################################################
#
# PrimerMapper
# Written by Damien O'Halloran
# The George Washington University
# Summer 2015
# readme.txt file for usage
#
#############################################################

use warnings;
use strict;

#### Load Modules
use Tkx;
use LWP::Simple;
use Cwd;
use Bio::SeqIO;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use List::Util 'max';
use base 'Exporter';
####

#########################

our @EXPORT =
  qw( get_sequences get_sequences_SNP loadfile loadfile_SNP calculate calculate_SNP match_positions 
  fuzzy_pattern make_approximate calclongTm calcdG calcTm selfie graphics_all_primers graphics_single_view 
  graphics_all_primers_SNP graphics_single_view_SNP primer_dimer primer_dimer_SNP clean_up web_blast_ncbi );

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

=pod
 
=head1 UTILITY
 
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

Author: Damien O'Halloran, The George Washington University, 2015

To run, execute as follows: 
>perl PrimerMapper_driver.pl 
  
=cut

####Primer Design Defaults
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
my $html          = "Y";
###########################

####Strings, Arrays, and Hashes
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

####SNP Primer Design Defaults
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

####SNP Strings and Hashes
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

####Load the Tkx interface
my $mw = Tkx::widget->new(".");
$mw->g_wm_title("PrimerMapper");
my $frm = $mw->new_ttk__frame( -padding => "3 3 12 12" );
$frm->g_grid( -column => 0, -row => 0, -sticky => "nwes" );

my $lf = $frm->new_ttk__labelframe( -text => "INPUT" );
$lf->new_ttk__frame( -padding => "3 3 12 12" );
$lf->g_grid(
    -column     => 0,
    -columnspan => 5,
    -rowspan    => 1,
    -row        => 0,
    -sticky     => "nwes"
);

my $ef = $lf->new_ttk__entry( -width => 40, -textvariable => \$fasta );
$ef->g_grid( -column => 1, -row => 0, -sticky => "we" );
my $em = $lf->new_ttk__button( -text => "Load File", -command => \&loadfile );
$em->g_grid( -column => 0, -row => 0, -sticky => "w" );
my $de = $lf->new_ttk__entry( -width => 40, -textvariable => \$seq_in );
$de->g_grid( -column => 4, -row => 0, -sticky => "w" );
my $bc = $lf->new_ttk__button(
    -text    => "Get Sequences",
    -command => sub { get_sequences(); }
);
$bc->g_grid( -column => 2, -row => 0, -sticky => "w" );

########## SNP TOP
my $elf = $frm->new_ttk__labelframe( -text => "SNP INPUT" );
$elf->new_ttk__frame( -padding => "3 3 12 12" );
$elf->g_grid(
    -column     => 0,
    -columnspan => 5,
    -rowspan    => 1,
    -row        => 14,
    -sticky     => "nwes",
    -pady       => 10
);

my $fef = $elf->new_ttk__entry( -width => 40, -textvariable => \$fasta_SNP );
$fef->g_grid( -column => 1, -row => 14, -sticky => "w" );
my $gem =
  $elf->new_ttk__button( -text => "Load SNP File", -command => \&loadfile_SNP );
$gem->g_grid( -column => 0, -row => 14, -sticky => "w" );
my $dej = $elf->new_ttk__entry( -width => 40, -textvariable => \$seq_in_SNP );
$dej->g_grid( -column => 4, -row => 14, -sticky => "w" );
my $gbc = $elf->new_ttk__button(
    -text    => "Get SNP Sequences",
    -command => sub { get_sequences_SNP(); }
);
$gbc->g_grid( -column => 2, -row => 14, -sticky => "w" );

my $djh = $frm->new_ttk__labelframe( -text => "PRIMER DESIGN" );
$djh->new_ttk__frame( -padding => "3 3 12 12" );
$djh->g_grid(
    -column     => 0,
    -columnspan => 2,
    -rowspan    => 1,
    -row        => 16,
    -sticky     => "nwes"
);

my $dgh = $djh->new_ttk__entry(
    -width        => 3,
    -textvariable => \$five_prime_SNP
);
$dgh->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $dhi = $djh->new_ttk__entry(
    -width        => 3,
    -textvariable => \$three_prime_SNP
);
$dhi->g_grid( -column => 5, -row => 2, -sticky => "we" );

my $djk = $djh->new_ttk__entry( -width => 3, -textvariable => \$kmer_max_SNP );
$djk->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $dkl = $djh->new_ttk__entry( -width => 3, -textvariable => \$kmer_min_SNP );
$dkl->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $dlm = $djh->new_ttk__entry( -width => 3, -textvariable => \$higher_gc_SNP );
$dlm->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $dmn = $djh->new_ttk__entry( -width => 3, -textvariable => \$lower_gc_SNP );
$dmn->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $dno = $djh->new_ttk__entry( -width => 3, -textvariable => \$upper_tm_SNP );
$dno->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $dop = $djh->new_ttk__entry( -width => 3, -textvariable => \$lower_tm_SNP );
$dop->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $dpq = $djh->new_ttk__entry( -width => 3, -textvariable => \$clamp_SNP );
$dpq->g_grid( -column => 2, -row => 6, -sticky => "e" );

my $dpz = $djh->new_ttk__entry( -width => 3, -textvariable => \$salt_SNP );
$dpz->g_grid( -column => 2, -row => 7, -sticky => "we" );

my $dqr = $djh->new_ttk__entry( -width => 3, -textvariable => \$spec_SNP );
$dqr->g_grid( -column => 5, -row => 6, -sticky => "w" );

my $drs = $djh->new_ttk__entry( -width => 3, -textvariable => \$mis_SNP );
$drs->g_grid( -column => 5, -row => 7, -sticky => "w" );

my $dpgz = $djh->new_ttk__entry( -width => 3, -textvariable => \$DNA_conc_SNP );
$dpgz->g_grid( -column => 7, -row => 2, -sticky => "we" );

my $dghy =
  $djh->new_ttk__entry( -width => 10, -textvariable => \$outputfile_SNP );
$dghy->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $dgjy =
  $djh->new_ttk__entry( -width => 10, -textvariable => \$selfie_cuttof_SNP );
$dgjy->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $fdgjy = $djh->new_ttk__entry( -width => 10, -textvariable => \$format );
$fdgjy->g_grid( -column => 7, -row => 6, -sticky => "we" );

my $dgjxy =
  $djh->new_ttk__entry( -width => 10, -textvariable => \$snp_distance );
$dgjxy->g_grid( -column => 7, -row => 5, -sticky => "we" );

my $rdgjxy =
  $djh->new_ttk__entry( -width => 10, -textvariable => \$min_size_seq );
$rdgjxy->g_grid( -column => 7, -row => 7, -sticky => "we" );

my $dabc = $djh->new_ttk__label( -text => "5' up from SNP" );
$dabc->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $dbca = $djh->new_ttk__label( -text => "3' down from SNP" );
$dbca->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $dcba = $djh->new_ttk__label( -text => "Primer length max" );
$dcba->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $dcab = $djh->new_ttk__label( -text => "Primer length min" );
$dcab->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $ddef = $djh->new_ttk__label( -text => "Upper GC%" );
$ddef->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $defg = $djh->new_ttk__label( -text => "Lower GC%" );
$defg->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $dhij = $djh->new_ttk__label( -text => "Upper Tm" );
$dhij->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $dijk = $djh->new_ttk__label( -text => "Lower Tm" );
$dijk->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $djkl = $djh->new_ttk__label( -text => "GC clamp (Y/N)" );
$djkl->g_grid( -column => 1, -row => 6, -sticky => "w" );

my $dklm = $djh->new_ttk__label( -text => "input specificity (Y/N)" );
$dklm->g_grid( -column => 4, -row => 6, -sticky => "w" );

my $dlmn = $djh->new_ttk__label( -text => "mis-matches" );
$dlmn->g_grid( -column => 4, -row => 7, -sticky => "w" );

my $dlzn = $djh->new_ttk__label( -text => "salt concentration (mM)" );
$dlzn->g_grid( -column => 1, -row => 7, -sticky => "w" );

my $dlnf = $djh->new_ttk__label( -text => "DNA concentration (nM)" );
$dlnf->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $dltf = $djh->new_ttk__label( -text => "SNP Primer file" );
$dltf->g_grid( -column => 6, -row => 3, -sticky => "w" );

my $dlqf = $djh->new_ttk__label( -text => "Self-complementarity" );
$dlqf->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $dlxqf = $djh->new_ttk__label( -text => "Minimum distance from SNP" );
$dlxqf->g_grid( -column => 6, -row => 5, -sticky => "w" );

my $adlxqf = $djh->new_ttk__label( -text => "Graphic format (.png or .gif)" );
$adlxqf->g_grid( -column => 6, -row => 6, -sticky => "w" );

my $ladlxqf = $djh->new_ttk__label( -text => "Minimum size sequence" );
$ladlxqf->g_grid( -column => 6, -row => 7, -sticky => "w" );

foreach ( Tkx::SplitList( $djh->g_winfo_children ) )
  {
    Tkx::grid_configure( $_, -padx => 5, -pady => 5 );
  }

my $dzs =
  $frm->new_ttk__labelframe( -text => "RUN", -width => 10, -height => 60 );
$dzs->new_ttk__frame( -padding => "3 3 12 12" );
$dzs->g_grid( -column => 0, -row => 28, -sticky => "nwes" );

my $dst = $dzs->new_ttk__button(
    -text    => "1: Design SNP Primers and Graphics",
    -command => sub { calculate_SNP(); graphics_all_primers_SNP(); graphics_single_view_SNP(); }
);
$dst->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $xduv = $dzs->new_ttk__button(
    -text    => "2: Multiplex PCR dimer scores",
    -command => sub { primer_dimer_SNP(); }
);
$xduv->g_grid( -column => 4, -row => 9, -sticky => "s" );

my $duv = $dzs->new_ttk__button(
    -text    => "3: Clean-up (run last)",
    -command => sub { clean_up(); }
);
$duv->g_grid( -column => 5, -row => 9, -sticky => "s" );

########## SNP BTM

my $doh = $frm->new_ttk__labelframe( -text => "PRIMER DESIGN" );
$doh->new_ttk__frame( -padding => "3 3 12 12" );
$doh->g_grid(
    -column     => 0,
    -columnspan => 2,
    -rowspan    => 1,
    -row        => 2,
    -sticky     => "nwes"
);

my $gh = $doh->new_ttk__entry(
    -width        => 3,
    -textvariable => \$five_prime_end
);
$gh->g_grid( -column => 2, -row => 2, -sticky => "we" );

my $hi = $doh->new_ttk__entry(
    -width        => 3,
    -textvariable => \$three_prime_end
);
$hi->g_grid( -column => 5, -row => 2, -sticky => "we" );

my $jk = $doh->new_ttk__entry( -width => 3, -textvariable => \$kmer_max );
$jk->g_grid( -column => 2, -row => 3, -sticky => "we" );

my $kl = $doh->new_ttk__entry( -width => 3, -textvariable => \$kmer_min );
$kl->g_grid( -column => 5, -row => 3, -sticky => "we" );

my $lm = $doh->new_ttk__entry( -width => 3, -textvariable => \$higher_gc );
$lm->g_grid( -column => 2, -row => 4, -sticky => "we" );

my $mn = $doh->new_ttk__entry( -width => 3, -textvariable => \$lower_gc );
$mn->g_grid( -column => 5, -row => 4, -sticky => "we" );

my $no = $doh->new_ttk__entry( -width => 3, -textvariable => \$upper_tm );
$no->g_grid( -column => 2, -row => 5, -sticky => "we" );

my $op = $doh->new_ttk__entry( -width => 3, -textvariable => \$lower_tm );
$op->g_grid( -column => 5, -row => 5, -sticky => "we" );

my $pq = $doh->new_ttk__entry( -width => 3, -textvariable => \$clamp );
$pq->g_grid( -column => 2, -row => 6, -sticky => "e" );

my $pz = $doh->new_ttk__entry( -width => 3, -textvariable => \$salt );
$pz->g_grid( -column => 2, -row => 7, -sticky => "we" );

my $qr = $doh->new_ttk__entry( -width => 3, -textvariable => \$spec );
$qr->g_grid( -column => 5, -row => 6, -sticky => "w" );

my $rs = $doh->new_ttk__entry( -width => 3, -textvariable => \$mis );
$rs->g_grid( -column => 5, -row => 7, -sticky => "w" );

my $pgz = $doh->new_ttk__entry( -width => 3, -textvariable => \$DNA_conc );
$pgz->g_grid( -column => 7, -row => 2, -sticky => "we" );

my $ghy = $doh->new_ttk__entry( -width => 10, -textvariable => \$outputfile );
$ghy->g_grid( -column => 7, -row => 3, -sticky => "we" );

my $gjy =
  $doh->new_ttk__entry( -width => 10, -textvariable => \$selfie_cuttof );
$gjy->g_grid( -column => 7, -row => 4, -sticky => "we" );

my $afdgjy = $doh->new_ttk__entry( -width => 10, -textvariable => \$format );
$afdgjy->g_grid( -column => 7, -row => 5, -sticky => "we" );

my $xafdgjy =
  $doh->new_ttk__entry( -width => 10, -textvariable => \$min_size_seq );
$xafdgjy->g_grid( -column => 7, -row => 6, -sticky => "we" );

my $htmlafdgjy =
  $doh->new_ttk__entry( -width => 10, -textvariable => \$html );
$htmlafdgjy->g_grid( -column => 7, -row => 7, -sticky => "we" );

my $zs =
  $frm->new_ttk__labelframe( -text => "RUN", -width => 10, -height => 40 );
$zs->new_ttk__frame( -padding => "3 3 12 12" );
$zs->g_grid( -column => 0, -row => 8, -sticky => "nwes" );

my $sp = $mw->new_ttk__labelframe(
    -text => "BLAST PRIMERS - insert primers here in fastA format" );
$sp->new_ttk__frame( -padding => "3 3 12 12" );
$sp->g_grid( -column => 0, -row => 10, -sticky => "nwes" );

my $text = $sp->new_tk__text( -width => 80, -height => 6, -wrap => "none" );
$text->g_grid;
my $thetext;

my $xxtuz = $mw->new_ttk__button(
    -text    => "EXIT",
    -command => sub { exit_program(); }
);
$xxtuz->g_grid( -column => 0, -row => 14, -sticky => "s" );

my $gfd = $mw->new_ttk__labelframe( -text => "SELECT DATABASE:" );
$gfd->g_grid( -column => 0, -row => 12, -sticky => "nwes" );

my $hme = $gfd->new_ttk__radiobutton(
    -text     => "Nucleotide collection (nt)",
    -variable => \$blaster,
    -value    => "nt"
);
my $ofz = $gfd->new_ttk__radiobutton(
    -text     => "Genomic Human",
    -variable => \$blaster,
    -value    => "human_genomic"
);
my $cez = $gfd->new_ttk__radiobutton(
    -text     => "Genomic others",
    -variable => \$blaster,
    -value    => "other_genomic"
);
my $ahme = $gfd->new_ttk__radiobutton(
    -text     => "EST others",
    -variable => \$blaster,
    -value    => "est_others"
);
my $aofz = $gfd->new_ttk__radiobutton(
    -text     => "EST Human",
    -variable => \$blaster,
    -value    => "est_human"
);
my $acez = $gfd->new_ttk__radiobutton(
    -text     => "EST Mouse",
    -variable => \$blaster,
    -value    => "est_mouse"
);

$hme->g_grid( -column => 0, -row => 13, -sticky => "w" );
$ofz->g_grid( -column => 1, -row => 13, -sticky => "w" );
$cez->g_grid( -column => 2, -row => 13, -sticky => "w" );
$ahme->g_grid( -column => 3, -row => 13, -sticky => "w" );
$aofz->g_grid( -column => 4, -row => 13, -sticky => "w" );
$acez->g_grid( -column => 5, -row => 13, -sticky => "w" );

my $tuz = $gfd->new_ttk__button(
    -text    => "BLAST",
    -command => sub { web_blast_ncbi(); }
);
$tuz->g_grid( -column => 0, -row => 14, -sticky => "w" );

my $st = $zs->new_ttk__button(
    -text    => "1: Design Primers and Graphics (run this first)",
    -command => sub { calculate(); graphics_all_primers(); graphics_single_view();}
);
$st->g_grid( -column => 1, -row => 9, -sticky => "e" );

my $ghuv = $zs->new_ttk__button(
    -text    => "2: Multiplex PCR dimer scores",
    -command => sub { primer_dimer(); }
);
$ghuv->g_grid( -column => 4, -row => 9, -sticky => "s" );

my $uv = $zs->new_ttk__button(
    -text    => "3: Clean-up (run last)",
    -command => sub { clean_up(); }
);
$uv->g_grid( -column => 5, -row => 9, -sticky => "s" );

my $abc = $doh->new_ttk__label( -text => "5' search area" );
$abc->g_grid( -column => 1, -row => 2, -sticky => "w" );

my $bca = $doh->new_ttk__label( -text => "3' search area" );
$bca->g_grid( -column => 4, -row => 2, -sticky => "w" );

my $cba = $doh->new_ttk__label( -text => "Primer length max" );
$cba->g_grid( -column => 1, -row => 3, -sticky => "w" );

my $cab = $doh->new_ttk__label( -text => "Primer length min" );
$cab->g_grid( -column => 4, -row => 3, -sticky => "w" );

my $def = $doh->new_ttk__label( -text => "Upper GC%" );
$def->g_grid( -column => 1, -row => 4, -sticky => "w" );

my $efg = $doh->new_ttk__label( -text => "Lower GC%" );
$efg->g_grid( -column => 4, -row => 4, -sticky => "w" );

my $hij = $doh->new_ttk__label( -text => "Upper Tm" );
$hij->g_grid( -column => 1, -row => 5, -sticky => "w" );

my $ijk = $doh->new_ttk__label( -text => "Lower Tm" );
$ijk->g_grid( -column => 4, -row => 5, -sticky => "w" );

my $jkl = $doh->new_ttk__label( -text => "GC clamp (Y/N)" );
$jkl->g_grid( -column => 1, -row => 6, -sticky => "w" );

my $klm = $doh->new_ttk__label( -text => "input specificity (Y/N)" );
$klm->g_grid( -column => 4, -row => 6, -sticky => "w" );

my $lmn = $doh->new_ttk__label( -text => "mis-matches" );
$lmn->g_grid( -column => 4, -row => 7, -sticky => "w" );

my $lzn = $doh->new_ttk__label( -text => "salt concentration (mM)" );
$lzn->g_grid( -column => 1, -row => 7, -sticky => "w" );

my $lnf = $doh->new_ttk__label( -text => "DNA concentration (nM)" );
$lnf->g_grid( -column => 6, -row => 2, -sticky => "w" );

my $ltf = $doh->new_ttk__label( -text => "Primer text file" );
$ltf->g_grid( -column => 6, -row => 3, -sticky => "w" );

my $lqf = $doh->new_ttk__label( -text => "Self-complementarity" );
$lqf->g_grid( -column => 6, -row => 4, -sticky => "w" );

my $alqf = $doh->new_ttk__label( -text => "Graphic format (.png or .gif)" );
$alqf->g_grid( -column => 6, -row => 5, -sticky => "w" );

my $galqf = $doh->new_ttk__label( -text => "Minimum size sequence" );
$galqf->g_grid( -column => 6, -row => 6, -sticky => "w" );

my $htmlgalqf = $doh->new_ttk__label( -text => "HTML5 Canvas (Y/N)" );
$htmlgalqf->g_grid( -column => 6, -row => 7, -sticky => "w" );

foreach ( Tkx::SplitList( $doh->g_winfo_children ) )
  {
    Tkx::grid_configure( $_, -padx => 5, -pady => 5 );
  }

sub exit_program
  {
    Tkx::tk___messageBox( -message => "Exiting PrimerMapper" );
    $mw->g_destroy;
  }

####Start subs
sub get_sequences
  {
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

sub get_sequences_SNP
  {
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

sub loadfile
  {
    $fasta = Tkx::tk___getOpenFile();

    $outfile = 'specificty.txt';
    open OUT, ">>$outfile";

    open IN, "<$fasta";
    while (<IN>)
      {
        print OUT if ( $_ !~ /(>.+)/ );
      }
    close IN;
    close OUT;

    #stringify the specificity file
    open OUT, $outfile or die "Couldn't open file: $!";
    $specificty = do { local $/; <OUT> };
    while (<OUT>)
      {
        $specificty .= $_;
      }
    close OUT;
    $specificty =~ s/\s//g;
    $specificty =~ s/[^agctn]//ig;

    print "\nSequences Loaded.\n";

  }

sub loadfile_SNP
  {
    $fasta_SNP = Tkx::tk___getOpenFile();

    $outfile_SNP = 'specificty_SNP.txt';
    open OUT_SNP, ">>$outfile_SNP ";

    open IN_SNP, "<$fasta_SNP";
    while (<IN_SNP>)
      {
        print OUT_SNP if ( $_ !~ />.+/ );
      }
    close IN_SNP;
    close OUT_SNP;

    #stringify the specificity file
    open OUT_SNP, $outfile_SNP or die "Couldn't open file: $!";
    $specificty_SNP = do { local $/; <OUT_SNP> };
    while (<OUT_SNP>)
      {
        $specificty_SNP .= $_;
      }
    close OUT_SNP;
    $specificty_SNP =~ s/\s//g;
    $specificty_SNP =~ s/[^agctn]//ig;

    print "\nSequences Loaded.\n";

  }

sub calculate
  {
    my $seqio = Bio::SeqIO->new(
        -file   => $fasta,
        -format => "fasta",
    );

    print "\nDesigning and validating primers....\n";

    while ( my $seqobj = $seqio->next_seq() )
      {
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
        $outputfile_html = "canvas_$id.txt";

    #skip to next if selected 5' or 3' length is longer than the sequence length
        if ( $five_prime_end > length $sequence )
          {
            next;
          }
        if ( $three_prime_end > length $sequence )
          {
            next;
          }

        #skip to next if sequence length is <100bps
        if ( length $sequence < $min_size_seq )
          {
            next;
          }

        open FUSIONFILE, ">>$outputfile";
        print FUSIONFILE "Header\tStart\tTm(degC)\tSequence\tSelf-comp\tGC%\n";
        close FUSIONFILE;
        
        
           if ($html eq "Y") {
              open CANVASFILE, ">>$outputfile_html";
               print CANVASFILE "FS\tLF\tRS\tLR\t>\t#\tSEQ\n\t\t\t\t$id\t$len_seq\t$sequence\n";
               close CANVASFILE;
           }

######## FORWARD PRIMER
        #start counting
        my $start = 1;
        for ( my $i = $start - 1 ; $i < $five_prime_end - 1 ; $i += 1 )
          {
            for ( my $a = $kmer_max ; $a >= $kmer_min ; $a-- )
              {
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
                if ( $a <= 36 )
                  {
                    $Tm = calcTm( $_, $DNA_conc, $salt );
                  }
                else
                  {
                    $Tm = calclongTm( $_, $DNA_conc, $salt );
                  }
                my $Tmrounded = sprintf( "%0.1f", $Tm );

                my $hairpin = calcdG($_);

                my $primer_end = $i + $a;

                #capture matches
                my $number_matches = 0;
                if ( $spec eq "Y" )
                  {
                    my $mis_mismatch = fuzzy_pattern( $_, $mis );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty );
                    $number_matches = @approximate_matches;
                  }
           

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                              if (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $spec eq "N" 
                    && $html eq "N")
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
                   

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "N" )
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
                   

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" 
                    && $html eq "N")
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
                  

                  }

                elsif (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y"
                    && $html eq "N" )
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
                   

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y" 
                    && $html eq "N")
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
                    
                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $number_matches < 2
                    && $spec eq "Y" 
                    && $html eq "N")
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
                   

                  }
            ############HTML
               
                elsif (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $spec eq "N" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                elsif (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/gc$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y"
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $_ =~ m/cg$/i
                    && $clamp eq "Y"
                    && $number_matches < 2
                    && $spec eq "Y" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tmrounded ge $lower_tm
                    && $Tmrounded le $upper_tm
                    && $percentGC ge $lower_gc
                    && $percentGC le $higher_gc
                    && $selfie_score < $selfie_cuttof
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $clamp eq "N"
                    && $number_matches < 2
                    && $spec eq "Y" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
                    print CANVASFILE "$i,\t$a,\n";
                    close(CANVASFILE);

                  }
              }
          }

######## REVERSE PRIMER

        #start counting for reverse primer
        for (
            my $j = length($sequence) - $three_prime_end ;
            $j < length($sequence) ;
            $j += 1
          )
          {
            for ( my $y = $kmer_max_SNP ; $y >= $kmer_min_SNP ; $y-- )
              {
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
                if ( $y <= 36 )
                  {
                    $Tm2 = calcTm( $_, $DNA_conc, $salt );
                  }
                else
                  {
                    $Tm2 = calclongTm( $_, $DNA_conc, $salt );
                  }
                my $Tm_rounded = sprintf( "%0.1f", $Tm2 );

                my $hairpin_r = calcdG($_);

                my $primer_start_R = $j + $y;

                #capture matches
                my $number_matches_R = 0;
                if ( $spec eq "Y" )
                  {
                    my $mis_mismatch = fuzzy_pattern( $_, $mis );
                    my @approximate_matches =
                      match_positions( $mis_mismatch, $specificty );
                    $number_matches_R = @approximate_matches;
                  }

                #define dinucleotide repeats and repetitive sequence
                #and print results if statements are matched
                if (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" 
                    && $html eq "N")
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
                     

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "N" )
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
                  

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "N" )
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
                   

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }

                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "N")
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
                    

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "N")
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
                    

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
                elsif (open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "N")
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
                   

                    print FUSIONFILE
"$id\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
                ############html
                elsif (   open( FUSIONFILE, ">>$outputfile" )
                    && $Tm_rounded ge $lower_tm
                    && $Tm_rounded le $upper_tm
                    && $percent_GC ge $lower_gc
                    && $percent_GC le $higher_gc
                    && $selfie_scoreR < $selfie_cuttof
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $spec eq "N" 
                    && $html eq "Y")
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
                     open( CANVASFILE, ">>$outputfile_html" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "Y" )
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $spec eq "N"
                    && $html eq "Y" )
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp eq "N"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "Y" )
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp eq "Y"
                    && $number_matches_R < 2
                    && $spec eq "Y" 
                    && $html eq "Y")
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
                    open( CANVASFILE, ">>$outputfile_html" ) or die;
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

sub calculate_SNP
  {

    $fasta_SNP_NEW = "processed_SNP.fasta";
    open INNER, '<', $fasta_SNP     or die "Can't read old file: $!";
    open OUTER, '>', $fasta_SNP_NEW or die "Can't write new file: $!";

    while (<INNER>)
      {
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

    while ( my $seqobj = $seqio->next_seq() )
      {
        my $sequence = $seqobj->seq();
        $id_SNP = $seqobj->id();

        if ( $id_SNP =~ m/pos=(.+).len/i )
          {
            $snp = $1;
          }

        ###################
        if ( $sequence =~ m/R/ )
          {
            $degen = "type_1";
          }
        elsif ( $sequence =~ m/Y/ )
          {
            $degen = "type_2";
          }
        elsif ( $sequence =~ m/S/ )
          {
            $degen = "type_3";
          }
        elsif ( $sequence =~ m/W/ )
          {
            $degen = "type_4";
          }
        elsif ( $sequence =~ m/K/ )
          {
            $degen = "type_5";
          }
        elsif ( $sequence =~ m/M/ )
          {
            $degen = "type_6";
          }
        elsif ( $sequence =~ m/B/ )
          {
            $degen = "type_7";
          }
        elsif ( $sequence =~ m/D/ )
          {
            $degen = "type_8";
          }
        elsif ( $sequence =~ m/H/ )
          {
            $degen = "type_9";
          }
        elsif ( $sequence =~ m/V/ )
          {
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
        if ( $five_prime_SNP > $snp )
          {
            next;
          }
        if ( $three_prime_SNP > ( $len_seq - $snp ) )
          {
            next;
          }

        #skip to next if sequence length is <100bps
        if ( length $sequence < $min_size_seq )
          {
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

######## FORWARD PRIMER
        #start counting
        my $upstream = $snp - $five_prime_SNP;

        for ( my $i = $upstream ; $i < $snp - $snp_distance ; $i += 1 )
          {
            for ( my $q = $kmer_max_SNP ; $q >= $kmer_min_SNP ; $q-- )
              {
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
                if ( $q <= 36 )
                  {
                    $Tm_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
                  }
                else
                  {
                    $Tm_SNP = calclongTm( $_, $DNA_conc_SNP, $salt_SNP );
                  }
                my $Tmrounded = sprintf( "%0.1f", $Tm_SNP );

                my $hairpin = calcdG($_);

                my $primer_end = $i + $q;

                #capture matches
                my $number_matches = 0;
                if ( $spec_SNP eq "Y" )
                  {
                    my $mis_mismatch = fuzzy_pattern( $_, $mis_SNP );
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
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
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
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                  }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
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
                    open( OLIGOS_SNP, ">>$out_image" ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                  }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
                    && $hairpin > "-9"
                    && $clamp_SNP eq "N"
                    && $spec_SNP eq "N" )
                  {
                    print FUSIONFILE_SNP
"$id_SNP\t$i\t$Tmrounded\tF:$_\t$selfie_score\t$percentGCrounded\n";
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$_";
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
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
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                  }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
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
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                  }
                elsif (open( FUSIONFILE_SNP, ">>$outputfile_SNP " )
                    && $Tmrounded ge $lower_tm_SNP
                    && $Tmrounded le $upper_tm_SNP
                    && $percentGC ge $lower_gc_SNP
                    && $percentGC le $higher_gc_SNP
                    && $selfie_score < $selfie_cuttof_SNP
                    && $_ !~ /AAAAA/i
                    && $_ !~ /TTTTT/i
                    && $_ !~ /GGGGG/i
                    && $_ !~ /CCCCC/i
                    && $_ !~ /ATATATAT/i
                    && $_ !~ /TATATATA/i
                    && $_ !~ /GCGCGCGC/i
                    && $_ !~ /CGCGCGCG/i
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
                    open( OLIGOS_SNP, ">>$out_image_SNP " ) or die;
                    print OLIGOS_SNP "$snp\t$selfie_score\t$i\t$primer_end\n";
                    close(OLIGOS_SNP);
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$i\t$foo_SNP \t$i\t$primer_end\n";
                    close(OUTS_SNP);

                  }
              }
          }

######## REVERSE PRIMER

        #start counting for reverse primer
        my $downstream = $snp + $three_prime_SNP;

        for ( my $j = $snp + $snp_distance ; $j < $downstream ; $j += 1 )
          {
            for ( my $c = $kmer_max_SNP ; $c >= $kmer_min_SNP ; $c-- )
              {
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
                if ( $c <= 36 )
                  {
                    $Tm2_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
                  }
                else
                  {
                    $Tm2_SNP = calclongTm( $_, $DNA_conc_SNP, $salt_SNP );
                  }
                my $Tm_rounded = sprintf( "%0.1f", $Tm2_SNP );

                my $hairpin_r = calcdG($_);

                my $primer_start_R = $j + $c;

                #capture matches
                my $number_matches_R = 0;
                if ( $spec_SNP eq "Y" )
                  {
                    my $mis_mismatch = fuzzy_pattern( $_, $mis_SNP );
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp_SNP eq "N"
                    && $spec_SNP eq "N" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp_SNP eq "Y"
                    && $spec_SNP eq "N" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $clamp_SNP eq "N"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^gc/i
                    && $clamp_SNP eq "Y"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
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
                    && $revR !~ /AAAAA/i
                    && $revR !~ /TTTTT/i
                    && $revR !~ /GGGGG/i
                    && $revR !~ /CCCCC/i
                    && $revR !~ /ATATATAT/i
                    && $revR !~ /TATATATA/i
                    && $revR !~ /GCGCGCGC/i
                    && $revR !~ /CGCGCGCG/i
                    && $hairpin_r > "-9"
                    && $_ =~ m/^cg/i
                    && $clamp_SNP eq "Y"
                    && $number_matches_R < 2
                    && $spec_SNP eq "Y" )
                  {
                    open( OLIGOS_SNP, ">>$out_image_SNP" ) or die;
                    print OLIGOS_SNP
                      "$snp\t$selfie_scoreR\t$primer_start_R\t$j\n";
                    close(OLIGOS_SNP);
                    push @array_length_SNP, $len_uniq_SNP;
                    push @array_name_SNP,   $id_uniq_SNP;
                    push @SNP_primers,      "$id_SNP:$revR";
                    open( OUTS_SNP, ">>$out_single_SNP" ) or die;
                    print OUTS_SNP "$j\t$foo_SNP \t$primer_start_R\t$j\n";
                    close(OUTS_SNP);

                    print FUSIONFILE_SNP
"$id_SNP\t$j\t$Tm_rounded\tR:$revR\t$selfie_scoreR\t$percentGC_rounded\n";

                  }
              }
          }

        ####################______Allele-Specific_BEGIN

        for ( my $p = $kmer_min_SNP ; $p <= $kmer_max_SNP ; $p++ )
          {
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
            if ( $p <= 36 )
              {
                $Tm3_SNP = calcTm( $revFA, $DNA_conc_SNP, $salt_SNP );
              }
            else
              {
                $Tm3_SNP = calclongTm( $revFA, $DNA_conc_SNP, $salt_SNP );
              }
            my $TmroundedA = sprintf( "%0.1f", $Tm3_SNP );

            my $hairpinA = calcdG($revFA);

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if (   open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.G\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.C\t$revFA.G\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.G\t$revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.C\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.A\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_10" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedA \tReverse: $revFA.T\t$revFA.G\t$revFA.C\t$selfie_scoreA\t$hairpinA\t$percentGCroundedA\n";

              }

          }

###################

        my $string_to_snp = $snp - 2;
        my $new_seq = substr( $sequence, 1, $string_to_snp );

        for ( my $u = $kmer_max_SNP ; $u >= $kmer_min_SNP ; $u-- )
          {
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
            if ( $u <= 36 )
              {
                $Tm4_SNP = calcTm( $_, $DNA_conc_SNP, $salt_SNP );
              }
            else
              {
                $Tm4_SNP = calclongTm( $_, $DNA_conc_SNP, $salt_SNP );
              }
            my $TmroundedAB = sprintf( "%0.1f", $Tm4_SNP );

            my $hairpinAB = calcdG($_);

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if (   open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_1" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.G\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_2" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.C\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_3" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.G\t$_.C\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_4" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_5" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_6" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.C\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_7" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.C\t$_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_8" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tReverse:$_.A\t$_.G\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
                && $degen eq "type_9" )
              {
                print AS_SNP
"$id_SNP\t$TmroundedAB \tForward: $_.A\t$_.C\t$_.T\t$selfie_scoreAB\t$hairpinAB\t$percentGCroundedAB\n";

              }

            elsif (open( AS_SNP, ">>$outputfile_SNP_AS" )
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
####SUBROUTINES####
####################################
use re qw(eval);
use vars qw($matchStart);

sub match_positions
  {
    my $pattern;
    local $_;
    ( $pattern, $_ ) = @_;
    my @results;
    local $matchStart;
    my $instrumentedPattern = qr/(?{ $matchStart = pos() })$pattern/;
    while (/$instrumentedPattern/g)
      {
        my $nextStart = pos();
        push @results, "[$matchStart..$nextStart)";
        pos() = $matchStart + 1;
      }
    return @results;
  }

####################################
sub fuzzy_pattern
  {
    my ( $original_pattern, $mismatches_allowed ) = @_;
    $mismatches_allowed >= 0
      or die "Number of mismatches must be greater than or equal to zero\n";
    my $new_pattern =
      make_approximate( $original_pattern, $mismatches_allowed );
    return qr/$new_pattern/;
  }

####################################
sub make_approximate
  {
    my ( $pattern, $mismatches_allowed ) = @_;
    if ( $mismatches_allowed == 0 ) { return $pattern }
    elsif ( length($pattern) <= $mismatches_allowed )
      {
        $pattern =~ tr/ACTG/./;
        return $pattern;
      }
    else
      {
        my ( $first, $rest ) = $pattern =~ /^(.)(.*)/;
        my $after_match = make_approximate( $rest, $mismatches_allowed );
        if ( $first =~ /[ACGT]/ )
          {
            my $after_miss = make_approximate( $rest, $mismatches_allowed - 1 );
            return "(?:$first$after_match|.$after_miss)";
          }
        else { return "$first$after_match" }
      }
  }

####################################
sub calclongTm
  {
    my $sequence = shift;
    my $DNA_nM   = shift;
    my $K_mM     = shift;
    return 81.5 + ( 16.6 * ( log( $K_mM / 1000.0 ) / log(10) ) ) +
      ( 41.0 * percentGC($sequence) ) - ( 600.0 / length($sequence) );
  }

####################################
sub calcdG
  {
    my $sequence = shift;
    my $sum;
    my @count;
    my $normalized_round;
    $sequence =~ tr/a-z/A-Z/;
    $sequence =~ s/\s//g;
    my $oligo_rev = reverse $sequence;
    $oligo_rev =~ tr/AGCT/TCGA/;
    my $seq_len = length($sequence);

    for my $i ( 0 .. $seq_len - 2 )
      {

        my $substr = substr $sequence, $i, 2;
        my $regex = $substr;
        if ( $oligo_rev =~ m/($regex)/ && $1 eq "AA" )
          {
            push( @count, -1 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "AT" )
          {
            push( @count, -0.88 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "TA" )
          {
            push( @count, -0.58 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "CA" )
          {
            push( @count, -1.45 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "GT" )
          {
            push( @count, -1.44 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "CT" )
          {
            push( @count, -1.28 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "GA" )
          {
            push( @count, -1.3 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "CG" )
          {
            push( @count, -2.17 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "GC" )
          {
            push( @count, -1.84 );
          }
        elsif ( $oligo_rev =~ m/($regex)/ && $1 eq "GG" )
          {
            push( @count, -1.42 );
          }

      }

    $sum = 0;
    for (@count)
      {
        $sum += $_;
      }

    my $begin        = 1;
    my $count_dimers = 0;

    for ( my $j = $begin - 1 ; $j < $seq_len - 2 ; $j += 1 )
      {

        my $checkF = substr( $sequence, $j, 2 );

        if ( $oligo_rev =~ m/$checkF/i )
          {

            $count_dimers = $count_dimers + 1;
          }
      }

    my $len_oligo_rev = length $oligo_rev;

    if ( $count_dimers != 0 )
      {
        my $multiplier       = $count_dimers / $len_oligo_rev;
        my $normalized_score = $multiplier * $sum;
        $normalized_round = sprintf( "%0.4f", $normalized_score );
      }
    else
      {
        $normalized_round = 0;
      }

    return $normalized_round;

  }

####################################
sub calcTm
  {
    my $sequence = shift;
    my $DNA_nM   = shift;
    my $K_mM     = shift;
    my $dH       = 0;
    my $dS       = 108;
    my $seq_len  = length($sequence);
    my $i;

    # Compute dH and dS
    for ( $i = 0 ; $i < ( $seq_len - 1 ) ; $i++ )
      {
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
sub selfie
  {
    my $primer_f  = shift;
    my $primer_rc = shift;
    my $FLEN      = length $primer_rc;
    my $RC_LEN    = length $primer_f;
    my $D         = [];
    for ( my $t = 0 ; $t <= $FLEN ; ++$t )
      {
        $D->[$t][0] = 0;
      }
    for ( my $p = 0 ; $p <= $RC_LEN ; ++$p )
      {
        $D->[0][$p] = $p;
      }
    for ( my $t = 1 ; $t <= $FLEN ; ++$t )
      {
        for ( my $p = 1 ; $p <= $RC_LEN ; ++$p )
          {
            $D->[$t][$p] =

              min3(
                substr( $primer_rc, $t - 1, 1 ) eq
                  substr( $primer_f, $p - 1, 1 )
                ? $D->[ $t - 1 ][ $p - 1 ]
                : $D->[ $t - 1 ][ $p - 1 ] + 1,

                $D->[ $t - 1 ][$p] + 1,

                $D->[$t][ $p - 1 ] + 1
              );
          }
      }

    my @matches   = ();
    my $bestscore = 10000000;

    for ( my $t = 1 ; $t <= $FLEN ; ++$t )
      {
        if ( $D->[$t][$RC_LEN] < $bestscore )
          {
            $bestscore = $D->[$t][$RC_LEN];
            @matches   = ($t);
          }
        elsif ( $D->[$t][$RC_LEN] == $bestscore )
          {
            push( @matches, $t );
          }
      }

    return $bestscore;
  }

####################################
sub min3
  {
    my ( $i, $j, $k ) = @_;
    my ($tmp);

    $tmp = ( $i < $j ? $i : $j );
    $tmp < $k ? $tmp : $k;
  }

####################################
sub graphics_all_primers
  {
    print "\nDesigning Grahical output files....\n";
    my $dir = cwd();
    my $orientation;
    my @gene_id;
    my @unique_length = uniq(@array_length);
    my @unique_name   = uniq(@array_name);
    my %hash;
    my $len;
    my $name_id;

    foreach my $unique_name (@unique_name)
      {
        if ( $unique_name =~ m/^\d+'(.+)?/i )
          {
            push @gene_id, $1;
          }
      }

    foreach my $unique_length (@unique_length)
      {
        if ( $unique_length =~ m/^.+'(\d+)?/i )
          {
            push @gene_length, $1;
          }
      }

    foreach my $fp ( glob("$dir/GRAPHIC_*.txt") )
      {
        open my $fh, "<", $fp or die;

        @hash{@gene_id} = @gene_length;
        while ( my ( $key, $val ) = each %hash )
          {
            $len = $val if $fp =~ m/$key/i;
          }
        while ( my ( $key, $val ) = each %hash )
          {
            $name_id = $key if $fp =~ m/$key/i;
          }

        my $outputfile = "$fp" . "$format";
        open( OUTFILE, ">$outputfile" ) or die;

        my $panel = Bio::Graphics::Panel->new(
            -length     => $len + 100,
            -width      => 800,
            -pad_left   => 30,
            -pad_right  => 30,
            -pad_top    => 20,
            -pad_bottom => 20,
        );

        my $full_length = Bio::SeqFeature::Generic->new(
            -start        => 1,
            -end          => $len + 100,
            -display_name => $name_id,
        );

        $panel->add_track(
            $full_length,
            -glyph        => 'arrow',
            -tick         => 2,
            -fgcolor      => 'black',
            -double       => 1,
            -label        => 1,
            -strand_arrow => 1,
        );

        my $track = $panel->add_track(
            -glyph        => 'transcript2',
            -label        => 1,
            -strand_arrow => 1,
            -bgcolor      => 'blue',
            -bump         => +1,
            -height       => 12,
        );

        while (<$fh>)
          {
            chomp;

            my ( $name, $score, $start, $end ) = split /\t/;
            if ( $start - $end < 1 )
              {
                $orientation = +1;
              }
            elsif ( $start - $end > 1 )
              {
                $orientation = -1;
              }
            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $start,
                -score        => $score,
                -start        => $start,
                -strand       => $orientation,
                -end          => $end
            );
            $track->add_feature($feature);

          }

        binmode OUTFILE;

        print OUTFILE $panel->png;
        close $fh;
      }

    close OUTFILE;

    print "\nGraphic Files Generated.\n";
  }

####################################
sub graphics_single_view
  {
    print "\nDesigning Single View Grahical output....\n";
    my $dir = cwd();
    my $orientation;
    my $lener     = eval join '+', @gene_length;
    my $first_seq = $gene_length[0];
    my $max       = max(@gene_length);

    my $outputfile = "single_view" . "$format";
    open( OUTFILE, ">>$outputfile" ) or die;

    open OUTS, "<", $out_single or die;

    my $panel = Bio::Graphics::Panel->new(
        -length     => $lener + $max,
        -width      => 4800,
        -pad_left   => 30,
        -pad_right  => 30,
        -pad_top    => 20,
        -pad_bottom => 20,
    );

    my $full_length = Bio::SeqFeature::Generic->new(
        -start        => 1,
        -end          => $lener + $max,
        -display_name => "single_view",
    );

    $panel->add_track(
        $full_length,
        -glyph        => 'arrow',
        -tick         => 2,
        -fgcolor      => 'black',
        -double       => 1,
        -label        => 1,
        -strand_arrow => 1,
    );

    my $track = $panel->add_track(
        -glyph        => 'transcript2',
        -label        => 1,
        -strand_arrow => 1,
        -bgcolor      => 'blue',
        -bump         => +1,
        -height       => 12,
    );

    while (<OUTS>)
      {
        chomp;

        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start - $end < 1 )
          {
            $orientation = +1;
          }
        elsif ( $start - $end > 1 )
          {
            $orientation = -1;
          }

        my $adder = $score - $first_seq;

        my $feature = Bio::SeqFeature::Generic->new(
            -display_name => $start,
            -score        => $score,
            -start        => $start + $adder,
            -strand       => $orientation,
            -end          => $end + $adder,
        );
        $track->add_feature($feature);

      }

    binmode OUTFILE;

    print OUTFILE $panel->png;

    close OUTS;
    close OUTFILE;

    print "\nSingle View Graphic File Generated.\n";
  }

####################################
sub graphics_all_primers_SNP
  {
    print "\nDesigning SNP Primer Grahical output files....\n";
    my $dir = cwd();
    my $orientation;
    my @gene_id;
    my @unique_length = uniq(@array_length_SNP);
    my @unique_name   = uniq(@array_name_SNP);
    my %hash;
    my $len;
    my $name_id;
    my ( $score, $start, $end );
    my $feature2;

    foreach my $unique_name (@unique_name)
      {
        if ( $unique_name =~ m/^\d+'(.+)?/i )
          {
            push @gene_id, $1;
          }
      }

    foreach my $unique_length (@unique_length)
      {
        if ( $unique_length =~ m/^.+'(\d+)?/i )
          {
            push @gene_length_SNP, $1;
          }
      }

    foreach my $fp ( glob("$dir/GRAPHIC_*.txt") )
      {
        open my $fh, "<", $fp or die;

        @hash{@gene_id} = @gene_length_SNP;
        while ( my ( $key, $val ) = each %hash )
          {
            $len = $val if $fp =~ m/$key/i;
          }
        while ( my ( $key, $val ) = each %hash )
          {
            $name_id = $key if $fp =~ m/$key/i;
          }

        my $outputfile_SNP = "$fp" . "$format";
        open( OUTFILE, ">$outputfile_SNP" ) or die;

        my $panel = Bio::Graphics::Panel->new(
            -length     => $len + 100,
            -width      => 800,
            -pad_left   => 30,
            -pad_right  => 30,
            -pad_top    => 20,
            -pad_bottom => 20,
        );

        my $full_length = Bio::SeqFeature::Generic->new(
            -start        => 1,
            -end          => $len + 100,
            -display_name => $name_id,
        );

        $panel->add_track(
            $full_length,
            -glyph        => 'arrow',
            -tick         => 2,
            -fgcolor      => 'black',
            -double       => 1,
            -label        => 1,
            -strand_arrow => 1,
        );

        my $track = $panel->add_track(
            -glyph        => 'transcript2',
            -label        => 1,
            -strand_arrow => 1,
            -bgcolor      => 'blue',
            -bump         => +1,
            -height       => 12,
        );
        my $track2 = $panel->add_track( -glyph => 'diamond', );

        while (<$fh>)
          {
            chomp;

            ( $snp, $score, $start, $end ) = split /\t/;
            if ( $start - $end < 1 )
              {
                $orientation = +1;
              }
            elsif ( $start - $end > 1 )
              {
                $orientation = -1;
              }
            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $start,
                -score        => $score,
                -start        => $start,
                -strand       => $orientation,
                -end          => $end
            );
            $feature2 = Bio::SeqFeature::Generic->new(
                -display_name => "SNP",
                -score        => $score,
                -start        => $snp,
                -end          => $snp
            );
            $track->add_feature($feature);

          }

        $track2->add_feature($feature2);

        binmode OUTFILE;

        print OUTFILE $panel->png;
        close $fh;
      }
    close OUTFILE;

    print "\nSNP Primer Graphic Files Generated.\n";
  }

####################################
sub graphics_single_view_SNP
  {
    print "\nDesigning Single View SNP Primer Grahical output file....\n";
    my $dir = cwd();
    my $orientation;
    my $lener     = eval join '+', @gene_length_SNP;
    my $first_seq = $gene_length_SNP[0];
    my $max       = max(@gene_length_SNP);

    my $outputfile_SNP = "single_view_SNP" . "$format";
    open( OUTFILE_SNP, ">>$outputfile_SNP" ) or die;

    open OUTS_SNP, "<", $out_single_SNP or die;

    my $panel = Bio::Graphics::Panel->new(
        -length     => $lener + $max,
        -width      => 4800,
        -pad_left   => 30,
        -pad_right  => 30,
        -pad_top    => 20,
        -pad_bottom => 20,
    );

    my $full_length = Bio::SeqFeature::Generic->new(
        -start        => 1,
        -end          => $lener + $max,
        -display_name => "single_view_SNP",
    );

    $panel->add_track(
        $full_length,
        -glyph        => 'arrow',
        -tick         => 2,
        -fgcolor      => 'black',
        -double       => 1,
        -label        => 1,
        -strand_arrow => 1,
    );

    my $track = $panel->add_track(
        -glyph        => 'transcript2',
        -label        => 1,
        -strand_arrow => 1,
        -bgcolor      => 'blue',
        -bump         => +1,
        -height       => 12,
    );

    while (<OUTS_SNP>)
      {
        chomp;

        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start - $end < 1 )
          {
            $orientation = +1;
          }
        elsif ( $start - $end > 1 )
          {
            $orientation = -1;
          }

        my $adder = $score - $first_seq;

        my $feature = Bio::SeqFeature::Generic->new(
            -display_name => $start,
            -score        => $score,
            -start        => $start + $adder,
            -strand       => $orientation,
            -end          => $end + $adder,
        );
        $track->add_feature($feature);

      }

    binmode OUTFILE_SNP;

    print OUTFILE_SNP $panel->png;

    close OUTS_SNP;
    close OUTFILE_SNP;

    print "\nSingle View SNP Primer Graphic File Generated.\n";
  }

####################################
sub primer_dimer
  {
    print
      "\nCalculating global primer dimer score file for multiplex PCR....\n";

    my $primer_dimer_input = "primers_combos.txt";
    my $primer_dimer;

    open IN, ">$primer_dimer_input" or die "Couldn't open file: $!";

    my $strings = \@primers;

    sub combine;

    print IN "@$_\n" for combine $strings, 2;
    close IN;

#######--GET_Primer_Dimer_Scores--#############
    my $primer_dimer_score = "primer_dimer_scores.tsv";
    open IN,    "<", "$primer_dimer_input" or die "Couldn't open file: $!";
    open OUTER, '>', $primer_dimer_score   or die "Can't write new file: $!";
    print OUTER
"BASED ON IN SILICO AND PCR TESTING, PRIMER DIMER SCORES SHOULD BE LESS THAN OR EQUAL TO 12 \n\n";
    while (<IN>)
      {

        if ( $_ =~ m/(.+):(.+)\s(.+):(.+)\n/i )
          {
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
sub primer_dimer_SNP
  {

    print
"\nCalculating global SNP primer dimer score file for multiplex PCR....\n";

    my $primer_dimer_input = "primers_combos.txt";
    my $primer_dimer;

    open IN, ">$primer_dimer_input" or die "Couldn't open file: $!";

    my $strings = \@SNP_primers;

    sub combine;

    print IN "@$_\n" for combine $strings, 2;
    close IN;

#######--GET_Primer_Dimer_Scores--#############
    my $primer_dimer_score = "primer_dimer_scores_SNP.tsv";
    open IN,    "<", "$primer_dimer_input" or die "Couldn't open file: $!";
    open OUTER, '>', $primer_dimer_score   or die "Can't write new file: $!";
    print OUTER
"BASED ON IN SILICO AND PCR TESTING, PRIMER DIMER SCORES SHOULD BE LESS THAN OR EQUAL TO 12 \n\n";
    while (<IN>)
      {

        if ( $_ =~ m/(.+):(.+)\s(.+):(.+)\n/i )
          {
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
####################################
sub combine
  {

    my ( $list, $n ) = @_;
    die "Insufficient list members for primer dimer analysis" if $n > @$list;

    return map [$_], @$list if $n <= 1;

    my @comb;

    for ( my $i = 0 ; $i + $n <= @$list ; ++$i )
      {
        my $val  = $list->[$i];
        my @rest = @$list[ $i + 1 .. $#$list ];
        push @comb, [ $val, @$_ ] for combine \@rest, $n - 1;
      }

    return @comb;
  }

####################################
sub uniq
  {
    my %seen;
    return grep { !$seen{$_}++ } @_;
  }

####################################
sub clean_up
  {
    my $dir = cwd();
    unlink glob "$dir/*.txt";

    print "\nDirectory has been cleaned.\n";
  }

####################################
=pod 
=WEB-BLAST
Note: must have web_blast.pl script in PATH for BLAST feature 
=cut
####################################
sub web_blast_ncbi
  {

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
=pod
 
=head1 DEFAULT INPUT PARAMETERS
 
PrimerMapper

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
   
=cut

