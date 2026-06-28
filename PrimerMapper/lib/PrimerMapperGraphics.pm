#
# Graphics module for PrimerMapper
#
# Please direct questions and support issues to <https://github.com/dohalloran/PrimerMapper/issues>
#
# Author: Damien O'Halloran, The George Washington University, 2015
#
# GNU GENERAL PUBLIC LICENSE
#
# POD documentation before the code

=head1 NAME
PrimerMapperGraphics - runs the graphical design side of things

=head2 SYNOPSIS
  use Cwd;
  use List::Util 'max';
  use Bio::Graphics;
  use Bio::SeqFeature::Generic;
  use base 'Exporter';
  
  for Usage, run by creating a PrimerMapper object within PrimerMapper_driver.pl
  which will call subroutines from PrimerMapperGraphics
  my $tmp = PrimerMapper->new();
  
  Default output as png 
  See $format to change output
  the default will be overwritten if changed

=head3 DESCRIPTION
The output from this package produces graphical primer::sequence
maps of primer distribution across each sequence and the entire input file

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
The rest of the documentation details each of the object methods.

=cut


package PrimerMapperGraphics;

use strict;
use warnings;
use Cwd;
use List::Util 'max';
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use base 'Exporter';

our @EXPORT =
  qw/ graphics_all_primers graphics_single_view graphics_all_primers_SNP graphics_single_view_SNP uniq /;

####################################

=head1 graphics_all_primers
 Title   :  graphics_all_primers
 Usage   :  graphics_all_primers( \@array_length, \@array_name, \@gene_length,
            $format );
 Function:  generates primer sequence maps for each inputted sequence 
 Returns :  GRAPHIC_<geneName>.png
=cut

sub graphics_all_primers {
    my $array_length = shift;
    my $array_name   = shift;
    my $gene_length  = shift;
    my $format       = shift;
    print "\nDesigning Graphical output files....\n";
    my $dir = cwd();
    my $orientation;
    my @gene_id;
    my @unique_length = uniq(@$array_length);
    my @unique_name   = uniq(@$array_name);
    my %hash;
    my $len;
    my $name_id;
    my @graphic_files = glob("$dir/GRAPHIC_*.txt");

    unless ( @unique_length && @unique_name && @graphic_files ) {
        warn "No primer graphics input was generated; skipping per-sequence graphics.\n";
        return 0;
    }

    foreach my $unique_name (@unique_name) {
        if ( $unique_name =~ m/^\d+'(.+)?/i ) {
            push @gene_id, $1;
        }
    }

    foreach my $unique_length (@unique_length) {
        if ( $unique_length =~ m/^.+'(\d+)?/i ) {
            push @$gene_length, $1;
        }
    }

    # BUG FIX: OUTFILE was a package-level bareword filehandle that was only
    # closed once outside the loop, meaning all but the last graphic file were
    # left open and could be corrupted. Now use a lexical $outfh per iteration.
    foreach my $fp (@graphic_files) {
        unless ( -s $fp ) {
            warn "$fp is empty; skipping.\n";
            next;
        }

        open my $fh, "<", $fp or do {
            warn "Cannot open $fp: $!\n";
            next;
        };

        @hash{@gene_id} = @$gene_length;
        while ( my ( $key, $val ) = each %hash ) {
            $len     = $val if $fp =~ m/$key/i;
            $name_id = $key if $fp =~ m/$key/i;
        }

        # Guard: skip if we couldn't resolve the sequence length
        unless ( defined $len && $len > 0 ) {
            warn "Could not determine length for $fp; skipping.\n";
            close $fh;
            next;
        }

        my $outputfile = "$fp$format";
        open my $outfh, ">", $outputfile or do {    # BUG FIX: was append >>
            warn "Cannot write $outputfile: $!\n";
            close $fh;
            next;
        };

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

        while (<$fh>) {
            chomp;
            my ( $name, $score, $start, $end ) = split /\t/;

            # BUG FIX: original used < 1 / > 1 which misses the == 0 case
            # (primer of length 1) and the == 2 case incorrectly.
            # Correct logic: left primer has start < end; right has start > end.
            if ( $start <= $end ) {
                $orientation = +1;
            }
            else {
                $orientation = -1;
            }

            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $start,
                -score        => $score,
                -start        => $start,
                -strand       => $orientation,
                -end          => $end,
            );
            $track->add_feature($feature);
        }

        binmode $outfh;
        print $outfh $panel->png;
        close $fh;
        close $outfh;
    }

    print "\nGraphic Files Generated.\n";
    return 1;
}

####################################

=head1 graphics_single_view
 Title   :  graphics_single_view
 Usage   :  graphics_single_view( \@gene_length, $format, $out_single );
 Function:  generates a single primer sequence map for all inputted sequences with their derived primers
 Returns :  single_view.png
=cut

sub graphics_single_view {
    my $gene_length = shift;
    my $format      = shift;
    my $out_single  = shift;
    print "\nDesigning Single View Graphical output....\n";
    my $dir = cwd();
    my $orientation;

    unless ( @$gene_length && -s $out_single ) {
        warn "No single-view primer graphics input was generated; skipping single view.\n";
        return 0;
    }

    my $lener     = eval join '+', @$gene_length;
    my $first_seq = $gene_length->[0];    # BUG FIX: was @$arr[0] (deprecated)
    my $max       = max(@$gene_length);

    my $outputfile = "single_view$format";

    # BUG FIX: was append (>>) which stacked multiple runs into one file,
    # making the graphic grow and eventually corrupt. Use > (overwrite).
    open my $outfh, ">", $outputfile or do {
        warn "Cannot write $outputfile: $!\n";
        return 0;
    };

    open my $outs_fh, "<", $out_single
      or do {
        warn "Cannot read $out_single: $!\n";
        close $outfh;
        return 0;
      };

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

    while (<$outs_fh>) {
        chomp;
        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start <= $end ) {    # BUG FIX: same orientation logic fix
            $orientation = +1;
        }
        else {
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

    binmode $outfh;
    print $outfh $panel->png;
    close $outs_fh;
    close $outfh;

    print "\nSingle View Graphic File Generated.\n";
    return 1;
}

####################################

=head1 graphics_all_primers_SNP
 Title   :  graphics_all_primers_SNP
 Usage   :  graphics_all_primers_SNP( \@array_length_SNP, \@array_name_SNP,
            \@gene_length_SNP, $format, $snp );
 Function:  generates a SNP primer sequence map for each inputted sequence 
 Returns :  GRAPHIC_<geneName>.png
=cut

sub graphics_all_primers_SNP {
    my $array_length_SNP = shift;
    my $array_name_SNP   = shift;
    my $gene_length_SNP  = shift;
    my $format           = shift;
    my $snp              = shift;
    print "\nDesigning SNP Primer Graphical output files....\n";
    my $dir = cwd();
    my $orientation;
    my @gene_id;
    my @unique_length = uniq(@$array_length_SNP);
    my @unique_name   = uniq(@$array_name_SNP);
    my %hash;
    my $len;
    my $name_id;
    my ( $score, $start, $end );
    my $feature2;
    my @graphic_files = glob("$dir/GRAPHIC_*.txt");

    unless ( @unique_length && @unique_name && @graphic_files ) {
        warn "No SNP primer graphics input was generated; skipping per-sequence SNP graphics.\n";
        return 0;
    }

    foreach my $unique_name (@unique_name) {
        if ( $unique_name =~ m/^\d+'(.+)?/i ) {
            push @gene_id, $1;
        }
    }

    foreach my $unique_length (@unique_length) {
        if ( $unique_length =~ m/^.+'(\d+)?/i ) {
            push @$gene_length_SNP, $1;
        }
    }

    foreach my $fp (@graphic_files) {
        unless ( -s $fp ) {
            warn "$fp is empty; skipping.\n";
            next;
        }

        open my $fh, "<", $fp or do {
            warn "Cannot open $fp: $!\n";
            next;
        };

        @hash{@gene_id} = @$gene_length_SNP;
        while ( my ( $key, $val ) = each %hash ) {
            $len     = $val if $fp =~ m/$key/i;
            $name_id = $key if $fp =~ m/$key/i;
        }

        unless ( defined $len && $len > 0 ) {
            warn "Could not determine length for $fp; skipping.\n";
            close $fh;
            next;
        }

        my $outputfile_SNP = "$fp$format";

        # BUG FIX: was bareword OUTFILE shared with non-SNP routine and
        # appended (>>) rather than overwritten. Now lexical + overwrite.
        open my $outfh, ">", $outputfile_SNP or do {
            warn "Cannot write $outputfile_SNP: $!\n";
            close $fh;
            next;
        };

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
        my $track2 = $panel->add_track( -glyph => 'diamond' );

        # BUG FIX: $feature2 was only added after the loop, meaning if the
        # file had zero primer lines $feature2 would be undef and
        # add_feature() would crash. Guard with defined check.
        $feature2 = undef;

        while (<$fh>) {
            chomp;
            ( $snp, $score, $start, $end ) = split /\t/;

            if ( $start <= $end ) {    # BUG FIX: same orientation logic
                $orientation = +1;
            }
            else {
                $orientation = -1;
            }

            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $start,
                -score        => $score,
                -start        => $start,
                -strand       => $orientation,
                -end          => $end,
            );

            # Recreate SNP diamond feature each line so it stays current
            $feature2 = Bio::SeqFeature::Generic->new(
                -display_name => "SNP",
                -score        => $score,
                -start        => $snp,
                -end          => $snp,
            );

            $track->add_feature($feature);
        }

        # Only add the SNP marker if we actually parsed at least one primer
        $track2->add_feature($feature2) if defined $feature2;

        binmode $outfh;
        print $outfh $panel->png;
        close $fh;
        close $outfh;
    }

    print "\nSNP Primer Graphic Files Generated.\n";
    return 1;
}

####################################

=head1 graphics_single_view_SNP
 Title   :  graphics_single_view_SNP
 Usage   :  graphics_single_view_SNP( \@gene_length_SNP, $format, $out_single_SNP );
 Function:  generates a single concatenated view of the SNP primer sequence maps for all inputted sequences
 Returns :  single_view_SNP.png
=cut

sub graphics_single_view_SNP {
    my $gene_length_SNP = shift;
    my $format          = shift;
    my $out_single_SNP  = shift;
    print "\nDesigning Single View SNP Primer Graphical output file....\n";
    my $dir = cwd();
    my $orientation;

    unless ( @$gene_length_SNP && -s $out_single_SNP ) {
        warn "No single-view SNP graphics input was generated; skipping SNP single view.\n";
        return 0;
    }

    my $lener     = eval join '+', @$gene_length_SNP;
    my $first_seq = $gene_length_SNP->[0];    # BUG FIX: deprecated array slice
    my $max       = max(@$gene_length_SNP);

    my $outputfile_SNP = "single_view_SNP$format";

    # BUG FIX: was append >> (same stacking problem as graphics_single_view)
    open my $outfh_snp, ">", $outputfile_SNP
      or do {
        warn "Cannot write $outputfile_SNP: $!\n";
        return 0;
      };

    open my $outs_snp_fh, "<", $out_single_SNP
      or do {
        warn "Cannot read $out_single_SNP: $!\n";
        close $outfh_snp;
        return 0;
      };

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

    while (<$outs_snp_fh>) {
        chomp;
        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start <= $end ) {    # BUG FIX: same orientation logic
            $orientation = +1;
        }
        else {
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

    binmode $outfh_snp;
    print $outfh_snp $panel->png;
    close $outs_snp_fh;
    close $outfh_snp;

    print "\nSingle View SNP Primer Graphic File Generated.\n";
    return 1;
}

####################################

=head1 uniq
 Title   :  uniq
 Usage   :  my @unique_length = uniq(@$array_length_SNP);
 Function:  used to reduce an array to unique entries only
 Returns :  new array 
=cut

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

####################################

1;
