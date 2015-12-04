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
The ouput from this package produces graphical primer::sequence
maps of primer distribution across each sequence and the entire input file

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
    print "\nDesigning Grahical output files....\n";
    my $dir = cwd();
    my $orientation;
    my @gene_id;
    my @unique_length = uniq(@$array_length);
    my @unique_name   = uniq(@$array_name);
    my %hash;
    my $len;
    my $name_id;

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

    foreach my $fp ( glob("$dir/GRAPHIC_*.txt") ) {
        open my $fh, "<", $fp or die;

        @hash{@gene_id} = @$gene_length;
        while ( my ( $key, $val ) = each %hash ) {
            $len = $val if $fp =~ m/$key/i;
        }
        while ( my ( $key, $val ) = each %hash ) {
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

        while (<$fh>) {
            chomp;

            my ( $name, $score, $start, $end ) = split /\t/;
            if ( $start - $end < 1 ) {
                $orientation = +1;
            }
            elsif ( $start - $end > 1 ) {
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
    print "\nDesigning Single View Grahical output....\n";
    my $dir = cwd();
    my $orientation;
    my $lener     = eval join '+', @$gene_length;
    my $first_seq = @$gene_length[0];
    my $max       = max(@$gene_length);

    my $outputfile = "single_view" . $format;
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

    while (<OUTS>) {
        chomp;

        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start - $end < 1 ) {
            $orientation = +1;
        }
        elsif ( $start - $end > 1 ) {
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
    print "\nDesigning SNP Primer Grahical output files....\n";
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

    foreach my $fp ( glob("$dir/GRAPHIC_*.txt") ) {
        open my $fh, "<", $fp or die;

        @hash{@gene_id} = @$gene_length_SNP;
        while ( my ( $key, $val ) = each %hash ) {
            $len = $val if $fp =~ m/$key/i;
        }
        while ( my ( $key, $val ) = each %hash ) {
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

        while (<$fh>) {
            chomp;

            ( $snp, $score, $start, $end ) = split /\t/;
            if ( $start - $end < 1 ) {
                $orientation = +1;
            }
            elsif ( $start - $end > 1 ) {
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
    print "\nDesigning Single View SNP Primer Grahical output file....\n";
    my $dir = cwd();
    my $orientation;
    my $lener     = eval join '+', @$gene_length_SNP;
    my $first_seq = @$gene_length_SNP[0];
    my $max       = max(@$gene_length_SNP);

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

    while (<OUTS_SNP>) {
        chomp;

        my ( $name, $score, $start, $end ) = split /\t/;

        if ( $start - $end < 1 ) {
            $orientation = +1;
        }
        elsif ( $start - $end > 1 ) {
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
