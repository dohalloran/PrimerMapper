#!/usr/bin/perl 

use strict;
use warnings;

my $out_dir = 'Browser_PrimerMapper';
if ( !-e $out_dir ) {
    mkdir($out_dir)
      or die "Failed to create the output directory ($out_dir): $!\n";
}

my $gene_name = $ARGV[0];
my $outfile   = "CANVAS.txt";
open CANVASFILE, ">$outfile" or die "Couldn't open file: $!";

my %data;
my @names;
while (<>) {
    chomp;
    my @list = split(/\t/);
    for ( my $i = 0 ; $i <= $#list ; $i++ ) {
        if ( $. == 1 ) {
            $names[$i] = $list[$i];
        }
        else {
            push @{ $data{ $names[$i] } }, $list[$i];
        }
    }
}
foreach (@names) {
    local $" = "\t";
    print CANVASFILE "$_\t@{$data{$_}}\n";
}

close CANVASFILE;

END;

my $output_canvas = "$gene_name.CANVAS.txt";
open CANVASFILE, "<$outfile" or die "Couldn't open file: $!";
open OUTER, ">$out_dir/$output_canvas" or die "Couldn't open file: $!";

while (<CANVASFILE>) {
    s/,$//;
    s/,\t\D//;
    print OUTER $_;
}

close OUTER;
close CANVASFILE;

unlink $outfile;


