#!/usr/bin/perl

use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

$ua = LWP::UserAgent->new;

$argc = $#ARGV + 1;

$program  = shift;
$database = shift;

if ( $program eq "megablast" ) {
    $program = "blastn&MEGABLAST=on";
}

if ( $program eq "rpsblast" ) {
    $program = "blastp&SERVICE=rpsblast";
}

# read and encode the queries
foreach $query (@ARGV) {
    open( QUERY, $query );
    while (<QUERY>) {
        $encoded_query = $encoded_query . uri_escape($_);
    }
}

# build the request
$args = "CMD=Put&PROGRAM=$program&DATABASE=$database&QUERY=" . $encoded_query;

$req = new HTTP::Request POST => 'https://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);

# get the response
$response = $ua->request($req);

# parse out the request id
$response->content =~ /^    RID = (.*$)/m;
$rid = $1;

# parse out the estimated time to completion
$response->content =~ /^    RTOE = (.*$)/m;
$rtoe = $1;

# wait for search to complete
sleep $rtoe;

# poll for results
while (true) {
    sleep 5;

    $req = new HTTP::Request GET =>
"https://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid&ALIGNMENTS=5&DESCRIPTIONS=5";
    $response = $ua->request($req);

    if ( $response->content =~ /\s+Status=WAITING/m ) {

        # print STDERR "Searching...\n";
        next;
    }

    if ( $response->content =~ /\s+Status=FAILED/m ) {
        print STDERR
"Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
    }

    if ( $response->content =~ /\s+Status=UNKNOWN/m ) {
        print STDERR "Search $rid expired.\n";
        exit 3;
    }

    if ( $response->content =~ /\s+Status=READY/m ) {
        if ( $response->content =~ /\s+ThereAreHits=yes/m ) {

            #  print STDERR "Search complete, retrieving results...\n";
            last;
        }
        else {
            print STDERR "No hits found.\n";
            exit 2;
        }
    }

    # if we get here, something unexpected happened.
    exit 5;
}    # end poll loop

# retrieve and display results
$req = new HTTP::Request GET =>
"https://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid&ALIGNMENTS=5&DESCRIPTIONS=5";
$response = $ua->request($req);

my $blastfile = "blast_primers.fasta";
my $outfile   = "blast_primers_results.fasta";
open( FILE, "<$blastfile" ) or print "unable to open file";
open( OUT,  ">$outfile" )   or print "unable to open file";
print OUT $response->content;
exit 0;

######################################################################
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
# This script is a modification of a public release provided by the
# National Center for Biotechnology Information.
# The original script is freely available to the public for use;
# The National Library of Medicine and the U.S. Government have not
# placed any restriction on its use or reproduction.
######################################################################

