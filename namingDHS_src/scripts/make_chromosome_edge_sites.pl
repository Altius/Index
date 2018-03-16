#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#          Generate ONE fake 1-bp DHS at each end for the namer, 
#              bedops -u with real input master list,
#              and removed by -n 100% on output, with names not even reserved?
#              but some namespace reserved.
#              Hmmm... just 1 placeholder at the very ends could be kept for naming but how to keep track?
#              WE keep track, because only WE assign the numbers.
#              So we reserve those numbers and their existence puts some elbow room around to avoid too-long names.
#
# equivalently, 
#       awk '{print $1"\t0\t1\n"$1"\t"($2-1)"\t"$2}' chrom.sizes
 

unless (@ARGV > 0) {
    die "Usage:  $0 genome.chrom.sizes\n";
}

my $chromSizesFile = $ARGV[0];
open (my $in, '<', $chromSizesFile) or die "$!";
while (<$in>) {
    chomp;
    my ($chrom, $size) = split /\t/;
    print join("\t", ($chrom, 0, 1)) . "\n";
    print join("\t", ($chrom, $size - 1, $size)) . "\n";
}
close $in;


