#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use List::Util qw( min );
 
#  Assign DHS identifiers for a semi-human-readable index.

############  Paths to required programs. ##########
my $thisScriptFile = $0;
chomp( my $bedops = `which bedops` );
chomp( my $bedmap = `which bedmap` );
chomp( my $sort_bed = `which sort-bed` );
chomp( my $closest_features = `which closest-features` );
unless (-x $bedops) {
    die "Failed to find bedops program\n";
}

chomp( my $scriptdir = `dirname $0` ); # a relative path is fine too
my $makeChromosomeEdgeSites = "$scriptdir/make_chromosome_edge_sites.pl";

#          Generate ONE fake 1-bp DHS at each end for the namer, 
#              bedops -u with real input master list,
#              and removed by -n 100% on output, with names not even reserved?
#              but some namespace reserved.
#              Hmmm... just 1 placeholder at the very ends could be kept for naming but how to keep track?
#              WE keep track, because only WE assign the numbers.
#              So we reserve those numbers and their existence puts some elbow room around to avoid too-long names.

##########  Decimal numbering, but could extend to anything else. ##########
my @symbols = ();
my $base = 10;
for (my $i = 0; $i < $base; ++$i) {
    push @symbols, $i;
}
my $symbolString = join("", @symbols );
my $zeroSymbol = $symbols[0]; 

##########  Usage in two distinct modes, new or updated label assignments. ##########
my $usageText =<<EOusage;

Usage:   $thisScriptFile  [options] | sort-bed - > labelled.txt

Options:
### Mode 1 :  denovo   Lay out numbers proportionally to a decent approximation of DHS space ### 
INPUT:
    --sizes    hg38.chrom.sizes
    --regions  bed/starch file of DHSs

    Prefix will be the chromosome name without leading 'chr', followed by a dot.

OUTPUT:
    STDOUT     bed4 with labels, requires sorting
    --edgefile Filename to write pseudo-DHS at chrom boundaries for future labelling (SAVE THIS FILE!)

### Mode 2 :  update   Make new labels consistent with previous ones ### 
INPUT:
    --regions  bed/starch file of new DHSs
    --previous bed/starch (BED4) of previously-labelled DHSs
    --retiredIDs optional text file of retired IDs that may or may not be in the previous BED file 
    --edgefile From first denovo pass, possibly edited if additional chromosomes in new genome build

    Prefix will be taken from adjacent DHSs or pseudo-DHSs via closest-features.

OUTPUT:
    STDOUT     bed4 with labels, requires sorting
EOusage

# TODO? [--overlap] bedmap overlap options to determine when to re-use an existing label (ONCE!)
#   (Currently requiring exact overlap).
# For master-list DHSs, at least be stringent enough that no two input regions will get the same old name!

##########  Read command-line options. ##########
my $chromSizesFile;
my $assignRegionsFile;
my $previousRegionsFile;
my $edgeCaseLabelsFile;
my $retiredIDFile;

GetOptions( "sizes=s" => \$chromSizesFile,
            "regions=s" => \$assignRegionsFile,
            "previous=s" => \$previousRegionsFile,
            "edgefile=s" => \$edgeCaseLabelsFile,
            "retiredIDs=s" => \$retiredIDFile )
or die $usageText;


# Two options always required, the regions to label and a special file for handling chromosome edge cases.
unless (defined($assignRegionsFile)) {
    die "Error:  No regions file specified for label assignment.\n$usageText";
}
unless (defined($edgeCaseLabelsFile)) {
    die "Error:  No edge file specified, required for input or output depending on assignment mode.\n$usageText";
}

#Read optional retired ids here before doing anything else
my %retiredIDs = ();

# Existence of a previous file determines which mode we run in.
if (defined( $previousRegionsFile )) {
    # Try to make new labels fit nicely in between the previous labels.
    unless (-e $previousRegionsFile) {
        die "Failed to find $previousRegionsFile\n";
    }
    # Labels to avoid re-assigning 
    readRetiredIDs();
    # Make new labels consistent with exact-matching and neighboring previous regions
    updateLabelAssignments();

} else {
    # Try to make new label numbers spread out evenly across the genome. 
    unless (defined($chromSizesFile)) {
        die "Error:  No genome.chrom.sizes file specified.\n$usageText\n";
    }
    # Need the fake edge markers sorted into the real input regions, but diverted on output.
    # Could do this more cryptically with awk instead of external script,
    # e.g. "awk '{print \$1\"\t0\t1\n\"\$1\"\t\"(\$2-1)\"\t\"\$2}' chrom.sizes";
    my %regionsByPrefix = ();
    my $cmd = "$makeChromosomeEdgeSites $chromSizesFile | sort-bed - | $bedops -u - $assignRegionsFile"; 
    open (my $in, "$cmd |") or die "$!";
    while (<$in>) {
        chomp;
        my $line = $_;
        my ($chrom,$min0,$max1,@etc) = split /\t/;
        my $labelPrefix = $chrom;
        if ($chrom =~ /^chr(.*)$/) {
            $labelPrefix = $1;
        }
        $labelPrefix .= ".";
        my $location = join("\t", ($chrom, $min0, $max1) );
        push @{$regionsByPrefix{$labelPrefix}}, $location;
    }
    close $in;
    # Do the thing
    if (-e $edgeCaseLabelsFile) {
        die "Warning:  $edgeCaseLabelsFile would be overwritten!\n"
            . "    Please move it out of the way first,\n"
            . "    because we don't want to clobber it by accident here.\n";
    }
    open (my $outEdges, "| $sort_bed - > $edgeCaseLabelsFile") or die "$!";
    for my $labelPrefix (sort keys %regionsByPrefix) {
        recursiveHierarchicalNumbering( $regionsByPrefix{$labelPrefix}, $symbolString, $labelPrefix, $outEdges, 1 );
    }
    close $outEdges;
}

######## Assign new labels consistent with existing labels. #########
sub updateLabelAssignments {
    my %fillIns = ();
    my $cmd = "$bedops -u $previousRegionsFile $edgeCaseLabelsFile "
        . "| $closest_features $assignRegionsFile - "
        . "| $bedmap --echo --echo-map --exact - $previousRegionsFile ";
    open (my $in, "$cmd |") or die "$!";
    while (<$in>) {
        chomp;
        my $mappedLine = $_;
        my ($element,$prev5,$prev3,$exactMatch) = split /\|/;
        #NOT USED:  my ($elementChrom,$elementMin0,$elementMax1) = split /\t/, $element;
        if (defined($exactMatch) and (length($exactMatch)>0)) {
            # If the match criteria was made more lax, we still want the old region for the old label
            print "$exactMatch\n";
        } else {
            # Assign a new name in between the nearest flanking old names
            my ($prev5Chrom,$prev5Min0,$prev5Max1,$prev5Name) = split /\t/, $prev5;
            my ($prev3Chrom,$prev3Min0,$prev3Max1,$prev3Name) = split /\t/, $prev3;
            unless (defined($prev5Name) and defined($prev3Name)) {
                die "Missing a named neighbor in ($mappedLine)\n";
            }
            my ($commonPrefix, $uniqueLeft, $uniqueRight) = compareFlankingIDs( $prev5Name, $prev3Name );
            my $key = "$commonPrefix|$uniqueLeft|$uniqueRight"; # just the ids would work too
            push @{$fillIns{$key}}, $element;
        }
    }
    close $in;

    foreach my $context (keys %fillIns) {
        my $listref = $fillIns{$context};
        my ($commonPrefix, $uniqueLeft, $uniqueRight) = split /\|/, $context;
        my @list = @$listref;
        my $numvals = scalar(@list);
        my @symsL = split //, $uniqueLeft;
        my @symsR = split //, $uniqueRight;
        # Zero-extend the boundary values to least common depth
        while (scalar(@symsL) < scalar(@symsR)) {
            push @symsL, $zeroSymbol;
        }
        while (scalar(@symsR) < scalar(@symsL)) {
            push @symsR, $zeroSymbol;
        }
        # Assuming decimal digits because it makes it really easy.
        # But can generalize if time permits, preferably in Python...
        my $limitL = int( join( '', @symsL ) );
        my $limitR = int( join( '', @symsR ) );
        if ($limitR < $limitL) {
            # Could happen if assembly changed
            ($limitL, $limitR) = ($limitR, $limitL); 
            # If these ids are far away or random, interpolation will still work just be kinda nonsensical.
        }
        my @middlestrings = ();
        for (my $i = $limitL + 1; $i < $limitR; ++$i) {
            push @middlestrings, "$i"; 
        }
        @middlestrings = excludeUsedIDs( $commonPrefix, \@middlestrings );
        while (scalar(@middlestrings) < $numvals) {
            @middlestrings = ();
            $limitL *= 10;
            $limitR *= 10;
            for (my $i = $limitL + 1; $i < $limitR; ++$i) {
                push @middlestrings, "$i"; 
            }
            @middlestrings = excludeUsedIDs( $commonPrefix, \@middlestrings );
        }
        my @indices = getIndicesForSplittingList( scalar(@list), scalar(@middlestrings) );
        my $prevIndex = 0;
        for (my $i = 0; $i < scalar(@list); ++$i) {
            my $region = $list[$i];
            my $nextIndex = $indices[$i];
            # If there's a lot of space, we don't want new IDs too close to old ones or each other
            # That way we limit the likelihood of having to make longer IDs.
            my $useIndex = int( ($prevIndex + $nextIndex) / 2 );
            my $newID = $commonPrefix . $middlestrings[$useIndex];
            print "$region\t$newID\n";
            $prevIndex = $nextIndex;
        }
    }
}

####### Remove retired labels from a list of candidates for assignment. #######
sub excludeUsedIDs {
    my ($prefix, $suffixlistref) = @_;
    my @result = ();
    foreach my $candidate (@$suffixlistref) {
        while ($candidate =~ /^(.*)$zeroSymbol$/) {
            $candidate = $1;
        }
        unless( defined( $retiredIDs{"$prefix$candidate"} ) ) {
            push @result, $candidate;
        }
    }
    return @result;
}

####### Find the maximal common prefix of neighboring identifiers. #######
sub compareFlankingIDs {
    my ( $prev5, $prev3 ) = @_;
    my $minlen = min( length($prev5), length($prev3) );
    my $common = 0;
    while (($common < $minlen) && (substr($prev5,$common,1) eq substr($prev3,$common,1))) {
        ++$common;
    }
    my $commonPrefix = substr( $prev5, 0, $common );
    die "assertion failed" unless ($commonPrefix eq substr( $prev3, 0, $common ));
    my $uniqueLeft = substr( $prev5, $common );
    my $uniqueRight = substr( $prev3, $common );
    #print "common( $prev5, $prev3 ) = ($commonPrefix, $uniqueLeft, $uniqueRight)\n";
    return ($commonPrefix, $uniqueLeft, $uniqueRight);
}

####### Method for spreading out namespace over a big region list the very first time #######
sub recursiveHierarchicalNumbering {
    my ($regionsRef, $symbolString, $prefix, $fileHandleForEdges, $startNumberingWithOne ) = @_;
    my @list = @$regionsRef;
    if (@list < 1) {
        # NOP
    } elsif (@list == 1) {
        my $line = $list[0];
        my $name = removeTrailingZeros( $prefix );
        my ($chrom,$min0,$max1) = split /\t/, $line;
        if ($max1 - $min0 == 1) {
            # Fake edge case, keeping out of main output
            print $fileHandleForEdges "$line\t$name\n";
        } else {
            # just writing to stdout for now
            print "$line\t$name\n";
        }
    } else {
        my @indices = getIndicesForSplittingList( length($symbolString) - $startNumberingWithOne, scalar(@list) );
        my $prevIndex = 0;
        for (my $i = 0; $i < scalar(@indices); ++$i) {
            my $nextIndex = $indices[$i];
            my @slice = @list[$prevIndex..($nextIndex-1)];
            my $symbol = substr( $symbolString, $i + $startNumberingWithOne, 1 );
            recursiveHierarchicalNumbering( \@slice, $symbolString, "$prefix$symbol", $fileHandleForEdges, 0 );
            $prevIndex = $nextIndex;
        }
    }
}

#######  Spread a small list more-or-less evenly across a larger list.  #######
sub getIndicesForSplittingList {
    my ($numGroups, $origLength) = @_;
    my $groups = $numGroups;
    my $remaining = $origLength;
    my $index = 0;
    my @result = ();
    while ($groups > 0) {
        my $size = int( 0.5 + ($remaining / $groups) );
        $index += $size;
        $remaining -= $size;
        push @result, $index;
        --$groups;
    }
    return @result;
}


#######  Avoid re-using IDs already assigned in any previous version.  #######
sub readRetiredIDs {
    # optional list of things that might not be in the old-id BED file
    if (defined( $retiredIDFile )) {
        open (my $ids, '<', $retiredIDFile) or die "$!";
        while (<$ids>) {
            chomp;
            $retiredIDs{removeTrailingZeros($_)} = 1;
        }
        close $ids;
    }
    # Also add ids from the previous file and the edge-case file
    for my $bedFile ($previousRegionsFile, $edgeCaseLabelsFile) {
        open (my $ids, "$bedops -u $bedFile |") or die "$!";
        while (<$ids>) {
            chomp;
            my ($chrom, $min0, $max1, $prevId) = split /\t/;
            if (defined($prevId)) {
                $retiredIDs{removeTrailingZeros($prevId)} = 1;
            }
        }
        close $ids;
    }
}

# Like representing the decimal fraction of real numbers, 
# adding zeroes doesn't change the ID number to be something different.
sub removeTrailingZeros {
    my $prevId = shift;
    while ($prevId =~ /^(.*)$zeroSymbol$/) {
        $prevId = $1;
    }
    return $prevId;
}

