#!/usr/bin/env perl

use warnings;
use strict;

open (FILE1, $ARGV[0]) or die "Cannot open input file: $!\n";

my %hash;
my @list;
my $key;
my $count;
my $diff;

while (<FILE1>) {
  chomp;
  @list = split(/\t/, $_);
  $hash{$list[0]} = $list[6];
}

open (FILE2, $ARGV[1]) or die "Cannot open input file: $!\n";
$count = do {local $/; <FILE2> };

close (FILE1);
close (FILE2);

open (OUTFILE, ">$ARGV[2]") or die "Cannot open output file: $!\n";

print OUTFILE "CHROM\tPos\tRef\tpresence\tabsence\n";

foreach $key (sort keys(%hash)) {
    $diff = $count - $hash{$key};
    if ($diff<0) {
    $diff = 0;
    }
    print OUTFILE "$key\t$key\tpresence\t$hash{$key}\t$diff\n";
}

close (OUTFILE);