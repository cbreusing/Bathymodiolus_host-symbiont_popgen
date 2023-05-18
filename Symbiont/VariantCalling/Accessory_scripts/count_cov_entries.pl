#!/usr/bin/env perl

use warnings;
use strict;
use List::Util qw(sum);

open (LIST, $ARGV[0]) or die "Cannot open input file: $!\n";

my @data;
my $presence = 0;
my $absence = 0;
my $i = 3;

print "CHROM\tPos\tRef\tpresence\tabsence\n";

while (<LIST>) {
  if ($. > 1) {
  chomp;
  @data = split(/\t/, $_);
  print "$data[0]\t$data[1]\t$data[2]\t";
        while ($i <= scalar(@data)-2) {
        $presence = $presence + $data[$i];
        $absence = $absence + $data[$i+1];
        $i = $i + 2;
  }
  print "$presence\t$absence\n";
  $presence = 0;
  $absence = 0;
  $i = 3;
  }
}

close (LIST);