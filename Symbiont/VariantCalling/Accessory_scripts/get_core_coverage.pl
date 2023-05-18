#!/usr/bin/env perl

use warnings;
use strict;

open (LIST1, $ARGV[0]) or die "Cannot open input file: $!\n";
open (LIST2, $ARGV[1]) or die "Cannot open input file: $!\n";

my %hash1;
my %hash2;
my @list1;
my @list2;
my $key;

while (<LIST1>) {
  chomp;
  @list1 = split(/\t/, $_);
  undef $hash1{$_};
}

while (<LIST2>) {
  chomp;
  @list2 = split(/\t/, $_);
  $hash2{$list2[0]} = join("\t", @list2[1..$#list2]);
}

close (LIST1);
close (LIST2);

open (OUTFILE, ">$ARGV[2]") or die "Cannot open output file: $!\n";

foreach $key (sort keys(%hash1)) {
  if (exists $hash2{$key}) {
    print OUTFILE "$key\t$hash2{$key}\n";
  }  
}