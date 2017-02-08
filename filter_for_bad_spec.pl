#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

open my $SPEC_LIST, "<confs/DR9Q_selection_specnames.lst";
open my $GOOD_SPECS, "<confs/good_categorized_spectra_through_code_subtracted_one.txt";

my @specs = ();

while(<$SPEC_LIST>){
	chomp;
	my $spec = $_; #=~ /([0-9]{4}\-[0-9]{5}\-[0-9]{4})/;
	push(@specs, $spec);
}

close $SPEC_LIST;

my @good_specs = ();

while(<$GOOD_SPECS>){
	push(@good_specs, $specs[$_]);
}

close $GOOD_SPECS;

my %good_specs = map{$_ => 1} @good_specs;
my @bad_specs = grep(!defined $good_specs{$_}, @specs);

print join("\n", @bad_specs);