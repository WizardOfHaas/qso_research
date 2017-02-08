#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;
use MongoDB;

open my $file, "<", $ARGV[0];
my @specs = <$file>;

my $client = MongoDB->connect('mongodb://localhost');
my $collection = $client->ns('qso_data.dr9q');

foreach my $spec (@specs){
	chomp $spec;
	my ($plate, $mjd, $fiberid) = $spec =~ /([0-9]{4})-([0-9]{5})-([0-9]{4})/;

	#print $plate."-".$mjd."-".$fiberid."\n";
	
	my $result = $collection->find({
		PLATE => 0+ $plate,
		MJD => 0+ $mjd,
		FIBERID => 0+ $fiberid
	});

	while(my $data = $result->next){
		print $spec.",".$data->{Z_VI}.",".$data->{SNR_1700}."\n";
	}
}