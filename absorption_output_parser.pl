#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;
use MongoDB;

open my $fh, "<", $ARGV[0];
my $raw_data = join("", <$fh>);
close $fh;

my @data = split("\n\n\n", $raw_data);

my $client = MongoDB->connect('mongodb://localhost');
my $collection = $client->ns('qso_data.raw_data');

for(my $i = 0; $i < scalar @data; $i++){
	my ($name_raw, $BI_raw, $vmins_raw, $vmaxs_raw, $BIs_raw, $EWs_raw, $depths_raw) = split("\n", $data[$i]);

	my ($plate, $mjd, $fiberid) = $name_raw =~ /([0-9]{4})-([0-9]{5})-([0-9]{4})/;
	my ($BI)  = $BI_raw =~ /\:\ ([0-9\.\-]*)/;

	my ($vmins_list) = $vmins_raw =~ /\:\ \[(.*)\]/;
	my ($vmaxs_list) = $vmaxs_raw =~ /\:\ \[(.*)\]/;
	my ($BIs_list) = $BIs_raw =~ /\:\ \[(.*)\]/;
	my ($EWs_list) = $EWs_raw =~ /\:\ \[(.*)\]/;
	my ($depths_list) = $depths_raw =~ /\:\ \[(.*)\]/;

	my @vmins = map {0+ $_} split(", ", $vmins_list);
	my @vmaxs = map {0+ $_} split(", ", $vmaxs_list);
	my @BI_individual = map {0+ $_} split(", ", $BIs_list);
	my @EW_individial = map {0+ $_} split(", ", $EWs_list);
	my @depth = map {0+ $_} split(", ", $depths_list);

	my $final_data = {
		SID => $plate."-".$mjd."-".(0+ $fiberid),
		PLATE => 0+ $plate,
		MJD => 0+ $mjd,
		FIBERID => 0+ $fiberid,
		BI => 0+ $BI,
		vmins => \@vmins,
		vmaxs => \@vmaxs,
		BI_individual => \@BI_individual,
		EW_individial => \@EW_individial,
		depth => \@depth,
		file => $name_raw
	};

	$collection->insert($final_data);
}