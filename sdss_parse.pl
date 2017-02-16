#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use MongoDB;

my $client = MongoDB->connect('mongodb://localhost');
my $collection = $client->ns('qso_data.lines');

open my $cfg_file, "<", $ARGV[0];

while(<$cfg_file>){
	my ($PLATEID, $MJD, $FIBER) = $_ =~ /spec\-([0-9]+)\-([0-9]+)\-([0-9]+)/;

	my $url = "http://dr9.sdss3.org/spectrumDetail?mjd=$MJD&fiber=$FIBER&plateid=$PLATEID";
	#print $url."\n";
	my $raw_data = `curl "$url"`;

	my @t = split(/<tbody>|<\/tbody>/, $raw_data);

	foreach my $entry (split(/<\/tr>/, $t[-2])){
	my ($name) = $entry =~ /<\!\-\- Line Name \-\-><td class\="name">(.*)<\/td>/;
	my ($wavelength) = $entry =~ /<\!\-\- Wavelength \-\-><td>(.*)<\/td>/;
	my ($z) = $entry =~ /<\!\-\- linez \-\-><td>-(.*)<\/td>/;
    my ($sigma) = $entry =~ /<\!\-\- linesigma \-\-><td>(.*)<\/td>/;
    my ($area) = $entry =~ /<\!\-\- linearea \-\-><td>(.*)<\/td>/;
    my ($contlevel) = $entry =~ /<\!\-\- linecontlevel \-\-><td>(.*)<\/td>/;

		if($name && defined $z){
			#Convert wavelength to v, for refrence purposes.
			my $CIV = 1549;
			my $abs_z = ($wavelength/$CIV) - 1;
			my $RC = (1 + $z)/(1 + $abs_z);
			my $betaC = (($RC**2) - 1)/(($RC**2) + 1);
			my $beta = $betaC * (-300000);

			my $data = {
				name => $name,
				wavelength => 0+ $wavelength,
				z => ($z? 0+ $z : 0),
				sigma => 0+ $sigma,
				area => 0+ $area,
				contlevel => 0+ $contlevel,
				SID => (0+ $PLATEID)."-".(0+ $MJD)."-".(0+ $FIBER),
				beta => 0+ $beta
			};

			$collection->insert($data);

			#print Dumper $data;
		}
	}
}
