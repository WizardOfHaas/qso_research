#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

my $MJD = 55475;
my $FIBER = 210;
my $PLATEID = 4207;

my $url = "http://dr9.sdss3.org/spectrumDetail?mjd=$MJD&fiber=$FIBER&plateid=$PLATEID";
print $url."\n";
my $raw_data = `curl "$url"`;

my @t = split(/<tbody>|<\/tbody>/, $raw_data);

foreach my $entry (split(/<\/tr>/, $t[-2])){
	my ($name) = $entry =~ /<\!\-\- Line Name \-\-><td class\="name">(.*)<\/td>/;
	my ($wavelength) = $entry =~ /<\!\-\- Wavelength \-\-><td>(.*)<\/td>/;
	my ($z) = $entry =~ /<\!\-\- linez \-\-><td>-(.*)<\/td>/;
        my ($sigma) = $entry =~ /<\!\-\- linesigma \-\-><td>(.*)<\/td>/;
        my ($area) = $entry =~ /<\!\-\- linearea \-\-><td>(.*)<\/td>/;
        my ($contlevel) = $entry =~ /<\!\-\- linecontlevel \-\-><td>(.*)<\/td>/;

	if($name){
		my $data = {
			name => $name,
			wavelength => 0+ $wavelength,
			z => ($z? 0+ $z : 0),
			sigma => 0+ $sigma,
			area => 0+ $area,
			ccontlevel => 0+ $contlevel
		};

		print Dumper $data;
	}
}
