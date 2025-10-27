#!/bin/env perl
# stat mapping rate
# only suitable for paired-end read mapped by bowtie2 log

use strict;
use warnings;

my $logfile = join(".", $ARGV[0], "log");
#print "$logfile\n";

open LOG, "$logfile";

my $total_read;
my $alg_conc_1;
my $unp_alg_1;
my $unmapped;

my $datestring = localtime();
print STDERR "$datestring: reading $ARGV[0].log ......\n";
while(<LOG>){
	chomp;
	if(/^(\d+) reads; of these:$/){
		$total_read = $1 * 2;
	}
	if(/(\d+) \(.*\) aligned concordantly exactly 1 time$/){
		$alg_conc_1 = $1 *2;
	}
	if(/(\d+) \(.+\) aligned exactly 1 time$/){
		$unp_alg_1 = $1;
	}
	if(/(\d+) \(.+\) aligned 0 times$/){
		$unmapped = $1;
	}
}

my $mapped = $total_read - $unmapped;
my $unique = $alg_conc_1 + $unp_alg_1;
my $mapped_rate = $mapped / $total_read;
my $unique_rate = $unique / $total_read;

print "$ARGV[0],$total_read,$mapped,$mapped_rate,$unique,$unique_rate\n";
