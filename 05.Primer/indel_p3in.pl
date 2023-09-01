#! /usr/bin/perl

use strict;
use warnings;
use List::Util qw/max/;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --input indelWithSEQ.filter.txt --output p3in --config primer3_config
#
# Required:
#
#	--input <string>               input filename.
#
#	--output <string>              output p3in filename.
#
#	--config <string>              the full path of primer3_config.
#
#	--minProdLen <integer>         min product length, default: 80.
#
#	--maxProdLen <integer>         max product length, default: 180.
#
#	--primerMinSize <integer>      min primer size, default: 18.
#
#	--primerOptSize <integer>      optimal primer size, default: 20.
#
#	--primerMaxSize <integer>      max primer size, default: 27.
#
#	--primerMinGC <floating>       min GC content of primer, default: 35.0.
#
#	--primerMaxGC <floating>       max GC content of primer, default: 65.0.
#
#	--primerMinTM <floating>       min TM of primer, default: 55.0.
#
#	--primerOptTM <floating>       optimal TM of primer, default: 60.0.
#
#	--primerMaxTM <floating>       max TM of primer, default: 65.0.
#
#	--primerNumReturn <integer>    number of primer pairs returned, default: 3.
#
############################################################

__EOUSAGE__

    ;

my $help_flag;
my $input;
my $output;
my $config;
my $minProdLen = 80;
my $maxProdLen = 180;
my $primerMinSize = 18;
my $primerOptSize = 20;
my $primerMaxSize = 27;
my $primerMinGC = 35.0;
my $primerMaxGC = 65.0;
my $primerMinTM = 55.0;
my $primerOptTM = 60.0;
my $primerMaxTM = 65.0;
my $primerNumReturn = 3;

&GetOptions('help|h'            => \$help_flag,
            'input|i=s'         => \$input,
            'output|o=s'        => \$output,
            'config|c=s'        => \$config,
            'minProdLen=i'      => \$minProdLen,
            'maxProdLen=i'      => \$maxProdLen,
            'primerMinSize=i'   => \$primerMinSize,
            'primerOptSize=i'   => \$primerOptSize,
            'primerMaxSize=i'   => \$primerMaxSize,
            'primerMinGC=s'     => \$primerMinGC,
            'primerMaxGC=s'     => \$primerMaxGC,
            'primerMinTM=s'     => \$primerMinTM,
            'primerOptTM=s'     => \$primerOptTM,
            'primerMaxTM=s'     => \$primerMaxTM,
            'primerNumReturn=i' => \$primerNumReturn);

unless ($input && $output && $config) {
	die $usage;
}

open IN,"$input";
open OUT,">$output";
my @len;

while(<IN>){
	chomp;
	next if /^SEQ_500bp/;
	my @a=split/\s+/,$_;
	my @b=split(/\,/,$a[5]);
	my @len;
	$len[0]=length($a[4]);
	my $i=1;
	foreach(@b){
		$len[$i]=length($_);
		$i+=1;
	}
	my $max=0;
	foreach(@len){
		if($_>$max){$max=$_}
	}
	if($max>$len[0]){$max=$len[0]}
	print OUT "SEQUENCE_ID=$a[3]\nSEQUENCE_TEMPLATE=$a[0]\nSEQUENCE_TARGET=249,$max\nPRIMER_TASK=generic\nPRIMER_PRODUCT_SIZE_RANGE=$minProdLen-$maxProdLen\nPRIMER_MIN_SIZE=$primerMinSize\nPRIMER_OPT_SIZE=$primerOptSize\nPRIMER_MAX_SIZE=$primerMaxSize\nPRIMER_MIN_GC=$primerMinGC\nPRIMER_MAX_GC=$primerMaxGC\nPRIMER_MIN_TM=$primerMinTM\nPRIMER_OPT_TM=$primerOptTM\nPRIMER_MAX_TM=$primerMaxTM\nPRIMER_NUM_RETURN=$primerNumReturn\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=$config\n=\n";
}
close IN;
close OUT;
