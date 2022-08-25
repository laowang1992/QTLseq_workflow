#!/bin/env perl
use strict;
use warnings;

open CSV, "$ARGV[0]";

print "<table>";
my $n = 1;
while(<CSV>){
	chomp;
	my @a = split /,/, $_;
	if($n == 1){
		s/,/<\/th><th>/g;
		print "<thead><tr><th>$_</th></tr></thead>";
	}elsif($n == 2){
		s/,/<\/td><td>/g;
		print "<tbody><tr><td>$_</td></tr>";
	}else{
		s/,/<\/td><td>/g;
		print "<tr><td>$_</td></tr>";
	}
	$n += 1;
}
if($n == 2){
	print "</table>";
}elsif($n > 2){
	print "</tbody></table>";
}

close CSV;
