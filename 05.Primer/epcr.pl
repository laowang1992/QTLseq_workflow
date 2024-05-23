#!/bin/env perl
# Date: 2024-05-23
# Usage: perl epcr.pl xxx.p3out output.txt

use strict;
use warnings;

open my $p3out_fh, '<', $ARGV[0] or die "Cannot open $ARGV[0]: $!";
open my $output, '>', $ARGV[1] or die "Cannot open $ARGV[1]: $!";

# 设置超时时间
my $timeout = 2;

print $output "ID	NUM	PrimerL	PrimerR	Hit	TmL	TmR	GCL	GCR	Length\n";
while (my %primers = get_next_primers($p3out_fh)) {
	for (my $i = 1; $i <= $primers{'number'}; ++$i) {
		my $id      = $primers{'id'};
		my $num     = $i -1;
		my $primerL = @{$primers{'p_left'}}[$num];
		my $primerR = @{$primers{'p_right'}}[$num];
		my $tmL     = @{$primers{'tm_left'}}[$num];
		my $tmR     = @{$primers{'tm_right'}}[$num];
		my $gcL     = @{$primers{'gc_left'}}[$num];
		my $gcR     = @{$primers{'gc_right'}}[$num];
		my $len     = @{$primers{'size'}}[$num];
		my @epcr;
		my $hit;
		eval {
			# 设置alarm，超时时间为$timeout秒
			local $SIG{ALRM} = sub {die "timeout\n"};
			alarm $timeout;
			
			# 执行外部命令
			@epcr = readpipe("re-PCR -s genome.hash -n 1 -g 1 $primerL $primerR 50-1000");
			
			# 取消alarm
			alarm 0;
		};
		
		if ($@) {
			if ($@ eq "timeout\n") {
				print "$id $num : e_PCR command timed out, skipping this iteration\n";
				$hit = "timeout";
				#next;  # 跳过当前循环，进入下一个循环
			} else {
				die $@;  # 处理其他可能的异常
			}
		} else {
			$hit = @epcr;
			$hit -= 2;
		}
		
		print $output "$id\t$num\t$primerL\t$primerR\t$hit\t$tmL\t$tmR\t$gcL\t$gcR\t$len\n"
		
	}
}

close $p3out_fh;
close $output;

sub get_next_primers {
	my ($fh) = @_;
	my %primers = (
		'p_left'    => [],
		'p_right'   => [],
		'tm_left'   => [],
		'tm_right'  => [],
		'gc_left'   => [],
		'gc_right'  => [],
		'size'      => [],
	);
	while (my $line = <$fh>) {
		chomp $line;
		# SEQUENCE_ID不能有空格
		if ($line =~ /^SEQUENCE_ID=(.+)$/) {
			$primers{'id'} = $1;
		}
		if ($line =~ /^PRIMER_LEFT_NUM_RETURNED=(\d+)$/) {
			$primers{'number'} = $1;
		}
		if ($line =~ /^PRIMER_LEFT_\d+_SEQUENCE=(.+)$/){
			push @{$primers{'p_left'}}, $1;
		}
		if ($line =~ /^PRIMER_RIGHT_\d+_SEQUENCE=(.+)$/){
			push @{$primers{'p_right'}}, $1;
		}
		if ($line =~ /^PRIMER_LEFT_\d+_TM=(.+)$/) {
			push @{$primers{'tm_left'}}, $1;
		}
		if ($line =~ /^PRIMER_RIGHT_\d+_TM=(.+)$/) {
			push @{$primers{'tm_right'}}, $1;
		}
		if ($line =~ /^PRIMER_LEFT_\d+_GC_PERCENT=(.+)$/) {
			push @{$primers{'gc_left'}}, $1;
		}
		if ($line =~ /^PRIMER_RIGHT_\d+_GC_PERCENT=(.+)$/) {
			push @{$primers{'gc_right'}}, $1;
		}
		if ($line =~ /^PRIMER_PAIR_\d+_PRODUCT_SIZE=(.+)$/) {
			push @{$primers{'size'}}, $1;
		}
		if ($line =~ /^=$/) {
			return %primers;
		}
	}
	return ();
}

