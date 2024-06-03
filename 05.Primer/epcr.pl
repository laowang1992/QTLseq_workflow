#!/bin/env perl
# Date: 2024-05-23
# Usage: perl epcr.pl --input p3out --output primer.txt [--timeout 2]

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --input p3out --output primer.txt [--timeout 2]
#
# Required:
#
#	--input <string>       input p3in filename.
#
#	--output <string>      output table filename.
#
# Optional:
#
#	--timeout <int>        timeout for each e-PCR command (default: 2 seconds)
#
############################################################

__EOUSAGE__

my $help_flag;
my $input;
my $output;
my $timeout = 2;  # 默认超时时间

&GetOptions('help|h' => \$help_flag,
            'input|i=s' => \$input,
            'output|o=s' => \$output,
            'timeout|T=i' => \$timeout,
            );

unless ($input && $output) {
    die $usage;
}

open my $p3out_fh, '<', $input or die "Cannot open $input: $!";
open my $output_fh, '>', $output or die "Cannot open $output: $!";

print $output_fh "ID\tNUM\tPrimerL\tPrimerR\tHit\tTmL\tTmR\tGCL\tGCR\tLength\n";
while (my %primers = get_next_primers($p3out_fh)) {
    for (my $i = 1; $i <= $primers{'number'}; ++$i) {
        my $id      = $primers{'id'};
        my $num     = $i - 1;
        my $primerL = @{$primers{'p_left'}}[$num];
        my $primerR = @{$primers{'p_right'}}[$num];
        my $tmL     = @{$primers{'tm_left'}}[$num];
        my $tmR     = @{$primers{'tm_right'}}[$num];
        my $gcL     = @{$primers{'gc_left'}}[$num];
        my $gcR     = @{$primers{'gc_right'}}[$num];
        my $len     = @{$primers{'size'}}[$num];
        my @epcr;
        my $hit;
        my $pid;
        eval {
            # 设置alarm，超时时间为$timeout秒
            local $SIG{ALRM} = sub {die "timeout\n"};
            alarm $timeout;

            # 执行外部命令，由readpipe改为打开句柄，并捕获re-PCR的pid，方便超时发生时杀死外部进程
            $pid = open my $cmd_fh, "-|", "re-PCR -s genome.hash -n 1 -g 1 $primerL $primerR 50-1000"
                or die "Cannot fork: $!";
            
            @epcr = <$cmd_fh>;
            close $cmd_fh;

            # 取消alarm
            alarm 0;
        };

        if ($@) {
            if ($@ eq "timeout\n") {
                print "$id $num : e_PCR command timed out, skipping this iteration\n";
                $hit = "timeout";
                # 杀死外部命令进程
                kill 9, $pid if $pid;
                # 原来发生超时会使用next跳过这个引物，现在继续处理，只是把$hit赋值"timeout"，注意这里会导致结果文件中Hit列数字和字符混合，后续步骤读取时要注意这个问题
            } else {
                die $@;  # 处理其他可能的异常
            }
        } else {
            $hit = @epcr;
            $hit -= 2;
        }

        print $output_fh "$id\t$num\t$primerL\t$primerR\t$hit\t$tmL\t$tmR\t$gcL\t$gcR\t$len\n"
    }
}

close $p3out_fh;
close $output_fh;

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
