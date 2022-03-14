#!/bin/env perl
use strict;
use warnings;

my $conf = $ARGV[0];
my $sample = $ARGV[1];
my $mem = $ARGV[2];

if(0){
	$sample = "HHH";
	$mem = 10;
	$conf = "../script/.conf"
}

my $outfile;
$outfile = join ".", $sample, "gatk.sh";
open OUT, ">$outfile";

#open OUT, ">$ARGV[2]";
print OUT "#PBS -N xx\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=${mem}gb\n\n";

print OUT ". $conf\n\n";

print OUT "cd \${work_dir}\/01.Mapping\n";

print OUT "java -Xmx${mem}g -jar \${picard} \\\n\tMarkDuplicates I=${sample}.sort.bam O=${sample}.dd.bam \\\n\tCREATE_INDEX=true REMOVE_DUPLICATES=true \\\n\tM=${sample}.dd.metics\n";

print OUT "samtools index ${sample}.dd.bam\n";

print OUT "genomeCoverageBed -ibam ${sample}.dd.bam -bga -g \${genome} | grep -w \"0\$\" > ${sample}.0cov.bedgraph\n\n";

print OUT "cd \${work_dir}/02.SNP_indel\n";

print OUT "java -Xmx${mem}g -jar \${gatk} \\\n\t-R \${genome} \\\n\t-T HaplotypeCaller -ERC GVCF \\\n\t-I ../01.Mapping/${sample}.dd.bam -o ${sample}.gatk.g.vcf \\\n\t\&> ${sample}.HaplotypeCaller.log\n";

close OUT;
