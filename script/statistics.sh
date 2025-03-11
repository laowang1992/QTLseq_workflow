
cd ${work_dir}/refseq

# Statistics
# fastQC
cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} `cut -f3,4 ../samples.txt`
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
# data stat
cd ${work_dir}/00.data/
Rscript dataStat.R -d 01.clean_data -o data_stat

cd ${work_dir}/01.Mapping
# coverage rate and depth
for sample in $(cut -f1 ${sampleInfo})
do
	pandepth -i $sample.dd.bam -w 100000 -t $thread -o $sample
done
Rscript CoverageStatistic.R --sampleInfo ../00.data/samples.txt --chrInfo ../refseq/chrom.txt --chrLen ../refseq/ref.len
# align rate
if [ $aligner = bowtie2 ];then
	echo -e "Sample,Total reads,Mapped reads,Mapped rate,Uniquely mapped reads,Uniquely mapped rate" > align_stat.csv
	for i in $(cut -f1 ${sampleInfo}); do perl alignStat.pl $i; done >> align_stat.csv
else
	echo "Sample,Unmapped reads,Uniquely mapped reads" > align_stat.csv
	for i in $(cut -f1 ${sampleInfo})
	do
		unmapped_reads=$(samtools idxstats $i.sort.bam | awk -F'\t' '{sum += $4} END {print sum}')
		unique_reads=$(samtools view -q 60 -c $i.sort.bam)
		echo "$i,$unmapped_reads,$unique_reads" >> align_stat.csv
	done
fi

cd ${work_dir}/02.SNP_indel
# VariantQC
java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf.gz --maxContigs `cat $work_dir/refseq/chrom.txt | wc -l` --threads ${thread}
