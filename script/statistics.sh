
cd ${work_dir}/refseq

# Statistics
# fastQC
cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} `cut -f3,4 ../samples.txt`
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
# data stat
cd ${work_dir}/00.data/
Rscript dataStat.R

cd ${work_dir}/01.Mapping
# coverage rate and depth
for sample in $(cut -f1 ${sampleInfo})
do
	pandepth -i $sample.dd.bam -w 100000 -t $thread -o $sample
done

# align rate
if [ $aligner = bowtie2 ];then
	echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
	for i in $(cut -f1 ${sampleInfo}); do perl alignStat.pl $i; done >> align_stat.csv
else
	echo "Sample,Mapping read,Unique mapping read" > align_stat.csv
	for i in $(cut -f1 ${sampleInfo})
	do
		mapping_reads=$(samtools view -c $i.sort.bam)
		unique_reads=$(samtools view -q 60 -c $i.sort.bam)
		echo "$i,$mapping_reads,$unique_reads" >> align_stat.csv
	done
fi

cd ${work_dir}/02.SNP_indel
# VariantQC
java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf.gz --maxContigs 19 --threads ${thread}

