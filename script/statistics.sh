
cd ${work_dir}/refseq
bedtools makewindows -w 200000 -s 100000 -g ref.len > ref.window.bed

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
# coverage rate
awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        genomeCoverageBed -ibam %.dd.bam -bga -g ${genome} "|" \
        grep -w "0$" ">" %.0cov.bedgraph
Rscript covStat.R ${genome}.fai

# depth
awk '{print $1}' ${sampleInfo} | \
	parallel -j ${thread} -I% --max-args 1 \
	bedtools bamtobed -i %.dd.bam ">" %.dd.bed
for i in $(cut -f1 ${sampleInfo})
do
bedtools coverage -b $i.dd.bed -a ${work_dir}/refseq/ref.window.bed -mean | \
	awk '{print $1"\t"$2+1"\t"$3"\t"$4}' > $i.dd.window.depth
rm $i.dd.bed
done
Rscript CoverageDepth.R --sampleInfo ${sampleInfo} --chrInfo ../refseq/chrom.txt --chrLen ../refseq/ref.len

# align rate
if [ $aligner = bowtie2 ];then
	echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
	for i in $(cut -f1 ${sampleInfo}); do perl alignStat.pl $i; done >> align_stat.csv
fi

cd ${work_dir}/02.SNP_indel
# VariantQC
java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf.gz --maxContigs 19 --threads ${thread}

