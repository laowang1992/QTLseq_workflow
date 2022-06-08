####################################################
# 数据准备
## 建基因组索引
cd ${work_dir}/refseq
# annovar建库
gffread ${gff} -T -o ${gtf}
gtfToGenePred  -genePredExt ${gtf} genome_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile ${genome} genome_refGene.txt --out genome_refGeneMrna.fa
#
samtools faidx ${genome}
java -jar ${picard} CreateSequenceDictionary R=${genome} O=${genome/fa/dict}
bowtie2-build ${genome} ${index}

awk '{print $1"\t"$2}' ${genome}.fai > ref.len
cp ref.len ${work_dir}/04.Analysis/
bedtools makewindows -w 200000 -s 100000 -g ref.len > ref.window.bed
####################################################


cat ${sampleInfo} | while read sample group fq1 fq2
do
# 质控 过滤
cd ${work_dir}/00.data/01.clean_data

fastp -i ${fq1} -o ./${sample}_1.clean.fastq.gz \
      -I ${fq2} -O ./${sample}_2.clean.fastq.gz \
      --json=./${sample}.json --html=${sample}.html --report_title="${sample} fastp report" \
      --thread=${thread} --length_required 50
# 比对
cd ${work_dir}/01.Mapping
bowtie2 --rg-id ${sample} --rg "PL:ILLUMINA" --rg "SM:${sample}" \
        -x ${index} \
        -1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
        -2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
        -p ${thread} \
        -S ${sample}.sam \
        2> ${sample}.log
# 排序
if [ $sort = sbb ];then
	sambamba view --format=bam --with-header --sam-input --nthreads=${thread} --output-filename ${sample}.bam ${sample}.sam
	sambamba sort --nthreads=${thread} --memory-limit=20GB --tmpdir=./ --out=${sample}.sort.bam ${sample}.bam
	#rm -rf sambamba-pid*
elif [ $sort = sts ];then
	samtools sort --threads ${thread} --output-fmt BAM -o ${sample}.sort.bam ${sample}.sam
fi
rm ${sample}.sam

done

cd ${work_dir}/01.Mapping
# 去重复
# sambamba markdup在batch队列使用28线程时会休眠，卡在这里不动
if [ $rmdup = picard ];then
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		java -Xmx20g -jar ${picard} \
			MarkDuplicates I=%.sort.bam O=%.dd.bam \
			CREATE_INDEX=true REMOVE_DUPLICATES=true \
			M=%.dd.metics
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		samtools index %.dd.bam
elif [ $rmdup = sbb ];then
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		sambamba markdup --remove-duplicates --nthreads=1 --tmpdir=./ %.sort.bam %.dd.bam
		#rm -rf sambamba-pid*
fi
rm *sort.bam *sort.bam.bai

awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        genomeCoverageBed -ibam %.dd.bam -bga -g ${genome} "|" \
        grep -w "0$" ">" %.0cov.bedgraph

awk '{print $1}' ${sampleInfo} | \
	parallel -j ${thread} -I% --max-args 1 \
	bedtools bamtobed -i %.dd.bam ">" %.dd.bed

cd ${work_dir}/02.SNP_indel

# GATK HaplotypeCaller多线程
for sample in $(awk '{print $1}' ${sampleInfo})
do
	java -Xmx20g -jar ${gatk} \
		-R ${genome} \
		-T HaplotypeCaller -ERC GVCF -nct ${thread} \
		-I ../01.Mapping/${sample}.dd.bam -o ${sample}.gatk.g.vcf.gz \
		&> ${sample}.HaplotypeCaller.log
done
#awk '{print $1}' ${sampleInfo} | \
#        parallel -j ${thread} -I% --max-args 1 \
#        java -Xmx20g -jar ${gatk} \
#        -R ${genome} \
#        -T HaplotypeCaller -ERC GVCF \
#        -I ../01.Mapping/%.dd.bam -o %.gatk.g.vcf.gz \
#        "&>" %.HaplotypeCaller.log

ls *g.vcf.gz > GVCFs.list
#cut -f1 ../00.data/samples.txt | sed 's/$/.gatk.g.vcf.gz/' > GVCFs.list

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T CombineGVCFs \
     -V GVCFs.list \
     -o ${filename}.gatk.g.vcf.gz

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T GenotypeGVCFs \
     -V ${filename}.gatk.g.vcf.gz -o ${filename}.gatk.vcf.gz

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T VariantFiltration \
     -o ${filename}.flt.vcf.gz -V ${filename}.gatk.vcf.gz \
     --filterName FilterQual --filterExpression "QUAL<30.0" \
     --filterName FilterQD --filterExpression "QD<13.0" \
     --filterName FilterMQ --filterExpression "MQ<20.0" \
     --filterName FilterFS --filterExpression "FS>20.0" \
     --filterName FilterMQRankSum --filterExpression "MQRankSum<-3.0" \
     --filterName FilterReadPosRankSum --filterExpression "ReadPosRankSum<-3.0" \
     2> ${filename}.flt.log

#grep -vP "\tFilter" ${filename}.flt.vcf > ${filename}.filter.vcf
# 改成bcftools
bcftools view -f PASS -o ${filename}.filter.vcf.gz -O z ${filename}.flt.vcf.gz
bcftools index -t ${filename}.filter.vcf.gz 

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType SNP \
     -V ${filename}.filter.vcf.gz -o ${filename}.filter.SNPs.vcf.gz

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType INDEL \
     -V ${filename}.filter.vcf.gz -o ${filename}.filter.INDELs.vcf.gz
# 保留双等位基因位点
bcftools view -o ${filename}.biallelic.SNPs.vcf.gz -O z -m 2 -M 2 ${filename}.filter.SNPs.vcf.gz

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T VariantsToTable \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
     -V ${filename}.biallelic.SNPs.vcf.gz -o ../04.Analysis/${filename}.biallelic.SNPs.table

#grep -v "##" ${filename}.biallelic.SNPs.vcf.gz | sed 's/^#CHROM/CHROM/' > ../04.Analysis/${filename}.biallelic.SNPs.txt

java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf

${work_dir}/04.Analysis
gzip ${filename}.biallelic.SNPs.table
#gzip ${filename}.biallelic.SNPs.txt

# annotation
cd ${work_dir}/03.Annotation
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.SNPs.vcf.gz --outfile ./${filename}.filter.SNPs.annovar.input
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.INDELs.vcf.gz --outfile ./${filename}.filter.INDELs.annovar.input
annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${filename}.filter.SNPs.anno --exonsort ./${filename}.filter.SNPs.annovar.input ../refseq
annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${filename}.filter.INDELs.anno --exonsort ./${filename}.filter.INDELs.annovar.input ../refseq

Rscript AnnoStat.R --snpvar ${filename}.filter.SNPs.anno.variant_function \
	--indelvar ${filename}.filter.INDELs.anno.variant_function \
	--snpex ${filename}.filter.SNPs.anno.exonic_variant_function \
	--indelex ${filename}.filter.INDELs.anno.exonic_variant_function

# 
cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript dataStat.R

cd ${work_dir}/01.Mapping

for i in $(cut -f1 ${sampleInfo})
do
bedtools coverage -b $i.dd.bed -a ${work_dir}/refseq/ref.window.bed -mean | \
	awk '{print $1"\t"$2+1"\t"$3"\t"$4}' > $i.dd.window.depth
rm $i.dd.bed
done

Rscript covStat.R ${genome}.fai
echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
for i in $(cut -f1 ${sampleInfo}); do perl alignStat.pl $i; done >> align_stat.csv

