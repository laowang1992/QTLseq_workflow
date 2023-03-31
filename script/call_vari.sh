####################################################
## 数据准备
## 建基因组索引
cd ${work_dir}/refseq
# annovar建库
gffread ${gff} -T -o ${gtf}
gtfToGenePred  -genePredExt ${gtf} genome_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile ${genome} genome_refGene.txt --out genome_refGeneMrna.fa
#
samtools faidx ${genome}
java -jar ${picard} CreateSequenceDictionary R=${genome} O=${genome/fa/dict}
if [ $aligner = bowtie2 ];then
	bowtie2-build ${genome} ${index}
elif [ $aligner = mem ];then
	bwa index ${genome} -p ${index}
elif [ $aligner = mem2 ];then
	bwa-mem2 index ${genome} -p ${index}
fi

awk '{print $1"\t"$2}' ${genome}.fai > ref.len
bedtools makewindows -w 200000 -s 100000 -g ref.len > ref.window.bed
cp ref.len ${work_dir}/03.Analysis/
cp chrom.txt ${work_dir}/03.Analysis/
cp ref.len ${work_dir}/04.Annotation/
cp chrom.txt ${work_dir}/04.Annotation/
####################################################

####################################################
## call variation
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
if [ $aligner = bowtie2 ];then
	bowtie2 --rg-id ${sample} --rg "PL:ILLUMINA" --rg "SM:${sample}" \
		-x ${index} \
		-1 ../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		-2 ../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		-p ${thread} \
		-S ${sample}.sam \
		2> ${sample}.log
elif [ $aligner = mem ];then
	bwa mem -t ${thread} -M -k 32 \
		-R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" \
		${index} \
		../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		-o ${sample}.sam \
		2> ${sample}.log
elif [ $aligner = mem2 ];then
	bwa-mem2 mem -t ${thread} -M -k 32 \
		-R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" \
		${index} \
		../00.data/01.clean_data/${sample}_1.clean.fastq.gz \
		../00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		-o ${sample}.sam \
		2> ${sample}.log
fi
# 排序
if [ $sort = sbb ];then
	sambamba view --format=bam --with-header --sam-input --nthreads=${thread} --output-filename ${sample}.bam ${sample}.sam
	sambamba sort --nthreads=${thread} --memory-limit=20GB --tmpdir=./ --out=${sample}.sort.bam ${sample}.bam
	rm ${sample}.bam
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

ls *g.vcf.gz > GVCFs.list
#cut -f1 ../00.data/samples.txt | sed 's/$/.gatk.g.vcf.gz/' > GVCFs.list

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T CombineGVCFs \
     -V GVCFs.list \
     -o ${filename}.gatk.g.vcf.gz \
     &> ${filename}.CombineGVCFs.log

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T GenotypeGVCFs \
     -V ${filename}.gatk.g.vcf.gz \
     -o ${filename}.gatk.vcf.gz \
     &> ${filename}.GenotypeGVCFs.log

# filter
#java -Xmx30g -jar ${gatk} \
#     -R ${genome} -T VariantFiltration \
#     -o ${filename}.flt.vcf.gz -V ${filename}.gatk.vcf.gz \
#     --filterName FilterQual --filterExpression "QUAL<30.0" \
#     --filterName FilterQD --filterExpression "QD<13.0" \
#     --filterName FilterMQ --filterExpression "MQ<20.0" \
#     --filterName FilterFS --filterExpression "FS>20.0" \
#     --filterName FilterMQRankSum --filterExpression "MQRankSum<-3.0" \
#     --filterName FilterReadPosRankSum --filterExpression "ReadPosRankSum<-3.0" \
#     &> ${filename}.VariantFiltration.log

# 20220813，修改过滤参数，以前的有点太严格，在202208_BnQTLseq项目中根据VariantQC结果，FS、MQRankSum和QD会过滤掉一半位点
java -Xmx30g -jar ${gatk} \
     -R ${genome} -T VariantFiltration \
     -o ${filename}.flt.vcf.gz -V ${filename}.gatk.vcf.gz \
     --filterName FilterQual --filterExpression "QUAL<30.0" \
     --filterName FilterQD --filterExpression "QD<2.0" \
     --filterName FilterMQ --filterExpression "MQ<20.0" \
     --filterName FilterFS --filterExpression "FS>60.0" \
     --filterName FilterMQRankSum --filterExpression "MQRankSum<-12.5" \
     --filterName FilterReadPosRankSum --filterExpression "ReadPosRankSum<-3.0" \
     &> ${filename}.VariantFiltration.log

# 改成bcftools，保留bi-allele
bcftools view -f PASS -o ${filename}.filter.SNPs.vcf.gz -O z -m 2 -M 2 -v snps --threads ${thread} ${filename}.flt.vcf.gz
bcftools view -f PASS -o ${filename}.filter.INDELs.vcf.gz -O z -m 2 -M 2 -v indels --threads ${thread} ${filename}.flt.vcf.gz
bcftools index -t ${filename}.filter.SNPs.vcf.gz
bcftools index -t ${filename}.filter.INDELs.vcf.gz

java -Xmx30g -jar ${gatk} \
     -R ${genome} -T VariantsToTable \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF GQ \
     -V ${filename}.filter.SNPs.vcf.gz -o ../03.Analysis/${filename}.filter.SNPs.table
java -Xmx30g -jar ${gatk} \
     -R ${genome} -T VariantsToTable \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF GQ \
     -V ${filename}.filter.INDELs.vcf.gz -o ../03.Analysis/${filename}.filter.INDELs.table

cd ${work_dir}/03.Analysis
gzip ${filename}.filter.SNPs.table
gzip ${filename}.filter.INDELs.table
# QTL-seq分析
Rscript QTLseq.R -i ${filename}.filter.SNPs.table.gz -p Parameter.csv

####################################################

####################################################
# annotation
cd ${work_dir}/04.Annotation
# 因为运行R脚本以后一些输出文件会被覆盖，所以先处理完SNP再用INDEL跑一遍
Rscript select_variation.R -i ../03.Analysis/${filename}.filter.SNPs.table.gz -p ../03.Analysis/Parameter.csv
for i in `cut -f1 -d',' ${work_dir}/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }'`
do
	cd ${work_dir}/04.Annotation/${i}
	cut -f1,2 ${i}.Depth_information.txt | awk '{if(NR!=1) print }' > ${i}.pos.txt
	bcftools view -R ${i}.pos.txt -Oz -o ./${i}.SNPs.vcf.gz --threads ${thread} ${work_dir}/02.SNP_indel/${filename}.filter.SNPs.vcf.gz
	convert2annovar.pl --format vcf4old ./${i}.SNPs.vcf.gz --outfile ./${i}.SNPs.annovar.input
	annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${i}.SNPs.anno --exonsort ./${i}.SNPs.annovar.input ${work_dir}/refseq
	rm ${i}.Depth_information.csv ${i}.Depth_information.txt ${i}.SNPs.vcf.gz ${i}.parameter.txt ${i}.pos.txt
done
cd ${work_dir}/04.Annotation
Rscript select_variation.R -i ../03.Analysis/${filename}.filter.INDELs.table.gz -p ../03.Analysis/Parameter.csv
for i in `cut -f1 -d',' ${work_dir}/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }'`
do
	cd ${work_dir}/04.Annotation/${i}
	cut -f1,2 ${i}.Depth_information.txt | awk '{if(NR!=1) print }' > ${i}.pos.txt
	bcftools view -R ${i}.pos.txt -Oz -o ./${i}.INDELs.vcf.gz --threads ${thread} ${work_dir}/02.SNP_indel/${filename}.filter.INDELs.vcf.gz
	convert2annovar.pl --format vcf4old ./${i}.INDELs.vcf.gz --outfile ./${i}.INDELs.annovar.input
	annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${i}.INDELs.anno --exonsort ./${i}.INDELs.annovar.input ${work_dir}/refseq
	rm ${i}.Depth_information.csv ${i}.Depth_information.txt ${i}.INDELs.vcf.gz ${i}.parameter.txt ${i}.pos.txt
done
cd ${work_dir}/04.Annotation
for i in `cut -f1 -d',' ${work_dir}/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }'`
do
	cd ${work_dir}/04.Annotation/${i}
	Rscript ../AnnoStat.R --snpvar ${i}.SNPs.anno.variant_function \
		--indelvar ${i}.INDELs.anno.variant_function \
		--snpex ${i}.SNPs.anno.exonic_variant_function \
		--indelex ${i}.INDELs.anno.exonic_variant_function
done

####################################################

####################################################
##  Statistics
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

