ulimit -n 1024000
# 建索引
cd ${work_dir}/refseq
samtools faidx ${genome}
awk '{print $1"\t"$2}' ${genome}.fai > ref.len
java -jar ${picard} CreateSequenceDictionary --REFERENCE ${genome} --OUTPUT ${genome%.*}.dict
if [ $aligner = bowtie2 ];then
	bowtie2-build ${genome} ${index}
elif [ $aligner = mem ];then
	bwa index ${genome} -p ${index}
elif [ $aligner = mem2 ];then
	bwa-mem2 index ${genome} -p ${index}
fi

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
	sambamba view --format=bam --with-header --sam-input --nthreads=${thread} --output-filename ${sample}.bam ${sample}.sam && rm ${sample}.sam || echo "sam转bam命令执行失败"
	sambamba sort --nthreads=${thread} --memory-limit=20GB --tmpdir=./ --out=${sample}.sort.bam ${sample}.bam && rm ${sample}.bam || echo "sort bam命令执行失败"
	#rm -rf sambamba-pid*
elif [ $sort = sts ];then
	samtools sort --threads ${thread} --output-fmt BAM -o ${sample}.sort.bam ${sample}.sam && rm ${sample}.sam || echo "sort bam命令执行失败"
fi
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
bcftools index --tbi --threads $thread ${filename}.filter.SNPs.vcf.gz
bcftools index --tbi --threads $thread ${filename}.filter.INDELs.vcf.gz
