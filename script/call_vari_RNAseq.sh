ulimit -n 1024000
# 建索引
cd ${work_dir}/refseq
samtools faidx ${genome}
awk '{print $1"\t"$2}' ${genome}.fai > ref.len
java -jar ${picard} CreateSequenceDictionary --REFERENCE ${genome} --OUTPUT ${genome%.*}.dict
gffread ${gff} -T -o ${gtf}
mkdir -p ${index}STARindex
STAR --runThreadN ${thread} \
	--runMode genomeGenerate \
	--genomeDir ${index}STARindex \
	--genomeFastaFiles ${genome} \
	--sjdbGTFfile ${gtf} \
	--sjdbOverhang ${sjdbOverhang}

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
	mkdir -p ${sample}
	STAR --twopassMode Basic \
		--runThreadN $thread --genomeDir ${index}STARindex \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA \
		--outSAMmapqUnique 60 \
		--outSAMmultNmax 1 \
		--outFilterMismatchNoverReadLmax 0.04 --outSJfilterReads Unique \
		--outFileNamePrefix ./$sample/$sample \
		--readFilesCommand gunzip -c \
		--readFilesIn $work_dir/00.data/01.clean_data/${sample}_1.clean.fastq.gz $work_dir/00.data/01.clean_data/${sample}_2.clean.fastq.gz
done

cd ${work_dir}/01.Mapping
# 去重复
# sambamba markdup在batch队列使用28线程时会休眠，卡在这里不动
if [ $rmdup = picard ];then
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		java -Xmx20g -jar ${picard} \
			MarkDuplicates \
			I=./%/"%"Aligned.sortedByCoord.out.bam \
			O=./%/"%".dedup.bam \
			CREATE_INDEX=true \
			REMOVE_DUPLICATES=true \
			VALIDATION_STRINGENCY=SILENT \
			M=./%/"%".dedup.metrics
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		samtools index ./%/"%".dedup.bam
elif [ $rmdup = sbb ];then
	awk '{print $1}' ${sampleInfo} | \
		parallel -j ${thread} -I% --max-args 1 \
		sambamba markdup --remove-duplicates --nthreads=1 --tmpdir=./% ./%/"%"Aligned.sortedByCoord.out.bam ./%/"%".dedup.bam
		#rm -rf sambamba-pid*
fi

# reads剪切
awk '{print $1}' ${sampleInfo} | \
	parallel -j ${thread} -I% --max-args 1 \
	java -jar ${gatk} -T SplitNCigarReads \
		-R ${genome} \
		-I ./%/"%".dedup.bam \
		-o ./%/"%".dedup_split.bam \
		-U ALLOW_N_CIGAR_READS

# 为$sample.dedup_split.bam建连接，方便整合到后续分析（如：pandepth分析）
for sample in $(awk '{print $1}' ${sampleInfo})
do
	ln -s ./$sample/$sample.dedup_split.bam ./$sample.dd.bam
	ln -s ./$sample/$sample.dedup_split.bai ./$sample.dd.bam.bai
done

# 变异检测
cd ${work_dir}/02.SNP_indel
# GATK HaplotypeCaller多线程
for sample in $(awk '{print $1}' ${sampleInfo})
do
	java -jar ${gatk} -T HaplotypeCaller \
		-R ${genome} \
		-ERC GVCF -nct ${thread} \
		-I ../01.Mapping/$sample/$sample.dedup_split.bam \
		-dontUseSoftClippedBases \
		-stand_call_conf 20.0 \
		-o ./$sample.gatk.g.vcf.gz \
		&> $sample.HaplotypeCaller.log
done

ls *g.vcf.gz > GVCFs.list

# 合并gvcf
java -jar ${gatk} \
     -R ${genome} -T CombineGVCFs \
     -V GVCFs.list \
     -o ${filename}.gatk.g.vcf.gz \
     &> ${filename}.CombineGVCFs.log

# gvcf转vcf
java -jar ${gatk} \
     -R ${genome} -T GenotypeGVCFs \
     -V ${filename}.gatk.g.vcf.gz \
     -o ${filename}.gatk.vcf.gz \
     &> ${filename}.GenotypeGVCFs.log

# 过滤低质量
java -jar ${gatk} \
     -T VariantFiltration \
     -R ${genome} \
     -V ${filename}.gatk.vcf.gz \
     -window 35 -cluster 3 \
     -filterName FilterFS -filter "FS > 30.0" \
     -filterName FilterQD -filter "QD < 2.0" \
     -o ./${filename}.flt.vcf.gz \
     &> ${filename}.VariantFiltration.log

# 改成bcftools，保留bi-allele
bcftools view -f PASS -o ${filename}.filter.SNPs.vcf.gz -O z -m 2 -M 2 -v snps --threads ${thread} ${filename}.flt.vcf.gz
bcftools view -f PASS -o ${filename}.filter.INDELs.vcf.gz -O z -m 2 -M 2 -v indels --threads ${thread} ${filename}.flt.vcf.gz
bcftools index --tbi --threads $thread ${filename}.filter.SNPs.vcf.gz
bcftools index --tbi --threads $thread ${filename}.filter.INDELs.vcf.gz
