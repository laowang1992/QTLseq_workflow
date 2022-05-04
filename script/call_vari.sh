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

samtools sort -@ ${thread} -O BAM -o ${sample}.sort.bam ${sample}.sam
# 使用sambamba代替samtools排序和Picard去重复
# sambamba markdup在batch队列使用28线程时会休眠，卡在这里不动，还是换回samtools和Picard吧

#sambamba view --format=bam --with-header --sam-input --nthreads=${thread} --output-filename ${i[0]}.bam ${i[0]}.sam
#sambamba sort -t ${thread} -m 20GB --tmpdir=./ -o ${i[0]}.sort.bam ${i[0]}.bam
#sambamba markdup -r -t ${thread} --tmpdir=./ ${i[0]}.sort.bam ${i[0]}.dd.bam
rm ${sample}.sam

done


awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        java -Xmx20g -jar ${picard} \
        MarkDuplicates I=%.sort.bam O=%.dd.bam \
        CREATE_INDEX=true REMOVE_DUPLICATES=true \
        M=%.dd.metics
rm *sort.bam

awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        samtools index %.dd.bam
#genomeCoverageBed -ibam test.sort.bam -bga -g ../refseq/TAIR10_genome.fa > test.cov.bedgraph
awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        genomeCoverageBed -ibam %.dd.bam -bga -g ${genome} "|" \
        grep -w "0$" ">" %.0cov.bedgraph

awk '{print $1}' ${sampleInfo} | \
	parallel -j ${thread} -I% --max-args 1 \
	bedtools bamtobed -i %.dd.bam ">" %.dd.bed


cd ${work_dir}/02.SNP_indel

awk '{print $1}' ${sampleInfo} | \
        parallel -j ${thread} -I% --max-args 1 \
        java -Xmx20g -jar ${gatk} \
        -R ${genome} \
        -T HaplotypeCaller -ERC GVCF \
        -I ../01.Mapping/%.dd.bam -o %.gatk.g.vcf \
        "&>" %.HaplotypeCaller.log

ls *g.vcf > GVCFs.list
#cut -f1 ../00.data/samples.txt | sed 's/$/.gatk.g.vcf/' > GVCFs.list

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T CombineGVCFs \
     -V GVCFs.list \
     -o ${filename}.gatk.g.vcf

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T GenotypeGVCFs \
     -V ${filename}.gatk.g.vcf -o ${filename}.gatk.vcf

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T VariantFiltration \
     -o ${filename}.flt.vcf -V ${filename}.gatk.vcf \
     --filterName FilterQual --filterExpression "QUAL<30.0" \
     --filterName FilterQD --filterExpression "QD<13.0" \
     --filterName FilterMQ --filterExpression "MQ<20.0" \
     --filterName FilterFS --filterExpression "FS>20.0" \
     --filterName FilterMQRankSum --filterExpression "MQRankSum<-3.0" \
     --filterName FilterReadPosRankSum --filterExpression "ReadPosRankSum<-3.0" \
     2> ${filename}.flt.log

grep -vP "\tFilter" ${filename}.flt.vcf > ${filename}.filter.vcf

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType SNP \
     -V ${filename}.filter.vcf -o ${filename}.filter.SNPs.vcf

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType INDEL \
     -V ${filename}.filter.vcf -o ${filename}.filter.INDELs.vcf

java -Xmx50g -jar ${gatk} \
     -R ${genome} -T VariantsToTable \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
     -V ${filename}.filter.SNPs.vcf -o ../04.Analysis/${filename}.filter.SNPs.table

grep -v "##" ${filename}.filter.SNPs.vcf | sed 's/^#CHROM/CHROM/' > ../04.Analysis/${filename}.filter.SNPs.txt

java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf

# annotation
cd ${work_dir}/03.Annotation
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.SNPs.vcf --outfile ./${filename}.filter.SNPs.annovar.input
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.INDELs.vcf --outfile ./${filename}.filter.INDELs.annovar.input
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

