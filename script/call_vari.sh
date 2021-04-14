####################################################
# 数据准备

# 建基因组索引
cd ${work_dir}/refseq
samtools faidx ${genome}
java -jar ${picard} CreateSequenceDictionary R=${genome} O=${genome/fa/dict}
bowtie2-build ${genome} ${index}
####################################################

IFS_OLD=$IFS
IFS=$'\n'

for i in $(cat ${sample})
do

IFS=$'\t'
i=($i)
IFS=$IFS_OLD

# 质控 过滤
cd ${work_dir}/00.data/01.clean_data

fastp -i ${i[2]} -o ./${i[0]}_1.clean.fastq.gz \
      -I ${i[3]} -O ./${i[0]}_2.clean.fastq.gz \
      --json=./${i[0]}.json --html=${i[0]}.html --report_title="${i[0]} fastp report" \
      --thread=8 --length_required 100

# 比对
cd ${work_dir}/01.Mapping

bowtie2 --rg-id ${i[0]} --rg "PL:ILLUMINA" --rg "SM:${i[0]}" \
        -x ${index} \
        -1 ../00.data/01.clean_data/${i[0]}_1.clean.fastq.gz \
        -2 ../00.data/01.clean_data/${i[0]}_2.clean.fastq.gz \
        -p ${thread} \
        -S ${i[0]}.sam \
        2> ${i[0]}.log

samtools sort -@ ${thread} -O BAM -o ${i[0]}.sort.bam ${i[0]}.sam
rm ${i[0]}.sam

done


awk '{print $1}' ${sample} | \
        parallel -j ${thread} -I% --max-args 1 \
        java -Xmx20g -jar ${picard} \
        MarkDuplicates I=%.sort.bam O=%.dd.bam \
        CREATE_INDEX=true REMOVE_DUPLICATES=true \
        M=%.dd.metics

awk '{print $1}' ${sample} | \
        parallel -j ${thread} -I% --max-args 1 \
        samtools index %.dd.bam
#genomeCoverageBed -ibam test.sort.bam -bga -g ../refseq/TAIR10_genome.fa > test.cov.bedgraph
awk '{print $1}' ${sample} | \
        parallel -j ${thread} -I% --max-args 1 \
        genomeCoverageBed -ibam %.dd.bam -bga -g ${genome} "|" \
        grep -w "0$" ">" %.0cov.bedgraph

cd ${work_dir}/02.SNP_indel

awk '{print $1}' ${sample} | \
        parallel -j ${thread} -I% --max-args 1 \
        java -Xmx20g -jar ${gatk} \
        -R ${genome} \
        -T HaplotypeCaller -ERC GVCF \
        -I ../01.Mapping/%.dd.bam -o %.gatk.g.vcf \
        "&>" %.HaplotypeCaller.log

#ls *g.vcf > GVCFs.list
cut -f1 ../00.data/samples.txt | sed 's/$/.g.vcf/' > GVCFs.list

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
     -V ${filename}.filter.SNPs.vcf -o ../03.Analysis/${filename}.filter.SNPs.table

grep -v "##" ${filename}.filter.SNPs.vcf | sed 's/^#CHROM/CHROM/' > ../03.Analysis/${filename}.filter.SNPs.txt

java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf



# 
cd ${work_dir}/00.data/00.raw_data
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript stat.R

cd ${work_dir}/01.Mapping
Rscript covStat.R ${genome}.fai
echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
for i in $(cut -f2 ${sample}); do perl alignStat.pl $i; done >> align_stat.csv

