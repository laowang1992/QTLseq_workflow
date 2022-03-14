#PBS -N xx
#PBS -l nods=1:ppn=1
#PBS -l mem=50gb
#PBS -q batch

cd ${work_dir}/02.SNP_indel

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
Rscript covStat.R ${genome}.fai
echo -e "Sample,Total read,Mapping read,Mapping rate,Unique mapping read,Unique mapping rate" > align_stat.csv
for i in $(cut -f1 ${sample}); do perl alignStat.pl $i; done >> align_stat.csv
