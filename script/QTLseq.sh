#!/bin/bash
cd ${work_dir}/refseq
cp ref.len ${work_dir}/03.Analysis/
cp chrom.txt ${work_dir}/03.Analysis/

cd ${work_dir}/02.SNP_indel
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
# QTL-seq·ÖÎö
Rscript QTLseq.R -i ${filename}.filter.SNPs.table.gz -p Parameter.csv
