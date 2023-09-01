#!/bin/bash
#. ../script/.conf
cd $work_dir/05.Primer

# lower case to UPPER case
seqkit seq -u $genome > genome.fa
# 建哈希数据库
famap -t N -b genome.famap genome.fa
fahash -b genome.hash -w 12 -f 3 genome.famap

# 
csvtk cut -f outPrefix,highParent,lowParent -D $'\t' $work_dir/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }' | while read outPrefix highParent lowParent
do
	mkdir -p $outPrefix
	cd $outPrefix
	
	ln -s ../genome.fa ./
	ln -s ../genome.famap ./
	ln -s ../genome.hash ./
	
	bcftools view -R $work_dir/04.Annotation/$i/$i.INDELs.pos.txt \
		-s $highParent,lowParent --force-samples \
		-Oz -o ./$i.INDELs.vcf.gz --threads \
		${thread} ${work_dir}/02.SNP_indel/$filename.filter.INDELs.vcf.gz
	
	bcftools annotate -I '%CHROM\_%POS' -x INFO,^FORMAT/GT ./$i.INDELs.vcf.gz | \
		grep -v "##" | \
		sed 's/^#CHROM/CHROM/' > $filename.indels.txt
	
	# 提取indels上下游各250bp序列
	perl extr_seq_for_indel_pri.pl --genome genome.fa --input $filename.indels.txt --output $filename.indelWithSEQ.txt
	# 整理成p3in格式
	perl indel_p3in.pl --input $filename.indelWithSEQ.txt --output $filename.p3in \
		--config $primerConfig \
		--minProdLen $minProdLen \
		--maxProdLen $maxProdLen \
		--primerMinSize $primerMinSize \
		--primerOptSize $primerOptSize \
		--primerMaxSize $primerMaxSize \
		--primerMinGC $primerMinGC \
		--primerMaxGC $primerMaxGC \
		--primerMinTM $primerMinTM \
		--primerOptTM $primerOptTM \
		--primerMaxTM $primerMaxTM \
		--primerNumReturn $primerNumReturn
	
	# 设计引物，每个位置设计3个
	primer3_core -strict_tags ${filename}.p3in > ${filename}.p3out
	# e-PCR
	perl epcr.pl ${filename}.p3out ${filename}.primer.txt
	# mkxls.pl已经整合到epcr.pl这一步中
	#perl mkxls.pl --input ${filename}.p3out --output ${filename}.primer.txt
	# merge
	Rscript merge.R ${filename}.indels.filter.txt ${filename}.primer.txt
	
	# 删除中间文件
	rm genome.fa genome.famap genome.hash ${filename}.p3* ${filename}.indelWithSEQ.txt ${filename}.indels.txt ${filename}.primer.txt
	cd ..
done
