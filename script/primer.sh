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
	
	bcftools view -R $work_dir/04.Annotation/$outPrefix/$outPrefix.INDELs.pos.txt \
		-s $highParent,$lowParent --force-samples \
		-Oz -o ./$outPrefix.INDELs.vcf.gz --threads \
		${thread} ${work_dir}/02.SNP_indel/$filename.filter.INDELs.vcf.gz
	
	bcftools annotate -I '%CHROM\_%POS' -x INFO,^FORMAT/GT ./$outPrefix.INDELs.vcf.gz | \
		grep -v "##" | \
		sed 's/^#CHROM/CHROM/' > $outPrefix.indels.txt
	
	# 提取indels上下游各250bp序列
	perl ../extr_seq_for_indel_pri.pl --genome genome.fa --input $outPrefix.indels.txt --output $outPrefix.indelWithSEQ.txt
	# 整理成p3in格式
	perl ../indel_p3in.pl --input $outPrefix.indelWithSEQ.txt --output $outPrefix.p3in \
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
	primer3_core -strict_tags $outPrefix.p3in > $outPrefix.p3out
	# e-PCR
	perl ../epcr.pl $outPrefix.p3out $outPrefix.primer.txt
	# merge
	Rscript ../merge.R $outPrefix.indels.filter.txt $outPrefix.primer.txt
	
	# 删除中间文件
	#rm genome.fa genome.famap genome.hash $outPrefix.p3* $outPrefix.indelWithSEQ.txt $outPrefix.indels.txt $outPrefix.primer.txt
	cd ..
done
