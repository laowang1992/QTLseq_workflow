#PBS -N xx
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb

. ./../script/.conf_bnHPC

cd ${work_dir}/01.Mapping
java -Xmx10g -jar ${picard} \
	MarkDuplicates I=HHHH.sort.bam O=HHHH.dd.bam \
	CREATE_INDEX=true REMOVE_DUPLICATES=true \
	M=HHHH.dd.metics
samtools index HHHH.dd.bam
genomeCoverageBed -ibam HHHH.dd.bam -bga -g ${genome} | grep -w "0$" > HHHH.0cov.bedgraph

cd ${work_dir}/02.SNP_indel
java -Xmx10g -jar ${gatk} \
	-R ${genome} \
	-T HaplotypeCaller -ERC GVCF \
	-I ../01.Mapping/HHHH.dd.bam -o HHHH.gatk.g.vcf \
	&> HHHH.HaplotypeCaller.log
