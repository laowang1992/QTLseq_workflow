#!/bin/bash
cd ${work_dir}/refseq
# annovar建库
gffread ${gff} -T -o ${gtf}
# 如果gff文件没有gene id信息，需要执行下面的命令，修改gtf文件
if false; then
awk -F'\t' 'BEGIN{OFS="\t"} {
    if($0 ~ /transcript_id/){
        match($9, /transcript_id "([^"]+)"/, arr);
        if(arr[1]!=""){
            gid=arr[1]"G";
            $9="gene_id \""gid"\"; gene_name \""gid"\"; " $9;
        }
    }
    print
}' ${gtf} > tmp.gtf && mv tmp.gtf ${gtf}
fi

gtfToGenePred  -genePredExt ${gtf} genome_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile ${genome} genome_refGene.txt --out genome_refGeneMrna.fa

cp ref.len ${work_dir}/04.Annotation/
cp chrom.txt ${work_dir}/04.Annotation/

# annotation
cd ${work_dir}/04.Annotation
# 
Rscript select_variation.R -i ../03.Analysis/${filename}.filter.SNPs.table.gz -p ../03.Analysis/Parameter.csv
for i in `cut -f1 -d',' ${work_dir}/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }'`
do
	cd ${work_dir}/04.Annotation/${i}
	cut -f1,2 ${i}.Depth_information.txt | awk '{if(NR!=1) print }' > ${i}.SNPs.pos.txt
	bcftools view -R ${i}.SNPs.pos.txt -Oz -o ./${i}.SNPs.vcf.gz --threads ${thread} ${work_dir}/02.SNP_indel/${filename}.filter.SNPs.vcf.gz
	convert2annovar.pl --format vcf4old ./${i}.SNPs.vcf.gz --outfile ./${i}.SNPs.annovar.input
	annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${i}.SNPs.anno --exonsort ./${i}.SNPs.annovar.input ${work_dir}/refseq
	rm ${i}.Depth_information.csv ${i}.Depth_information.txt ${i}.SNPs.vcf.gz ${i}.parameter.txt
done
cd ${work_dir}/04.Annotation
Rscript select_variation.R -i ../03.Analysis/${filename}.filter.INDELs.table.gz -p ../03.Analysis/Parameter.csv
for i in `cut -f1 -d',' ${work_dir}/03.Analysis/Parameter.csv | awk '{if(NR!=1) print }'`
do
	cd ${work_dir}/04.Annotation/${i}
	cut -f1,2 ${i}.Depth_information.txt | awk '{if(NR!=1) print }' > ${i}.INDELs.pos.txt
	bcftools view -R ${i}.INDELs.pos.txt -Oz -o ./${i}.INDELs.vcf.gz --threads ${thread} ${work_dir}/02.SNP_indel/${filename}.filter.INDELs.vcf.gz
	convert2annovar.pl --format vcf4old ./${i}.INDELs.vcf.gz --outfile ./${i}.INDELs.annovar.input
	annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${i}.INDELs.anno --exonsort ./${i}.INDELs.annovar.input ${work_dir}/refseq
	rm ${i}.Depth_information.csv ${i}.Depth_information.txt ${i}.INDELs.vcf.gz ${i}.parameter.txt
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
