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
      --thread=${thread} --length_required 50

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
