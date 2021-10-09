# QTLseq_workflow
from NGS data to QTL-seq result
# 一 QTL-seq原理
QTL-seq<sup>[[1](#ref)]</sup>是一种将Bulked‐segregant analysis (BSA)<sup>[[2](#ref),[3](#ref)]</sup>和高通量测序相结合，快速定位QTL的方法。目标性状有差异的双亲构建的分离群体中分别选取极端表型个体进行等量混合构建两个极端表型bulked DNA pool并进行测序。随后进行变异分析筛选出双亲间SNP位点并分别计算两个bulked DNA pool中每个SNP位点上某一亲本基因型read覆盖深度占该位点总read深度的比值，即SNP index，通过两个bulked DNA pool的SNP index相减即得到ΔSNP index。在基因组所有区域中，目标基因及其连锁的区域由于根据表型受到相反的选择在两个bulked DNA pool中表现出不同的趋势，因此ΔSNP index会显著偏离0附近；另一方面，于目标形状无关的区域则两个bulked DNA pool则表现为相似的变化趋势，因此ΔSNP index会在0附近波动（图1）。
<br/>
![图1 QTL-seq原理示意图](https://user-images.githubusercontent.com/35584208/121320629-e9595300-c93f-11eb-8745-e1f153dd59c6.png)
# 二 QTL-seq流程图
![image](https://user-images.githubusercontent.com/35584208/121321053-5967d900-c940-11eb-9d22-be6614ce71c0.png)
# 三 QTL-seq分析流程
## 1	取样、提取DNA、建库和测序
选取分离群体中极端表型个体各30-50株以及双亲取幼嫩组织，分别提取DNA并将分离群体极端表型个体按照等量原则混合构建bulked DNA pools，然后进行建库和高通量测序（混池DNA深度宜与单株个数相同）。（具体流程应参考实验设计以及公司测序报告）
## 2	数据过滤
使用fastp<sup>[[4](#ref)]</sup>（version: 0.20.0）对raw data进行过滤得到clean data。统计过滤前后total bases、total reads、Q30、Q20、GC content以及有效数据比率（data_stat.csv/txt），同时使用FastQC（version: 0.11.9）对过滤前后的数据进行质量评估（QC/sample_fastqc.html）。
## 3	比对到参考基因组
使用Bowtie2<sup>[[5](#ref)]</sup>（version: 2.4.1）软件将clean read比对到甘蓝型油菜ZS11参考基因组<sup>[[6](#ref)]</sup>，得到SAM（Sequence Alignment/Map）格式文件并对比对结果进行统计（01.Mapping/align_stat.csv），随后使用SAMtools（version: 1.9）对比对结果（SAM文件）按照染色体和位置进行排序并转换为BAM（Binary Alignment/Map）格式文件<sup>[[7](#ref)]</sup>，使用Picard tools<sup>[[8](#ref)]</sup>（version: 2.23.2）去除建库过程中产生的PCR重复并统计样本基因组覆盖率（01.Mapping/cov_stat.txt）。
## 4	变异分析
变异分析使用Genome Analysis Toolkit，GATK<sup>[[9](#ref)]</sup>（version: 3.8-0-ge9d806836）完成，首先使用GATK的HaplotypeCaller功能对样本单独分析再使用CombineGVCFs功能合并，随后使用GenotypeGVCFs功能得到SNP和INDEL信息，最后使用VariantFiltration功能过滤原始的变异位点的到可靠的变异信息。
## 5	QTL-seq分析
在变异分析得到的可靠SNP位点中，筛选亲本内纯和并且亲本间不同的位点，并且进一步过滤掉双亲和两各混池中低覆盖深度的位点，此时统计剩余SNP位点在基因组上的分布情况并绘图。利用最终筛选出的SNP进行下一步的分析，首先计算每个混池的SNP index值，随后计算ΔSNP index值并绘图，其中点图按照500kb窗口大小、250kb步长进行滑窗统计，折线图是由R软件（version: 4.0.2）扩展包QTLseqr<sup>[[10](#ref)]</sup>（version: 0.7.5.2）按照2Mb窗口大小统计得到，置信区间按照Takagi 等（2013）<sup>[[1](#ref)]</sup>描述方法计算得到，超出95或99%置信区间的区域视为目标性状候选QTL。
## 6	引物设计（若需要请单独备注）
为进行QTL验证和进一步精细定位，根据变异分析中得到的可靠INDEL位点进行全基因组范围的引物设计。

我们首先需要确定引物设计的位点。在vcf文件的第10及以后的列为各样本信息，从中可以得到各样本基因型，对应第9列GT位置，两个数字中间用“ / ”分开，这两个数字表示二倍体sample的基因型。0表示样本中有ref的allele；1表示样本中有variant的allele；2表示有第二个variant的allele。

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	D	P	R	R3	W	Y	qW	qY
scaffoldA01	193227	.	GAATGAATGA	G	361.22	PASS	AC=3;AF=0.167;AN=18;BaseQRankSum=0.493;ClippingRankSum=0.00;DP=67;ExcessHet=0.4012;FS=0.000;MLEAC=3;MLEAF=0.167;MQ=41.18;MQRankSum=0.253;QD=20.07;ReadPosRankSum=-7.280e-01;SOR=1.022	GT:AD:DP:GQ:PL	0/0:7,0:7:21:0,21,232	0/0:8,0:8:9:0,9,243	0/0:5,0:5:15:0,15,170	0/0:5,0:5:15:0,15,175	1/1:0,8:8:24:360,24,0	0/0:8,0:8:24:0,24,308	0/0:4,0:4:9:0,9,135	0/1:8,2:10:60:60,0,847	0/0:12,0:12:33:0,33,495
scaffoldA01	276029	.	CA	CAA,C	1809.04	PASS	AC=10,4;AF=0.556,0.222;AN=18;BaseQRankSum=-3.850e-01;ClippingRankSum=0.00;DP=107;ExcessHet=4.8226;FS=15.247;MLEAC=10,4;MLEAF=0.556,0.222;MQ=41.00;MQRankSum=0.00;QD=19.04;ReadPosRankSum=0.00;SOR=1.026	GT:AD:DP:GQ:PL	1/2:0,5,2:7:55:167,55,55,130,0,162	1/2:1,3,4:10:74:157,85,124,74,0,121	0/1:1,2,0:3:17:44,0,17,47,23,70	1/1:0,12,0:12:36:317,36,0,317,36,317	2/2:0,0,10:10:30:278,278,278,30,30,0	1/1:1,12,0:13:12:291,12,0,294,36,317	0/1:3,8,0:11:30:179,0,30,187,54,242	0/1:3,6,0:9:37:116,0,37,125,54,179	0/1:7,15,0:22:95:316,0,95,336,140,476

```

在鉴定的INDELs中筛选两个亲本内纯和且在亲本间不同的位点，即两亲本基因型分别为0/0、1/1或者1/1、0/0的位点，再根据INDEL位置提取上下游各250 bp基因组序列，随后使用primer3<sup>[[11](#ref)]</sup>（version: 2.5.0）进行引物设计（引物退火温度、长度、GC含量和扩增片段长短等偏好需提前说明，同时需根据实验室电泳仪器能区分的INDEL大小挑选合适的位点，引物根据参考基因组设计，其真实位置和特异性需根据序列比对和实验结果进行验证）。
 <br>
 <br>

---
 <div id="ref"></div>

 # 参考文献
 - [1]	Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi. QTL ‐seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations[J]. The Plant Journal, 2013, 74(1).
- [2]	Giovannoni James J., Wing Rod A., Ganal Martin W., Tanksley Steven D..  Isolation of molecular markers from specific chromosomal intervals using DNA pools from existing mapping populations[J]. Narnia, 1991, 19(23).
- [3]	R. W. Michelmore, I. Paran,R. V. Kesseli. Identification of Markers Linked to Disease-Resistance Genes by Bulked Segregant Analysis: A Rapid Method to Detect Markers in Specific Genomic Regions by Using Segregating Populations[J]. Proceedings of the National Academy of Sciences of the United States of America, 1991, 88(21).
- [4]	Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ preprocessor[J]. Bioinformatics, 2018, 34(17).
- [5]	Ben Langmead, Steven L Salzberg. Fast gapped-read alignment with Bowtie 2[J]. Nature Methods: Techniques for life scientists and chemists, 2012, 9(4).
- [6]	Jia-Ming Song, Zhilin Guan, Jianlin Hu, Chaocheng Guo, Zhiquan Yang, Shuo Wang, Dongxu Liu, Bo Wang, Shaoping Lu, Run Zhou, Wen-Zhao Xie, Yuanfang Cheng, Yuting Zhang, Kede Liu, Qing-Yong Yang, Ling-Ling Chen, Liang Guo. Eight high-quality genomes reveal pan-genome architecture and ecotype differentiation of Brassica napus[J]. Nature Plants, 2020, 6(1).
- [7]	Li Heng, Handsaker Bob, Wysoker Alec, Fennell Tim, Ruan Jue, Homer Nils, Marth Gabor, Abecasis Goncalo, Durbin Richard. The Sequence Alignment/Map format and SAMtools.[J]. Bioinformatics (Oxford, England), 2009, 25(16).
- [8]	“Picard Toolkit.” 2019. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/; Broad Institute
- [9]	Aaron McKenna, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian Cibulskis, Andrew Kernytsky, Kiran Garimella, David Altshuler, Stacey Gabriel, Mark Daly, Mark A. DePristo. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data[J]. Cold Spring Harbor Laboratory Press, 2010, 20(9).
- [10]	Ben N. Mansfeld, Rebecca Grumet. QTLseqr: An R Package for Bulk Segregant Analysis with Next‐Generation Sequencing[J]. The Plant Genome, 2018, 11(2).
- [11]	Untergasser Andreas, Cutcutache Ioana, Koressaar Triinu, Ye Jian, Faircloth Brant C., Remm Maido, Rozen Steven G.. Primer3—new capabilities and interfaces[J]. Narnia, 2012, 40(15).
