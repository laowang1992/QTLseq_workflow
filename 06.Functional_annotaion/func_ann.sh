


cd ${work_dir}/06.Functional_annotation
gff=../refseq/zs11.v0.gff3
bed=../refseq/zs11.v0.bed
ann=
chr=
start=
end=

python -m jcvi.formats.gff bed --type=gene --key=ID ${gff} -o ${bed}

Rscript extr_func_ann.R $bed $ann $chr $start $end
