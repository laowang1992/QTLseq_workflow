#!/bin/bash
# parse_config.sh
# 从 config.json 读取配置并导出为环境变量

config=${1:-config.json}

# 全局参数
export work_dir=$(jq -r '.global.work_dir' "$config")
export thread=$(jq -r '.global.thread' "$config")
export sampleInfo=$(jq -r '.global.sampleInfo' "$config")
export genome=$(jq -r '.global.genome' "$config")
export gff=$(jq -r '.global.gff' "$config")
export gtf=$(jq -r '.global.gtf' "$config")

# 软件
export picard=$(jq -r '.tools.picard' "$config")
export gatk=$(jq -r '.tools.gatk' "$config")
export DISCVRSeq=$(jq -r '.tools.DISCVRSeq' "$config")

# 变异检测
export index=$(jq -r '.call_vari.index' "$config")
export aligner=$(jq -r '.call_vari.aligner' "$config")            # 比对软件，bwotie2、mem2（bwa-mem2）、mem（bwa-mem）或star（STAR）
export sjdbOverhang=$(jq -r '.call_vari.sjdbOverhang' "$config")  # 该参数用于STAR建库，通常为read length-1
export sort=$(jq -r '.call_vari.sort' "$config")                  # 排序软件，sbb指sambamba，sts指samtools
export rmdup=$(jq -r '.call_vari.rmdup' "$config")                # 去重软件，sbb指sambamba，picard指picard
export filename=$(jq -r '.call_vari.filename' "$config")

# 引物设计
export primerConfig=$(jq -r '.primer_design.primerConfig' "$config")
export minProdLen=$(jq -r '.primer_design.minProdLen' "$config")
export maxProdLen=$(jq -r '.primer_design.maxProdLen' "$config")
export primerMinSize=$(jq -r '.primer_design.primerMinSize' "$config")
export primerOptSize=$(jq -r '.primer_design.primerOptSize' "$config")
export primerMaxSize=$(jq -r '.primer_design.primerMaxSize' "$config")
export primerMinGC=$(jq -r '.primer_design.primerMinGC' "$config")
export primerMaxGC=$(jq -r '.primer_design.primerMaxGC' "$config")
export primerMinTM=$(jq -r '.primer_design.primerMinTM' "$config")
export primerOptTM=$(jq -r '.primer_design.primerOptTM' "$config")
export primerMaxTM=$(jq -r '.primer_design.primerMaxTM' "$config")
export primerNumReturn=$(jq -r '.primer_design.primerNumReturn' "$config")

# 指控、统计
export DPwinSize=$(jq -r '.statistics.win_size' "$config")

# 有需要的话添加上RNA seq分析的部分
# geneEXPR 模块参数
#export GENEEXPR_ENABLED=$(jq -r '.geneEXPR.enabled' "$config")
#export GENEEXPR_METHOD=$(jq -r '.geneEXPR.method' "$config")
#export GENEEXPR_MIN=$(jq -r '.geneEXPR.min_expression' "$config")

# DEG 模块参数
#export DEG_ENABLED=$(jq -r '.DEG.enabled' "$config")
#export DEG_METHOD=$(jq -r '.DEG.method' "$config")
#export DEG_P=$(jq -r '.DEG.P_value' "$config")
#export DEG_LOG2FC=$(jq -r '.DEG.log2FC' "$config")

# 也可以检查参数是否为空
#if [[ -z "$THREADS" || -z "$INPUT_DIR" ]]; then
#    echo "[Error] Missing essential config in $config" >&2
#    exit 1
#fi

# 遍历所有导出的变量，自动展开 ${...}
# 将原来的.conf改成了json格式，但是jq只是返回字符串，不能解析嵌套的变量名
# 改成json完全就是看着标准，但是下面的`compgen -v`会返回所有环境变量感觉不太安全
# 以后不行再改回.conf文件吧
for var in $(compgen -v); do
    [[ "$var" == BASH_* ]] && continue
    eval "value=\${$var}"
    if [[ "$value" == *'${'* ]]; then
        eval "export $var=$(eval echo \"$value\")"
    fi
done
