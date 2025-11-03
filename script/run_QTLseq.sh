#!/bin/bash
#  run_BSA.sh
#  
#  Copyright 2021 WangPF <wangpf0608@126.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#!/bin/bash
# Author:
#	WangPF
# Program:
# 	BSA 
# History:
# 	2021-02-20	First release  
# 	2023-03-31	Second release
set -euo pipefail

# 1. 加载软件环境
. ./env.sh

# 2. 解析配置文件（只解析一次）
. ./parse_config.sh config.json

cd ${work_dir}/script
# 这里为了省事，用两个脚本分开处理DNA测序和RNA测序，可以考虑后续整合成一个脚本
if [ "$aligner" = "star" ]; then
    sh call_vari_RNAseq.sh
elif [ "$aligner" = "mem" ] || [ "$aligner" = "mem2" ] || [ "$aligner" = "bowtie2" ]; then
    sh call_vari.sh
else
    echo "Error: unknown aligner '$aligner'"
    exit 1
fi
sh QTLseq.sh
sh annotation.sh
sh primer.sh
sh statistics.sh
