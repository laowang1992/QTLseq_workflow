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
# 	2021/02/20	First release  

## 加载配置文件
. ./.conf

## call variation
sh call_vari.sh

##
cd ${work_dir}/03.Analysis
Rscript bsa.R \
	--input ${filename}.filter.SNPs.txt --out ${out} \
	--highP ${highP} --lowP ${lowP} --highB ${highB} --lowB ${lowB} \
	--minQ ${minQ} \
	--minHPdp ${minHPdp} --maxHPdp ${maxHPdp} \
	--minLPdp ${minHPdp} --maxLPdp ${maxLPdp} \
	--minHBdp ${minHBdp} --maxHBdp ${maxHBdp} \
	--minLBdp ${minHBdp} --maxLBdp ${maxLBdp} \
	--winSize ${winSize1} --winStep ${winStep1} \
	--minN ${minN} \
	--width ${width} --height ${height}

Rscript QTLseqr.R \
	--input ${filename}.filter.SNPs.table --out ${out} \
	--highP ${highP} --lowP ${lowP} --highB ${highB} --lowB ${lowB} \
	--bulkSuieH ${bulkSuieH} --bulkSuieL ${bulkSuieL} \
	--minQ ${minQ} \
	--minHPdp ${minHPdp} --maxHPdp ${maxHPdp} \
	--minLPdp ${minHPdp} --maxLPdp ${maxLPdp} \
	--minHBdp ${minHBdp} --maxHBdp ${maxHBdp} \
	--minLBdp ${minHBdp} --maxLBdp ${maxLBdp} \
	--popType ${popType} \
	--winSize ${winSize2} \
	--minN ${minN} \
	--width ${width} --height ${height}

## primer
cd ${work_dir}/04.Primer
#sh primer.sh
