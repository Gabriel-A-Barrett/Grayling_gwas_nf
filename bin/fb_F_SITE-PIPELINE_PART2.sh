#!/bin/bash
cd /labs/Urban/Urban_RAD_ArcticGrayling/freebayes/test
module load vcftools

FILTERFILE5="../results/fb.F.mac3.minQ20.minDP3.maxDP127.miss65.ab.bs.mq.ps.pp.snps.ri.vcf"
FILTERFILE6=$(echo $FILTERFILE5 | sed 's/vcf/hwe.vcf/')
FILTERFILE7=$(echo $FILTERFILE6 | sed 's/mac3/maf.01/')

awk '{print $1"\t"$2}' filtered.hwe > badloci.hwe
vcftools --vcf $FILTERFILE5 --exclude-positions badloci.hwe --recode --recode-INFO-all --stdout > $FILTERFILE6
vcftools --vcf $FILTERFILE6 --maf 0.01 --recode --recode-INFO-all --stdout > $FILTERFILE7

