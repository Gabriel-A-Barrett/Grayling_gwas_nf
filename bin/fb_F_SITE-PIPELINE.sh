#!/bin/bash
#SBATCH --job-name=fb_dedup_F_filter
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=6G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
module load vcftools
module load bcftools
module load vcflib
#cd /labs/Urban/Urban_RAD_ArcticGrayling/LG_PIPELINE
FILE=$1
PREFIX=$(basename $1 .vcf.gz)
POPMAP=$2
BAD_SAMPLES=$3
SUBSET=${PREFIX}.F.vcf.gz
FILTERFILE=${PREFIX}.F.mac.6.minQ20.minDP3.maxDP67.vcf.gz
FILTERFILE2=$(echo $FILTERFILE | sed 's/vcf.gz/miss65.vcf/')
FILTERFILE3=$(echo $FILTERFILE2 | sed 's/vcf/ab.bs.mq.ps.pp.vcf/')
FILTERFILE4=$(echo $FILTERFILE3 | sed 's/vcf/snps.vcf/')
FILTERFILE5=$(echo $FILTERFILE4 | sed 's/vcf/ri.vcf/')
FILTERFILE6=$(echo $FILTERFILE6 | sed 's/vcf/hwe.vcf')
# Subset Individuals
#########################
vcftools --gzvcf $FILE --remove $BAD_SAMPLES \
--recode --recode-INFO-all --stdout | gzip -c > $SUBSET
# Filter Low Quality Sites
###############################
vcftools --gzvcf $SUBSET --mac 6 --minQ 20 \
--minDP 3 --min-meanDP 3 \
--maxDP 67 --max-meanDP 67 \
--recode --recode-INFO-all --stdout | gzip -c > $FILTERFILE
# PER-POP MISSINGNESS
########################## 
POPMAP=REGMAP.txt
sh filter_miss_by_pop.sh $FILTERFILE $POPMAP .65 1 $FILTERFILE2
rm *Kup* *Sag* *Col*
# INFO
########
bcftools filter -i'AB>0.25 && AB<0.75 | AB < 0.01' $FILTERFILE2 | \
bcftools filter -i'SAF / SAR > 100 && SRF / SRR > 100 | SAR / SAF > 100 && SRR / SRF > 100' | \
bcftools filter -i'MQM / MQMR > 0.9 && MQM / MQMR < 1.05' | \
bcftools filter -i'PAIRED > 0.05 && PAIREDR > 0.05 && PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05' | \
bcftools filter -i'QUAL/ INFO/DP > .25' > $FILTERFILE3
# HWE !Decompose into SNPs!
############################
vcfallelicprimitives $FILTERFILE3 --keep-info --keep-geno > $FILTERFILE4
vcftools --vcf $FILTERFILE4 --remove-indels --recode --recode-INFO-all --stdout > $FILTERFILE5
# Having trouble... outputs *hwe.recode.vcf
perl ${baseDir}/bin/filter_hwe_by_pop.pl -v $FILTERFILE5 -p $POPMAP -h 0.001 -c .5 -o consensus.vcf
rm *Kup* *Sag* *Col*
#vcftools --vcf $FILTERFILE5 --exclude-positions exclude.hwe --recode --recode-INFO-all --stdout | gzip -c > $FILTERFILE6