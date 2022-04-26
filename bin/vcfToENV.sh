#!/bin/bash
VCF=$1
GROUP=$2 # BayPass_FullGroup_list.txt
LFMM=$3
ENV=$(echo $LFMM | sed -e 's/lfmm_//' -e 's/.env//')
PREFIX=$(echo $VCF | sed 's/.vcf.gz//')

awk '!/NA/ {print $1}' $LFMM > ${ENV}.indv
# retain samples with environmental data
bcftools view -I --force-samples -S ${ENV}.indv ${VCF} > ${ENV}_${VCF} # force samples not present in vcf header, -I = don't recaculate INFO fields
# Get final individual list 
bcftools query -l ${ENV}_${VCF} > ${ENV}_vcf.indv
plink --vcf ${ENV}_${VCF} --recode --allow-extra-chr --out "${ENV}_consensus" 

# Remove absent individuals in group_list
awk 'NR==FNR{a[$1];next} $1 in a' ${ENV}_vcf.indv $GROUP > joined${ENV}_Group_list.txt
# Remove absent individuals in LEA.env
awk -F'[ ]' 'NR==FNR{a[$1];next} $1 in a' ${ENV}_vcf.indv $LFMM > ${ENV}_m_lfmm.env

# convert vcf to baypass
cat ${ENV}_${VCF} | perl /home/gbarrett/nf_ArcticGrayling_LandscapeGenomics/bin/vcf2baypass.pl $GROUP ${ENV}.frq
