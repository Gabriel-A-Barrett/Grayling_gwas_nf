#!/bin/bash

VCF=$1
ENV=$2
GROUP=$3 # BayPass_FullGroup_list.txt
PREFIX=$(echo $VCF | sed 's/.vcf.gz//')

# Remove Na's and queue individuals from meta into env. 

bcftools query -l $VCF > consensus.indv

if [[ $ENV == 'temp' || $ENV == 'pit' ]]; then
    
    # remove samples with NA's and query list of ind. ID's
    awk '!/NA/ {print $1}' ../../meta/lfmm_${ENV}.env > ${ENV}.indv
    # retain samples with environmental data
    bcftools view -I --force-samples -S ${ENV}.indv ${VCF} > ${ENV}_${VCF} # force samples not present in vcf header, -I = don't recaculate INFO fields
    # Get individual list 
    bcftools query -l ${ENV}_${VCF} > ${ENV}_vcf.indv
    # Remove absent individuals in group_list
    awk 'NR==FNR{a[$1];next} !/NA/ $1 in a' ${ENV}_vcf.indv $GROUP > joined${ENV}_Group_list.txt
    # Remove absent individuals in LEA.env
    awk -F'[ ]' 'NR==FNR{a[$1];next} $1 in a' ${ENV}_vcf.indv lfmm_${ENV}.env > m_${ENV}_lfmm.env

    cat $VCF | perl vcf2baypass.pl $GROUP ${ENV}_${VCF}

elif [ $ENV == 'lat' ]; then

    awk 'NR==FNR{a[$1];next} $1 in a' consensus.indv $GROUP > joinedFull_Group_list.txt
    # Remove absent individuals in LEA.env
    awk -F'[ ]' 'NR==FNR{a[$1];next} $1 in a' consensus.indv lfmm_${ENV}.env > m_${ENV}_lfmm.env

    cat $VCF | perl vcftobaypass.pl ${GROUP} full_${VCF}

else # this applies to discrete and elevation gwas

    # Remove absent individuals in LEA.env
    awk -F'[ ]' 'NR==FNR{a[$1];next} $1 in a' consensus.indv lfmm_${ENV}.env > m_${ENV}_lfmm.env

fi

exit
