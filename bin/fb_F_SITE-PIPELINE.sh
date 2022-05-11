
popmap="\$(awk -F',' '{print \$4"\t"\$7}' ${meta})"


# Subset Individuals
#########################
vcftools --gzvcf ${vcf} --remove ${probSamples} \
--recode --recode-INFO-all --stdout |

# Site-Based Filters
########################
vcftools --gzvcf - --mac 6 --minQ 20 \
--minDP 3 --min-meanDP 3 \
--maxDP 67 --max-meanDP 67 \
--recode --recode-INFO-all --stdout |

# PER-POP MISSINGNESS
########################## 
sh filter_miss_by_pop.sh - "\${popmap}" .65 1 |

# INFO
########
bcftools filter -i'AB>0.25 && AB<0.75 | AB < 0.01' | \
bcftools filter -i'SAF / SAR > 100 && SRF / SRR > 100 | SAR / SAF > 100 && SRR / SRF > 100' | \
bcftools filter -i'MQM / MQMR > 0.9 && MQM / MQMR < 1.05' | \
bcftools filter -i'PAIRED > 0.05 && PAIREDR > 0.05 && PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05' | \
bcftools filter -i'QUAL/ INFO/DP > .25' |

# HWE !Decompose into SNPs!
############################
vcfallelicprimitives --keep-info --keep-geno |
vcftools --vcf - --remove-indels --recode --recode-INFO-all --stdout |

# Having trouble... outputs *hwe.recode.vcf
# Does not find file even w/ full written
perl ${baseDir}/bin/filter_hwe_by_pop.pl -v - -p "\${popmap}" -h 0.001 -c .5 -o consensus
rm *Kup* *Sag* *Col*
#vcftools --vcf \$FILTERFILE5 --exclude-positions exclude.hwe --recode --recode-INFO-all --stdout | gzip -c > \$FILTERFILE6