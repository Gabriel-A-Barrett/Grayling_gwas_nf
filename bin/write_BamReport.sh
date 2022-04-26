#!/bin/bash
# ISSUE writes one individual to tsv files and does not output float for mapping rate
#if [[ ! -e bamstats.tsv ]]; then # if file doesn't exist create it with header
#    echo -e "ind\traw_reads\taligned_reads\t0_mapping_quality\tmissmatch\terror_rate" > bamStats.tsv
#fi

STAT=$1

IND=$(basename $STAT .stats)
SEQ=$(cat ${STAT} | grep -P $'SN\tsequences:' | awk '{print $3}')
ALIG=$(cat ${STAT} | grep -P $'SN\treads mapped:' | awk '{print $4}')
MQ0=$(cat ${STAT} | grep -P $'SN\treads MQ0:' | awk '{print $4}')
MIS=$(cat ${STAT} | grep -P $'SN\tmismatches:' | awk '{print $3}')
ER=$(cat ${STAT} | grep -P $'SN\terror rate:' | awk '{print $4}')
MRATE=$(printf "%2.3f" $(($ALIG / $SEQ)))

echo -e "$IND\t$SEQ\t$ALIG\t$MQ0\t$MIS\t$ER" >> bamStats.tsv

echo -e "${IND}\t${MRATE}" >> mapRate.tsv

exit