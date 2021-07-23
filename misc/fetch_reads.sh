#! /bin/bash

#######################################################################
#
#   This script is used for extracting reads using cluster index
#
#######################################################################

if [ $# -eq 0 ]; then
    echo "bash fetch_reads.sh [reads_resultfile] [cluster_index]";exit 0
fi

while [[ $1 == -* ]]; do
    case "$1" in
      -h|--help|-\?) echo "bash fetch_reads.sh [reads_resultfile] [cluster_index]" ; exit 0;; 
    esac
done

readsfile=$1
clusterindex=$2


cat ${readsfile} | awk -v var="${clusterindex}" 'NR == var'