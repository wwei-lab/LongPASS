#! /bin/bash

#################################################################
#
#   script for split of the final result
#
#################################################################

#author : yanhong hong

if [ $# -eq 0 ]; then
    echo "bash result_sort_split.sh [original_resultfile] [outdir]";exit 0
fi

while [[ $1 == -* ]]; do
    case "$1" in
      -h|--help|-\?) echo "bash result_sort_split.sh [original_resultfile] [outdir]" ; exit 0;; 
    esac
done

resultfile=$1

outdir=$2

result_tss_sort=${outdir}"sort_tss_"`basename ${resultfile}`

result_tes_sort=${outdir}"sort_tes_"`basename ${resultfile}`

resultfile_tss_trim=${outdir}"trim_"`basename ${result_tss_sort}`

resultfile_tes_trim=${outdir}"trim_"`basename ${result_tes_sort}`

resultfile_tss_readsid=${outdir}"readsid_"`basename ${result_tss_sort}`

resultfile_tes_readsid=${outdir}"readsid_"`basename ${result_tes_sort}`

cat ${resultfile} | awk '$11 == "tss"' |sort -k1,1 -k3n,3 -k4n,4 -k2,2 > ${result_tss_sort}

cat ${resultfile} | awk '$11 == "tes"' |sort -k1,1 -k3n,3 -k4n,4 -k2,2 > ${result_tes_sort}

cat ${result_tss_sort} | rev | cut -f2- | rev | awk '{print NR"\t"$0}' > ${resultfile_tss_trim}

cat ${result_tes_sort} | rev | cut -f2- | rev | awk '{print NR"\t"$0}' > ${resultfile_tes_trim}

cat ${result_tss_sort} | rev | cut -f1 | rev |awk '{print NR"\t"$0}' > ${resultfile_tss_readsid}

cat ${result_tes_sort} | rev | cut -f1 | rev |awk '{print NR"\t"$0}' > ${resultfile_tes_readsid}