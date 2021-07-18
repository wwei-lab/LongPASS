#######################################
#
#   script for combine replicates
#
#######################################

file="$1"

awk -F " "  '
    NF {a[$1"\t"$2"\t"$3]+=$4; b[$1"\t"$2"\t"$3]+=$5} 
    END {for(i in a)print i"\t"a[i]"\t"b[i]}
' ${file}|sort