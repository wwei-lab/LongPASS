#######################################
#
#   script for combine replicates
#
#######################################

#This script was used for filtering, but is now deprecated.

file="$1"
threshold="$2"

awk -v var="$threshold" -F " "  '
    NF {a[$1"\t"$2"\t"$3]+=$4; b[$1"\t"$2"\t"$3]+=$5} 
    END {for(i in a)
        {
            if (a[i] < var) a[i]=0
            if (b[i] < var) b[i]=0 
        };
        for(i in a)
        {
            if(a[i] != 0 || b[i] != 0)   print i"\t"a[i]"\t"b[i]
        };
    }
' ${file}|sort