import pysam;
import sys;
import os;
import re;

def is_spliced(thisRead):
    if(re.search('(\d+)N',read.cigarstring)):
        return True
    else:
        return False

infile = sys.argv[1]
outfile = sys.argv[2]

bamFP = pysam.AlignmentFile(infile, "rb")
otfile = pysam.AlignmentFile(outfile, "wb", template=bamFP)

for read in bamFP:
    if( not( read.is_unmapped ) ):
        if( read.has_tag("ts") ):
            ts = read.get_tag("ts")
        else:
            ts = "*"
        read.qname = "%s:%s:%s:%s:%s" % (read.get_tag("RG"),read.get_tag("tp"),
                                ts, read.get_tag("NM"), read.qname)

        if( read.has_tag("SA") ):
            continue
        else:
            if( is_spliced(read)):
                # spliced reads
                otfile.write(read)
