import sys
import pysam

bamfilepath = sys.argv[1]

tbsfilepath = sys.argv[2]

bamfile = pysam.AlignmentFile(bamfilepath,"rb")

# bam to tbs
tbs_dict = dict()
for read in bamfile:
    if read.flag >= 2048:
        continue
    readname = read.qname
    chrom = read.reference_name
    read.reference_start = read.reference_start + 1
    cigarlist = read.cigar
    if not read.has_tag("ts"):
        continue
    if read.get_tag("ts") == "+":
        if read.is_reverse:
            if cigarlist[0][0] == 4:
                if cigarlist[0][1] <= 10:
                    read.reference_start = int(read.reference_start) - cigarlist[0][1]
            if (chrom,read.reference_start,"-") not in tbs_dict:
                tbs_dict[(chrom,read.reference_start,"-")] = [0,1,[],[]]
            else:
                tbs_dict[(chrom,read.reference_start,"-")][1] += 1
            
            tbs_dict[(chrom,read.reference_start,"-")][3].append(readname)
                
            if (chrom,read.reference_end,"-") not in tbs_dict:
                tbs_dict[(chrom,read.reference_end,"-")] = [1,0,[],[]]
            else:
                tbs_dict[(chrom,read.reference_end,"-")][0] += 1
                
            tbs_dict[(chrom,read.reference_end,"-")][2].append(readname)
        else:
            reference_end = read.reference_end
            if cigarlist[-1][0] == 4:
                if cigarlist[-1][1] <= 10:
                    reference_end = int(read.reference_end) + cigarlist[-1][1]
            if (chrom,read.reference_start,"+") not in tbs_dict:
                tbs_dict[(chrom,read.reference_start,"+")] = [1,0,[],[]]
            else:
                tbs_dict[(chrom,read.reference_start,"+")][0] += 1
                
            tbs_dict[(chrom,read.reference_start,"+")][2].append(readname)
                
            if (chrom,reference_end,"+") not in tbs_dict:
                tbs_dict[(chrom,reference_end,"+")] = [0,1,[],[]]
            else:
                tbs_dict[(chrom,reference_end,"+")][1] += 1
                
            tbs_dict[(chrom,reference_end,"+")][3].append(readname)
    else:
        if read.is_reverse:
            reference_end = read.reference_end
            if cigarlist[-1][0] == 4:
                if cigarlist[-1][1] <= 10:
                    reference_end = int(read.reference_end) + cigarlist[-1][1]
            if (chrom,read.reference_start,"+") not in tbs_dict:
                tbs_dict[(chrom,read.reference_start,"+")] = [1,0,[],[]]
            else:
                tbs_dict[(chrom,read.reference_start,"+")][0] += 1
                
            tbs_dict[(chrom,read.reference_start,"+")][2].append(readname)
                
            if (chrom,reference_end,"+") not in tbs_dict:
                tbs_dict[(chrom,reference_end,"+")] = [0,1,[],[]]
            else:
                tbs_dict[(chrom,reference_end,"+")][1] += 1
                
            tbs_dict[(chrom,reference_end,"+")][3].append(readname)
        else:
            if cigarlist[0][0] == 4:
                if cigarlist[0][1] <= 10:
                    read.reference_start = int(read.reference_start) - cigarlist[0][1]
            if (chrom,read.reference_start,"-") not in tbs_dict:
                tbs_dict[(chrom,read.reference_start,"-")] = [0,1,[],[]]
            else:
                tbs_dict[(chrom,read.reference_start,"-")][1] += 1
                
            tbs_dict[(chrom,read.reference_start,"-")][3].append(readname)
                
            if (chrom,read.reference_end,"-") not in tbs_dict:
                tbs_dict[(chrom,read.reference_end,"-")] = [1,0,[],[]]
            else:
                tbs_dict[(chrom,read.reference_end,"-")][0] += 1
                
            tbs_dict[(chrom,read.reference_end,"-")][2].append(readname)

with open(tbsfilepath,"w") as f:
    for chrom_pos_strand,counts in tbs_dict.items():
        chrom = chrom_pos_strand[0]
        pos = chrom_pos_strand[1]
        strand = chrom_pos_strand[2]
        tsscount = counts[0]
        tescount = counts[1]
        tssreadsnames = counts[2]
        tesreadsnames = counts[3]
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,strand,tsscount,tescount,",".join(tssreadsnames) + ",",",".join(tesreadsnames) + ","))