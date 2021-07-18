file1 = open('out_cluster2.bed','w')

with open("out_cluster2.txt",'r') as file2:
    for line in file2:
        strand = line.strip('\n').split('\t')[1]
        chro = line.strip('\n').split('\t')[0]
        start = line.strip('\n').split('\t')[2]
        end = line.strip('\n').split('\t')[3]
        name = line.strip('\n').split('\t')[5]

        writeline = "chr%s\t%s\t%s\t%s\t%s\t%s\n" %(chro,start,end,name,"0",strand)

        file1.write(writeline)

file1.close()
