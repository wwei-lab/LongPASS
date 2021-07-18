import subprocess
prefix = 'mixA_outcluster5'
file1 = open(prefix+'_tss.bed','w')
file2 = open(prefix+'_tes.bed','w')

try:
    command1 = "cat " + prefix + ".txt|grep tss > " + prefix + "_tss.txt"
    command2 = "cat " + prefix + ".txt|grep tes > " + prefix + "_tes.txt"
    subprocess.check_call(command1,shell = True)
    subprocess.check_call(command2,shell = True)
except:
    pass

with open(prefix+"_tss.txt",'r') as f:
    for line in f:
        strand = line.strip('\n').split('\t')[1]
        chro = line.strip('\n').split('\t')[0]
        start = line.strip('\n').split('\t')[2]
        end = line.strip('\n').split('\t')[3]
        name = line.strip('\n').split('\t')[5]
        score = line.strip('\n').split('\t')[7]

        writeline = "%s\t%s\t%s\tcluster_dominant:%s\t%s\t%s\n" %(chro,int(start)-1,end,name,score,strand)
        file1.write(writeline)

f.close()

with open(prefix+"_tes.txt",'r') as f:
    for line in f:
        strand = line.strip('\n').split('\t')[1]
        chro = line.strip('\n').split('\t')[0]
        start = line.strip('\n').split('\t')[2]
        end = line.strip('\n').split('\t')[3]
        name = line.strip('\n').split('\t')[5]
        score = line.strip('\n').split('\t')[7]

        writeline = "%s\t%s\t%s\tcluster_dominant:%s\t%s\t%s\n" %(chro,int(start)-1,end,name,score,strand)
        file2.write(writeline)

f.close()   
