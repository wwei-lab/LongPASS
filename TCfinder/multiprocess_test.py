# import numpy as np
# import multiprocessing as mp
# import os



# class Tester:

#     num = 0.0
#     name = 'none'
#     def __init__(self,tnum=num, tname=name):
#         self.num  = tnum
#         self.name = tname

#     def __str__(self):
#         return '%f %s' % (self.num, self.name)

# def mod(test, nn, out_queue):
#     print (test.num)
#     test.num = np.random.randn()
#     print (test.num)
#     test.name = nn
#     out_queue.put(test)

# def test2(l):
#     print(l,os.getpid())


# if __name__ == '__main__':       
#     # num = 1
#     # out_queue = mp.Queue()
#     # tests = np.empty(num, dtype=object)
#     # for it in range(num):
#     #     tests[it] = Tester(tnum=it*1.0)


#     # print ('\n')
#     # workers = [ mp.Process(target=mod, args=(test, 'some', out_queue) ) for test in tests ]

#     # for work in workers: work.start()

    

#     # for work in workers: work.join()

#     # res_lst = []
#     # for j in range(len(workers)):
#     #     res_lst.append(out_queue.get())

#     # for test in res_lst: print (test)
#     list1 = [[1,2,3],[3,4,5]]

#     p = mp.Pool(5)

#     p.map(test2,list1)

#     p.close()
#     p.join()

f = open('/home/hongyanhong/TSS_TES_project/start_end.ctss5','w')

with open('/home/hongyanhong/TSS_TES_project/start_end.ctss','r') as file:
    for line in file:
        chro = line.split('\t')[0]
        pos = line.split('\t')[1]
        strand = line.split('\t')[2]
        start = line.split('\t')[3]
        end = line.split('\t')[4].strip()
        if strand == "-":
            f.write("%s\t%s\t%s\t%s\t%s\n"%(chro,pos,strand,end,start))
        else:
            f.write(line)

f.close()
        
    