import unittest
import LongPass.cluster_functions as cluster_functions
from LongPass.cluster import para_cluster


class TestClusterFunction(unittest.TestCase):
    '''
    test for cluster functions.
    '''
    def load_postestfile(self,filepath):
        inps = []
        outs = []
        max_dists = []
        with open(filepath,'r') as testfile:
            for line in testfile:
                if line.startswith("in:"):
                    inps.append([int(inp) for inp in line.split(':')[1].split(',')])
                elif line.startswith("max_dist:"):
                    max_dists.append(int(line.split(":")[1]))
                elif line.startswith("out:"):
                    outs.append([int(out) for out in line.split(':')[1].split(',')])
                else:
                    pass
        return [inps,max_dists,outs]


    
    def load_postpmtestfile(self,filepath):
        inps = []
        outs = []
        tpms = []
        with open(filepath,'r') as testfile:
            for line in testfile:
                if line.startswith("in:"):
                    inps.append([int(inp) for inp in line.split(':')[1].split(',')])
                elif line.startswith("tpm:"):
                    tpms.append([float(tpm) for tpm in line.split(':')[1].split(',')]) 
                elif line.startswith("outs:"):
                    outs.append([string2num(out) for out in line.split(':')[1].strip().split(',')])
                else:
                    pass
            
        return [inps,tpms,outs]

    def load_clusterlistfile(self,filepath):
        pos = []
        tpms = []
        starts = []
        ends = []
        num_sites = []
        dominant_sites = []
        total_tpms = []
        max_tpms = []
        minimal_densities = []
        maximal_densities = []
        with open(filepath ,'r') as testfile:
            for line in testfile:
                if line.startswith("pos:"):
                    pos.append([int(pos) for pos in line.split(":")[1].split(",")])
                elif line.startswith("tpms:"):
                    tpms.append([float(tpm) for tpm in line.split(':')[1].split(',')])
                elif line.startswith("starts:"):
                    starts.append([int(start) for start in line.split(':')[1].split(',')])
                elif line.startswith("ends:"):
                    ends.append([int(end) for end in line.split(':')[1].split(',')])
                elif line.startswith("num_sites:"):
                    num_sites.append([int(num_site) for num_site in line.split(':')[1].split(',')])
                elif line.startswith("dominant_sites:"):
                    dominant_sites.append([int(dominant_site) for dominant_site in line.split(":")[1].split(",")])
                elif line.startswith("total_tpms:"):
                    total_tpms.append([float(total_tpm) for total_tpm in line.split(":")[1].split(",")])
                elif line.startswith("max_tpms:"):
                    max_tpms.append([float(max_tpm) for max_tpm in line.split(":")[1].split(",")])
                elif line.startswith("minimal_densities:"):
                    minimal_densities.append([round(float(minimal_density),5) for minimal_density in line.split(":")[1].split(",")])
                elif line.startswith("maximal_densities:"):
                    maximal_densities.append([round(float(maximal_density),5) for maximal_density in line.split(":")[1].split(",")])
                else:
                    raise Exception("unexpected test content!")
            
        return [pos,tpms,starts,ends,num_sites,dominant_sites,total_tpms,max_tpms,minimal_densities,maximal_densities]



    def test_dist_cluster(self):
        postestlist = self.load_postestfile('/home/hongyanhong/TCfinder/test/test_pos_list.txt')
        inps = postestlist[0]
        max_dists = postestlist[1]
        outs = postestlist[2]
        testouts = []
        for inp,max_dist in zip(inps,max_dists):
            clusters = cluster_functions.distclu(inp,max_dist)
            starts = [cluster.start for cluster in clusters]
            ends = [cluster.end for cluster in clusters]
            testout = starts + ends
            testouts.append(testout)
        
        self.assertEqual(testouts,outs)

    def test_paraclu_findbreak(self):
        postpmtestlist = self.load_postpmtestfile('/home/hongyanhong/TCfinder/test/test_postpm_list.txt')
        inps = postpmtestlist[0]
        tpms = postpmtestlist[1]
        outs = postpmtestlist[2]
        testouts = []
        for inp,tpm in zip(inps,tpms):
            paramlist = cluster_functions.paraclu_findbreak(inp,tpm)
            testouts.append(paramlist)

        self.assertEqual(testouts,outs)

    def test_paracluster(self):
        testparaclusterlist = self.load_clusterlistfile("/home/hongyanhong/TCfinder/test/test_cluster_list.txt")
        pos = testparaclusterlist[0]
        tpms = testparaclusterlist[1]
        starts = testparaclusterlist[2]
        ends = testparaclusterlist[3]
        num_sites = testparaclusterlist[4]
        dominant_sites = testparaclusterlist[5]
        total_tpms = testparaclusterlist[6]
        max_tpms = testparaclusterlist[7]
        minimal_densities = testparaclusterlist[8]
        maximal_densities = testparaclusterlist[9]
        test_cluster_list = [] # one elemnet for one test case.
        test_outs_list = []
        for i in range(len(pos)):
            para_cluster_list = []
            
            for start,end,num_site,dominant_site,total_tpm,max_tpm,minimal_density, maximal_density in zip(starts[i],ends[i],num_sites[i],dominant_sites[i],total_tpms[i],max_tpms[i],minimal_densities[i],maximal_densities[i]):
                paracluster = para_cluster(start,end,num_site,dominant_site,total_tpm,max_tpm,minimal_density,maximal_density)
                para_cluster_list.append(paracluster)
            
            test_cluster_list.append(para_cluster_list)

        for p,tpm in zip(pos,tpms):
            cluster_functions.paraclu(p,tpm)
            test_outs_list.append(para_cluster_list)
            para_cluster_list = [] #清空全局变量

        self.assertEqual(test_cluster_list,test_outs_list)
        

def string2num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
