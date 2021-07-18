#class of transcription boundaries sites(TSS OR TES)(TBS)
import normalization_functions
from cluster_functions import paraclu,distclu
import multiprocessing
import math
import itertools
import operator
from time import *
from functools import partial

class TBS_all_collection(object):
    '''
    class of transcription boundary sites collection.
    '''
    def __init__(self,TBS_all_list,TBS_strand_collection_list):
        self.TBS_strand_collection_list = TBS_strand_collection_list
        self.TBSlist = TBS_all_list
        self.reads_counts = [int(TBS.reads_count) for TBS in self.TBSlist]
        self.total_count = sum(self.reads_counts)
        # self.chrs = [TBS.chr for TBS in self.TBSlist]
        # self.pos = [TBS.pos for TBS in self.TBSlist]
        # self.strands = [TBS.strand for TBS in self.TBSlist]

    def TBS_filter(self,num):
        pass

    def normalization(self,method):
        if method == "simplecpm":
            self.normalized_counts = normalization_functions.simpleTpm(self.reads_counts, self.total_count)
        elif method == "powerlaw":
            slope,intercept = normalization_functions.fit_power_law_to_values(self.reads_counts)
            coefficients = [slope,intercept]
            self.normalized_counts = normalization_functions.fit_values_to_referent_power_law(self.reads_counts, coefficients)
        else:
            self.normalized_counts = self.reads_counts

        for i,TBS in enumerate(self.TBSlist):
            self.TBSlist[i].normalized_count = self.normalized_counts[i]
        
        g = 0
        for i,tbs_strand_collection in enumerate(self.TBS_strand_collection_list):
            for j,tbs_chr_collection in enumerate(tbs_strand_collection.tbs_chr_collection_list):
                for k,tbs in enumerate(tbs_chr_collection.tbs_chr_list):
                    self.TBS_strand_collection_list[i].tbs_chr_collection_list[j].tbs_chr_list[k].normalized_count = self.normalized_counts[g]
                    g += 1


    def clustering(self,cpu,params_dict,method = "paraclu"):
        all_tbs_chr_list = []
        for i,tbs_strand_collection in enumerate(self.TBS_strand_collection_list):
            strand = tbs_strand_collection.strand
            for j,tbs_chr_collection in enumerate(tbs_strand_collection.tbs_chr_collection_list):
                chro = tbs_chr_collection.chro
                cluster_list = []
                tss_pos_list = [int(tbs.pos) for tbs in tbs_chr_collection.tbs_chr_list if tbs.eventtype == "tss"]
                tes_pos_list = [int(tbs.pos) for tbs in tbs_chr_collection.tbs_chr_list if tbs.eventtype == "tes"]
                tss_cpm_list = [float(tbs.normalized_count) for tbs in tbs_chr_collection.tbs_chr_list if tbs.eventtype == "tss"]
                tes_cpm_list = [float(tbs.normalized_count) for tbs in tbs_chr_collection.tbs_chr_list if tbs.eventtype == "tes"]
                all_tbs_chr_list.append([strand,chro,tss_pos_list,tss_cpm_list,cluster_list,"tss"]) # three dimensional list
                all_tbs_chr_list.append([strand,chro,tes_pos_list,tes_cpm_list,cluster_list,"tes"])

        if method == "paraclu":
            p = multiprocessing.Pool(int(cpu))
            all_cluster_list = p.map(paraclu,all_tbs_chr_list)
            p.close()
            p.join()
        else:
            p = multiprocessing.Pool(int(cpu))
            maxdist = params_dict.get("maxdist")
            all_cluster_list = p.map(partial(distclu,max_dist = maxdist),all_tbs_chr_list)
            p.close()
            p.join()

        self.all_cluster_list = all_cluster_list
    
    @staticmethod
    def cluster_filter1(cluster_list_params):
        cluster_list = cluster_list_params[0]
        params = cluster_list_params[1]
        minStability,maxLength,removeSingletons,keepSingletonAbove = params
        cluster_filter1_list = []
        cluster_id = 1
        for j,cluster in enumerate(cluster_list):
            start = cluster.start
            end = cluster.end
            num_sites = cluster.num_sites
            dominant_site = cluster.dominant_site
            total_cpm = cluster.total_cpm
            dominant_site_cpm = cluster.dominant_site_cpm
            min_d = cluster.min_d
            max_d = cluster.max_d
            chro = cluster.chro
            strand = cluster.strand
            eventtype = cluster.eventtype
            if (max_d/min_d > minStability and end - start + 1 < maxLength):
                if (removeSingletons == True):
                    if (start == end and total_cpm > keepSingletonAbove) or (end != start):
                        cluster.cluster_id = cluster_id
                        cluster_filter1_list.append(cluster)
                        cluster_id += 1
                else:
                    cluster.cluster_id = cluster_id
                    cluster_filter1_list.append(cluster)
                    cluster_id += 1
        
        return cluster_filter1_list

    @staticmethod
    def cluster_filter2_step3(cluster_filter2_tmp_selected_list_cluster_filter1_list):
        cluster_filter2_tmp_selected_list,cluster_filter1_list = cluster_filter2_tmp_selected_list_cluster_filter1_list
        cluster_filter2_list = []
        cluster_filter2_tmp_selected_list_all_ids = [cluster[0] for cluster in cluster_filter2_tmp_selected_list]
        for cluster in cluster_filter1_list:
            if cluster.cluster_id in cluster_filter2_tmp_selected_list_all_ids:
                cluster_filter2_list.append(cluster)

        return cluster_filter2_list


    def cluster_filter(self,params_dict):
        cpu = params_dict.get("cpu")
        minStability = float(params_dict.get("minStability"))
        maxLength = int(params_dict.get("maxLength"))
        removeSingletons = bool(params_dict.get("removeSingletons"))
        keepSingletonAbove = float(params_dict.get("keepSingletonAbove"))
        reducetoNoneoverlap = bool(params_dict.get("reducetoNoneoverlap"))

        all_cluster_list = self.all_cluster_list
        params = [minStability,maxLength,removeSingletons,keepSingletonAbove]

        '''
        filter1
        '''

        all_cluster_list_params = []
        for i,cluster_list in enumerate(all_cluster_list):
            cluster_list_params = [cluster_list,params]
            all_cluster_list_params.append(cluster_list_params)

        p = multiprocessing.Pool(cpu)
        all_cluster_filter1_list = p.map(self.cluster_filter1,all_cluster_list_params)
        p.close()
        p.join()

        '''
        filter1 done
        '''
        if not reducetoNoneoverlap:
            self.all_cluster_list = all_cluster_filter1_list
        else:
            '''
            filter2
            step1 : Sort in ascending order of start, descending order of end
            '''

            all_cluster_filter2_tmp_list = []

            for cluster_filter1_list in all_cluster_filter1_list:
                cluster_filter2_tmp_list = [[cluster.cluster_id,cluster.start,cluster.end] for cluster in cluster_filter1_list]
                cluster_filter2_tmp_list = sorted(cluster_filter2_tmp_list, key = lambda x: (x[1], -x[2]))
                all_cluster_filter2_tmp_list.append(cluster_filter2_tmp_list)

            all_cluster_filter2_tmp_selected_list = []

            '''
            filter2
            step1 done
            '''

            '''
            fllter2
            step 2 :remove enveloped intervals.
            '''

            for cluster_filter2_tmp_list in all_cluster_filter2_tmp_list:
                cluster_filter2_tmp_selected_list = []
                prev_end = 0
                for cluster_info in cluster_filter2_tmp_list:
                    if cluster_info[2] > prev_end:
                        cluster_filter2_tmp_selected_list.append(cluster_info)
                        prev_end = cluster_info[2]
                all_cluster_filter2_tmp_selected_list.append(cluster_filter2_tmp_selected_list)
            
            '''
            filter 2
            step2 : done
            '''

            '''
            filter2
            step 3:time consuming,select clusters according to cluster_ids.
            '''

            p = multiprocessing.Pool(cpu)
            all_cluster_filter2_list = p.map(self.cluster_filter2_step3,list(zip(all_cluster_filter2_tmp_selected_list,all_cluster_filter1_list)))
            p.close()
            p.join()

            self.all_cluster_list = all_cluster_filter2_list

class TBS_strand_collection(object):
    '''
    class of transcription boundary sites for a specific strand collection.
    '''
    def __init__(self,strand,TBS_chr_collection_list):
        self.strand = strand
        self.tbs_chr_collection_list = TBS_chr_collection_list

class TBS_chr_collection(object):
    '''
    class of transcription boundary sites for a specific chromosome collection.
    '''
    def __init__(self,chro,TBS_chr_list):
        self.chro = chro
        self.tbs_chr_list = TBS_chr_list  #Elements of this list are TBS object.
        

class TBS(object):
    '''
    class of transcription boundary sites
    '''
    def __init__(self,chro,pos,strand,reads_count,boundary_type):
        self.chr = chro
        self.pos = pos      #should be integer.
        self.strand = strand #should be "+" or "-"
        self.reads_count = reads_count #should be integer.
        self.eventtype = boundary_type  #should be TSS or TES.
    

class RBS(object):
    '''
    class of reads boundary site.
    '''
    def __init__(self,chro,pos,strand,reads_count,boundary_type):
        self.chr = chro
        self.pos = pos      #should be integer.
        self.strand = strand #should be "+" or "-"
        self.reads_count = reads_count #should be integer.
        self.type = boundary_type  #should be start or end.

    def convert2TBS(self,):
        if (self.strand == "+"):
            if (self.type.lower() == "start"):
                return TBS(self.chr,self.pos,self.strand,self.reads_count,"TSS")
            elif (self.type.lower() == "end"):
                return TBS(self.chr,self.pos,self.strand,self.reads_count,"TES")
            else:
                raise Exception("sorry!Non standard input for the boundary column.")
        elif (self.strand == "-"):
            if (self.type.lower() == "start"):
                return TBS(self.chr,self.pos,self.strand,self.reads_count,"TES")
            elif (self.type.lower() == "end"):
                return TBS(self.chr,self.pos,self.strand,self.reads_count,"TSS")
            else:
                raise Exception("sorry!Non standard input for the boundary column.")
        else:
            raise Exception("sorry!Non standard input for the strand column.")

    def reads_filter(self,):
        pass

    
