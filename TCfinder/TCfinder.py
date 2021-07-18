"""
Copyright 2021 Yanhong Hong (a524664992@gmail.com)

This module is the main entry of TCfinder.It is executed when user runs `run_TCfinder.py` 
(directly from the source directory).
"""

from misc import load_lrtsp_objects
import logging
import misc
import sys
import os
from multiprocessing import Pool,Queue,Process,Manager
import multiprocessing as mp
import TBS
from time import *
from functools import partial
import tempfile
import subprocess

# def combine_replicate_file(all_combine_file_path,boundary_type):
#     prev = []
#     prev_start = 0
#     prev_end = 0
#     prev_site = 0
#     unique_all_combine_file_path = "unique_" + all_combine_file_path
#     f = open(unique_all_combine_file_path,'w')
#     with open(all_combine_file_path,'r') as all_combine_file:
#         for line in all_combine_file:
#             try:
#                 chro_pos_strand = line.split('\t')[:3]
#             except Exception as e:
#                 print("check your input!")
#             if boundary_type == "tbs":
#                 try:
#                     start = line.split('\t')[3]
#                     end = line.split('\t')[4].strip()
#                 except Exception as e:
#                     print("check your input!")
#             else:
#                 try:
#                     site = line.split('\t')[3].strip()
#                 except Exception as e:
#                     print("check your input!")
#             if chro_pos_strand == prev:
#                 if boundary_type == "tbs":
#                     prev_start = prev_start + start
#                     prev_end = prev_end + end
#                 else:
#                     prev_site = prev_site + site
#             else:
#                 if prev == []: #initial type
#                     prev = chro_pos_strand
#                     if boundary_type == "tbs":
#                         prev_start = start
#                         prev_end = end
#                     else:
#                         prev_site = site
#                     continue
#                 if boundary_type == "tbs":
#                     newline = "%s\t%s\t%s\t%s\t%s\n" % (prev[0],prev[1],prev[2],prev_start,prev_end)
#                     f.write(newline)
#                     prev = chro_pos_strand
#                     if boundary_type == "tbs":
#                         prev_start = start
#                         prev_end = end
#                     else:
#                         prev_site = site
        
#         newline = "%s\t%s\t%s\t%s\t%s\n" % (prev[0],prev[1],prev[2],prev_start,prev_end)
#         f.write(newline)





def load_tbs_objects(logger,args):    
    logger.info("loading files...\n")


    tssfiles = args.tssfile
    tesfiles = args.tesfile
    tbsfiles = args.tbsfile

    tss_collections_dict = {}
    tes_collections_dict = {}
    tbs_collections_dict = {}

    if len(tssfiles) > 0:
        if not args.replicate:
            for tssfilepath in tssfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tssfilepath))
                tss_collection = misc.load_lrtsp_objects(tssfilepath, "tss")
                tss_collections_dict[os.path.basename(tssfilepath)] = tss_collection
        else:
            combine_command = "cat "
            for tssfilepath in tssfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tssfilepath))
                combine_command = combine_command + tssfilepath + " "
            combine_command += "|sort -k1,1 -k2n,2 > all_combine_tssfiles.lrtsp"
            try:
                combine_status = subprocess.check_call(combine_command,shell=True)
            except:
                pass
            script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/misc/combine_tbsfile.sh"
            combine_command2 = "bash " + script_path + " all_combine_tssfiles.lrtsp |awk  '{print $1,$2,$3,$4}'|sed 's/ /\t/g'|sort -k1,1 -k2n,2 > unique_all_combine_tssfiles.lrtsp"
            try:
                combine2_status = subprocess.check_call(combine_command2,shell=True)
            except:
                pass
            tss_collection = misc.load_lrtsp_objects("unique_all_combine_tssfiles.lrtsp", "tss")
            tss_collections_dict["unique_all_combine_tssfiles.lrtsp"] = tss_collection
            try:
                subprocess.check_call("rm unique_all_combine_tssfiles.lrtsp",shell=True)
                subprocess.check_call("rm all_combine_tssfiles.lrtsp",shell=True)
            except:
                pass

    if len(tesfiles) > 0:
        if not args.replicate:
            for tesfilepath in tesfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tesfilepath))
                tes_collection = misc.load_lrtsp_objects(tesfilepath, "tes")
                tes_collections_dict[os.path.basename(tesfilepath)] = tes_collection
        else:
            combine_command = "cat "
            for tesfilepath in tesfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tesfilepath))
                combine_command = combine_command + tesfilepath + " "
            combine_command += "|sort -k1,1 -k2n,2 > all_combine_tesfiles.lrtsp"
            try:
                combine_status = subprocess.check_call(combine_command,shell=True)
            except:
                pass
            script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/misc/combine_tbsfile.sh"
            combine_command2 = "bash " + script_path + " all_combine_tesfiles.lrtsp |awk  '{print $1,$2,$3,$4}'|sed 's/ /\t/g'|sort -k1,1 -k2n,2 > unique_all_combine_tesfiles.lrtsp"
            try:
                combine2_status = subprocess.check_call(combine_command2,shell=True)
            except:
                pass
            tes_collection = misc.load_lrtsp_objects("unique_all_combine_tesfiles.lrtsp", "tes")
            tes_collections_dict["unique_all_combine_tesfiles.lrtsp"] = tes_collection
            try:
                subprocess.check_call("rm unique_all_combine_tesfiles.lrtsp",shell=True)
                subprocess.check_call("rm all_combine_tesfiles.lrtsp",shell=True)
            except:
                pass

    if len(tbsfiles) > 0:
        if not args.replicate:
            '''
            TBS file should look like this:
            chr pos strand startcount endcount
            1   14000   +   2   0
            1   15361   -   0   1
            '''
            for tbsfilepath in tbsfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tbsfilepath))
                tbs_collection = misc.load_lrtsp_objects(tbsfilepath, "tbs")
                tbs_collections_dict[os.path.basename(tbsfilepath)] = tbs_collection

        else:
            combine_command = "cat "
            for tbsfilepath in tbsfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tbsfilepath))
                combine_command = combine_command + tbsfilepath + " "
            combine_command += "|sort -k1,1 -k2n,2 > all_combine_tbsfiles.lrtsp > all_combine_tbsfiles.lrtsp"
            try:
                combine_status = subprocess.check_call(combine_command,shell=True)
            except:
                pass
            script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/misc/combine_tbsfile.sh"
            combine_command2 = "bash " + script_path + " all_combine_tbsfiles.lrtsp |sort -k1,1 -k2n,2 > unique_all_combine_tbsfiles.lrtsp"
            print (combine_command2)
            try:
                combine2_status = subprocess.check_call(combine_command2,shell=True)
            except:
                pass
            tbs_collection = misc.load_lrtsp_objects("unique_all_combine_tbsfiles.lrtsp", "tbs")
            tbs_collections_dict["unique_all_combine_tbsfiles.lrtsp"] = tbs_collection
            try:
                pass
                # subprocess.check_call("rm all_combine_tbsfiles.lrtsp",shell=True)
                subprocess.check_call("rm unique_all_combine_tbsfiles.lrtsp",shell=True)
            except:
                pass

    return [tss_collections_dict,tes_collections_dict,tbs_collections_dict]

def normalization(ts_collection,method):
    ts_collection.normalization(method)
    return ts_collection

def read_params_dict(args):
    if args.params == "":
        params_dict = {
            "value_range" : [10,1000],
            "maxdist" : 20,
            "minStability" : 1,
            "maxLength" : 500,
            "removeSingletons" : True,
            "keepSingletonAbove" : 1.1,
            "reducetoNoneoverlap" : True,
        }
    else:
        params_dict = {}
        with open(args.params,'r') as paramsfile:
            for line in paramsfile:
                if ":" in line:
                    pass
                else:
                    key = line.split("=")[0].strip()
                    value = line.split("=")[1].strip()
                    params_dict[key] = value

    return params_dict

def writing(args,ts_tmp_collections_list):
    if args.clustering == "paraclu":
        with open(args.o,'w') as file:
            ts_collection_all_clusters_list = ts_tmp_collections_list[0].all_cluster_list
            for ts_collection_clusters in ts_collection_all_clusters_list:
                for cluster in ts_collection_clusters:
                    start = cluster.start
                    end = cluster.end
                    strand = cluster.strand
                    chro = cluster.chro
                    num_sites = cluster.num_sites
                    dominant_site = cluster.dominant_site
                    total_cpm = cluster.total_cpm
                    dominant_site_cpm = cluster.dominant_site_cpm
                    max_d = cluster.max_d
                    min_d = cluster.min_d
                    event_type = cluster.eventtype
                    line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chro,strand,start,end,num_sites,dominant_site,dominant_site_cpm,round(total_cpm,5),round(max_d,5),round(min_d,5),event_type)
                    file.write(line)
    else:
        with open(args.o,'w') as file:
            ts_collection_all_clusters_list = ts_tmp_collections_list[0].all_cluster_list
            for ts_collection_clusters in ts_collection_all_clusters_list:
                for cluster in ts_collection_clusters:
                    start = cluster.start
                    end = cluster.end
                    strand = cluster.strand
                    chro = cluster.chro
                    event_type = cluster.eventtype
                    line = "%s\t%s\t%s\t%s\t%s\n" % (chro,strand,start,end,event_type)
                    file.write(line)

def main():
    logger = logging.getLogger('TCfinder')
    logger.setLevel(level=logging.INFO)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    args = misc.get_args()
    #read params:
    params_dict = read_params_dict(args)
    params_dict["cpu"] = args.cpu
    #handler:
    if args.log:
        if args.log == True:
            log_index = 1
            while os.path.exists("TBSfinder%s.log" % log_index):
                log_index += 1
            args.log = "TBSfinder%s.log"% log_index
        handler = logging.FileHandler(args.log)
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    ts_collections_list = load_tbs_objects(logger,args)

    #### normalization ####

    logger.info("start normalization...")
    logger.info("using core numbers:%s" %args.cpu)

    if args.normalization not in ["raw","simplecpm","powerlaw"]:
        logger.error("Unknown normalization method!Normalization methods supported are raw,simplecpm and powerlaw")
        raise Exception()
    elif args.normalization == "simplecpm":
        logger.info("using simplecpm to normalize the data.")
    elif args.normalization == "powerlaw":
        logger.info("using powerlaw to normalize the data.")
    else:
        logger.info("keep the raw counts of the data.")

    p = Pool(int(args.cpu))

    ts_tmp_collections_list = []

    for ts_collection_dict in ts_collections_list:
        if bool(ts_collection_dict):
            for ts_type,ts_collection in ts_collection_dict.items():
                ts_tmp_collections_list.append(ts_collection)

    ts_tmp_collections_list = p.map(partial(normalization,method = args.normalization),ts_tmp_collections_list)

    p.close()
    p.join()

    #update the ts_collections_list
    for i,ts_collection_dict in enumerate(ts_collections_list):
        if bool(ts_collection_dict):
            for ts_type,ts_collection in ts_collection_dict.items():
                ts_collections_list[i][ts_type] = ts_tmp_collections_list.pop(0)

    logger.info("\nNormalization done!\n")

    ### normalization done ###

    ### start clustering ###
    logger.info("start clustering...\n")
    logger.info("\nusing core numbers:%s\n" %args.cpu)
    if args.clustering not in ["paraclu","distclu"]:
        logger.error("Unknown clustering method!\ncluster methods supported are paraclu and distclu")
        raise Exception()
    elif args.clustering == "paraclu":
        logger.info("using paraclu to cluster the data.")
    elif args.clustering == "distclu":
        logger.info("using distclu to cluster the data.")

    ts_tmp_collections_list = []

    for ts_collection_dict in ts_collections_list:
        if bool(ts_collection_dict):
            for ts_type,ts_collection in ts_collection_dict.items():
                ts_tmp_collections_list.append(ts_collection)

    for i,ts_collection in enumerate(ts_tmp_collections_list):
        ts_collection.clustering(args.cpu,params_dict,args.clustering)
        ts_tmp_collections_list[i] = ts_collection

    logger.info("clustering done!")
    ### clustering done ###

    ### cluster filter ###
    if args.clustering == "paraclu":
        logger.info("start clustering filter...\n")

        for i,ts_collection in enumerate(ts_tmp_collections_list):
            ts_collection.cluster_filter(params_dict)
            ts_tmp_collections_list[i] = ts_collection

        logger.info("clustering filter done...\n")
    ### cluster filter done ###

    ### writing results ###
    logger.info("writing...\n")
    writing(args,ts_tmp_collections_list)
    logger.info("writing done!")
    ### writing done ###

if __name__ == "__main__":
	main()
