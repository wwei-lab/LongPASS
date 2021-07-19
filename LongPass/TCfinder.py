"""
Copyright 2021 Yanhong Hong (a524664992@gmail.com)

This module is the main entry of TCfinder.It is executed when user runs `run_TCfinder.py` 
(directly from the source directory).
"""

from misc import load_lrtsp_objects,bedtotbs,get_genestrand
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
import pysam
import pybedtools
import csv

def loadbamfile(bamfilepath):
    #bam file filter:keep lines that have ts tags.
    readbamfile = pysam.AlignmentFile(bamfilepath,"rb")
    writebamfilepath = "./tagfilter_" + os.path.basename(bamfilepath)
    writebamfilepath2 = "./flagfilter_tagfilter_" + os.path.basename(bamfilepath)
    writebamfile = pysam.AlignmentFile(writebamfilepath,"wb",template  = readbamfile)
    writebamfile2 = pysam.AlignmentFile(writebamfilepath2,"wb",template  = readbamfile)
    for read in readbamfile:
        if read.has_tag("ts"):
            writebamfile.write(read)
            if int(read.flag) < 2048:
                writebamfile2.write(read)

    readbamfile.close()
    writebamfile.close()
    writebamfile2.close()

    #bam to bed and extract ts tag information
    if writebamfilepath2.endswith(".bam"):
        writebedfilepath = writebamfilepath2.replace(".bam",".bed")
    else:
        writebedfilepath = writebamfilepath2 + ".bed"
    bedfile = pybedtools.BedTool(writebamfilepath2).bamtobed(tag = "ts").saveas(writebedfilepath)

    #bed file to lrtsp file
    writelrtsppath = writebedfilepath.replace(".bed", ".lrtsp")
    writebedfilepath2 = "genestrand_" + os.path.basename(writebedfilepath)  # path for writing bed file that convert to genestrand
    bedtotbs(writebedfilepath, writelrtsppath, writebedfilepath2)

    return [writebedfilepath,writelrtsppath,writebedfilepath2]

def noreplicate_peak_filter(filepath,peak_threshold,logger,filetype):
    # script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/misc/combine_tbsfile.sh"
    # filter_command = "cat " + filepath + "|sort -k1,1 -k2n,2 > " + filepathtmp

    # try:
    #     filter_status = subprocess.check_call(filter_command,shell=True)
    # except:
    #     pass
    filepathtmp = filepath + ".tmp"
    finalfilepath = "final_" + os.path.basename(filepath)

    combine_files([filepath], filepathtmp, logger, columns_num=2)

    # if filetype == "tss" or filetype == "tes":
    #     filter_command2 = "bash " + script_path + " " + filepathtmp + " " + str(peak_threshold) + " |awk  '{print $1,$2,$3,$4}'|sed 's/ /\t/g'|sort -k1,1 -k2n,2 > " + finalfilepath
    # else:
    #     filter_command2 = "bash " + script_path + " " + filepathtmp + " " + str(peak_threshold) + " |sort -k1,1 -k2n,2 > " + finalfilepath
    # try:
    #     filter_status2 = subprocess.check_call(filter_command2,shell = True)
    # except:
    #     pass

    peak_filter(filepathtmp, filetype, finalfilepath, peak_threshold)
    os.remove(filepathtmp)

    return finalfilepath

def peak_filter(all_combine_tbsfiles,tbs_type,write_file,peak_threshold):
    peak_threshold = int(peak_threshold)
    with open(all_combine_tbsfiles) as f:
        tss_dict = {}
        tes_dict = {}
        for line in f:
            line = line.strip()
            chro = line.split('\t')[0]
            pos = line.split('\t')[1]
            strand = line.split('\t')[2]
            chro_pos_strand = chro + "\t" + pos + "\t" + strand
            if tbs_type == "tss":
                tss = int(line.split('\t')[3])
                if chro_pos_strand not in tss_dict:
                    tss_dict[chro_pos_strand] = tss
                else:
                    tss_dict[chro_pos_strand] += tss
            elif tbs_type == "tes":
                tes = int(line.split('\t')[3])
                if chro_pos_strand not in tes_dict:
                    tes_dict[chro_pos_strand] = tes
                else:
                    tes_dict[chro_pos_strand] += tes
            elif tbs_type == "tbs":
                tss = int(line.split('\t')[3])
                tes = int(line.split('\t')[4])
                if chro_pos_strand not in tss_dict:
                    tss_dict[chro_pos_strand] = tss
                else:
                    tss_dict[chro_pos_strand] += tss

                if chro_pos_strand not in tes_dict:
                    tes_dict[chro_pos_strand] = tes
                else:
                    tes_dict[chro_pos_strand] += tes

        if tbs_type == "tss":
            for chro_pos_strand,count in tss_dict.items():
                if count < peak_threshold:
                    tss_dict[chro_pos_strand] = 0
            
            with open(write_file , "w") as f:
                for chro_pos_strand,count in tss_dict.items():
                    if count != 0:
                        f.write("%s\t%s\n"%(chro_pos_strand,count))

        if tbs_type == "tes":
            for chro_pos_strand,count in tes_dict.items():
                if count < peak_threshold:
                    tes_dict[chro_pos_strand] = 0

            with open(write_file , "w") as f:
                for chro_pos_strand,count in tes_dict.items():
                    if count != 0:
                        f.write("%s\t%s\n"%(chro_pos_strand,count))

        if tbs_type == "tbs":
            for chro_pos_strand,count in tss_dict.items():
                if count < peak_threshold:
                    tss_dict[chro_pos_strand] = 0

            for chro_pos_strand,count in tes_dict.items():
                if count < peak_threshold:
                    tes_dict[chro_pos_strand] = 0

            with open(write_file , "w") as f:
                for chro_pos_strand in tss_dict:
                    tss = tss_dict[chro_pos_strand]
                    tes = tes_dict[chro_pos_strand]
                    if (tss != 0 or tes != 0):
                        f.write("%s\t%s\t%s\n"%(chro_pos_strand,tss,tes))


def load_tbs_objects(peak_threshold,logger,args):
    logger.info("loading files...\n")

    tssfiles = args.tssfile
    tesfiles = args.tesfile
    tbsfiles = args.tbsfile
    rangefiles = args.rangefile
    bamfiles = args.bamfile

    tss_collections_dict = {}
    tes_collections_dict = {}
    tbs_collections_dict = {}

    #input file type:rangefile
    # range file needs annotation file to find the transcript strand.
    '''
    rangefile:
    chrIS   38705   96443
    chrIS   44260   96459
    '''
    if len(rangefiles) > 0:
        all_transcript_list = []
        if (not args.gtf) or len(args.gtf) > 1:
            logger.error("no gtf file or more than one gtf file provided!") # range file must provide gtf file.
            raise Exception()
        else:
            all_transcript_list = misc.read_gtf_transcript(args.gtf[0])

        for rangefile in rangefiles:
            rangefile_basename = os.path.basename(rangefile)
            write_tbs_filename = os.path.splitext(rangefile_basename)[0] + ".lrtsp"
            logger.info("\t>>>reading %s\n" % os.path.basename(rangefile))
            if args.gtf:
                all_ranges = misc.read_range_withoutstrand_file(rangefile)
                all_ranges_filterBylength = misc.filter_by_cover_radio(args,float(args.length_radio_threshold),all_ranges,all_transcript_list)
                all_ranges_ToPeak = misc.range_to_peak(all_ranges_filterBylength)
            else:
                all_ranges = misc.read_range_withstrand_file(rangefile)
                all_ranges_ToPeak = misc.range_to_peak(all_ranges)
            with open(write_tbs_filename,'w') as f:
                for chro_strand_pos,counts in all_ranges_ToPeak.items():
                    chro = chro_strand_pos.split('_')[0]
                    strand = chro_strand_pos.split('_')[1]
                    pos = chro_strand_pos.split('_')[2]
                    if strand == "+":
                        write_line = "%s\t%s\t%s\t%s\t%s\n" % (chro,pos,strand,counts[0],counts[1])
                    elif strand == "-":
                        write_line = "%s\t%s\t%s\t%s\t%s\n" % (chro,pos,strand,counts[1],counts[0])
                    f.write(write_line)
            tbsfiles.append(write_tbs_filename)

    #input file type:bamfile
    if len(bamfiles) > 0:
        logger.info("coverting bam file to lrtsp file")
        p = Pool(int(args.cpu))
        bed_tbsfiles = p.map(loadbamfile,bamfiles)
        bedfiles = [bed_tbsfile[0] for bed_tbsfile in bed_tbsfiles]
        tbsfiles = [bed_tbsfile[1] for bed_tbsfile in bed_tbsfiles]
        genestrandbedfiles = [bed_tbsfile[2] for bed_tbsfile in bed_tbsfiles]
        p.close()
        p.join()

    #input file type:tssfile
    if len(tssfiles) > 0:
        if not args.replicate:
            for tssfilepath in tssfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tssfilepath))
                final_tssfile_path = noreplicate_peak_filter(tssfilepath,peak_threshold,logger,"tss")
                tss_collection = misc.load_lrtsp_objects(final_tssfile_path, "tss")
                tss_collections_dict[os.path.basename(final_tssfile_path)] = tss_collection
        else:
            #replicate case,peak filter
            all_combine_tssfiles = "all_combine_tssfiles.lrtsp"
            unique_all_combine_tssfiles = "unique_all_combine_tssfiles.lrtsp"
            tss_collection = replicate_case_getcollection(tssfiles, all_combine_tssfiles, unique_all_combine_tssfiles, logger,"tss",peak_threshold)
            tss_collections_dict[unique_all_combine_tssfiles] = tss_collection
            try:
                os.remove(all_combine_tssfiles)
                os.remove(unique_all_combine_tssfiles)
            except:
                pass

    #input file type:tesfile
    if len(tesfiles) > 0:
        if not args.replicate:
            for tesfilepath in tesfiles:
                logger.info("\t>>>reading %s\n" % os.path.basename(tesfilepath))
                final_tesfile_path = noreplicate_peak_filter(tesfilepath,peak_threshold,logger,"tes")
                tes_collection = misc.load_lrtsp_objects(final_tesfile_path, "tes")
                tes_collections_dict[os.path.basename(final_tesfile_path)] = tes_collection
        else:
            #replicate case,peak filter
            all_combine_tesfiles = "all_combine_tesfiles.lrtsp"
            unique_all_combine_tesfiles = "unique_all_combine_tesfiles.lrtsp"
            tes_collection = replicate_case_getcollection(tesfiles, all_combine_tesfiles, unique_all_combine_tesfiles, logger,"tes",peak_threshold)
            tes_collections_dict[unique_all_combine_tesfiles] = tes_collection
            try:
                os.remove(all_combine_tesfiles)
                os.remove(unique_all_combine_tesfiles)
            except:
                pass

    #input file type:tbs file
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
                final_tbsfile_path = noreplicate_peak_filter(tbsfilepath,peak_threshold,logger,"tbs")
                tbs_collection = misc.load_lrtsp_objects(final_tbsfile_path, "tbs")
                tbs_collections_dict[os.path.basename(final_tbsfile_path)] = tbs_collection
        else:
            all_combine_tbsfiles = "all_combine_tbsfiles.lrtsp"
            unique_all_combine_tbsfiles = "unique_all_combine_tbsfiles.lrtsp"
            tbs_collection = replicate_case_getcollection(tbsfiles, all_combine_tbsfiles, unique_all_combine_tbsfiles, logger,"tbs",peak_threshold)
            tbs_collections_dict[unique_all_combine_tbsfiles] = tbs_collection
            try:
                os.remove(all_combine_tbsfiles)
                os.remove(unique_all_combine_tbsfiles)
            except:
                pass

    return [tss_collections_dict,tes_collections_dict,tbs_collections_dict,bedfiles,genestrandbedfiles]

def replicate_case_getcollection(tbsfiles,all_combine_tbsfiles,unique_all_combine_tbsfiles,logger,tbs_type,peak_threshold):
    combine_files(tbsfiles,all_combine_tbsfiles,logger,columns_num=2)
    peak_filter(all_combine_tbsfiles, tbs_type, unique_all_combine_tbsfiles, peak_threshold)
    tbs_collection = misc.load_lrtsp_objects(unique_all_combine_tbsfiles, tbs_type)
    return tbs_collection

def non_replicate_case_getcollection():
    pass

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
            "peak_threshold" : 10,
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

def reads_assign_to_cluster(cluster_args_combinegenestrand):
    cluster, args, combine_genestrand = cluster_args_combinegenestrand
    pos_bed = ""
    start = cluster.start - 1
    end = cluster.end
    strand = cluster.strand
    chro = cluster.chro
    num_sites = cluster.num_sites
    dominant_site = cluster.dominant_site - 1
    total_cpm = cluster.total_cpm
    dominant_site_cpm = cluster.dominant_site_cpm
    max_d = cluster.max_d
    min_d = cluster.min_d
    event_type = cluster.eventtype
    all_pos = sorted(cluster.all_pos)
    if len(args.bamfile) > 0:

        combine_bed_file = pybedtools.BedTool(combine_genestrand)
        for pos in all_pos:
            pos_bed += chro + "\t" + str(pos - 1) + "\t" + str(pos) + "\t0\t0\t" + strand + "\n" # 1-base to 0-base
        pos_bed = pos_bed.strip('\n')
        pos_bed = pybedtools.BedTool(pos_bed,from_string= True).sort()

        overlaps = combine_bed_file.intersect(pos_bed,sorted = True,s = True)
        print("done")
        readsids = [overlap.name for overlap in overlaps]

        readsidsstr = ",".join(readsids)

    line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chro,strand,start,end,num_sites,dominant_site,dominant_site_cpm,round(total_cpm,5),round(max_d,5),round(min_d,5),event_type,readsidsstr)
    
    with open(args.o,"a") as f:
        f.write(line)

def combine_files(files_list,write_combine_file,logger,columns_num):
    lists = []
    for file_ in files_list:
        logger.info("\t>>>reading %s\n" % os.path.basename(file_))
        with open(file_,"r") as f:
            for line in f:
                lists.append(line.rstrip().split('\t'))
    
    if columns_num == 1:
        sorted_results = sorted(lists, key=lambda x:(x[0]))
    
    if columns_num == 2:
        sorted_results = sorted(lists, key=lambda x:(x[0], int(x[1])))

    if columns_num == 3:
        sorted_results = sorted(lists, key=lambda x:(x[0], int(x[1]), int(x[2])))

    with open(write_combine_file,"w") as f2:
        writer = csv.writer(f2,delimiter = '\t')
        writer.writerows(sorted_results)

def writing(args,ts_tmp_collections_list,genestrandbedfiles,logger):
    if args.clustering == "paraclu":
        combine_genestrand = "combine_genestrand.bed"
        combine_files(genestrandbedfiles,combine_genestrand,logger,3)

        ts_collection_all_clusters_list = ts_tmp_collections_list[0].all_cluster_list
        for ts_collection_clusters in ts_collection_all_clusters_list:
            p = Pool(int(args.cpu))
            all_cluster_args = []
            for cluster in ts_collection_clusters:
                cluster_args_combinegenestrand = [cluster,args,combine_genestrand]
                all_cluster_args.append(cluster_args_combinegenestrand)

            p.map(reads_assign_to_cluster,all_cluster_args)
            p.close()
            p.join()

    else:
        with open(args.o,'w') as file:
            ts_collection_all_clusters_list = ts_tmp_collections_list[0].all_cluster_list
            for ts_collection_clusters in ts_collection_all_clusters_list:
                for cluster in ts_collection_clusters:
                    start = cluster.start - 1
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

    #load lrtsp file
    peak_threshold = params_dict["peak_threshold"]
    ts_collections_list = load_tbs_objects(peak_threshold,logger,args)
    genestrandbedfiles = ts_collections_list[-1]
    bedfiles = ts_collections_list[-2]
    ts_collections_list = ts_collections_list[:-2]

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
    writing(args,ts_tmp_collections_list,genestrandbedfiles,logger)
    logger.info("writing done!")
    ### writing done ###


if __name__ == "__main__":
	main()