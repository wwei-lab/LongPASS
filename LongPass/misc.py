import argparse
from version import __version__
import multiprocessing
import os
import logging
from TPS import TPS_all_collection,TPS_chr_collection,TPS_strand_collection

__doc__ = ''' 
   ,--,                        ,----,                                       
,---.'|                      ,/   .`|              ,-.----.                 
|   | :    ,-.----.        ,`   .'  :   .--.--.    \    /  \     .--.--.    
:   : |    \    /  \     ;    ;     /  /  /    '.  |   :    \   /  /    '.  
|   ' :    ;   :    \  .'___,/    ,'  |  :  /`. /  |   |  .\ : |  :  /`. /  
;   ; '    |   | .\ :  |    :     |   ;  |  |--`   .   :  |: | ;  |  |--`   
'   | |__  .   : |: |  ;    |.';  ;   |  :  ;_     |   |   \ : |  :  ;_     
|   | :.'| |   |  \ :  `----'  |  |    \  \    `.  |   : .   /  \  \    `.  
'   :    ; |   : .  /      '   :  ;     `----.   \ ;   | |`-'    `----.   \ 
|   |  ./  ;   | |  \      |   |  '     __ \  \  | |   | ;       __ \  \  | 
;   : ;    |   | ;\  \     '   :  |    /  /`--'  / :   ' |      /  /`--'  / 
|   ,/     :   ' | \.'     ;   |.'    '--'.     /  :   : :     '--'.     /  
'---'      :   : :-'       '---'        `--'---'   |   | :       `--'---'   
           |   |.'                                 `---'.|                  
           `---'                                     `---`                  
'''

logger = logging.getLogger('TCfinder.misc')

def get_genestrand(ts_tag,read_strand):
    if ts_tag == "43":
        genestrand = read_strand
    elif ts_tag == "45":
        if read_strand == "+":
            genestrand = "-"
        elif read_strand == "-":
            genestrand = "+"

    return genestrand

def bedtotps(bedfile,writefile,writefile2):
    f2 = open(writefile,"w")
    f3 = open(writefile2,"w") # this is for pos bed file.
    chr_pos_strand_dict = {}
    with open(bedfile,"r") as f:
        for line in f:
            # sites convert to 1-base.
            # tps file is 1-base
            chrname = line.split('\t')[0].strip()
            startsite = int(line.split('\t')[1].strip()) + 1
            endsite = int(line.split('\t')[2].strip())
            readid = line.split('\t')[3].strip()
            ts_tag = line.split('\t')[4].strip()
            read_strand = line.split('\t')[5].strip()

            genestrand = get_genestrand(ts_tag,read_strand)

            if genestrand == "+":
                chr_pos_strand1 = "%s\t%s\t%s" %(chrname,startsite,genestrand)
                chr_pos_strand2 = "%s\t%s\t%s" %(chrname,endsite,genestrand)
                
                if chr_pos_strand1 not in chr_pos_strand_dict:
                    chr_pos_strand_dict[chr_pos_strand1] = [1,0]
                else:
                    chr_pos_strand_dict[chr_pos_strand1][0] += 1
                if chr_pos_strand2 not in chr_pos_strand_dict:
                    chr_pos_strand_dict[chr_pos_strand2] = [0,1]
                else:
                    chr_pos_strand_dict[chr_pos_strand2][1] += 1
            else:
                chr_pos_strand1 = "%s\t%s\t%s" %(chrname,endsite,genestrand)
                chr_pos_strand2 = "%s\t%s\t%s" %(chrname,startsite,genestrand)

                if chr_pos_strand1 not in chr_pos_strand_dict:
                    chr_pos_strand_dict[chr_pos_strand1] = [1,0]
                else:
                    chr_pos_strand_dict[chr_pos_strand1][0] += 1
         
                if chr_pos_strand2 not in chr_pos_strand_dict:
                    chr_pos_strand_dict[chr_pos_strand2] = [0,1]
                else:
                    chr_pos_strand_dict[chr_pos_strand2][1] += 1

            f3.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrname,startsite - 1,startsite,readid,"0",genestrand))
            f3.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrname,endsite,endsite + 1,readid,"0",genestrand))

    for chr_pos_strand,counts in chr_pos_strand_dict.items():
        writeline = "%s\t%s\n" %(chr_pos_strand,'\t'.join([str(count) for count in counts]))
        f2.write(writeline)

    f2.close()
    f3.close()

def read_gtf_transcript(gtf_path):
  all_transcript_list = []
  with open(gtf_path,"r") as gtf_file:
    for line in gtf_file:
      if line.split('\t')[2] == "transcript":
        gene_id = line.split('\t')[8].split(';')[0]
        chro = line.split('\t')[0]
        strand = line.split('\t')[6]
        start = int(line.split('\t')[3])
        end = int(line.split('\t')[4])
        width = end - start + 1
        all_transcript_list.append((chro,strand,start,end,width,gene_id))
  return all_transcript_list


def read_range_withoutstrand_file(range_file):
  all_ranges = []
  with open(range_file,'r') as rangefile:
    for line in rangefile:
      chro = line.split('\t')[0]
      start = int(line.split('\t')[1])
      end = int(line.split('\t')[2])
      width = end - start + 1
      all_ranges.append((chro,start,end,width))

  return all_ranges

def read_range_withstrand_file(range_file):
  all_ranges = []
  with open(range_file,'r') as rangefile:
    for line in rangefile:
      chro = line.split('\t')[0]
      start = int(line.split('\t')[1])
      end = int(line.split('\t')[2])
      strand = line.split(('\t')[3])

      all_ranges.append((chro,strand,start,end))

def filter_by_cover_radio(args,radio_threshold,all_ranges,all_gtf_transcript):
  #need a buffer length?
  #filter the ranges and correct the strand to the gtf file.
  all_ranges_filterBylength = []
  for range_ in all_ranges:
    range_chro,range_start,range_end,range_width = range_
    for gtf_transcript in all_gtf_transcript:
      gtf_chro,gtf_strand,gtf_start,gtf_end,gtf_width,gene_id = gtf_transcript
      buffer_length = 100
      if range_chro == gtf_chro and range_start >= gtf_start - buffer_length and range_end <= gtf_end + buffer_length and range_width / gtf_width >= radio_threshold:
        all_ranges_filterBylength.append((range_chro,gtf_strand,range_start,range_end))

  return all_ranges_filterBylength

def range_to_peak(all_ranges_filterBylength):
  #dict to store position count
  chro_strand_pos_dict = {}
  for range_ in all_ranges_filterBylength:
    chro,strand,start,end = range_
    chro_strand_start = chro + "_" + strand + "_" + str(start)
    chro_strand_end = chro + "_" + strand + "_" + str(end)
    if chro_strand_start in chro_strand_pos_dict:
      chro_strand_pos_dict[chro_strand_start][0] += 1
    else:
      chro_strand_pos_dict[chro_strand_start] = [1,0]
    if chro_strand_end in chro_strand_pos_dict:
      chro_strand_pos_dict[chro_strand_end][1] += 1
    else:
      chro_strand_pos_dict[chro_strand_end] = [0,1]

  # for pos,counts_list in chro_strand_pos_dict.items():
  #   if counts_list[0] < peak_threshold:
  #       chro_strand_pos_dict[pos][0] = 0
  #   if counts_list[1] < peak_threshold:
  #       chro_strand_pos_dict[pos][1] = 0

  # chro_strand_pos_dict = {i:chro_strand_pos_dict[i] for i in chro_strand_pos_dict if chro_strand_pos_dict[i] != [0,0]}

  return chro_strand_pos_dict


def get_args():
    parser = argparse.ArgumentParser(description="",
    )
    
    logger.info(__doc__)

    parser.add_argument("-v", "--version",
                         help="Print version and exit.",
                         action="version",
                         version='TPSfinder {}'.format(__version__),
                         default="")

    parser.add_argument("--params",
                          default="",
                          metavar="file",
    )

    parser.add_argument("--length_radio_threshold",
                          help="length radio for screen.",
                          default="0.85"
    )

    parser.add_argument("--peak_threshold",
                        help = "peak threshold.",
                        default = "5",
    )

    parser.add_argument("--clustering",
                          help="method used for clustering,default is paraclu",
                          default="paraclu",
    )
    
    parser.add_argument("--rangefile",
                        help="Data is range file.",
                        nargs="+",
                        metavar="file",
                        default=[]
    )    
    parser.add_argument("--tssfile",
                         help="Data is transcription start site file.",
                         nargs='+',
                         metavar="file",
                         default=[])

    parser.add_argument("--tesfile",
                         help="Data is transcription end site file.",
                         nargs='+',
                         metavar="file",
                         default=[])

    parser.add_argument("--tpsfile",
                         help="Data is transcription start and end site file.",
                         nargs='+',
                         metavar="file",
                         default=[])

    parser.add_argument("--bamfile",
                      help="mapping result of minimap",
                      nargs='+',
                      metavar="file",
                      default=[])

    parser.add_argument("--gtf",
                          help="annotation file",
                          nargs='+',
                          metavar="file",
                          default= []
    )

    parser.add_argument("--cpu",
                        help = "the number of cpu that will be used to process the data.",
                        default=multiprocessing.cpu_count(),
                        )

    parser.add_argument("--normalization",
                        help = "the methods used to normalize the count.",
                        default="raw",)

    parser.add_argument("--log",
                        help= "path for the output log file.",
                        action="store_true",
                        # default = "TPSfinder%s.log" % log_index,
                        )
              
    parser.add_argument("--replicate",
                        help="option for whether combine the samples."
                        ,action="store_true")

    out_index = 1
    while os.path.exists("out_cluster%s.txt" % out_index):
        out_index += 1
    parser.add_argument("-o",
                    help= "path for the output cluster file.",
                    metavar="file",
                    default = "out_cluster%s.txt" % out_index,
                    )
    
    args = parser.parse_args()
    return args

def load_lrtsp_objects(filepath,boundary_type = "tps"):
    from TPS import TPS
    logger = logging.getLogger('logger_main.logging_core')
    '''
    expected input file:
    chrname pos strand  count
    '''
    lrtsp_objects_positive_strand = {}
    lrtsp_objects_minus_strand = {}
    lrtsp_objects = {
      "+" : lrtsp_objects_positive_strand,
      "-" : lrtsp_objects_minus_strand,
    }

    TPS_all_list = []
    if boundary_type != "tps":
      with open(filepath,'r') as lrtspfile:
          for line in lrtspfile:
              if len(line.split('\t')) > 4:
                logger.warn("extra columns(>4) found in the file,use the fourth column as the count for tss or tes.")
              elif len(line.split('\t')) <= 3:
                logger.error("uncomplete file!please check your input!")
                raise Exception()
              chrname = line.split('\t')[0]
              pos = line.split('\t')[1]
              strand = line.split('\t')[2]
              count = line.split('\t')[3]
              lrtsp_object = TPS(chrname, pos, strand, count, boundary_type)
              TPS_all_list.append(lrtsp_object)
              if strand == '+':
                if chrname not in lrtsp_objects_positive_strand:
                  lrtsp_objects_positive_strand[chrname] = []
                  lrtsp_objects_positive_strand[chrname].append(lrtsp_object)
                else:
                  lrtsp_objects_positive_strand[chrname].append(lrtsp_object)
              elif strand == '-':
                if chrname not in lrtsp_objects_minus_strand:
                  lrtsp_objects_minus_strand[chrname] = []
                  lrtsp_objects_minus_strand[chrname].append(lrtsp_object)
                else:
                  lrtsp_objects_minus_strand[chrname].append(lrtsp_object)
              else:
                logger.warning("non-standard strand name found in the file,strand name is %s .still loading..." %strand)
                if strand not in lrtsp_objects:
                  lrtsp_objects[strand] = {}
                  lrtsp_objects[strand][chrname] = []
                else:
                  if chrname in lrtsp_objects[strand]:
                    lrtsp_objects[strand][chrname].append(lrtsp_object)
                  else:
                    lrtsp_objects[strand][chrname] = []
                    lrtsp_objects[strand][chrname].append(lrtsp_object)

    else:
      with open(filepath,'r') as tpsfile:
        lrtsp_object_tss = None
        lrtsp_object_tes = None
        for line in tpsfile:
          if len(line.split('\t')) <5:
            logger.error("uncomplete file!please check!")
            raise Exception()
          elif len(line.split('\t')) >5:
            logger.warn("extra columns(>5) found in the file,use the fourth and fifth columns for the tss and tes counts.")
          chrname = line.split('\t')[0]
          pos = line.split('\t')[1]
          strand = line.split('\t')[2]
          startcount = int(line.split('\t')[3])
          endcount = int(line.split('\t')[4])

          if startcount != 0 and endcount != 0:
            lrtsp_object_tss = TPS(chrname, pos, strand, startcount, "tss")
            lrtsp_object_tes = TPS(chrname, pos, strand, endcount, "tes")
          elif startcount != 0:
            lrtsp_object_tss = TPS(chrname, pos, strand, startcount, "tss")
          elif endcount != 0:
            lrtsp_object_tes = TPS(chrname, pos, strand, endcount, "tes")
          else:
            logger.warn("start site and end site counts are zero at position %s" %pos)
          
          lrtsp_object = [lrtsp_object for lrtsp_object in [lrtsp_object_tss,lrtsp_object_tes] if lrtsp_object is not None]


          if strand == '+':
            if chrname not in lrtsp_objects_positive_strand:
              lrtsp_objects_positive_strand[chrname] = []
              lrtsp_objects_positive_strand[chrname].extend(lrtsp_object)
            else:
              lrtsp_objects_positive_strand[chrname].extend(lrtsp_object)
          elif strand == '-':
            if chrname not in lrtsp_objects_minus_strand:
              lrtsp_objects_minus_strand[chrname] = []
              lrtsp_objects_minus_strand[chrname].extend(lrtsp_object)
            else:
              lrtsp_objects_minus_strand[chrname].extend(lrtsp_object)
          else:
            logger.warning("non-standard strand name found in the file,strand name is %s .still loading..." %strand)
            if strand not in lrtsp_objects:
              lrtsp_objects[strand] = {}
              lrtsp_objects[strand][chrname] = []
            else:
              if chrname in lrtsp_objects[strand]:
                lrtsp_objects[strand][chrname].extend(lrtsp_object)
              else:
                lrtsp_objects[strand][chrname] = []
                lrtsp_objects[strand][chrname].extend(lrtsp_object)
          
          lrtsp_object_tss = None
          lrtsp_object_tes = None

    TPS_strand_collection_list = []

    for strand,lrtsp_objects_strand in lrtsp_objects.items():
      TPS_chr_collection_list = []
      for chro,TPS_chr_list in lrtsp_objects_strand.items():
        TPS_chr_collection_list.append(TPS_chr_collection(chro, TPS_chr_list))
        TPS_all_list.extend(TPS_chr_list)
      TPS_strand_collection_list.append(TPS_strand_collection(strand,TPS_chr_collection_list))
    
    return TPS_all_collection(TPS_all_list,TPS_strand_collection_list)

  