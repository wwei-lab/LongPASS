import argparse
from version import __version__
import multiprocessing
import os
import logging
from TBS import TBS_all_collection,TBS_chr_collection,TBS_strand_collection

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

def get_args():
    parser = argparse.ArgumentParser(description="",
    )
    
    logger.info(__doc__)

    parser.add_argument("-v", "--version",
                         help="Print version and exit.",
                         action="version",
                         version='TBSfinder {}'.format(__version__),
                         default="")

    parser.add_argument("--params",
                          default="",
                          metavar="file",
    )

    parser.add_argument("--clustering",
                          help="method used for clustering,default is paraclu",
                          default="paraclu",
    )
    
    parser.add_argument("--tssfile",
                         help="Data is transcription start site file.",
                         nargs='+',
                         metavar="file",
                         default="")

    parser.add_argument("--tesfile",
                         help="Data is transcription end site file.",
                         nargs='+',
                         metavar="file",
                         default="")

    parser.add_argument("--tbsfile",
                         help="Data is transcription start and end site file.",
                         nargs='+',
                         metavar="file",
                         default="")

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
                        # default = "TBSfinder%s.log" % log_index,
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

def load_lrtsp_objects(filepath,boundary_type = "tbs"):
    from TBS import TBS
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

    TBS_all_list = []
    if boundary_type != "tbs":
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
              lrtsp_object = TBS(chrname, pos, strand, count, boundary_type)
              TBS_all_list.append(lrtsp_object)
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
      with open(filepath,'r') as tbsfile:
        lrtsp_object_tss = None
        lrtsp_object_tes = None
        for line in tbsfile:
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
            lrtsp_object_tss = TBS(chrname, pos, strand, startcount, "tss")
            lrtsp_object_tes = TBS(chrname, pos, strand, endcount, "tes")
          elif startcount != 0:
            lrtsp_object_tss = TBS(chrname, pos, strand, startcount, "tss")
          elif endcount != 0:
            lrtsp_object_tes = TBS(chrname, pos, strand, endcount, "tes")
          else:
            logger.warn("start site and end site counts are zero at position %s" %pos)
          
          lrtsp_object = [lrtsp_object for lrtsp_object in [lrtsp_object_tss,lrtsp_object_tes] if lrtsp_object is not None]
          TBS_all_list.extend(lrtsp_object)

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

    TBS_strand_collection_list = []

    for strand,lrtsp_objects_strand in lrtsp_objects.items():
      TBS_chr_collection_list = []
      for chro,TBS_chr_list in lrtsp_objects_strand.items():
        TBS_chr_collection_list.append(TBS_chr_collection(chro, TBS_chr_list))
      TBS_strand_collection_list.append(TBS_strand_collection(strand,TBS_chr_collection_list))
    
    return TBS_all_collection(TBS_all_list,TBS_strand_collection_list)