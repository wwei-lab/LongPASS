######################################
#  classes about cluster
######################################


class dist_cluster(object):
    def __init__(self,start,end,chro,strand,eventtype):
        self.chro = chro
        self.strand = strand
        self.start = start
        self.end = end
        self.eventtype = eventtype      

class para_cluster(object):
    def __init__(self,start,end,num_sites,dominant_site,total_cpm,dominant_site_cpm,min_d,max_d,chro,strand,eventtype):
        self.chro = chro
        self.strand = strand
        self.start = start
        self.end = end
        self.num_sites = num_sites
        self.dominant_site = dominant_site
        self.total_cpm = total_cpm
        self.dominant_site_cpm = dominant_site_cpm
        self.min_d = min_d
        self.max_d = max_d
        self.eventtype = eventtype
