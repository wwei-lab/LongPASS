#function used to cluster all boundary sites.
from cluster import dist_cluster
from cluster import para_cluster
import math

def distclu(strand_chro_pos_cpm_para_cluster_list_eventtype,max_dist):
    '''
    param:
    pos:a list of positions of sites
    max_dist:The max distance for the cluster algorithum
    return:
    clusters:list of cluster of sites.
    '''
    strand,chro,pos,cpm,para_cluster_list,event_type = strand_chro_pos_cpm_para_cluster_list_eventtype
    pos.sort()
    distances = list(map(lambda x,y:x-y-1,pos[1:],pos[:-1]))
    #dist indexed that are more than the max distances
    dist_indexes = [x[0] for x in enumerate(distances) if x[1] >= max_dist]
    ends = [pos[dist_index] for dist_index in dist_indexes]
    ends.append(pos[-1])
    starts = [pos[dist_index + 1]-1 for dist_index in dist_indexes]
    starts.insert(0,pos[0]-1)

    clusters = [dist_cluster(start, end,chro, strand,  event_type) for start,end in zip(starts,ends)]

    return clusters

def paraclu_findbreak(pos,cpm):
    '''
    param:
    pos:a list of positions of sites.
    cpm:a list of cpm for different sites.
    return:
    
    '''
    num_sites = len(pos)
    total_cpm = math.fsum(cpm)
    if num_sites == 1:
        minimal_density = math.inf
        break_site = None
        start_site = pos[0]
        end_site = start_site
        return (break_site,minimal_density,num_sites,total_cpm,start_site,end_site)
    pos_cpm_list = [list(a) for a in zip(pos, cpm)]
    pos_cpm_list = sorted(pos_cpm_list, key=lambda x: x[0])
    start_site = pos_cpm_list[0][0]
    end_site = pos_cpm_list[-1][0]
    prefix_cpm_cum = 0
    prefix_cpm_cum_list = []
    for pos_cpm in pos_cpm_list:
        cpm_single = pos_cpm[1]
        prefix_cpm_cum = prefix_cpm_cum + cpm_single
        prefix_cpm_cum_list.append(prefix_cpm_cum)
    prefix_cpm_cum_list = prefix_cpm_cum_list[:-1]
    prefix_segment_lengths_list = list(map(lambda x:x-pos_cpm_list[0][0],[i[0] for i in pos_cpm_list][1:]))
    prefix_segment_densities_list = [x/y for x, y in zip(prefix_cpm_cum_list, prefix_segment_lengths_list)]

    suffix_cpm_cum = 0
    suffix_cpm_cum_list = []
    for pos_cpm in list(reversed(pos_cpm_list)):
        cpm_single = pos_cpm[1]
        suffix_cpm_cum = suffix_cpm_cum + cpm_single
        suffix_cpm_cum_list.append(suffix_cpm_cum)
    suffix_cpm_cum_list = suffix_cpm_cum_list[:-1]
    suffix_segment_lengths_list = list(reversed(list(map(lambda x:pos_cpm_list[-1][0]-x,[pos[0] for pos in pos_cpm_list][:-1]  ) ) ) )
    suffix_segment_densities_list = [x/y for x, y in zip(suffix_cpm_cum_list, suffix_segment_lengths_list)]

    prefix_segment_minimal_density = min(prefix_segment_densities_list)
    prefix_segment_minimal_density_index = prefix_segment_densities_list.index(min(prefix_segment_densities_list))
    suffix_segment_minimal_density = min(suffix_segment_densities_list)
    suffix_segment_minimal_density_index = suffix_segment_densities_list.index(min(suffix_segment_densities_list))

    if prefix_segment_minimal_density >= suffix_segment_minimal_density:
        minimal_density = suffix_segment_minimal_density
        break_site = len(pos) - 1 - suffix_segment_minimal_density_index
    elif prefix_segment_minimal_density < suffix_segment_minimal_density:
        minimal_density = prefix_segment_minimal_density
        break_site = prefix_segment_minimal_density_index + 1
    else:
        raise Exception("sorry!something went wrong while finding the break site! probably because of the data type")

    return [break_site,minimal_density,num_sites,total_cpm,start_site,end_site]

def paraclu(strand_chro_pos_cpm_para_cluster_list_eventtype,minimal_density = -math.inf):
    strand,chro,pos,cpm,para_cluster_list,event_type = strand_chro_pos_cpm_para_cluster_list_eventtype
    if (len(pos) > 0):
        pos_cpm_list = [list(a) for a in zip(pos, cpm)]
        pos_cpm_list = sorted(pos_cpm_list, key=lambda x: x[0])
        pos = [i[0] for i in pos_cpm_list[:]]
        cpm = [i[1] for i in pos_cpm_list[:]]
        params = paraclu_findbreak(pos,cpm)
        break_site = params[0]
        maximal_density = params[1]
        num_sites = params[2]
        total_cpm = params[3]
        start_site = params[4]
        end_site = params[5]
        max_cpm = max(cpm)
        max_cpm_index = cpm.index(max_cpm)
        dominant_site = pos[max_cpm_index]

        if (maximal_density != math.inf):
            next_min = max(minimal_density,maximal_density)
            paraclu((strand,chro,pos[:break_site],cpm[:break_site],para_cluster_list,event_type),next_min)
            paraclu((strand,chro,pos[break_site:],cpm[break_site:],para_cluster_list,event_type),next_min)

        para_cluster_list.append(para_cluster(start_site,end_site,num_sites,dominant_site,total_cpm,max_cpm,minimal_density,maximal_density,chro,strand,event_type,pos))
        # the assignment of the chromesome and strand information is placed in the entry of the pipeline.
 
        return para_cluster_list
    else:
        return