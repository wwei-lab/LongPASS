#functions used to normalize the sites counts.
from scipy import stats
from scipy.special import zeta
import math

def simpleTpm(site_counts,total_count):
    '''
    params:
    site_count:count of a site.
    total_count:total number of site count of different TBS object.
    return:
    the normalized result of the site count:cpm.
    '''
    try:
        cpm = [round((float(site_count)/total_count * (10**6)),5) for site_count in site_counts]
    except Exception as e:
        print (e,"Something went wrong while running the simpleTpm normalization,please check the reads count column of the input file.")
    return cpm

def fit_power_law_to_values(sites_counts,value_range = [10,1000]):
    '''
    param:
    sites_counts:list of dirrerent TBS obejct's site count.
    value_range:two number,one for the start of the range,one for the end of the range.This range is used to remove tags counts 
    that are out of the range of the input.
    return:
    the two coefficients describing the log form of the power law distribution that fitted to values.
    '''

    sites_counts_nums = [[x,sites_counts.count(x)] for x in sorted(list(set(sites_counts)) ,reverse = True)]
    '''
    sites_counts_nums:[[sites_counts:sites_counts_nums],...]
    '''
    cum = 0
    for sites_counts_item in sites_counts_nums:
        sites_counts_item[1] = sites_counts_item[1] + cum
        cum = sites_counts_item[1]
    
    # now drop the range of sites_counts we dont want.
    min_range = value_range[0]
    max_range = value_range[1]
    sites_counts_nums = [sites_counts_item for sites_counts_item in sites_counts_nums if max_range >= sites_counts_item[0] >= min_range]
    #fit to linear model
    x = [math.log(sites_counts_item[0]) for sites_counts_item in sites_counts_nums]
    y = [math.log(sites_counts_item[1]) for sites_counts_item in sites_counts_nums]
    slope, intercept, r, p, std_err = stats.linregress(x, y)

    return slope,intercept

def fit_values_to_referent_power_law(sites_counts,coefficients,alpha = 1.25,T = 10**6):
    '''
    param:
    sites_counts:list of dirrerent TBS obejct's site count.
    coefficients:output of function fit_power_law_to_values.
    slope of the referent power-law distribution.
    total number of sites_counts in the referent power-law distribution.
    return:
    the normalized sites counts.
    '''
    slope = coefficients[0]
    intercept = coefficients[1]
    para_lambda = (T/(zeta(alpha) * math.exp(intercept)))**(1/alpha)
    Beta = -1 * slope/alpha
    normalized_sites_counts = list(map(lambda x:round(para_lambda * (x**Beta),5),sites_counts))

    return normalized_sites_counts





