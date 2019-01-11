from utils import *
import glob
import os
import math
import operator

fhinput = open('../data/Census_allTue_May_23_12-08-15_2017.tsv')
cosmic_genes = []


def prep_file_paths(key):
    paths_raw = glob.glob(key)

    dict_paths = {}
    for filename in paths_raw:
          file_core = filename[filename.rindex('/')+1:]
          if file_core.startswith('cc'):
              num_genes = int(file_core.split('_')[1][1:])
              dict_paths[num_genes] = filename

    sorted_dict = sorted(dict_paths.items(), key=operator.itemgetter(0))

    paths = []
    for i in range(len(sorted_dict)):
        paths.append(sorted_dict[i][1])
    return paths

def cosmic_overlap(filename):
    with open(filename) as f:
        fhmodule = f.readlines()
        module_genes = []
        for line in fhmodule:
           words = line.strip().split()
           for word in words:
               if word not in module_genes:
                  module_genes.append(word)

    count = 0
    for gene in module_genes:
       if gene in cosmic_genes:
           count += 1
    return (count, len(module_genes))

def calculate_tpr_fpr(filename):
    (count_overlap, count) = cosmic_overlap(filename)
    tpr = count_overlap / float(len(cosmic_genes))
    sth_sth = 9859
    fpr = float(count - count_overlap) / ( sth_sth - len(cosmic_genes))
    return (tpr, fpr)

def calculate_all(hotnet_paths, our_paths):
    line = fhinput.readline()
    for line in fhinput:
        cosmic_genes.append(line.split()[0])

    with open("tpr_fpr_hotnet2.txt", "w+") as f:
        for path in hotnet_paths:
            (tpr, fpr) = calculate_tpr_fpr(path)
            f.write(str(tpr) + " " + str(fpr) + "\n")

    with open("tpr_fpr_our.txt", "w+") as f:
        for path in our_paths:
            (tpr, fpr) = calculate_tpr_fpr(path)
            f.write(str(tpr) + " " + str(fpr) + "\n")

    with open("count_overlap_hotnet2.txt", "w+") as f:
        for path in hotnet_paths:
            (cnt1, cnt2) = cosmic_overlap(path)
            f.write(str(cnt1) + "\n")

    with open("count_overlap_our.txt", "w+") as f:
        for path in our_paths:
            (cnt1, cnt2) = cosmic_overlap(path)
            f.write(str(cnt1) + "\n")


our_path_subdirs = ["m*g1cov*g2cov", "mutex_m*sqrt_g1cov*g2cov", "mutex_m*sqrt_g1cov*g2cov+d1", "mutex_m*sqrt_g1cov*g2cov+d2", "mutex_m*sqrt_g1cov*g2cov+d3", "log_mutex_m*sqrt_g1cov*g2cov", "log_mutex_m*sqrt_g1cov*g2cov+logd1"]

our_paths = prep_file_paths('../hint/out/our_subnetworks/' + our_path_subdirs[6] +'/*.txt')
hotnet_paths = prep_file_paths('../hotnet2_subnetworks_direction/*.txt')

calculate_all(hotnet_paths, our_paths)
