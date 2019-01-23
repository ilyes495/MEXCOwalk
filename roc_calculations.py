from utils import *
from glob import glob
import os
import math
import operator
from sklearn.metrics import auc
from matplotlib import pyplot as plt
import operator
from pylab import *

plt.xlabel('FP rate')
plt.ylabel('TP rate')
axes(frameon=0)
save_path = '../results/roc_old_cosmic/'
# new cosmic file: cancer_gene_census_2019_01_11_genes.txt
# old cosmic file: Census_allTue_May_23_12-08-15_2017.tsv
with open('../data/Census_allTue_May_23_12-08-15_2017.tsv', 'r') as f:

    # for old comsic file
    cosmic_genes = [line.split('\t')[0] for line in f.readlines()[1:]]

    # for new comsic file
    #cosmic_genes = [line.rstrip() for line in f.readlines()[:]]
    print('There are {} cosmins genes'.format(len(cosmic_genes)))
model2area = {}
# fhinput = open('../data/Census_allTue_May_23_12-08-15_2017.tsv')
# cosmic_genes = []
# lines = fhinput.readlines()
# for line in fhinput:
#     cosmic_genes.append(line.split()[0])

def prep_file_paths(key):
    paths_raw = glob(key)

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
    #print('iside cosmic: ', filename)

    with open(filename, 'r') as f:
        model_gene = []
        for line in f.readlines():
            model_gene.extend(line.rstrip().split())
        model_gene =set(model_gene)
        #print('len model_gene ', len(model_gene))
        #print('some genes: ', model_gene[:10])

    count_overlap = set(cosmic_genes).intersection(model_gene)

    return len(count_overlap), len(model_gene)

def calculate_tpr_fpr(filename):
    (count_overlap, count) = cosmic_overlap(filename)
    tpr = count_overlap / float(len(cosmic_genes))
    sth_sth = 9859
    fpr = float(count - count_overlap) / ( sth_sth - len(cosmic_genes))
    return (tpr, fpr)

def calculate_all(our_paths, key):
    tprs, fprs = [],[]
    #tprs_hotnet2, fprs_hotnet2 = [],[]

    # with open("../results/roc/tpr_fpr_hotnet2.txt", "w+") as f:
    #     for path in hotnet_paths[:len(our_paths)]:
    #         (tpr, fpr) = calculate_tpr_fpr(path)
    #         f.write(str(tpr) + " " + str(fpr) + "\n")
    #         tprs_hotnet2.append(tpr)
    #         fprs_hotnet2.append(fpr)

    with open(save_path+"fps_tps/tpr_fpr_our_"+key +".txt", "w+") as f:
        for path in our_paths:
            (tpr, fpr) = calculate_tpr_fpr(path)
            f.write(str(tpr) + " " + str(fpr) + "\n")
            tprs.append(tpr)
            fprs.append(fpr)

    # with open("../results/roc/count_overlap_hotnet2.txt", "w+") as f:
    #     for path in hotnet_paths:
    #         (cnt1, cnt2) = cosmic_overlap(path)
    #         f.write(str(cnt1) + "\n")

    with open(save_path+"fps_tps/count_overlap_our.txt", "w+") as f:
        for path in our_paths:
            (cnt1, cnt2) = cosmic_overlap(path)
            f.write(str(cnt1) + "\n")


    area =auc(fprs,tprs )
    model2area[key] = area
    print('\nArea under ROC {}: {}'.format(key, area))

    # plt.xlabel('FP rate')
    # plt.ylabel('TP rate')
    # plt.title('ROC for hotnet2 ')
    # plt.plot(fprs_hotnet2, tprs_hotnet2)
    # plt.show()
    # plt.savefig('../results/roc/roc_hotnet2.png')
    #plt.title('ROC for {}'.format(key))
    plt.plot(fprs, tprs, '-o')


cc_path = '../hint/out/connected_components_isolarge_n2500_whh/'
#our_path_subdirs = glob('../hint/out/connected_components_isolarge/*')
models  = ["hotnet2","memcover_v1","memcover_v2","mutex_t07_nsep_cov_nsep","mutex_t07_nsep_cov", ]

for model in models:
    key = model
    #if not key in models: continue
    # if key in ['hotnet2', 'memcover0.08','memcover0.2', 'memcover', 'memcover_basic', 'mutex_t07_nsep_cov_nsep' ]:
    #     continue
    #model = '../hint/out/connectrd_components_isolarge_v2_n2500/{}'.format(key)
    our_paths = prep_file_paths(cc_path+model +'/*.txt')
    #hotnet_paths = prep_file_paths('../hint/out/connected_components_isolarge/hotnet2/*.txt')

    # print('path is ', our_path_subdirs[0])
    # print('Filename is ', our_paths[0])
    calculate_all( our_paths, key)

model2area = sorted(model2area.items(), key=operator.itemgetter(1), reverse=True)
with open(save_path+'areas.txt', 'w') as f:
    for c in model2area:
        model, area = c
        print('{} {}'.format(model, area))
        f.write('{} {}\n'.format(model, area))

art = []
legend = plt.legend(models, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                    edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.4))
art.append(legend)
frame = legend.get_frame()
grid()
plt.savefig(save_path+'plots/roc_of_interest.png',additional_artists=art,
            bbox_inches="tight")
plt.close()
