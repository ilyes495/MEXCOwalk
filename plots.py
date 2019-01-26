from utils import *
import glob
import os
import math
import operator
from tqdm import trange, tqdm
from matplotlib import pyplot as plt
#from matplotlib import rcParams
from pylab import *

ts= [5, 6, 7, 8,9]
models = ["mutex_t0{}_nsep_cov".format(t) for t in ts]
colors = ['C0', 'C1', 'C4', 'C2', 'C3']

def diff_mutex_thresh_plot():
    ps = ['iwavg_cov','wavg_mutex','wavg_covmutex']
    labels = ['Coverage Score (CS)', 'Mutual Exclusion score ( MS )',  'Driver Module Set Score (DMSS)' ]
    for p,l in zip(ps,labels):
        Ns= []
        with open('../hint/out/evaluation_tab/{}.txt'.format(p)) as f:
            lines = f.readlines()
            models_ = lines[0].rstrip().split('\t')[1:]
            model2w = {s:[] for s in models_}
            for line in lines[1:]:
                line = line.rstrip().split('\t')
                Ns.append(int(line[0]))
                for m,w in zip(models_, line[1:]):
                    if float(w) == 0: continue
                    try:
                        model2w[m].append(float(w))
                    except:
                        print(m,w)

        axes(frameon=0)
        for m, c in zip(models, colors):
            if m == 'hier_hotnet2':
                plot([554,806], model2w[m], 'k*', markersize=12)
            elif m == 'memcover_v3':
                plot(Ns[:-9], model2w[m], '-o', color = c)
            else:
                try:
                    plot(Ns, model2w[m], '-o', color = c)
                except:
                    print(m, len(Ns), len(model2w[m]))

        art = []
        legend_ = ["MEXCOwalk_{}0{}".format(r'$\theta$',t) for t in ts]
        legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                            edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.3))
        art.append(legend)
        frame = legend.get_frame()
        plt.xlabel('total_genes')
        plt.ylabel(l)

        xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
        xticks(xtick, fontsize='x-small')
        # plt.ylim(0.79, 1)
        # axes(frameon=0)
        grid()
        plt.savefig('../results/roc_MEXCOwalk/plots/{}.pdf'.format(p), format= 'pdf',additional_artists=art,
                    bbox_inches="tight", dpi = 800)
        plt.close()


def Logpval_per_comp():
    models  = ["hotnet2","memcover_v1","memcover_v2","memcover_v3","mutex_t07_nsep_cov", 'hier_hotnet']
    labels = [
        "Hotnet2",
        "MEMCover_v1",
        "MEMCover_v2",
        "MEMCover_v3",
        "MEXCOwalk",
        "Hierarchical Hotnet",
        ]
    model2label = {m:l for m,l in zip(models,labels )}

    model2w = {}
    with open('../hint/out/cancer_subtype_test_results/Logpval_per_comp.txt') as f:
        lines = f.readlines()
        Ns = [int(n) for n in lines[0].rstrip().split('\t')]
        for line in lines[1:]:
            line = line.rstrip().split('\t')
            models_ = line[0]
            W = [float(w)*-1 for w in line[1:] if w != 'N/A']
            model2w[models_] = W

    axes(frameon=0)
    legend_ = []
    for m in models:
        legend_.append(model2label[m])
        if m == 'hier_hotnet':
            plot([554,806], model2w[m], 'k*', markersize=12)
        elif m == 'memcover_v3':
            plot(Ns[:-9], model2w[m], '-o')
        else:
            try:
                plot(Ns, model2w[m], '-o')
            except:
                print(m, len(Ns), len(model2w[m]))

    art = []


    legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                        edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.3))
    art.append(legend)
    frame = legend.get_frame()
    plt.xlabel('total_genes')
    plt.ylabel('Combined p-value (-log)')

    xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
    xticks(xtick, fontsize='x-small')
    plt.ylim(-0.01, 31)
    # axes(frameon=0)
    grid()
    plt.savefig('../hint/out/cancer_subtype_test_results/Logpval.pdf',format = 'pdf', additional_artists=art,
                bbox_inches="tight", dpi = 800)
    plt.close()


# diff_mutex_thresh_plot()


Logpval_per_comp()
