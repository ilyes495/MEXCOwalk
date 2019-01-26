import matplotlib.pyplot as plt
from  pylab import *
import numpy as np
import seaborn as sns

sns.set(color_codes=True)

source_folder = '../hint/out/edge_weights/pre/'
outfile = '../hint/figures/hist_data.txt'

models = [ "mutex_t07_nsep_cov", "mutex_t07_nsep_cov_nsep"]

hist_list = []
with open(outfile, 'w') as f_out:
    count = 0
    #for model in models:

    #infile = source_folder+model+ ".txt"
    infile = source_folder + "mutex_nsep.txt"
    list_edges = []

    with open(infile, 'r') as f_in:
        lines = f_in.readlines()

        for line in lines:
            line = line.strip().split()
            list_edges.append(float(line[2]))
            if float(line[2])>1:
                print (float(line[2]))

        f_in.close()
    # print "len", len(list_edges)
    # print "plotting..."

    # plt.figure(count)
    # axes(frameon=0)
    # plt.hist(list_edges, bins = 10, range =(0,1))
    sns.distplot(list_edges, bins = 10, kde = True,hist= True,  color ='r')
    plt.xlabel("MEX{} scores".format(r'$_n$'))
    plt.ylabel("Count")
    count += 1
    # sns.grid()
    plt.savefig('../hint/figures/histo.pdf', format = 'pdf',
                bbox_inches="tight", dpi = 800)
    plt.close()
