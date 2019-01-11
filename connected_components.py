import networkx as nx
from utils import *
import numpy as np
import scipy.sparse as sp
import operator
import os


# Given the total number of genes N in the connected components, this script
# finds a threshold for hotnet2 and our algorithm that generates connected
# components with total number of genes N
models = [
    # "mutex", "cov", "mutex_cov", \
    # "mutex_wesme", "mutex_wesme_cov", \
    # "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
    # "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \
    #"cov_nsep",
      "cov_ncomb"

    # "mutex_t10_cov", \
    # "mutex_t05_ncomb_cov", "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov", "mutex_t05_nsep_cov_nsep", \
    # "mutex_t06_ncomb_cov", "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov", "mutex_t06_nsep_cov_nsep", \
    # "mutex_t07_ncomb_cov", "mutex_t07_ncomb_cov_ncomb", "mutex_t07_nsep_cov", "mutex_t07_nsep_cov_nsep", \
    # "mutex_t08_ncomb_cov", "mutex_t08_ncomb_cov_ncomb", "mutex_t08_nsep_cov", "mutex_t08_nsep_cov_nsep", \
    # "mutex_t09_ncomb_cov", "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov", "mutex_t1_ncomb_cov", \
    # "mutex_t1_ncomb_cov_ncomb", "mutex_t1_nsep_cov", "mutex_t1_nsep_cov_nsep", "mutex_t09_nsep_cov_nsep", \
    # "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", "mutex_t06_ncomb_cov_ncomb","mutex_t06_nsep_cov_nsep",\
    # "mutex_t07_ncomb_cov_ncomb","mutex_t07_nsep_cov_nsep", "mutex_t08_ncomb_cov_ncomb", "mutex_t08_nsep_cov_nsep",\
    # "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", "mutex_t10_ncomb_cov", "mutex_t10_nsep_cov",\
    # "mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep"\
    #"mutex_t09_nsep_cov_ncomb",\
    #"mutex_nsep_t06_a07_cov_nsep_d1", \
    # "mutex_nsep_t06_a07_cov_nsep_d2", "mutex_nsep_t06_a07_cov_nsep_d3",\
    # "mutex_nsep_t07_a07_cov_nsep_d1", "mutex_nsep_t07_a07_cov_nsep_d2", "mutex_nsep_t07_a07_cov_nsep_d3",\
    # "mutex_nsep_t08_a07_cov_nsep_d1", "mutex_nsep_t08_a07_cov_nsep_d2", "mutex_nsep_t08_a07_cov_nsep_d3"
     ]
for key in models:

    #key = "mutex_nsep_cov"
    threshold_start = 0.000214

    our_path = "../" + network_name + "/out/connected_components/" + key + "/"
    path = our_path
    print path

    if not os.path.exists(our_path):
        os.mkdir(our_path)

    id_to_gene = load_id_to_gene()
    LARGE_NUM = 100000
    k = 3
    our_E = sp.load_npz("../" + network_name + "/out/random_walk/"+ key+"_sparse_matrix_e.npz")
    E = our_E.toarray()

    # find threshold

    num_start = 2500    # max n
    num_end = 100       # min n

    # create the initial graph, ommit the edges smaller than the threshold
    N = len(id_to_gene)
    G = nx.DiGraph()

    for i in range(N):
        G.add_node(i)

    for i in range(N):
        for j in range(N):
            if E[i][j] > threshold_start:
                G.add_edge(j, i)


    # find the list of the connected components
    list_graphs = []
    num_nodes = 0
    count_comp = 0
    smallest_weight = LARGE_NUM
    smallest_edge_str = ''
    for comp in nx.strongly_connected_components(G):
        if len(comp) >= k:
            subG = G.subgraph(comp).copy()
            list_graphs.append(subG)
            num_nodes += len(comp)

            for e in subG.edges():
                u = e[0]
                v = e[1]
                key = str(count_comp) + '_' + str(u) + '_' + str(v)
                if E[v][u] < smallest_weight:
                    smallest_edge_str = key
                    smallest_weight = E[v][u]
            count_comp += 1
            print smallest_edge_str
    print 'number of nodes in the beginning is : ', num_nodes
    print smallest_edge_str


    while num_start >= num_end:
        print("a" + str(num_start))
        iteration = 0
        threshold_found = 0
        num_target_genes = num_start
        while num_nodes > num_target_genes:
             print(str(num_target_genes) + " " + str(num_nodes) + " " + str(iteration))

             #remove smallest edge
             print "smallest weight", smallest_weight

             smallest_edge = smallest_edge_str.split('_')
             threshold_found = smallest_weight
             graph_index = int(smallest_edge[0])

             subG = list_graphs[graph_index]
             num_nodes -= len(subG.nodes())

             subG.remove_edge(int(smallest_edge[1]), int(smallest_edge[2]))
             subG_comps = nx.strongly_connected_components(subG)

             del list_graphs[graph_index]
             largestCompSize = 20
             for comp in subG_comps:
                 if len(comp) >= k:
                     list_graphs.append(subG.subgraph(comp).copy())
                     num_nodes += len(comp)

             num_comps = len(list_graphs)
             print "numcomps", num_comps
             smallest_weight = LARGE_NUM
             smallest_edge_str = ''
             for i in range(num_comps):
                 graph = list_graphs[i]

                 if len(graph.nodes())>=largestCompSize:
                     gmlpath = path + "gml_n" + str(num_target_genes) + "largestComp.gml"
                     print "printing inside if largest cc", len(graph.nodes)
                     nx.write_gml(graph, gmlpath)
                 if len(graph.nodes()) >= k:
                     for e in graph.edges():
                         u = e[0]
                         v = e[1]
                         if E[v][u] < smallest_weight:
                             key = str(i) + '_' + str(u) + '_' + str(v)
                             smallest_edge_str = key
                             smallest_weight = E[v][u]
             iteration += 1


        outfilename = path + "cc_n" + str(num_target_genes) + "_k" + str(k) + "_d" + str(threshold_found) + ".txt"


        # with open(outfilename, "w") as f:
        #     for i in range(len(list_graphs)):
        #         for node_index in list_graphs[i]:
        #             f.write(id_to_gene[node_index] + " ")
        #         f.write("\n")
        #     f.close()
        num_start -= 100

