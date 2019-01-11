import networkx as nx
from utils import *
import numpy as np
import scipy.sparse as sp
import operator
import os
from collections import deque
import glob
import tqdm
#change in 08.01.2019 corrected the large split part.


# Given the total number of genes N in the connected components, this script
# finds a threshold for hotnet2 and our algorithm that generates connected
# components with total number of genes N
models = [
    #"cov"
    # "hotnet2",
    "mutex", "mutex_cov", #"cov",\
    "mutex_wesme", "mutex_wesme_cov", \
    "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
    "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \



    #these have threshold >0.00007
    "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", \
    "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov_nsep", \
    "mutex_t07_ncomb_cov_ncomb",  "mutex_t07_nsep_cov_nsep", \
    "mutex_t08_ncomb_cov_ncomb",  "mutex_t08_nsep_cov_nsep", \
    "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", \

    #these have common threshold > 0.0002
    "mutex_t05_ncomb_cov", "mutex_t05_nsep_cov", "mutex_t06_ncomb_cov",
    "mutex_t06_nsep_cov", \
    "mutex_t07_ncomb_cov", "mutex_t07_nsep_cov", "mutex_t08_ncomb_cov", "mutex_t08_nsep_cov", \
    "mutex_t09_ncomb_cov", "mutex_t09_nsep_cov",\
    "mutex_t10_cov", \
    "mutex_t05_ncomb_cov_nsep", "mutex_t05_nsep_cov_ncomb", "mutex_t06_ncomb_cov_nsep", "mutex_t06_nsep_cov_ncomb", \
    "mutex_t07_ncomb_cov_nsep", "mutex_t07_nsep_cov_ncomb", "mutex_t08_ncomb_cov_nsep", "mutex_t08_nsep_cov_ncomb", \
    "mutex_t09_ncomb_cov_nsep", "mutex_t09_nsep_cov_ncomb"  \

     ]
id_to_gene = load_id_to_gene()
gene_to_id = load_gene_to_id()

def find_largest_comp(list_comp_graph):
    large = list_comp_graph[0]
    large_index = 0
    for i in range(len(list_comp_graph)):
        if len(large) < len(list_comp_graph[i]):
            large = list_comp_graph[i]
            large_index = i
    print "large", large
    return (large,large_index)

def list_graph_to_list_comp(list_graph):
    list_comp = []
    for i in range (len(list_graph)):
        list_comp.extend([list_graph[i].nodes])

    return list_comp


def max_comp_length(list_graph):
    max_length = len(list_graph[0])

    for i in range (len(list_graph)):
        if len(list_graph[i])> max_length:
            max_length = len(list_graph[i])

    return max_length

def find_max_outdegree(comp):
    max_out_degree = 0
    largest_node_id = 0
    largest_gene_id = ''

    for n in comp.nodes():
        out_n = comp.out_degree(n)
        if max_out_degree < out_n:
            max_out_degree = out_n
            largest_node_id = n
            largest_gene_id = id_to_gene[n]
    print largest_gene_id

    return (max_out_degree, largest_node_id, largest_gene_id)

#star construction by finding all interconnected genes
def star_construction(large_comp, largest_node_id):
    all_neighbor_of_largest = set()
    star_nodes_set = set([largest_node_id])

    for e in large_comp.edges:

        if largest_node_id in e:

            all_neighbor_of_largest.update(e)

    print "neighbors", all_neighbor_of_largest

    all_neighbor_of_largest.remove(largest_node_id)

    for u in all_neighbor_of_largest:
        in_comp = True
        for e in large_comp.edges:
            if e[0] == u:
                v = e[1]

                if v != largest_node_id and (v not in all_neighbor_of_largest):

                    in_comp = False
                    print "false", e
                    break
            elif e[1] == u:
                v = e[0]

                if v != largest_node_id and (v not in all_neighbor_of_largest):
                    in_comp = False
                    print "false", e
                    break

        if in_comp == True:
            star_nodes_set.update([u])

    print "\n\n", star_nodes_set, "\n\n"

    return list(star_nodes_set)


#passes a list of component graph and returns connected commponents graph?


def split_large_components(list_comp_graph, large_threshold = 0):


    list_comp = list_graph_to_list_comp(list_comp_graph)
    set_all_genes_before = set([n for comp in list_comp for n in comp])

    threshold_small_comp = 3


    gene_set = set()
    for i in range(len(list_comp)):
        gene_set = gene_set.union(set(list_comp[i]))
    #print "gene_set", gene_set

    max_out_degree_all = 0
    #finding large graph threshold, take the max of max-out-degree of each component
    for comp in list_comp_graph:
         max_out_degree_pre, largest_node_id_pre, largest_gene_id_pre = find_max_outdegree(comp)
         print "maxoutdegree", max_out_degree_pre
         if max_out_degree_all < max_out_degree_pre:
               max_out_degree_all = max_out_degree_pre


    if large_threshold == 0:
        threshold_large_comp = max_out_degree_all
    else:
        threshold_large_comp = large_threshold


    list_graph_leftover = []
    list_large_components = []
    for comp in list_comp_graph:
        if len(comp)>= threshold_large_comp:
            list_large_components.append(comp)
            print "\nanother large component added\n"
        else:
            list_graph_leftover.append(comp)

    #going through all large components
    all_modified_component_list = []
    print "\nlenght of large comps\n",len(list_large_components)
    for lc in list_large_components:
        main_comp_list = []
        small_comp_list = []
        #large_graph_queue
        large_comp_queue = deque()
        large_comp_queue.append(lc)
        while len(large_comp_queue) >0:
            print "\nnew element in large queue\n"
            largest_comp = large_comp_queue.popleft()
            print "len", len(largest_comp)
            max_out_degree, largest_node_id, largest_gene_id = find_max_outdegree(largest_comp)

            reduced_comps = largest_comp.copy()

            removable_nodes_list = star_construction(largest_comp,largest_node_id)
            print "\nremovable nodes list", removable_nodes_list

            print "\nnodes before removal", len(reduced_comps)
            print "edges before removal", len(reduced_comps.edges)

            #reduced comps -> largest graph - large gene and neighbors
            reduced_comps.remove_nodes_from(removable_nodes_list)

            print "\nnodes after removal", len(reduced_comps)
            print "edges after removal", len(reduced_comps.edges)

            #adding to LOM or SCL
            #checking if star construction is less than k
            largest_gene_graph = largest_comp.subgraph(removable_nodes_list).copy()
            if len(largest_gene_graph) < threshold_small_comp:
                small_comp_list.append(largest_gene_graph)
            else:
                main_comp_list.append(largest_gene_graph)

            scc = nx.strongly_connected_components(reduced_comps)

            for comp in scc:
                print "comp", comp
                if len(comp)> threshold_large_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    large_comp_queue.append(g)
                    print "\nlarge_comp_queue"
                elif len(comp)< threshold_small_comp:
                    g = reduced_comps.subgraph(comp).copy()
                    small_comp_list.append(g)
                    print "\nsmall"
                else:
                    g = reduced_comps.subgraph(comp).copy()
                    main_comp_list.append(g)
                    print "\nLOM"

        print "len list main b4", len(main_comp_list )

        temp_list = []
        for scm in range(len(main_comp_list)):

            if len(main_comp_list[scm]) < threshold_small_comp:
                small_comp_list.append(main_comp_list[scm])

            else:
                temp_list.append(main_comp_list[scm])

        main_comp_list = temp_list[:]
        print "len list main aftr", len(main_comp_list)

        for scomp in range(len(small_comp_list)):
            max_comp_score = 0
            max_comp_index = 0
            for comp_index in  range(len(main_comp_list)):
                comp_score = 0
                for gene_m in main_comp_list[comp_index]:

                    for gene_s in small_comp_list[scomp]:

                        if (gene_s, gene_m) in lc.edges:
                            comp_score +=1
                        if (gene_m, gene_s) in lc.edges:
                            comp_score += 1

                if comp_score > max_comp_score:
                    max_comp_score = comp_score
                    max_comp_index = comp_index

            print "\ncomponent", small_comp_list[scomp].nodes, "maxcomp score+index", max_comp_score, max_comp_index
            print "list of comps nodes before", len (main_comp_list[max_comp_index]),main_comp_list[max_comp_index].nodes
            main_comp_list[max_comp_index].add_nodes_from(small_comp_list[scomp])

            temp_subgraph_nodes = main_comp_list[max_comp_index].nodes
            print "tempsubg", temp_subgraph_nodes
            # also add the edges
            main_comp_list[max_comp_index] = lc.subgraph(temp_subgraph_nodes).copy()
            print "list of comps nodes after", len(main_comp_list[max_comp_index]), main_comp_list[max_comp_index].nodes

        all_modified_component_list.extend(main_comp_list[:])

    for x in all_modified_component_list:
         print x.nodes

    set_all_genes_after = set([n for comp in list_graph_leftover[:] + all_modified_component_list[:] for n in comp.nodes])

    print('Set before: ', set_all_genes_before )
    print('Set after: ', set_all_genes_after )

    assert set_all_genes_after == set_all_genes_before
    return list_graph_leftover[:] + all_modified_component_list[:]



for key in tqdm.tqdm(models):

    #key = "mutex_nsep_cov"
    #threshold_start = 0.00007  # for 1000ncomb
    #threshold_start = 0.00009 #for 1000ncomb
    #threshold_start = 0.00029 #for 100ncomb

    #input file path
    filepath = "../" + network_name + "/out/connected_components_original/" + key + "/"
    our_path = "../" + network_name + "/out/connected_components_isolarge/" + key + "/"

    try:
        # Get the starting threshold from the file name at n =1000
        threshold_start = float(glob.glob(filepath+'cc_n1000_*')[0].split('/')[-1].split('d')[1][:-4])
    except:
        print(glob.glob(filepath+'cc_n1000_*'))

    path = our_path
    print path

    if not os.path.exists(our_path):
        os.mkdir(our_path)

    LARGE_NUM = 100000
    k = 3
    our_E = sp.load_npz("/Volumes/My Passport/random_walk/"+ key+"_sparse_matrix_e.npz")
    E = our_E.toarray()
    # find threshold
    num_start = 371  # max n
    num_end = 371  # min n

    # create the initial graph, omit the edges smaller than the threshold
    print('DiGraph...')
    N = len(id_to_gene)
    G = nx.DiGraph()
    print('DiGraph...')
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

            # remove smallest edge

            smallest_edge = smallest_edge_str.split('_')
            threshold_found = smallest_weight
            graph_index = int(smallest_edge[0])

            subG = list_graphs[graph_index]
            num_nodes -= len(subG.nodes())


            subG.remove_edge(int(smallest_edge[1]), int(smallest_edge[2]))
            subG_comps = nx.strongly_connected_components(subG)

            del list_graphs[graph_index]
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

                if len(graph.nodes()) >= k:
                    for e in graph.edges():
                        u = e[0]
                        v = e[1]
                        if E[v][u] < smallest_weight:
                            key = str(i) + '_' + str(u) + '_' + str(v)
                            smallest_edge_str = key
                            smallest_weight = E[v][u]
            iteration += 1

        ################################################################################################################
        #finding largest component
        #list_graphs contain all comps
        #largest_comp contains the largest comp

        outfilename = filepath + "cc_n" + str(num_target_genes) + "_d" + str(threshold_found) + ".txt"


        #prints the top [largest component size - difference] genes (or target genes)
        with open(outfilename, "w") as f:
            for i in range(len(list_graphs)):
                f.write(" ".join([id_to_gene[idx] for idx in list_graphs[i]])+"\n")


        largest_comp, largest_comp_index = find_largest_comp(list_graphs)
        component_list = list_graphs[:]
        #print len(largest_comp)
        largeset = set([id_to_gene[idx] for idx in largest_comp])
        print len(largest_comp), " largeset", largeset

        largest_gene_comp_list = []
        count = 0

        ##list_comp_graphs, largest_gene_graph = isolate_large(component_list)
        # 1: identifying the largest component graph
        # 2: identifying the gene with most connections
        # 3: removing that gene and its neighbors (largest_gene_graph)from the main graph (input graph)
        # 4: Run SCC on the remaining graph and list the connected components (list_comp_graphs)
        # 5: the function was run on only the largest gene, but other genes remain from the connected components
        # calculated earlier, save those components in (list_graphs_leftover).
        # 6: The first instance we only work with the largest component


        #inputs = component list, threshold (0 means largest component out degree)
        final_comp_list_graphs = split_large_components(component_list)

        # write gml from deleted genes
        # gmlpath = path + "del_genes_gml_n" + str(num_target_genes) + "largestComp.gml"
        # nx.write_gml(largest_gene_graph, gmlpath)



        outfilename = path + "cc_n" + str(num_target_genes) + "_"+ str(len(largest_comp))+ "_"+\
                      str(len(final_comp_list_graphs)) + "_d" + str(threshold_found) + ".txt"


        #prints the top [largest component size - difference] genes (or target genes)
        with open(outfilename, "w") as f:
            for i in range(len(final_comp_list_graphs)):
                for node_index in final_comp_list_graphs[i]:
                    f.write(id_to_gene[node_index] + " ")
                f.write("\n")
            f.close()



        num_start -= 100
