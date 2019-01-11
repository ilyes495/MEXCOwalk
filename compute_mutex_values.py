from utils import *
from hist import plot_histogram_list
import sys
import math
import networkx as nx

path = '../' + network_name + '/out/edge_weights/pre/'
out_file = path + 'mutex.txt'
mod_out_file = path + 'mutex_ncomb.txt'
avgmod_out_file = path + 'mutex_nsep.txt'

def compute_mutex():
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id() # string to indices
    edge_list = load_edge_list()
    print data
    G = nx.Graph()

    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())

    N = len(genes)
    num_samples = len(load_patients_to_indices())  # number of patients

    weight_out_file = out_file
    fhout = open(weight_out_file, 'w+')

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_patients = []
        gene2_patients = []

        if gene1 in data:
            gene1_patients = data[gene1]
        if gene2 in data:
            gene2_patients = data[gene2]

        gene1_count = len(gene1_patients)
        gene2_count = len(gene2_patients)
        union = len(list(set(gene1_patients).union(gene2_patients)))

        if gene1_count + gene2_count != 0:
            mutex = float(union) / (gene1_count + gene2_count)
            print >>fhout, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(mutex)

    #plot_histogram_list(final_scores, 100)
    #fhoutdensity.close()
    fhout.close()

def compute_mutex_ncomb():
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id() # string to indices
    edge_list = load_edge_list()

    G = nx.Graph()

    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())

    N = len(genes)
    num_samples = len(load_patients_to_indices())  # number of patients

    weight_out_file = mod_out_file
    fhout = open(weight_out_file, 'w+')

    count_1s = 0
    count_0s = 0

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_neighbours = G[gene1_index]
        gene2_neighbours = G[gene2_index]

        gene_set = set(gene1_neighbours).union(set(gene2_neighbours))

        union_set = set()

        for gene in gene_set:
            if id_to_gene[gene] in data:
                union_set = union_set.union(data[id_to_gene[gene]])

        union = len(union_set)

        count = 0
        for gene in gene_set:
            if id_to_gene[gene] in data:
                count += len(data[id_to_gene[gene]])

        # print(str(union) + " " + str(count))
        if count != 0:
            mutex = float(union) / count

            if mutex == 1:
               count_1s += 1

            print >>fhout, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(mutex)
        else:
            count_0s += 1


    print (count_0s)
    #plot_histogram_list(final_scores, 100)
    #fhoutdensity.close()
    fhout.close()

def compute_mutex_nsep():
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id() # string to indices
    edge_list = load_edge_list()

    G = nx.Graph()

    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    G.remove_edges_from(G.selfloop_edges())

    N = len(genes)
    num_samples = len(load_patients_to_indices())  # number of patients

    weight_out_file = avgmod_out_file
    fhout = open(weight_out_file, 'w+')

    count_1s = 0
    count_0s = 0

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_neighbours = G[gene1_index]
        gene2_neighbours = G[gene2_index]

        union_set1 = set()
        union_set2 = set()
        count1 = 0
        count2 = 0

        if gene1 in data:
            union_set1 = set(data[gene1])
            count1 = len(data[gene1])
        if gene2 in data:
            union_set2 = set(data[gene2])
            count2 = len(data[gene2])

        for gene in gene1_neighbours:
            if id_to_gene[gene] in data and id_to_gene[gene] != gene2_index:  # exclude gene2 ?
                union_set1 = union_set1.union(data[id_to_gene[gene]])
                count1 += len(data[id_to_gene[gene]])
        for gene in gene2_neighbours:
            if id_to_gene[gene] in data and id_to_gene[gene] != gene1_index:
                union_set2 = union_set2.union(data[id_to_gene[gene]])
                count2 += len(data[id_to_gene[gene]])

        # print(str(union) + " " + str(count))

        union_set1 = len(union_set1)
        union_set2 = len(union_set2)

        m1 = 0
        m2 = 0
        if count1 != 0:
            m1 = float(union_set1) / count1
        if count2 != 0:
            m2 = float(union_set2) / count2

        mutex = (m1+m2)/2

        if mutex == 1:
           count_1s += 1

        print >>fhout, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(mutex)

    print (count_0s)
    #plot_histogram_list(final_scores, 100)
    #fhoutdensity.close()
    fhout.close()


compute_mutex()
#compute_mutex_ncomb()
#compute_mutex_nsep()
