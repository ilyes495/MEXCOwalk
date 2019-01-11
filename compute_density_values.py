from utils import *
import sys
import math
import networkx as nx

# look at the subgraph induced by u, v and their neighbors, calculate the number of existing edges, divide by # of possible edges
def calc_edge_density_v1(G, index1, index2):
    index1_neighbors = G[index1].keys()
    index2_neighbors = G[index2].keys()
    all_nodes = list(set(index1_neighbors).union( set(index2_neighbors) ))
    all_nodes.append(index1)
    all_nodes.append(index2)
    num_nodes = len(all_nodes)
    Gsub = G.subgraph(all_nodes)
    density = float(Gsub.number_of_edges()) / (num_nodes * (num_nodes-1) / 2)
    return density

# look at the subgraph induced by u, v and their neighbors, calculate the number of existing edges, divide by # of possible edges
def calc_edge_density_v2(G, index1, index2):
    index1_neighbors = G[index1].keys()
    index2_neighbors = G[index2].keys()
    core_nodes = list(set(index1_neighbors).union( set(index2_neighbors) ))
    core_nodes.append(index1)
    core_nodes.append(index2)
    outer_layer = []
    for node in core_nodes:
        neighbors = G[node].keys()
        for n in neighbors:
            if n not in core_nodes:
                outer_layer.append(n)

    all_nodes = list(core_nodes)
    all_nodes.extend(outer_layer)
    Gsub = G.subgraph(all_nodes)
    w_in = 0
    w_out = 0
    for edge in Gsub.edges():
        node1_core = edge[0] in core_nodes
        node2_core = edge[1] in core_nodes
        if node1_core and node2_core:
            w_in += 1
        elif node1_core or node2_core:
            w_out += 1
    return float(w_in) / (w_in + w_out)

def calc_gene_density(G, genes, edge_list):

    dict_degree = {}
    dict_density = {}
    for gene in genes:
        index = gene_to_id[gene]
        neighbors = G[index].keys() #[n for n in G.neighbors_iter(index)]
        dict_degree[index] = len(neighbors)
        neighbors.append(index)
        num_neighbors = len(neighbors)
        Gsub = G.subgraph(neighbors)
        density = float(Gsub.number_of_edges()) / (num_neighbors * (num_neighbors-1) / 2)
        dict_density[index] = density
    return (dict_degree, dict_density)


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

d1_file = open('../' + network_name + '/out/edge_weights/pre/density1.txt', 'w')
d2_file = open('../' + network_name + '/out/edge_weights/pre/density2.txt', 'w')
d3_file = open('../' + network_name + '/out/edge_weights/pre/density3.txt', 'w')

# gene density
(dict_gene_degree, dict_gene_density) = calc_gene_density(G, genes, edge_list)

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

    if gene1_count + gene2_count != 0:
        density_v1 = calc_edge_density_v1(G, gene1_index, gene2_index)
        density_v2 = calc_edge_density_v2(G, gene1_index, gene2_index)
        density_v3 = (dict_gene_density[gene1_index] + dict_gene_density[gene2_index]) / 2
        print >>d1_file, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(density_v1)
        print >>d2_file, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(density_v2)
        print >>d3_file, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(density_v3)

d1_file.close()
d2_file.close()
d3_file.close()
