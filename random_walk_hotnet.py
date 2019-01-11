from utils import *
import scipy.sparse as sp
import numpy as np

# v2: fixed the F * H multiplication so that H only contains those genes that are expressed and in the network
# now, the size of H is 6930 * 6930, in v1 it was 9859 * 9859
def normalizeW(A):
    n = len(A)
    W = np.zeros((n, n))
    for j in range(n):
        d_j = float(A[:,j].sum())
        if(d_j == 0):    # added this check
            continue
        for i in range(n):
            W[i,j] = A[i,j]/d_j
    return W

def computeF(W, beta):
    n = len(W)
    return beta*np.linalg.inv(np.eye(n)-(1-beta)*W)

def computeH_():
    # some genes in patient_data are not in the graph, so we have to ignore those
    heat = load_gene_vs_patient_data_ready()
    gene_to_id = load_gene_to_id()  # gene symbols to integer indices
    gene_list = load_gene_list() # gene symbols
    # find the common genes between the heats and the network genes
    common_gene_list = set(heat.keys()).intersect(set(gene_list))
    N = len(gene_list)
    h = np.zeros((N, N))
    for key in heat:
        #  there are genes in heat data that do not exist in the network
        if key not in gene_to_id:
            continue
        gene_id = gene_to_id[key]
        h[gene_id][gene_id] = heat[key]
    return h, common_gene_list # size 9859 * 9859


def computeH():
    gene_to_id = load_gene_to_id()
    lines = []
    with open("../data/pan12gene2freq.txt") as f:
        lines = f.readlines()
    N = len(gene_to_id)
    h = np.zeros((N, N))
    yes = 0
    no = 0
    for line in lines:
        line = line.split()
        if line[0] not in gene_to_id.keys():
            no += 1
            continue
        yes += 1
        gene_id = gene_to_id[line[0]]
        h[gene_id][gene_id] = line[1]
    print(yes)
    print(no)
    return h

def computeW():
    gene_to_id = load_gene_to_id()
    edge_list = load_edge_list()

    N = len(gene_to_id)
    w = np.zeros((N, N))

    for edge in edge_list:
        w[edge[0]][edge[1]] = 1
        w[edge[1]][edge[0]] = 1
    return normalizeW(w)


path = "../" + network_name + "/out/random_walk/"

W = computeW()
sp_w = sp.csc_matrix(W)
sp.save_npz(path + "hotnet_sparse_matrix_w.npz", sp_w)
print("W")
print(W)

F = computeF(W, network_beta)
sp_f = sp.csc_matrix(F)
sp.save_npz(path + "hotnet_sparse_matrix_f.npz", sp_f)
print("F")
print(F)

H = computeH()
sp_h = sp.csc_matrix(H)
sp.save_npz(path + "hotnet_sparse_matrix_h.npz", sp_h)
print("H")
print(H)

E = np.dot(F, H)

sp_e = sp.csc_matrix(E)
sp.save_npz(path + "hotnet_sparse_matrix_e.npz", sp_e)
print("E")
print(E)
