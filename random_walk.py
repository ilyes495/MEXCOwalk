from utils import *
import scipy.sparse as sp
import numpy as np
import os

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

def computeH_pan():
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

def computeH():
    # some genes in patient_data are not in the graph, so we have to ignore those
    data = load_gene_vs_patient_data()
    N = len(load_patients_to_indices())  # number of patients
    gene_list = load_gene_list()
    genes = {}

    for key in gene_list:
        if key in data:
            genes[key] = len(data[key]) / float(N)

    N = len(gene_list) # number of genes
    h = np.zeros((N, N))
    gene_to_id = load_gene_to_id()
    no = 0
    for key in genes:
        # TODO there are genes in patiend data that does not exist in the graph
        if key not in gene_to_id:
            no += 1
            continue
        gene_id = gene_to_id[key]
        h[gene_id][gene_id] = genes[key]
    print(no)
    return h

def computeW():
    gene_to_id = load_gene_to_id()

    N = len(gene_to_id)
    w = np.zeros((N, N))

    with open(weight_out_file) as f:
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = line.split()
            w[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
    return normalizeW(w)

def computeW_combine(alpha1, alpha2):
    gene_to_id = load_gene_to_id()

    N = len(gene_to_id)
    w1 = np.zeros((N, N))
    w2 = np.zeros((N, N))
    with open(weight_out_file) as f:
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = line.split()
            #print line
            w1[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w1[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
            w2[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[3]
            w2[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[3]


    w1_norm = normalizeW(w1)
    w2_norm = normalizeW(w1)
    return alpha1 * w1_norm + alpha2 * w2_norm

if __name__ == "__main__":
    network_file = 'hint'
    keys = [
    # "mutex", "cov", "mutex_cov", \
    # "mutex_wesme", "mutex_wesme_cov", \
    # "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
    # "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \
    #
    # "mutex_t10_cov", \
    # "mutex_t05_ncomb_cov", "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov", "mutex_t05_nsep_cov_nsep", \
    # "mutex_t06_ncomb_cov", "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov", "mutex_t06_nsep_cov_nsep", \
    # "mutex_t07_ncomb_cov", "mutex_t07_ncomb_cov_ncomb", "mutex_t07_nsep_cov", "mutex_t07_nsep_cov_nsep", \
    # "mutex_t08_ncomb_cov", "mutex_t08_ncomb_cov_ncomb", "mutex_t08_nsep_cov", "mutex_t08_nsep_cov_nsep", \
    # "mutex_t09_ncomb_cov", "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov", "mutex_t09_nsep_cov_nsep", \
    # "mutex_t1_ncomb_cov", "mutex_t1_ncomb_cov_ncomb", "mutex_t1_nsep_cov", "mutex_t1_nsep_cov_nsep"\
    # "mutex_nsep_t06_a07_cov_nsep_d1", "mutex_nsep_t06_a07_cov_nsep_d2", "mutex_nsep_t06_a07_cov_nsep_d3", \
    # "mutex_nsep_t07_a07_cov_nsep_d1", "mutex_nsep_t07_a07_cov_nsep_d2", "mutex_nsep_t07_a07_cov_nsep_d3", \
    # "mutex_nsep_t08_a07_cov_nsep_d1", "mutex_nsep_t08_a07_cov_nsep_d2", "mutex_nsep_t08_a07_cov_nsep_d3"
    # "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov_nsep", \
    # "mutex_t07_ncomb_cov_ncomb", "mutex_t07_nsep_cov_nsep", "mutex_t08_ncomb_cov_ncomb", "mutex_t08_nsep_cov_nsep", \
    # "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", "mutex_t10_ncomb_cov", "mutex_t10_nsep_cov",\
    # "mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep" \
    #     "mutex_nsep_t06_a07_cov_nsep_d1", "mutex_nsep_t06_a07_cov_nsep_d2", "mutex_nsep_t06_a07_cov_nsep_d3", \
    #     "mutex_nsep_t07_a07_cov_nsep_d1", "mutex_nsep_t07_a07_cov_nsep_d2", "mutex_nsep_t07_a07_cov_nsep_d3", \
    #     "mutex_nsep_t08_a07_cov_nsep_d1", "mutex_nsep_t08_a07_cov_nsep_d2", "mutex_nsep_t08_a07_cov_nsep_d3", \
    #
    #     "mutex_nsep_t06_a07_cov_ncomb_d1", "mutex_nsep_t06_a07_cov_ncomb_d2", "mutex_nsep_t06_a07_cov_ncomb_d3", \
    #     "mutex_nsep_t07_a07_cov_ncomb_d1", "mutex_nsep_t07_a07_cov_ncomb_d2", "mutex_nsep_t07_a07_cov_ncomb_d3", \
    #     "mutex_nsep_t08_a07_cov_ncomb_d1", "mutex_nsep_t08_a07_cov_ncomb_d2", "mutex_nsep_t08_a07_cov_ncomb_d3", \

        "mutex_ncomb_t06_a07_cov_ncomb_d1", "mutex_ncomb_t06_a07_cov_ncomb_d2", "mutex_ncomb_t06_a07_cov_ncomb_d3", \
        "mutex_ncomb_t07_a07_cov_ncomb_d1", "mutex_ncomb_t07_a07_cov_ncomb_d2", "mutex_ncomb_t07_a07_cov_ncomb_d3", \
        "mutex_ncomb_t08_a07_cov_ncomb_d1", "mutex_ncomb_t08_a07_cov_ncomb_d2", "mutex_ncomb_t08_a07_cov_ncomb_d3", \

        "mutex_ncomb_t06_a07_cov_nsep_d1", "mutex_ncomb_t06_a07_cov_nsep_d2", "mutex_ncomb_t06_a07_cov_nsep_d3", \
        "mutex_ncomb_t07_a07_cov_nsep_d1", "mutex_ncomb_t07_a07_cov_nsep_d2", "mutex_ncomb_t07_a07_cov_nsep_d3", \
        "mutex_ncomb_t08_a07_cov_nsep_d1", "mutex_ncomb_t08_a07_cov_nsep_d2", "mutex_ncomb_t08_a07_cov_nsep_d3", \
        ]

    for key in keys:

        path = "../" + network_file + "/out/random_walk/"
        weight_out_file = "../" + network_file + "/out/edge_weights/" + key + '.txt'

        W = computeW()
        sp_w = sp.csc_matrix(W)
        sp.save_npz(path  + key + "_sparse_matrix_w.npz", sp_w)
        print("W")
        print(W)

        F = computeF(W, network_beta)
        sp_f = sp.csc_matrix(F)
        sp.save_npz(path + key + "_sparse_matrix_f.npz", sp_f)
        print("F")
        print(F)

        H = computeH_pan()
        sp_h = sp.csc_matrix(H)
        sp.save_npz(path  + key + "_sparse_matrix_h.npz", sp_h)
        print("H")
        print(H)

        E = np.dot(F, H)
        sp_e = sp.csc_matrix(E)
        sp.save_npz(path  + key +"_sparse_matrix_e.npz", sp_e)
        print("E")
        print(E)
