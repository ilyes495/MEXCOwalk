from utils import *
import sys
import math
import networkx as nx

def compute_edge_weights(key):
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id() # string to indices
    edge_list = load_edge_list()

    N = len(genes)
    num_samples = len(load_patients_to_indices())  # number of patients
    alpha = 0.7

    mutex_scores = {}
    prefix = '../' + network_name + '/out/edge_weights/pre/'
    with open(prefix + 'mutex.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            mutex_scores[line[0]+" "+line[1]] = float(line[2])

    mutex_ncomb_scores = {}
    with open(prefix + 'mutex_ncomb.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            mutex_ncomb_scores[line[0]+" "+line[1]] = float(line[2])

    mutex_nsep_scores = {}
    with open(prefix + 'mutex_nsep.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            mutex_nsep_scores[line[0]+" "+line[1]] = float(line[2])

    mutex_wesme_scores = {}
    with open(prefix + 'mutex_wesme.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            mutex_wesme_scores[line[0]+" "+line[1]] = float(line[2])

    cov_scores = {}
    with open(prefix + 'cov.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            cov_scores[line[0]+" "+line[1]] = float(line[2])

    cov_ncomb_scores = {}
    with open(prefix + 'cov_ncomb.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            cov_ncomb_scores[line[0]+" "+line[1]] = float(line[2])

    cov_nsep_scores = {}
    with open(prefix + 'cov_nsep.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            cov_nsep_scores[line[0]+" "+line[1]] = float(line[2])

    density1_scores = {}
    with open(prefix + 'density1.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            density1_scores[line[0]+" "+line[1]] = float(line[2])

    density2_scores = {}
    with open(prefix + 'density2.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            density2_scores[line[0]+" "+line[1]] = float(line[2])

    density3_scores = {}
    with open(prefix + 'density3.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            density3_scores[line[0]+" "+line[1]] = float(line[2])

    weight_out_file = '../' + network_name + '/out/edge_weights/' + key + '.txt'
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

        if gene1_count + gene2_count != 0:
            if key == 'mutex':
                mutex = mutex_scores[gene1 + " " + gene2]
                res = mutex
            elif key == 'cov':
                cov = cov_scores[gene1 + " " + gene2]
                res = cov
            elif key == 'mutex_cov':
                mutex = mutex_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex * cov
            elif key == 'mutex_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                res = mutex
            elif key == 'mutex_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                res = mutex
            elif key == 'mutex_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex * cov
            elif key == 'mutex_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex * cov
            elif key == 'cov_ncomb':
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                res = cov
            elif key == 'cov_nsep':
                cov = cov_nsep_scores[gene1 + " " + gene2]
                res = cov
            elif key == 'mutex_wesme':
                mutex = 1
                if gene1 + " " + gene2 not in mutex_wesme_scores:
                    if gene2 + " " + gene1 not in mutex_wesme_scores:
                        mutex = 1
                    else:
                        mutex = mutex_wesme_scores[gene2 + " " + gene1]
                else:
                    mutex = mutex_wesme_scores[gene1 + " " + gene2]
                res = mutex
            elif key == 'mutex_wesme_cov':
                mutex = 1
                if gene1 + " " + gene2 not in mutex_wesme_scores:
                    if gene2 + " " + gene1 not in mutex_wesme_scores:
                        mutex = 1
                    else:
                        mutex = mutex_wesme_scores[gene2 + " " + gene1]
                else:
                    mutex = mutex_wesme_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex * cov
            elif key == 'mutex_t10_cov':
                mutex = mutex_scores[gene1 + " " + gene2]
                if mutex < 1.0:
                    mutex = 0
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex * cov

#threshold 0.5

            elif key == 'mutex_t05_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex*cov
            elif key == 'mutex_t05_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                cov = cov_scores[gene1 + " " + gene2]
                res = mutex*cov
            elif key == 'mutex_t05_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t05_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t05_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t05_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.5:
                    mutex = 0
                res = mutex*cov

# threshold 0.6

            elif key == 'mutex_t06_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t06_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t06_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t06_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t06_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex * cov
            elif key == 'mutex_t06_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = mutex * cov

# threshold 0.7

            elif key == 'mutex_t07_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t07_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t07_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t07_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t07_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex * cov
            elif key == 'mutex_t07_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = mutex * cov

# threshold 0.8

            elif key == 'mutex_t08_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t08_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t08_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t08_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t08_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex * cov
            elif key == 'mutex_t08_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = mutex * cov

# threshold 0.9

            elif key == 'mutex_t09_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t09_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t09_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t09_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t09_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex * cov
            elif key == 'mutex_t09_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 0.9:
                    mutex = 0
                res = mutex * cov

# threshold 1

            elif key == 'mutex_t10_ncomb_cov':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t10_nsep_cov':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t10_ncomb_cov_ncomb':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t10_ncomb_cov_nsep':
                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t10_nsep_cov_nsep':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov
            elif key == 'mutex_t10_nsep_cov_ncomb':
                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                if mutex < 1:
                    mutex = 0
                res = mutex*cov

#Additional | Alpha set to 0.7
            elif key == 'mutex_nsep_t06_a07_cov_nsep_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t06_a07_cov_nsep_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t06_a07_cov_nsep_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t07_a07_cov_nsep_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t07_a07_cov_nsep_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t07_a07_cov_nsep_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_nsep_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_nsep_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_nsep_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            #phase2

            elif key == 'mutex_nsep_t06_a07_cov_ncomb_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t06_a07_cov_ncomb_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t06_a07_cov_ncomb_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_nsep_t07_a07_cov_ncomb_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t07_a07_cov_ncomb_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t07_a07_cov_ncomb_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_ncomb_d1':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_ncomb_d2':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_nsep_t08_a07_cov_ncomb_d3':

                mutex = mutex_nsep_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            #phase 3 mutex ncomb
            elif key == 'mutex_ncomb_t06_a07_cov_ncomb_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t06_a07_cov_ncomb_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t06_a07_cov_ncomb_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t07_a07_cov_ncomb_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t07_a07_cov_ncomb_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t07_a07_cov_ncomb_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_ncomb_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_ncomb_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_ncomb_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_ncomb_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            #mutex ncomb cov nsep
            elif key == 'mutex_ncomb_t06_a07_cov_nsep_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t06_a07_cov_nsep_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t06_a07_cov_nsep_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.6:
                    mutex = 0
                res = alpha*(mutex * cov) + (1-alpha)*density

            elif key == 'mutex_ncomb_t07_a07_cov_nsep_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t07_a07_cov_nsep_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t07_a07_cov_nsep_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.7:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_nsep_d1':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density1_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_nsep_d2':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density2_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            elif key == 'mutex_ncomb_t08_a07_cov_nsep_d3':

                mutex = mutex_ncomb_scores[gene1 + " " + gene2]
                cov = cov_nsep_scores[gene1 + " " + gene2]
                density = density3_scores[gene1 + " " + gene2]
                if mutex < 0.8:
                    mutex = 0
                res = alpha * (mutex * cov) + (1 - alpha) * density

            print >>fhout, id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(res)

    fhout.close()

#naming requires cautions
models = [
# "mutex", "cov", "mutex_cov", \
# "mutex_wesme", "mutex_wesme_cov", \
# "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
# "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \
#
# "mutex_t10_cov", \
# "mutex_t05_ncomb_cov", "mutex_t05_ncomb_cov_nsep", "mutex_t05_nsep_cov", "mutex_t05_nsep_cov_ncomb", \
# "mutex_t06_ncomb_cov", "mutex_t06_ncomb_cov_nsep", "mutex_t06_nsep_cov", "mutex_t06_nsep_cov_ncomb", \
# "mutex_t07_ncomb_cov", "mutex_t07_ncomb_cov_nsep", "mutex_t07_nsep_cov", "mutex_t07_nsep_cov_ncomb", \
# "mutex_t08_ncomb_cov", "mutex_t08_ncomb_cov_nsep", "mutex_t08_nsep_cov", "mutex_t08_nsep_cov_ncomb", \
# "mutex_t09_ncomb_cov", "mutex_t09_ncomb_cov_nsep", "mutex_t09_nsep_cov", "mutex_t09_nsep_cov_ncomb", \


# "mutex_t05_ncomb_cov_ncomb","mutex_t05_nsep_cov_nsep", "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov_nsep",\
# "mutex_t07_ncomb_cov_ncomb", "mutex_t07_nsep_cov_nsep", "mutex_t08_ncomb_cov_ncomb","mutex_t08_nsep_cov_nsep",\
# "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", "mutex_t10_ncomb_cov", "mutex_t10_nsep_cov",
# "mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep"\

"mutex_nsep_t06_a07_cov_nsep_d1", "mutex_nsep_t06_a07_cov_nsep_d2", "mutex_nsep_t06_a07_cov_nsep_d3",\
"mutex_nsep_t07_a07_cov_nsep_d1", "mutex_nsep_t07_a07_cov_nsep_d2", "mutex_nsep_t07_a07_cov_nsep_d3",\
"mutex_nsep_t08_a07_cov_nsep_d1", "mutex_nsep_t08_a07_cov_nsep_d2", "mutex_nsep_t08_a07_cov_nsep_d3",\

"mutex_nsep_t06_a07_cov_ncomb_d1", "mutex_nsep_t06_a07_cov_ncomb_d2", "mutex_nsep_t06_a07_cov_ncomb_d3",\
"mutex_nsep_t07_a07_cov_ncomb_d1", "mutex_nsep_t07_a07_cov_ncomb_d2", "mutex_nsep_t07_a07_cov_ncomb_d3",\
"mutex_nsep_t08_a07_cov_ncomb_d1", "mutex_nsep_t08_a07_cov_ncomb_d2", "mutex_nsep_t08_a07_cov_ncomb_d3",\

"mutex_ncomb_t06_a07_cov_ncomb_d1", "mutex_ncomb_t06_a07_cov_ncomb_d2", "mutex_ncomb_t06_a07_cov_ncomb_d3",\
"mutex_ncomb_t07_a07_cov_ncomb_d1", "mutex_ncomb_t07_a07_cov_ncomb_d2", "mutex_ncomb_t07_a07_cov_ncomb_d3",\
"mutex_ncomb_t08_a07_cov_ncomb_d1", "mutex_ncomb_t08_a07_cov_ncomb_d2", "mutex_ncomb_t08_a07_cov_ncomb_d3",\

"mutex_ncomb_t06_a07_cov_nsep_d1", "mutex_ncomb_t06_a07_cov_nsep_d2", "mutex_ncomb_t06_a07_cov_nsep_d3",\
"mutex_ncomb_t07_a07_cov_nsep_d1", "mutex_ncomb_t07_a07_cov_nsep_d2", "mutex_ncomb_t07_a07_cov_nsep_d3",\
"mutex_ncomb_t08_a07_cov_nsep_d1", "mutex_ncomb_t08_a07_cov_nsep_d2", "mutex_ncomb_t08_a07_cov_nsep_d3",\
]

for m in models:
    compute_edge_weights(m)
