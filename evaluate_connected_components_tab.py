from utils import *
import glob
import os
import math
import operator


def geo_mean_underflow(iterable):
    a = np.log(iterable)
    return np.exp(a.sum()/len(a))

def calculate_scores_for_subnetworks(filename, edge_list, outfile, num_genes):
    fhout = open(outfile, 'a+')
    num_samples = len(load_patients_to_indices())  # number of patients
    gene_to_id = load_gene_to_id()
    id_to_gene = load_id_to_gene()
    data = load_gene_vs_patient_data()

    scores = []
    with open(filename) as f:
        connected_components = f.readlines()
        num_components = float(len(connected_components))
        t_coverage = 0.0
        t_mutex = 0.0
        t_density = 0.0
        t_covmutex = 0.0


        for component_line in connected_components:
            component = component_line.split()
            coverage_list = []
            union_set = set([])
            individual_count = 0
            mult_count = 1
            for i in component:
                if i not in data:
                    print(i + " not found" + filename)
                    continue
                union_set = set(data[i]).union(union_set)
                individual_count = individual_count + len(data[i])
                mult_count *= (float(len(data[i]))/float(num_samples))
                coverage_list.append(len(data[i]) / float(num_samples))
            union_set = len(list(union_set))

            num_in = 0
            num_out = 0
            for edge in edge_list:
                a = id_to_gene[edge[0]] in component
                b = id_to_gene[edge[1]] in component
                if a and b:
                    num_in = num_in + 1
                elif a or b:
                    num_out = num_out + 1


            density = 0.0
            if num_in + num_out != 0:
                density = round(float(num_in) / (num_in + num_out),3)

            mutex = round(float(union_set) / individual_count,3)
            #coverage = float(geo_mean_underflow(coverage_list))#float(mult_count)**(1.0/len(component))
            coverage = round(float(union_set)/num_samples, 3)


            cov_mutex = mutex*coverage

            t_density = t_density + density
            t_coverage = t_coverage + coverage
            t_mutex = t_mutex + mutex
            t_covmutex = t_covmutex + cov_mutex

    avg_density = t_density / num_components
    avg_cov = t_coverage / num_components
    avg_mutex = t_mutex / num_components
    avg_covmutex = t_covmutex / num_components

    print >>fhout, str(num_genes) + '\t' + str(num_components) + \
    '\t' + str(t_density) + '\t' + str(avg_density) + \
    '\t' + str(t_coverage) + '\t' + str(avg_cov) + \
    '\t' + str(t_mutex) + '\t' + str(avg_mutex) + \
    '\t' + str(t_covmutex) + '\t' + str(avg_covmutex)

    return

def read_all_genes(filename):
    l = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            for i in line:
                l.append(i)
    return l

#****** COSMIC RELATED*******
def read_subnetworks(filename):
    fhmodule = open(filename)
    module_genes = []
    module_list = []
    huge_module = []
    for line in fhmodule:
       words = line.strip().split()
       module_list.append(words)
       for word in words:
           if word not in module_genes:
              module_genes.append(word)
       if 'TP53' in words:
           huge_module = words

    return (module_genes, module_list, huge_module)

def find_overlap(module, cosmic_genes):
    count = 0
    for gene in module:
       if gene in cosmic_genes:
           count += 1
    return (count, len(module))

def read_cosmic_genes():
    fhinput = open('../analysis_files/Census_allTue_May_23_12-08-15_2017.tsv')
    cosmic_genes = []
    line = fhinput.readline()
    for line in fhinput:
        cosmic_genes.append(line.split()[0])
    return cosmic_genes

def cosmic_overlap_analysis(our_file, hotnet2_file, cosmic_output_file):
    cosmic_genes = read_cosmic_genes()
    l = []
    hotnet_list = []
    our_list = []

    with open(cosmic_output_file, "a") as f:
        (our_module_genes, our_module_list, our_huge_module) = read_subnetworks(our_file)
        (hotnet_module_genes, hotnet_module_list, hotnet_huge_module) = read_subnetworks(hotnet2_file)
        #print(our_module_genes)
        overlap = len(list(set(our_module_genes).intersection(set(hotnet_module_genes))))
        overlap_huge = len(list(set(our_huge_module).intersection(set(hotnet_huge_module))))
        (count_overlap_ours, count) =  find_overlap(our_module_genes, cosmic_genes)
        (count_overlap_hotnet2, count) =  find_overlap(hotnet_module_genes, cosmic_genes)
        # number of genes, our overlap, hotnet overlap
        f.write(str(i*100 + 100) + "\t" + str(count_overlap_ours) + "\t" + str(count_overlap_hotnet2) + "\n")

    return [count_overlap_ours, count_overlap_hotnet2]

def prep_file_paths(key):
    paths_raw = glob.glob(key)

    dict_paths = {}
    for filename in paths_raw:
          #file_core = filename[filename.rindex('/')+1:]
          directory, file_core = os.path.split(filename)
          print file_core
          if file_core.startswith('cc'):
              num_genes = int(file_core.split('_')[1][1:])
              dict_paths[num_genes] = filename

    sorted_dict = sorted(dict_paths.items(), key=operator.itemgetter(0))

    paths = []
    for i in range(len(sorted_dict)):
        paths.append(sorted_dict[i][1])
    return paths



type_index = 2

models = [
    "hotnet2",\
    #"cluster_one",\
    "mutex", "cov", "mutex_cov", \
    "mutex_wesme", "mutex_wesme_cov", \
    "mutex_ncomb", "cov_nsep", "mutex_ncomb_cov", \
    "mutex_nsep", "cov_ncomb", "mutex_nsep_cov", \

    "mutex_t10_cov", \
    "mutex_t05_ncomb_cov", "mutex_t05_ncomb_cov_nsep", "mutex_t05_nsep_cov", "mutex_t05_nsep_cov_ncomb", \
    "mutex_t06_ncomb_cov", "mutex_t06_ncomb_cov_nsep", "mutex_t06_nsep_cov", "mutex_t06_nsep_cov_ncomb", \
    "mutex_t07_ncomb_cov", "mutex_t07_ncomb_cov_nsep", "mutex_t07_nsep_cov", "mutex_t07_nsep_cov_ncomb", \
    "mutex_t08_ncomb_cov", "mutex_t08_ncomb_cov_nsep", "mutex_t08_nsep_cov", "mutex_t08_nsep_cov_ncomb", \
    "mutex_t09_ncomb_cov", "mutex_t09_ncomb_cov_nsep", "mutex_t09_nsep_cov", "mutex_t09_nsep_cov_ncomb",\

    "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", "mutex_t06_ncomb_cov_ncomb","mutex_t06_nsep_cov_nsep",\
    "mutex_t07_ncomb_cov_ncomb","mutex_t07_nsep_cov_nsep", "mutex_t08_ncomb_cov_ncomb", "mutex_t08_nsep_cov_nsep",\
    "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", #"mutex_t10_ncomb_cov", "mutex_t10_nsep_cov",\
    # "mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep"\


    # "mutex_nsep_t06_a07_cov_ncomb_d1", "mutex_nsep_t06_a07_cov_ncomb_d2", "mutex_nsep_t06_a07_cov_ncomb_d3",\
    # "mutex_nsep_t07_a07_cov_ncomb_d1", "mutex_nsep_t07_a07_cov_ncomb_d2", "mutex_nsep_t07_a07_cov_ncomb_d3",\
    # "mutex_nsep_t08_a07_cov_ncomb_d1", "mutex_nsep_t08_a07_cov_ncomb_d2", "mutex_nsep_t08_a07_cov_ncomb_d3"
     ]
for key in models:
    our_path_key = key #'mutex_t05_ncomb_cov'

    hotnet_paths = prep_file_paths('../hint/out/connected_components/hotnet2/*.txt')
    our_paths = prep_file_paths('../hint/out/connected_components/' + our_path_key +'/*.txt')

    #os.system('mkdir ' + '../hint/out/evaluation_tab/'+ our_path_key)
    newpath = '../hint/out/evaluation_tab/'+ our_path_key
    if not os.path.exists(newpath):
        os.mkdir(newpath)

    hotnet_list = []
    our_list = []
    edge_list_original = load_edge_list()
    for i in range(min(25, len(our_paths))):
        num_genes = str(i*100 + 100)

        # COSMIC analysis
        cosmic_output_file = '../hint/out/evaluation_tab/' +  our_path_key + '/cosmic_' + our_path_key  + '.txt'
        cosmic_output_huge_file = '../hint/out/evaluation_tab/' + our_path_key  + '/cosmic_' + our_path_key  + '_huge.txt'
        [our_overlap, hotnet2_overlap] = cosmic_overlap_analysis(our_paths[i], hotnet_paths[i], cosmic_output_file)

        # GO TERM analysis
        go_output_file = '../hint/out/evaluation_tab/' +  our_path_key + '/GO_ne_' + our_path_key  + '.txt'
        os.system('python check_go_assoc_ne_tab.py ' + our_paths[i] + ' ' +  hotnet_paths[i] + ' ' + go_output_file )

        # Our evaluation
        outfile = '../hint/out/evaluation_tab/'+ our_path_key + '/our_eval_' + our_path_key + '.txt'
        outfilehotnet2 = '../hint/out/evaluation_tab/'+ our_path_key + '/our_eval_hotnet2_' + our_path_key + '.txt'
        calculate_scores_for_subnetworks(our_paths[i], edge_list_original, outfile, num_genes)
        calculate_scores_for_subnetworks(hotnet_paths[i], edge_list_original, outfilehotnet2, num_genes)
