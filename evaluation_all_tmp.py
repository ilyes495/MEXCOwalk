from utils import *
import glob
import os
import math
import operator
from tqdm import trange, tqdm

def calculate_scores_for_subnetworks(filename, edge_list, outfile,outfile_tab, num_genes):

    fhout = open(outfile, 'w')
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

        print >>fhout, 'num_comp\tcoverage\tmutex\tcov_mutex\tdensity'
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
            union_set = len(list(union_set))#?

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



            #coverage = float(geo_mean_underflow(coverage_list))#**(1.0/len(component))#round(float(union_set) / num_samples,3)
            try:
                coverage = round(float(union_set) / num_samples,3)
                mutex = round(float(union_set) / individual_count,3) if individual_count != 0 else 0
            except:
                print('file: ',outfile)
            cov_mutex = mutex*coverage

            t_density = t_density + density
            t_coverage = t_coverage + coverage
            t_mutex = t_mutex + mutex
            t_covmutex = t_covmutex + cov_mutex

            print >>fhout, str(len(component)) + '\t' + str(coverage) +\
            '\t' + str(mutex) + '\t' +  str(cov_mutex) + '\t' +\
            str(density)
            # scores.append(a[0]*coverage + a[1]*mutex + a[2]*density)
            scores.append([len(component), coverage, mutex, cov_mutex, density])
    avg_density = t_density / num_components
    avg_cov = t_coverage / num_components
    avg_mutex = t_mutex / num_components
    avg_covmutex = t_covmutex / num_components

    print >>fhout, '****AVERAGE ACROSS ALL MODULES*****'
    print >>fhout, str(num_genes) + '\t' + str(num_components) + \
    '\t' + str(t_density) + '\t' + str(avg_density) + \
    '\t' + str(t_coverage) + '\t' + str(avg_cov) + \
    '\t' + str(t_mutex) + '\t' + str(avg_mutex) + \
    '\t' + str(t_covmutex) + '\t' + str(avg_covmutex)
    fhout.close()

    fhout = open(outfile_tab, 'a+')
    print >>fhout, str(num_genes) + '\t' + str(num_components) + \
    '\t' + str(t_density) + '\t' + str(avg_density) + \
    '\t' + str(t_coverage) + '\t' + str(avg_cov) + \
    '\t' + str(t_mutex) + '\t' + str(avg_mutex) + \
    '\t' + str(t_covmutex) + '\t' + str(avg_covmutex)
    return [t_density, avg_density, t_coverage, avg_cov, t_mutex, avg_mutex, t_covmutex, avg_covmutex]   #[t_coverage, t_mutex, t_density, t_density+t_coverage+t_mutex]

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

def cosmic_overlap_analysis(our_file, hotnet2_file,cosmic_genes, module_name, num_genes):
    #print our_file

    cosmic_file = '../hint/out/tmp/evaluation/{}.txt'.format(module_name)
    cosmic_huge_file = '../hint/out/tmp/evaluation/{}_huge.txt'.format(module_name)
    cosmic_tab_file = '../hint/out/tmp/evaluation_tab/{}.txt'.format(module_name)

    #cosmic_genes = read_cosmic_genes()
    l = []
    hotnet_list = []
    our_list = []

    our_largest_module = ''
    hotnet_largest_module = ''
    our_largest_module_cosmic = ''
    hotnet_largest_module_cosmic = ''
    x_axes = ''
    (our_module_genes, our_module_list, our_huge_module) = read_subnetworks(our_file)
    (hotnet_module_genes, hotnet_module_list, hotnet_huge_module) = read_subnetworks(hotnet2_file)
    #print(our_module_genes)
    overlap = len(list(set(our_module_genes).intersection(set(hotnet_module_genes))))
    overlap_huge = len(list(set(our_huge_module).intersection(set(hotnet_huge_module))))
    (count_overlap_ours, count) =  find_overlap(our_module_genes, cosmic_genes)
    (count_overlap_hotnet2, count) =  find_overlap(hotnet_module_genes, cosmic_genes)

    with open(cosmic_file, "a") as f:
        f.write("***Top " + str(num_genes) + " genes***\n")
        f.write("Ours vs hotnet2 subnetw. total overlap:  " + str(overlap) + "\n")
        f.write("Huge module sizes, ours: " +  str(len(our_huge_module))+   " Hotnet2: " +  str(len(hotnet_huge_module)) + " overlap: " + str(overlap_huge) + "\n")
        f.write("Overlap with cosmic, ours: " + str(count_overlap_ours) + " hotnet2:  " + str(count_overlap_hotnet2) +  " diff: " + str(count_overlap_ours - count_overlap_hotnet2) +"\n")
        count_overlap_huge_ours = find_overlap(our_huge_module, cosmic_genes)[0]
        count_overlap_huge_hotnet = find_overlap(hotnet_huge_module, cosmic_genes)[0]
        f.write("Overlap with cosmic only huge module, ours: "  + str(count_overlap_huge_ours) + " hotnet2: " + str(count_overlap_huge_hotnet) +  " diff:" + str(count_overlap_huge_ours - count_overlap_huge_hotnet) +  "\n\n" )
        x_axes += str(num_genes) + '\t'
        our_largest_module += str(len(our_huge_module)) + '\t'
        hotnet_largest_module += str(len(hotnet_huge_module)) + '\t'
        our_largest_module_cosmic += str(count_overlap_huge_ours) + '\t'
        hotnet_largest_module_cosmic += str(count_overlap_huge_hotnet) + '\t'

    with open(cosmic_huge_file, "w+") as f:
      f.write(x_axes + '\n')
      f.write(our_largest_module + '\n')
      f.write(hotnet_largest_module + '\n')
      f.write(our_largest_module_cosmic + '\n')
      f.write(hotnet_largest_module_cosmic + '\n'	)

    with open(cosmic_tab_file, "a") as f:
      # number of genes, our overlap, hotnet overlap
      f.write(str(num_genes) + "\t" + str(count_overlap_ours) + "\t" + str(count_overlap_hotnet2) + "\n")

    return [count_overlap_ours, count_overlap_hotnet2]

def generate_cosmic_analysis_file(path_pre = "../hint/out/tmp/evaluation_tab/"):
    d = {}
    for key in tqdm(models, desc='running cosmic analysis'):
        #print key + "\n"
        if key == "hotnet2":
            continue
        with open(path_pre + key + "/cosmic_" + key + ".txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                N = line[0].strip()
                our = line[1].strip()
                d[str(N) + key] = our

    key = "hotnet2"
    with open(path_pre + key + "/cosmic_" + key + ".txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            N = line[0].strip()
            hotnet = line[2].strip()
            d[str(N) + key] = hotnet

    with open(path_pre + "cosmic" + ".txt", "w") as ff:
        ff.write('N\t')
        for k in models:
            ff.write(k + '\t')
        ff.write('\n')
        num_genes_list = range(100,1100,100)#+[371]+range(400,1100,100)
        for i in num_genes_list:
            ff.write(str(i))
            for key in models:
                if str(i) + key in d:
                    ff.write('\t' + d[str(i) + key])
                else:
                    ff.write('\t' + ' ')
            ff.write('\n')




def generate_our_eval_files(path_pre = "../hint/out/tmp/evaluation_tab/"):
    eval_list = ['num_components', \
    't_density', 'avg_density', \
    't_coverage', 'avg_cov', \
    't_mutex', 'avg_mutex', \
    't_covmutex', 'avg_covmutex']

    for i in trange(len(eval_list), desc='running our evaluation'):
        eval = eval_list[i]
        d = {}
        for key in models:
            #print(key)
            # if key == "hotnet2":
            #     continue
            with open(path_pre + key + "/our_eval_" + key + ".txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.split()
                    N = line[0].strip()
                    e = line[i+1].strip()
                    d[str(N) + key] = e

        # key = "hotnet2"
        # with open(path_pre + 'hotnet2' + "/our_eval_hotnet2.txt", "r") as f:
        #     lines = f.readlines()
        #     for line in lines:
        #         line = line.split()
        #         N = line[0].strip()
        #         e = line[i+1].strip()
        #         d[str(N) + key] = e
        with open(path_pre + eval + ".txt", "w") as f:
            f.write('N\t')
            for k in models:
                f.write(k + '\t')
            f.write('\n')
            num_genes_list = range(100,1100,100)#+[371]+range(400,1100,100)
            for i in num_genes_list:
                f.write(str(i))
                for key in models:
                    if str(i) + key in d:
                        f.write('\t' + d[str(i) + key])
                    else:
                        f.write('\t' + ' ')
                f.write('\n')

def calculate_weighted_scores():
    evals = [ 'iwavg_cov', 'wavg_mutex', 'wavg_density']
    evals_pos = [1, 2, 4]

    dir_pre = "../hint/out/tmp/evaluation/"
    dir_post = "/optimized_function_comparison/"


    for i in trange(len(evals), desc='Calculating weighte scores'):
        eval = evals[i]
        pos = evals_pos[i]
        with open("../hint/out/tmp/evaluation_tab/"+ eval + ".txt", "w") as eval_file:
            eval_file.write("N\t")
            for model in models:
                eval_file.write(model + "\t")

            eval_file.write("\n")
            num_genes_list = range(100,1100,100)#+[371]+range(400,1100,100)

            for n in num_genes_list:
                eval_file.write(str(n) + "\t")
                for model in models:
                    score = 0.0
                    filepath = dir_pre + model + dir_post + model + "_" + str(n) + ".txt"

                    if os.path.exists(filepath):
                        with open(dir_pre + model + dir_post + model + "_" + str(n) + ".txt") as n_file:
                            lines = n_file.readlines()[1:-2]
                            if pos == 1:
                                score = 0.0
                                genes = 0
                                weight_sum = 0
                                for line in lines:
                                    line = line.strip().split()
                                    weight_sum += (1-(float(line[0])/n))
                                    score += ((1-(float(line[0])/n))) * float(line[pos])
                                    #score += float(line[pos])/float(line[0])
                                    genes += float(line[0])
                                score /= weight_sum
                            else:
                                score = 0.0
                                genes = 0
                                for line in lines:
                                    line = line.strip().split()
                                    score += float(line[0])*float(line[pos])
                                    genes += float(line[0])
                                score /= genes
                        eval_file.write(str(score) + '\t')
                    else:
                        eval_file.write(' ' + '\t')

                eval_file.write("\n")

def prep_file_paths(key):
    paths_raw = glob.glob(key)
    #print ('printing inside prep')

    dict_paths = {}
    for filename in paths_raw:

          directory, file_core = os.path.split(filename)
          #file_core = filename[filename.rindex('/')+1:]
          #print file_core
          if file_core.startswith('cc'):
              num_genes = int(file_core.split('_')[1][1:])
              dict_paths[num_genes] = filename

    sorted_dict = sorted(dict_paths.items(), key=operator.itemgetter(0))

    paths = []
    for i in range(len(sorted_dict)):
        paths.append(sorted_dict[i][1])
    return paths

models = [
    "hotnet2",
    "mutex", "cov", "mutex_cov", \
    'mutex_t07_nsep_cov_nsep',\
    'mutex_t07_nsep_cov_nsep_k6',\
    'mutex_t07_nsep_cov_nsep_k9',\
    'mutex_t07_nsep_cov_nsep_k12'\
    # "mutex_wesme", "mutex_wesme_cov", \
    # "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
    # "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \
    # "mutex_t10_cov", \
    #
    # #these have common threshold > 0.0002
    # "mutex_t05_ncomb_cov", "mutex_t05_nsep_cov", "mutex_t06_ncomb_cov",
    # "mutex_t06_nsep_cov", \
    # "mutex_t07_ncomb_cov", "mutex_t07_nsep_cov", "mutex_t08_ncomb_cov", "mutex_t08_nsep_cov", \
    # "mutex_t09_ncomb_cov", "mutex_t09_nsep_cov",\
    #
    # #these have threshold >0.00007
    # "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", \
    # "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov_nsep", \
    #  "mutex_t07_ncomb_cov_ncomb",  "mutex_t07_nsep_cov_nsep", \
    #  "mutex_t08_ncomb_cov_ncomb",  "mutex_t08_nsep_cov_nsep", \
    #  "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", \
    #
    # "mutex_t05_ncomb_cov_nsep", "mutex_t05_nsep_cov_ncomb", "mutex_t06_ncomb_cov_nsep", "mutex_t06_nsep_cov_ncomb", \
    # "mutex_t07_ncomb_cov_nsep", "mutex_t07_nsep_cov_ncomb", "mutex_t08_ncomb_cov_nsep", "mutex_t08_nsep_cov_ncomb", \
    # "mutex_t09_ncomb_cov_nsep", "mutex_t09_nsep_cov_ncomb",  \
    #
    # #Memcover modules
    # "memcover0.2",
    # "memcover0.08"
     ]
# Read the file once, instead of reading it in a loop many times
cosmic_genes = read_cosmic_genes()
for key in tqdm([]):

    our_path_key = key #'mutex_t05_ncomb_cov'
    #if key == 'hotnet2': continue

    hotnet_paths = prep_file_paths('../hint/out/connected_components_isolarge/hotnet2/*.txt')
    our_paths = prep_file_paths('../hint/out/tmp/connected_components_isolarge/' + our_path_key +'/*.txt')
    #print ('done printing paths')
    print('our paths: ', len(our_paths))
    # print('our paths: ', our_paths)


    newpath = '../hint/out/tmp/evaluation/{0}'.format(our_path_key)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
        newpathext = newpath + '/optimized_function_comparison'
        if not os.path.exists(newpathext):
            os.mkdir(newpathext)

    newpath = '../hint/out/tmp/evaluation_tab/{0}'.format(our_path_key)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
    # os.system('mkdir ' + '../hint/out/evaluation/'+ our_path_key)
    # os.system('mkdir ' + '../hint/out/evaluation/'+ our_path_key + '/optimized_function_comparison')

    hotnet_list = []
    our_list = []
    edge_list_original = load_edge_list()
    fhout = open('../hint/out/tmp/evaluation/'+ our_path_key + '/optimized_function_comparison/summary_' + our_path_key  +'.txt', 'w')
    num_genes_list = range(100,1100,100)#+[371]+range(400,1100,100)
    #print(len(num_genes_list), len(our_paths)+1)

    skip_count = 0 # this is a trick to handle the invariance in the number of files
    for i in trange(min(10, len(our_paths)), desc='running main evaluation'):

        num_genes = str(num_genes_list[i])
        if (key == "hotnet2" or key == 'memcover0.08') and (num_genes == '371'):
            skip_count +=1
            print('key is {} and num_genes is {}'.format(key, num_genes))
            continue
        # COSMIC analysis
        module_name = our_path_key + '/cosmic_' + our_path_key
        [our_overlap, hotnet2_overlap] = cosmic_overlap_analysis(our_paths[i-skip_count], hotnet_paths[i], cosmic_genes, module_name, num_genes)
        print >>fhout, our_overlap, hotnet2_overlap,

        outfile = '../hint/out/tmp/evaluation/'+ our_path_key + '/optimized_function_comparison/' + our_path_key + '_' + str(num_genes) + '.txt'
        outfile_tab = '../hint/out/tmp/evaluation_tab/'+ our_path_key + '/our_eval_' + our_path_key + '.txt'

        our_scores = calculate_scores_for_subnetworks(our_paths[i-skip_count], edge_list_original, outfile, outfile_tab, num_genes)

        for j in range(len(our_scores)):
            print >>fhout, our_scores[j],# hotnet_scores[j],
        print >>fhout, '\n',


    fhout.close()

generate_cosmic_analysis_file()
generate_our_eval_files()
calculate_weighted_scores()
