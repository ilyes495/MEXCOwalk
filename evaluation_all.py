from utils import *
import glob
import os
import math
import operator
from tqdm import trange, tqdm
from matplotlib import pyplot as plt
#from matplotlib import rcParams
from pylab import *

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

        total_genes = 0
        print >>fhout, 'num_comp\tcoverage\tmutex\tcov_mutex\tdensity'
        for component_line in connected_components:
            component = component_line.split()
            coverage_list = []
            union_set = set([])
            individual_count = 0
            mult_count = 1
            for i in component:
                total_genes+=1
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
    print >>fhout, str(total_genes) + '\t' + str(num_components) + \
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

    cosmic_file = '../hint/out/evaluation/{}.txt'.format(module_name)
    cosmic_huge_file = '../hint/out/evaluation/{}_huge.txt'.format(module_name)
    cosmic_tab_file = '../hint/out/evaluation_tab/{}.txt'.format(module_name)

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

def generate_cosmic_analysis_file(path_pre = "../hint/out/evaluation_tab/"):
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
        num_genes_list = range(100,600,100)+[554]+range(600,900,100)+[806]+range(900,2600,100)
        for i in num_genes_list:
            ff.write(str(i))
            for key in models:
                if str(i) + key in d:
                    ff.write('\t' + d[str(i) + key])
                else:
                    ff.write('\t' + ' ')
            ff.write('\n')




def generate_our_eval_files(path_pre = "../hint/out/evaluation_tab/"):
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
            num_genes_list = range(100,600,100)+[554]+range(600,900,100)+[806]+range(900,2600,100)
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

    dir_pre = "../hint/out/evaluation/"
    dir_post = "/optimized_function_comparison/"


    for i in trange(len(evals), desc='Calculating weighte scores'):
        eval = evals[i]
        pos = evals_pos[i]
        with open("../hint/out/evaluation_tab/"+ eval + ".txt", "w") as eval_file:
            eval_file.write("N\t")
            for model in models:
                eval_file.write(model + "\t")

            eval_file.write("\n")
            num_genes_list =range(100,600,100)+[554]+range(600,900,100)+[806]+range(900,2600,100)

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

models  = ["hotnet2","memcover_v1","memcover_v2","memcover_v3","mutex_t07_nsep_cov", 'hier_hotnet2']

# Read the file once, instead of reading it in a loop many times
cosmic_genes = read_cosmic_genes()
for key in tqdm([]):

    our_path_key = key #'mutex_t05_ncomb_cov'
    #if key == 'hotnet2': continue

    hotnet_paths = prep_file_paths('../hint/out/connected_components_isolarge_n2500_whh/hotnet2/*.txt')
    our_paths = prep_file_paths('../hint/out/connected_components_isolarge_n2500_whh/' + our_path_key +'/*.txt')
    #print ('done printing paths')
    # print('our paths: ', len(our_paths))
    # print('our paths: ', our_paths)


    newpath = '../hint/out/evaluation/{0}'.format(our_path_key)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
        newpathext = newpath + '/optimized_function_comparison'
        if not os.path.exists(newpathext):
            os.mkdir(newpathext)

    newpath = '../hint/out/evaluation_tab/{0}'.format(our_path_key)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
    # os.system('mkdir ' + '../hint/out/evaluation/'+ our_path_key)
    # os.system('mkdir ' + '../hint/out/evaluation/'+ our_path_key + '/optimized_function_comparison')

    hotnet_list = []
    our_list = []
    edge_list_original = load_edge_list()
    fhout = open('../hint/out/evaluation/'+ our_path_key + '/optimized_function_comparison/summary_' + our_path_key  +'.txt', 'w')
    num_genes_list = range(100,600,100)+[554]+range(600,900,100)+[806]+range(900,2600,100)
    #num_genes_list_hh = [554]+[806]
    # print(len(num_genes_list), len(our_paths))

    skip_count = 0 # this is a trick to handle the invariance in the number of files

    for i in trange(min(27, len(our_paths)), desc='running main evaluation'):

        num_genes = 806 if key == 'hier_hotnet2_k2' else(554 if  key == 'hier_hotnet2_k3' else  str(num_genes_list[i]))
        #indices file index and hotnet file index

        #if len(our_paths) == len()

        #print "#",num_genes
        # if (key == "hotnet2" or key == 'hier_hotnet_k2' or key == 'hier_hotnet_k3') and (num_genes == '554' or num_genes == '806'):
        #     skip_count +=1
        #     print('key is {} and num_genes is {}'.format(key, num_genes))
        #     continue
        # # COSMIC analysis
        # print "skip_count", skip_count

        hotnet_path = glob.glob('../hint/out/connected_components_isolarge_n2500_whh/hotnet2/cc_n{}_*'.format(num_genes))[0]
        our_path = glob.glob('../hint/out/connected_components_isolarge_n2500_whh/{}/cc_n{}_*'.format(key,num_genes))[0]


        module_name = our_path_key + '/cosmic_' + our_path_key
        [our_overlap, hotnet2_overlap] = cosmic_overlap_analysis(our_path, hotnet_path, cosmic_genes, module_name, num_genes)
        print >>fhout, our_overlap, hotnet2_overlap,

        outfile = '../hint/out/evaluation/'+ our_path_key + '/optimized_function_comparison/' + our_path_key + '_' + str(num_genes) + '.txt'
        outfile_tab = '../hint/out/evaluation_tab/'+ our_path_key + '/our_eval_' + our_path_key + '.txt'

        our_scores = calculate_scores_for_subnetworks(our_path, edge_list_original, outfile, outfile_tab, num_genes)

        for j in range(len(our_scores)):
            print >>fhout, our_scores[j],# hotnet_scores[j],
        print >>fhout, '\n',


    fhout.close()

# models  = ["hotnet2","memcover_v1","memcover_v2","memcover_v3","mutex_t07_nsep_cov", 'hier_hotnet2']

# generate_cosmic_analysis_file()
# generate_our_eval_files()
# calculate_weighted_scores()

# models  = ["hotnet2","memcover_v1","memcover_v2","memcover_v3","mutex_t07_nsep_cov", 'hier_hotnet2']
ps = ['iwavg_cov','wavg_mutex','wavg_covmutex']
labels = ['Coverage Score (CS)', 'Mutual Exclusion Score ( MS )',  'Driver Module Set Score (DMSS)' ]
for p,l in zip(ps,labels):
    Ns= []
    with open('../hint/out/evaluation_tab/{}.txt'.format(p)) as f:
        lines = f.readlines()
        models_ = lines[0].rstrip().split('\t')[1:]
        model2w = {s:[] for s in models_}
        for line in lines[1:]:
            line = line.rstrip().split('\t')
            Ns.append(int(line[0]))
            for m,w in zip(models_, line[1:]):
                if float(w) == 0: continue
                try:
                    model2w[m].append(float(w))
                except:
                    print(m,w)
    print(l)
    for m in models:
        print(m,'\n',model2w[m])

    axes(frameon=0)
    for m in models:
        if m == 'hier_hotnet2':
            plot([554,806], model2w[m], 'k*', markersize=12)
        elif m == 'memcover_v3':
            plot(Ns[:-9], model2w[m], '-o')
        else:
            try:
                plot(Ns, model2w[m], '-o')
            except:
                print(m, len(Ns), len(model2w[m]))

    art = []
    legend_ = [
        "Hotnet2",
        "MEMCover_v1",
        "MEMCover_v2",
        "MEMCover_v3",
        "MEXCOwalk",
        "Hierarchical Hotnet",
        ]

    legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                        edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.3))
    art.append(legend)
    frame = legend.get_frame()
    plt.xlabel('total_genes')
    plt.ylabel(l)

    xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
    xticks(xtick, fontsize='x-small')
    if p == 'iwavg_cov':
        plt.ylim(-0.01, 0.21)
    elif p == 'wavg_mutex':
        # plt.ylim(0.01, 1.)
        pass
        # yticks(list(range(1, 10, 1))*0.1)
    else:
        plt.ylim(-0.01, 0.18)
    # axes(frameon=0)
    grid()
    plt.savefig('../hint/out/evaluation_tab/plots/{}.pdf'.format(p), format= 'pdf',transparent = True, additional_artists=art,
                bbox_inches="tight", dpi = 800)
    plt.close()
