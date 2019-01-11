# models with hotnet
keys = [
    "hotnet2",
    #"cluster_one", \
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
    #"mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep"\


    # "mutex_nsep_t06_a07_cov_ncomb_d1", "mutex_nsep_t06_a07_cov_ncomb_d2", "mutex_nsep_t06_a07_cov_ncomb_d3",\
    # "mutex_nsep_t07_a07_cov_ncomb_d1", "mutex_nsep_t07_a07_cov_ncomb_d2", "mutex_nsep_t07_a07_cov_ncomb_d3",\
    # "mutex_nsep_t08_a07_cov_ncomb_d1", "mutex_nsep_t08_a07_cov_ncomb_d2", "mutex_nsep_t08_a07_cov_ncomb_d3"
     ]


path_pre = "../hint/out/evaluation_tab/"

# GO
def generate_go_analysis_file():
    d = {}
    for key in keys:
        if key == "hotnet2":
            continue
        with open(path_pre + key + "/GO_ne_" + key + ".txt", "r") as f:
            lines = f.readlines()
            for i in range(0, len(lines), 6):
                N = lines[i].strip()
                our = lines[i+1].split()[0].strip()
                d[str(N) + key] = our
    print "d", d
    key = "hotnet2"
    with open(path_pre + keys[1] + "/GO_ne_" + keys[1] + ".txt", "r") as f:
        lines = f.readlines()
        for i in range(0, len(lines), 6):
            N = lines[i].strip()
            hotnet = lines[i+4].split()[0].strip()
            d[str(N) + key] = hotnet

    with open(path_pre + "go" + ".txt", "w") as ff:
        ff.write('N\t')
        for k in keys:
            ff.write(k + '\t')
        ff.write('\n')
        for i in range(100, 2600, 100):
            ff.write(str(i))
            for key in keys:
                if str(i) + key in d:
                    ff.write('\t' + d[str(i) + key])
                else:
                    ff.write('\t' + '0')
            ff.write('\n')

# COSMIC
def generate_cosmic_analysis_file():
    d = {}
    for key in keys:
        print key + "\n"
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
    with open(path_pre + keys[1] + "/cosmic_" + keys[1] + ".txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            N = line[0].strip()
            hotnet = line[2].strip()
            d[str(N) + key] = hotnet

    with open(path_pre + "cosmic" + ".txt", "w") as ff:
        ff.write('N\t')
        for k in keys:
            ff.write(k + '\t')
        ff.write('\n')
        for i in range(100, 2600, 100):
            ff.write(str(i))
            for key in keys:
                if str(i) + key in d:
                    ff.write('\t' + d[str(i) + key])
                else:
                    ff.write('\t' + '0')
            ff.write('\n')


eval_list = ['num_components', \
't_density', 'avg_density', \
't_coverage', 'avg_cov', \
't_mutex', 'avg_mutex', \
't_covmutex', 'avg_covmutex']

def generate_our_eval_files():
    for i in range(len(eval_list)):
        eval = eval_list[i]
        d = {}
        for key in keys:
            print(key)
            if key == "hotnet2":
                continue
            with open(path_pre + key + "/our_eval_" + key + ".txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.split()
                    N = line[0].strip()
                    e = line[i+1].strip()
                    d[str(N) + key] = e
        key = "hotnet2"
        with open(path_pre + keys[1] + "/our_eval_hotnet2.txt", "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                N = line[0].strip()
                e = line[i+1].strip()
                d[str(N) + key] = e
        with open(path_pre + eval + ".txt", "w") as f:
            f.write('N\t')
            for k in keys:
                f.write(k + '\t')
            f.write('\n')
            for i in range(100, 2600, 100):
                f.write(str(i))
                for key in keys:
                    if str(i) + key in d:
                        f.write('\t' + d[str(i) + key])
                    else:
                        f.write('\t' + '0')
                f.write('\n')

# generate_go_analysis_file()
generate_cosmic_analysis_file()
generate_our_eval_files()
