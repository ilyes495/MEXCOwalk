import os.path
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
    "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", "mutex_t10_ncomb_cov", "mutex_t10_nsep_cov",\
    "mutex_t10_ncomb_cov_ncomb", "mutex_t10_nsep_cov_nsep", "mutex_t10_nsep_cov_ncomb", "mutex_t10_ncomb_cov_nsep"\


    # "mutex_nsep_t06_a07_cov_ncomb_d1", "mutex_nsep_t06_a07_cov_ncomb_d2", "mutex_nsep_t06_a07_cov_ncomb_d3",\
    # "mutex_nsep_t07_a07_cov_ncomb_d1", "mutex_nsep_t07_a07_cov_ncomb_d2", "mutex_nsep_t07_a07_cov_ncomb_d3",\
    # "mutex_nsep_t08_a07_cov_ncomb_d1", "mutex_nsep_t08_a07_cov_ncomb_d2", "mutex_nsep_t08_a07_cov_ncomb_d3"
     ]
# models = [
# "hotnet2", \
# "mutex", "cov", "mutex_cov", \
# "mutex_wesme", "mutex_wesme_cov", \
# "mutex_mod", "cov_mod", "mutex_mod_cov", \
# "mutex_avgmod", "cov_avgmod", "mutex_avgmod_cov", \
#
# "mutex_t10_cov", \
# "mutex_t05_mod_cov", "mutex_t05_mod_cov_mod", "mutex_t05_avgmod_cov", "mutex_t05_avgmod_cov_avgmod", \
# "mutex_t06_mod_cov", "mutex_t06_mod_cov_mod", "mutex_t06_avgmod_cov", "mutex_t06_avgmod_cov_avgmod", \
# "mutex_t07_mod_cov", "mutex_t07_mod_cov_mod", "mutex_t07_avgmod_cov", "mutex_t07_avgmod_cov_avgmod", \
# "mutex_t08_mod_cov", "mutex_t08_mod_cov_mod", "mutex_t08_avgmod_cov", "mutex_t08_avgmod_cov_avgmod", \
# "mutex_t09_mod_cov", "mutex_t09_mod_cov_mod", "mutex_t09_avgmod_cov", "mutex_t09_avgmod_cov_avgmod"
# ]

evals = [ 'wavg_cov', 'wavg_mutex', 'wavg_density']
evals_pos = [1, 2, 4]

dir_pre = "../hint/out/evaluation/"
dir_post = "/optimized_function_comparison/"


for i in range(len(evals)):
    eval = evals[i]
    pos = evals_pos[i]
    with open("../hint/out/evaluation_tab/"+ eval + ".txt", "w") as eval_file:
        eval_file.write("n\t")
        for model in models:
            eval_file.write(model + "\t")
        eval_file.write("\n")
        for n in range(100, 2600, 100):
            eval_file.write(str(n) + "\t")
            for model in models:
                score = 0.0
                filepath = dir_pre + model + dir_post + model + "_" + str(n) + ".txt"
                if os.path.exists(filepath):
                    with open(dir_pre + model + dir_post + model + "_" + str(n) + ".txt") as n_file:
                        lines = n_file.readlines()[1:-2]
                        score = 0.0
                        genes = 0
                        for line in lines:
                            line = line.strip().split()
                            score += float(line[0])*float(line[pos])
                            genes += float(line[0])
                        score /= genes

                eval_file.write(str(score) + '\t')
            eval_file.write("\n")
