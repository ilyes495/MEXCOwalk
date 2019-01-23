#
#
# list of models
# for each modele I have a list of files n=100,200,...
# for each file there are components

from glob import glob
import pandas as pd
import numpy as np
from tqdm import tqdm, trange
random_state = 1234
from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import accuracy_score
from matplotlib import pyplot as plt
#from matplotlib import rcParams
from pylab import *
import operator
import os

# params = {
#    'axes.labelsize': 8,
#    # 'text.fontsize': 8,
#    'legend.fontsize': 10,
#    'xtick.labelsize': 10,
#    'ytick.labelsize': 10,
#    'text.usetex': False,
#    'figure.figsize': [7, 4] # instead of 4.5, 4.5
#    }
# rcParams.update(params)

class Knearst ():
    def __init__(self, models_path, list_models, csv_path, gene_file, cv=10,  n_neighbors = 1,):
        # list_modules: list, list_models
        # models_path: string, the path to the root folder of connected_components_isolarge
        # csv_path: string, path to the csv file contianing all the cancers

        self.models = list_models
        self.models_path = models_path
        self.csv_path = csv_path
        self.gene_file = gene_file
        self.cv = cv
        self.n_neighbors = n_neighbors
        self.mappings(self.gene_file)

    def get_file(self,path_model, n):
        list = glob(path_model+'/cc_n{}_*.txt'.format(n))
        assert len(list) > 0, 'file was not found {}/cc_n{}_*.txt'.format(path_model, n)
        return list[0]

    def get_modules(self, file):
        delemiter = '\t' if 'memcover' in file or 'hier' in file else ' '
        with open(file, 'r') as f:
            modules = [l.rstrip().split(delemiter) for l in f.readlines() ]
        # print('len of modules', len(modules))
        # print(modules[:10])
        # genes = []
        # for m in modules:
        #     genes.extend(m)
        # print(len(genes))

        return modules

    def mappings(self, gene_file):
        with open(gene_file, 'r') as f:
            self.gene2id = {g.rstrip():id_ for id_, g in enumerate(f.readlines())}
            self.id2gene = {id_:g for id_, g in self.gene2id.items()}
        # print('Mapping DONE!')

    def run_knearest(self,model, n):

        file  = self.get_file(self.models_path+model, n)
        modules = self.get_modules(file)
        #print('There are {} modules in {}_{}'.format(len(modules), model, n))
        modules = self.check_module(modules)
        #print('modules checked, there are {} modules'.format(len(modules)))
        self.set_data(modules)

        scores_per_module = []
        for i in trange(len(modules)):
            module = modules[i]
            if len(module) == 0:
                print('a module was skipped')
                continue
            score = self.knearst(module)
            scores_per_module.append(score)
        assert i == len(modules)-1, 'not all modules are processed'
        #print('The median of model {} is {}\n'.format(model,np.median(scores_per_module)))
        self.write_log(scores_per_module, model, n)
        return np.mean(scores_per_module)

    def set_data(self,modules):
        # laod the data with all the genes that are present in our files
        # for a given n (number of genes)
        ms = []
        for m in modules:
            ms.extend(m)
        self.data = pd.read_csv(self.csv_path, usecols =ms+['y'])
        #print('Data loaded!')

    def check_module(self, module):
        # it take a list of genes, or list[list[]] of genes,
        # and return check if they exisit in pancancer dataset or not,
        # if not, the gene is removed from the list
        module_ = []
        if isinstance(module[0], list):
            ms= []
            for m in module:
                m_ = self.check_module(m)
                if m_ != '':
                    ms.append(m_)
            return ms
        if isinstance(module, str):
                print(module)
                return [module] if module in self.gene2id.keys() else ''
        for m in module:
            if m in self.gene2id.keys():
                module_.append(m)

        if len(module_) != len(module):
            print('Some genes are not present in the dataset')
            #print('New module: ', module_)
        return module_

    def get_data(self, module):
        # return the data for the genes given in module
        data = self.data[module+['y']]
        return data.sample(frac = 1.0) # This is for shuflling


    def knearst(self, module):

        data = self.get_data(module)
        y = data['y'].values
        X = data[module].values

        kf = StratifiedKFold(n_splits=self.cv, random_state=random_state, shuffle=True)
        score_per_fold = []
        for train_index, test_index in kf.split(X,y):
               X_train, X_test = X[train_index], X[test_index]
               y_train, y_test = y[train_index], y[test_index]
               cls  = NearestNeighbors(self.n_neighbors)
               #print(cls)
               # fit the model
               cls.fit(X_train)
               # print('X_train shape: ', X_train.shape)
               # print('X_test shape: ', X_test.shape)

               # return a list of the indices of nearest neighbor for each point in the test set
               neighbors = cls.kneighbors(X_test, return_distance=False)
               # print('neighbors shape: ', neighbors.shape)
               # print('neighbors: ', neighbors[:10])

               # get the label correspoding to the closest point from train set

               ## TODO:  make it work for any number of neighbors
               predict = np.array([y_train[idx[0]] for idx in neighbors])
               # print ('predict shape: ', predict.shape)
               # print ('y_test shape: ', y_test.shape)

               score  = accuracy_score(y_test.reshape(-1,1), predict)
               score_per_fold.append(score)

        return np.mean(score_per_fold)

    def run(self, min_n, max_n):

        for model in self.models:
            score_per_n = []
            t = trange(min_n, max_n, 100)
            for n in t:
                #print('n: ', str(n)+'\n')
                t.set_description('# max gene: {}'.format(n))
                t.refresh()

                score = self.run_knearest(model, n)
                score_per_n.append(score)
            #self.plot(score_per_n,min_n, max_n, model)


    def write_log(self, scores, model, n):

        if not os.path.exists('../results/knearst/{}/'.format(model)):
            os.mkdir('../results/knearst/{}/'.format(model))

        with open('../results/knearst/{}/{}.txt'.format(model, n), 'w') as f:
            f.write('median {}\n'.format(np.median(scores)))
            f.write('max {}\n'.format(np.max(scores)))
            f.write('min {}\n'.format(np.min(scores)))
            f.write('mean {}\n'.format(np.mean(scores)))

    def plot(self, scores,min_n, max_n, model):
        plt.xlabel('N')
        plt.ylabel('Accuracy')
        #plt.title('Median of average cv accuracy per module for {}'.format(model))
        plt.plot(range(min_n, max_n, 100),scores)


if __name__ == '__main__':

    models = [
        "hotnet2",
        "memcover_v1",
        "memcover_v2",
        "mutex_t07_nsep_cov_nsep",
        "mutex_t07_nsep_cov",
        "memcover_v3",
        "hier_hotnet"
        # "mutex", "mutex_cov", #"cov",\
        # "mutex_wesme", "mutex_wesme_cov", \
        # "mutex_ncomb", "cov_ncomb", "mutex_ncomb_cov", \
        # "mutex_nsep", "cov_nsep", "mutex_nsep_cov", \
        # #these have threshold >0.00007
        # "mutex_t05_ncomb_cov_ncomb", "mutex_t05_nsep_cov_nsep", \
        # "mutex_t06_ncomb_cov_ncomb", "mutex_t06_nsep_cov_nsep", \
        # "mutex_t07_ncomb_cov_ncomb",  "mutex_t07_nsep_cov_nsep", \
        # "mutex_t08_ncomb_cov_ncomb",  "mutex_t08_nsep_cov_nsep", \
        # "mutex_t09_ncomb_cov_ncomb", "mutex_t09_nsep_cov_nsep", \
        #
        # #these have common threshold > 0.0002
        # "mutex_t05_ncomb_cov", "mutex_t05_nsep_cov", "mutex_t06_ncomb_cov",
        # "mutex_t06_nsep_cov", \
        # "mutex_t07_ncomb_cov", "mutex_t07_nsep_cov", "mutex_t08_ncomb_cov", "mutex_t08_nsep_cov", \
        # "mutex_t09_ncomb_cov", "mutex_t09_nsep_cov",\
        # "mutex_t10_cov", \
        # "mutex_t05_ncomb_cov_nsep", "mutex_t05_nsep_cov_ncomb", "mutex_t06_ncomb_cov_nsep", "mutex_t06_nsep_cov_ncomb", \
        # "mutex_t07_ncomb_cov_nsep", "mutex_t07_nsep_cov_ncomb", "mutex_t08_ncomb_cov_nsep", "mutex_t08_nsep_cov_ncomb", \
        # "mutex_t09_ncomb_cov_nsep", "mutex_t09_nsep_cov_ncomb"  \

         ]
    # models_path = '../hint/out/connected_components_isolarge_n2500_whh/'
    # csv_path = '../data/pancancer_all_subtypes.csv'
    # gene_file = '../data/hint_inters_pan.txt'
    # Knearst = Knearst(models_path, models, csv_path, gene_file, cv=10)
    # # special Ns, 554, 806
    # min_n = 806
    # max_n = 806+10
    # Knearst.run(min_n, max_n)
    plt.xlabel('N')
    plt.ylabel('Accuracy')
    axes(frameon=0)
    median_idx = 0
    mean_idx = 3
    for model in models:
        if model == 'hier_hotnet': continue
        median_scores = []
        mean_scores = []
        a = list(range(100,600,100))+[554]+ list(range(600,900, 100))+[806]
        N = a + list(range(900, 1700, 100)) if model  == 'memcover_v3' else a+list(range(900, 2600, 100))
        for n in N:
            results_file = '../results/knearst/{}/{}.txt'.format(model, n)
            with open(results_file, 'r') as f:
                lines= f.readlines()
                mean = float(lines[mean_idx].rstrip().split()[-1])
                mean_scores.append(mean)
                median = float(lines[median_idx].rstrip().split()[-1])
                median_scores.append(median)
        # print(mean_scores)
        plt.plot(N,median_scores, '-o')


median_scores = []
mean_scores = []
for n in [554,806]:
    results_file = '../results/knearst/{}/{}.txt'.format('hier_hotnet', n)
    with open(results_file, 'r') as f:
        lines= f.readlines()
        mean = float(lines[mean_idx].rstrip().split()[-1])
        mean_scores.append(mean)
        median = float(lines[median_idx].rstrip().split()[-1])
        median_scores.append(median)
# print(mean_scores)
plt.plot([554,806],median_scores, 'k*',)

art = []
legend = plt.legend(models, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                    edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.4))
art.append(legend)
frame = legend.get_frame()
# frame.set_facecolor('0.9')
# frame.set_edgecolor('0.9')
xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
xticks(xtick)
plt.ylim(0.79, 1)
# axes(frameon=0)
grid()
plt.savefig('../results/knearst/plots/Median_Accuracy.png',additional_artists=art,
            bbox_inches="tight")
plt.close()
