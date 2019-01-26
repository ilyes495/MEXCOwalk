from scipy.stats import norm
from lifelines import CoxPHFitter
# %matplotlib inline
import os
from os import path
from glob import glob
import numpy as np
import pandas as pd

from utils import *
import sys, io
import warnings

# for each module
# for significant cancer:
#     get clinical dataset
#     get gene expression [moule]
#     run coxph


cancer_types = ['BLCA','BRCA','CRC', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC']
#cancer_types= ['GBM']
out_path = '../hint/out/'
data_path = '../data/'
data_files = 'pancancer/'
clinic_path = data_path+'clinical_data/'
clinic_file = '{}_clinical_data.txt'
files = glob(data_path+data_files+'*/[!MANIFEST]*.txt')
cancer2file = {c.split('/')[-1].split('.')[0]:c for c in files}

with open('../data/exp_final_gene_list.txt', 'r') as f:
    final_genes= [l.rstrip() for l in f.readlines()]
with open('../data/hint_inters_pan.txt', 'r') as f:
    hint_u_pan= [l.rstrip() for l in f.readlines()]

def get_modules_cancer(model ='mutex_t07_nsep_cov_nsep' , n = 100):
    comps = []
    significant_cancers =[]
    with open('../hint/out/cancer_subtype_test_results/{}/cc_n{}_cancer_subtype_tests.txt'.format(model, n), 'r') as f:
        lines = f.readlines()[:]
        for i in range(0,len(lines), 12):
            comp = lines[i].strip().split('[')[-1].split(',')
            comp = [c.replace(']', '').replace('\'', '').strip() for c in comp]
            comps.append(comp)
            type_interest = [j-1 for j in range(1,12) if len(lines[i+j].split()) ==8]
            significant_cancers.append(type_interest)
    print('get_modules_cancer Done!')
    return  comps, significant_cancers

def get_data(modules, cancer ):
    clinical_data = clinic_file.format(cancer)
    df_clinic = pd.read_csv(clinic_path+clinical_data, usecols = ['patient_id','age', 'gender','event','time_to_event'], delimiter = '\t', low_memory=False)
    # patients = df_clinic['patient_id'].values
    df_clinic = df_clinic.set_index('patient_id')
    df_clinic.index.names = ['patients']
    df_clinic =df_clinic[~df_clinic.isnull().any(axis=1)]

    df = pd.read_csv(cancer2file[cancer], delimiter = '\t', low_memory=False)
    df.drop(df.index[0], inplace=True)
    gene_and_id = df.iloc[:, 0].values
    df.index = final_genes
    df.drop(df.columns[0], axis=1, inplace= True)
    df = df.loc[list(hint_u_pan)]
    df = df.T[modules].astype('float')
    # df.index.names = ['patients']
    df.columns.name = 'genes'

    if cancer != 'LAML':
        keep_patient = [x.split('-')[3][:2] == '01' for x in df.index ]
    else:
        keep_patient = [x.split('-')[3][:2] == '03' for x in df.index]

    df = df.loc[keep_patient]
    df.index = ['-'.join(x.split('-')[:3]) for x in df.index]
    df.index.names = ['patients']

    df = df.apply(invers_norm, axis = 0)

    df2 = df_clinic.join(df)
    df2 = df2.loc[~df2[df2.columns[4:]].isnull().any(axis=1)]

    print('Data loaded!')
    return df2

def invers_norm(data):
    rank = np.argsort(data)
    q = [0 for _ in range(len(rank))]
    for i,k in enumerate(rank):
        q[k] = i+1
    q = (np.array(q)-0.5) /len(q)
    return norm.ppf(q)

results_path = '../results/CoxPh/'
comps, significant_cancers = get_modules_cancer()
results = {}

cancer2modules = {}
for i,m in enumerate(comps):
    for t in significant_cancers[i]:
        if cancer_types[t] in cancer2modules.keys():
            cancer2modules[cancer_types[t]].append(m)
        else:
            cancer2modules[cancer_types[t]] = [m]

module2id = {str(m):i for i,m in enumerate(comps)}
id2module = {i:str(m) for i,m in enumerate(comps)}

all_genes = [g for c in comps for g in c]

for cancer in cancer_types:
    if not cancer in cancer2modules.keys() :
        print('{} skipped'.format(cancer))
        continue

    print('Cancer {}'.format(cancer))
    if cancer == 'CRC':
        df = get_data(all_genes, 'COADREAD')
    else:
        df = get_data(all_genes, cancer)


    for module in cancer2modules[cancer]:

        module_id = module2id[str(module)]
        if not os.path.exists(results_path+'{}/'.format(module_id)):
            os.mkdir(results_path+'{}/'.format(module_id))
        csv_path = results_path+'{}/{}.csv'.format(module_id, cancer)

        if not cancer in ['OV', 'UCEC']:
            tmp_df  = df[module+['age', 'gender','event','time_to_event']]
            cph = CoxPHFitter()
            cph.fit(tmp_df, duration_col='time_to_event', event_col='event', show_progress=False)
        else:
            tmp_df  = df[module+['age','event','time_to_event']]
            cph = CoxPHFitter()
            cph.fit(tmp_df, duration_col='time_to_event', event_col='event', show_progress=False)
        s = cph.summary

        # capture the output of print_summary()
        stdout = sys.stdout
        sys.stdout = io.StringIO()
        cph.print_summary()
        output = sys.stdout.getvalue()
        sys.stdout = stdout

        if str(module) in results.keys():
            results[str(module)].append('Cancer: {} \n'.format(cancer)+str(output))
        else:
            results[str(module)] = ['Cancer: {} \n'.format(cancer)+str(output)]

        s.to_csv(csv_path)

for m in module2id.keys():
    if not m in results.keys(): continue
    module_id = module2id[str(m)]
    txt_path = results_path+'{}/summary.txt'.format(module_id)
    with open(txt_path, 'w') as f:
        f.write('\n\n'.join(results[m]))

# for cancer in cancer_types:
#     if not cancer in cancer2modules.keys(): continue
#     for module in cancer2modules[cancer]:
#
