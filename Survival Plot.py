
# coding: utf-8

# In[6]:

import numpy as np
import pandas as pd
import lifelines as ll

get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font = fm.FontProperties(size=25)

#larger
from pylab import rcParams
rcParams['figure.figsize'] = 15, 10

from lifelines.estimation import KaplanMeierFitter
from lifelines import CoxPHFitter


# In[7]:

# Breast is run 10 and K = 3
# Colon is run 1 and K = 5
# Kidney is run 5 and K = 4
# Lung run is 4 and K = 8 | run is also 4 for K = 3


# In[ ]:




# In[23]:

cancer_type = 'Kidney'
if cancer_type == 'Lung':
    cancer_name = 'LSCC'
    rand_run = 4
    cluster_num = 3
elif cancer_type == 'Breast':
    cancer_name = 'BIC'
    rand_run = 10
    cluster_num = 3
elif cancer_type == 'Colon':
    cancer_name = 'COAD'
    rand_run = 1
    cluster_num = 5
elif cancer_type == 'Kidney':
    cancer_name = 'KRCC'
    rand_run = 5
    cluster_num = 4
    

clinical_filename = '/home/mllab/Dropbox/Tunde/analysis/new_datasets/' + cancer_type + '/' + cancer_type +  '_clinical_data.txt'
drug_filename = '/home/mllab/Dropbox/Tunde/analysis/new_datasets/' + cancer_type + '/' + cancer_type +  '_drug_data.txt'
cluster_label = '/home/mllab/Dropbox/Tunde/analysis/new_results/MVKKM_newdata/labels/' + cancer_type + '_mvkkm_labels_K' +  str(cluster_num) + '_random_init_run' + str(rand_run) + '_rbf_kernel_old_gamma_with_shankar_distance_second_derv.txt'
pam50_file = '/home/mllab/Dropbox/Tunde/analysis/PAM50/bioclassifier_R/Breastcancerpan50subtyperesult_pam50scores.txt'


# In[24]:

cluster_l = open(cluster_label).readlines()
cluster_list = [s.rstrip() for s in cluster_l]
np.array(cluster_list)
cluster_list



# In[25]:

df = pd.read_csv(clinical_filename, delimiter='\t')
df['group'] = cluster_list
df.head()


# In[26]:

kmf = KaplanMeierFitter()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'Arial'

for group in sorted(df['group'].unique()):
    g = df.group == group
    T = df[g]['days_to_last_followup']
    C = df[g]['event']
    kmf.fit(T, event_observed=C, label='Cluster - ' + group + ' (' + str(len(T)) + ')')
    kmf.survival_function_.plot(ax=ax,  linewidth=4.0)
kmf2 = plt.gcf()
plt.title(cancer_name,fontsize=30)
plt.xlabel('Time in Days',fontsize=30)
plt.ylabel('Survival Rate',fontsize=30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.xlim(xmin=0)
plt.legend(prop=font,frameon=False)
plt.savefig(cancer_type + '_survival_plot_improved.eps',format='eps')








# In[3]:

cox_data = df[['time_to_event','event','gender','age','group']]
cph = CoxPHFitter()
cph.fit(cox_data, 'time_to_event', event_col='event')
cph.confidence_intervals_


# In[125]:

cph.print_summary()


# In[126]:

cph.summary


# In[127]:

cox_data


# In[ ]:




# In[13]:

drug_df = pd.read_csv(drug_filename, delimiter='\t')
drug_df.head()


# In[14]:

drug_df[' Measure_of_Response ']


# In[15]:

df['therapy_type'] = drug_df[' Therapy_Type ']
df['drug_name'] = drug_df[' Drug_name ']
df.head()


# In[16]:

c = df[df['group'] == '1']
len(c[c['gender'] == 0]) > 0


# In[17]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        # ensure that there exist gender data in each cluster
        if len(clust[clust['gender'] == 1]) > 0 and len (clust[clust['gender'] == 1]) > 0 : 
            cluster_male = clust[clust['gender'] == 0]
            cluster_female = clust[clust['gender'] == 1]
            T = cluster_male['days_to_last_followup']
            C = cluster_male['event']
            T2 = cluster_female['days_to_last_followup']
            C2 = cluster_female['event']
        plt.figure()
        ax = plt.subplot(111)

        kmf.fit(T, event_observed=C, label=['Male'])
        kmf.survival_function_.plot(ax=ax)
        kmf.fit(T2, event_observed=C2, label=['Female'])
        kmf.survival_function_.plot(ax=ax)

        plt.title('Lifespans of different gender in Cluster ' + clusters)
        kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass
        


# In[18]:

clust = df[df['group'] == '1']
cluster_chemo = clust[clust['therapy_type'].str.match('chemotherapy')]
cluster_chemo


# In[19]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for therapy in clust['therapy_type'].unique():
            data = clust[clust['therapy_type'].str.match(therapy)]
            T = data['days_to_last_followup']
            C = data['event']
            ax = plt.subplot(111)
            kmf.fit(T, event_observed=C, label=[therapy + '(' + str(len(data)) + ')'])
            kmf.survival_function_.plot(ax=ax)
            plt.title('Lifespans in Cluster ' + clusters + ' Based on Therapy Type')
            kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass


# In[20]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for therapy in ['chemotherapy', 'hormone therapy']:
            data = clust[clust['therapy_type'].str.match(therapy)]
            T = data['days_to_last_followup']
            C = data['event']
            ax = plt.subplot(111)
            kmf.fit(T, event_observed=C, label=[therapy + '(' + str(len(data)) + ')'])
            kmf.survival_function_.plot(ax=ax)
            plt.title('Lifespans in Cluster ' + clusters + ' Based on Therapy Type')
            kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass


# In[21]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for drug in clust['drug_name'].unique():
            data = clust[clust['drug_name'].str.match(drug)]
            T = data['days_to_last_followup']
            C = data['event']
            ax = plt.subplot(111)
            kmf.fit(T, event_observed=C, label=[drug + '(' + str(len(data)) + ')'])
            kmf.survival_function_.plot(ax=ax)
            plt.title('Lifespans in Cluster ' + clusters + ' Based on Drug Used')
            kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass


# In[ ]:




# In[22]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for drug in ['tamoxifen', 'adriamycin', 'doxorubicin', 'NA']:
            data = clust[clust['drug_name'].str.match(drug)]
            T = data['days_to_last_followup']
            C = data['event']
            ax = plt.subplot(111)
            kmf.fit(T, event_observed=C, label=[drug + ' (' + str(len(data)) + ')'])
            kmf.survival_function_.plot(ax=ax)
            plt.title('Lifespans in Cluster ' + clusters + ' Based on Drug Used')
            kmf2 = plt.gcf()
            
except ValueError:  #raised if one gender is empty.
    pass


# In[23]:

df


# In[24]:

clust2 = df[df['group'] == '2']


# In[25]:

clust2


# In[26]:

data = clust2[clust2['drug_name'].str.match('fluorouracil')]
data2 = clust2[~clust2['drug_name'].str.match('fluorouracil')]


# In[27]:

T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['fluorouracil ' + str(len(data))])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others '  + str(len(data2)) ])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 2 Based on Drug Used')
kmf2 = plt.gcf()


# In[28]:

clust2 = df[df['group'] == '1']
data = clust2[clust2['drug_name'].str.match('fluorouracil')]
data2 = clust2[~clust2['drug_name'].str.match('fluorouracil')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['fluorouracil'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 1 Based on Drug Used')
kmf2 = plt.gcf()


# In[29]:

clust2 = df[df['group'] == '3']
data = clust2[clust2['drug_name'].str.match('fluorouracil')]
data2 = clust2[~clust2['drug_name'].str.match('fluorouracil')  & ~clust2['drug_name'].str.match('NA')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['fluorouracil'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 1 Based on Drug Used')
kmf2 = plt.gcf()   


# In[30]:

data2


# In[31]:

clust2 = df[df['group'] == '2']
data = clust2[clust2['drug_name'].str.match('xeloda')]
data2 = clust2[~clust2['drug_name'].str.match('xeloda')  & ~clust2['drug_name'].str.match('NA')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['xeloda'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 1 Based on Drug Used')
kmf2 = plt.gcf()   


# In[32]:

clust2 = df[df['group'] == '1']
data = clust2[clust2['therapy_type'].str.match('chemotherapy')]
data2 = clust2[~clust2['therapy_type'].str.match('chemotherapy')  & ~clust2['therapy_type'].str.match('NA')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['chemotherapy'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 1 Based on Theraphy')
kmf2 = plt.gcf()   


# In[33]:

clust2 = df[df['group'] == '2']
data = clust2[clust2['therapy_type'].str.match('chemotherapy')]
data2 = clust2[~clust2['therapy_type'].str.match('chemotherapy')  & ~clust2['therapy_type'].str.match('NA')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['chemotherapy'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 2 Based on Theraphy')
kmf2 = plt.gcf()  


# In[34]:

clust2 = df[df['group'] == '3']
data = clust2[clust2['therapy_type'].str.match('chemotherapy')]
data2 = clust2[~clust2['therapy_type'].str.match('chemotherapy')  & ~clust2['therapy_type'].str.match('NA')]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=['chemotherapy'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=['others'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans in Cluster 3 Based on Theraphy')
kmf2 = plt.gcf()  


# In[28]:

pam50_data = pd.read_csv(pam50_file,delimiter='\t')
pam50_data.head()


# In[29]:

df['pam_subtype'] = pam50_data['Call']


# In[37]:

try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for therapy in ['chemotherapy', 'hormone therapy']:
            for subtype in clust['pam_subtype'].unique():
                data = clust[(clust['therapy_type'].str.match(therapy)) & (clust['pam_subtype'].str.match(subtype))]
                T = data['days_to_last_followup']
                C = data['event']
                ax = plt.subplot(111)
                kmf.fit(T, event_observed=C, label=[therapy + ' ' + subtype + ' (' + str(len(data)) + ')'])
                kmf.survival_function_.plot(ax=ax)
                plt.title('Lifespans in Cluster ' + clusters + ' Based on Therapy Type')
                kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass


# In[38]:

subtype = 'LumA'
try:
    for clusters in df['group'].unique():
        clust = df[df['group'] == clusters]
        plt.figure()
        for therapy in ['chemotherapy', 'hormone therapy']:
            data = clust[(clust['therapy_type'].str.match(therapy)) & (clust['pam_subtype'].str.match(subtype))]
            T = data['days_to_last_followup']
            C = data['event']
            ax = plt.subplot(111)
            kmf.fit(T, event_observed=C, label=[therapy + ' ' + subtype + ' (' + str(len(data)) + ')'])
            kmf.survival_function_.plot(ax=ax)
            plt.title('Lifespans in Cluster ' + clusters + ' Based on Therapy Type')
            kmf2 = plt.gcf()
except ValueError:  #raised if one gender is empty.
    pass


# In[43]:

subtype = 'LumA'
data = df[(df['group'] == '1') & (df['pam_subtype'].str.match(subtype))]
data2 = df[(df['group'] == '3') & (df['pam_subtype'].str.match(subtype))]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=[subtype + ' Cluster 1 (' + str(len(data)) + ')'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=[subtype + ' Cluster 3 (' + str(len(data2)) + ')'])
kmf.survival_function_.plot(ax=ax)
plt.title('BIC LumA subtype',fontsize=25)
plt.xlabel('Time in Days',fontsize=30)
plt.ylabel('Survival Rate',fontsize=30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.legend(prop=font,frameon=False)
kmf2 = plt.gcf()
plt.savefig(cancer_type + '_' + subtype + '_survival_plot.eps',format='eps')


# In[45]:

data2


# In[40]:

kmf = KaplanMeierFitter()
ax = plt.subplot(111)
plt.rcParams['font.family'] = 'Arial'

for group in sorted(df['pam_subtype'].unique()):
    g = df.pam_subtype == group
    T = df[g]['days_to_last_followup']
    C = df[g]['event']
    kmf.fit(T, event_observed=C, label= group + ' (' + str(len(T)) + ')')
    kmf.survival_function_.plot(ax=ax,  linewidth=4.0)
kmf2 = plt.gcf()
plt.title("PAM50 labels",fontsize=30)
plt.xlabel('Time in Days',fontsize=30)
plt.ylabel('Survival Rate',fontsize=30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.xlim(xmin=0)
plt.legend(prop=font,frameon=False)
plt.savefig('PAM50_survival_plot_improved.eps',format='eps')





# In[41]:

from lifelines.statistics import multivariate_logrank_test
results = multivariate_logrank_test(np.array(df['days_to_last_followup']),np.array(df['pam_subtype']), np.array(df['event']) )
results.print_summary()


# In[64]:

from lifelines.statistics import logrank_test
results = logrank_test(T, T2, C, C2)
results.print_summary()


# In[41]:

data = df[(df['group'] == '1') & (df['pam_subtype'].str.match('LumA'))]
data


# In[42]:

number_of_patients = len(df)


# In[43]:

horm = np.zeros(number_of_patients)
both_type = np.zeros(number_of_patients)
chemo = np.zeros(number_of_patients)
drug_filename2 = '/home/tunthey/Dropbox/Tunde/analysis/new_datasets/' + cancer_type + '/' + cancer_type +  '_drug_data_v3.txt'
fh = open(drug_filename2)
for ctr,line in enumerate(fh.readlines()):
    if ctr==0: continue # exclude header
    dat = line.split(',')
    if 'chemotherapy' in dat and 'hormone therapy' in dat :
        both_type[ctr] = 1
    if 'chemotherapy' in dat:
        chemo[ctr] = 1
    if 'hormone therapy' in dat:
        horm[ctr] = 1


# In[44]:

df['both_type'] = both_type
df['chemo'] = chemo
df['horm'] = horm


# In[45]:

df


# In[46]:

subtype = 'LumA'
data = df[(df['group'] == '1') & (df['pam_subtype'].str.match(subtype)) & (df['both_type'] == 1)]
data2 = df[(df['group'] == '1') & (df['pam_subtype'].str.match(subtype))  & (df['therapy_type'].str.match('chemotherapy'))]
data3 = df[(df['group'] == '1') & (df['pam_subtype'].str.match(subtype))  & (df['therapy_type'].str.match('hormone therapy'))]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
T3 = data3['days_to_last_followup']
C3 = data3['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=[subtype + ' Cluster 1 Both (' + str(len(data)) + ')'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=[subtype + ' Cluster 1 Chemo only (' + str(len(data2)) + ')'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T3, event_observed=C3, label=[subtype + ' Cluster 1 Hormone only (' + str(len(data3)) + ')'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans base on LumA subtype Cluster 1')
kmf2 = plt.gcf()


# In[47]:

subtype = 'LumA'
data = df[(df['group'] == '3') & (df['pam_subtype'].str.match(subtype)) & (df['both_type'] == 1)]
data2 = df[(df['group'] == '3') & (df['pam_subtype'].str.match(subtype))  & (df['therapy_type'].str.match('chemotherapy'))]
data3 = df[(df['group'] == '3') & (df['pam_subtype'].str.match(subtype))  & (df['therapy_type'].str.match('hormone therapy'))]
T = data['days_to_last_followup']
C = data['event']
T2 = data2['days_to_last_followup']
C2 = data2['event']
T3 = data3['days_to_last_followup']
C3 = data3['event']
ax = plt.subplot(111)
kmf.fit(T, event_observed=C, label=[subtype + ' Cluster 3 Both (' + str(len(data)) + ')'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T2, event_observed=C2, label=[subtype + ' Cluster 3 Chemo only (' + str(len(data2)) + ')'])
kmf.survival_function_.plot(ax=ax)
kmf.fit(T3, event_observed=C3, label=[subtype + ' Cluster 3 Hormone only (' + str(len(data3)) + ')'])
kmf.survival_function_.plot(ax=ax)
plt.title('Lifespans base on LumA subtype Cluster 3')
kmf2 = plt.gcf()


# In[48]:

len(df[(df['group'] == '3') & (df['therapy_type'].str.match('chemotherapy')) & (df['pam_subtype'].str.match(subtype))])


# In[49]:

len(df[(df['group'] == '3') & (df['therapy_type'].str.match('hormone therapy')) & (df['pam_subtype'].str.match(subtype))])


# In[50]:

len(df[(df['group'] == '3') & (df['both_type'] == 1) & (df['pam_subtype'].str.match(subtype))])


# In[51]:

df[(df['group'] == '3') & (df['both_type'] == 1) & (df['pam_subtype'].str.match(subtype)) & (df['vital_status'].str.match('dead')) & (df['time_to_event'] > 1000)]


# In[52]:

df[(df['group'] == '1') & (df['both_type'] == 1) & (df['pam_subtype'].str.match(subtype))]


# In[ ]:



