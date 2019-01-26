
## A script for finding every cox coefficient and pvalue for every mRNA in GBM Tier 3 data downloaded Jan. 5th, 2016


## Load necessary modules
# from rpy2 import robjects as ro
import numpy as np
import os
# ro.r('library(survival)')
import re

##This call will only work if you are running python from the command line.
##If you are not running from the command line manually type in your paths.
# BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

def load_clinical():
    f=open('../data/clinical_data/nationwidechildrens.org_clinical_follow_up_v1.0_gbm.txt')
    ##get the column indexes needed
    columns=f.readline().split('\t')
    patient_column=columns.index('bcr_patient_barcode')
    alive_column=columns.index('last_contact_days_to')
    death_column=columns.index('death_days_to')
    f.readline()
    f.readline()
    data=[i.split('\t') for i in f]
    ## A patient can be listed multiple times in the file. The most recent listing (furthest down in the file), contains the most recent
    ## follow up data.  This code checks if the patient has already been loaded into the list, and if so, takes the more recent data.
    ## This required an empty value in the list initialization.
    ## Data is: [[Patient ID, time(days), Vital status],[Patient ID, time(days), Vital status],...]
    clinical1=[['','','']]
    for i in data:
        if clinical1[-1][0]==i[patient_column]:
            if re.search('^[0-9]+$',i[death_column]):
                clinical1[-1]=[i[patient_column],int(i[death_column]),'Dead']
            elif re.search('^[0-9]+$',i[alive_column]):
                clinical1[-1]=[i[patient_column],int(i[alive_column]),'Alive']
            else:
                pass
        else:
            if re.search('^[0-9]+$',i[death_column]):
                clinical1.append([i[patient_column],int(i[death_column]),'Dead'])
            elif re.search('^[0-9]+$',i[alive_column]):
                clinical1.append([i[patient_column],int(i[alive_column]),'Alive'])
            else:
                pass

    ## Removing the empty value.
    clinical=clinical1[1:]


    ## Sex and age information were taken from the "clinical_patient" file.  A dictionary was created for sex.
    more_clinical={}


    sex_dict={}
    sex_dict['MALE']=0
    sex_dict['FEMALE']=1



    ## The "clinical_patient" file can also contain patients not listed in the follow_up files.
    ## In these cases the clinical data for these patients gets appended to a new clinical list.

    f=open('../data/clinical_data/nationwidechildrens.org_clinical_patient_gbm.txt')
    ##get the column indexes needed
    columns=f.readline().split('\t')
    sex_column=columns.index('gender')
    age_column=columns.index('age_at_initial_pathologic_diagnosis')
    patient_column=columns.index('bcr_patient_barcode')
    alive_column=columns.index('last_contact_days_to')
    death_column=columns.index('death_days_to')
    f.readline()
    f.readline()
    clinical4=[]
    data=[i.split('\t') for i in f]
    for i in data:
        try:
            more_clinical[i[patient_column]]=[0,sex_dict[i[sex_column]],int(i[age_column])]
            if re.search('^[0-9]+$',i[death_column]):
                clinical4.append([i[patient_column],int(i[death_column]),'Dead'])
            elif re.search('^[0-9]+$',i[alive_column]):
                clinical4.append([i[patient_column],int(i[alive_column]),'Alive'])
            else:
                pass
        except:
            pass

    new_clinical=[]


    ##It is possible that the clinical data in the clinical_patient file is more up to date than the follow_up files
    ##All the clinical data is merged checking which data is the most up to date
    for i in clinical4:
        if i[0] not in [j[0] for j in clinical]:
            new_clinical.append(i)
        else:
            if i[1]<=clinical[[j[0] for j in clinical].index(i[0])][1]:
                new_clinical.append(clinical[[j[0] for j in clinical].index(i[0])])
            else:
                new_clinical.append(i)


    ##also do the reverse since clinical can contain patients not included in clinical4
    for i in clinical:
        if i[0] not in [j[0] for j in new_clinical]:
            new_clinical.append(i)




    ## only patients who had a follow up time greater than 0 days are included in the analysis
    clinical=[i for i in new_clinical if i[1]>0]



    final_clinical=[]

    ## A new list containing both follow up times and sex and age is constructed.
    ## Only patients with sex and age information are included.
    ## Data is [[Patient ID, time (days), vital status, 0, sex, age at diagnosis],...]

    for i in clinical:
        if i[0] in more_clinical:
            final_clinical.append(i+more_clinical[i[0]])
    return final_clinical
