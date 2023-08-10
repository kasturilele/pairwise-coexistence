#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:12:20 2023

@author: kasturilele

further testing of model 5 - taking the parameters from rejection sampling and seeing if that changes anything about the model predictions
"""

import pandas as pd
import random

param_subset = pd.read_csv("parameters_subset.csv")
param_subset = param_subset.to_dict('records')

strain_names = ['F.sanfranciscensis', 'L.brevis', 'L.plantarum', 'C.paralimentarius','A.malorum',
                'S.cerevisiae','W.anomalus','K.humilis','K.servazzii']

#initialize dictionary for parameter subset replicates
param_subset_dict = {x: [] for x in strain_names}

#populate parameter subset dictionary
for i in range(0,len(param_subset)):
     tempDict = param_subset[i]
     for j in range(0, len(strain_names)): 
         tempRep = tempDict[strain_names[j]]
         param_subset_dict[strain_names[j]].append(tempRep)
         

seed_start = 1998
for olN in range(10):
    random.seed(seed_start)
    testList1 = []
    for lN in range (9):
        param_subList = param_subset_dict[strain_names[lN]]
        testList1.append(random.choice(param_subList))
    print(testList1)
    seed_start = seed_start + 1 