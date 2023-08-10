#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:41:03 2023

@author: kasturilele
change from 4.4
 - changing the appender functions so that they can pick a random replicate and not a fixed replicate for each set of parameters
 - also, make models with fewer species in each model - leave one out
 - start with only the five species we know coexist in the community
 
change from 5 - 
 - changing the pairwise appender function so that it picks randomly between replicates as well
          
"""

import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import csv
import random


single_data = pd.read_csv("single_parameters.csv")
single_data = single_data.to_dict('records')

paired_data = pd.read_csv("paired_parameters_1.csv")
paired_data = paired_data.to_dict('records')

#function that converts a list to a tuple - for the output of the next function
def convert(list):
    return (*list, )

#final function - constructs the equations as it goes
def testFun5(t, y, r, a_s, a_p, N):
    outputs = []
    for j in range(0,N):
        pairs = []
        for k in range(0,N):
            if j == k:
                pairs.append(0)
            else:
                pairs.append(a_p[j][k]*y[k])
        outputs_calc = y[j]*(r[j] + a_s[j]*y[j] + sum(pairs))
        outputs.append(outputs_calc)
    outputs_tuple = convert(outputs)
    return(outputs_tuple)
              
t_int = np.linspace(0,48,100) #defines a vector of timepoints to solve over

strain_names = ['F.sanfranciscensis', 'L.brevis', 'L.plantarum', 'C.paralimentarius','A.malorum',
                'S.cerevisiae','W.anomalus','K.humilis','K.servazzii']
all_strains = strain_names

#all_strains = ['F.sanfranciscensis', 'L.brevis', 'C.paralimentarius','A.malorum',
#               'S.cerevisiae','W.anomalus','K.humilis','K.servazzii']


strain_color_main = ['#d9d9d9','#882255','#ddcc77','#aa4499','#cc6677','#117733','#88CCEE','#332288','#44aa99']
#strain_color = ['#d9d9d9','#882255', '#aa4499','#cc6677','#117733','#88CCEE','#332288','#44aa99' ]
strain_color = strain_color_main

#making empty lists to save the model results and replicates picked for each run of the model for later
y_endpoint = []
repsPickedS = []
repsPickedP = []

#initialize single parameter dictionary
single_dict = {x: [] for x in strain_names}

#populate single parameter dictionary
for i in range(0,len(single_data)):
     tempStrain = single_data[i]
     tempKey = tempStrain['Strains']
     tempPars = [tempStrain['N0'],tempStrain['r'],tempStrain['a11']]
     single_dict[tempKey].append(tempPars)

#initialize paired parameter dictionary
paired_dict = {x: {} for x in strain_names}

for i in range(0,len(strain_names)):
    outerVar = strain_names[i]
    delStrain = [x for x in strain_names if x != outerVar]
    paired_dict[outerVar] = {x: [] for x in delStrain}

#populate paired parameter dictionary
for i in range(0,len(paired_data)):
    tempStrain = paired_data[i]
    tempKey1 = tempStrain['Strain_1']
    tempKey2 = tempStrain['Strain_2']
    paired_dict[tempKey1][tempKey2].append(tempStrain['a12'])
    paired_dict[tempKey2][tempKey1].append(tempStrain['a21'])
    

#define a function to make the input list per replicate for single parameters
def appender(strains,array,iList,j,n):
    for a in range(0,n):
        S = strains[a]
        iTemp = iList[a]
        array.append(single_dict[S][iTemp][j])
    return array


# define a similar function for paired parameters
def appenderPairsNew(strains,iList,n):
    pairMatrix = [[0] * n for _ in range(n)]
    idx = 0
    for j in range(0,n):
        for k in range(j+1,n):
            S1 = strains[j]
            S2 = strains[k]
            i = iList[idx]
            pairMatrix[j][k] = paired_dict[S1][S2][i]
            pairMatrix[k][j] = paired_dict[S2][S1][i]
            idx += 1
    return pairMatrix

#function that counts how many replicates are there in a given pair of strains
def varCounter(strains,x,y):
    vS1 = strains[x]
    vS2 = strains[y]
    varLen = len(paired_dict[vS1][vS2])
    return varLen


# 9-species equation solver - with more randomized sampling from parameter distributions
seed_start = 1998

for olN in range(1000):
    
    #the fixed parameters are set here for each run
    num_species = 9
    #N_init = appender(all_strains, [], i, 0, num_species)
    N_init = [1600,3,199,272,171,530,34,390,169]
    
    random.seed(seed_start)
    testList1 = []
    testList2 = []
    
    for lN in range (num_species):
        testList1.append(random.randint(0,5))
    
    for xlN in range(0,num_species):
        for ylN in range(xlN+1,num_species):
            uB = varCounter(all_strains, xlN, ylN) - 1
            testList2.append(random.randint(0,uB))

    #set the parameters that vary - these are now randomized for both single species as well as pairwise parameters
    g_rate = appender(all_strains, [], testList1, 1, num_species)
    alpha_single = appender(all_strains, [], testList1, 2, num_species)
    
    alpha_pair = appenderPairsNew(all_strains, testList2, num_species)
    sol = solve_ivp(testFun5, [0, 48], N_init, args = (g_rate, alpha_single, alpha_pair, num_species),dense_output=True)

    y_inter = sol.sol(t_int) #stores values for all timepoints being solved over
    y_endpoint.append(y_inter[:,99])
    repsPickedS.append(testList1)
    repsPickedP.append(testList2)
    
    # for pltnum in range(0,num_species):
    #     plt.plot(t_int, y_inter[pltnum], color = strain_color[pltnum], linewidth = 1, label = all_strains[pltnum])
    
    # plt.title('logistic growth model (with interaction)')
    # plt.legend(loc = "upper left")
    # plt.show()
        
    seed_start = seed_start + 1 
    print(seed_start)

  
with open('Outnew_1000.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')
      
    write.writerow(all_strains)
    write.writerows(y_endpoint)

with open('Outnew_1000_reps_single.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')
      
    write.writerow(all_strains)
    write.writerows(repsPickedS)
    
with open('Outnew_1000_reps_paired.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')
    
    write.writerows(repsPickedP)





