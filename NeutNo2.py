#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 10:56:22 2022

@author: robertoking

In this file I will perfect the permutation test and try using classes to organise my work
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import time

def Gaussian(x):
    out = np.exp((-D**2)/2)
    return out
#%%
Nd = 100  #Data Sample
Nm = 100 #Model Sample
Nperm = 5
NT = 5

D = pd.read_csv('gof-master/neut_events_4D.csv')
D = np.array(D)

M = pd.read_csv('gof-master/nuwro_events_4D.csv')
M = np.array(M)

DSamp = D[0:Nd]
MSamp = M[0:Nm]

#%%
T0 = mf.T_Function_4d(DSamp,MSamp,Psi=mf.DRecip)

BigList = np.concatenate((DSamp,MSamp))
T_Perm = []
for i in range(NT):
    np.random.shuffle(BigList)
    MNew = BigList[0:Nm]
    DNew = BigList[Nm:Nm+Nd]
    T_Perm.append(mf.T_Function_4d(DNew,MNew,Psi = mf.DRecip))
    print("Finished TPerm number",i)

#%%
T_Below = [i for i in T_Perm if i < T0]
p = len(T_Below)/len(T_Perm)
print(p)

#%%
start = time.time()
P_List = []
#%%

BigList = np.concatenate((DSamp,MSamp))
for j in range(Nperm):
    MSamp = M[Nm*i:Nm*(i+1)]
    T0 = mf.T_Function_4d(DSamp,MSamp,Psi=mf.DRecip)
    
    T_Perm = []
    for i in range(NT):
        np.random.shuffle(BigList)
        MNew = BigList[0:Nm]
        DNew = BigList[Nm:Nm+Nd]
        T_Perm.append(mf.T_Function_4d(DNew,MNew,Psi = mf.DRecip))
        print("Currently of p value", j, "on TPerm number",i)
    
    T_Above = [i for i in T_Perm if i > T0]
    p = len(T_Above)/len(T_Perm)
    P_List.append(p)
    
print("Time Elapsed =",time.time()-start,"s")
#%%
# P1_List = [i for i in P_List]

plt.hist(P_List,bins=10,range = [0,1],density = False,histtype = 'step')
plt.xlabel('P-Value')
plt.ylabel('Frequency')
plt.title('Permutation Test P-Values, Hist Events = 200 \n Sims per PVal = 200, Nd = 50,Nm = 50')














    