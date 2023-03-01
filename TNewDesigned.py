#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 18:20:46 2022

@author: robertoking

New T-Function designed from scratch

Takes in two lists of 4d points
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import itertools as it
import time


def FindDisp(PointA,PointB):
    return abs(PointB - PointA)

def Term1(Points,Psi):
    N = len(Points)
    output = 0
    
    for i in range(0,N):
        for j in range(i+1,N):
            output += Psi(abs(Points[i] - Points[j]))
    
    return (1/(N*(N-1)))*output

def VTerm1(Points,Psi):
    N = len(Points)
    Points = np.array(Points)
    
    Displacements = np.array([pair[1]-pair[0] for pair in it.combinations(Points,r=2)])
    # print('d')
    # print(Displacements)
    
    
    PsiVals = [Psi(Disp) for Disp in Displacements]
    # print(PsiVals)
    
    return (1/(N*(N-1)))*np.sum(PsiVals)
    
    
    
    
def Term3(DataPoints,MCPoints,Psi):
    Nd , Nmc = len(DataPoints) , len(MCPoints)
    output = 0
    
    for i in range(Nd):
        for j in range(Nmc):
            output += Psi(abs(DataPoints[i] - MCPoints[j]))
    return (1/(Nd * Nmc))*output


def T_func(DataPoints,MCPoints,Psi):
    t1 = Term1(DataPoints,Psi)
    t2 = Term1(MCPoints,Psi)
    t3 = Term3(DataPoints,MCPoints,Psi)
    
    return t1 + t2 + t3

def VT_func(DataPoints,MCPoints,Psi):
    t1 = VTerm1(DataPoints,Psi)
    t2 = VTerm1(MCPoints,Psi)
    t3 = Term3(DataPoints,MCPoints,Psi)
    
    return t1 + t2 + t3

def P(List,threshold=0.05):
    Listbelow = [i for i in List if i < threshold]
    p = len(Listbelow)/len(List)
    return p

#%%
D = pd.read_csv('gof-master/nuwro_events_4D.csv')
# D = pd.read_csv('gof-master/neut_events_4D.csv')
D = np.array(D)

M = pd.read_csv('gof-master/nuwro_events_4D.csv')
# M = pd.read_csv('gof-master/neut_events_4D.csv')
M = np.array(M)

D1 = [i[0] for i in D]
D2 = [i[1] for i in D]
D3 = [i[2] for i in D]
D4 = [i[3] for i in D]

M1 = [i[0] for i in M]
M2 = [i[1] for i in M]
M3 = [i[2] for i in M]
M4 = [i[3] for i in M]

sigma1 = mf.FWHM(M1)
sigma2 = mf.FWHM(M2)
sigma3 = mf.FWHM(M3)
sigma4 = mf.FWHM(M4)

def VectorGaussian(Vector):
    power = -0.5*((Vector[0]/sigma1)**2 + (Vector[1]/sigma2)**2 + (Vector[2]/sigma3)**2 + (Vector[3]/sigma4)**2)
    return np.exp(power)


#%%
Ld = 200
Lmc = Ld

DSamp = D[0:Ld]
MSamp = M[0:Lmc]

start_time = time.time()

myval = T_func(DSamp,MSamp,VectorGaussian)

print('T =',myval)
print('Time Elapsed =',time.time()-start_time,'s')

#%%
start_time = time.time()

NoofT = 100
NoofP = 100

MSize = NoofP * Lmc
Marray = M[0:MSize]
np.random.shuffle(Marray)

DSample = D[0:Ld] #M[0:Lmc]

P_values = []
#%%

for i in range(NoofP):
    
    MSample = Marray[Lmc*i:Lmc*(i+1)]
    BigList = np.concatenate((MSample,DSample))
    
    T0 = T_func(DSample,MSample,VectorGaussian)
    T_Perm = []
    for j in range(NoofT):
        np.random.shuffle(BigList)
        MNew = BigList[0:Lmc]
        DNew = BigList[Lmc:Lmc+Ld]
        T_Perm.append(T_func(DNew,MNew,VectorGaussian))
        print("Currently of p value", i, "on TPerm number",j)
    
    T_above = [t for t in T_Perm if t > T0]
    p = len(T_above)/len(T_Perm)
    P_values.append(p)
        
print('Time Elapsed =',time.time()-start_time,'s')
    


#%%
fig = plt.figure(figsize=(5,5))
gs = plt.GridSpec(1,1)
ax = plt.subplot(gs[0])

ax.hist(P_values,bins=10,range = [0,1],density = True,histtype = 'step',lw=2)
ax.set_xlim([0,1])
ax.set_xlabel('P-Value')
ax.set_ylabel('Frequency Density')
fig.suptitle('$N_{d} = %s,N_{m} = %s$' %(Ld,Lmc))

df = pd.DataFrame(P_values)
df.to_csv('P_Vals_100_100_300_100')






















#%%

y = np.array([4,3,2])
x = np.array([1,7,9])


#%%

Points = np.random.randint(1,6,size=(200,4))
mPoints = np.random.randint(11,16,size=(200,4))
# Points = np.zeros(shape=(4,4))
Points[-1]=[1,0,0,0]
# print(Points)
N = len(Points)
Points = np.array(Points)

t1 = time.time()
# VTerm1(Points,VectorGaussian)
# Displacements = np.zeros(shape=(int(len(Points)*(len(Points)-1)/2),4))
# Displacements[:] = [Points[i+1:]-Points[i] for i in range(N)]
# Displacements = np.array([Points[i+1:]-Points[i] for i in range(N)])
# Displacements = np.array([pair[1]-pair[0] for pair in it.combinations(Points,r=2)])
# Displacements.flatten()

# for i in range(0,N):
#     Disp = Points[i+1:] - Points[i]
#     Displacements = np.concatenate((Displacements,Disp))
    # print(i,Disp)
    

# PsiVals = np.zeros(shape = (len(Displacements)))
# PsiVals[:] = Psi(Displacements[:])
T = VT_func(Points,mPoints,VectorGaussian)
print('T',T)
print(time.time()-t1,'s')













