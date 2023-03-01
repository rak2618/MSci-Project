#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 17:51:07 2022

@author: robertoking
"""

import NeutNo3 as NN
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import MyFunctions as mf

#%%
Nd = 50  #Data Sample
Nm = 50 #Model Sample
NoOfTSims = 100
NoOfPSims = 1

D = pd.read_csv('gof-master/neut_events_4D.csv')
D = np.array(D)

M = pd.read_csv('gof-master/nuwro_events_4D.csv')
M = np.array(M)

D1 = [i[0] for i in D]
D2 = [i[1] for i in D]
D3 = [i[2] for i in D]
D4 = [i[3] for i in D]

M1 = [i[0] for i in M]
M2 = [i[1] for i in M]
M3 = [i[2] for i in M]
M4 = [i[3] for i in M]

DSamp = D[0:Nd]
MSamp = M[0:Nm]

DSamp2 = D[Nd:2*Nd]

sigM1 = np.std(M1)
sigM2 = np.std(M2)
sigM3 = np.std(M3)
sigM4 = np.std(M4)

SigM1 = mf.FWHM(M1)
SigM2 = mf.FWHM(M2)
SigM3 = mf.FWHM(M3)
SigM4 = mf.FWHM(M4)


globalSigma = np.mean([SigM1,SigM2,SigM3,SigM4])
globalSigm  = (SigM1+SigM2+SigM3+SigM4)/4

sigD1 = np.std(D1)
sigD2 = np.std(D2)
sigD3 = np.std(D3)
sigD4 = np.std(D4)

SigD1 = mf.FWHM(D1)
SigD2 = mf.FWHM(D2)
SigD3 = mf.FWHM(D3)
SigD4 = mf.FWHM(D4)

def VectGaussianTerm1(vector):
    out = np.exp(-0.5*(vector[0]**2/sigD1**2 + vector[1]**2/sigD2**2 + vector[2]**2/sigD3**2 + vector[3]**2/sigD4**2))
    return out

def VectGaussianTerm2(vector):
    out = np.exp(-0.5*(vector[0]**2/sigM1**2 + vector[1]**2/sigM2**2 + vector[2]**2/sigM3**2 + vector[3]**2/sigM4**2))
    return out

def VectGaussianTerm3(vector):
    out = np.exp(-0.5*(vector[0]**2/sigD1*sigM1 + vector[1]**2/sigD2*sigM2 + vector[2]**2/sigD3*sigM3 + vector[3]**2/sigD4*sigM4))
    return out



def vectgauss1(vector):
    out = np.exp(-0.5*(vector[0]**2/SigD1**2 + vector[1]**2/SigD2**2 + vector[2]**2/SigD3**2 + vector[3]**2/SigD4**2))
    return out

def vectgauss2(vector):
    out = np.exp(-0.5*(vector[0]**2/SigM1**2 + vector[1]**2/SigM2**2 + vector[2]**2/SigM3**2 + vector[3]**2/SigM4**2))
    return out

def vectgauss3(vector):
    out = np.exp(-0.5*(vector[0]**2/SigD1*SigM1 + vector[1]**2/SigD2*SigM2 + vector[2]**2/SigD3*SigM3 + vector[3]**2/SigD4*SigM4))
    return out

#%%
T_Vals_disp = []
for i in range(NoOfTSims):
    print(i)
    np.random.shuffle(M)
    MNew = M[0:Nm]
    # T_Vals_disp.append(NN.T_FUnction(DSamp,MNew,NN.VectGaussianTerm2,NN.VectGaussianTerm2,NN.VectGaussianTerm2))
    T_Vals_disp.append(NN.T_FUnction(DSamp,MNew,vectgauss1,vectgauss2,vectgauss3))


plt.hist(T_Vals_disp)
print('sigma = ',np.std(T_Vals_disp))

#%%
T_Vals_dist = []

def Gaussian(x):
    out = np.exp((-x**2)/2*globalSigma**2)
    return out

for i in range(NoOfTSims):
    print(i)
    np.random.shuffle(M)
    MNew = M[0:Nm]
    T_Vals_dist.append(mf.T_Function_4d(DSamp,MNew,Psi = Gaussian))

plt.hist(T_Vals_dist)
print('sigma = ',np.std(T_Vals_dist))
    




