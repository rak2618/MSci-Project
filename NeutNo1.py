#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 19:48:13 2022

@author: robertoking
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import time

#%%
Nd = 50  #Data Sample
Nm = 500 #Model Sample
NoOfSims = 200

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

D34 = [[i[2],i[3]] for i in D]
M34 = [[i[2],i[3]] for i in M]
#%%
mf.histogram2d(D,'D',2,D,'D',3,NBins = 200,mycolours = 'plasma')
mf.scatterplot(M,'M',1,M,'M',2,N=100)
#%%
Data = D[0:50]
MCPoints = M[0:500]

T = mf.T_Function_4d(Data,MCPoints,Psi = mf.D)
print(T)




#%%

T1 = mf.T_Function_2d(D34[:50],M34[:50],Psi = mf.DSquared)
print(T1)

#%%
start = time.time()
T_Vals = []
for i in range(0,NoOfSims):
    T_Vals.append(mf.T_Function_4d(D[(i)*50:(i+1)*50],M[(i)*50:(i+1)*50],Psi = mf.DRecip))


#%%
plt.hist(T_Vals,bins='auto',histtype='step',density = True ,label='$\u03C8(x) = 1/x $')   #e^{-x^{2}/2\u03C3^{2}}
plt.xlabel('T-Value')
plt.ylabel('Frequency Density')
plt.title('T-value histogram for Neut vs Nuwro (4 dimensions)')
plt.legend()

#%%

plt.title('$cos(\u03B8_{\u03BC})$ histogram',fontweight='light')
plt.hist(D2,200,color='black')
plt.xlim(0,1)
plt.xlabel('$cos(\u03B8_{\u03BC})$')
plt.ylabel('Frequency, f',fontsize=15)
# num = '$f_{max}$'
# denom = '2'
# label = '$\frac{%s}{%s}$' %(num,denom)
label = r'$\frac{f_{max}}{2}$'
plt.yticks([10408,20815],labels=[label,'$f_{max}$'],fontsize=20,fontweight='bold')



















