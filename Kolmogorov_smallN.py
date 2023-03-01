#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:10:42 2022

@author: robertoking
"""

import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import Kolmogorov as kg
from scipy.stats import ks_2samp
import time

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

def P(List,threshold=1.36):
    Listabove = [i for i in List if i > threshold]
    p = len(Listabove)/len(List)
    return p

#%%
start_time = time.time()
# N = 100
Sims = 2000
Nvalues = [100,200,200,200]
scaleMC = 10

p_vals1, p_vals2, p_vals3, p_vals4 = [], [], [], []
for j in Nvalues:
    N=j
    dcur1, dcur2, dcur3, dcur4 = [], [], [], []
    Dpoints1 = random.sample(D1,N)
    Dpoints2 = random.sample(D2,N)
    Dpoints3 = random.sample(D3,N)
    Dpoints4 = random.sample(D4,N)
    for i in range(Sims):
        print(j,i+1,'/',Sims)
        Mpoints1 = random.sample(M1,N*scaleMC)
        Mpoints2 = random.sample(M2,N*scaleMC)
        Mpoints3 = random.sample(M3,N*scaleMC)
        Mpoints4 = random.sample(M4,N*scaleMC)
        
        # Dpoints = np.array(D1[N*i:N*(i+1)])
        # Mpoints = np.array(M1[N*i:N*(i+1)])
        
        # dcur.append(kg.Kolm(Dpoints,Mpoints))
        dcur1.append(ks_2samp(Dpoints1,Mpoints1)[1])
        dcur2.append(ks_2samp(Dpoints2,Mpoints2)[1])
        dcur3.append(ks_2samp(Dpoints3,Mpoints3)[1])
        dcur4.append(ks_2samp(Dpoints4,Mpoints4)[1])
    p_vals1.append(dcur1)
    p_vals2.append(dcur2)
    p_vals3.append(dcur3)
    p_vals4.append(dcur4)

print('Time taken',time.time()-start_time,'s')
#%%
nobins = 20
maxheight = 550
dens = True
linwid = 1.5

def plot4():
    fig = plt.figure(figsize=(10,10))
    gs = plt.GridSpec(2,2)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])
    
    ax1.hist(p_vals1,bins = nobins,density = dens,histtype = 'step',lw=linwid)
    ax2.hist(p_vals2,bins = nobins,density = dens,histtype = 'step',lw=linwid,label=[i for i in Nvalues])
    ax3.hist(p_vals3,bins = nobins,density = dens,histtype = 'step',lw=linwid)
    ax4.hist(p_vals4,bins = nobins,density = dens,histtype = 'step',lw=linwid)
    ax2.legend(fontsize=12)
    
    ax1.set_xlim(0,1)
    ax2.set_xlim(0,1)
    ax3.set_xlim(0,1)
    ax4.set_xlim(0,1)
    
    ax1.set_ylabel('Freq Density')
    ax2.set_ylabel('Freq Density')
    ax3.set_ylabel('Freq Density')
    ax4.set_ylabel('Freq Density')
    
    ax1.set_title('$E_{\u03BC}$')
    ax2.set_title('$cos \u03B8_{\u03BC}$')
    ax3.set_title('$E_{p}$')
    ax4.set_title('$cos \u03B8_{p}$')
    
    fig.suptitle('P-distributions for different N, for different dimensions',weight='bold')
    # fig.legend()
    
plot4()
# for i in range(len(d_vals)):
#     rootn = np.sqrt(Nvalues[i])
#     for j in range(len(d_vals[i])):
#         d_vals[i][j] = d_vals[i][j]/rootn

#%%
d_vals = []
for j in Nvalues:
    N=j
    dcur = []
    Dpoints = random.sample(D1,N)
    for i in range(Sims):
        print(j,i+1,'/',Sims)
        Mpoints = random.sample(M1,N*scaleMC)
        
        # Dpoints = np.array(D1[N*i:N*(i+1)])
        # Mpoints = np.array(M1[N*i:N*(i+1)])
        
        # dcur.append(kg.Kolm(Dpoints,Mpoints))
        dcur.append(ks_2samp(Dpoints,Mpoints)[1])
    d_vals.append(dcur)

def plot5():
    plt.hist(d_vals[0],bins = nobins,density = dens,histtype = 'step',lw=linwid,label='N = %s' %( Nvalues[0]))
    # plt.hist(d_vals[1],bins = nobins,histtype = 'step',label='N = 500')
    plt.hist(d_vals[1],bins = nobins,density = dens,histtype = 'step',lw=linwid,label='N = %s' %( Nvalues[1]))
    plt.hist(d_vals[2],bins = nobins,density = dens,histtype = 'step',lw=linwid,label='N = %s' %( Nvalues[2]))
    plt.hist(d_vals[3],bins = nobins,density = dens,histtype = 'step',lw=linwid,label='N = %s' %( Nvalues[3]))
    
    # plt.plot([1.358,1.358],[0,maxheight],color='black',lw = 1,ls='-',label='5% significance line')
    # plt.plot([1.628,1.628],[0,maxheight],color='black',lw = 1,ls='dashed',label='1% significance line')
    # plt.plot([1.949,1.949],[0,maxheight],color='black',lw = 1,ls=':',label='0.1% significance line')
    
    plt.title('Proton Angle: P-value histogram for different N, Simulations = %s' %(Sims))
    plt.xlabel('P-value')
    plt.ylabel('Frequency Density')
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.legend()

plot5()

 #%%
plt.title('d value histogram for N = %s , Simulations = %s' %(N,Sims))
plt.hist(d_vals,bins = 100)
#%%


print(kg.Kolm(D1[:N],M1[:N]))