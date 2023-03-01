#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:38:24 2022

@author: robertoking
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf

from scipy.spatial import distance as dt
import time
import random


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


def FNN(k,point,List):
    '''

    Parameters
    ----------
    k : TYPE
        DESCRIPTION.
    point : TYPE
        DESCRIPTION.
    List : list of vectors
        DESCRIPTION.

    Returns
    -------
    K Indexes of points in List that are closest to point

    '''
    List = np.array(List)
    point = np.array(point)
    g = List - point
    dists = np.linalg.norm(g,axis=1)
    sortdists = np.sort(dists)
    distsneeded = sortdists[1:k+1]
    inds = [np.where(dists == val)[0][0] for val in distsneeded]
    
    return inds

def MixedT(List1,List2,k):
    T = 0
    LongList = np.concatenate((List1,List2))
    N1 = len(List1)
    N2 = len(List2)
    for i in range(len(List1)):
        inds = FNN(k,List1[i],LongList)
        for l in inds:
            if l + 1 <= N1:
                T += 1
    
    for j in range(len(List2)):
        inds = FNN(k,List2[j],LongList)
        for l in inds:
            if l + 1 > N1:
                T += 1
                
    return (1/(k*(N1+N2)))*T
        
    
def PullT(ListA, ListB,K):
    Tval = MixedT(ListA,ListB,K)
    na = len(ListA)
    nb = len(ListB)
    n = na + nb
    
    muT = (na*(na-1) + nb*(nb-1))/(n*(n-1))
    sigTsquared = (1/(n*K))*(na*nb/(n**2) + 4*na**2 *nb**2/(n**4))
    sigT = np.sqrt(sigTsquared)
    
    ullT = (Tval - muT ) / sigT
    return ullT

def P(List,threshold=1.64):
    Listabove = [i for i in List if i > threshold]
    p = len(Listabove)/len(List)
    return p

#%%
Nd  = 1000          #Data sample size
Nmc = 1000          #MCarlo sample size
myK = 10        #No of nearest neighbours checked

MSamp = M[0:Nmc]
DSamp = D[0:Nd ]

start_time = time.time()

# print('Tvalue =', MixedT(MSamp,DSamp,myK))
print('PullT =', PullT(MSamp,DSamp,myK))

print('Time Elapsed',time.time()-start_time,'s')
            



#%%
'Creating a basic T histogram'
start_time = time.time()
T_Vals = []

Nd  = 100          #Data sample size
Nmc = 1000          #MCarlo sample size
myK = 10        #No of nearest neighbours checked
NoofTs = 2000

n = Nd + Nmc
muT = (Nd*(Nd-1) + Nmc*(Nmc-1))/(n*(n-1))
sigTsquared = (1/(n*myK))*(Nd*Nmc/(n**2) + 4*Nd**2 *Nmc**2/(n**4))
sigT = np.sqrt(sigTsquared)

random.seed(1)
Dind = random.sample(range(len(D)), Nd)
DSample = np.array([D[i] for i in Dind])

random.seed(26)
print('loop time')
for i in range(NoofTs):
    # start_time1 = time.time()
    # print(i/NoofTs * 100,'%')
    print(i)
    'Create Monte-Carlo Random Sample'
    Mind = random.sample(range(len(M)),Nmc)
    MSample = np.array([M[i] for i in Mind])
    # print(time.time()-start_time1,'s')
    
    T_Vals.append(MixedT(DSample,MSample,myK))
    # T = MixedT(DSample,MSample,myK)
    # T_Vals.append(PullT(DSample,MSample,myK))
    

print('loop done')
T_Vals = np.array(T_Vals)
T_Vals = (T_Vals - muT )/sigT

print('Time Elapsed',time.time()-start_time,'s')
    

#%%
Tabove = [i for i in T_Vals if i > 1.64]
rejectionpower = len(Tabove)/len(T_Vals)
print('Rejection Power =',rejectionpower)

plt.hist(T_Vals,bins='auto',density = False,histtype = 'step',lw=2)
plt.xlabel('Pull-T')
plt.ylabel('Frequency')
plt.title('Mixed Sampled Method T-Distribution, Neut v Nuwro,\n Hist Events = %s,$ N_{d} = %s,N_{mc} = %s, N_{k} = %s $' %(NoofTs,Nd,Nmc,myK))
# plt.savefig('MixedSamp1')



#%%
'Testing the effect of K on the T value'
noTs = 500
Nd = 100
Nmc = 10*Nd

start_time = time.time()

T_Vals = []
for Kval in np.arange(5,30,5):
    Dind = random.sample(range(len(D)), Nd)
    DSample = np.array([D[i] for i in Dind])
    
    PermTvals = []
    for i in range(noTs):
        'Create Random Sample'
        Mind = random.sample(range(len(M)),Nmc)
        MSample = np.array([M[i] for i in Mind])
        
        # PermTvals.append(MixedT(DSample,MSample,Kval))
        PermTvals.append(PullT(DSample,MSample,Kval))
        
        print('Kval =',Kval,'Tval =',i)
    
    T_Vals.append(PermTvals)
    
print('Time Elapsed',time.time()-start_time,'s')


#%%
plt.hist(T_Vals,bins='auto',density=True,histtype = 'step',lw=1.5,label=[i for i in np.arange(5,30,5)])
plt.plot([1.64,1.64],[0,0.7],c='black',label =  '5% significance line')
plt.plot([1.28,1.28],[0,0.7],':',c='black',label = '10% significance line')
plt.legend()
plt.title('T Histograms for different $n_{k}$ ,%s events per hist, $n_{d}$ = %s, $n_{mc}$ = %s' %(noTs,Nd,Nmc))
plt.xlabel('Pull-T')
plt.ylabel('Frequency')
plt.savefig('PullTHistograms',format='pdf')

#%%
T_Avg = []
for i in range(len(T_Vals)):
    T_Avg.append(np.mean(T_Vals[i]))
    
plt.plot(np.arange(5,Nd+Nmc,5)/(Nd+Nmc),T_Avg,'o')
plt.xlabel('$ n_{k}/n$')
plt.ylabel('T value (averaged over ten)')
plt.title('T values for different $n_{k}$')


#%%
'Investigating the effect of Nmc on the T histogram and p value'

myK = 10
noTs = 500
Nd = 100
NMCVals = [200,500,1000,2000]

start_time = time.time()
T_Vals = []
for N_MC in NMCVals:
    Dind = random.sample(range(len(D)),  Nd)
    DSample = np.array([D[i] for i in Dind])
    
    tempTVals = []
    for i in range(noTs):
        'Create Random Sample'
        Mind = random.sample(range(len(M)),N_MC)
        MSample = np.array([M[i] for i in Mind])
        
        tempTVals.append(PullT(DSample,MSample,myK))
        
        print('Testing Nmc =',N_MC,'Tval =',i+1)
    
    T_Vals.append(tempTVals)

timetaken = mf.minutes(time.time()-start_time)
print('Time Elapsed',timetaken[0],'mins',timetaken[1],'s')
#%%
'Investigating the effect of Nd and Nmc on the T-Hist and p value'
myK = 10
noTs = 500
Ndvals = [100,200,500,1000]
Nmcvals = Ndvals
Nvals = Ndvals

start_time = time.time()
T_Vals = []
timevals = []
for N in Nvals:
    t1 = time.time()
    Dind = random.sample(range(len(D)),  N)
    DSample = np.array([D[i] for i in Dind])
    
    tempTVals = []
    for i in range(noTs):
        'Create Random Sample'
        Mind = random.sample(range(len(M)),N)
        MSample = np.array([M[i] for i in Mind])
        
        tempTVals.append(PullT(DSample,MSample,myK))
        
        print('Testing N =',N,'Tval =',i+1)
    
    T_Vals.append(tempTVals)
    print(time.time()-t1)
    timevals.append(time.time()-t1)
    
timetaken = mf.minutes(time.time()-start_time)
print('Time Elapsed',timetaken[0],'mins',timetaken[1],'s')

#%%
def plot4():
    plt.hist(T_Vals[0],bins='auto',density=True,histtype = 'step',lw=1.5,label=Nvals[0])
    plt.hist(T_Vals[1],bins='auto',density=True,histtype = 'step',lw=1.5,label=Nvals[1])
    plt.hist(T_Vals[2],bins='auto',density=True,histtype = 'step',lw=1.5,label=Nvals[2])
    plt.hist(T_Vals[3],bins='auto',density=True,histtype = 'step',lw=1.5,label=Nvals[3])
    plt.plot([1.64,1.64],[0,0.7],c='black',label =  '5% significance line')
    plt.plot([1.28,1.28],[0,0.7],':',c='black',label = '10% significance line')
    
    plt.legend(title='$N_{d}$ ,$n_{mc}$ = 10$N_{d}$',loc='upper left',fontsize=8)
    plt.title('$Pull_{T}$ Histograms for different $n_{mc},N_{d}$,\n %s events per hist, $n_{k}$ = %s' %(noTs,myK))
    plt.xlabel('$Pull_{T}$')
    plt.ylabel('Frequency Density')
    
plot4()
pd.DataFrame(T_Vals).to_csv("LongTvals.csv")
#%%
# plt.hist(T_Vals,bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals)
plt.hist(T_Vals[0],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[0])
plt.hist(T_Vals[1],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[1])
plt.hist(T_Vals[2],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[2])
plt.hist(T_Vals[3],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[3])
plt.plot([1.64,1.64],[0,0.7],c='black',label =  '5% significance line')
plt.plot([1.28,1.28],[0,0.7],':',c='black',label = '10% significance line')

plt.legend()
plt.title('Pull T Histograms for different $n_{mc}$,\n %s events per hist, $n_{k}$ = %s,  $n_{d}$ = %s' %(noTs,myK,Nd))
plt.xlabel('Pull-T')
plt.ylabel('Frequency Density')
plt.savefig('PullTHistograms',format='pdf')





#%%
# k=3
# List = np.random.randint(low=1,high=10,size=(10,4))
# List = np.array(List)
# point = List[np.random.randint(0,len(List))]
# point = np.array(point)
# print('p',point)
# print('list',List)
# g = List - point

# dists = np.linalg.norm(g,axis=1)
# print('dis',dists)
# sortdists = np.sort(dists)
# print(sortdists)
# distsneeded = sortdists[1:k+1]
# print('dn',distsneeded)
# inds = []
# # for i in range(len(distsneeded)):
# #     inds.append(np.where(dists == distsneeded[i])[0][0])
    
# inds = [np.where(dists == val)[0][0] for val in distsneeded]
# print('indexes',inds)

