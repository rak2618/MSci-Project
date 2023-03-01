#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 12:19:31 2022

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

#%%



def I(point1,point2,List1,List2):
    
    # print('p1',point1)
    # print(List1)
    if point1 in List1:
        if point2 in List1:
            out = 1
        else:
            out = 0
    
    if point1 in List2:
        if point2 in List2:
            out = 1
        else:
            out = 0
    
    return out

def findnearestneighbour(k,Point,List1,List2):
    '''
    

    Parameters
    ----------
    k : integer
        The kth closest point to Point will be found.
    Point : VEctor in the form of a list
        The target point.
    List1 : Array of Vectors
        The function will be checking the points in the lists for the closest one.
    List2 : Array of Vectors
        The function will be checking the points in the lists for the closest one.

    Returns
    -------
    The closest point and which list it was in.

    '''
    
    Distlist1 = []
    Distlist2 = []
    
    for i in range(len(List1)):
        Distlist1.append(mf.dist_4d(Point,List1[i]))
    for i in range(len(List2)):
        Distlist2.append(mf.dist_4d(Point,List2[i]))
    
    BigList = np.concatenate((Distlist1,Distlist2))
    BigList.sort()
    
    DistNeeded = BigList[k+1]                           #Need k+1 th place since the zeroeth index is dist to point itself
    LargeList = np.concatenate((Distlist1,Distlist2))
    index = np.where(LargeList == DistNeeded)
    LongList = np.concatenate((List1,List2))
    
    
    return LongList[index]

def FNN(k,Point,List1,List2):
    '''
    Parameters
    ----------
    k : integer
        The kth closest point to Point will be found.
    Point : VEctor in the form of a list
        The target point.
    List1 : Array of Vectors
        The function will be checking the points in the lists for the closest one.
    List2 : Array of Vectors
        The function will be checking the points in the lists for the closest one.

    Returns
    -------
    The closest point and which list it was in.

    '''
    
    Distlist1 = []
    Distlist2 = []
    
    for i in range(len(List1)):
        Distlist1.append(mf.dist_4d(Point,List1[i]))
    for i in range(len(List2)):
        Distlist2.append(mf.dist_4d(Point,List2[i]))
    
    # print(Distlist1[0])
    
    BigList = np.concatenate((Distlist1,Distlist2))
    BigList.sort()
    
    DistNeeded = BigList[k+1]                           #Need k+1 th place since the zeroeth index is dist to point itself
    LargeList = np.concatenate((Distlist1,Distlist2))
    index = np.where(LargeList == DistNeeded)
    LongList = np.concatenate((List1,List2))
    
    
    return LongList[index]

def findnearestneighbour(k,Point,List1,List2):
    '''
    Parameters
    ----------
    k : integer
        The kth closest point to Point will be found.
    Point : VEctor in the form of a list
        The target point.
    List1 : Array of Vectors
        The function will be checking the points in the lists for the closest one.
    List2 : Array of Vectors
        The function will be checking the points in the lists for the closest one.

    Returns
    -------
    The closest point and which list it was in.

    '''
    
    Distlist1 = [np.sqrt((Point[0]-s[0])**2+(Point[1]-s[1])**2+(Point[2]-s[2])**2+(Point[3]-s[3])**2) for s in List1]
    Distlist2 = [np.sqrt((Point[0]-s[0])**2+(Point[1]-s[1])**2+(Point[2]-s[2])**2+(Point[3]-s[3])**2) for s in List2]
    # print(Distlist1[0])
    
    BigList = np.concatenate((Distlist1,Distlist2))
    BigList.sort()
    
    DistNeeded = BigList[k+1]                           #Need k+1 th place since the zeroeth index is dist to point itself
    LargeList = np.concatenate((Distlist1,Distlist2))
    index = np.where(LargeList == DistNeeded)
    LongList = np.concatenate((List1,List2))
    
    
    return LongList[index]

def FNN_rescaled(k,Point,List1,List2):
    
    Distlist1 = []
    Distlist2 = []
    
    for i in range(len(List1)):
        Distlist1.append(mf.dist_4d(Point,List1[i]))
    for i in range(len(List2)):
        Distlist2.append(mf.dist_4d(Point,List2[i]))
    
    BigList = np.concatenate((Distlist1,Distlist2))
    BigList.sort()
    
    DistNeeded = BigList[k+1]                           #Need k+1 th place since the zeroeth index is dist to point itself
    LargeList = np.concatenate((Distlist1,Distlist2))
    index = np.where(LargeList == DistNeeded)
    LongList = np.concatenate((List1,List2))
    
    return 0

def MixedT(ListA,ListB,K):
    Tvalue = 0
    for i in np.concatenate((ListA,ListB)):
        for j in range(K):
            P1 = i
            P2 = findnearestneighbour(j,P1,ListA,ListB)
            
            # print('p1',P1)
            # print('p2',P2)
            
            Ival = I(P1,P2,ListA,ListB)
            # print(Ival)
            
            Tvalue += Ival
    
    na = len(ListA)
    nb = len(ListB)
    
    return (1/(K*(na+nb)))*Tvalue

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
    




#%%


Nd  = 100          #Data sample size
Nmc = 100          #MCarlo sample size
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

NoofTs = 0
for i in range(NoofTs):
    # start_time1 = time.time()
    print(i/NoofTs * 100,'%')
    'Create Random Sample'
    Dind = random.sample(range(len(D)), Nd)
    Mind = random.sample(range(len(M)),Nmc)
    
    DSample = np.array([D[i] for i in Dind])
    MSample = np.array([M[i] for i in Mind])
    # print(time.time()-start_time1,'s')
    
    # T_Vals.append(MixedT(DSample,MSample,myK))
    T_Vals.append(PullT(DSample,MSample,myK))
    
    
print('Time Elapsed',time.time()-start_time,'s')
    

#%%
plt.hist(T_Vals,bins='auto',density = False,histtype = 'step',lw=2)
plt.xlabel('Pull-T')
plt.ylabel('Frequency')
plt.title('Mixed Sampled Method T-Distribution, Neut v Nuwro,\n Hist Events = %s,$ N_{d} = %s,N_{mc} = %s, N_{k} = %s $' %(NoofTs,Nd,Nmc,myK))
# plt.savefig('MixedSamp1')


#%%
'Testing the effect of K on the T value'

start_time = time.time()

T_Vals = []
for Kval in np.arange(5,30,5):
    PermTvals = []
    for i in range(0):
        'Create Random Sample'
        Dind = random.sample(range(len(D)), Nd)
        Mind = random.sample(range(len(M)),Nmc)
        
        DSample = np.array([D[i] for i in Dind])
        MSample = np.array([M[i] for i in Mind])
        
        # PermTvals.append(MixedT(DSample,MSample,Kval))
        PermTvals.append(PullT(DSample,MSample,Kval))
        
        print('Kval =',Kval,'Tval =',i)
    
    T_Vals.append(PermTvals)
    
print('Time Elapsed',time.time()-start_time,'s')


#%%
plt.hist(T_Vals,bins=10,histtype = 'step',lw=1.5,label=[i for i in np.arange(5,30,5)])
plt.legend()
plt.title('T Histograms for different $n_{k}$ ,200 events per hist')
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

myK = 2
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
# plt.hist(T_Vals,bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals)
plt.hist(T_Vals[0],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[0])
plt.hist(T_Vals[1],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[1])
plt.hist(T_Vals[2],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[2])
plt.hist(T_Vals[3],bins='auto',density=True,histtype = 'step',lw=1.5,label=NMCVals[3])
plt.legend()
plt.title('T Histograms for different $n_{mc}$, %s events per hist, $n_{k}$ = %s,  $n_{d}$ = %s' %(noTs,myK,Nd))
plt.xlabel('Pull-T')
plt.ylabel('Frequency Density')
plt.savefig('PullTHistograms',format='pdf')





































