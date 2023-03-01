#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:30:54 2022

@author: robertoking
"""

import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import time


def cumsum1(List1d):
    List1d.sort()
    x = np.repeat(List1d,2)
    x = np.append(x,x[-1]+1)
    y = np.arange(0,len(List1d)+1)             #Y counts from 0 to number of data points
    y = np.repeat(y,2)
    y = np.delete(y,0)
    # y = np.delete(y,-1)
    return x , y

def cumsum2(List1d):
    List1d.sort()
    x = np.repeat(List1d,1)
    y = np.arange(0,len(List1d))             #Y counts from 0 to number of data points
    return x , y

def cumsum(List1d):
    List1d.sort()
    x = np.repeat(List1d,1)
    y = np.arange(0,len(List1d))/len(List1d)             #Y counts from 0 to 1 and broadcast to data points
    return x , y

def nearestabove(array,value):
    if value > array.max():
        raise Exception('Value is greater than max array value')
    indexs = np.where(array >= value)
    index = indexs[0][0]
    return index

def Kolm(DataPoints,MCPoints):
    if len(DataPoints) != len(MCPoints):
        raise Exception('Datapoints and MC Points need to have equal length')
    
    N = len(DataPoints)
    Dxvals , Dyvals = cumsum(DataPoints)
    Mxvals , Myvals = cumsum(MCPoints)
    
    Dxvals = Dxvals[Dxvals < Mxvals.max()]
    Dyvals = Dyvals[:len(Dxvals)]
    
    NewYvals = []
    for i in range(len(Dxvals)-1):                  #Finds y values in the MC list in line with the Data x values
        
        xVal = Dxvals[i]
        if xVal > Mxvals.max():
            raise Exception('Value is greater than max array value')
        Mindexs = np.where(Mxvals >= xVal)
        Mindex = Mindexs[0][0]
        
        # Mindex = nearestabove(Mxvals,xVal)
        
        x0 = Mxvals[Mindex-1]
        x1 = Mxvals[Mindex]
        y0 = Myvals[Mindex-1]
        y1 = Myvals[Mindex]
        
        yprime = (y1-y0)*(xVal-x0)/(x1-x0)+y0
        NewYvals.append(yprime)
    
    Nprime = len(NewYvals)
    diff = np.array(NewYvals)-Dyvals[:Nprime]
    absdiff = abs(diff)
    D = max(absdiff)
    d = np.sqrt(N)*D
    
    return d
    
    
    
    
#%%


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
N=100
hd = plt.hist(D1[0:N],cumulative=True,histtype='step',density=True,bins=N)
hm = plt.hist(M1[0:N],cumulative=True,histtype='step',density=True,bins=N)
plt.plot(D1[0:N],hd[0] - hm[0])
#%%
N = 1000  #Sample size

# plt.hist(D1,bins=400)
# plt.hist(M1,bins=400)
M1 = M1[0:len(D1)]


Dxvals , Dyvals = cumsum(random.sample(D1,N))
Mxvals , Myvals = cumsum(random.sample(M1,N))

Dxvals = Dxvals[Dxvals < Mxvals.max()]
Dyvals = Dyvals[:len(Dxvals)]

# plt.plot(Dxvals,Dyvals,label='Data Events')
# plt.plot(Mxvals,Myvals,label='MC Events')
# plt.xlabel('Energy Value (GeV)')
# plt.ylabel('Cumulative Events')
# plt.title('Number of Events at or below a Muon Energy')
# plt.legend()

#%%
start = time.time()
NewYvals = []
for i in range(len(Dxvals)-1):
    print(100*i/N,'% finished')#,len(Dxvals)-1)
    xVal = Dxvals[i]
    Mindexs = np.where(Mxvals >= xVal)
    Mindex = Mindexs[0][0]
    
    x0 = Mxvals[Mindex-1]
    x1 = Mxvals[Mindex]
    y0 = Myvals[Mindex-1]
    y1 = Myvals[Mindex]
    
    yprime = (y1-y0)*(xVal-x0)/(x1-x0)+y0
    NewYvals.append(yprime)

print("Time Elapsed =",time.time()-start,"s")
#%%
Nprime = len(NewYvals)
diff = np.array(NewYvals)-Dyvals[:Nprime]
absdiff = abs(diff)
D = max(absdiff)
print('D',D)
d = np.sqrt(N)*D
print('d =',format(d,'.4f'))

plt.plot(Dxvals,Dyvals,label='Neut')
plt.plot(Dxvals[:Nprime],np.array(NewYvals),label='Nuwro')
# plt.plot(Mxvals,Myvals,label='MC Events')
plt.plot(Dxvals[:Nprime],diff,label='Difference')
plt.plot([-1,10],[0,0],c='black',ls=':')
plt.xlim(-0.5,8)
plt.xlabel('$ cos(\u03B8_{p}) $')
plt.ylabel('Cumulative Probability')
plt.title('Probability of event at or below $ cos(\u03B8_{p}) $ \n Sample Size = 1000')
plt.legend()


#%%
NValues4Emu = [100,1000,10000,100000,50000,20000,30000,40,30,80000,60000,70000,40000,90000,5000]
dvalues4Emu = [1.311358,1.138524,2.994089,9.952921,7.139264,4.208383,5.762638,0.570025,0.758920,9.278446,7.725843,8.781403,6.4025,9.584634,1.884489]


NValues4thetamu = [600,700,800,1000,2000,3000,4000,5000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000]
dvalues4thetamu = [1.1542,0.9495,0.9236,0.8903,1.1065,1.5087,1.2702,1.2492,1.2226,1.5233,2.1031,2.4350,2.4839,2.8844,3.1869,3.5326,3.6869,3.6596]

NValues4Ep = [80    ,90    ,100   ,200   ,300   ,400   ,500   ,600   ,700   ,800   ,1000  ,2000  ,3000  ,4000  ,5000  ,10000 ,20000 ,30000 ,40000 ,50000 ,60000 ,70000  ,80000  ,90000  ,100000 ]
dvalues4Ep = [1.2357,1.0596,1.4052,1.0719,2.0299,2.4098,1.7226,1.9296,1.8621,1.6834,2.0271,2.4168,2.5079,2.8927,3.0540,4.1903,6.2514,6.9691,8.3366,8.8957,9.6003,10.1203,10.9799,11.5084,12.3725]

NValues4Ptheta = [100   ,200   ,300   ,400   ,500   ,1000  ,2000  ,3000  ,4000  , 5000 ,6000  ,7000  ,8000  ,9000  ,10000 ,20000 ,30000 ,40000 ,50000 ,60000 ,70000 ,80000 ,90000 ,100000]
dvalues4Ptheta = [3.5555,2.6958,3.3072,1.7140,1.3063,0.9669,1.1965,1.3782,1.2042,1.3188,1.4220,1.7174,1.9384,2.0278,2.0631,3.0917,4.2938,5.3485,5.7703,6.2654,6.9761,7.2305,7.5369,8.0666]

linxvals = np.linspace(0,100000,200)

plt.plot(NValues4thetamu,dvalues4thetamu,'o')
plt.plot(linxvals,0.012*np.sqrt(linxvals),label = '$ d = 0.012\sqrt{N} $')
plt.plot([-10000,120000],[1.358,1.358],'--',label='5% Significance level')
plt.grid()
plt.xlim([-5000,110000])
plt.ylim(0)
plt.ylabel('Value of d')
plt.xlabel('Sample Size')
plt.legend()
plt.title('Certainty of having different data against Sample Size \n for Muon Angle')
    
    
    


#%%

 # maxMx = Mxvals.max()
        # maxDx = Dxvals.max()
        # if MaxDx > MaxMx:
        #     popindices = np.where(Dxvals > maxMx)
        
            






