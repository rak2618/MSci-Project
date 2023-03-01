#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 12:04:49 2022

@author: robertoking
"""
import numpy as np
cimport numpy as cnp
import random
# cimport random as cr
import pandas as pd

import time

cdef unsigned long long int total
cdef int k
cdef float t1, t2, t

t1 = time.time()

for k in range(1000000000):
    total = total + k
print("Total =", total)

t2 = time.time()
t = t2-t1
print("%.100f" % t)



D = pd.read_csv('gof-master/neut_events_4D.csv')
M = pd.read_csv('gof-master/nuwro_events_4D.csv')

# cdef list Da
# cdef cnp.ndarray Mo
Da = list(D)
Mo = list(M)

D1 = [i[0] for i in Da]
D2 = [i[1] for i in Da]
D3 = [i[2] for i in Da]
D4 = [i[3] for i in Da]

M1 = [i[0] for i in Mo]
M2 = [i[1] for i in Mo]
M3 = [i[2] for i in Mo]
M4 = [i[3] for i in Mo]


cdef int Nd,Nmc
Nd  = 50
Nmc = 50




# cdef list  T_Vals 
# for i in range(200):
#     # start_time1 = time.time()
#     # print(i/NoofTs * 100,'%')
#     # 'Create Random Sample'
#     Dind = random.sample(range(len(D)), Nd)
#     Mind = random.sample(range(len(M)),Nmc)
    
#     DSample = np.array([Da[i] for i in Dind])
#     MSample = np.array([Mo[i] for i in Mind])
    
#     for j in range(len(DSample)):
#         DSample
#     # print(time.time()-start_time1,'s')
    
#     # T_Vals.append(MixedT(DSample,MSample,myK))
    

cpdef float dist(list p1,list p2):
    return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2+(p1[3]-p2[3])**2)
    
    
    
cpdef FNN(int k,list Point,list List1,list List2):
    
    cdef list Distlist1, Distlist2
    cdef int N1,N2 
    N1,N2 = len(List1),len(List2)
    Distlist1, Distlist2 = [] , []
    
    cdef int i
    for i in range(N1):
        Distlist1.append(dist(Point,List1[i]))
    for i in range(N2):
        Distlist2.append(dist(Point,List2[i]))
    
    BigList = np.concatenate((Distlist1,Distlist2))
    BigList.sort()
    
    DistNeeded = BigList[k+1]                           #Need k+1 th place since the zeroeth index is dist to point itself
    LargeList = np.concatenate((Distlist1,Distlist2))
    index = np.where(LargeList == DistNeeded)
    LongList = np.concatenate((List1,List2))
    
    return LongList[index]




















