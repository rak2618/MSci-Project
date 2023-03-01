#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:46:00 2022

@author: robertoking
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MyFunctions as mf
import itertools as it
from bisect import bisect
import time

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

#%%
Ld = 100000
Lmc= 100

Dsamp = D[:Ld ]
Msamp = M[:Lmc]

def bin_4d_datapoints(points,dimbins = 5):
    wvals = [i[0] for i in points]
    xvals = [i[1] for i in points]
    yvals = [i[2] for i in points]
    zvals = [i[3] for i in points]
    
    wbinedges = np.linspace(min(wvals),max(wvals),dimbins+1)
    xbinedges = np.linspace(min(xvals),max(xvals),dimbins+1)
    ybinedges = np.linspace(min(yvals),max(yvals),dimbins+1)
    zbinedges = np.linspace(min(zvals),max(zvals),dimbins+1)
    
    bin_edges = np.concatenate((wbinedges,xbinedges,ybinedges,zbinedges))
    
    binheights = np.zeros(shape=(dimbins,dimbins,dimbins,dimbins))
    for i in range(len(points)):
        p = points[i]
        # print(p[0])
        # print(wbinedges.max())
        wbin = bisect(wbinedges,p[0]) - 1
        if wbin == dimbins:
            wbin -= 1
        xbin = bisect(xbinedges,p[1]) - 1
        if xbin == dimbins:
            xbin -= 1
        ybin = bisect(ybinedges,p[2]) - 1
        if ybin == dimbins:
            ybin -= 1
        zbin = bisect(zbinedges,p[3]) - 1
        if zbin == dimbins:
            zbin -= 1
        # print(wbin)
        
        binheights[wbin,xbin,ybin,zbin] += 1
    
    return binheights,bin_edges


bh = bin_4d_datapoints(Dsamp)[0]

def listbins(points,dimbins = 10):
    wvals = [i[0] for i in points]
    xvals = [i[1] for i in points]
    yvals = [i[2] for i in points]
    zvals = [i[3] for i in points]
    
    wbinedges = np.linspace(min(wvals),max(wvals),dimbins+1)
    xbinedges = np.linspace(min(xvals),max(xvals),dimbins+1)
    ybinedges = np.linspace(min(yvals),max(yvals),dimbins+1)
    zbinedges = np.linspace(min(zvals),max(zvals),dimbins+1)
    
    bin_edges = np.concatenate((wbinedges,xbinedges,ybinedges,zbinedges))
    
    binheights = np.zeros(shape=(dimbins,dimbins,dimbins,dimbins))
    for i in range(len(points)):
        p = points[i]
        # print(p[0])
        # print(wbinedges.max())
        wbin = bisect(wbinedges,p[0]) - 1
        if wbin == dimbins:
            wbin -= 1
        xbin = bisect(xbinedges,p[1]) - 1
        if xbin == dimbins:
            xbin -= 1
        ybin = bisect(ybinedges,p[2]) - 1
        if ybin == dimbins:
            ybin -= 1
        zbin = bisect(zbinedges,p[3]) - 1
        if zbin == dimbins:
            zbin -= 1
        # print(wbin)
        
        binheights[wbin,xbin,ybin,zbin] += 1
    
    bh = np.ravel(binheights)
    g = np.zeros(shape=(dimbins**4,2))
    g[:,0] = range(0,dimbins**4,1)
    g[:,1] = bh[:]
    
    return g

gh = listbins(Dsamp)


def mergebins(array2d):
    'Takes in a 2d array of indices and bin heights'
    H = np.ones(shape = (len(array2d),3))
    'add extra colum for bin size'
    H[:,0] = array2d[:,0]
    H[:,1] = array2d[:,1]
    
    '''
    To continue:
    Sort H by column index 1 so by bin size. 
    Merge first n bins so they add up to at least 10. 
    Resort
    Merge next bins so they add to 10
    Keep going until no bins have less than 10 bin size
        '''
    Hsorted = H[H[:,1].argsort()]
    
    return Hsorted

hs = mergebins(gh)


