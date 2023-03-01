#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 18:28:46 2021

@author: robertoking
"""

'Plotting Functions from packages'

import numpy as np
import matplotlib.pyplot as plt
import time


from scipy.stats import norm
from matplotlib import cm 
import matplotlib.colors as colors

G = [13.0,10.0, 12.0, 13.0, 11.0, 14.0, 21.0, 32.0, 40.0, 16.0, 10.0, 19.0, 20.0, 18.0, 11.0, 23.0, 11.0, 11.0, 28.0, 40.0, 55.0, 34.0, 10.0,
 49.0, 132.0, 195.0, 151.0, 79.0, 25.0, 61.0, 108.0, 88.0, 49.0, 10.0, 64.0, 89.0, 94.0, 67.0, 31.0, 121.0, 211.0, 232.0, 151.0, 329.0, 552.0,
 455.0, 221.0, 423.0, 357.0, 158.0, 329.0, 371.0, 244.0, 495.0, 587.0, 867.0, 696.0, 663.0, 774.0]

coords = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], [22, 23, 24, 25, 26, 27, 28, 29, 30],
 [31, 32], [33, 34], [35, 36, 37], [38], [39], [40], [41], [42, 43, 44], [45], [46], [47], [48, 49, 50, 51], [52, 53], [54], [55, 56],
 [57, 58, 59], [60], [61], [62], [63], [64], [65], [66], [67], [68], [69], [70], [71], [72], [73], [74], [75], [76], [77], [78],
 [79], [80], [81], [82], [83], [84], [85], [86], [87], [88], [89], [90], [91], [92], [93], [94], [95], [96], [97], [98], [99]]


#%%
def minutes(seconds):
    mins = np.floor(seconds/60)
    secs = (seconds/60 - mins)*60
    return [mins,secs]



#%%


def FWHM(List1d, nobins = 500):
    BinHeights,BinEdges = np.histogram(List1d,bins = nobins)
    HM = 0.5*max(BinHeights)  #Half Max value
    aboveHM = [i for i in BinHeights if i > HM]
    BinsaboveHM = len(aboveHM)  #Find number of bins that have values above half max
    BinWidth = BinEdges[1]-BinEdges[0]
    FWHM = BinsaboveHM*BinWidth
    return FWHM

def dist_2d_temp(point1,point2):     #takes 2 coordinates and finds their euclidean distance
    rsquared = (point1[0] - point2[0])**2 + (point1[1]-point2[1])**2
    return np.sqrt(rsquared)

def dist_2d(x1,y1,x2,y2):
    rsquared = (x1-x2)**2 + (y1-y2)**2
    return np.sqrt(rsquared)
    
def PoissonTest(f_obs,f_exp):
    
    Fpoisson = 0
    
    if len(f_obs) - len(f_exp) != 0:
        raise Exception('Lengths of Observed data and Expected data do not match')
        
    else:
        for i in range(len(f_obs)):
            Fpoisson += f_exp[i] - f_obs[i] + f_obs[i]*np.log(f_obs[i] / f_exp[i])
    
    return 2*Fpoisson

def norm_bins(x1,x2,Nbins):
    mu = (x2+x1)/2
    FMin = norm.cdf(x1-mu)
    FMax = norm.cdf(x2-mu)
    FValues = np.linspace(FMin,FMax,Nbins+1)
    # print('FValues are',FValues)
    xvals = norm.ppf(FValues)+mu
    # print('xvals are',xvals)
    return list(xvals)

def plot(func,x1,x2):
    x = np.linspace(x1,x2,10000)
    y = func(x)
    plt.plot(x,y,color='black')
    
def plotlines(func,xvals):
    yvals = func(xvals)
    for i in range(len(xvals)):
        X = [xvals[i],xvals[i]]
        Y = [0       ,yvals[i]]
        plt.plot(X,Y,color='red',ls=':')

def plotareas(func,x1,x2,Nbins):
    plot(func,x1,x2)
    plotlines(func,norm_bins(x1,x2,Nbins))

#%%
'2d spiral bins section'
def plotshapes(points,Heights):                     #takes list of lists of coordinates
    Heights = np.array(Heights)
    offset = Heights + np.abs(Heights.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.jet(norm(fracs.tolist()))
    print(color_values)
    
    plt.figure()
    for i in range(len(points)):
        points[i].append(points[i][0]) #repeat the first point to create a 'closed loop'
        
        xs, ys = zip(*points[i]) #create lists of x and y values
        plt.fill(xs, ys, color=color_values[i]) 
    plt.grid()
    plt.show()

def plotsquare(points,colour,ax):                     #takes list of coordinates
    points = list(points)
    points.append(points[0]) #repeat the first point to create a 'closed loop'
        
    xs, ys = zip(*points) #create lists of x and y values
    ax.fill(xs, ys, color=colour) 
    


def flatten_list(_2d_list):
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list   
    
def get_cell_vertices(cell):          #Input one cell coordinate

    Y  = cell[0]                       #First cell coordinate gives y value
    X  = cell[1]
    
    topleft     = [Y  ,X  ]
    topright    = [Y  ,X+1]
    bottomleft  = [Y+1,X  ]         #Y goes downwards
    bottomright = [Y+1,X+1]
    
    return [topleft,topright,bottomright,bottomleft]

def get_cells_vertices(cells):
    corners = []
    for i in cells:
        corners.append(get_cell_vertices(i))
    return get_list_of_lists(list(set(map(tuple,flatten_list(corners)))))

def get_list_of_lists(list_of_tuples):
    list_of_lists = []                                                          
    for tuple in list_of_tuples:
        list_of_lists.append(list(tuple))

    return list_of_lists

#%%
'Spiral section'
def spiral_cw(A):
    A = np.array(A)
    out = []
    while(A.size):
        out.append(A[0])        # take first row
        A = A[1:].T[::-1]       # cut off first row and rotate counterclockwise
    return np.concatenate(out)


def spiral_coords(pos,wid):
    
    if pos + 1 > wid**2 :
        raise Exception('Position index out of range of spiral')
    
    arr = []
    for i in range(wid):
        sublist = []
        for j in range(wid):
            sublist.append(wid*i+j)
        arr.append(sublist)
    
    arr_spiral = spiral_cw(arr)
    arr = np.array(arr)
    
    return list(np.ravel(np.where(arr == arr_spiral[pos])))

def spiral_bins(Heights):          #works on the unit grid
    X = spiral_cw(Heights)
    Bins = len(Heights)
    G = []
    coords = []
    coord = []
    
    
    val = 0
    for i in range(Bins**2):
        if val < 10:
            val += X[i]
            coord.append(i)
    
        else:
            G.append(val)
            val = X[i]
                
            coords.append(coord)
            coord = [i]
    
    if val != 0:            
        G.append(val)
        coords.append(coord)
        
    return G , coords

def transform(spiralcoords):                        
    N = int(np.sqrt(len(flatten_list(spiralcoords))))
    XYcoords = []
    for i in range(len(spiralcoords)):
        sublist = []
        for j in range(len(spiralcoords[i])):
            sublist.append(spiral_coords(spiralcoords[i][j],N))
        XYcoords.append(sublist)
    return XYcoords

def translate(coords,Xoffset,Yoffset,Xscalar,Yscalar):
    out = []
    for i in range(len(coords)):
        out.append([coords[i][0]/Xscalar - Xoffset, coords[i][1]/Yscalar - Yoffset])
            
    return out
            
def makegrid(Xl,Yl):
    fig = plt.figure(figsize = (5,5))
    gs = plt.GridSpec(1,1)
    
    ax1 = plt.subplot(gs[0])
    
    ax1.set_xlim(0,Xl)
    ax1.set_ylim(0,Yl)
    ax1.set_aspect('equal')
    ax1.set_yticks(np.arange(0,Yl,1))
    ax1.set_xticks(np.arange(0,Xl,1))



    
# makegrid(10,10)

Heights = [1,2,3,4]
Squares = [[0,0],[0,1],[1,1],[1,0]]
def fillgrid(Heights,Squares,Length):
    '''
    Heights is a 1d, Length^2 list with every bin colour
    Squares is a list of coordinates for each square, corresponding with height value
    Length is side length of the grid
    '''
    makegrid(Length,Length)
    
    Heights = np.array(Heights)
    offset = Heights + np.abs(Heights.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.viridis(norm(fracs.tolist()))
    # print(color_values)
    
    for i in range(Length**2):
        vertices = get_cell_vertices(Squares[i])
        plotsquare(vertices,colour = color_values[i])
    
        
    
def BinChilling(G,coords,density):
    hh = []
    pp = []
    for i in range(len(coords)):
        
        for j in range(len(coords[i])):
            if density:
                h = G[i]/len(coords[i])
            else:
                h = G[i]
            hh.append(h)
            pp.append(coords[i][j])
        
    return hh , pp

#%%
def DSquared(D):
    out = D**2
    return out

def D(D):
    out = D
    return out

def DRecip(D):
    out = 1/D
    return out

def DGauss(D):
    out = np.exp((-D**2)/2)
    return out


def scatterplot(PointsA,PointsAType,indexA,PointsB,PointsBType,indexB,N=50):
    '''
    Takes a list of vectors, only works in 2D
    Takes in data A at index A
    Takes in data B at index B
    Takes in N which is the number of sample points from the data
    '''
    
    if PointsAType == 'M':
        if indexA == 0:
            LabelA = 'MC Muon Energy'
        if indexA == 1:
            LabelA = 'MC Muon Angle'
        if indexA == 2:
            LabelA = 'MC Proton Energy'
        if indexA == 3:
            LabelA = 'MC Proton Angle'
            
    if PointsAType == 'D':
        if indexA == 0:
            LabelA = 'Data Muon Energy'
        if indexA == 1:
            LabelA = 'Data Muon Angle'
        if indexA == 2:
            LabelA = 'Data Proton Energy'
        if indexA == 3:
            LabelA = 'Data Proton Angle'
    
    # else:
    #     LabelA = 'No Label Found'
        
    if PointsBType == 'M':
        if indexB == 0:
            LabelB = 'MC Muon Energy'
        if indexB == 1:
            LabelB = 'MC Muon Angle'
        if indexB == 2:
            LabelB = 'MC Proton Energy'
        if indexB == 3:
            LabelB = 'MC Proton Angle'
            
    if PointsBType == 'D':
        if indexB == 0:
            LabelB = 'Data Muon Energy'
        if indexB == 1:
            LabelB = 'Data Muon Angle'
        if indexB == 2:
            LabelB = 'Data Proton Energy'
        if indexB == 3:
            LabelB = 'Data Proton Angle'
    
    # else:
    #     LabelB = 'No Label Found'
        
    
    fig = plt.figure(figsize = (6,5))
    fig.suptitle(LabelB + ' vs ' + LabelA)
    gs = plt.GridSpec(1,1)
    ax1 = plt.subplot(gs[0])
    
    x_points = [i[indexA] for i in PointsA[0:N]]
    y_points = [i[indexB] for i in PointsB[0:N]]
    
    ax1.scatter(x_points  ,y_points   ,color='blue',marker = 'x'  ,label = '%s Points' %(len(x_points)))
    ax1.set_xlabel(LabelA)
    ax1.set_ylabel(LabelB)
    ax1.grid()
    ax1.legend()
    
    
    
def histogram2d(PointsA,PointsAType,indexA,PointsB,PointsBType,indexB,NBins = 100,mycolours = 'jet'):
    '''
    Takes a list of vectors, only works in 2D
    Takes in data A at index A
    Takes in data B at index B
    Takes in N which is the number of sample points from the data
    '''
    
    if PointsAType == 'M':
        if indexA == 0:
            LabelA = 'MC $E_{\u03BC}$'             #muon energy
        if indexA == 1:
            LabelA = 'MC $cos(\u03F4_{\u03BC})$'   #muon angle
        if indexA == 2:
            LabelA = 'MC $E_{p}$'                  #proton energy
        if indexA == 3:
            LabelA = 'MC $cos(\u03F4_{p})$'        #proton angle
            
    if PointsAType == 'D':
        if indexA == 0:
            LabelA = 'Data $E_{\u03BC}$'
        if indexA == 1:
            LabelA = 'Data $cos(\u03F4_{\u03BC})$'
        if indexA == 2:
            LabelA = 'Data $E_{p}$'
        if indexA == 3:
            LabelA = 'Data $cos(\u03F4_{p})$'
    
    # else:
    #     LabelA = 'No Label Found'
        
    if PointsBType == 'M':
        if indexB == 0:
            LabelB = 'MC $E_{\u03BC}$'
        if indexB == 1:
            LabelB = 'MC $cos(\u03F4_{\u03BC})$'
        if indexB == 2:
            LabelB = 'MC $E_{p}$'  
        if indexB == 3:
            LabelB = 'MC $cos(\u03F4_{p})$'
            
    if PointsBType == 'D':
        if indexB == 0:
            LabelB = 'Data $E_{\u03BC}$'
        if indexB == 1:
            LabelB = 'Data $cos(\u03F4_{\u03BC})$'
        if indexB == 2:
            LabelB = 'Data $E_{p}$'
        if indexB == 3:
            LabelB = 'Data $cos(\u03F4_{p})$'
    
    # else:
    #     LabelB = 'No Label Found'
        
    
    # fig = plt.figure(figsize = (6,5))
    # gs = plt.GridSpec(1,1)
    # ax1 = plt.subplot(gs[0])
    
    x_points = [i[indexA] for i in PointsA]
    y_points = [i[indexB] for i in PointsB]
    
    plt.hist2d(x_points  ,y_points ,bins = NBins  ,cmap = mycolours)
    plt.title('2D Histogram of ' + LabelB + ' against ' + LabelA )#+ '\n Bins = (%s,%s)' %(NBins,NBins))
    plt.colorbar()
    plt.xlabel(LabelA)
    plt.ylabel(LabelB)
    
    
    
    
    

def Term1(Points,Psi):
    N = len(Points)
    out = 0
    
    for i in range(0,N):
        for j in range(i+1,N):
            dist = abs(Points[i]-Points[j])
            out += Psi(dist)
    
    return (1/(N*(N-1)))*out

def Term3(Data,MCPoints,Psi):
    N_Data = len(Data)
    N_MC   = len(MCPoints)
    out = 0
    
    for i in range(N_MC):
       for j in range(N_Data):
           dist = abs(Data[i]-MCPoints[j])
           out += Psi(dist)
    
    return (1/(N_MC*N_Data))*out

def T_Function(Data,MCPoints,Psi):
    term1 = Term1(Data,Psi)
    term2 = Term1(MCPoints,Psi)
    term3 = Term3(Data,MCPoints,Psi)
    
    return term1 + term2 + term3

#%%
'''
Need to design a 4d T-function that can be imported into any script

Will be designed to take in vectors for every point, rather than  
'''

def dist_Ndim(Point1,Point2):
    dim = len(Point1)
    
    Rsquared = 0
    
    for i in range(dim):
        Rsquared += (Point1[i]-Point2[i])**2
    
    return np.sqrt(Rsquared)

def dist_2_d(Point1,Point2):
    
    
    Rsquared = (Point1[0]-Point2[0])**2 + (Point1[1]-Point2[1])**2
    
    return np.sqrt(Rsquared)

def Term1_2d(Points,Psi):
    N = len(Points)
    out = 0
    
    for i in range(0,N):
        for j in range(i+1,N):
            out += Psi(dist_2_d(Points[i],Points[j]))
    
    return (1/(N*(N-1)))*out

def Term3_2d(DataPoints,MCPoints,Psi):
    N_MC = len(MCPoints)
    N_Data = len(DataPoints)
    
    out = 0
    for i in range(N_MC):
        for j in range(N_Data):
            out += Psi(dist_2_d(MCPoints[i],DataPoints[j]))
    
    return (1/(N_MC*N_Data))*out

def T_Function_2d(DataPoints,MCPoints,Psi):
    term1 = Term1_2d(DataPoints,Psi)
    term2 = Term1_2d(MCPoints,Psi)
    term3 = Term3_2d(DataPoints,MCPoints,Psi)
    
    return term1 + term2 - term3

#%%
def dist_4d(Point1,Point2):
    
    Rsquared = (Point1[0]-Point2[0])**2 + (Point1[1]-Point2[1])**2 + (Point1[2]-Point2[2])**2 + (Point1[3]-Point2[3])**2
    
    return np.sqrt(Rsquared)


def Term1_4d(Points,Psi):
    N = len(Points)
    out = 0
    
    for i in range(0,N):
        for j in range(i+1,N):
            out += Psi(dist_4d(Points[i],Points[j]))
    
    return (1/(N*(N-1)))*out

def Term3_4d(DataPoints,MCPoints,Psi):
    N_MC = len(MCPoints)
    N_Data = len(DataPoints)
    
    out = 0
    for i in range(N_MC):
        for j in range(N_Data):
            out += Psi(dist_4d(MCPoints[i],DataPoints[j]))
    
    return (1/(N_MC*N_Data))*out

def T_Function_4d(DataPoints,MCPoints,Psi):
    term1 = Term1_4d(DataPoints,Psi)
    term2 = Term1_4d(MCPoints,Psi)
    term3 = Term3_4d(DataPoints,MCPoints,Psi)
    
    return term1 + term2 - term3
        
        


# plot(norm.pdf,-3,3)
# plotlines(norm.pdf,createbins(-3,3,10))



# X = get_cells_vertices(transform(coords)[0])
# # plotshapes([X],[1])


# shapes = []

# for i in range(len(coords)):
#     shapes.append(get_cells_vertices(transform(coords)[i]))

# heights = np.arange(0,len(shapes),1)

# plotshapes([shapes[0]],heights)


# V = [3,7,2,9]
# U = [32,94,23,37]
# start = time.time()
# dist_Ndim(V,U)
# print(time.time()-start)


