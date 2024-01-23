#!/usr/bin/env python3
#
# Jun 2024.1.16
#

import sys,os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy import interpolate as ip

from ioxsf import *
from numerical import *
from extract import *

aB = 0.52917720859
Ha = 27.21138505

def integrate(fname,l,m):

    nn = 401

    primvec,convvec,grids,org,span,symbol,natom,atomcoord,\
            wanval = read_xsf(fname)

#    wanvalsum = np.sum(np.abs(wanval))**2
#    wannorm = np.sqrt(wanvalsum)
#    wanval = wanval / wannorm

    val = extract(fname,l,m)
#    valsum = np.sum(np.abs(val))**2
#    norm = np.sqrt(valsum)
#    val = val / norm

    nx = grids[0]
    ny = grids[1]
    nz = grids[2]

    origin = np.zeros((3))
    stepvec = np.zeros((3,3))
    for i in range(3):
        stepvec[i,i] = span[i,i] / (grids[i]-1)
        origin[i] = -stepvec[i,i]

    ww = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ww = ww + np.abs(val[i,j,k])

# densities in xsf are already Bohr^-3 ?
#    vv = volume(stepvec)
    vv = volume(stepvec) / (aB*aB*aB)
    wfnorm = ww * vv
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                val[i,j,k] = val[i,j,k] / wfnorm

    cx = int(nx/2) + 1
    cy = int(ny/2) + 1
    cz = int(nz/2) + 1
    r = nx*stepvec[0,0] - cx*stepvec[0,0]
#    r = nx*stepvec[0,0]*span[0,0] -cx*stepvec[0,0]*span[0,0]
    step = np.linspace(0,r,cx-2)
    x = np.linspace(0,r,nn)
    pnl = np.zeros((26,nn))

#
# 26
#
#  6 : x,-x,y,-y,z,-z
# 12 : xy*4,yz*4,xz*4
#  8 : xyz,-x-y-z,-xyz,x-yz,xy-z
#    : -x-yz,-xy-z,x-y-z
#    : 

# 6
#100
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[i+cx,cy,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[0] = fit(x) 

#-100
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[cx-i,cy,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[1] = fit(x) 

#001
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[cx,cy,i+cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[2] = fit(x) 

#00-1
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[cx,cy,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[3] = fit(x) 


#010
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[cx,cy+i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[4] = fit(x) 

#0-10
    tmp = []
    n = 0
    for i in range(int(nx/2)-1):
        tmp.append([step[n],val[cx,cy-i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[5] = fit(x) 

    nr = int(int(nx/2)/np.sqrt(2))
    step = np.linspace(0,r,nr)
    x = np.linspace(0,r,nn)
#
# xy 
#110
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy+i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[6] = fit(x) 

#-1-10
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy-i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[7] = fit(x) 

#-110
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy+i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[8] = fit(x) 

#1-10
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy-i,cz]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[9] = fit(x) 

# yz
#011
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx,cy+i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[10] = fit(x) 

#0-1-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx,cy-i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[11] = fit(x) 

#0-11
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx,cy-i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[12] = fit(x) 

#01-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx,cy+i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[13] = fit(x) 

# xz
#101
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[14] = fit(x) 

#-10-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[15] = fit(x) 

#-101
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[16] = fit(x) 

#10-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[17] = fit(x) 

#
# xyz
#
#  8 : xyz,-x-y-z,-xyz,x-yz,xy-z
#    : -x-yz,-xy-z,x-y-z

    nr = int(int(nx/2)/np.sqrt(3))
    step = np.linspace(0,r,nr)
    x = np.linspace(0,r,nn)

#111
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy+i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[18] = fit(x) 

#-1-1-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy-i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[19] = fit(x) 

#-111
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy+i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[20] = fit(x) 

#1-11
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy-i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[21] = fit(x) 

#11-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy+i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[22] = fit(x) 

#-1-11
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy-i,cz+i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[23] = fit(x) 

#-11-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx-i,cy+i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[24] = fit(x) 

#1-1-1
    tmp = []
    n = 0
    for i in range(nr):
        tmp.append([step[n],val[cx+i,cy-i,cz-i]])
        n += 1
    tmp = np.array(tmp).astype(float)
    fit = ip.interp1d(tmp[:,0],tmp[:,1],kind='cubic')

#    plt.scatter(tmp[:,0],tmp[:,1])
#    plt.plot(x,fit(x))
    pnl[25] = fit(x) 

#    plt.imshow(val[14,:,:])

    pnl = pnl * Ha

    sumpnl = np.zeros((nn))
    for j in range(26):
        n = 0
        for i in range(nn):
            sumpnl[n] += pnl[j,n]
            n += 1
            
    avgpnl = sumpnl / 26 
#    np.savetxt("avgpnl.dat",avgpnl)
#    plt.scatter(x,avgpnl,color="red")
    plt.plot(x,avgpnl,color="red")

    for i in range(26):
#        plt.scatter(x,pnl[i],color="black")
        plt.plot(x,pnl[i],color="black")


#    plt.show()
    plt.savefig(fname.split(".")[0]+".jpg")

    return x,pnl



if __name__=="__main__":


    print("integration.")
    fname = sys.argv[1]
    l = int(sys.argv[2])
    m = int(sys.argv[3])
    x,pnl = integrate(sys.argv[1],l,m)

    for i in range(26):
        np.savetxt("Pnl_"+str(i)+".dat",pnl[i])
    np.savetxt("r1.dat",x/aB)


