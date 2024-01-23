#!/usr/bin/env python3
#
# Jun 2024.1.12
#

import sys,os
import numpy as np
from ioxsf import *
from numerical import *

aB = 0.52917720859


__all__=['extract']

def extract(fname,l,m):

    primvec,convvec,grids,org,span,symbol,natom,atomcoord,\
            wanval = read_xsf(fname)

    origin = np.zeros((3))
    stepvec = np.zeros((3,3))
    for i in range(3):
        stepvec[i,i] = 1.0 / grids[i]
        origin[i] = -stepvec[i,i]

    val = np.zeros((grids[0],grids[1],grids[2]))
    mx = -(grids[0])/2 
    my = -(grids[1])/2 
    mz = -(grids[2])/2 

    ix = mx
    iy = my
    iz = mz

    n = 0
    for k in range(grids[2]):
        if (iz == grids[2]/2):
            iz = mz
        for j in range(grids[1]):
            if (iy == float(grids[1]/2)):
                iy = my
            for i in range(grids[0]):
                if (ix == float(grids[0]/2)):
                    ix = mx
                x = float(ix) * stepvec[0,0] + origin[0]
                y = float(iy) * stepvec[1,1] + origin[1]
                z = float(iz) * stepvec[2,2] + origin[2]
#origin.                
#                if ix == 0 and iy == 0 and iz == 0:
#                    print("origin ",i*stepvec[0,0],j*stepvec[1,1],k*stepvec[2,2])
#                    print("origin ",i,j,k)

                r,theta,phi = cart2sph(x,y,z)
                zlm = tesseral(l,m,theta,phi)
                val[i,j,k] = np.round(zlm * wanval[i,j,k],12)
                n += 1
                ix += 1
            iy += 1
        iz += 1

    write_xsf(fname,primvec,convvec,natom,symbol,atomcoord,\
            grids,org,span,val)

    return val


if __name__=="__main__":


    print("extraction.")
    l = int(sys.argv[2])
    m = int(sys.argv[3])
    extract(sys.argv[1],l,m)
