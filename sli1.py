#!/usr/bin/env python3
#
# Jun 2023.12.29
#

import sys, os, glob
import numpy as np
import math 

# already angstrom in VESTA.
# aB = 0.52917721092
EPS9 = 1.E-9

def sli1(tmpr1,tmppnl1,tmpr2,tmppnl2,k):

    nn = len(tmpr1)
    r1 = np.zeros((nn))
    r2 = np.zeros((nn))
    pnl1 = np.zeros((nn))
    pnl2 = np.zeros((nn))
    for i in range(len(tmpr1)):
        r1[i] = float(tmpr1[i])
    for i in range(len(tmpr2)):
        r2[i] = float(tmpr2[i])
    for i in range(len(tmppnl1)):
        pnl1[i] = float(tmppnl1[i])
    for i in range(len(tmppnl2)):
        pnl2[i] = float(tmppnl2[i])

#    if r1 != r2:
#        print("We expect r1 = r2 at present.")
#        sys.exit()

    r1[0] = 0.001

    direct = True
    f1 = np.zeros((nn))
    f2 = np.zeros((nn))
    # test. compute on the same electron.
    if (direct):
        for i in range(nn):
            f1[i] = pnl1[i] ** 2
            f2[i] = pnl2[i] ** 2

#
# expression:
#
#  f(i) - f(i-1) = h/12 (5 a0 + a1 - a2)
#
# W.G. Bickley:
# 
#  y1 - y0 = h/12 (5 q0 + 8 q1 - q2) 
# 

#
# compute r1
#

    xi = np.zeros((nn))
    xj = np.zeros((nn))
    xa= np.zeros((nn))
    xb = np.zeros((nn))

    ho12 = r1[1] / 12.0 # h over 12. set for second pont.

# no need ni, 
    #ni = 40 # numbers of point per block?

# first point is already set as zero in above.
# initialize coefficients.

    a0 = 0.0
    b0 = 0.0

#
# calc r1 integral
#

    for i in range(2,nn,2):

        eras = 8.0 * f1[i-1]
        a2 = a0
        b2 = b0 
    # r^k+1 / r^1 
        a1 = eras * r1[i-1]**(k+1) / r1[i-1]
        b1 = eras / r1[i-1]**(k+1) 
        a2 = f1[i] * r1[i]**(k+1) / r1[i]
        b2 = f1[i] / r1[i]**(k+1)

        eras = ho12 * (5.0*a0+a1-a2)
        xi[i-1] = xi[i-2] + eras 
        xa[i-1] = xa[i-2] + np.abs(eras)

        eras = ho12 * (5.0*b0+b1-b2)
        xj[i-1] = xj[i-2] + eras 
        xb[i-1] = xb[i-2] + np.abs(eras)
        
        eras = ho12 * (-a0 + a1 + 5.0*a2)
        xi[i] = xi[i-1] + eras 
        xa[i] = xa[i-1] + np.abs(eras)

        eras = ho12 * (-b0 + b1 + 5.0*b2)
        xj[i] = xj[i-1] + eras 
        xb[i] = xb[i-1] + np.abs(eras)

    #    ni = ni - 2
    #    if ni <= 0:
    #        ho12 = ho12 + ho12
    #        ni = 40

    for i in range(2,nn):
        eras = r1[i]**(k+1) / r1[i]
        xi[i] = xi[i] / r1[i]**(k+1) + (xj[nn-1] - xj[i]) * eras
        xa[i] = xa[i] / r1[i]**(k+1) + (xb[nn-1] - xb[i]) * eras

    #
    # calc r2 integral
    #

        a0 = 0.0
        b0 = 0.0
        a1 = 0.0
        b1 = 0.0
        ho3 = r1[1] / 1.5
    #    ni = 40

        for i in range(2,nn,2):
            a2 = f2[i] * xi[i]
            b2 = np.abs(f2[i])* xa[i]

            eras = 4.0 * f2[i-2]
            a1 = a1 + ho3 * (a0 + eras*xi[i-1] + a2)
            b1 = b1 + ho3 * (b0 + np.abs(eras) * xa[i-1] + b2)
            a0 = a2
            b0 = b2 
    #        ni = ni -2
    #        if ni < 0:
    #            ho3 = ho3 + ho3 
    #            ni = 40

    return a1

if __name__=="__main__":

    r1fname = open(sys.argv[1],"r")
    pnl1fname = open(sys.argv[2],"r")
    pnl2fname = open(sys.argv[3],"r")
    k = int(sys.argv[4])
    rdat = r1fname.readlines()
    pnl1dat = pnl1fname.readlines()
    pnl2dat = pnl2fname.readlines()

#    print(rdat)
#    print(pnldat)

    a1 = sli1(rdat,pnl1dat,rdat,pnl2dat,k)
    print(a1)
