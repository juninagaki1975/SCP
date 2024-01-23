#!/usr/bin/env python3
#
# Jun 2024.1.22
#

import sys,os
import numpy as np
import subprocess as sp

cmd = "rm ./out"
try:
    sp.check_output(cmd,shell=True)
except sp.CalledProcessError as e:
    print(e)


dorb = ["dxy","dyz","dz2","dxz","dx2my2"]
avg = 0.0
for k in range(0,5,2):
    print("k = ",k)
    n = 4
    for l in range(4):
        for j in range(l+1,n+1):
            enesum = 0.0
            for i in range(26):
                cmd = "./sli1.py NiO/unitcell/"+dorb[l]+"/r1.dat NiO/unitcell/"+dorb[l]+\
                    "/Pnl_"+str(i)+".dat NiO/unitcell/"+dorb[j]+"/Pnl_"+str(i)+".dat "+str(k)
                try:
                    ene = sp.check_output(cmd,shell=True)
                except sp.CalledProcessError as e:
                    print("error...")
                    sys.exit()

                enesum = enesum + float(ene)
            enesum = np.round(enesum / 26,6)
            avg = avg + enesum
            print(dorb[l],dorb[j],enesum)
    print(np.around(avg/10,6))
    print()
    avg = 0.0




