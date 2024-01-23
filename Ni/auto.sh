#!/bin/bash

./integrate.py Ni_up_dxy.xsf 2 -2
wait
mv *dat ./Ni/unitcell/dxy/
wait

./integrate.py Ni_up_dyz.xsf 2 -1
wait
mv *dat ./Ni/unitcell/dyz/
wait

./integrate.py Ni_up_dz2.xsf 2 0
wait
mv *dat ./Ni/unitcell/dz2/
wait

./integrate.py Ni_up_dxz.xsf 2 1
wait
mv *dat ./Ni/unitcell/dxz/
wait

./integrate.py Ni_up_dx2my2.xsf 2 2
wait
mv *dat ./Ni/unitcell/dx2my2/
wait
