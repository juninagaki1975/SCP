#!/bin/bash

./integrate.py NiO_up_dxy.xsf 2 -2
wait
mv *dat ./NiO/unitcell/dxy/
wait

./integrate.py NiO_up_dyz.xsf 2 -1
wait
mv *dat ./NiO/unitcell/dyz/
wait

./integrate.py NiO_up_dz2.xsf 2 0
wait
mv *dat ./NiO/unitcell/dz2/
wait

./integrate.py NiO_up_dxz.xsf 2 1
wait
mv *dat ./NiO/unitcell/dxz/
wait

./integrate.py NiO_up_dx2my2.xsf 2 2
wait
mv *dat ./NiO/unitcell/dx2my2/
wait
