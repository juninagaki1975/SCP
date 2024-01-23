#!/usr/bin/env python3
#
# Jun 2024.1.12
#

import sys,os
import numpy as np
from ioxsf import read_xsf
import math as mt

__all__ = ['volume','sph2cart','cart2sph','tesseral']

EPS9 = 1.0E-9

def volume(vec):

    vv = np.abs(vec[0][0] * (vec[1][1]*vec[2][2] - vec[2][1]*vec[1][2]) - \
                vec[0][1] * (vec[1][0]*vec[2][2] - vec[2][0]*vec[1][2]) + \
                vec[0][2] * (vec[1][0]*vec[2][1] - vec[2][0]*vec[1][1]) )

    return vv

def sph2cart(r,theta,phi):

    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    return x,y,z

def cart2sph(x,y,z):

    r = np.sqrt(x**2+y**2+z**2)

    if x == 0:
        if y > 0:
            phi = np.pi/2
        elif y < 0:
            phi = np.pi*3/2
        else:
            phi = 0
    elif (x > 0 and y > 0):
        phi = np.arctan(y/x)
    elif (x < 0 and y > 0) :
        phi = np.arctan(y/x) + np.pi
    elif (x < 0 and y <= 0) :
        phi = np.arctan(y/x) + np.pi
    elif (x > 0 and y <= 0) :
        phi = np.arctan(y/x) + (2*np.pi)
    if z == 0.0:
        theta = np.pi / 2
    elif z > 0:
        theta = np.arctan(np.sqrt(x**2+y**2)/z)
    elif z < 0:
        theta = np.arctan(np.sqrt(x**2+y**2)/z) + np.pi
#debug        
#    print("theta",np.rad2deg(theta))
#    print("phi",np.rad2deg(phi))


# original algorithm. studied from QE.
#
#    r = np.sqrt(x**2 + y**2 + z**2)
#    if (r < EPS9):
#        theta = np.pi / 2.0
#    else:
#        theta = np.arccos(z/r)
#    if (x > EPS9 or x < EPS9):
#        phi = mt.atan2(y,x)
#    else:
#        phi = np.sin(y) * np.pi / 2.0

    return r, theta, phi


def tesseral(l,m,theta,phi):

    if (l == 0 and m == 0): # s
        Zlm = 1.0 / 2.0 * np.sqrt(np.pi)

    elif (l == 1 and m == -1): # py
        Zlm = 0.5*np.sqrt(3.0/np.pi)*np.sin(theta)*np.sin(phi)

    elif (l == 1 and m == 0): # pz
        Zlm = 0.5*np.sqrt(3.0/np.pi)*np.cos(theta)

    elif (l == 1 and m == 1): # px
        Zlm = 0.5*np.sqrt(3.0/np.pi)*np.sin(theta)*np.cos(phi)

    elif (l == 2 and m == -2): # dxy
        Zlm = 0.25*np.sqrt(15.0/np.pi)*np.sin(theta)**2 *np.sin(2*phi)

    elif (l == 2 and m == -1): # dyz
        Zlm = 0.5*np.sqrt(15.0/np.pi)*np.sin(theta)*np.cos(theta)*np.sin(phi)

    elif (l == 2 and m == 0): # dz2
#        Zlm = 0.5*np.sqrt(5.0/(2.0*np.pi))*(3.0*np.cos(theta)**2+1.0)
        Zlm = 0.125*np.sqrt(5.0/np.pi)*(3.0*np.cos(2.0*theta)+1.0)

    elif (l == 2 and m == 1): # dxz
        Zlm = 0.5*np.sqrt(15.0/np.pi)*(np.sin(theta)*np.cos(theta)*np.cos(phi))

    elif (l == 2 and m == 2): # dx2-y2
        Zlm = 0.25*np.sqrt(15.0/np.pi)*(np.sin(theta)**2 *np.cos(2*theta))

    return Zlm


def axis_rotation(theta, axis_name):

    c = np.cos(theta)
    s = np.sin(theta)
    if axis_name =='x':
        rotation_matrix = np.array([[1, 0, 0],
                                    [0, c, -s],
                                    [0, s, c]])
    if axis_name =='y':
        rotation_matrix = np.array([[c, 0, s],
                                    [0, 1, 0],
                                    [-s, 0, c]])
    elif axis_name =='z':
        rotation_matrix = np.array([[c, -s, 0],
                                    [s, c, 0],
                                    [0, 0, 1]])
    return rotation_matrix


if __name__=="__main__":

    print("numerical.")
