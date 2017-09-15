# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 15:45:23 2017

@author: sajid
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import njit

def repeat_pattern(X,Y,Z):
    flag_ = np.where((X>0)&(Y>0))
    flag1 = np.where((X>0)&(Y<0))
    flag1 = tuple((flag1[0][::-1],flag1[1]))
    Z[flag1] = Z[flag_]
    flag2 = np.where((X<0)&(Y>0))
    flag2 = tuple((flag2[0],flag2[1][::-1]))
    Z[flag2] = Z[flag_]
    flag3 = np.where((X<0)&(Y<0))
    flag3 = tuple((flag3[0][::-1],flag3[1][::-1]))
    Z[flag3] = Z[flag_]
    return Z

@njit   # equivalent to "jit(nopython=True)".
def partial_fill(x,y,step,r1,r2,n):
    x_ = np.linspace(x-step/2,x+step/2,n)
    y_ = np.linspace(y-step/2,y+step/2,n)
    cnts = 0
    for i in range(n):
        for j in range(n):
            z = (x_[i] * x_[i] + y_[j] * y_[j])
            if r1*r1 < z < r2*r2:
                cnts += 1
    fill_factor = cnts/(n*n)
    return fill_factor

#find the radius of the nth zone
def zone_radius(n,f,wavel):
    return np.sqrt(n*wavel*f + ((n*wavel)/2)**2)

def make_ring(X,Y,flag,r1,r2,n):
    z = np.zeros(np.shape(X))
    Z = np.sqrt(X**2+Y**2)
    for l in range(len(flag[0])):
        i = flag[0][l]
        j = flag[1][l]
        if 0.75*r1< Z[i][j] < 1.25*r2:
            x1 = X[i][j]
            y1 = Y[i][j]
            z[i][j] = partial_fill(x1,y1,step,r1,r2,n)
    z[tuple((flag[1],flag[0]))] = z[tuple((flag[0],flag[1]))]
    return z

#based upon page 215 of Chris's book

f = 6*10**(-3)             #focal length in meters 
e = 5000                   #Energy in EV
wavel = (1240/e)*10**(-9)   #Wavelength in meters

#Gold(Au) @ 5kEV, from Henke et.al
delta = 1.2143*10**(-4)
beta  = 2.6233*10**(-5) 

n = 7500
zp = np.zeros((n,n))
x =  np.linspace(-50e-6,50e-6,n)
X,Y = np.meshgrid(x,x)
step = x[-1]-x[-2]
flag = np.where((X>0)&(Y>0)&(X>=Y))

N = 833 #number of zones
radius = np.zeros(N)
for k in range(N):
    radius[k] = zone_radius(k,f,wavel)

#2D ZP
zp = np.zeros(np.shape(X))
for i in range(N):
    if i%2 == 1 :
        print(i)
        r1 = radius[i-1]
        r2 = radius[i]
        ring = make_ring(X,Y,flag,r1,r2,250)
        zp+=ring
zp = repeat_pattern(X,Y,zp)