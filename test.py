# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:28:29 2017

@author: sajid
"""

import numpy as np
from numba import jit,njit

def func_1(x1,y1,x2,y2,r1,r2,n):
     x11,y11 = np.meshgrid(np.linspace(x1,x2,n),np.linspace(y1,y2,n))
     z1 = np.sqrt(x11**2+y11**2)
     a = np.where((z1>(r1)) & (z1<(r2)))
     fill_factor = len(a[0])/(n*n)
     return fill_factor
 
@jit(nopython=True,parallel=True)
def func_2(x1,y1,x2,y2,r1,r2,n):
    x_ = np.linspace(x1,x2,n)
    y_ = np.linspace(y1,y2,n)
    cnts = 0
    for i in range(n):
        for j in range(n):
            z = (x_[i] * x_[i] + y_[j] * y_[j])
            if r1*r1 < z < r2*r2:
                cnts += 1
    fill_factor = cnts/(n*n)
    return fill_factor

@njit   # equivalent to "jit(nopython=True)".
def func_3(x1,y1,x2,y2,r1,r2,n):
    x_ = np.linspace(x1,x2,n)
    y_ = np.linspace(y1,y2,n)
    cnts = 0
    for i in range(n):
        for j in range(n):
            z = (x_[i] * x_[i] + y_[j] * y_[j])
            if r1*r1 < z < r2*r2:
                cnts += 1
    fill_factor = cnts/(n*n)
    return fill_factor


x1 = 1.0
x2 = -1.0
y1 = 1.0
y2 = -1.0
r1 = 0.5
r2 = 0.75
n = 2500