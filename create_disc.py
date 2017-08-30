# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 11:17:12 2017

@author: sajid
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_region(sides_x,sides_y,r,step,n):
     plt.plot([sides_x[0],sides_x[0]],[sides_y[0],sides_y[1]],'g-')
     plt.plot([sides_x[1],sides_x[1]],[sides_y[0],sides_y[1]],'g-')
     plt.plot([sides_x[0],sides_x[1]],[sides_y[0],sides_y[0]],'g-')
     plt.plot([sides_x[0],sides_x[1]],[sides_y[1],sides_y[1]],'g-')
     x1 = np.linspace(sides_x[0],sides_x[1],n)
     y1 = np.linspace(sides_y[0],sides_y[1],n)
     x11,y11 = np.meshgrid(x1,y1)
     z1 = np.sqrt(x11**2+y11**2)
     #a = np.where((z1>(r-step/20)) & (z1<(r+step/20)) )
     a = np.where(z1<(r)) 
     plt.scatter(x1[a[0]],y1[a[1]])     
    

def get_partial_fill(sides_x,sides_y,r1,r2,step,n):
     x1 = np.linspace(sides_x[0],sides_x[1],n)
     y1 = np.linspace(sides_y[0],sides_y[1],n)
     x11,y11 = np.meshgrid(x1,y1)
     z1 = np.sqrt(x11**2+y11**2)
     a = np.where((z1>(r1)) & (z1<(r2)))
     fill_factor = len(a[0])/(n*n)
     return fill_factor

def create_disc(X,Y,step,r1,r2,n):
    Z = np.sqrt(X**2+Y**2)
    Z1 = np.zeros(np.shape(Z))
    l = np.shape(Z)[0]
    for i in range(l):
        for j in range(l):
            #the function get_partial_fill gets the value of the pixel at each point
            #by calcuating how much of the structure is within the box represented
            #by the pixel
            x1 = X[i][j]
            y1 = Y[i][j]
            #get boundaries of bounding cell
            sides_x = np.array([x1-step/2,x1+step/2])
            sides_y = np.array([y1-step/2,y1+step/2])
            #Note that r1 must be less than r2
            Z1[i][j] = get_partial_fill(sides_x,sides_y,r1,r2,step,n)
    return Z1