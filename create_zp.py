# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:17:50 2017

@author: sajid
"""

import numpy as np
import matplotlib.pyplot as plt
import urllib
import create_disc as c
 
def zone_radius(n,f,wavel):
    return np.sqrt(n*wavel*f + ((n*wavel)/2)**2)

def get_property(mat,energy):
    url = "http://henke.lbl.gov/cgi-bin/pert_cgi.pl"
    data = {'Element':str(mat), 'Energy':str(energy), 'submit':'Submit Query'}
    data = urllib.parse.urlencode(data)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    resp = urllib.request.urlopen(req)
    respDat = resp.read()
    response = respDat.split()
    d =  b'g/cm^3<li>Delta'
    i = response.index(d)
    delta = str(response[i+2])[:str(response[i+2]).index('<li>Beta')][2:]
    beta = str(response[i+4])[2:-1]
    return float(delta),float(beta)


mat = 'Au'
energy = 5000                   #Energy in EV
f = 6*10**(-3)                  #focal length in meters 
wavel = (1240/energy)*10**(-9)  #Wavelength in meters
delta,beta = get_property(mat,energy)

N = 100
n = 50
zp = np.zeros((N,N))
x =  np.linspace(-25e-6,25e-6,N)
X,Y = np.meshgrid(x,x)
step = x[-1]-x[-2]

N = 25 #number of zones
radius = np.zeros(N)

#2D ZP
for i in range(N):
    radius[i] = zone_radius(i,f,wavel)
    if i%2 == 1 :
        disc1 = c.create_disc(X,Y,step,radius[i-1],radius[i],n)
        zp = zp + disc1
    #print(i)
