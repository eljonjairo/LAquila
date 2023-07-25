#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 12:04:09 2022

@author: jon
"""

import os 
import warnings
    
import numpy as np
from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")

os.system('clear')
# m to Km
m = 1000

# outVelo = 'Outputs/Velo/ItalyCenterBianchi_1000m.txt'
outVelo = '../Outputs/Velo/ItalyCenterHerrmannCIA_1000m.txt'

# X and Y UTM limits (Km) for Model
xmin = 290000
xmax = 450000
ymin = 4550000
ymax = 4830000
zmin = -60000

dh = 4000
dz = 500

print("  ")
print(" START PROGRAM ")
print("  ")

xmax = xmax+dh
ymax = ymax+dh
x = np.arange(xmin,xmax,dh)
y = np.arange(ymin,ymax,dh)
z = np.arange(zmin,dz,dz)

xMat, yMat = np.meshgrid(x, y, indexing='xy')

xVec = xMat.flatten(order='F')
yVec = yMat.flatten(order='F')

# # Bianchi* (Tomado de Ameri 2012)
# # Velocity layers in z: 1.0, 2.0, 5.0 27.0 42.0
# H   = np.array([ 0,    1.0,  2.0,  5.0,  27.0, 42.0 ])
# Vp  = np.array([ 3.16, 4.83, 5.76, 6.51, 7.00, 7.80 ])
# Vs  = np.array([ 1.70, 2.60, 3.10, 3.50, 3.80, 4.20 ])
# Rho = np.array([ 2500, 2840, 2940, 3150, 3260, 3500 ])

# Herrmann initial Model
# Velocity layers in z: 1.5, 4.5, 7.5, 14.5
# H   = np.array([ 0.00, 1.50, 4.50, 7.50, 14.50, 29.50, 35.50, 43.50 ])
# Vp  = np.array([ 3.75, 5.00, 6.00, 6.30,  6.00,  6.70,  7.20,  7.90 ])
# Vs  = np.array([ 2.14, 2.86, 3.43, 3.57,  3.43,  3.78,  3.94,  4.40 ])
# Rho = np.array([ 2275, 2515, 2687, 2754,  2687,  2850,  2956,  3212 ])

# Herrmann CIA Model
# Velocity layers in z: 1.5, 4.5, 7.5, 14.5, 29.5, 35.5, 43.5

H   = np.array([ 0.00, 1.50, 4.50, 7.50, 14.50, 29.50, 35.50, 43.50 ])
Vp  = np.array([ 3.75, 4.94, 6.01, 5.55,  5.88,  7.11,  7.10,  7.90 ])
Vs  = np.array([ 2.14, 2.82, 3.43, 3.15,  3.36,  4.01,  3.99,  4.40 ])
Rho = np.array([ 2275, 2485, 2706, 2609,  2677,  3010,  3012,  3276 ])

Zh = -H*m
nZh = int(Zh.size)-1

nx = int(x.size)
ny = int(y.size)
nz = int(z.size)
n = nx*ny*nz

with open (outVelo,'w') as f:
       f.write("%d %d %d %d \n" %(n,nx,ny,nz))
       for iz in range(0,nz):
           print(f' {iz+1} of {nz}')
           for izh in range(0,nZh):
               if( z[iz] <= Zh[izh] and z[iz] > Zh[izh+1] ):
                   m1 = Vp[izh]*1000;
                   m2 = Vs[izh]*1000;
                   m3 = Rho[izh];
               if( z[iz] <= Zh[nZh] ):
                   m1 = Vp[nZh]*1000;
                   m2 = Vs[nZh]*1000;
                   m3 = Rho[nZh];
                   
           for iy in range(0,ny):
               for ix in range(0,nx):
                   f.write("%8.3f %8.3f %8.3f %5.1f %5.1f %5.1f \n" %(x[ix],y[iy],z[iz],m1,m2,m3))
                   
                

                   
#       f.write("%d \n" %(ntriF))
#       for itri in range(0,ntriF):
#           f.write("%6d %6d %12.3f #Set Maximum area on Facets (1) \n"\
#                   %(itri+1,itri+1,area))


f.close()
print()
print(f' {outVelo} file created ...')



