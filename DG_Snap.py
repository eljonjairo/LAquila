#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate Snapshoots of the surface propagation
#
# John Diaz August 2022

# To Do List:
#   Filtering 
#   Same colorbar
#   Save movie

import os
import warnings

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from pathlib import Path

from scipy import signal
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")
plt.close('all')

os.system('clear')

# DG folder
DGFolder = Path('../DGrun/Ref_13Stats/LAquilaCirella03_dhF500m_ID_1/')

# Space coordinates (km)
# LAquila 13 Stations
xini = 310.0
xend = 430.0
yini = 4620.0
yend = 4740.0

# Maximum Velocity for colorbar (m/s)
Vmax = 0.2

# Simulation Time (s), Spatial step (Km) and Time (s) step
time = 40.0
dx = 0.5
dt = 0.2

# Filtering parameters
lowcut  = 1.0                              # srate low pass cut frecuency
highcut = 0.01                             # srate cut high pass frecuency
fs = 1/dt                                  # Sample rate

print("  ")
print(" START PROGRAM ")
print("  ")

xend = xend+dx
yend = yend+dx
X = np.arange(xini,xend,dx)
Y = np.arange(yini,yend,dx)

nx = X.size
ny = Y.size

xMat, yMat = np.meshgrid(X,Y)
zMat = np.zeros((nx,ny))

# Read binary files from fortran
fnameDGx = DGFolder.joinpath('VX.snap')
#fnameDGy = DGFolder.joinpath('VY.snap')

print(" Loading snap files:\n ")
print(f" {fnameDGx}")
#print(f" {fnameDGy}")

Vx = np.fromfile(fnameDGx, dtype=np.float32)
#Vy = np.fromfile(fnameDGy, dtype=np.float32)

nsnap = int(Vx.size/(nx*ny))

Vx = np.reshape(Vx,(nx,ny,nsnap),order = 'F')
#Vy = np.reshape(Vy,(nx,ny,nsnap),order = 'F')

# Filtering Snapshots
# Coefs fol filtering
w = lowcut/(fs/2)                         # Normalize the frequency
b, a = signal.butter(4, w, 'low')
VxF = np.zeros((nx,ny,nsnap))
VyF = np.zeros((nx,ny,nsnap))

for ix in range (0,nx):
    for iy in range (0,ny):
        VxF[ix,iy,:] = signal.filtfilt(b, a, Vx[ix,iy,:])
#        VyF[ix,iy,:] = signal.filtfilt(b, a, Vy[ix,iy,:])

for isnap in range(0,nsnap):
    print(f" Processing snap {isnap} of {nsnap} ")
    fig = plt.figure(figsize=(12,12))
    ax = fig.subplots(1,1)
    V = VxF[:,:,isnap]
    Vm=ax.pcolormesh(xMat,yMat,V,cmap=cm.seismic,vmin=-Vmax,vmax=Vmax)
    plt.colorbar(Vm,label='Vx (m/s)',orientation="horizontal",shrink=.68)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel(' X (Km) ')
    ax.set_ylabel(' Y (Km) ')


























