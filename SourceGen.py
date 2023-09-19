#!/usr/bin/env python3

#
# Generate sliprate files for DG
#
# John Diaz January 2023
#
# To Do List:
#   
#    Sliprate animation
#   

import os
import warnings

import numpy as np
import pickle

from scipy import signal
from scipy.io import FortranFile

import Source
import plotly.io as io

io.renderers.default='svg'

# Load Fault Object
inFault = '../Outputs/3DFaults/LAquilaCirella03_dhF500m.pickle'

# Output files directory
outDir = "Outputs/Source/"

# Output files directory
outDir = "Outputs/Source/"

# ID of simulation
ID = 0

# Source Inputs
tmax = 10.0;
dt = 0.01;
nt = int(tmax/dt)+1
lowcut  = 1.0                              # srate low pass cut frecuency
highcut = 0.01                             # srate cut high pass frecuency
fn = 1/dt/2                                # Nyquist frecuency (Sample rate/2)

# Magnitud Moment
Mw = 6.0

# Rake angle Cirella 2012
rake = -90

# Number of broadband escenarios
nscen = 4
# Characteristic asperity size (Km)
asp_size = 5.0
#Corner wave number
kc = 5.0

if __name__ == '__main__':

    warnings.filterwarnings("ignore")
 
    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Source Program      * ")
    print(" ************************************************ ")
    print()


    # Load the Fault Dict
    with open(inFault, 'rb') as handle:
        Fault = pickle.load(handle)

    LAquilaSource = Source.Source(Fault, Mw, asp_size, kc, rake)
  
    print(LAquilaSource)

    # Velocity Model H(Km),Vp(Km/s),Vs(Km/s) and Rho(Kg/m3) Bianchi, Ameri*
    h  =  np.array([ 1.00, 1.00, 3.00, 22.0, 15.0, 99.0 ])
    vp =  np.array([ 3.16, 4.83, 5.76, 6.51, 7.00, 7.80 ])
    vs =  np.array([ 1.70, 2.60, 3.10, 3.50, 3.80, 4.20 ])
    rho = np.array([ 2500, 2840, 2940, 3150, 3260, 3500 ])

    LAquilaSource.set_mu_1d(vs, rho, h)
    LAquilaSource.adjust_Mw()
    
    
# nstk = Fault['nstk']
# ndip = Fault['ndip']

# nflt = nstk*ndip
# SlipMat  = Fault['SlipMat']
# RiseTMat = Fault['RiseTMat']
# RupTMat = Fault['RupTMat']
# Slipmean = np.mean(SlipMat)

# stkMat = Fault['stkMat']
# dipMat = Fault['dipMat']

# it_start = np.round(RupTMat/dt)
# it_start = it_start.astype(int)
# it_end   = np.round((RupTMat+RiseTMat)/dt)
# it_end   = it_end.astype(int)

# # Sliprate calculation
# SRate = np.zeros([ndip,nstk,nt])
# SRateFilt = np.zeros([ndip,nstk,nt])
# maxsr = np.zeros([ndip,nstk])
# # Filtered Sliprate
# b, a = signal.butter(4, lowcut/fn)

# for idip in range (0,ndip):
#     for istk in range (0,nstk):
#         iti = it_start[idip,istk]
#         itf = it_end[idip,istk]
#         SRate[idip,istk,iti:itf] = SlipMat[idip,istk]/RiseTMat[idip,istk]
#         SRatetmp = SRate[idip,istk,:]
#         SRateFilt[idip,istk,:] = signal.filtfilt(b, a, SRatetmp)
#         maxsr[idip,istk] = np.amax(SRate[idip,istk,:])

# tfig = np.array([0,1,2,3,4,5,7])
# ntfig = int(tfig.size)
# itfig = np.zeros([ntfig,])
# for it in range (0,ntfig):
#     itfig[it] = int(tfig[it]/dt)


# fig, axes = plt.subplots(nrows=8, ncols=1,sharex='col', sharey='row')
# for i,ax in enumerate(axes.flat):
#     it = int(i/dt)
#     itime = ("%3.1f" %(it*dt) )
#     SRatefig = SRateFilt[:,:,it]
#     im=ax.pcolormesh(stkMat,dipMat,SRatefig,cmap="hot_r",vmin=0,vmax=0.3)
#     ax.text(32,7,itime,fontsize=8,fontweight='bold',color='black')
# fig.colorbar(im, ax=axes.ravel().tolist())

# # Slip positive in the direction of dip and in the direction of strike
# SRdip = (-SRate.flatten(order='F'))*np.sin(np.deg2rad(rake))
# SRstk = (SRate.flatten(order='F'))*np.cos(np.deg2rad(rake))
# SRdip = np.float32(SRdip)
# SRstk = np.float32(SRstk)
# Nsr = np.arange(0,SRdip.size,1)

# outname = Fault['outname']

# SRdipName = '../Outputs/Source/srate_dip_'+outname+'_ID_'+str(ID)
# SRstkName = '../Outputs/Source/srate_str_'+outname+'_ID_'+str(ID)

# fsrd = FortranFile(SRdipName, 'w')
# fsrd.write_record(SRdip)
# fsrs = FortranFile(SRstkName, 'w')
# fsrs.write_record(SRstk)

# # Write fcoor file
# fcoorHeader = "%d  %d %4.2f " %(nflt, nt, dt)
# fcoor = np.array(Fault['fcoor'])*1000
# fcoorName = '../Outputs/3DFaults/fcoor_'+outname+'_ID_'+str(ID)+'.in'

# with open(fcoorName,'wb') as f:
#     np.savetxt(f, fcoor, header=fcoorHeader, comments=' ',fmt = '%9.4f')

# print(f" Coordinates saved in file: {fcoorName}" )
# print(f" SlipRate dip saved in file: {SRdipName}" )
# print(f" SlipRate stk saved in file: {SRstkName}" )

# # Write slip file
# slipHeader = "%d %d %d" %(nstk,ndip,nflt)

# slip = SlipMat.flatten(order='F')
# rise = RiseTMat.flatten(order='F')
# rupt = RupTMat.flatten(order='F')

# total = np.array([slip,rise,rupt]).transpose()

# print(slip[0])
# print(rise[0])
# print(rupt[0])
# slipName = '../Outputs/3DFaults/slip_'+outname+'_ID_'+str(ID)+'.in'

# with open(slipName,'wb') as f:
#     np.savetxt(f, total, header=slipHeader, comments=' ',fmt = '%9.4f')

# print(f" slip saved in file: {slipName}" )



