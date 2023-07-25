#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Fault class
# 

import utm
from scipy.io import loadmat

import numpy as np
import matplotlib.pyplot as plt
import scipy 
import matplotlib.colors

import plotly.graph_objs as go
import plotly.tools as tls

import plotly.io as io
from plotly.offline import plot

from matplotlib import cm


class Fault():
    def __init__(self,name, dhF):
        self.name = name
        self.dhF = dhF

    def __str__(self):
        return " Fault name: " + str(self.name) + " dhF: " + str(self.dhF) + " Km" 

    def LoadMat(self,infile):
        print()
        print(f" Loading matlab file from {infile} ")
        print()

        # Load Matlab Finite Fault input file
        inFault = loadmat(infile)
        Fault = inFault['s2009LAQUIL03CIRE']
        
        # Load original fault coordinates, slip, rise time and rupture time
        self.ZFinMat = -Fault['geoZ'][0][0]            # Input Z coord (Km)  
        self.LatFinMat = Fault['geoLAT'][0][0]         # Input Lat coord (°)
        self.LonFinMat = Fault['geoLON'][0][0]         # Input Lon coord (°)
        ndip_in, nstk_in = self.ZFinMat.shape          # Input ndip, nstk
        self.SlipinMat  = Fault['slipSPL'][0][0]/100   # Input Slip in meters.
        self.RiseTinMat = Fault['riseSPL'][0][0]       # Input risetime (s)
        self.RupTinMat  = Fault['timeSPL'][0][0]       # Input rupture time (s)
        
        # Xutm, Yutm coord Km
        self.XFinMat, self.YFinMat, tmp1, tmp2 = \
        utm.from_latlon(self.LatFinMat, self.LonFinMat,33,'T')
        self.XFinMat = self.XFinMat/1000
        self.YFinMat = self.YFinMat/1000
        self.XFin3D = self.XFinMat.flatten(order='F').transpose()
        self.YFin3D = self.YFinMat.flatten(order='F').transpose()
        self.ZFin3D = self.ZFinMat.flatten(order='F').transpose()
        self.Slipin3D = self.SlipinMat.flatten(order='F').transpose()
        self.RiseTin3D = self.RiseTinMat.flatten(order='F').transpose()
        self.RupTin3D = self.RiseTinMat.flatten(order='F').transpose()
        
        # Hypocenter coordinates (Km) 
        self.hypolon = Fault['evLON'][0][0][0][0]   # Hypocenter lon coord (°)
        self.hypolat = Fault['evLAT'][0][0][0][0]   # Hypocenter lat coord (°)
        self.hypoz = -Fault['evDPT'][0][0][0][0]    # Hypocenter lon coord (Km)
        
        # Hypocenter Xutm, Yutm coord (Km)
        self.hypox, self.hypoy, tmp1, tmp2 = utm.from_latlon(self.hypolat, 
                                                             self.hypolon,33,'T')
        self.hypox = self.hypox/1000
        self.hypoy = self.hypoy/1000
        
        # Calculate magnitude and unitary vector of delta in strike direction 
        vec_dstk = np.array([ self.XFinMat[0,1]-self.XFinMat[0,0], 
                              self.YFinMat[0,1]-self.YFinMat[0,0],
                              self.ZFinMat[0,1]-self.ZFinMat[0,0] ])
        dstk_in = np.linalg.norm(vec_dstk)
        self.univec_dstk = vec_dstk/dstk_in

        # Calculate magnitude and unitary vector of delta in dip direction 
        vec_ddip = np.array([ self.XFinMat[1,0]-self.XFinMat[0,0], 
                              self.YFinMat[1,0]-self.YFinMat[0,0],
                              self.ZFinMat[1,0]-self.ZFinMat[0,0] ])
        ddip_in = np.linalg.norm(vec_ddip)
        self.univec_ddip = vec_ddip/ddip_in

        # Calculate output delta vector in strike and dip directions
        self.dstkVec = self.univec_dstk*self.dhF
        self.ddipVec = self.univec_ddip*self.dhF
        self.dstk = np.linalg.norm(self.dstkVec)
        self.ddip = np.linalg.norm(self.ddipVec)

        # Calculate length and output number of points in strike an dip direction        
        stk_len = round((nstk_in-1)*dstk_in)
        dip_len = round((ndip_in-1)*ddip_in)
        self.nstk = int(stk_len/self.dhF)+1
        self.ndip = int(dip_len/self.dhF)+1
        
        # Define input and output strike and dip vectors
        self.dipinVec = np.linspace(0, dip_len, ndip_in)
        self.stkinVec = np.linspace(0, stk_len, nstk_in)
        self.stkVec = np.linspace(0, stk_len, self.nstk)
        self.dipVec = np.linspace(0, dip_len, self.ndip)
        
        stk_len = round((self.nstk-1)*self.dstk)
        dip_len = round((self.ndip-1)*self.ddip)

        print(" Original Fault Dimensions:")
        print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " 
              %(stk_len, nstk_in, dstk_in) )
        print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" 
              %(dip_len, ndip_in, ddip_in) )
        
        print(" Hypocenter Coordinates x, y and z (Km): %6.2f %6.2f %6.2f " 
              %(self.hypox, self.hypoy, self.hypoz) )
        print()
        print(" Output Fault Dimensions:")
        print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " 
              %(stk_len, self.nstk, self.dstk) )
        print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" 
              %(dip_len, self.ndip, self.ddip) )

    def PlotXYSlipin(self):

        fig, ax = plt.subplots()
        ax.set_xlabel(" X (Km)")
        ax.set_ylabel(" Y (Km)")
        ax.set_aspect('equal',adjustable='box')
        fp=ax.pcolormesh(self.XFinMat/1000, self.YFinMat/1000, self.SlipinMat, 
                         cmap=cm.viridis)
        plt.colorbar(fp,location='right', label=" Slip (m) ", shrink=.6)
        ax.set_title(" Input Slip ")
        plt.show()

    def PlotXYZSlipin(self, azim, dist, elev):
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface( self.XFinMat, self.YFinMat, self.ZFinMat,  
                          facecolors=cm.viridis(self.SlipinMat), linewidth=0, 
                          antialiased=True )
        ax.scatter(self.hypox, self.hypoy, self.hypoz, color='red', marker='*')
        ax.set_xlabel(" X (Km)")
        ax.set_ylabel(" Y (Km)")
        ax.set_zlabel(" Z (Km)")
        ax.set_aspect('equal',adjustable='box') 
        ax.set_title(" Input Slip ")
        # Azimuth angle, distance an elevation of camera 
        ax.azim = azim                                           
        ax.dist = dist 
        ax.elev = elev 
        plt.show()
       
    def InterpolateXYZCoords(self):
        # Coordinates of the first fault point
        inivec = np.array([ self.XFinMat[0,0], self.YFinMat[0,0],
                            self.ZFinMat[0,0] ])

        # Creta arrays for interpolated coordinates
        XFMat = np.zeros((self.ndip, self.nstk))
        YFMat = np.zeros((self.ndip, self.nstk))
        ZFMat = np.zeros((self.ndip, self.nstk))

        # Initiate arrays to calculate the index of hypocenter in interpolated
        # XYZ coordinates
        hypod = np.zeros((self.ndip, self.nstk))       # Hypocenter distance
        hypoxy = [self.hypox, self.hypoy]

        for istk in range (0, self.nstk):
            delta_stk = istk*self.dstkVec
            for idip in range (0, self.ndip):
                delta_dip = idip*self.ddipVec
                XFMat[idip,istk] = inivec[0] + delta_stk[0] + delta_dip[0]
                YFMat[idip,istk] = inivec[1] + delta_stk[1] + delta_dip[1]
                ZFMat[idip,istk] = inivec[2] + delta_stk[2] + delta_dip[2]

                # Calculate hypocenter distance in XY coords.
                XYvec = [XFMat[idip,istk],YFMat[idip,istk]]
                hypod[idip,istk] = scipy.spatial.distance.euclidean(XYvec,hypoxy)


        # Find the indexes of the hypocenter
        self.hypoidip, self.hypoistk = np.where(hypod == np.min(hypod))

        # Calculate the rigth z fault coords using hypocenter as reference,
        # the hypocenter have to be over the fault
        zmov = self.hypoz - ZFMat[self.hypoidip, self.hypoistk]
        ZFMat = ZFMat + zmov
        
        # From matrix to column vector following fortran 
        XF3D = XFMat.flatten(order='F').transpose()
        YF3D = YFMat.flatten(order='F').transpose()
        ZF3D = ZFMat.flatten(order='F').transpose()

        self.XFMat = XFMat
        self.YFMat = YFMat
        self.ZFMat = ZFMat
        self.XF3D = XF3D
        self.YF3D = YF3D
        self.ZF3D = ZF3D
        self.fcoor = np.array((XF3D,YF3D,ZF3D)).transpose()
        
    def InterpolateSlip(self):
        # Slip Interpolation
        SlipF = scipy.interpolate.interp2d(self.stkinVec, self.dipinVec, 
                                     self.SlipinMat, kind = "cubic")
        self.SlipMat = SlipF(self.stkVec, self.dipVec)

    def PlotXYZSlip(self, azim, dist, elev):
       
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface( self.XFMat, self.YFMat, self.ZFMat,  
                         facecolors=cm.viridis(self.SlipMat), linewidth=0, 
                         antialiased=False )
        ax.scatter(self.hypox, self.hypoy, self.hypoz, color='red', marker='*')
        ax.set_xlabel(" X (Km)")
        ax.set_ylabel(" Y (Km)")
        ax.set_zlabel(" Z (Km)")
        ax.set_aspect('equal',adjustable='box') 
        ax.set_title(" Interpolated Slip ")
        # Azimuth angle, distance an elevation of camera 
        ax.azim = azim                                           
        ax.dist = dist 
        ax.elev = elev 
        plt.show()
        
        
    def PlotlySlip(self):
        io.renderers.default='svg'
        
        
        #fig = go.Surface(z=self.SlipMat, colorscale='Viridis',showscale=True, 
        #                 lighting=dict(ambient=0.9))
        #fig.show()
        fig = go.Figure(data=[go.Surface(x=self.XFMat, y=self.YFMat, 
                                         z=self.ZFMat, colorscale="Viridis", 
                                         surfacecolor=self.SlipMat, lighting=dict(ambient=0.7))])
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.25, y=-1.25, z=0.2))

        fig.update_layout(scene_camera=camera, title="Interpolated Slip")
        fig.show()
        
    def PlotXYSlip(self):

        fig, ax = plt.subplots()
        ax.set_xlabel(" X (Km)")
        ax.set_ylabel(" Y (Km)")
        ax.set_aspect('equal',adjustable='box')
        fp=ax.pcolormesh(self.XFMat/1000, self.YFMat/1000, self.SlipMat, 
                         cmap=cm.viridis)
        plt.colorbar(fp,location='right', label=" Slip (m) ", shrink=.6)
        ax.set_title(" Interpolated Slip ")
        plt.show()   