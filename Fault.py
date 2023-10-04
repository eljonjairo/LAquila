#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Fault class
# 

import utm
from scipy.io import loadmat

import numpy as np
import scipy 

import plotly.graph_objs as go
import plotly.io as io
import plotly.express as px
import plotly.tools as tls
import plotly.figure_factory as ff


import pickle

class Fault():
    def __init__(self,name, dhF):
        self.dhF = dhF
        self.name = name + "_dhF" + str(int(self.dhF)) + "m"

    def __str__(self):
        return " Fault name: " + str(self.name) + " dhF: " + str(self.dhF) + " Km" 

    def load_mat_file(self,infile):
        print()
        print(f" Loading matlab file from {infile} ")
        print()

        # Load Matlab Finite Fault input file
        inFault = loadmat(infile)
        Fault = inFault['s2009LAQUIL03CIRE']
        
        # Load original fault coordinates, slip, rise time and rupture time
        self.ZFinMat = -Fault['geoZ'][0][0]*1000        # Input Z coord (m)  
        self.LatFinMat = Fault['geoLAT'][0][0]          # Input Lat coord (°)
        self.LonFinMat = Fault['geoLON'][0][0]          # Input Lon coord (°)
        self.ndip_in, self.nstk_in = self.ZFinMat.shape # Input ndip, nstk
        self.SlipinMat  = Fault['slipSPL'][0][0]/100    # Input Slip in meters.
        self.RiseTinMat = Fault['riseSPL'][0][0]        # Input risetime (s)
        self.RupTinMat  = Fault['timeSPL'][0][0]        # Input rupture time (s)
   
        # Hypocenter coordinates (Km) 
        self.hypolon = Fault['evLON'][0][0][0][0]     # Hypocenter lon coord (°)
        self.hypolat = Fault['evLAT'][0][0][0][0]     # Hypocenter lat coord (°)
        self.hypoz = -Fault['evDPT'][0][0][0][0]*1000 # Hypocenter lon coord (m)
   
        # Hypocenter Xutm, Yutm coord (Km)
        self.hypox, self.hypoy, tmp1, tmp2 = utm.from_latlon(self.hypolat, 
                                                        self.hypolon,33,'T')
        self.hypox = self.hypox
        self.hypoy = self.hypoy 
   
        # Xutm, Yutm coord Km
        self.XFinMat, self.YFinMat, tmp1, tmp2 = \
        utm.from_latlon(self.LatFinMat, self.LonFinMat,33,'T')
        
        # Calculate magnitude and unitary vector of delta in strike direction 
        vec_dstk = np.array([ self.XFinMat[0,1]-self.XFinMat[0,0], 
                              self.YFinMat[0,1]-self.YFinMat[0,0],
                              self.ZFinMat[0,1]-self.ZFinMat[0,0] ])
        self.dstk_in = np.linalg.norm(vec_dstk)
        
        # Calculate magnitude and unitary vector of delta in dip direction 
        vec_ddip = np.array([ self.XFinMat[1,0]-self.XFinMat[0,0], 
                              self.YFinMat[1,0]-self.YFinMat[0,0],
                              self.ZFinMat[1,0]-self.ZFinMat[0,0] ])
        self.ddip_in = np.linalg.norm(vec_ddip)
        
        # Calculate length and output number of points in strike an dip direction        
        self.stk_len_in = round((self.nstk_in-1)*self.dstk_in)
        self.dip_len_in = round((self.ndip_in-1)*self.ddip_in)
          
        print(" Original Fault Dimensions:")
        print(" Strike (Km): %6.3f nstk: %d dstk (Km): %6.3f " 
              %(self.stk_len_in/1000, self.nstk_in, self.dstk_in/1000) )
        print(" Dip (Km): %6.3f ndip: %d ddip (Km): %6.3f" 
              %(self.dip_len_in/1000, self.ndip_in, self.ddip_in/1000) )
        print(f" Subfaults: {len(self.ZFinMat.flatten(order='F'))}")
    
        print(" Hypocenter Coordinates x, y and z (Km): %6.3f %6.3f %6.3f " 
              %(self.hypox, self.hypoy, self.hypoz) )
        print()     
   
    def set_full_fault(self):   
         
        # Xutm, Yutm coord Km
        self.XFinMat, self.YFinMat, tmp1, tmp2 = \
        utm.from_latlon(self.LatFinMat, self.LonFinMat,33,'T')
        self.XFin3D = self.XFinMat.flatten(order='F').transpose()
        self.YFin3D = self.YFinMat.flatten(order='F').transpose()
        self.ZFin3D = self.ZFinMat.flatten(order='F').transpose()
        self.Slipin3D = self.SlipinMat.flatten(order='F').transpose()
        self.RiseTin3D = self.RiseTinMat.flatten(order='F').transpose()
        self.RupTin3D = self.RiseTinMat.flatten(order='F').transpose()
        
        
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

        # Calculate output fault's size
        self.nstk = int(self.stk_len_in/self.dhF)+1
        self.ndip = int(self.dip_len_in/self.dhF)+1
        
        # Define input and output strike and dip vectors
        self.dipinVec = np.linspace(0, self.dip_len_in, self.ndip_in)
        self.stkinVec = np.linspace(0, self.stk_len_in, self.nstk_in)
        self.stkVec = np.linspace(0, self.stk_len_in, self.nstk)
        self.dipVec = np.linspace(0, self.dip_len_in, self.ndip)
        
        stk_len = round((self.nstk-1)*self.dstk)
        dip_len = round((self.ndip-1)*self.ddip)
       
        print(" Output Fault Dimensions:")
        print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " 
              %(stk_len/1000, self.nstk, self.dstk/1000) )
        print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" 
              %(dip_len/1000, self.ndip, self.ddip/1000) )
      
        
    def effective_size(self):
        
        slip = self.SlipinMat
       
        # Cumulative slip in strike and dip direction
        sum_slip_stk = np.sum(slip, axis=0)
        sum_slip_dip = np.sum(slip, axis=1)
       
        # Calculate the integral, and the effective length following
        # Source Scaling Properties from Finite-Fault-Rupture Models (Mai & Beroza 2000)
        acorr_stk = np.correlate(sum_slip_stk, sum_slip_stk, mode = 'full')
        x_stk = np.arange(0, len(acorr_stk), 1, dtype=int)*self.dstk_in
        Wacf_stk = scipy.integrate.simpson(acorr_stk, x_stk)/max(acorr_stk)
         
        acorr_dip = np.correlate(sum_slip_dip,sum_slip_dip, mode = 'full')
        x_dip = np.arange(0, len(acorr_dip), 1, dtype=int)*self.ddip_in
        Wacf_dip = scipy.integrate.simpson(acorr_dip, x_dip)/max(acorr_dip)
        
        print(" Effective length in strike direction: %5.2f Km." % (Wacf_stk) )
        print(" Effective length in dip direction: %5.2f Km." % (Wacf_dip) )
     
        # Calculate effective number of nodes
        nstk_eff = round(Wacf_stk/self.dstk_in)+1
        ndip_eff = round(Wacf_dip/self.ddip_in)+1
           
        # Select the interval with length Leff in the strike dir which
        # maximize the slip sum
        istk = 0
        jstk = nstk_eff
        istk_eff = istk
        sum_slip_eff = np.sum(sum_slip_stk[istk:jstk])
        while( jstk < self.nstk_in ):
             istk+=1
             jstk+=1
             sum_slip_tmp = np.sum(sum_slip_stk[istk:jstk])
             if ( sum_slip_tmp > sum_slip_eff ):
                   istk_eff = istk
                   sum_slip_eff = sum_slip_tmp
      
        # Select the interval with length Leff in the dip dir which
        # maximize the slip sum
        idip = 0
        jdip = ndip_eff
           
        idip_eff = idip
        sum_slip_eff = np.sum(sum_slip_dip[idip:jdip])
        while( jdip < self.ndip_in):
               idip+=1
               jdip+=1
               print(jdip)
               sum_slip_tmp = np.sum(sum_slip_dip[idip:jdip])
               if ( sum_slip_tmp > sum_slip_eff ):
                   idip_eff = idip
                   sum_slip_eff = sum_slip_tmp
           
        slip_eff = slip[idip_eff:idip_eff+ndip_eff, istk_eff:istk_eff+nstk_eff]
        print()
        print( f" Original slip matrix dimensions: {slip.shape} ")
        print( f" Effective slip matrix dimensions: {slip_eff.shape} ")
        print()
        return np.array([idip_eff,idip_eff+ndip_eff,istk_eff,istk_eff+nstk_eff])
       
    def set_effective_fault(self):   
               
        # Xutm, Yutm coord Km
        self.XFinMat, self.YFinMat, tmp1, tmp2 = \
        utm.from_latlon(self.LatFinMat, self.LonFinMat,33,'T')
     
        # Calculate effective dimensions
        
        eff_dim = self.effective_size()
    
        # Extracting effective part of the fault. 
        self.XFinMat = self.XFinMat[eff_dim[0]:eff_dim[1], 
                                    eff_dim[2]:eff_dim[3]]
        self.YFinMat = self.YFinMat[eff_dim[0]:eff_dim[1], 
                                    eff_dim[2]:eff_dim[3]]
        self.ZFinMat = self.ZFinMat[eff_dim[0]:eff_dim[1], 
                                    eff_dim[2]:eff_dim[3]]
        self.SlipinMat = self.SlipinMat[eff_dim[0]:eff_dim[1], 
                                        eff_dim[2]:eff_dim[3]]
        self.RiseTinMat = self.RiseTinMat[eff_dim[0]:eff_dim[1], 
                                       eff_dim[2]:eff_dim[3]]
        self.RupTinMat = self.RupTinMat[eff_dim[0]:eff_dim[1], 
                                      eff_dim[2]:eff_dim[3]]
       
        self.ndip_in, self.nstk_in = self.XFinMat.shape
        # Calculate length and output number of points in strike an dip direction        
        self.stk_len_in = round((self.nstk_in-1)*self.dstk_in)
        self.dip_len_in = round((self.ndip_in-1)*self.ddip_in)
       
        self.XFin3D = self.XFinMat.flatten(order='F').transpose()
        self.YFin3D = self.YFinMat.flatten(order='F').transpose()
        self.ZFin3D = self.ZFinMat.flatten(order='F').transpose()
        self.Slipin3D = self.SlipinMat.flatten(order='F').transpose()
        self.RiseTin3D = self.RiseTinMat.flatten(order='F').transpose()
        self.RupTin3D = self.RiseTinMat.flatten(order='F').transpose()
       
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

        # Calculate output fault's size
        self.nstk = int(self.stk_len_in/self.dhF)+1
        self.ndip = int(self.dip_len_in/self.dhF)+1
       
        # Define input and output strike and dip vectors
        self.dipinVec = np.linspace(0, self.dip_len_in, self.ndip_in)
        self.stkinVec = np.linspace(0, self.stk_len_in, self.nstk_in)
        
        self.stkVec = np.linspace(0, self.stk_len_in, self.nstk)
        self.dipVec = np.linspace(0, self.dip_len_in, self.ndip)
       
        stk_len = round((self.nstk-1)*self.dstk)
        dip_len = round((self.ndip-1)*self.ddip)
      
        print(" Output Fault Dimensions:")
        print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " 
              %(stk_len/1000, self.nstk, self.dstk/1000) )
        print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" 
              %(dip_len/1000, self.ndip, self.ddip/1000) )
        
  
     
    # *************************************************************************
    # *                                                                       *
    # *              Interpolation methods section                            *
    # *                                                                       *
    # *************************************************************************
    def interpolate_xyz_coords(self):
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
        self.n_sub_faults = self.nstk*self.ndip
        self.fcoor = np.array((XF3D,YF3D,ZF3D)).transpose()
        print(f" Subfaults: {len(self.XF3D)}")
        
    def interpolate_slip(self):
        # Slip Interpolation
        SlipF = scipy.interpolate.interp2d(self.stkinVec, self.dipinVec, 
                                     self.SlipinMat, kind = "cubic")
        self.SlipMat = SlipF(self.stkVec, self.dipVec)
 
    def interpolate_rise_time(self):
        # Slip Interpolation
        riset_fun = scipy.interpolate.interp2d(self.stkinVec, self.dipinVec, 
                                     self.RiseTinMat, kind = "cubic")
        self.riset_mat = riset_fun(self.stkVec, self.dipVec)
        
    
    def interpolate_rupt_time(self):
        # Slip Interpolation
        rupt_time_fun = scipy.interpolate.interp2d(self.stkinVec, self.dipinVec, 
                                     self.RupTinMat, kind = "cubic")
        self.rupt_time = rupt_time_fun(self.stkVec, self.dipVec)
 
    # *************************************************************************
    # *                                                                       *
    # *         Triangulatiopn methods section                                *
    # *                                                                       *
    # *************************************************************************
    def triangulate_fault(self):
        print()
        print(" Fault Triangulation")
        index_fault = np.arange(0,self.XF3D.size).reshape((self.ndip,self.nstk)
                                                          ,order='F')
        ntri  = (self.nstk-1)*(self.ndip-1)*2
        self.tri   = np.zeros([ntri,3],dtype=int)
        #xy_3D = np.array((self.XF3D,self.YF3D)).transpose()

        # Delaunay triangulation
        # tri = Delaunay(xy_3D).simplices
        self.ntri = int(self.tri.size/3)
        jtri = -1
        for istk in range (0,self.nstk-1):
            for idip in range (0,self.ndip-1):
                jtri += 1
                self.tri[jtri,0] = index_fault[idip,istk]
                self.tri[jtri,1] = index_fault[idip,istk+1]
                self.tri[jtri,2] = index_fault[idip+1,istk+1]
                jtri += 1
                self.tri[jtri,0] = index_fault[idip,istk]
                self.tri[jtri,1] = index_fault[idip+1,istk+1]
                self.tri[jtri,2] = index_fault[idip+1,istk]
        
        print(f" Number of facets: {ntri}")
        self.trib_marker = np.ones(self.ntri,)
        # Calculate unitary normal, strike and dip vector at each facet
        self.univector = np.zeros((self.ntri,9))
        # Vector normal to earth surface
        nsurf = np.array([0,0,-1])

        for itri in range(0,self.ntri):
            iv0 = self.tri[itri,0]
            iv1 = self.tri[itri,1]
            iv2 = self.tri[itri,2]
            v0 = np.array([ self.XF3D[iv0], self.YF3D[iv0], self.ZF3D[iv0]])
            v1 = np.array([ self.XF3D[iv1], self.YF3D[iv1], self.ZF3D[iv1]])
            v2 = np.array([ self.XF3D[iv2], self.YF3D[iv2], self.ZF3D[iv2]])
            self.vec_normal = np.cross(v1-v0,v2-v0)
            self.vec_normal = self.vec_normal/np.linalg.norm(self.vec_normal)
            self.vec_strike = np.cross(self.vec_normal,nsurf)
            self.vec_strike = self.vec_strike/np.linalg.norm(self.vec_strike)
            self.vec_dip = np.cross(self.vec_strike,self.vec_normal)
            self.vec_dip = self.vec_dip/np.linalg.norm(self.vec_dip)
            
            self.univector[itri,0:3] = self.vec_normal
            self.univector[itri,3:6] = self.vec_strike
            self.univector[itri,6:9] = self.vec_dip

    def add_nodes_above_below(self):
    # Add nodes above an below to the fault
        x_above = self.XF3D + self.vec_normal[0]*self.dhF
        y_above = self.YF3D + self.vec_normal[1]*self.dhF
        z_above = self.ZF3D + self.vec_normal[2]*self.dhF
        x_below = self.XF3D - self.vec_normal[0]*self.dhF
        y_below = self.YF3D - self.vec_normal[1]*self.dhF
        z_below = self.ZF3D - self.vec_normal[2]*self.dhF
        print()
        print(f" {len(x_above)} nodes added above the fault ")
        print(f" {len(x_below)} nodes added below the fault ")
        print()
        self.XF3D_add = np.concatenate((x_above, x_below), axis=None)
        self.YF3D_add = np.concatenate((y_above, y_below), axis=None)
        self.ZF3D_add = np.concatenate((z_above, z_below), axis=None)
     
    # *************************************************************************
    # *                                                                       *
    # *                 Plots methods section                                 *
    # *                                                                       *
    # *************************************************************************
    def plot_xyz_slipin(self):
        io.renderers.default='svg'
            
        colorbar=dict(lenmode='fraction', len=0.75, thickness=20, 
                      bordercolor="black", title="<b> slip(m) </b>", x=0.2)
        
        data = [go.Surface(x=self.XFinMat/1000, y=self.YFinMat/1000,
                           z=self.ZFinMat/1000, surfacecolor=self.SlipinMat,
                           colorscale=px.colors.sequential.Viridis, 
                           colorbar=colorbar, showscale=True,
                           lighting=dict(ambient=0.7))]
        
        tickfont = dict(color="black", size=16, family="Arial Black")
        
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.5, y=-2.0, z=1.2))
        
        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[360,385],
                     tickfont=tickfont, showbackground=True)
        yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[4670,4720], 
                     tickfont=tickfont, showbackground=True)
        zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[-18,2], tickfont=tickfont,
                     showbackground=True)
  
        margin = dict(r=30, l=10, b=30, t=20)

        title = dict(text="<b>Input Slip </b>", font_family="Arial Blak", 
                     font_color="black", x=0.5, y=0.85)
        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')
        
        layout = go.Layout(scene = scene, margin=margin, width=800, height=350, 
                           title=title)
        
        fig = go.Figure(data=data, layout=layout)
        fig.show()
        
    def plot_xyz_slip(self):
        io.renderers.default='svg'
            
        colorbar=dict(lenmode='fraction', len=0.75, thickness=20, bordercolor="black",
                      title="<b> slip(m) </b>", x=0.2)
        
        data = [go.Surface(x=self.XFMat/1000, y=self.YFMat/1000, z=self.ZFMat/1000, 
                          surfacecolor=self.SlipMat,
                          colorscale=px.colors.sequential.Viridis, colorbar=colorbar, showscale=True,
                          lighting=dict(ambient=0.7))]
           
        tickfont = dict(color="black", size=16, family="Arial Black")
      
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.5, y=-2.0, z=1.2))
        
        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[360,385],
                     tickfont=tickfont, showbackground=True)
        yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[4670,4720], 
                     tickfont=tickfont, showbackground=True)
        zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[-18,2], tickfont=tickfont,
                     showbackground=True)

        margin = dict(r=30, l=10, b=30, t=20)

        title = dict(text="<b>Interpolated Slip </b>", font_family="Arial Blak", 
                     font_color="black", x=0.5, y=0.85)
  
        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')
        
        layout = go.Layout(scene = scene, margin=margin, width=800, height=350, 
                           title=title)
        
        fig = go.Figure(data=data, layout=layout)
        fig.show()
    
    def plot_xyz_model_slip(self):
        io.renderers.default='svg'
                
        colorbar=dict(lenmode='fraction', len=0.75, thickness=20, bordercolor="black",
                          title="<b> slip(m) </b>", x=0.2)
            
        data = [go.Surface(x=self.XFMat/1000, y=self.YFMat/1000, z=self.ZFMat/1000, 
                          surfacecolor=self.SlipMat,
                          colorscale="Viridis", colorbar=colorbar, showscale=True,
                          lighting=dict(ambient=0.7))]
               
        tickfont = dict(color="black", size=16, family="Arial Black")
          
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                          eye=dict(x=-0.5, y=-2.0, z=1.5))
            
        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=10, range=[300,420],
                     tickfont=tickfont, showbackground=True)
        yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=10, range=[4610,4750], 
                     tickfont=tickfont, showbackground=True)
        zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=6, range=[-60,2], tickfont=tickfont,
                     showbackground=True)

        margin = dict(r=30, l=50, b=20, t=20)

        title = dict(text="<b>Interpolated Slip </b>", font_family="Arial Blak", 
                         font_color="black", x=0.5, y=0.85)
      
        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')
            
        layout = go.Layout(scene = scene, margin=margin, width=800, height=350, 
                           title=title)
            
        fig = go.Figure(data=data, layout=layout)
        fig.show()
        
    def compare_xyz_slip(self):
        io.renderers.default='svg'
            
        colorbar=dict(lenmode='fraction', len=0.9, thickness=10, 
                      bordercolor="black", title="<b> slip(m) </b>")
        
        data_a = go.Surface(x=self.XFinMat/1000, y=self.YFinMat/1000, z=self.ZFinMat/1000, 
                            surfacecolor=self.SlipinMat,
                            colorbar=colorbar, 
                            showscale=False, lighting=dict(ambient=0.9))
           
        data_b = go.Surface(x=self.XFMat/1000, y=self.YFMat/1000, z=self.ZFMat/1000, 
                            surfacecolor=self.SlipMat, colorbar=colorbar,
                            showscale=True, lighting=dict(ambient=0.9))
        
        tickfont = dict(color="black", size=12, family="Arial Black")
      
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.5, y=-2.0, z=1.2))
        
        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, gridcolor="white",
                     showticklabels=True, nticks=6, range=[360,385],
                     tickfont=tickfont, showbackground=True)
        yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=6, range=[4670,4720], 
                     tickfont=tickfont, showbackground=True)
        zaxis = dict(title="<b> z (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=6, range=[-18,2], tickfont=tickfont,
                     showbackground=True)
        margin = dict(r=5, l=5, b=10, t=20)
        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')
               
        fig = tls.make_subplots(rows=1, cols=2,specs=[[{'is_3d': True},
                                                       {'is_3d': True} ]],
                                horizontal_spacing = 0,
                                subplot_titles=["<b>Input Slip </b>",
                                                "<b>Interpolated Slip </b>"])
        fig.append_trace(data_a, 1, 1)
        fig.append_trace(data_b, 1, 2)
        fig.update_scenes(scene, row=1, col=2)
        fig.update_scenes(scene, row=1, col=1)
        fig.update_layout(height=310, width=800, margin=margin)
        fig.show()
    
    def plot_triangulation(self):
        tickfont = dict(color="black", size=19, family="Arial Black")
      
        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.0, y=-1.0, z=0.7))
        label_z = dict(text="<b>z (Km)</b>", font_family="Arial Blak", 
                         font_color="black", font_size=26)
        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, gridcolor="white",
                     showticklabels=True, nticks=2, range=[360,390],
                     tickfont=tickfont, showbackground=True)
        yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=2, range=[4670,4710], 
                     tickfont=tickfont, showbackground=True)
        zaxis = dict(title=label_z, showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=4, range=[-15,-1], tickfont=tickfont,
                     showbackground=True)
        margin = dict(r=10, l=20, b=20, t=20)
        
        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')
        
        title = dict(text="<b>Fault Triangulation</b>", font_family="Arial Blak", 
                         font_color="black", font_size=24, x=0.5, y=0.97)
        
        layout = go.Layout(scene = scene, margin=margin, width=1600, 
                           height=1200, title=title)
     
        fig = ff.create_trisurf(x=self.XF3D/1000, y=self.YF3D/1000, z=self.ZF3D/1000, 
                                colormap=[(1.0, 1.0, 0.6), (0.95, 0.95, 0.6)],
                                simplices=self.tri, plot_edges=True)
        fig.update_layout(layout)
        fig.update(layout_coloraxis_showscale=False)
           
        fig.show()   
   
  
    # *************************************************************************
    # *                                                                       *
    # *                     save files methods section                        *                     
    # *                                                                       *
    # *************************************************************************
    def save_fault(self,dir):
        print()
        out_file = dir + self.name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(self, object_file)
        object_file.close()
        print(f" Fault object saved in: {out_file} ")

    def write_univector(self,dir):
        print()
        # Write .vector file
        fvector_header = "%d" %(self.ntri)
        fvector = dir + self.name + ".vector"
        with open(fvector,'wb') as f:
            np.savetxt(f, self.univector,header=fvector_header, 
                       comments=' ',fmt='%10.6f')
        f.close()

        print(f" vector file saved in: {fvector} ")
        
    def write_fcoor(self,dir,id,nt,dt):
        print()
        # Write fcoor file
        fcoorHeader = "%d  %d %4.2f " %(self.n_sub_faults, nt, dt)

        fcoorName = dir + self.name + "_ID_" + id +".in"

        with open(fcoorName,'wb') as f:
            np.savetxt(f, self.fcoor, header=fcoorHeader, comments=' ',fmt = '%9.4f')

        print(f" Coordinates file saved in: {fcoorName} " )
        print(f" Number of subfaults: {self.n_sub_faults} " )
        print(f" Number of time steps: {nt} " )
        print(f" time step: {dt} " )
        
        f.close()

               