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

import pickle

class Fault():
    def __init__(self,name, dhF):
        self.name = name
        self.dhF = dhF

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
        self.fcoor = np.array((XF3D,YF3D,ZF3D)).transpose()
        
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
 
    def triangulate_fault(self):
        index_fault = np.arange(0,self.XF3D.size).reshape((self.ndip,self.nstk)
                                                          ,order='F')
        ntri  = (self.nstk-1)*(self.ndip-1)*2
        tri   = np.zeros([ntri,3],dtype=int)
        #xy_3D = np.array((self.XF3D,self.YF3D)).transpose()

        # Delaunay triangulation
        # tri = Delaunay(xy_3D).simplices
        self.ntri = int(tri.size/3)
        jtri = -1
        for istk in range (0,self.nstk-1):
            for idip in range (0,self.ndip-1):
                jtri += 1
                tri[jtri,0] = index_fault[idip,istk]
                tri[jtri,1] = index_fault[idip,istk+1]
                tri[jtri,2] = index_fault[idip+1,istk+1]
                jtri += 1
                tri[jtri,0] = index_fault[idip,istk]
                tri[jtri,1] = index_fault[idip+1,istk+1]
                tri[jtri,2] = index_fault[idip+1,istk]

        self.trib_marker = np.ones(self.ntri,)
        # Calculate unitary normal, strike and dip vector at each facet
        self.univector = np.zeros((self.ntri,9))
        # Vector normal to earth surface
        nsurf = np.array([0,0,-1])

        for itri in range(0,self.ntri):
            iv0 = tri[itri,0]
            iv1 = tri[itri,1]
            iv2 = tri[itri,2]
            v0 = np.array([ self.XF3D[iv0], self.YF3D[iv0], self.ZF3D[iv0]])
            v1 = np.array([ self.XF3D[iv1], self.YF3D[iv1], self.ZF3D[iv1]])
            v2 = np.array([ self.XF3D[iv2], self.YF3D[iv2], self.ZF3D[iv2]])
            self.vec_normal = np.cross(v1-v0,v2-v0)
            self.vec_strike = np.cross(self.vec_normal,nsurf)
            self.vec_vdip = np.cross(self.vec_strike,self.vec_normal)
            self.univector[itri,0:3] \
            = self.vec_normal/np.linalg.norm(self.vec_normal)
            self.univector[itri,3:6] \
            = self.vec_strike/np.linalg.norm(self.vec_strike)
            self.univector[itri,6:9] \
            = self.vec_dip/np.linalg.norm(self.vec_dip)

    def add_nodes_above_below(self):
    # Add nodes above an below to the fault
        x_above = self.XF3D + self.vec_normal[0]*self.dhF
        y_above = self.YF3D + self.vec_normal[1]*self.dhF
        z_above = self.ZF3D + self.vec_normal[2]*self.dhF
        x_below = self.XF3D - self.vec_normal[0]*self.dhF
        y_below = self.YF3D - self.vec_normal[1]*self.dhF
        z_below = self.ZF3D - self.vec_normal[3]*self.dhF
        self.XF3D_add =np.concatenate((x_above, x_below), axis=None)
        self.YF3D_add =np.concatenate((y_above, y_below), axis=None)
        self.ZF3D_add =np.concatenate((z_above, z_below), axis=None)
     
    def plot_xyz_slipin(self):
        io.renderers.default='svg'
            
        colorbar=dict(lenmode='fraction', len=0.75, thickness=20, bordercolor="black",
                      title="<b> slip(m) </b>", x=0.2)
        
        data = [go.Surface(x=self.XFinMat, y=self.YFinMat, z=self.ZFinMat, 
                          surfacecolor=self.SlipinMat,
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
        
        data = [go.Surface(x=self.XFMat, y=self.YFMat, z=self.ZFMat, 
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
            
        data = [go.Surface(x=self.XFMat, y=self.YFMat, z=self.ZFMat, 
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
        
        data_a = go.Surface(x=self.XFinMat, y=self.YFinMat, z=self.ZFinMat, 
                            surfacecolor=self.SlipinMat,
                            colorbar=colorbar, 
                            showscale=False, lighting=dict(ambient=0.9))
           
        data_b = go.Surface(x=self.XFMat, y=self.YFMat, z=self.ZFMat, 
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
        
        
    def save_fault(self,dir):
        out_file = dir + self.name + "_dhF" + str(self.dhF*1000) + ".pickle"

        print()
        object_file = open(out_file, 'wb')
        pickle.dump(self, object_file)
        object_file.close()
        print(f" Fault object saved in: {out_file} ")

    def write_univector(self,dir):
        print()
        # Write .vector file
        fvector_header = "%d" %(self.ntri)
        fvector = dir + self.name + "_dhF" + str(self.dhF*1000) + ".vector"
        with open(fvector,'wb') as f:
            np.savetxt(f, self.univector,header=fvector_header, 
                       comments=' ',fmt='%10.6f')
        f.close()

        print(f" vector file saved in: {fvector} ")
        