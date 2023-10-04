#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#       
# Source class
# 

import numpy as np
import statsmodels.api as sm
from scipy import integrate
from math import log10
from scipy.interpolate import splev, splrep

import plotly.graph_objs as go
import plotly.io as io
import plotly.express as px

import plotly.tools as tls

import matplotlib.pyplot as plt
import Fault


class Source():
    def __init__(self, fault_object, Mw, asp_size, kc, rake):
        
        self.target_Mw = Mw
        self.asp_size = asp_size
        self.kc = kc
        self.rake = rake
        # copy all atributes form Fault object.
        for attribute_key in fault_object.__dict__.keys():
            self.__dict__[attribute_key] = fault_object.__dict__[attribute_key]
     
 
    def __str__(self):
        source = " Fault info \n"
        source += " name: " + str(self.name) + "\n"
        source += " dhF : " + str(self.dhF) + "\n"
        source += " ndip : " + str(self.ndip) + "\n"
        source += " nstk : " + str(self.nstk) + "\n \n"
        source += " Source info \n"
        source += " asp_size : " + str(self.asp_size) + "\n"
        source += " kc : " + str(self.kc) + "\n"
        source += " Target Mw: " + str(self.target_Mw) + "\n"
        return source
    
    
        
    
    def set_mu_1d(self, vs, rho, h):
        z = np.zeros(h.size)
        z[0] = 0
        for ih in range(0,h.size-1):
            z[ih+1] = z[ih]-h[ih]          
       
        self.mu = np.zeros(self.ZF3D.size)
        
        for iz in range(self.ZF3D.size):
            for ih in range(0,z.size-1):
                if ( self.ZF3D[iz] <= z[ih] and self.ZF3D[iz] > z[ih+1]):
                    self.mu[iz] = vs[ih]*vs[ih]*rho[ih]*1e6
     
        
    def adjust_Mw(self):
                 
        M0_target = pow(10,(self.target_Mw+10.75)/1.5)
   
        moment = 0
        slip = self.SlipMat.flatten(order='F')
        area = pow(self.dhF*1000,2)
        
        slip_umbral = 0.05*np.mean(slip)
      
        for iz in range(self.mu.size):
            if (slip[iz] > slip_umbral):
                moment += area*slip[iz]*self.mu[iz]
            
        moment = moment*1e-7   # from N-m to dyn-cm    
        self.original_Mw = -10.75 + 1.5*log10(moment)
        print(" Original Mw: %5.2f" % (self.original_Mw))
        
        self.adjust_slip = self.SlipMat*M0_target/moment
        slip = self.adjust_slip.flatten(order='F')
        
        for iz in range(self.mu.size):
            if (slip[iz] > slip_umbral):
                moment += area*slip[iz]*self.mu[iz]
            
        moment = moment*1e-7   # from N-m to dyn-cm   
        self.adjust_Mw = -10.75 + 1.5*log10(moment)
        
        print(" Adjusted Mw: %5.2f" % (self.adjust_Mw))
        
    def compare_xyz_slip(self):
        io.renderers.default='svg'
            
        colorbar=dict(lenmode='fraction', len=0.9, thickness=10, 
                      bordercolor="black", title="<b> slip(m) </b>")
        
        data_a = go.Surface(x=self.XFMat, y=self.YFMat, z=self.ZFMat, 
                            surfacecolor=self.SlipMat,
                            colorbar=colorbar, 
                            showscale=False, lighting=dict(ambient=0.9))
           
        data_b = go.Surface(x=self.XFMat, y=self.YFMat, z=self.ZFMat, 
                            surfacecolor=self.adjust_slip, colorbar=colorbar,
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

    def compare_slip_contour_2d(self):

        tickfont = dict(color="black", size=12, family="Arial Black")
        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont, 
                     autorange="reversed")
        figa = go.Figure(data=go.Heatmap( z=self.SlipMat, x= self.stkVec, 
                                         y=self.dipVec, hoverongaps = False))
        
        figa.update_layout(xaxis=xaxis, yaxis=yaxis)
        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours=dict(coloring='lines',showlabels=True)
        figc = go.Figure(data=go.Contour(z=self.rupt_time, x=self.stkVec, 
                                         y=self.dipVec, contours=contours, 
                                         line_width=2, showscale=False,
                                         colorscale=colorscale))
        figa.add_trace(figc.data[0])
        figb = go.Figure(data=go.Heatmap( z=self.adjust_slip, x= self.stkVec, 
                                         y=self.dipVec, hoverongaps = False))
        scene = dict(xaxis=xaxis, yaxis=yaxis)
        figb.update_layout(xaxis=xaxis, yaxis=yaxis)
        figb.add_trace(figc.data[0])
        fig = tls.make_subplots(rows=1, cols=2, horizontal_spacing = 0.1,
                                subplot_titles=["<b> Original Slip </b>",
                                                "<b> Adjusted Slip </b>"])
        fig.append_trace(figa.data[0], 1, 1)
        fig.add_trace(figc.data[0], 1, 1)
        fig.append_trace(figb.data[0], 1, 2)
        fig.add_trace(figc.data[0], 1, 2)
        fig["layout"]["xaxis1"].update(title="<b> strike (Km) </b>", 
                                       tickfont=tickfont) 
        fig["layout"]["yaxis1"].update(title="<b> dip (Km) </b>", 
                                       tickfont=tickfont, autorange="reversed")
        fig["layout"]["xaxis2"].update(title="<b> strike (Km) </b>", 
                                       tickfont=tickfont) 
        fig["layout"]["yaxis2"].update(title="<b> dip (Km) </b>", 
                                       tickfont=tickfont, autorange="reversed")
        margin = dict(r=5, l=5, b=10, t=20)
        fig.update_layout(height=310, width=800, margin=margin)
      
        fig.show()
        
        
    def original_slip_psd(self):
    
        # Taking the fourier transform centered at the origin
        slip_fft = np.fft.fft2(self.adjust_slip)
        slip_fft_shift = np.fft.fftshift(slip_fft)
        # Get power spectral density
        slip_psd_amplitude = np.abs(slip_fft)**2
        slip_shift_psd_amplitude = np.abs(slip_fft_shift)**2
        slip_psd_amplitude_log = np.log(slip_psd_amplitude)
        slip_shift_psd_amplitude_log = np.log(slip_shift_psd_amplitude)
        self.npad = 512
        # Prepare wavenumbers centered at the origin
        kstk = (2*np.pi/(self.npad*self.dhF))*np.arange(-self.npad//2,self.npad//2-1, 1.)
        kdip = (2*np.pi/(self.npad*self.dhF))*np.arange(-self.npad//2, self.npad//2-1, 1.)
        kave = 0.5 * (kstk[1:] + kdip[:-1])
        
        fig = go.Figure(data=go.Heatmap( z=slip_shift_psd_amplitude_log, x=kstk,
                                        y=kdip,hoverongaps = False, showscale=True))
        #fig = px.scatter(slip_psd_log, log_x=True)
        fig.show()
        pass
    
    def low_wavenumber_scenario(self):
        npad = 512
        
        pass
        
