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

class Source:
    def __init__(self, Fault, Mw, asp_size, kc, rake):
        self.name = Fault.name
        self.dhF = Fault.dhF 
        self.ndip = Fault.ndip
        self.nstk = Fault.nstk
        self.target_Mw = Mw
        self.asp_size = asp_size
        self.kc = kc
        self.rake = rake
        self.slip = Fault.SlipMat
        self.stk = Fault.stkVec*1000
        self.dip = Fault.dipVec*1000
        self.XF3D = Fault.XF3D
        self.YF3D = Fault.YF3D
        self.ZF3D = Fault.ZF3D

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
    
    def effective_size(self):
        # Cumulative slip in strike and dip direction
        sum_slip_stk = np.sum(self.slip, axis=0)
        sum_slip_dip = np.sum(self.slip, axis=1)
        
        acorr_stk = sm.tsa.acf(sum_slip_stk, nlags = len(sum_slip_stk)-1)
        acorr_dip = sm.tsa.acf(sum_slip_dip, nlags = len(sum_slip_dip)-1)
        
        Wacf_stk = integrate.simpson(acorr_stk, self.stk)/max(acorr_stk)
        Wacf_dip = integrate.simpson(acorr_dip, self.dip)/max(acorr_dip)
        
        pass
        
    
    def set_mu_1d(self, vs, rho, h):
        z = np.zeros(h.size)
        z[0] = 0
        for ih in range(0,h.size-1):
            z[ih+1] = z[ih]-h[ih]          
       
        self.mu = np.zeros(self.ZF3D.size)
        
        for iz in range(self.ZF3D.size):
            #print(iz)
            for ih in range(0,z.size-1):
                if ( self.ZF3D[iz] <= z[ih] and self.ZF3D[iz] > z[ih+1]):
                    self.mu[iz] = vs[ih]*vs[ih]*rho[ih]*1e6
     
        
    def adjust_Mw(self):
                 
        M0_target = pow(10,(self.target_Mw+10.75)/1.5)
   
        moment = 0
        slip =self.slip.flatten(order='F')
        area =pow(self.dhF*1000,2)
        
        for iz in range(self.mu.size):
           moment += area*slip[iz]*self.mu[iz]
            
        moment = moment*1e-7   # from N-m to dyn-cm    
        self.original_Mw = -10.75 + 1.5*log10(moment)
        print(f" Original Mw: %5.2f" % (self.original_Mw))
        self.adjust_slip = self.slip*M0_target/moment
        slip =self.adjust_slip.flatten(order='F')
        for iz in range(self.mu.size):
           moment += area*slip[iz]*self.mu[iz]
            
        moment = moment*1e-7   # from N-m to dyn-cm   
        self.adjust_Mw = -10.75 + 1.5*log10(moment)
        
        print(f" Adjusted Mw: %5.2f" % (self.adjust_Mw))
        
        
        
        pass




