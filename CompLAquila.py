#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Compare velocity for Amatrice 2017
#
# John Diaz October 2022

# To Do List:
# Add Data
    
import os 
import warnings

import obspy
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
plt.close('all') 
os.system('clear')

# DG folder
DGFolder = "../DGrun/LAquilaCirella03_dhF500m_ID_2/"

# Data folder
DataFolder = "../Data/Streams/"

#Stats = [ "AQK","AQU","AQV","AQA","AQG","GSA","MTR","FMG","ANT","AVZ","CSO1","LSS","SUL"]
Stats  = [ "AQK","AQU","AQV","AQA","AQG","GSA","MTR","FMG","ANT","AVZ","CSO1","LSS","SUL"]


# Filtering parameters
lowcut  = 0.5                              # low pass cut frecuency
highcut = 0.05                             # cut high pass frecuency

# initial time for plotting components (s)
tiy = 30;  tix = 80;  tiz = 130;
xmin = 0; xmax = 205
ymin = -40; ymax = 2

print("  ")
print(" START PROGRAM ")
print("  ")

di=0

for stat in Stats:
    # DG Synthetic Data Reading
    DGfile = DGFolder + stat + "_DGVEL.pickle"
    print(f" Reading DG file {DGfile}")
    # Read the obspy Stream with three components of Synthetic velocity
    DGst = obspy.read(DGfile)    
    DGvx = DGst[0].copy()
    DGvy = DGst[1].copy()
    DGvz = DGst[2].copy()
    
    # Integrate to get the Synthetic displacement
    DGdx = DGst[0].integrate(method='spline')
    DGdy = DGst[1].integrate(method='spline')
    DGdz = DGst[2].integrate(method='spline')
    
    Datafile = DataFolder + stat + "_VEL.pickle"
    print(f" Reading Data file {Datafile}")
    # Reading the Obspy Stream with three components of Data velocity
    Datast = obspy.read(Datafile)
    # Take the starttime of the Synthetic which is the event time
    t = DGvx.stats.starttime
    # Take the data from the event time plus 40 seconds (Synthetic simulation time)
    Datavx = Datast[0].copy().slice(t,t+40)
    Datavy = Datast[1].copy().slice(t,t+40)
    Datavz = Datast[2].copy().slice(t,t+40)
    # Integrate to get the Data displacement
    Datadx = Datast[0].integrate(method='spline').slice(t,t+40)
    Datady = Datast[1].integrate(method='spline').slice(t,t+40)
    Datadz = Datast[2].integrate(method='spline').slice(t,t+40)
        
    OBStimex = Datavx.times(reftime=DGvx.stats.starttime)
    SYNtimex = DGvx.times()
  
    # Bandpass filtering of velocity
    DGvx.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    DGvy.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    DGvz.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datavx.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datavy.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datavz.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
  
    # Calculate the maximum value of the velocity 
    maxv = max([ abs(Datavx.max()), abs(Datavy.max()), abs(Datavz.max()) ])
    
    # Normalice Synthetic and Data with the maximum velocity of the Data 
    Datavx.normalize(norm=maxv)
    DGvx.normalize(norm=maxv)
    Datavy.normalize(norm=maxv)
    DGvy.normalize(norm=maxv)
    Datavz.normalize(norm=maxv)
    DGvz.normalize(norm=maxv)
    
    
    print(DGvy.stats.starttime)
    # Get the time for Synthetic and Data
    OBStimex = Datavx.times(reftime=DGvy.stats.starttime)
    SYNtimex = DGvx.times()
    OBStimey = Datavy.times(reftime=DGvy.stats.starttime)
    SYNtimey = DGvy.times()
    OBStimez = Datavz.times(reftime=DGvy.stats.starttime)
    SYNtimez = DGvz.times()
   
    # Velocity Comparison
    di -= 2.2
    fig = plt.figure(1)
    plt.title(" LAquila ")    
    plt.plot(OBStimey+tiy,Datavy.data+di,color='k')
    plt.plot(SYNtimey+tiy,DGvy.data+di,color='r')
    plt.plot(OBStimex+tix,Datavx.data+di,color='k')
    plt.plot(SYNtimex+tix,DGvx.data+di,color='r')
    plt.plot(OBStimez+tiz,Datavz.data+di,color='k')
    plt.plot(SYNtimez+tiz,DGvz.data+di,color='r')
    plt.text(4, di, stat,fontsize=10,fontweight='bold')
    plt.text(2, 0, 'Station',fontsize=8,fontweight='bold')
    plt.text(40, 0, ' Vy',fontsize=8,fontweight='bold')
    plt.text(90, 0, ' Vx',fontsize=8,fontweight='bold')
    plt.text(145, 0, ' Vz',fontsize=8,fontweight='bold')
    plt.text(175, 0.8, ' max vel ',fontsize=8,fontweight='bold')
    plt.text(177, 0, ' (cm/s)',fontsize=8,fontweight='bold')
    plt.text(177, di, '{:5.2f}'.format(maxv),fontsize=10,fontweight='bold')
    plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
    plt.legend(['Data','DG'],loc='lower center')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    
    # Bandpass filtering of displacement
    DGdx.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    DGdy.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    DGdz.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datadx.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datady.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
    Datadz.filter("bandpass",freqmin=highcut,freqmax=lowcut,zerophase=True)
 
    # Calculate the maximum value of the velocity 
    maxd = max([ abs(Datadx.max()), abs(Datady.max()), abs(Datadz.max()) ])
    
    # Normalice Synthetic and Data with the maximum velocity of the Data 
    Datadx.normalize(norm=maxd)
    DGdx.normalize(norm=maxd)
    Datady.normalize(norm=maxd)
    DGdy.normalize(norm=maxd)
    Datadz.normalize(norm=maxd)
    DGdz.normalize(norm=maxd)
    
    # Velocity Comparison
  
    fig = plt.figure(2)
    plt.title(" LAquila ")    
    plt.plot(OBStimey+tiy,Datady.data+di,color='k')
    plt.plot(SYNtimey+tiy,DGdy.data+di,color='r')
    plt.plot(OBStimex+tix,Datadx.data+di,color='k')
    plt.plot(SYNtimex+tix,DGdx.data+di,color='r')
    plt.plot(OBStimez+tiz,Datadz.data+di,color='k')
    plt.plot(SYNtimez+tiz,DGdz.data+di,color='r')
    plt.text(4, di, stat,fontsize=10,fontweight='bold')
    plt.text(2, 0, 'Station',fontsize=8,fontweight='bold')
    plt.text(40, 0, ' Vy',fontsize=8,fontweight='bold')
    plt.text(90, 0, ' Vx',fontsize=8,fontweight='bold')
    plt.text(145, 0, ' Vz',fontsize=8,fontweight='bold')
    plt.text(175, 0.8, ' max disp ',fontsize=8,fontweight='bold')
    plt.text(177, 0, ' (cm)',fontsize=8,fontweight='bold')
    plt.text(177, di, '{:5.2f}'.format(maxd),fontsize=10,fontweight='bold')
    plt.tick_params(left = False, right = False , labelleft = False ,
               labelbottom = False, bottom = False)
    plt.legend(['Data','DG'],loc='lower center')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
  

   








