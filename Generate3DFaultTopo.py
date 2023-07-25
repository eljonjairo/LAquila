
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Compare velocity for Amatrice 2017
#
# John Diaz January 2023

# To Do List:
# Load Stations
# Gray Scale topo?


import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from scipy.spatial import distance
import pickle
from pathlib import Path
import pygmt
import utm

import matplotlib.colors
from matplotlib import cm
from scipy.io import loadmat


warnings.filterwarnings("ignore")
os.system('clear')

dhF = 0.5  # Output subfaults size in Km

inDir  = Path('../Input/')
outFaultDir = Path('../Outputs/3DFaults/')
outTopoDir = Path('../Outputs/Topo/')
name = "LAquilaCirella03_dhF"  # Output files name
# Input matlab file name
inFmat = 's2009LAQUIL03CIRE.mat'

# Lon-Lat limits for Topo Download
Latmin = 41.0
Latmax = 44.0
Lonmin = 12.5
Lonmax = 14.5

# # X and Y UTM limits (Km) for Model 23 Stations
# xmin = 300.0
# xmax = 440.0
# ymin = 4560.0
# ymax = 4820.0
# zmin = -60.0

# X and Y UTM limits (Km) for Model 13 Stations
xmin = 300.0
xmax = 420.0
ymin = 4610.0
ymax = 4750.0
zmin = -60.0

# m to Km
m = 1000

print("  ")
print(" START PROGRAM ")
print("  ")

# Load Earth relief data for the entire globe and a subset region
region = [Lonmin,Lonmax,Latmin,Latmax]
grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

zTopoMat = grid.data
TopoLat = grid.lat.data
TopoLon = grid.lon.data

zmax = np.max(zTopoMat)/m

TopoLonMat, TopoLatMat = np.meshgrid(TopoLon,TopoLat)
xTopoMat,yTopoMat, tmp1, tmp2 = utm.from_latlon(TopoLatMat, TopoLonMat,33,'N')
xTopoMat = xTopoMat/m
yTopoMat = yTopoMat/m

fig = pygmt.Figure()
pygmt.makecpt(cmap="geo", series=[-10000, 10000])
fig.grdimage(grid=grid, region=region,
             frame=['WSrt+t" Topografia Italia Central"', "xa0.4", "ya0.4"])
fig.colorbar(frame=["xa2000f500+lElevación ", "y+lm"])

fig.show()

cvals  = [ 0, 100, 1000, 2000, 3000, 4000, ]
colors = ["darkgreen","green","forestgreen","yellowgreen",
          "orange","maroon","sienna","brown","white"]

norm=plt.Normalize(min(cvals),max(cvals))
tuples = list(zip(map(norm,cvals), colors))
cTopo = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

#Load Matlab Finite Fault input file
infile = inDir.joinpath(inFmat)

Fault = loadmat(infile)

Fault = Fault['s2009LAQUIL03CIRE']

ZFinMat = -Fault['geoZ'][0][0]
LatFinMat = Fault['geoLAT'][0][0]
LonFinMat = Fault['geoLON'][0][0]
ndipin, nstkin = ZFinMat.shape
SlipinMat  = Fault['slipSPL'][0][0]/100     #Slip in meters.
RiseTinMat = Fault['riseSPL'][0][0]
RupTinMat  = Fault['timeSPL'][0][0]

# Hypocenter coordinates (Km) 
hypolon = Fault['evLON'][0][0][0][0]
hypolat = Fault['evLAT'][0][0][0][0]
hypoz = -Fault['evDPT'][0][0][0][0]

XFinMat,YFinMat, tmp1, tmp2 = utm.from_latlon(LatFinMat, LonFinMat,33,'T')
XFinMat = XFinMat/m
YFinMat = YFinMat/m

fig = plt.figure(figsize = (10,10))
ax = fig.subplots(1,1)
mp=ax.pcolormesh(xTopoMat,yTopoMat,zTopoMat*m, cmap=cTopo)
plt.colorbar(mp,location='bottom',label="Elevación (m)",shrink=.6)
fp=ax.pcolormesh(XFinMat,YFinMat,SlipinMat, cmap=cm.viridis)
plt.colorbar(fp,location='right',label="Slip (m)",shrink=.6)

ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_aspect('equal',adjustable='box')
ax.set_title('Topografía Italia Central Input Slip')


# Interpolation of fault plane
dd = Fault['invDzDx'][0][0]
X = Fault['geoX'][0][0]
Y = Fault['geoY'][0][0]

dstk = np.array([ XFinMat[0,1]-XFinMat[0,0], YFinMat[0,1]-YFinMat[0,0],
                  ZFinMat[0,1]-ZFinMat[0,0] ])
dstkin = np.linalg.norm(dstk)

ddip = np.array([ XFinMat[1,0]-XFinMat[0,0], YFinMat[1,0]-YFinMat[0,0],
                  ZFinMat[1,0]-ZFinMat[0,0] ])
ddipin = np.linalg.norm(ddip)

dstk = dstk*dhF
ddip = ddip*dhF

# Calculate the strike and dip unitary vetors
univec_stk = np.linalg.norm(dstk)
univec_dip = np.linalg.norm(ddip)

stk = round((nstkin-1)*dstkin)
dip = round((ndipin-1)*ddipin)

hypox, hypoy, tmp1, tmp2 = utm.from_latlon(hypolat,hypolon,33,'T')
hypox = hypox/m
hypoy = hypoy/m

print()
print(" Original Fault Dimensions:")
print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " %(stk,nstkin,dstkin) )
print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" %(dip,ndipin,ddipin) )
print(" Hypocenter Coordinates x, y and z (Km): %6.2f %6.2f %6.2f " %(hypox,hypoy,hypoz) )

dipinVec = np.linspace(0, dip, ndipin)
stkinVec = np.linspace(0, stk, nstkin)
stkinMat, dipinMat = np.meshgrid(stkinVec,dipinVec)

#interpolation
nstk = int(stk/dhF)+1
ndip = int(dip/dhF)+1
stkVec = np.linspace(0, stk, nstk)
dipVec = np.linspace(0, dip, ndip)
stkMat, dipMat = np.meshgrid(stkVec,dipVec)

# Slip Interpolation
SlipF = interpolate.interp2d(stkinVec,dipinVec,SlipinMat, kind = "linear")
SlipMat = SlipF(stkVec, dipVec)
# RiseTime Interpolation
RiseTF = interpolate.interp2d(stkinVec,dipinVec,RiseTinMat, kind = "linear")
RiseTMat = RiseTF(stkVec, dipVec)
# Rupture Time Interpolation
RupTF = interpolate.interp2d(stkinVec,dipinVec,RupTinMat, kind = "linear")
RupTMat = RupTF(stkVec, dipVec)

#Coordinates Interpolation
inivec = np.array([ XFinMat[0,0], YFinMat[0,0], ZFinMat[0,0] ])

XFMat = np.zeros((ndip,nstk))
YFMat = np.zeros((ndip,nstk))
ZFMat = np.zeros((ndip,nstk))

hypod = np.zeros((ndip,nstk))       # Hypocenter distance
hypoxy = [hypox, hypoy]

for istk in range (0,nstk):
    delta_stk = istk*dstk
    for idip in range (0,ndip ):
        delta_dip = idip*ddip
        XFMat[idip,istk] = inivec[0] + delta_stk[0] + delta_dip[0]
        YFMat[idip,istk] = inivec[1] + delta_stk[1] + delta_dip[1]
        ZFMat[idip,istk] = inivec[2] + delta_stk[2] + delta_dip[2]

        # Calculate hypocenter distance in XY coords.
        XYvec = [XFMat[idip,istk],YFMat[idip,istk]]
        hypod[idip,istk] = distance.euclidean(XYvec,hypoxy)

# Find the indexes of the hypocenter
hypoidip, hypoistk = np.where(hypod == np.min(hypod))

# Calculate the rigth z fault coords using hypocenter as reference,
# the hypocenter have to be over the fault
zmov = hypoz - ZFMat[hypoidip,hypoistk]
ZFMat = ZFMat + zmov

# From matrix to column vector following fortran 
XF3D = XFMat.flatten(order='F').transpose()
YF3D = YFMat.flatten(order='F').transpose()
ZF3D = ZFMat.flatten(order='F').transpose()

fcoor = np.array((XF3D,YF3D,ZF3D)).transpose()

fig = plt.figure(figsize = (10,10))
ax = fig.subplots(1,1)
mp=ax.pcolormesh(SlipMat, cmap=cm.viridis)
plt.colorbar(mp,location='bottom',label="Slip (m)",shrink=.6)

levels = np.linspace(0, np.max(RupTMat),9)


fig = plt.figure(constrained_layout=True)
ax1 = fig.add_subplot(121)
ax1.pcolormesh(stkinMat,dipinMat,SlipinMat)
cs=ax1.contour(stkinMat,dipinMat,RupTinMat,levels,colors=('w',),linewidths=(0.3,),origin='lower')
ax1.clabel(cs, fmt='%2.1f', colors='w', fontsize=10)
#fig.colorbar(mslip,ax=ax1,label='Slip (m)', orientation="horizontal",shrink=.95,ticks=Slipticks)
ax1.set_xlabel(' Strike (Km) ')
ax1.set_ylabel(' Dip (Km) ')
ax1.set_xlim([0,stk])
ax1.set_ylim([0,dip])
ax1.set_aspect('equal',adjustable='box')
ax1.set_title(" Input Slip ")

plt.gca().invert_yaxis()
ax2 = fig.add_subplot(122)
ax2.pcolormesh(stkMat,dipMat,SlipMat)
cs=ax2.contour(stkMat,dipMat,RupTMat,levels,colors=('w',),linewidths=(0.3,),origin='lower')
ax2.clabel(cs, fmt='%2.1f', colors='w', fontsize=10)
#fig.colorbar(mslip,ax=ax2,label='Slip (m)', orientation="horizontal",shrink=.95,ticks=Slipticks)
ax2.set_xlabel(' Strike (Km) ')
ax2.set_xlim([0,stk])
ax2.set_ylim([0,dip])
ax2.set_aspect('equal',adjustable='box')
plt.gca().invert_yaxis()
ax2.set_title(" Interpolate Slip ")

fig = plt.figure()
ax3 = fig.add_subplot(121, projection='3d')
surf = ax3.plot_surface( XFMat, YFMat, ZFMat, facecolors=cm.hsv(SlipMat), linewidth=0,
                        antialiased=False )
ax3.scatter(hypox,hypoy,hypoz,color='black', marker='*')
ax3.set_xlabel(" X (Km)")
ax3.set_ylabel(" Y (Km)")
ax3.set_xlim(xmin,xmax)
ax3.set_ylim(ymin,ymax)
ax3.set_zlim(-15.0,15.0)
ax3.set_aspect('equal',adjustable='box')
ax3.azim = -120
ax3.dist = 10
ax3.elev = 10

ax4 = fig.add_subplot(122, projection='3d')
surf = ax4.plot_surface( XFMat, YFMat, ZFMat, facecolors=cm.hsv(SlipMat), linewidth=0,
                        antialiased=False )
ax4.scatter(hypox,hypoy,10,color='black', marker='*')
ax4.set_xlabel(" X (Km)")
ax4.set_ylabel(" Y (Km)")
ax4.set_xlim(xmin,xmax)
ax4.set_ylim(ymin,ymax)
ax4.set_zlim(-15.0,15.0)
ax4.set_aspect('equal',adjustable='box')
ax4.azim = -90
ax4.dist = 10
ax4.elev = 90
plt.show()


IFMat = np.arange(0,XF3D.size).reshape((ndip,nstk),order='F')
ntri  = (nstk-1)*(ndip-1)*2
tri   = np.zeros([ntri,3],dtype=int)
XY3D  = np.array((XF3D,YF3D)).transpose()

# Delaunay triangulation
# tri = Delaunay(XY3D).simplices
ntri = int(tri.size/3)
jtri = -1
for istk in range (0,nstk-1):
    for idip in range (0,ndip-1):
        jtri += 1
        tri[jtri,0] = IFMat[idip,istk]
        tri[jtri,1] = IFMat[idip,istk+1]
        tri[jtri,2] = IFMat[idip+1,istk+1]
        jtri += 1
        tri[jtri,0] = IFMat[idip,istk]
        tri[jtri,1] = IFMat[idip+1,istk+1]
        tri[jtri,2] = IFMat[idip+1,istk]

triBmarker = np.ones(ntri,)
# Calculate unitary normal, strike and dip vector at each facet
univector = np.zeros((ntri,9))
# Vector normal to earth surface
nsurf = np.array([0,0,-1])

for itri in range(0,ntri):
    iv0 = tri[itri,0]
    iv1 = tri[itri,1]
    iv2 = tri[itri,2]
    v0 = np.array([ XF3D[iv0], YF3D[iv0], ZF3D[iv0]])
    v1 = np.array([ XF3D[iv1], YF3D[iv1], ZF3D[iv1]])
    v2 = np.array([ XF3D[iv2], YF3D[iv2], ZF3D[iv2]])
    vnormal = np.cross(v1-v0,v2-v0)
    vstrike = np.cross(vnormal,nsurf)
    vdip    = np.cross(vstrike,vnormal)
    univector[itri,0:3] = vnormal/np.linalg.norm(vnormal)
    univector[itri,3:6] = vstrike/np.linalg.norm(vstrike)
    univector[itri,6:9] = vdip/np.linalg.norm(vdip)


# Add nodes above an below to the fault
XF3Dplus = XF3D+vnormal[0]*dhF
YF3Dplus = YF3D+vnormal[1]*dhF
ZF3Dplus = ZF3D+vnormal[2]*dhF
XF3Dminus = XF3D-vnormal[0]*dhF
YF3Dminus = YF3D-vnormal[1]*dhF
ZF3Dminus = ZF3D-vnormal[2]*dhF

XF3Dadd = np.concatenate((XF3Dplus,XF3Dminus),axis=None)
YF3Dadd = np.concatenate((YF3Dplus,YF3Dminus),axis=None)
ZF3Dadd = np.concatenate((XF3Dplus,ZF3Dminus),axis=None)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(XF3D,YF3D,ZF3D, marker ='.', color='b',label ='Fault Nodes')
ax.scatter(XF3Dplus,YF3Dplus,ZF3Dplus, marker ='.', color='r', label='ExtraNodes above')
ax.scatter(XF3Dminus,YF3Dminus,ZF3Dminus, marker ='.', color='g', label="ExtraNodes below")
ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_zlabel(" Z (Km)")
ax.legend(loc ="upper left")


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_trisurf(XF3D,YF3D,ZF3D, triangles=tri)
ax.azim = -100
ax.dist = 10
ax.elev = 10
ax.set_title(' Delaunay Fault Triangulation ')

name = name+str(int(dhF*1000))+'m'

# Output Fault Dictionary
Fault = {}

Fault['dhF'] = dhF

Fault['nstk'] = nstk
Fault['ndip'] = ndip

Fault['XF3D'] = XF3D
Fault['YF3D'] = YF3D
Fault['ZF3D'] = ZF3D

Fault['XF3Dadd'] = XF3Dadd
Fault['YF3Dadd'] = YF3Dadd
Fault['ZF3Dadd'] = ZF3Dadd/m

Fault['stkVec'] = stkVec
Fault['dipVec'] = dipVec
Fault['stkinVec'] = stkinVec
Fault['dipinVec'] = dipinVec
Fault['stkMat'] = stkMat
Fault['dipMat'] = dipMat
Fault['stkinMat'] = stkinMat
Fault['dipinMat'] = dipinMat

Fault['ntri'] = ntri
Fault['tri']  = tri
Fault['triBmarker'] = triBmarker

Fault['SlipMat']  = SlipMat
Fault['RupTMat']  = RupTMat
Fault['RiseTMat'] = RiseTMat

Fault['fcoor']   = fcoor
Fault['outname'] = name

Fault['hypox'] = hypox
Fault['hypoy'] = hypoy
Fault['hypoz'] = hypoz
Fault['hypoistk'] = hypoistk
Fault['hypoidip'] = hypoidip


print()
print(" Output Fault Dimensions:")
print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " %(stk,nstk,np.linalg.norm(dstk)) )
print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f " %(dip,ndip,np.linalg.norm(ddip)) )

outfile = outFaultDir.joinpath(name+'.pickle')

fileObj = open(outfile, 'wb')
pickle.dump(Fault, fileObj)
fileObj.close()

print()
print(f" Fault info saved in:  {outfile} ")

# Write .vector file
fvectorHeader = "%d" %(ntri)
fvector = outFaultDir.joinpath(name+'.vector')
with open(fvector,'wb') as f:
    np.savetxt(f, univector,header=fvectorHeader, comments=' ',fmt='%10.6f')
f.close()

print()
print(f" Fault vector file written in:  {fvector} ")

# Output Topo Dictionary
Topo = {}
Topo['zTopoMat'] = zTopoMat/m
Topo['TopoLon'] = TopoLon
Topo['TopoLat'] = TopoLat
Topo['TopoLonMat'] = TopoLonMat
Topo['TopoLatMat'] = TopoLatMat
Topo['xTopoMat'] = xTopoMat
Topo['yTopoMat'] = yTopoMat

name ='Topo_'+str(Lonmin)+'a'+str(Lonmax)+'_'+str(Latmin)+'a'+str(Latmax)

outTopofile = outTopoDir.joinpath(name+'.pickle')

TopoObj = open(outTopofile, 'wb')
pickle.dump(Topo, TopoObj)
TopoObj.close()

print()
print(f" Topo info saved in:  {outTopofile} ")


print("  ")
print(" END PROGRAM ")
print("  ")




