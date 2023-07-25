#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Model class
# 

# TO DO LIST:
# Generate Topo
# Generate Velo full
# Generate velo fault
#
# John Diaz June 2023
#

import utm
import pygmt
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
import pickle
import numpy as np

class Model:
    def __init__(self, name, xmin, xmax, ymin, ymax ):
        self.name = name
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __str__(self):
        model = " Model's name: " + str(self.name) + "\n"
        model += " XUTM limits (km) [" + str(self.xmin) + ',' + str(self.xmax) + "]\n"
        model += " YUTM limits (km) [" + str(self.ymin) + ',' + str(self.ymax) + "]\n"

        return model

    def GenerateTopoItaly(self):

        print(" Generate Topo Italy ")

        # Transfrom from UTM to lat lon coordinates
        lat_min, lon_min = utm.to_latlon(self.xmin*1000, self.ymin*1000, 33, 'N')
        lat_max, lon_max = utm.to_latlon(self.xmax*1000, self.ymax*1000, 33, 'N')
        lat_min = round(lat_min*10-2)/10
        lat_max = round(lat_max*10+2)/10
        lon_min = round(lon_min*10-2)/10
        lon_max = round(lon_max*10+2)/10

        # Load Earth relief data for the entire globe and a subset region
        region = [lon_min,lon_max,lat_min,lat_max]
        grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

        zTopoMat = grid.data
        TopoLat = grid.lat.data
        TopoLon = grid.lon.data

        print(" Topo downloaded with limits: ")
        print(f" lon: [{lon_min},{lon_max}] and lat: [{lat_min},{lat_max}] ")
        print()

        fig = pygmt.Figure()
        pygmt.makecpt(cmap="geo", series=[-6000, 6000])
        fig.grdimage(grid=grid, region=region,
        frame=['WSrt+t" Topografia Italia Central"', "xa0.4", "ya0.4"])
        fig.colorbar(frame=["xa2000f500+lElevación ", "y+lm"])
        #fig.show()

        TopoLonMat, TopoLatMat = np.meshgrid(TopoLon,TopoLat)
        self.xTopoMat, self.yTopoMat, t1, t2 = utm.from_latlon(TopoLatMat, TopoLonMat,33,'N')

        self.out_name =  "Topo_" + str(lon_min) + 'a' + str(lon_max)
        self.out_name += "_" + str(lat_min) + 'a' + str(lat_max)

    def save(self,out_dir):

        outTopofile = out_dir + self.out_name + ".pickle"
        TopoObj = open(outTopofile, 'wb')
        pickle.dump(self, TopoObj)
        TopoObj.close()

        print(f" Topo info saved in:  {outTopofile} ")
        print()



        #cvals  = [ -100, 0, 100, 1000, 2000, 3000, 4000, ]
        #colors = ["blue","darkgreen","green","forestgreen","yellowgreen",
        #          "orange","maroon","sienna","brown","white"]
#
#        norm = plt.Normalize(min(cvals),max(cvals))
#        tuples = list(zip(map(norm,cvals), colors))
#        cTopo = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)




