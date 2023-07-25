#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate model's output files
#
#
# John Diaz June 2023
#

import warnings

import Model


if __name__ == '__main__':

    warnings.filterwarnings("ignore")
    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Model Program       * ")
    print(" ************************************************ ")
    print()

    # Model Name
    name = 'LAquila'

    # Output dir for topo file
    out_dir = '../Outputs/Topo/'

    # X and Y UTM limits (Km) for Central Italy Model 13 Stations
    xmin = 300.0
    xmax = 420.0
    ymin = 4610.0
    ymax = 4750.0
    zmin = -60.0

    # Creates a new instance of the Model class
    LAquilaModel = Model.Model(name,xmin,xmax,ymin,ymax)
    print(LAquilaModel)
    LAquilaModel.GenerateTopoItaly()
    LAquilaModel.save(out_dir)






