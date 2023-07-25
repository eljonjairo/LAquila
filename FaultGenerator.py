#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate fault's output files
#
#
# John Diaz June 2023
#
#from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
import plotly.graph_objects as go
import warnings
import Fault
#init_notebook_mode()
import plotly.io as io
io.renderers.default='svg'
if __name__ == '__main__':

    warnings.filterwarnings("ignore")
 
    
    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Fault Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03'

    # Output dir for topo file
    out_dir = '../Outputs/3DFaults/'

    dhF = 0.5  # Output subfaults size in Km

    # Creates a new instance of the Fault class
    LAquilaFault = Fault.Fault(name,dhF)
    print(LAquilaFault)

    LAquilaFault.LoadMat('../Input/s2009LAQUIL03CIRE.mat')
    LAquilaFault.PlotXYSlipin()
    LAquilaFault.PlotXYZSlipin(-120,10,10)   
    LAquilaFault.InterpolateXYZCoords()
    LAquilaFault.InterpolateSlip()
    LAquilaFault.PlotXYSlip()
    LAquilaFault.PlotXYZSlip(-120,10,10)   
    
    LAquilaFault.PlotlySlip()

















