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

    LAquilaFault.load_mat_file('../Input/s2009LAQUIL03CIRE.mat')
   
    #LAquilaFault.plot_xyz_slipin()
    LAquilaFault.interpolate_xyz_coords()
    LAquilaFault.interpolate_slip()
    LAquilaFault.interpolate_rise_time()
    LAquilaFault.interpolate_rupt_time()
  
   # LAquilaFault.plot_xyz_slip()
   # LAquilaFault.plot_xyz_model_slip()
    LAquilaFault.compare_xyz_slip()

    LAquilaFault.triangulate_fault()
    LAquilaFault.plot_triangulation()
    
    # plot own triangulation with scatter 3d and then lines
    #LAquilaFault.write_univector(out_dir)
    #LAquilaFault.save_fault(out_dir)















