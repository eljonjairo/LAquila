#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate fault's output files
#
#
# John Diaz June 2023
#

import warnings
import Fault

import plotly.express as px

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
    name = 'LAquilaCirella03_eff'

    # Output dir for topo file
    out_dir = '../Outputs/3DFaults/'

    dhF = 500  # Output subfaults size in m

    # Creates a new instance of the Fault class
    LAquilaFault = Fault.Fault(name,dhF)
    print(LAquilaFault)

    LAquilaFault.load_mat_file('../Input/s2009LAQUIL03CIRE.mat')
   
    #LAquilaFault.set_full_fault()
    LAquilaFault.set_effective_fault()
    LAquilaFault.plot_xyz_slipin()
    LAquilaFault.interpolate_xyz_coords()
    LAquilaFault.interpolate_slip()
    LAquilaFault.interpolate_rise_time()
    LAquilaFault.interpolate_rupt_time()
  
    LAquilaFault.plot_xyz_slip()
    LAquilaFault.plot_xyz_model_slip()
    LAquilaFault.compare_xyz_slip()
    
    # df = px.data.tips()
    # fig = px.scatter(df, x="total_bill", y="tip", color="size",
    #              title="Numeric 'size' values mean continuous color")

    # fig.show()
    
    LAquilaFault.triangulate_fault()
    LAquilaFault.add_nodes_above_below()
    LAquilaFault.plot_triangulation()
  
    LAquilaFault.write_univector(out_dir)

    LAquilaFault.save_fault(out_dir)















