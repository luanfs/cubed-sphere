#-----------------------------------------------------------------------
# Python script to plot scalar field outputs
# Luan Santos - 2022
#-----------------------------------------------------------------------

# Source code directory
srcdir = "src/"

import sys
import os.path
sys.path.append(srcdir)


import numpy as np
from configuration      import get_parameters, gridname, loadgrid
from plot               import plot_fields_list
from constants          import datadir, Nlat, Nlon


# Parameters
colormap = 'jet'
map_projection = 'mercator'

def main():
    # Get the parameters
    N, kind, midpoint, resolution = get_parameters()

    # Grid name
    grid_name = gridname(N, kind, midpoint, resolution)

    fields = ('area', 'length', 'sinc')
    plot_fields_list(fields, grid_name, colormap, map_projection)

main()

