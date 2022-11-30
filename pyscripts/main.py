#########################################################
#
# Cubed-sphere generation code
#
# Luan da Fonseca Santos (luan.santos@usp.br)
#
##########################################################

# Source code directory
srcdir = "src/"

import sys
import os.path
sys.path.append(srcdir)

# Imports
from configuration      import get_parameters, gridname, loadgrid
from constants          import Nlat, Nlon, griddir, datadir, graphdir
from plot               import plot_grid
def main():
    # Get the parameters
    N, kind, midpoint, resolution = get_parameters()

    # Grid name
    filename = gridname(N, kind, midpoint, resolution)

    # Load the grid data
    vert_lat, vert_lon, center_lat, center_lon, \
    midu_lat, midu_lon, midv_lat, midv_lon = loadgrid(filename, int(N))

    # Plot the grid
    map_projection = "mercator"
    plot_grid(vert_lat, vert_lon, center_lat, center_lon, \
              midu_lat, midu_lon, midv_lat, midv_lon, \
              filename, int(N), map_projection)

main()


