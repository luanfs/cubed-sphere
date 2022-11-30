####################################################################################
#
# This module contains all the constants needed for the other modules
#
# Luan da Fonseca Santos - January 2022
####################################################################################

import math

# Physical parameters
day2sec        = 24.0*3600.0    # Days in seconds
sec2day        = 1.0/day2sec    # Seconds in days
erad           = 6371.0*10**3   # Earth mean radius (m)
eradi          = 1.0/erad       # Inverse of erad
rotation_speed = 7.29212e-5     # Angular speed of rotation of the Earth (radians/s)
grav           = 9.80616        # Gravitational acceleration of the Earth (m s^-2)
gravi          = 1.0/grav       # Inverse of grav

# Mathematical constants
pi      = math.pi
pio2    = pi*0.5
pio4    = pio2*0.5
rad2deg = 180.0/pi              # Radians to degrees conversion
deg2rad = 1.0/rad2deg           # Degrees to radians conversion

# Grid parameters
nbfaces = 6                     # Number of faces

# Directory parameters
griddir  = "../grid/"              # Output grid data directory
datadir  = "../data/"              # Output data directory
graphdir = "../graphs/"            # Graphs directory
pardir   = "../par/"               # Parameter files directory

# Plotting parameters
Nlat = 1000                     # Number of points in the latlon grid used for plotting
Nlon = 2*Nlat
