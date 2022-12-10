#-----------------------------------------------------------------------
# Python script to run the discrete divergence test case
# for all reconstruction method.
# This script analyses the convergence with respect to the
# grid size.
# Luan Santos - 2022
#-----------------------------------------------------------------------

# Source code directory
srcdir = "src/"

import sys
import os.path
sys.path.append(srcdir)


import numpy as np
from configuration      import get_parameters, get_adv_parameters, gridname, replace_line
from plot               import plot_fields_list
from constants          import datadir, pardir, graphdir
from errors             import plot_errors_loglog, plot_convergence_rate
import subprocess

# Parameters
N = (16, 32, 64, 128, 256, 512, 1024) # Values of N
reconmethods = ('ppm', 'hyppm')       # reconstruction methods

# Program to be run
program = "./main"
run = True # Run the simulation?

# Plotting parameters
colormap = 'seismic'
map_projection = 'mercator'

def main():
    # Get the parameters
    _, kind, midpoint, resolution = get_parameters()

    # Define divergence test in mesh.par'
    replace_line(pardir+'mesh.par', '2', 13)

    # Get adv parameters
    ic, vf, recon = get_adv_parameters()

    # Error arrays
    error_linf = np.zeros((len(N),len(reconmethods)))
    error_l1   = np.zeros((len(N),len(reconmethods)))
    error_l2   = np.zeros((len(N),len(reconmethods)))

    # compile the code
    subprocess.run('cd .. ; make', shell=True)

    r = 0
    for recon in reconmethods:
        k = 0
        # Update reconstruction method in advection.par
        replace_line(pardir+'advection.par', recon, 7)

        for n in N:
            # Grid name
            grid_name = gridname(n, kind, midpoint, resolution)

            # Div error name
            div_name = "div_vf"+vf+"_"+recon

            # File to be opened
            filename = datadir+div_name+"_"+grid_name+"_errors.txt"
            print(filename)

            # Update N in mesh.par
            replace_line(pardir+'mesh.par', str(n), 3)

            # Run the program
            if (run):
                subprocess.run('cd .. ; ./main', shell=True)

            errors = np.loadtxt(filename)
            error_linf[k,r] = errors[0]
            error_l1[k,r] = errors[1]
            error_l2[k,r] = errors[2]
            k = k+1

            fields = (div_name, div_name+'_error')
            plot_fields_list(fields, grid_name, colormap, map_projection)

        r = r+1

    # Plot errors
    r = 0
    for recon in reconmethods:
        div_name = "div_vf"+vf+"_"+recon+"_"+ grid_name
        filename = graphdir+div_name+'_errors'
        title = div_name
        plot_errors_loglog(N, error_linf[:,r], error_l1[:,r], error_l2[:,r], filename, title)
        filename = graphdir+div_name+'_CR'
        plot_convergence_rate(N, error_linf[:,r], error_l1[:,r], error_l2[:,r], filename, title)
        r = r+1
main()

