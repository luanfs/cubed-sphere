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
N = (16, 32, 64, 128, 256, 512) # Values of N

#N = (16,32,64) # Values of N

reconmethods = ('ppm',) # reconstruction methods
#reconmethods = ('hyppm',)       # reconstruction methods

splitmethods =('lr96',)
#splitmethods =('lr96', 'l04', 'pl07')

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
    error_linf = np.zeros((len(N),len(reconmethods),len(splitmethods)))
    error_l1   = np.zeros((len(N),len(reconmethods),len(splitmethods)))
    error_l2   = np.zeros((len(N),len(reconmethods),len(splitmethods)))

    # compile the code
    subprocess.run('cd .. ; make', shell=True)

    r = 0
    o = 0
    for recon in reconmethods:
        o = 0
        for opsplit in splitmethods:
            # Update reconstruction method in advection.par
            replace_line(pardir+'advection.par', recon, 7)

            # Update splitting method in advection.par
            replace_line(pardir+'advection.par', opsplit, 9)

            k = 0
            for n in N:
                # Grid name
                grid_name = gridname(n, kind, midpoint, resolution)

                # Div error name
                div_name = "div_vf"+vf+"_"+recon+"_"+opsplit

                # File to be opened
                filename = datadir+div_name+"_"+grid_name+"_errors.txt"
                print(filename)

                # Update N in mesh.par
                replace_line(pardir+'mesh.par', str(n), 3)

                # Run the program
                if (run):
                    subprocess.run('cd .. ; ./main', shell=True)

                #print(k,r,o)
                errors = np.loadtxt(filename)
                error_linf[k,r,o] = errors[0]
                error_l1[k,r,o] = errors[1]
                error_l2[k,r,o] = errors[2]
                k = k+1

                fields = (div_name, div_name+'_error')
                plot_fields_list(fields, grid_name, colormap, map_projection)

            o = o+1
        r = r+1

    # Plot errors
    r = 0
    o = 0
    for recon in reconmethods:
        o = 0
        for opsplit in splitmethods:
            div_name = "div_vf"+vf+"_"+recon+"_"+opsplit+"_"+ grid_name
            filename = graphdir+div_name+'_errors'
            title = div_name
            plot_errors_loglog(N, error_linf[:,r,o], error_l1[:,r,o], error_l2[:,r,o], filename, title)
            filename = graphdir+div_name+'_CR'
            plot_convergence_rate(N, error_linf[:,r,o], error_l1[:,r,o], error_l2[:,r,o], filename, title)
            o = o+1
        r = r+1

main()
