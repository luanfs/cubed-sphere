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
from plot               import plot_scalar_field
from constants          import datadir, pardir, graphdir, Nlat, Nlon
from errors             import plot_errors_loglog, plot_convergence_rate
import subprocess

# Parameters
#N = (16, )
N = (16,32,64,128,256,512) # Values of N
reconmethods = ('ppm','ppm') # reconstruction methods
#splitmethods = ('avlt','avlt') # splitting
splitmethods = ('avlt','pl07') # splitting
mtmethods = ('mt0', 'pl07') # metric tensor formulation
dpmethods = ('rk1', 'rk1') # departure point formulation

# Program to be run
program = "./main"
run = True # Run the simulation?

# Plotting parameters
colormap = 'seismic'
map_projection = 'mercator'

def main():
    # Get the parameters
    _, kind  = get_parameters()

    # Define divergence test in mesh.par'
    replace_line(pardir+'mesh.par', '2', 11)

    # Define velocity
    vf = '1'
    replace_line(pardir+'advection.par', vf, 5)

    # Get adv parameters
    ic, _, recon = get_adv_parameters()

    # time steps
    dt = np.zeros(len(N))

    # Initial time step
    if vf=='1':
        dt[0] = 0.025
    elif vf=='2' or vf=='3' or vf=='4':
        dt[0] = 0.0125
    else:
        print('Error - invalid vf')
        exit()

    for k in range(1, len(N)):
        dt[k] = 0.5*dt[k-1]

    # Error arrays
    error_linf = np.zeros((len(N),len(reconmethods)))
    error_l1   = np.zeros((len(N),len(reconmethods)))
    error_l2   = np.zeros((len(N),len(reconmethods)))

    # Lat/lon aux vars
    lats = np.linspace(-90.0, 90.0, Nlat+1)
    lons = np.linspace(-180.0, 180.0, Nlon+1)
    lats, lons = np.meshgrid(lats, lons)

    # compile the code
    subprocess.run('cd .. ; make', shell=True)

    for i in range(0, len(reconmethods)):
        # Method parameters
        opsplit = splitmethods[i]
        recon = reconmethods[i]
        mt = mtmethods[i]
        dp = dpmethods[i]

        # Update reconstruction method in advection.par
        replace_line(pardir+'advection.par', recon, 9)

        # Update splitting method in advection.par
        replace_line(pardir+'advection.par', opsplit, 11)

        # Update metric tensor method in advection.par
        replace_line(pardir+'advection.par', mt, 13)

        # Update departure point  method in advection.par
        replace_line(pardir+'advection.par', dp, 15)


        k = 0
        for n in N:
            # Update time step
            replace_line(pardir+'advection.par', str(dt[k]), 7)

            # Grid name
            grid_name = gridname(n, kind)

            # Div error name
            div_name = "div_vf"+vf+"_"+recon+"_"+opsplit+"_"+mt

            # File to be opened
            filename = datadir+div_name+"_"+grid_name+"_errors.txt"
            print(filename)

            # Update N in mesh.par
            replace_line(pardir+'mesh.par', str(n), 3)

            # Run the program
            if (run):
                subprocess.run('cd .. ; ./main', shell=True)

            errors = np.loadtxt(filename)
            error_linf[k,i] = errors[0]
            error_l1[k,i] = errors[1]
            error_l2[k,i] = errors[2]
            k = k+1

            #--------------------------------------------------------
            # Plot error
            # Open the file and reshape
            fname = div_name+'_error_'+grid_name
            f = open(datadir+fname+'.dat', 'rb')
            data = np.fromfile(f, dtype=np.float64)
            data = np.reshape(data, (Nlat+1, Nlon+1))
            data = np.transpose(data)
            qabs_max = np.amax(abs(data))
            qmin, qmax = -qabs_max, qabs_max
            plot_scalar_field(data, lats, lons, \
                             colormap, map_projection, fname, qmin, qmax)
            if vf == '4':
                # Plot divergence
                # Open the file and reshape
                fname = div_name+'_'+grid_name
                f = open(datadir+fname+'.dat', 'rb')
                data = np.fromfile(f, dtype=np.float64)
                data = np.reshape(data, (Nlat+1, Nlon+1))
                data = np.transpose(data)
                plot_scalar_field(data, lats, lons, \
                                 colormap, map_projection, fname)
            #--------------------------------------------------------

    # plot errors for different all schemes in  different norms
    error_list = [error_linf, error_l1, error_l2]
    norm_list  = ['linf','l1','l2']
    norm_title  = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']

if __name__ == '__main__':
    main()
