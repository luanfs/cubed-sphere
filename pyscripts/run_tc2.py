#-----------------------------------------------------------------------
# Python script to run the interpolation at duogrid test case
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
N = (16, 32, 64, 128, 256,) # Values of N
degrees = (1, 2, 3, 4,) # interpolation degrees

# Program to be run
program = "./main"
run = True # Run the simulation?


def main():
    # Get the parameters
    _, kind  = get_parameters()

    # Define interpolation test in mesh.par'
    replace_line(pardir+'mesh.par', '2', 11)

    # Define velocity and scalar field
    vf = '1'
    ic = '5'

    # update parameters
    replace_line(pardir+'interpolation.par', ic, 3)
    replace_line(pardir+'interpolation.par', vf, 5)

    # Error arrays
    error_q = np.zeros((len(N),len(degrees)))
    error_ucontra = np.zeros((len(N),len(degrees)))
    error_ucovari = np.zeros((len(N),len(degrees)))

    # compile the code
    subprocess.run('cd .. ; make', shell=True)

    for i in range(0, len(degrees)):
        # Method parameters
        d = degrees[i]

        # Update interpolation order in interpolation.par
        replace_line(pardir+'interpolation.par', str(d), 7)

        k = 0
        for n in N:
            # Grid name
            grid_name = gridname(n, kind)

            # interp error name
            interp_name = "interp_ic"+ic+"_vf"+vf+"_id"+str(d)

            # File to be opened
            filename = datadir+interp_name+"_"+grid_name+"_errors.txt"
            print(filename)

            # Update N in mesh.par
            replace_line(pardir+'mesh.par', str(n), 3)

            # Run the program
            if (run):
                subprocess.run('cd .. ; ./main', shell=True)

            errors = np.loadtxt(filename)
            error_q[k,i] = errors[0]
            error_ucontra[k,i] = errors[1]
            error_ucovari[k,i] = errors[2]

            k = k + 1

    # plot errors for different all schemes in  different norms
    error_list = [ error_q, error_ucontra, error_ucovari]
    name = ['ic'+str(ic), 'ucontra_vf'+str(vf), 'ucovari_vf'+str(vf)]
    titles = ['Ghost cells of scalar field' , 'Ghost cells of vector field', \
              'Ghost cells of vector field']
    e = 0
    for error in error_list:
        emin, emax = np.amin(error[:]), np.amax(error[:])

        # convergence rate min/max
        n = len(error)
        CR = np.abs(np.log(error[1:n])-np.log(error[0:n-1]))/np.log(2.0)
        CRmin, CRmax = np.amin(CR), np.amax(CR)
        errors = []
        fname = []

        for k in range(0, len(degrees)):
            errors.append(error[:,k])
            fname.append('Degree - '+str(degrees[k]))

        title = 'Interpolation error - ' + titles[e]
        filename = graphdir+'interp_'+name[e]+'_parabola_errors.pdf'
        plot_errors_loglog(N, errors, fname, filename, title, emin, emax)

        # Plot the convergence rate
        title = 'Convergence rate - ' + titles[e]
        filename = graphdir+'interp_'+name[e]+'_convergence_rate.pdf'
        plot_convergence_rate(N, errors, fname, filename, title, 0, 5)
        e = e+1

if __name__ == '__main__':
    main()
