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
N = (48, ) # Values of N
reconmethods = ('ppm', ) # reconstruction methods
splitmethods = ('avlt', ) # splitting
mtmethods    = ('mt0', ) # metric tensor formulation
dpmethods    = ('rk1'  ,  'rk1' ) # departure point formulation
mfixers      = ('none' ,  'none') # mass fixers 
edgetreat    = ('duogrid', 'duogrid') # edge treatments


# Program to be run
program = "./main"
run = True # Run the simulation?
#run = False # Run the simulation?

# Plotting parameters
colormap = 'seismic'
map_projection = 'mercator'
#map_projection = 'sphere'

def main():
    # Get the parameters
    _, kind  = get_parameters()

    # Define divergence test in mesh.par'
    replace_line(pardir+'mesh.par', '3', 11)

    # Define velocity
    vf = '1'
    replace_line(pardir+'advection.par', vf, 5)
    ic = '1'
    replace_line(pardir+'advection.par', ic, 3)

    # remove output on screen
    replace_line(pardir+'mesh.par', '0', 9)

    # number of plots
    Nplots = 2
    replace_line(pardir+'advection.par', str(Nplots), 9)
 
    # interpolation degree
    interpd = '3'
    replace_line(pardir+'advection.par', interpd, 23)

    # Get adv parameters
    _, _, recon = get_adv_parameters()

    # time steps
    dt = np.zeros(len(N))

    # Initial time step
    if vf=='1':
        dt[0] = 3600
    if vf=='2':
        dt[0] = 2400
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
        mf = mfixers[i]
        et = edgetreat[i]

        # Update reconstruction method in advection.par
        replace_line(pardir+'advection.par', recon, 11)

        # Update splitting method in advection.par
        replace_line(pardir+'advection.par', opsplit, 13)

        # Update metric tensor method in advection.par
        replace_line(pardir+'advection.par', mt, 15)

        # Update departure point method in advection.par
        replace_line(pardir+'advection.par', dp, 17)

        # Update mass fixer in advection.par
        replace_line(pardir+'advection.par', mf, 19)

        # Update edge treatment advection.par
        replace_line(pardir+'advection.par', et, 21)

        k = 0
        for n in N:
            # Update time step
            replace_line(pardir+'advection.par', str(dt[k]), 7)

            # Grid name
            grid_name = gridname(n, kind)

            # Div error name
            div_name = "div_ic"+str(ic)+"_vf"+vf+"_"+opsplit+"_"+recon+"_mt"+mt+"_"+dp+\
            "_mf"+mf+"_et"+et+"_id"+interpd

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
            cfl = errors[3]
            massvar = errors[4]
            cfl = str("{:.2e}".format(cfl))
            massvar = str("{:.2e}".format(massvar))
            k = k+1

            #--------------------------------------------------------
            # Plot error
            # Open the file and reshape
            fname = div_name+'_div_error_'+grid_name
            f = open(datadir+fname+'.dat', 'rb')
            data = np.fromfile(f, dtype=np.float64)
            data = np.reshape(data, (Nlat+1, Nlon+1))
            data = np.transpose(data)
            qabs_max = np.amax(abs(data))
            qmin, qmax = -qabs_max, qabs_max
            title = 'N = '+str(n)+', CFL = '+str(cfl)+', ic = '+str(ic)+', vf = '+str(vf)\
            +', sp = '+str(opsplit)+'\n recon = '+ str(recon)+', dp = '+str(dp)+', mt = '\
            +str(mt)+', mf = '+str(mf) +', et = '+str(et)+\
            ', mass variation = '+massvar+'\n \n'

            plot_scalar_field(data, lats, lons, \
                             colormap, map_projection, fname, title, qmin, qmax)
            if vf == '4':
                # Plot divergence
                # Open the file and reshape
                fname = div_name+'_div_'+grid_name
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
    scheme_names = ['PL07', 'PL07-PR', 'AVLT-AF', 'AVLT-PR']
    e = 0
    for error in error_list:
        emin, emax = np.amin(error[:]), np.amax(error[:])

        # convergence rate min/max
        n = len(error)
        CR = np.abs(np.log(error[1:n])-np.log(error[0:n-1]))/np.log(2.0)
        CRmin, CRmax = np.amin(CR), np.amax(CR)
        errors = []
        fname = []
        for k in range(0, len(reconmethods)):
            errors.append(error[:,k])
            opsplit = splitmethods[k]
            recon = reconmethods[k]
            mt = mtmethods[k]
            dp = dpmethods[k]
            mf = mfixers[k]
            et = edgetreat[k]
            #name.append(opsplit+'/'+recon+'/'+mt+'/'+dp+'/'+mf+'/'+et)
            fname.append(scheme_names[k])

        title = 'Divergence error, vf='+ str(vf)+', norm='+norm_title[e]
        filename = graphdir+'cs_div_vf'+str(vf)+'_norm'+norm_list[e]+'_parabola_errors.pdf'

        plot_errors_loglog(N, errors, fname, filename, title, emin, emax)

        # Plot the convergence rate
        title = 'Convergence rate, vf=' + str(vf)+', norm='+norm_title[e]
        filename = graphdir+'cs_div_vf'+str(vf)+'_norm'+norm_list[e]+'_convergence_rate.pdf'
        plot_convergence_rate(N, errors, fname, filename, title, CRmin, CRmax)
        e = e+1

if __name__ == '__main__':
    main()
