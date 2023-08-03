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
N = (16, 32, 64, 128,) # Values of N
reconmethods = ('hyppm', 'hyppm', 'hyppm', 'hyppm') # reconstruction methods
splitmethods = ('pl07' ,  'pl07', 'avlt', 'avlt' ) # splitting
mtmethods    = ('pl07' ,  'pl07', 'mt0' , 'mt0') # metric tensor formulation
dpmethods    = ('rk1'  ,  'rk1' , 'rk2' , 'rk2') # departure point formulation
mfixers      = ('none' ,  'gpr'  , 'af'  , 'gpr') # mass fixers 
edgetreat    = ('pl07' , 'duogrid', 'duogrid', 'duogrid') # edge treatments

# Program to be run
program = "./main"
run = True # Run the simulation?

# Plotting parameters
#map_projection = 'sphere'
map_projection = 'mercator'

def main():
    # Get the parameters
    _, kind  = get_parameters()

    # Define advection test in mesh.par'
    replace_line(pardir+'mesh.par', '4', 11)

    # remove output on screen
    replace_line(pardir+'mesh.par', '0', 9)

    # Define ic
    ic = '2'
    replace_line(pardir+'advection.par', ic, 3)

    # Define velocity
    vf = '1'
    replace_line(pardir+'advection.par', vf, 5)

    # interpolation degree
    interpd = '3'
    replace_line(pardir+'advection.par', interpd, 23)

    # time steps
    dt = np.zeros(len(N))

    # number of plots
    Nplots = 2
    replace_line(pardir+'advection.par', str(Nplots), 9)
    tplots = np.linspace(0, Nplots-1, Nplots, dtype=np.uint8)
    timeplots = np.linspace(0.0, 5.0, Nplots)

    # Initial time step
    if vf=='1':
        dt[0] = 0.025
    elif vf=='2'or vf=='4':
        dt[0] = 0.0125
    elif vf=='3':
        dt[0] = 0.00625
    else:
        print('Error - invalid vf')
        exit()

    # min/max in plot
    if ic=='1' or ic=='2' or ic=='3' or ic=='4':
       qmin, qmax = -0.2, 1.2
    else:
        print('Error - invalid ic')
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

            # Advection error name
            adv_name = "adv_ic"+ic+"_vf"+vf+"_"+opsplit+"_"+recon+"_mt"+mt+"_"+dp\
            +"_mf"+mf+"_et"+et+"_id"+interpd

            # File to be opened
            filename = datadir+adv_name+"_"+grid_name+"_errors.txt"
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
            # Plot
            # Open the file and reshape
            for t in tplots:
                # scalar field
                colormap = 'jet'
                fname = adv_name+'_Q_t'+str(t)+'_'+grid_name
                f = open(datadir+fname+'.dat', 'rb')
                data = np.fromfile(f, dtype=np.float64)
                data = np.reshape(data, (Nlat+1, Nlon+1))
                data = np.transpose(data)
                dmin, dmax = np.amin(data), np.amax(data)
                dmin = str("{:.2e}".format(dmin))
                dmax = str("{:.2e}".format(dmax))
                time = str("{:.2e}".format(timeplots[t]))
                title = 'N = '+str(n)+', Time = '+time+', CFL = '+str(cfl)+', ic = '+str(ic)+', vf = '+str(vf)+\
                '\n sp = '+str(opsplit)+', recon = '+ str(recon)+', dp = '+str(dp)+', mt = '\
                +str(mt)+', mf = '+str(mf) +', et = '+str(et)+\
                '\n min = '+dmin+', max = '+dmax+', mass variation = '+massvar+'\n \n'

                plot_scalar_field(data, lats, lons, \
                             colormap, map_projection, fname, title, qmin, qmax)

                # plot the error
                colormap = 'seismic'
                fname = adv_name+'_Q_error_t'+str(t)+'_'+grid_name
                if os.path.exists(datadir+fname+'.dat'): # The file exists
                    f = open(datadir+fname+'.dat', 'rb')
                    data = np.fromfile(f, dtype=np.float64)
                    data = np.reshape(data, (Nlat+1, Nlon+1))
                    data = np.transpose(data)
                    dabs_max = np.amax(abs(data))
                    dmin, dmax = -dabs_max, dabs_max
 
                    time = str("{:.2e}".format(timeplots[t]))
                    title = 'Error - N = '+str(n)+', Time = '+time+', CFL = '+str(cfl)+', ic = '+str(ic)+', vf = '+str(vf)+\
                    '\n sp = '+str(opsplit)+', recon = '+ str(recon)+', dp = '+str(dp)+\
                    ', mt = '+str(mt)+', mf = '+str(mf)+', et = '+str(et)+'\n \n'
                    plot_scalar_field(data, lats, lons, \
                                 colormap, map_projection, fname, title, dmin, dmax)

     
            #--------------------------------------------------------

    # plot errors for different all schemes in  different norms
    error_list = [error_linf, error_l1, error_l2]
    norm_list  = ['linf','l1','l2']
    norm_title  = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']

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
            fname.append(opsplit+'/'+recon+'/'+mt+'/'+dp+'/'+mf+'/'+et)

        title = 'Advection error, ic = '+str(ic)+', vf = '+ str(vf)+', norm='+norm_title[e]
        filename = graphdir+'cs_adv_ic'+str(ic)+'_vf'+str(vf)+'_norm'+norm_list[e]+'_parabola_errors.pdf'

        plot_errors_loglog(N, errors, fname, filename, title, emin, emax)

        # Plot the convergence rate
        title = 'Convergence rate, ic = '+str(ic)+', vf=' + str(vf)+', norm='+norm_title[e]
        filename = graphdir+'cs_adv_ic'+str(ic)+'_vf'+str(vf)+'_norm'+norm_list[e]+'_convergence_rate.pdf'
        plot_convergence_rate(N, errors, fname, filename, title, CRmin, CRmax)
        e = e+1

if __name__ == '__main__':
    main()
