#-----------------------------------------------------------------------
# This script analyses the convergence with respect to the
# grid size of the shallow water model
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
N = (16, 32, 64, 128, 256) # Values of N
reconmethods = ('hyppm', 'hyppm', 'hyppm') # reconstruction methods
splitmethods = ( 'pl07', 'avlt', 'avlt' ) # splitting
mtmethods    = ( 'pl07', 'mt0' , 'mt0') # metric tensor formulation
dpmethods    = ( 'rk1' , 'rk2' , 'rk2') # departure point formulation
mfixers      = ( 'gpr'  , 'af'  , 'gpr') # mass fixers 
edgetreat    = ('duogrid', 'duogrid', 'duogrid') # edge treatments


#reconmethods = ('hyppm',) # reconstruction methods
#splitmethods = ('avlt',  ) # splitting
#mtmethods    = ('mt0' , ) # metric tensor formulation
#dpmethods    = ('rk2' , ) # departure point formulation
#mfixers      = ('gpr'  , ) # mass fixers 
#edgetreat    = ('duogrid', )# edge treatments


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
    replace_line(pardir+'mesh.par', '5', 11)

    # remove output on screen
    replace_line(pardir+'mesh.par', '0', 9)

    # Define ic
    ic = '0'
    replace_line(pardir+'swm.par', ic, 3)

    # final integration time (days)
    tf = '12'
    replace_line(pardir+'swm.par', tf, 5)

    # interpolation degree
    interpd = '3'
    replace_line(pardir+'swm.par', interpd, 23)

    # average process degree (linear or cubic)
    avd = '3'
    replace_line(pardir+'swm.par', avd, 25)

    # time steps
    dt = np.zeros(len(N))

    # number of plots
    Nplots = 2
    replace_line(pardir+'swm.par', str(Nplots), 9)
    tplots = np.linspace(0, Nplots-1, Nplots, dtype=np.uint8)
    timeplots = np.linspace(0.0, 5.0, Nplots)

    # Initial time step
    if ic=='0' or ic=='1' or ic=='2':
        dt[0] = 8000
        #dt[0] = 0.025 #3000
    else:
        print('Error - invalid ic')
        exit()

    # min/max in plot
    if ic=='1':
       qmin, qmax = 400, 1200
    elif ic=='2':
       qmin, qmax = 800.0, 3400.0
    else:
        if ic != '0':
            print('Error - invalid ic')
            exit()


    for k in range(1, len(N)):
        dt[k] = 0.5*dt[k-1]

    # Error arrays
    error_linf = np.zeros((len(N),len(reconmethods)))
    error_l1   = np.zeros((len(N),len(reconmethods)))
    error_l2   = np.zeros((len(N),len(reconmethods)))

    # consistency error
    if ic=='0':
        error_linf_div = np.zeros((len(N),len(reconmethods)))
        error_l1_div   = np.zeros((len(N),len(reconmethods)))
        error_l2_div   = np.zeros((len(N),len(reconmethods)))
        error_linf_rv = np.zeros((len(N),len(reconmethods)))
        error_l1_rv   = np.zeros((len(N),len(reconmethods)))
        error_l2_rv   = np.zeros((len(N),len(reconmethods)))
        error_linf_av = np.zeros((len(N),len(reconmethods)))
        error_l1_av   = np.zeros((len(N),len(reconmethods)))
        error_l2_av   = np.zeros((len(N),len(reconmethods)))
        error_linf_h_po = np.zeros((len(N),len(reconmethods)))
        error_linf_av_pu = np.zeros((len(N),len(reconmethods)))
        error_linf_av_pv = np.zeros((len(N),len(reconmethods)))
        error_linf_u_po = np.zeros((len(N),len(reconmethods)))
        error_linf_v_po = np.zeros((len(N),len(reconmethods)))


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

        # Update reconstruction method in swm.par
        replace_line(pardir+'swm.par', recon, 11)

        # Update splitting method in swm.par
        replace_line(pardir+'swm.par', opsplit, 13)

        # Update metric tensor method in swm.par
        replace_line(pardir+'swm.par', mt, 15)

        # Update departure point method in swm.par
        replace_line(pardir+'swm.par', dp, 17)

        # Update mass fixer in swm.par
        replace_line(pardir+'swm.par', mf, 19)

        # Update edge treatment swm.par
        replace_line(pardir+'swm.par', et, 21)

        k = 0
        for n in N:
            # Update time step
            replace_line(pardir+'swm.par', str(dt[k]), 7)

            # Grid name
            grid_name = gridname(n, kind)

            # swm error name
            swm_name = "swm_ic"+ic+"_"+opsplit+"_"+recon+"_mt"+mt+"_"+dp\
            +"_mf"+mf+"_et"+et+"_id"+interpd+"_av"+avd

            # File to be opened
            filename = datadir+swm_name+"_"+grid_name+"_errors.txt"
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

            if ic=='0':
                error_linf_div[k,i] = errors[5]
                error_l1_div[k,i]   = errors[6]
                error_l2_div[k,i]   = errors[7]
                error_linf_rv[k,i] = errors[8]
                error_l1_rv[k,i]   = errors[9]
                error_l2_rv[k,i]   = errors[10]
                error_linf_av[k,i] = errors[11]
                error_l1_av[k,i]   = errors[12]
                error_l2_av[k,i]   = errors[13]
                error_linf_h_po[k,i] = errors[14]
                error_linf_av_pu[k,i] = errors[15]
                error_linf_av_pv[k,i] = errors[16]
                error_linf_u_po[k,i] = errors[17]
                error_linf_v_po[k,i] = errors[18]


            k = k+1

            #--------------------------------------------------------
            # Plot
            # Open the file and reshape
            for t in tplots:
                # scalar field
                colormap = 'jet'
                fname = swm_name+'_H_t'+str(t)+'_'+grid_name
                f = open(datadir+fname+'.dat', 'rb')
                data = np.fromfile(f, dtype=np.float64)
                data = np.reshape(data, (Nlat+1, Nlon+1))
                data = np.transpose(data)
                dmin, dmax = np.amin(data), np.amax(data)
                dmin = str("{:.2e}".format(dmin))
                dmax = str("{:.2e}".format(dmax))
                time = str("{:.2e}".format(timeplots[t]))
                title = 'N = '+str(n)+', Time = '+time+', CFL = '+str(cfl)+', ic = '+str(ic)+\
                '\n sp = '+str(opsplit)+', recon = '+ str(recon)+', dp = '+str(dp)+', mt = '\
                +str(mt)+', mf = '+str(mf) +', et = '+str(et)+\
                '\n min = '+dmin+', max = '+dmax+', mass variation = '+massvar+'\n \n'

                if ic != '0':
                    plot_scalar_field(data, lats, lons, \
                             colormap, map_projection, fname, title, qmin, qmax)

                # plot the error
                colormap = 'seismic'
                fields = ['H',]
                for fd in fields:
                    fname = swm_name+'_'+fd+'_error_t'+str(t)+'_'+grid_name
                    if os.path.exists(datadir+fname+'.dat'): # The file exists
                        f = open(datadir+fname+'.dat', 'rb')
                        data = np.fromfile(f, dtype=np.float64)
                        data = np.reshape(data, (Nlat+1, Nlon+1))
                        data = np.transpose(data)
                        dabs_max = np.amax(abs(data))
                        dmin, dmax = -dabs_max, dabs_max

                        time = str("{:.2e}".format(timeplots[t]))
                        title = 'Error - N = '+str(n)+', Time = '+time+', CFL = '+str(cfl)+', ic = '+str(ic)+\
                        '\n sp = '+str(opsplit)+', recon = '+ str(recon)+', dp = '+str(dp)+\
                        ', mt = '+str(mt)+', mf = '+str(mf)+', et = '+str(et)+'\n \n'
                        plot_scalar_field(data, lats, lons, \
                                colormap, map_projection, fname, title, dmin, dmax)

                if ic=='0':
                    # plot the error
                    colormap = 'seismic'
                    #colormap = 'jet'
                    fields = ['div', 'rel_vort', 'abs_vort']
                    #fields = ['div',]
                    names = ['Divergence', 'Relative vorticity', 'Absolute vorticity',]

                    for fd in range(0,len(fields)):
                        fname = swm_name+'_'+fields[fd]+'_error_t'+str(t)+'_'+grid_name
                        if os.path.exists(datadir+fname+'.dat'): # The file exists
                            f = open(datadir+fname+'.dat', 'rb')
                            data = np.fromfile(f, dtype=np.float64)
                            data = np.reshape(data, (Nlat+1, Nlon+1))
                            data = np.transpose(data)
                            dabs_max = np.amax(abs(data))
                            dmin, dmax = -dabs_max, dabs_max
                            #dmin, dmax = 0.0, 1.0

                            time = str("{:.2e}".format(timeplots[t]))
                            title = names[fd]+\
                            ' error - N = '+str(n)+', Time = '+time+', CFL = '+str(cfl)+', ic = '+str(ic)+\
                            '\n sp = '+str(opsplit)+', recon = '+ str(recon)+', dp = '+str(dp)+\
                            ', mt = '+str(mt)+', mf = '+str(mf)+', et = '+str(et)+'\n \n'
                            plot_scalar_field(data, lats, lons, \
                                    colormap, map_projection, fname, title, dmin, dmax)


            #--------------------------------------------------------

    # plot errors for different all schemes in  different norms
    error_list = [error_linf, ]#error_l1, error_l2]
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

        title = 'SWM error, ic = '+str(ic)+', norm='+norm_title[e]
        filename = graphdir+'cs_swm_ic'+str(ic)+'_norm'+norm_list[e]+'_parabola_errors.pdf'

        plot_errors_loglog(N, errors, fname, filename, title, emin, emax)

        # Plot the convergence rate
        title = 'Convergence rate, ic = '+str(ic)+', norm='+norm_title[e]
        filename = graphdir+'cs_swm_ic'+str(ic)+'_norm'+norm_list[e]+'_convergence_rate.pdf'
        plot_convergence_rate(N, errors, fname, filename, title, CRmin, CRmax)
        e = e+1

    # consistency error
    if ic=='0':
        # plot the error
        fields = ['div', 'rel_vort', 'abs_vort', 'h_po', 'av_pu', 'av_pv', 'u_po', 'v_po']
        names = ['Divergence', 'Relative vorticity', 'Absolute vorticity',\
        'Height field at B grid', 'Absolute vorticity at pu', 'Absolute vorticity at pv',\
        'Covariant u at B grid', 'Covariant v at B grid']

        error_list_div   = [error_linf_div, ]#error_l1_div, error_l2_div]
        error_list_rv    = [error_linf_rv , ]#error_l1_rv , error_l2_rv]
        error_list_av    = [error_linf_av , ]#error_l1_av , error_l2_av]
        error_list_h_po  = [error_linf_h_po, ]
        error_list_av_pu  = [error_linf_av_pu, ]
        error_list_av_pv  = [error_linf_av_pv, ]
        error_list_u_po  = [error_linf_u_po, ]
        error_list_v_po  = [error_linf_v_po, ] 
        error_lists = [error_list_div, error_list_rv, error_list_av, error_list_h_po,\
        error_list_av_pu, error_list_av_pv, error_list_u_po, error_list_v_po]

        for fd in range(0,len(fields)):
            error_list = error_lists[fd]
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

                title = names[fd]+' error, ic = '+str(ic)+', norm='+norm_title[e]
                filename = graphdir+'cs_swm_'+fields[fd]+'_ic'+str(ic)+'_norm'+norm_list[e]+'_parabola_errors.pdf'

                plot_errors_loglog(N, errors, fname, filename, title, emin, emax)

                # Plot the convergence rate
                title = names[fd]+' - convergence rate, ic = '+str(ic)+', norm='+norm_title[e]
                filename = graphdir+'cs_swm_'+fields[fd]+'_ic'+str(ic)+'_norm'+norm_list[e]+'_convergence_rate.pdf'
                plot_convergence_rate(N, errors, fname, filename, title, CRmin, CRmax)
                e = e+1


if __name__ == '__main__':
    main()
