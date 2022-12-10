####################################################################################
#
# This module includes routines that compute the error in L_inf, L_1 and L_2 norms
# and plot the error.
# Luan da Fonseca Santos - May 2022
#
####################################################################################

import matplotlib.pyplot as plt
import numpy as np

####################################################################################
# Routine that plot the errors in different norms in log-log scale
####################################################################################
def plot_errors_loglog(N, error_linf, error_l1, error_l2, filename, title):
    # Plot the error graph
    plt.loglog(N, error_linf, color='green', marker='x', label = '$L_\infty$')
    plt.loglog(N, error_l1 , color='blue',  marker='o', label = '$L_1$')
    plt.loglog(N, error_l2 , color='red',   marker='D', label = '$L_2$')
    #plt.ylim(10.0**(-3),np.amax(10*error_linf))

    # Reference lines
    nref   = len(N)
    Norder = N[0:nref]
    order1 = np.zeros(nref)
    order2 = np.zeros(nref)
    order3 = np.zeros(nref)
    ref = np.amax(error_linf)
    order1[0], order2[0], order3[0] = ref, ref, ref

    order4 = np.zeros(nref)
    order4[0] = ref
    for i in range(1, nref):
        order1[i] = order1[i-1]/2.0
        order2[i] = order2[i-1]/4.0
        order3[i] = order3[i-1]/8.0
        order4[i] = order4[i-1]/16.0
    plt.loglog(Norder, order1 , ':' , color='black', label = '1st order')
    #plt.loglog(Norder, order2 , '--', color='black', label = '2nd order')
    #plt.loglog(Norder, order3 , '-.', color='black', label = '3rd order')
    #plt.loglog(Norder, order4 , '--', color='black', label = '4rd order')
    plt.xlabel('N (number of cells)')
    plt.ylabel('Error')
    plt.legend()
    plt.grid(True, which="both")
    plt.title(title)
    plt.savefig(filename, format='png')
    plt.close()

####################################################################################
# Compute and plot the convergence rate
####################################################################################
def plot_convergence_rate(N, error_linf, error_l1, error_l2, filename, title):
    n = len(N)

    CR_linf = np.abs(np.log(error_linf[1:n])-np.log(error_linf[0:n-1]))/np.log(2.0)
    CR_l1   = np.abs(np.log(error_l1[1:n])  -np.log(error_l1[0:n-1]))/np.log(2.0)
    CR_l2   = np.abs(np.log(error_l2[1:n])  -np.log(error_l2[0:n-1]))/np.log(2.0)
    #plt.ylim(0,2)
    plt.xscale('log')
    plt.plot(N[1:n], CR_linf, color='green', marker='x', label = '$L_\infty$')
    plt.plot(N[1:n], CR_l1, color='blue',  marker='o', label = '$L_1$')
    plt.plot(N[1:n], CR_l2, color='red',   marker='D', label = '$L_2$')

    plt.xlabel('N (number of cells)')
    plt.ylabel('Convergence rate')
    plt.legend()
    plt.grid(True, which="both")
    plt.title(title)
    plt.savefig(filename, format='png')
    plt.close()

####################################################################################
# Print the errors of the ith simulation
####################################################################################
def print_errors_simul(error_linf, error_l1, error_l2, i):
    # Output
    if i > 0:
        print('Norms          (Linf,    L1,       L2) ')
        print('Error E_'+str(i)+'    :',"{:.2e}".format(error_linf[i]), "{:.2e}".format(error_l1[i]), "{:.2e}".format(error_l2[i]))
        print('Ratio E_'+str(i)+'/E_'+str(i-1)+':',"{:.2e}".format(error_linf[i-1]/error_linf[i]), "{:.2e}".format(error_l1[i-1]/error_l1[i]), "{:.2e}".format(error_l2[i-1]/error_l2[i]))
    else:
        print('Error(Linf, L1, L2) :',"{:.2e}".format(error_linf[i]), "{:.2e}".format(error_l1[i]), "{:.2e}".format(error_l2[i]))

####################################################################################
# Returns the L_inf, L_1 and L_2 errors
####################################################################################
def compute_errors(Q, Qref):
    # Relative errors in different metrics
    E = abs(Qref - Q)
    n = E.size

    # L_inf error
    error_inf = np.amax(abs(E))

    # L_1 error
    error_1 = np.sum(E)/n

    # L_2 error
    error_2 = np.sqrt(np.sum(E*E)/n)

    return error_inf, error_1, error_2

####################################################################################
# Plot evolution with time of the values given in the list values
####################################################################################
def plot_time_evolution(values, Tf, vlabels, ylabel, filename, title):
    # Plot the error graph
    n = len(values)
    Nsteps = len(values[0])
    times = np.linspace(0, Tf, Nsteps)
    colors = ('black', 'blue', 'green', 'red', 'purple')
    for k in range(0, n):
        plt.semilogy(times, values[k], color=colors[k], label = vlabels[k])

    plt.ylabel(ylabel)
    plt.xlabel('Time (seconds)')
    plt.legend()
    plt.grid(True, which="both")
    plt.title(title)
    plt.savefig(filename)
    plt.close()
