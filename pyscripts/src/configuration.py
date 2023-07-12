####################################################################################
#
# This module contains all the routines that get the needed
# parameters from the /par directory.
#
# Luan da Fonseca Santos - January 2022
# (luan.santos@usp.br)
####################################################################################

from constants import*
import os.path
import numpy as np

def get_parameters():
   # The standard grid file mesh.par must exist in par/ directory
   file_path = pardir+"mesh.par"

   if os.path.exists(file_path): # The file exists
      # Open the grid file
      confpar = open(file_path, "r")

      # Read the grid file lines
      confpar.readline()
      confpar.readline()
      N = confpar.readline()
      confpar.readline()
      kind = confpar.readline()
      confpar.readline()
      gridload = confpar.readline()
      confpar.readline()
      showonscreen = confpar.readline()
      confpar.readline()

      # Close the file
      confpar.close()

      # rm \n
      N = rmspace(N)
      kind = rmspace(kind)

   else:   # The file does not exist
      print("ERROR in get_grid_parameters: file mesh.par not found in /par.")
      exit()

   return N, kind

def get_adv_parameters():
   # The standard adv file advection.par must exist in par/ directory
   file_path = pardir+"advection.par"

   if os.path.exists(file_path): # The file exists
      # Open the grid file
      confpar = open(file_path, "r")

      # Read the grid file lines
      confpar.readline()
      confpar.readline()
      ic = confpar.readline()
      confpar.readline()
      vf = confpar.readline()
      confpar.readline()
      flux = confpar.readline()

      # Close the file
      confpar.close()

      # rm \n
      ic = rmspace(ic)
      vf = rmspace(vf)
      flux = rmspace(flux)

   else:   # The file does not exist
      print("ERROR in get_adv_parameters: file advection.par not found in /par.")
      exit()

   return ic, vf, flux



'''-----------------------------------
 Give a name to the grid
-------------------------------------'''
def gridname(N, kind):
    filename = kind+"_"+str(N)
    return filename

'''-----------------------------------
 Remove the \n in a string
-------------------------------------'''
def rmspace(value):
    return ''.join(value.splitlines())

def loadgrid(filename, N):
    ng = 4
    # Open the vertices file
    f = open(griddir+filename+"_po.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+ng+1, N+ng+1,2))
    vert_lat = grid[:,:,:,0]*rad2deg
    vert_lon = grid[:,:,:,1]*rad2deg
    #print(np.amin(vert_lat), np.amax(vert_lat))
    #print(np.amin(vert_lon), np.amax(vert_lon))

    # Open the centers file
    f = open(griddir+filename+"_pc.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+ng, N+ng, 2))
    center_lat = grid[:,:,:,0]*rad2deg
    center_lon = grid[:,:,:,1]*rad2deg

    # Open the midu file
    f = open(griddir+filename+"_pu.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+ng+1, N+ng, 2))
    midu_lat = grid[:,:,:,0]*rad2deg
    midu_lon = grid[:,:,:,1]*rad2deg

    # Open the midv file
    f = open(griddir+filename+"_pv.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+ng, N+ng+1, 2))
    midv_lat = grid[:,:,:,0]*rad2deg
    midv_lon = grid[:,:,:,1]*rad2deg


    return vert_lat, vert_lon, center_lat, center_lon, \
              midu_lat, midu_lon, midv_lat, midv_lon


def replace_line(filename, content, line_number):
    import re
    if os.path.exists(filename): # The file exists
        # Open the grid file
        file  = open(filename, "r")
        lines = file.readlines()

        # new content
        lines[line_number-1] = content+'\n'

        # Close the file
        file.close()

        # Write new file
        with open(filename, 'w') as file:
            for line in lines:
                file.write(line)

    else:   # The file does not exist
        print("ERROR in edit_file_line: file"+filename+" not found in /par.")
        exit()

