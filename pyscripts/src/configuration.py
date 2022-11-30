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
      midpoint = confpar.readline()
      confpar.readline()
      gridload = confpar.readline()
      confpar.readline()
      showonscreen = confpar.readline()
      confpar.readline()

      # Close the file
      confpar.close()

      resolution = 'unif'

      # rm \n
      N = rmspace(N)
      kind = rmspace(kind)
      midpoint = rmspace(midpoint)
      resolution = rmspace(resolution)

   else:   # The file does not exist
      print("ERROR in get_grid_parameters: file mesh.par not found in /par.")
      exit()

   return N, kind, midpoint, resolution

'''-----------------------------------
 Give a name to the grid
-------------------------------------'''
def gridname(N, kind, midpoint, resolution):
    filename = kind+"_"+str(N)+"_"+midpoint+"_"+resolution
    return filename

'''-----------------------------------
 Remove the \n in a string
-------------------------------------'''
def rmspace(value):
    return ''.join(value.splitlines())

def loadgrid(filename, N):
    # Open the vertices file
    f = open(griddir+filename+"_vert.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+1, N+1,2))
    vert_lat = grid[:,:,:,0]*rad2deg
    vert_lon = grid[:,:,:,1]*rad2deg
    #print(np.amin(vert_lat), np.amax(vert_lat))
    #print(np.amin(vert_lon), np.amax(vert_lon))

    # Open the centers file
    f = open(griddir+filename+"_center.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N, N, 2))
    center_lat = grid[:,:,:,0]*rad2deg
    center_lon = grid[:,:,:,1]*rad2deg

    # Open the midu file
    f = open(griddir+filename+"_midu.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N+1, N, 2))
    midu_lat = grid[:,:,:,0]*rad2deg
    midu_lon = grid[:,:,:,1]*rad2deg

    # Open the midu file
    f = open(griddir+filename+"_midv.dat",'rb')
    grid = np.fromfile(f, dtype=np.float64)
    grid = np.reshape(grid,(nbfaces, N, N+1, 2))
    midv_lat = grid[:,:,:,0]*rad2deg
    midv_lon = grid[:,:,:,1]*rad2deg


    return vert_lat, vert_lon, center_lat, center_lon, \
              midu_lat, midu_lon, midv_lat, midv_lon
