####################################################################################
#
# Python script to plot cubed-sphere grids
#
# Luan da Fonseca Santos - November 2022
# (luan.santos@usp.br)
####################################################################################

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from constants import *
####################################################################################
# This routine plots the grid.
####################################################################################
fig_format = 'png' # Figure format
def plot_grid(vert_lat, vert_lon, centers_lat, centers_lon, \
              midu_lat, midu_lon, midv_lat, midv_lon, \
              filename, N, map_projection):

   # Figure resolution
   dpi = 100


   print("--------------------------------------------------------")
   print('Plotting latlon grid using '+map_projection+' projection...')
   if map_projection == "mercator":
      plateCr = ccrs.PlateCarree()
      plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
   elif map_projection == 'sphere':
      plateCr = ccrs.Orthographic(central_longitude=-60.0, central_latitude=40.0)
      plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi)
   else:
      print('ERROR: Invalid map projection.')
      exit()

   plateCr._threshold = plateCr._threshold/10.
   ax = plt.axes(projection=plateCr)
   ax.stock_img()

   # Plot vertices
   for p in range(0, nbfaces):
       for i in range(0, N):
           for j in range(0, N):
               A_lon, A_lat = vert_lon[p,i,j], vert_lat[p,i,j]
               B_lon, B_lat = vert_lon[p,i+1,j], vert_lat[p,i+1,j]
               C_lon, C_lat = vert_lon[p,i+1,j+1], vert_lat[p,i+1,j+1]
               D_lon, D_lat = vert_lon[p,i,j+1], vert_lat[p,i,j+1]

               plt.plot([A_lon, B_lon], [A_lat, B_lat],linewidth=1, color='blue', transform=ccrs.Geodetic())
               plt.plot([B_lon, C_lon], [B_lat, C_lat],linewidth=1, color='blue', transform=ccrs.Geodetic())
               plt.plot([C_lon, D_lon], [C_lat, D_lat],linewidth=1, color='blue', transform=ccrs.Geodetic())
               plt.plot([D_lon, A_lon], [D_lat, A_lat],linewidth=1, color='blue', transform=ccrs.Geodetic())

   # Plot centers
   for p in range(0, nbfaces):
       for i in range(0, N):
           for j in range(0, N):
               A_lon, A_lat = centers_lon[p,i,j], centers_lat[p,i,j]

               plt.plot(A_lon, A_lat, color='blue', linewidth = 0.1, marker = '*', transform=ccrs.Geodetic())

   # Plot midu
   for p in range(0, nbfaces):
       for i in range(0, N+1):
           for j in range(0, N):
               A_lon, A_lat = midu_lon[p,i,j], midu_lat[p,i,j]

               plt.plot(A_lon, A_lat, color='blue', linewidth = 0.1, marker = '+', transform=ccrs.Geodetic())

   # Plot midv
   for p in range(0, nbfaces):
       for i in range(0, N):
           for j in range(0, N+1):
               A_lon, A_lat = midv_lon[p,i,j], midv_lat[p,i,j]

               plt.plot(A_lon, A_lat, color='blue', linewidth = 0.1, marker = '+', transform=ccrs.Geodetic())


   if map_projection == 'mercator':
      ax.gridlines(draw_labels=True)

   # Save the figure
   plt.savefig(graphdir+filename+"_"+map_projection+'.'+fig_format, format=fig_format)
   print('Figure has been saved in ../graphs/'+filename+'_'+map_projection+'.'+fig_format)
   print("--------------------------------------------------------\n")
   plt.close()
