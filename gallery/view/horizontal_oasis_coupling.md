## Plot 2

![horizontal_oasis_coupling.png](../figures/horizontal_oasis_coupling.png)

````python
#!/bin/python3
# --------------------------------------------------------
#
#                 Author  (    date    ) :
#             J. Pianezze ( 29.09.2023 )
#
#                    ~~~~~~~~~~~~~~~
#       Script used to verify OASIS exchanges between
#                Meso-NH and toy models
#                    ~~~~~~~~~~~~~~~
#
# --------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os, glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.colors import BoundaryNorm
curdir_path = os.getcwd()+'/'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #########################################################
# ###           To be defined by user                   ###
# #########################################################

name_file_send_mnh = glob.glob('*mesonh_01.nc')[0]
name_file_recv_mnh = glob.glob('*mesonh_02.nc')[0]
name_file_send_toy = glob.glob('*toyexe_02.nc')[0]
name_file_recv_toy = glob.glob('*toyexe_01.nc')[0]
name_var01         = name_file_send_mnh[5:8]
name_var02         = name_file_recv_mnh[5:8]

# #########################################################

# ---------------------------------------
#   Create directory to store figures
# ---------------------------------------
try :
  os.mkdir(curdir_path+name_var01+'_'+name_var02+'/')
except OSError:
  print('Directory already created')
else:
  print('Making directory')

# ---------------------------------------
#   Open files
# ---------------------------------------
file_send_toy = netCDF4.Dataset(name_file_send_toy)
file_recv_mnh = netCDF4.Dataset(name_file_recv_mnh)
file_recv_toy = netCDF4.Dataset(name_file_recv_toy)
file_send_mnh = netCDF4.Dataset(name_file_send_mnh)
file_grids    = netCDF4.Dataset('grids.nc')

# ---------------------------------------
#   Read lon/lat
# ---------------------------------------
lon_toy = file_grids.variables['toyt.lon'] ; nlon_toy = np.shape(lon_toy)[1]
lat_toy = file_grids.variables['toyt.lat'] ; nlat_toy = np.shape(lat_toy)[0]
lon_mnh = file_grids.variables['ssea.lon'] ; nlon_mnh = np.shape(lon_mnh)[1]
lat_mnh = file_grids.variables['ssea.lat'] ; nlat_mnh = np.shape(lat_mnh)[0]

# ---------------------------------------
#   Read variables
# ---------------------------------------
var_send_toy = file_send_toy.variables[name_file_send_toy[0:8]][-2,:,:]*1000.0
var_recv_toy = file_recv_toy.variables[name_file_recv_toy[0:8]][-1,:,:]
var_send_mnh = file_send_mnh.variables[name_file_send_mnh[0:8]][-2,:,:]
var_recv_mnh = file_recv_mnh.variables[name_file_recv_mnh[0:8]][-1,:,:]*1000.0

mask_mnh     = (var_send_mnh[:,:] > 1E10)
var_send_mnh = np.ma.MaskedArray(var_send_mnh, mask=mask_mnh)
var_recv_mnh = np.ma.MaskedArray(var_recv_mnh, mask=mask_mnh)

mask_toy    = (var_recv_toy[:,:] == 0.0)
var_send_toy = np.ma.MaskedArray(var_send_toy, mask=mask_toy)
var_recv_toy = np.ma.MaskedArray(var_recv_toy, mask=mask_toy)

# -----------------------------------------------------------
#   Create figure
# -----------------------------------------------------------
fig = plt.figure()

# -----------------------------------------------------------
#   Define colormap and norm
# -----------------------------------------------------------
cmap_wnd   = plt.cm.RdBu_r
cmap_toy   = plt.cm.RdBu_r

levels_wnd = np.arange(  0.0, 0.26, 0.01)
levels_toy = np.arange(-10.0, 10.1,  0.1)

norm_wnd   = BoundaryNorm(levels_wnd, ncolors=cmap_wnd.N, clip=True)
norm_toy   = BoundaryNorm(levels_toy, ncolors=cmap_toy.N, clip=True)

#----------------------
ax   = fig.add_subplot(221)
plt.title('(a) Send by MNH')
cs   = plt.pcolormesh(lon_mnh[:,:],lat_mnh[:,:],var_send_mnh[:,:],cmap=cmap_wnd,norm=norm_wnd)
cbar = plt.colorbar(cs,orientation='vertical',format='%.2f')
plt.tick_params(axis='x',which='both',labelbottom=False)
ax.set_xlim(( max(np.min(lon_mnh[1:-1,1:-1]),np.min(lon_toy[1:-1,1:-1])), min(np.max(lon_mnh[1:-1,1:-1]),np.max(lon_toy[1:-1,1:-1])) ))
ax.set_ylim(( max(np.min(lat_mnh[1:-1,1:-1]),np.min(lat_toy[1:-1,1:-1])), min(np.max(lat_mnh[1:-1,1:-1]),np.max(lat_toy[1:-1,1:-1])) ))

#----------------------
ax   = fig.add_subplot(222)
plt.title('(b) Received by TOY')
cs   = plt.pcolormesh(lon_toy[:,:],lat_toy[:,:],var_recv_toy[:,:],cmap=cmap_wnd,norm=norm_wnd)
cbar = plt.colorbar(cs,orientation='vertical',format='%.2f')
plt.tick_params(axis='x',which='both',labelbottom=False)
plt.tick_params(axis='y',which='both',labelleft  =False)
ax.set_xlim(( max(np.min(lon_mnh[1:-1,1:-1]),np.min(lon_toy[1:-1,1:-1])), min(np.max(lon_mnh[1:-1,1:-1]),np.max(lon_toy[1:-1,1:-1])) ))
ax.set_ylim(( max(np.min(lat_mnh[1:-1,1:-1]),np.min(lat_toy[1:-1,1:-1])), min(np.max(lat_mnh[1:-1,1:-1]),np.max(lat_toy[1:-1,1:-1])) ))

#----------------------
ax   = fig.add_subplot(223)
plt.title('(c) Send by TOY')
cs   = plt.pcolormesh(lon_toy[:,:],lat_toy[:,:],var_send_toy[:,:],cmap=plt.cm.RdBu_r,vmin=np.min(var_send_toy), vmax=np.max(var_send_toy))  
cbar = plt.colorbar(cs,orientation='vertical',format='%.1f')
ax.set_xlim(( max(np.min(lon_mnh[1:-1,1:-1]),np.min(lon_toy[1:-1,1:-1])), min(np.max(lon_mnh[1:-1,1:-1]),np.max(lon_toy[1:-1,1:-1])) ))
ax.set_ylim(( max(np.min(lat_mnh[1:-1,1:-1]),np.min(lat_toy[1:-1,1:-1])), min(np.max(lat_mnh[1:-1,1:-1]),np.max(lat_toy[1:-1,1:-1])) ))
   
#----------------------
ax   = fig.add_subplot(224)
plt.title('(d) Received by MNH')
cs   = plt.pcolormesh(lon_mnh[:,:],lat_mnh[:,:],var_recv_mnh[:,:],cmap=plt.cm.RdBu_r,vmin=np.min(var_send_toy), vmax=np.max(var_send_toy))  
cbar = plt.colorbar(cs,orientation='vertical',format='%.1f')
plt.tick_params(axis='y',which='both',labelleft=False)
ax.set_xlim(( max(np.min(lon_mnh[1:-1,1:-1]),np.min(lon_toy[1:-1,1:-1])), min(np.max(lon_mnh[1:-1,1:-1]),np.max(lon_toy[1:-1,1:-1])) ))
ax.set_ylim(( max(np.min(lat_mnh[1:-1,1:-1]),np.min(lat_toy[1:-1,1:-1])), min(np.max(lat_mnh[1:-1,1:-1]),np.max(lat_toy[1:-1,1:-1])) ))

#------------------------
plt.savefig(curdir_path+name_var01+"_"+name_var02+"/"+name_var01+"_"+name_var02+".png", bbox_inches='tight', dpi=400)
plt.close()
#------------------------
````
