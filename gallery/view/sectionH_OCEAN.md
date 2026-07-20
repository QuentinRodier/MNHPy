## Plot 2

![sectionH_OCEAN.png](../figures/sectionH_OCEAN.png)

````python
#!/usr/bin/env python3
"""
@author: Quentin Rodier
Creation : 23/06/2021

Last modifications
"""
import matplotlib as mpl
mpl.use('Agg')
from read_MNHfile import read_netcdf
from Panel_Plot import PanelPlot
from misc_functions import *
import cartopy.crs as ccrs
import numpy as np
import os
#
#  User's parameter / Namelist
#
path=""
LnameFiles = ['SPWAN.1.25m00.003.nc','SPWAN.2.25m00.003.nc']

Dvar_input = {'f1':['WT','TKET','THT','level_w','ni','nj'],
              'f2':['WT','TKET','THT','level_w','ni','nj']}

#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path=path, removeHALO=False)

################################################################
#########          PANEL 1
###############################################################
Panel = PanelPlot(2,3, [25,15],'', titlepad=20, minmaxpad=1.03, timepad=-0.10, colorbarpad=0.03, labelcolorbarpad = 13, colorbaraspect=35, lateralminmaxpad=1.02)

Lplot = [ Dvar['f1']['WT'] ,Dvar['f1']['TKET'], Dvar['f1']['THT']-273.15,
          Dvar['f2']['WT'] ,Dvar['f2']['TKET'], Dvar['f2']['THT']-273.15 ]

lon = [Dvar['f1']['ni'], Dvar['f1']['ni'], Dvar['f1']['ni'],
       Dvar['f2']['ni'], Dvar['f2']['ni'], Dvar['f2']['ni']]
lat = [Dvar['f1']['nj'], Dvar['f1']['nj'], Dvar['f1']['nj'],
       Dvar['f2']['nj'], Dvar['f2']['nj'], Dvar['f2']['nj']]
Llevel = [97]*len(Lplot) 
Ltitle = ['Vertical velocity D1', 'Subgrid TKE D1', 'Temperature D1','Vertical velocity D2', 'Subgrid TKE D2', 'Temperature D2']
Lcbarlabel = ['cm/s','m2/s2','°C']*2
Lxlab = ['X (m)']*len(Lplot)
Lylab = ['Y (m)']*len(Lplot)
Lminval = [-7., 0, 10.31]*2
Lmaxval = [7., 2E-4, 10.3625]*2
Lstep = [0.1, 5E-6, 1E-5]*2
Lstepticks = [1, 2.5E-5, 1E-2]*2
Lcolormap = ['seismic','gist_rainbow_r','gist_rainbow_r']*2
Lfacconv = [100.0,1,1]*2
LaddWhite = [False,True,False]*2
Lcbformatlabel=[False,True,False]*2

fig = Panel.psectionH(lon=lon, lat=lat, Llevel=Llevel, Lvar=Lplot, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Lminval=Lminval, Lmaxval=Lmaxval, 
                                Lstep=Lstep, Lstepticks=Lstepticks, Lcolormap=Lcolormap, Lcbarlabel=Lcbarlabel, 
                                LaddWhite_cm=LaddWhite, Lfacconv=Lfacconv, Lcbformatlabel=Lcbformatlabel)
fig.tight_layout()
Panel.save_graph(1,fig)
````
