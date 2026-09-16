## Plot 2

![sectionH_vectors_007janvier_full.png](../figures/sectionH_vectors_007janvier_full.png)

````python
#!/usr/bin/env python3
"""
@author: Quentin Rodier
Creation : 07/01/2021

Last modifications
"""
import matplotlib as mpl
mpl.use('Agg')
import cartopy.crs as ccrs
from read_MNHfile import read_netcdf
from Panel_Plot import PanelPlot
import os

os.system('rm -f tempgraph*')
#
#  User's parameter / Namelist
#
#
path=""
LnameFiles = ['16JAN.1.12B18.001dg.nc', '16JAN.2.12B18.001dg.nc']

Dvar_input = {
'f1':['MRV700HPA','THT850HPA','UT850HPA','VT850HPA','UT700HPA','VT700HPA', 'ALT_PRESSURE','ALT_U','ALT_V', 'ZS', 'latitude', 'longitude'],
'f2':['MRV700HPA','THT850HPA','UT850HPA','VT850HPA','UT700HPA','VT700HPA', 'ALT_PRESSURE', 'ZS', 'ALT_U','ALT_V','latitude', 'longitude']
}

#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path=path, removeHALO=True)

################################################################
#########          PANEL 1  
###############################################################
Panel1 = PanelPlot(2,2, [20,20],'007_janvier domaine 1 16JAN.1.12B18.001dg.nc', minmaxpad=1.05)

Lplot = [ Dvar['f1']['ZS'],Dvar['f1']['THT850HPA'], Dvar['f1']['MRV700HPA'],Dvar['f1']['ALT_PRESSURE']]
lon = [Dvar['f1']['longitude']]*len(Lplot)
lat = [Dvar['f1']['latitude']]*len(Lplot)
Ltitle = ['Orography', 'Potential Temperature at 850hPa', 'Water vapor mixing at 700hPa','Pressure at z = 9000m']
Lcbarlabel = ['m','K', 'g/kg', 'hPa']
Lxlab = ['longitude']*len(Lplot)
Lylab = ['latitude']*len(Lplot)
Lminval = [0, 285, 0.9, 286]
Lmaxval = [300, 289, 2.6, 294]
Lstep = [10, 0.25, 0.1, 0.4]
Lstepticks = [50, 1, 0.2, 0.4]
Lfacconv = [1.0, 1.0, 1.0, 1./100.0]
Lcolormap = ['terrain', 'gist_rainbow_r', 'gist_rainbow_r', 'gist_rainbow_r']
Lprojection = [ccrs.PlateCarree()]*len(Lplot)
Llvl = [0, 0, 0, 0]
fig1 = Panel1.psectionH(lon=lon, lat=lat, Lvar=Lplot, Lcarte=[], Llevel=Llvl, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Lminval=Lminval, Lmaxval=Lmaxval, 
                                Lstep=Lstep, Lstepticks=Lstepticks, Lcolormap=Lcolormap, Lcbarlabel=Lcbarlabel, Lproj=Lprojection, Lfacconv=Lfacconv)

Lplot1 = [ Dvar['f1']['UT850HPA'], Dvar['f1']['UT700HPA'], Dvar['f1']['ALT_U']]
Lplot2 = [ Dvar['f1']['VT850HPA'], Dvar['f1']['VT700HPA'], Dvar['f1']['ALT_V']]
Ltitle = ['Wind at 850hPa', 'Wind at 700hPa', 'Wind at 9000m']
Lxlab = ['longitude']*len(Lplot1)
Lylab = ['latitude']*len(Lplot1)
Llegendval = [20,20,40]
Lcbarlabel = ['(m/s)']*len(Lplot1)
Larrowstep = [2]*len(Lplot1)
Lwidth = [0.002]*len(Lplot1)
Lcolor = ['black']*len(Lplot1)
Lprojection = [ccrs.PlateCarree()]*len(Lplot1)
Llvl = [0]*len(Lplot1)
fig2 = Panel1.pvector(Lxx=lon, Lyy=lat, Lvar1=Lplot1, Lvar2=Lplot2, Lcarte=[], Llevel=Llvl, Lxlab=Lxlab, Lylab=Lylab, 
                      Ltitle=Ltitle, Lwidth=Lwidth, Larrowstep=Larrowstep, Lproj=Lprojection,
                      Lcolor=Lcolor, Llegendval=Llegendval, Lcbarlabel=Lcbarlabel, Lid_overlap=[2,4,6], ax=fig1.axes)

Panel1.save_graph(1,fig2)

################################################################
#########          PANEL 2 
###############################################################
Panel2 = PanelPlot(2,2, [20,20],'007_janvier domaine 2 16JAN.1.12B18.001dg.nc', minmaxpad=1.05)

Lplot = [ Dvar['f2']['ZS'],Dvar['f2']['THT850HPA'], Dvar['f2']['MRV700HPA'],Dvar['f2']['ALT_PRESSURE']]
lon = [Dvar['f2']['longitude']]*len(Lplot)
lat = [Dvar['f2']['latitude']]*len(Lplot)
Ltitle = ['Orography', 'Potential Temperature at 850hPa', 'Water vapor mixing at 700hPa','Pressure at z = 9000m']
Lcbarlabel = ['m','K', 'g/kg', 'hPa']
Lxlab = ['longitude']*len(Lplot)
Lylab = ['latitude']*len(Lplot)
Lminval = [0, 285, 0.9, 286]
Lmaxval = [300, 289, 2.6, 294]
Lstep = [10, 0.25, 0.1, 0.4]
Lstepticks = [50, 1, 0.2, 0.4]
Lfacconv = [1.0, 1.0, 1.0, 1./100.0]
Lcolormap = ['terrain', 'gist_rainbow_r', 'gist_rainbow_r', 'gist_rainbow_r']
Lprojection = [ccrs.PlateCarree()]*len(Lplot)
Llvl = [0]*len(Lplot)
fig1 = Panel2.psectionH(lon=lon, lat=lat, Lvar=Lplot, Lcarte=[], Llevel=Llvl, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Lminval=Lminval, Lmaxval=Lmaxval, 
                                Lstep=Lstep, Lstepticks=Lstepticks, Lcolormap=Lcolormap, Lcbarlabel=Lcbarlabel, Lproj=Lprojection, Lfacconv=Lfacconv)

Lplot1 = [ Dvar['f2']['UT850HPA'], Dvar['f2']['UT700HPA'], Dvar['f2']['ALT_U']]
Lplot2 = [ Dvar['f2']['VT850HPA'], Dvar['f2']['VT700HPA'], Dvar['f2']['ALT_V']]
Ltitle = ['Wind at 850hPa', 'Wind at 700hPa', 'Wind at 9000m']
Llegendval = [20,20,40]
Lxlab = ['longitude']*len(Lplot1)
Lylab = ['latitude']*len(Lplot1)
Lcbarlabel = ['(m/s)']*len(Lplot1)
Larrowstep = [2]*len(Lplot1)
Lwidth = [0.002]*len(Lplot1)
Lcolor = ['black']*len(Lplot1)
Lprojection = [ccrs.PlateCarree()]*len(Lplot1)
Llvl = [0]*len(Lplot1)
fig2 = Panel2.pvector(Lxx=lon, Lyy=lat, Lvar1=Lplot1, Lvar2=Lplot2, Lcarte=[], Llevel=Llvl, Lxlab=Lxlab, Lylab=Lylab, 
                      Ltitle=Ltitle, Lwidth=Lwidth, Larrowstep=Larrowstep, Lproj=Lprojection,
                      Lcolor=Lcolor, Llegendval=Llegendval, Lcbarlabel=Lcbarlabel, Lid_overlap=[2,4,6], ax=fig1.axes)

Panel2.save_graph(2,fig2)
````
