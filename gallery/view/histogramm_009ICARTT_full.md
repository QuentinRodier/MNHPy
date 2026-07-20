## Plot 2

![histogramm_009ICARTT_full.png](../figures/histogramm_009ICARTT_full.png)

````python
#!/usr/bin/env python3
"""
@author: Quentin Rodier
Creation : 07/01/2021

Last modifications
"""
import matplotlib as mpl
mpl.use('Agg')
from read_MNHfile import read_netcdf
from Panel_Plot import PanelPlot
import cartopy.crs as ccrs
import os
import numpy as np

os.system('rm -f tempgraph*')
#
#  User's parameter / Namelist
#
path=""

LnameFiles = ['ICART.1.SEG01.001dg.nc', 'ICART.1.SEG01.002dg.nc']

Dvar_input = {
'f1':['MRC','COT','O3T','O3_PROD','O3_LOSS','CO_PROD','CO_LOSS','level','ZTOP', 'longitude','latitude','level_w','time',
      'CO_BUDGET','O3_BUDGET','O3_CHREACLIST','CO_CHREACLIST'],
'f2':['MRC','COT','O3T','O3_PROD','O3_LOSS','CO_PROD','CO_LOSS','level','ZTOP', 'longitude','latitude','level_w','time',
      'CO_BUDGET','O3_BUDGET','O3_CHREACLIST','CO_CHREACLIST']
}
#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path=path, removeHALO=True)

################################################################
#########          PANEL 1  # Horizontal cross-section
###############################################################
Panel1 = PanelPlot(3,3, [25,17],'Horizontal section at 1150m, 19h', titlepad=25, minmaxpad=1.04, timepad=-0.07, colorbarpad=0.01, colorbaraspect=40, labelcolorbarpad = 13)

Lplot = [ Dvar['f1']['MRC'][:,:,:], Dvar['f1']['COT'][:,:,:], Dvar['f1']['O3T'][:,:,:], Dvar['f1']['O3_PROD'][:,:,:],
         Dvar['f1']['O3_LOSS'][:,:,:], Dvar['f1']['CO_PROD'][:,:,:], Dvar['f1']['CO_LOSS'][:,:,:]]

LaxeX = [Dvar['f1']['longitude']]*len(Lplot)
LaxeY = [Dvar['f1']['latitude']]*len(Lplot)
Ltitle = ['Cloud mixing ratio', 'Carbon monoxyde CO ','Ozone O3', 'Ozone production', 'Ozone destruction','Carbon monoxyde production','Carbon monoxyde destruction']
Lcbarlabel = ['g/kg', 'ppbv','ppbv','ppbv/h','ppbv/h','ppbv/h','ppbv/h']
Lylab = ['latitude']*len(Lplot)
Lminval = [ 0, 107.5, 0, 70, 70, 0.5, 0.5 ]
Lmaxval = [ 0.2, 137.5, 70, 130, 130, 1.7, 1.7 ]
Lstep = [ 0.01, 2.5, 5, 5, 5, 0.1, 0.1]
Lstepticks = Lstep
Lfacconv = [ 1, 1, 1, 1e9*3600, 1e9*3600,1e9*3600, 1e9*3600]
Lcolormap = ['gist_ncar']*len(Lplot)
Llvl = [14]*len(Lplot)
LaddWhite_cm = [True, False, False, False, False, False, False]
Lprojection = [ccrs.PlateCarree()]*len(Lplot)

fig1 = Panel1.psectionH(lon=LaxeX, lat=LaxeY, Lvar=Lplot, Llevel=Llvl,Lylab=Lylab, Ltitle=Ltitle, Lminval=Lminval, Lmaxval=Lmaxval, 
                                Lstep=Lstep, Lstepticks=Lstepticks, Lcolormap=Lcolormap, Lcbarlabel=Lcbarlabel, Lfacconv=Lfacconv, 
                                colorbar=True, LaddWhite_cm=LaddWhite_cm, Lproj=Lprojection)

Panel1.save_graph(1,fig1)

################################################################
#########          PANEL 2  # Horizontal cross-section
###############################################################
Panel2 = PanelPlot(3,3, [25,17],'Horizontal section at 1150m, 20h', titlepad=25, minmaxpad=1.04, timepad=-0.07, colorbarpad=0.01, colorbaraspect=40, labelcolorbarpad = 13)

Lplot = [ Dvar['f2']['MRC'][:,:,:], Dvar['f2']['COT'][:,:,:], Dvar['f2']['O3T'][:,:,:], Dvar['f2']['O3_PROD'][:,:,:],
         Dvar['f2']['O3_LOSS'][:,:,:], Dvar['f2']['CO_PROD'][:,:,:], Dvar['f2']['CO_LOSS'][:,:,:]]

fig2 = Panel2.psectionH(lon=LaxeX, lat=LaxeY, Lvar=Lplot, Llevel=Llvl,Lylab=Lylab, Ltitle=Ltitle, Lminval=Lminval, Lmaxval=Lmaxval, 
                                Lstep=Lstep, Lstepticks=Lstepticks, Lcolormap=Lcolormap, Lcbarlabel=Lcbarlabel, Lfacconv=Lfacconv, 
                                colorbar=True, LaddWhite_cm=LaddWhite_cm, Lproj=Lprojection)

Panel2.save_graph(2,fig2)

################################################################
#########          PANEL 3  # Bar plots Budget chemical reactions
###############################################################
Dvar['f1']['CO_BUDGET_mean'] = np.mean(Dvar['f1']['CO_BUDGET'][:,13,:,:],axis=(1,2))  #  {x,y} Average on height = 1150m
Dvar['f1']['O3_BUDGET_mean'] = np.mean(Dvar['f1']['O3_BUDGET'][:,13,:,:],axis=(1,2))  #  {x,y} Average on height = 1150m
Dvar['f2']['CO_BUDGET_mean'] = np.mean(Dvar['f2']['CO_BUDGET'][:,13,:,:],axis=(1,2))  #  {x,y} Average on height = 1150m
Dvar['f2']['O3_BUDGET_mean'] = np.mean(Dvar['f2']['O3_BUDGET'][:,13,:,:],axis=(1,2))  #  {x,y} Average on height = 1150m

Panel3 = PanelPlot(2,2, [20,20],'Chemical budgets', titlepad=25, minmaxpad=1.04, timepad=-0.07, colorbarpad=0.01, colorbaraspect=40, labelcolorbarpad = 13)

Lplot = [Dvar['f1']['CO_BUDGET_mean'], Dvar['f1']['O3_BUDGET_mean'], Dvar['f2']['CO_BUDGET_mean'], Dvar['f2']['O3_BUDGET_mean'] ]
Lbins=[Dvar['f1']['CO_CHREACLIST'], Dvar['f1']['O3_CHREACLIST'],Dvar['f2']['CO_CHREACLIST'], Dvar['f2']['O3_CHREACLIST']]
Ltitle = ['Carbon monoxyde CO chemical reactions', 'Ozone O3 chemical reactions']*2
Lylab = ['Budget (ppbv/h)']*len(Lplot)
Ltime = [Dvar['f1']['time'], Dvar['f1']['time'], Dvar['f2']['time'], Dvar['f2']['time']]
Lylim=[(-1,1), (-100, 100)]*2
Lfacconv = [1E9*3600]*len(Lplot)
Lwidth=[0.95]*len(Lplot)
Lcolors=[]
for var in Lplot:
    cc=['']*len(var)
    for n,val in enumerate(var):
        if val<0:
            cc[n]='blue'
        elif val>=0:
            cc[n]='red'
    Lcolors.append(cc)

fig3 = Panel3.pXY_bar(Lbins=Lbins, Lvar=Lplot, Lylim=Lylim, Lfacconv=Lfacconv, Ltitle=Ltitle, Lylab=Lylab, Lcolor=Lcolors, Lwidth=Lwidth, Ltime=Ltime)

#  Handle a new axis at y=0 for each graphs
for i,var in enumerate(Lplot):
    ax2 = fig3.axes[i].twinx() #  Clone the existing axis
    ax2_x = ax2.get_xaxis()
    ax2_x.set_label('Chemical reactions')
    ax2_y = ax2.get_yaxis() #  Get the new Y axe and hide it
    ax2_y.set_visible(False)
    fig3.axes[i].spines['bottom'].set_position('center') #  Move the original axis to the center

Panel3.save_graph(3,fig3)
````
