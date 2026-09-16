## Plot 2

![XY_multisimple_GABLS1.png](../figures/XY_multisimple_GABLS1.png)

````python
Panel = PanelPlot(3,3, [25,25],'8-9h time averaged vertical profiles', titlepad=11, minmaxpad=1.04, timepad=-0.07, colorbarpad=0.03, labelcolorbarpad = 13, colorbaraspect=40)

Lplot = [np.mean(Dvar['f1'][(LG_MEAN,'MEAN_TH')][:,t_beg:t_end],axis=1), np.mean(Dvar['f2'][(LG_MEAN,'MEAN_TH')][:,t_beg:t_end],axis=1),
         np.mean(Dvar['f3'][(LG_MEAN,'MEAN_TH')][:,t_beg:t_end],axis=1)]
LaxeZ = [Dvar['f1']['level_les'], Dvar['f2']['level_les'], Dvar['f3']['level_les']]
Ltitle = ['MEAN_TH']*len(Lplot)
Llinelabel = ['1D BL89','1D RM17', 'LES']
Lxlab = ['theta (K)']*len(Lplot)
Lylab = ['altitude (m)']*len(Lplot)
Lylim = [(0,400.1)]*len(Lplot)
Lxlim = [(262, 268)]*len(Lplot)
Llinecolor = ['red','blue', 'black']
LaxisColor = ['black']*len(Lplot)
Llinewidth = [3]*len(Lplot)
fig = Panel.pXY_lines(Lxx=Lplot, Lyy=LaxeZ, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Llinewidth=Llinewidth,
                        Lylim=Lylim, Lxlim=Lxlim, Llinelabel=Llinelabel, Llinecolor=Llinecolor,LaxisColor=LaxisColor)

Lplot = [np.mean(Dvar['f1']['WIND'][:,t_beg:t_end],axis=1), np.mean(Dvar['f2']['WIND'][:,t_beg:t_end],axis=1), np.mean(Dvar['f3']['WIND'][:,t_beg:t_end],axis=1)]
Ltitle = ['Wind speed']*len(Lplot)
Lxlab = ['Wind speed (m/s)']*len(Lplot)
Lxlim = [(0, 11)]*len(Lplot)
fig = Panel.pXY_lines(Lxx=Lplot, Lyy=LaxeZ, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, ax=fig.axes, Llinewidth=Llinewidth,
                        Lylim=Lylim, Lxlim=Lxlim, Llinelabel=Llinelabel, Llinecolor=Llinecolor,LaxisColor=LaxisColor)

Lplot = [Dvar['f1']['SBL_H'][:], Dvar['f2']['SBL_H'][:], Dvar['f3']['SBL_H'][:]]

Ltitle = ['Boundary layer height']*len(Lplot)
LaxeTime = [Dvar['f1']['time_les']/3600.0, Dvar['f2']['time_les']/3600.0, Dvar['f3']['time_les']/3600.0]
Lxlab = ['Time (h)']*len(Lplot)
Lxlim = [(0, 9)]*len(Lplot)
Lylim = [(0, 300.1)]*len(Lplot)
fig = Panel.pXY_lines(Lyy=Lplot, Lxx=LaxeTime, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, ax=fig.axes,Llinewidth=Llinewidth,
                        Lylim=Lylim, Lxlim=Lxlim, Llinelabel=Llinelabel, Llinecolor=Llinecolor,LaxisColor=LaxisColor)
Panel.save_graph(1,fig)
````
