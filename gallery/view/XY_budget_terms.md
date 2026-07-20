## Plot 2

![XY_budget_terms.png](../figures/XY_budget_terms.png)

````python
Panel = PanelPlot(2,2, [20,20],'COPT81 avec Mask', titlepad=11, minmaxpad=1.04, timepad=-0.07, colorbarpad=0.03, labelcolorbarpad = 13, colorbaraspect=40)

# Budget of potential temperature
nmask=0 #Convective mask, criteria in set_mask.f90
ntime=7 # 8th hour

Lplot = [Dvar['f1'][('/Budgets/TH','SFR')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','DEPS')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','DEPG')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','REVA')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','DEPI')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','IMLT')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','GMLT')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','DRYG')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','ACC')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','RIM')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','BERFI')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','CFRZ')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','WETG')][nmask,ntime,:]
         ]
LaxeZ = [Dvar['f1']['mask_level']]*len(Lplot)
Ltitle = ['Bilan température potentielle : partie convective - MASK1']*len(Lplot)
Llinelabel = ['SFR','DEPS','DEPG','REVA','DEPI','IMLT','GMLT','DRYG','ACC','RIM','BERFI','CFRZ','WETG']
Lxlim = [(-0.7E-2, 0.7E-2)]*len(Lplot)
Lxlab = ['Terme du bilan (K)']*len(Lplot)
Lylab = ['altitude (m)']*len(Lplot)
Lylim = [(0,12000.0)]*len(Lplot)
LaxisColor = ['black']*len(Lplot)
Llinewidth = [3]*len(Lplot)
LfacconvX=[1]*len(Lplot)
Llinecolor = ['red','green','blue','cyan','indigo','bisque','brown','orange','yellow',
              'magenta','gray','lightblue','black']

fig = Panel.pXY_lines(Lxx=Lplot, Lyy=LaxeZ, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Llinewidth=Llinewidth, LfacconvX=LfacconvX,
                        Lylim=Lylim, Lxlim=Lxlim, Llinelabel=Llinelabel, Llinecolor=Llinecolor,LaxisColor=LaxisColor,
                        ax=fig.axes)

nmask=1 #Convective mask, criteria in set_mask.f90
Lplot = [Dvar['f1'][('/Budgets/TH','SFR')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','DEPS')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','DEPG')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','REVA')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','DEPI')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','IMLT')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','GMLT')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','DRYG')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','ACC')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','RIM')][nmask,ntime,:],
         Dvar['f1'][('/Budgets/TH','BERFI')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','CFRZ')][nmask,ntime,:], Dvar['f1'][('/Budgets/TH','WETG')][nmask,ntime,:]
         ]
Lxlim = [(-0.7E-3, 0.7E-3)]*len(Lplot)
Ltitle = ['Bilan température potentielle : partie stratiforme - MASK2']*len(Lplot)

fig = Panel.pXY_lines(Lxx=Lplot, Lyy=LaxeZ, Lxlab=Lxlab, Lylab=Lylab, Ltitle=Ltitle, Llinewidth=Llinewidth, LfacconvX=LfacconvX,
                        Lylim=Lylim, Lxlim=Lxlim, Llinelabel=Llinelabel, Llinecolor=Llinecolor,LaxisColor=LaxisColor,
                        ax=fig.axes)
Panel.save_graph(1,fig)
````
