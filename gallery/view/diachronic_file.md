## Read diachronic file

````python
LnameFiles = ['RBL89.1.ECH01.000.nc','GABL1.1.ECH01.000.nc']
LG_MEAN = '/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart/'
LG_SBG = '/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart/'
LG_RES = '/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart/'

Dvar_input = {
'f1':[(LG_SBG,'SBG_TKE'),(LG_SBG,'SBG_WU'),(LG_SBG,'SBG_WV'),(LG_SBG,'SBG_KM'),(LG_SBG,'SBG_KH'),(LG_SBG,'SBG_WTHL'),(LG_SBG,'SBG_THL2'),
      (LG_MEAN,'MEAN_U'),(LG_MEAN,'MEAN_V'),(LG_MEAN,'MEAN_TH'),
      'time_les','level_les'],
'f2':[(LG_SBG,'SBG_TKE'),(LG_SBG,'SBG_WU'),(LG_SBG,'SBG_WV'),(LG_SBG,'SBG_KM'),(LG_SBG,'SBG_KH'),(LG_SBG,'SBG_WTHL'),(LG_SBG,'SBG_THL2'),
      (LG_RES,'RES_KE'),(LG_RES,'RES_WU'),(LG_RES,'RES_WV'),(LG_RES,'RES_WTH'),(LG_RES,'RES_TH2'),
      (LG_MEAN,'MEAN_U'),(LG_MEAN,'MEAN_V'),(LG_MEAN,'MEAN_TH'),
      'time_les','level_les']
}

#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path=path, removeHALO=False)
````
