## Read multiple file

````python
LnameFiles = ['ICART.1.SEG01.001dg.nc', 'ICART.1.SEG01.002dg.nc']
Lvariables = ['MRC','COT','O3T','O3_PROD','O3_LOSS','CO_PROD','CO_LOSS',
              'level','ZTOP', 'longitude','latitude','level_w','time',
              'CO_BUDGET','O3_BUDGET','O3_CHREACLIST','CO_CHREACLIST']

Dvar_input = {'f1':Lvariables, 'f2':Lvariables}

#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path="", removeHALO=True)
````
