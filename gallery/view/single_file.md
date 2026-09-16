## Read single file

````python
from read_MNHfile import read_netcdf

#  List of Meso-NH files present in the path
path=""
LnameFiles = ['DUST7.1.SEG02.004.nc']

Dvar_input = {'f1':['ZS', 'UT','VT', 'WT','THT',
      'DSTM03T','DSTM33T','DSTM02T','DSTM32T','DSTM01T','DSTM31T','F_DST001P1','F_DST002P1','F_DST003P1',
      'latitude','longitude','level',
      'INPRR','ACPRR','PABST','RCT','RVT','RRT','LSTHM']}

#  Read the variables in the files
Dvar = {}
Dvar = read_netcdf(LnameFiles, Dvar_input, path=path)
````
