## Full documentation

````python
def read_netcdf(LnameFiles, Dvar_input, path='.', get_data_only=True, del_empty_dim=True, removeHALO=True):
    """Read a netCDF4 Meso-NH file
    For each file, call functions to read diachronic or synchronous file

    Parameters
    ----------
    LnameFiles : list of str
        list of Meso-NH netCDF4 file (diachronic or synchronous)

    Dvar_input : Dict{'fileNumber' : 'var_name',('group_name','var_name')}
        where
        'fileNumber' is a str corresponding to 'f' + the file number in LnameFiles (by order)
        'var_name' is the exact str of the netCDF4 variable name
        ('group_name','var_name') is the exact tuple of the (sub-)groups name and the netCDF4 variable name
        e.g. : {'f1':['ZS', 'WT','ni', 'level'],
                'f2':[('/LES_budgets/Cartesian/Not_time_averaged/Not_normalized/cart/',MEAN_TH'),('/Budgets/RI','AVEF')]
                }

    path : str
        unique path of the files

    get_data_only : bool, default: True
        if True,  the function returns Dvar as masked_array (only data)
        if False, the function returns Dvar as netCDF4._netCDF4.Variable

    del_empty_dim : bool, default: True
        if get_data_only=True and del_empty_dim=True, returns Dvar as an array without dimensions with size 1 and 0
        e.g. : an array of dimensions (time_budget, cart_level, cart_nj, cart_ni) with shape (180,1,50,1) is returned (180,50)

    removeHALO : bool, default: True
        if True, remove first and last (NHALO=1) point [1:-1] if get_data_only=True on each
        level, level_w, ni, ni_u, ni_v, nj, nj_u, nj_v dimensions

    Returns
    -------
    Dvar : Dict
        Dvar[ifile]['var_name']                if the group contains only one variable
        Dvar[ifile][('group_name','var_name')] if the group contains more than one variable
    """
````
