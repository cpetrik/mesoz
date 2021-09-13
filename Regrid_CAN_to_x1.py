#!/usr/bin/env python 

import os
import numpy as np
import xarray as xr
import xesmf as xe
from glob import glob

# input dataset
in_files = glob('/home/cmp/Public/CMIP6_output_native_grid/CAN/*os*.nc')

grid_out = xr.Dataset({'lat': (['lat'], np.arange(89.5, -90, -1)),
                     'lon': (['lon'], np.arange(-179.5, 180, 1))})

for f in in_files:
    ds = xr.open_dataset(f)
    
    file=str.replace(f,'/home/cmp/Public/CMIP6_output_native_grid/CAN/','')
    print('reading file: ' + file)
    var=str.split(file,'_')[0]
    dvar = ds[var]
    
    ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
    
    # generate regridder
    regridder = xe.Regridder(ds, grid_out, 'bilinear', periodic=True, ignore_degenerate=True, reuse_weights=True)
    dvar_out = regridder(dvar)
    
    # generate dataset
    dvar_out.attrs.update(dvar.attrs)
    dout = dvar_out.to_dataset()
    
    # write file
    file_out = str.replace(file,'.nc','_onedeg.nc')
    print('writing file: ' + file_out)
    dout.to_netcdf('output/'+file_out, mode='w')
    