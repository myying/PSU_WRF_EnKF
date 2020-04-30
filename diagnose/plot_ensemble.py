#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import util
import sys
import wrf_functions as wrf
import tropical_cyclone as tc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

workdir = '/glade/scratch/mying/Patricia/'
casename = sys.argv[1] #'multiscale'
tstr = sys.argv[2] #'201510211500'
nens = 60
# varname = 'slp'
# clevels = np.arange(97000, 101000, 200)
varname = 'wind'
clevels = np.arange(0, 80, 5)

plt.switch_backend('Agg')
plt.figure(figsize=(15, 9))

###observation: samuri retrieval
# f = Dataset('/glade/work/mying/data/Patricia/samurai_analyses/patricia_1800UTC_22_Oct_pass2_cartesian_wind.nc')
# lon_obs = np.array(f['longitude'])
# lat_obs = np.array(f['latitude'])
# wspd_obs = np.array(f['WSPD'])[0, 3, :, :]
# ax = plt.subplot(3,7,1,projection=ccrs.PlateCarree())
# ax.contourf(lon_obs, lat_obs, wspd_obs, clevels, cmap='jet')

###from model forecast
for m in range(nens):
  # filename = workdir+casename+'/forecast_201510212100/{:03d}'.format(m+1)+'/wrfout_d02_'+util.wrf_time_string(tstr)
  filename = workdir+casename+'/fc/'+tstr+'/wrfinput_d02_{:03d}'.format(m+1)
  dx = 3
  hsize = int(300/dx)
  lv = 7
  lat = wrf.getvar(filename, 'XLAT')[0, 180-hsize:180+hsize, 180-hsize:180+hsize]
  lon = wrf.getvar(filename, 'XLONG')[0, 180-hsize:180+hsize, 180-hsize:180+hsize]
  tmp = wrf.getvar(filename, varname)
  if (tmp.ndim == 3):
    var = tmp[0, 180-hsize:180+hsize, 180-hsize:180+hsize]
  if (tmp.ndim == 4):
    var = tmp[0, lv, 180-hsize:180+hsize, 180-hsize:180+hsize]
  ax = plt.subplot(6,10,m+1,projection=ccrs.PlateCarree())
  ax.contourf(lon, lat, var, clevels, cmap='jet')

plt.tight_layout()
plt.savefig('/glade/work/mying/data/Patricia/diag/'+casename+'_ens_'+varname+'_'+tstr+'.png', dpi=100)
