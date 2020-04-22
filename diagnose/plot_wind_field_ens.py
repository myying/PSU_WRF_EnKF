#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset
import util
import os
import wrf_functions as wrf
import tropical_cyclone as tc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

workdir = '/glade/scratch/mying/Patricia/'
tstr = '201510220500'
casename = 'control'
clevels = np.arange(0, 50, 5)
nens = 48

plt.switch_backend('Agg')
plt.figure(figsize=(20, 20))

###observation: samuri retrieval
# f = Dataset('/glade/work/mying/data/Patricia/samurai_analyses/patricia_1800UTC_22_Oct_pass1_cartesian_wind.nc')
# lon_obs = np.array(f['longitude'])
# lat_obs = np.array(f['latitude'])
# wspd_obs = np.array(f['WSPD'])[0, 0, :, :]
# ax = plt.subplot(7,7,1,projection=ccrs.PlateCarree())
# ax.contourf(lon_obs, lat_obs, wspd_obs, clevels, cmap='jet')

# f = open('/glade/work/mying/data/Patricia/airborne_radar/gridded/20151022I1radar.dat', 'r')
# print(f.read())

###from model forecast
for m in range(nens):
  # filename = workdir+casename+'/forecast_201510212100/{:03d}'.format(m+1)+'/wrfout_d03_'+util.wrf_time_string(tstr)
  filename = workdir+casename+'/fc/'+tstr+'/wrfinput_d02_{:03d}'.format(m+1)
  hsize = int(200/3)
  lv = 0
  lat = wrf.getvar(filename, 'XLAT')[0, 180-hsize:180+hsize, 180-hsize:180+hsize]
  lon = wrf.getvar(filename, 'XLONG')[0, 180-hsize:180+hsize, 180-hsize:180+hsize]
  wspd = wrf.getvar(filename, 'wind')[0, lv, 180-hsize:180+hsize, 180-hsize:180+hsize]
  ax = plt.subplot(7,7,m+2,projection=ccrs.PlateCarree())
  ax.contourf(lon, lat, wspd, clevels, cmap='jet')

plt.tight_layout()
plt.savefig(casename+'_windspeed.png', dpi=100)
