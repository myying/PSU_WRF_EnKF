#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt
import wrf_functions as wrf
import cartopy.crs as ccrs

workdir = '/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'
varname = 'U10'
m = 26
clevel = np.arange(-75, 75, 5)
center_lat = 16.5
center_lon = -105.4
lat = wrf.getvar(workdir+'fort.80011', 'XLAT')[0, :, :]
lon = wrf.getvar(workdir+'fort.80011', 'XLONG')[0, :, :]

var1 = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)[0, :, :]

plt.switch_backend('Agg')
plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

c = ax.contourf(lon, lat, var1, clevel, cmap='seismic')
plt.colorbar(c)
ax.set_aspect('equal', 'box')
ax.set_extent([center_lon-4, center_lon+4, center_lat-4, center_lat+4])
ax.coastlines('50m')

plt.savefig("1.png", dpi=200)
plt.close()
