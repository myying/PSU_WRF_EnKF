#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt
import wrf_functions as wrf
import cartopy.crs as ccrs

workdir = '/glade/scratch/mying/Patricia/test_201510230600/enkf_multiscale/'
varname = 'U10'
m = 1
clevel = np.arange(-30, 32, 2)
center_lat = 16.5
center_lon = -105.4
lat = wrf.getvar(workdir+'fort.80011', 'XLAT')[0, :, :]
lon = wrf.getvar(workdir+'fort.80011', 'XLONG')[0, :, :]

var1 = wrf.getvar(workdir+'scale1/fort.{}'.format(m+50010), varname)[0, :, :]
var1[np.where(var1>30)] = 30
var1[np.where(var1<-30)] = -30

plt.switch_backend('Agg')
plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

c = ax.contourf(lon, lat, var1, clevel, cmap='seismic')
cbar = plt.colorbar(c)
cbar.ax.tick_params(labelsize=20)
ax.plot(center_lon, center_lat, marker='+', markersize=20, color='k')
ax.set_aspect('equal', 'box')
ax.set_extent([center_lon-2.5, center_lon+2.5, center_lat-2.5, center_lat+2.5])
ax.coastlines('50m')

plt.savefig("1.png", dpi=200)
plt.close()
