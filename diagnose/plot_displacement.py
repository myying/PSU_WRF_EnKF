#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt
import wrf_functions as wrf
import cartopy.crs as ccrs

workdir = '/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'
varname = 'U10'
m = 26
clevel = np.arange(-40, 40, 2)
center_lat = 16.5
center_lon = -105.4
lat = wrf.getvar(workdir+'fort.80011', 'XLAT')[0, :, :]
lon = wrf.getvar(workdir+'fort.80011', 'XLONG')[0, :, :]

var1 = wrf.getvar(workdir+'scale1/fort.{}'.format(m+50010), varname)[0, :, :]
var2 = wrf.getvar(workdir+'scale1/fort.{}'.format(m+70010), varname)[0, :, :]
u, v = util.optical_flow_HS(var1, var2, 5)

for factor in np.arange(0.0, 1.1, 0.1):
  plt.switch_backend('Agg')
  plt.figure(figsize=(8, 6))
  ax = plt.axes(projection=ccrs.PlateCarree())

  var3 = util.warp(var1, -factor*u, -factor*v)

  c = ax.contourf(lon, lat, var3, clevel, cmap='seismic')
  plt.colorbar(c)
  ax.quiver(lon[::10, ::10], lat[::10, ::10], factor*v[::10, ::10], factor*u[::10, ::10], scale=300, headwidth=6, headlength=8, headaxislength=6, linewidths=3)
  ax.set_aspect('equal', 'box')
  ax.set_extent([center_lon-4, center_lon+4, center_lat-4, center_lat+4])
  ax.coastlines('50m')

  plt.savefig("{:4.2f}.png".format(factor), dpi=200)
  plt.close()
