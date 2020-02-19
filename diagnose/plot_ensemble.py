#!/usr/bin/env python3
import numpy as np
import util
# import matplotlib.pyplot as plt
import wrf_functions as wrf
# import sys
# import cartopy.crs as ccrs

workdir = '/glade/scratch/mying/Patricia/control/run/201510212100/wrf_ens/'
filename = 'wrfout_d02_2015-10-21_21:00:00'
#workdir1 = '/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'
varname = 'wind'
# m = int(sys.argv[1])
# clevel = np.arange(10, 100, 10)
# center_lat = 12.9
# center_lon = -96.9
# lat = wrf.getvar(workdir+'001/'+filename, 'XLAT')[0, :, :]
# lon = wrf.getvar(workdir+'001/'+filename, 'XLONG')[0, :, :]
# plt.switch_backend('Agg')
# plt.figure(figsize=(6, 5))

nens = 60
varout = np.zeros((nens, 360, 360))
for m in range(nens):
  varout[m, :, :] = wrf.getvar(workdir+'{:03d}/'.format(m+1)+filename, varname)[0, :, :]
  print(np.max(varout[m, :, :]))
wrf.ncwrite(varname+'_ens.nc', varname, varout)

# for i in range(1):
  # if(i==0):
    # var = wrf.getvar(workdir+'{:03d}/'.format(m)+filename, varname)[0, :, :]
  # if(i==1):
  #   var = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)[0, :, :]
  # if(i==2):
  #   var = wrf.getvar(workdir1+'fort.{}'.format(m+90010), varname)[0, :, :]
  # var[np.where(var>90)] = 90.

#   ax = plt.subplot(1,1,i+1,projection=ccrs.PlateCarree())
#   c = ax.contourf(lon, lat, var, clevel, cmap='Reds')
#   plt.colorbar(c)
#   ax.set_extent([center_lon-4, center_lon+4, center_lat-4, center_lat+4])
#   ax.coastlines('50m')

# plt.tight_layout()
# plt.savefig("{:03d}.png".format(m), dpi=200)
