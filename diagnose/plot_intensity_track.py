#!/usr/bin/env python3
import numpy as np
import util
import os
import wrf_functions as wrf
import tropical_cyclone as tc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

diag_dir = '/glade/work/mying/data/Patricia/diag/'
work_dir = '/glade/scratch/mying/Patricia/'
t_start = '201510211200'
nt = 28 #32
casename = 'multiscale'

###observation: tcvitals best track
obs_t_inv = 6
t_obs = np.arange(0, 73, obs_t_inv)
vmax_obs = np.array([40, 50, 60, 75, 90, 115, 150, 180, 185, 180, 110, 50, 25])*0.51444
pmin_obs = np.array([1001, 997, 991, 981, 969, 957, 920, 886, 872, 878, 946, 985, 1000])
# vmax_obs = np.array([  18,   21,   28,    36,    41,    59,    72,    85,    90,    90,    72,    44,    21])
# pmin_obs = np.array([1004,  999,  994,   985,   975,   958,   924,   886,   880,   879,   924,   970,  1000])
lat_obs = np.array([ 12.9, 13.0, 13.3,  14.0,  14.7,  15.1,  15.8,  16.5,  17.2,  18.2,  19.6,  21.4,  23.2])
lon_obs = np.array([-96.9,-98.6,-99.9,-101.7,-103.2,-104.1,-104.9,-105.4,-105.6,-105.3,-104.9,-104.0,-102.3])

tc_lat = np.zeros(nt)
tc_lon = np.zeros(nt)
tc_vmax = np.zeros(nt)
tc_pmin = np.zeros(nt)
tc_lat[:] = np.nan
tc_lon[:] = np.nan
tc_vmax[:] = np.nan
tc_pmin[:] = np.nan
for t in range(0, nt):
  t_str = util.advance_time(t_start, t*60)
  print(t_str)
  filename = work_dir+casename+'/fc/'+t_str+'/wrfinput_d02_mean'
  # if(os.path.exists(filename)):
    #skip
  lat = wrf.getvar(filename, 'XLAT')[0, :, :]
  lon = wrf.getvar(filename, 'XLONG')[0, :, :]
  wind_speed = util.smooth(wrf.getvar(filename, 'wind')[0, 0, :, :], 0)
  p_pert = wrf.getvar(filename, 'P')[0, 0, :, :]
  slp = wrf.getvar(filename, 'PSFC')[0, :, :]/100
  j, i = tc.find_center(p_pert)
  tc_lat[t] = lat[j, i]
  tc_lon[t] = lon[j, i]
  tc_pmin[t] = slp[j, i]
  tc_vmax[t] = tc.maximum_wind(wind_speed, j, i)

plt.switch_backend('Agg')
plt.figure(figsize=(15, 5))

###track plots
ax = plt.subplot(1,3,1,projection=ccrs.PlateCarree())
ax.plot(lon_obs, lat_obs, color='k', marker='o')
ax.plot(tc_lon, tc_lat, 'r')
ax.plot(tc_lon[::obs_t_inv], tc_lat[::obs_t_inv], 'ro')
ax.set_extent([-108, -96, 12, 24])
ax.coastlines(resolution='50m')
ax.set_xticks(np.arange(-108, -95, 2), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(12, 25, 2), crs=ccrs.PlateCarree())

###intensity: wind max and pres min
ax = plt.subplot(1,3,2)
ax.plot(t_obs, vmax_obs, color='k', marker='o')
ax.plot(np.arange(0, nt), tc_vmax, 'r')
ax.plot(np.arange(0, nt, obs_t_inv), tc_vmax[::obs_t_inv], 'ro')
ax.set_ylim(0, 110)
ax.set_xticks(np.arange(0, 73, 12))
ax.set_xticklabels(np.array(['21/12', '22/00', '22/12', '23/00', '23/12', '24/00', '24/12']))

ax = plt.subplot(1,3,3)
ax.plot(t_obs, pmin_obs, color='k', marker='o')
ax.plot(np.arange(0, nt), tc_pmin, 'r')
ax.plot(np.arange(0, nt, obs_t_inv), tc_pmin[::obs_t_inv], 'ro')
ax.set_ylim(860, 1020)
ax.set_xticks(np.arange(0, 73, 12))
ax.set_xticklabels(np.array(['21/12', '22/00', '22/12', '23/00', '23/12', '24/00', '24/12']))

plt.tight_layout()
plt.savefig(casename+'_intensity_track.png', dpi=100)
