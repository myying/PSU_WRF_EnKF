#!/usr/bin/env python3
import numpy as np
import util
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

workdir = '/glade/scratch/mying/Patricia/'
casename = 'control'
t_start = '201510211200'
nt = 37
nens = 60

###observation: tcvitals best track
obs_t_inv = 6
t_obs = np.arange(0, 73, obs_t_inv)
vmax_obs = np.array([40, 50, 60, 75, 90, 115, 150, 180, 185, 180, 110, 50, 25])*0.51444
pmin_obs = np.array([1001, 997, 991, 981, 969, 957, 920, 886, 872, 878, 946, 985, 1000])
lat_obs = np.array([ 12.9, 13.0, 13.3,  14.0,  14.7,  15.1,  15.8,  16.5,  17.2,  18.2,  19.6,  21.4,  23.2])
lon_obs = np.array([-96.9,-98.6,-99.9,-101.7,-103.2,-104.1,-104.9,-105.4,-105.6,-105.3,-104.9,-104.0,-102.3])

tc_lat = np.zeros((nt, nens))
tc_lon = np.zeros((nt, nens))
tc_vmax = np.zeros((nt, nens))
tc_pmin = np.zeros((nt, nens))
tc_lat[:, :] = np.nan
tc_lon[:, :] = np.nan
tc_vmax[:, :] = np.nan
tc_pmin[:, :] = np.nan
for t in range(0, nt):
  t_str = util.advance_time(t_start, t*60)
  tdir = workdir+casename+'/fc/'+t_str
  if (os.path.exists(tdir+'/tc_center_ens.npy')):
    latlon = np.load(tdir+'/tc_latlon_ens.npy')
    tc_lat[t, :] = latlon[0, :]
    tc_lon[t, :] = latlon[1, :]
    tc_vmax[t, :] = np.load(tdir+'/tc_vmax_ens.npy')
    tc_pmin[t, :] = np.load(tdir+'/tc_pmin_ens.npy')

plt.switch_backend('Agg')
plt.figure(figsize=(15, 5))
cmap = [plt.cm.jet(m) for m in np.linspace(0, 1, nens)]

###track plots
ax = plt.subplot(1,3,1,projection=ccrs.PlateCarree())
for m in range(nens):
  ax.plot(tc_lon[:, m], tc_lat[:, m], color=cmap[m][0:3]) #, marker='.')
ax.plot(lon_obs, lat_obs, color='k', marker='o')
ax.set_extent([-108, -96, 12, 24])
ax.coastlines(resolution='50m')
ax.set_xticks(np.arange(-108, -95, 2), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(12, 25, 2), crs=ccrs.PlateCarree())

###intensity: wind max and pres min
ax = plt.subplot(1,3,2)
for m in range(nens):
  ax.plot(np.arange(0, nt), tc_vmax[:, m], color=cmap[m][0:3]) #, marker='.')
ax.plot(t_obs, vmax_obs, color='k', marker='o')
ax.set_ylim(0, 110)
ax.set_xticks(np.arange(0, 73, 12))
ax.set_xticklabels(np.array(['21/12', '22/00', '22/12', '23/00', '23/12', '24/00', '24/12']))

ax = plt.subplot(1,3,3)
for m in range(nens):
  ax.plot(np.arange(0, nt), tc_pmin[:, m], color=cmap[m][0:3]) #, marker='.')
ax.plot(t_obs, pmin_obs, color='k', marker='o')
ax.set_ylim(860, 1020)
ax.set_xticks(np.arange(0, 73, 12))
ax.set_xticklabels(np.array(['21/12', '22/00', '22/12', '23/00', '23/12', '24/00', '24/12']))

plt.tight_layout()
plt.savefig(casename+'_intensity_track_ens.png', dpi=100)
