#!/usr/bin/env python3
import numpy as np
import sys
import util
import wrf_functions as wrf
import tropical_cyclone as tc

workdir = '/glade/scratch/mying/Patricia/'
casename = sys.argv[1] #'multiscale/fc'
tstr = sys.argv[2] #'201510211200'
nens = int(sys.argv[3]) #60
domain_id = 2
print(tstr)

tc_center = np.zeros((2, nens))
tc_latlon = np.zeros((2, nens))
tc_vmax = np.zeros(nens)
tc_pmin = np.zeros(nens)
for m in range(nens):
  filename = workdir+casename+'/'+tstr+'/wrfinput_d{:02d}_{:03d}'.format(domain_id, m+1)
  lat = wrf.getvar(filename, 'XLAT')[0, :, :]
  lon = wrf.getvar(filename, 'XLONG')[0, :, :]
  wind_speed = wrf.getvar(filename, 'wind')[0, 0, :, :]
  p_pert = wrf.getvar(filename, 'MU')[0, :, :]
  slp = wrf.getvar(filename, 'PSFC')[0, :, :]/100
  j, i = tc.find_center(p_pert)
  tc_center[0, m] = j
  tc_center[1, m] = i
  tc_latlon[0, m] = lat[j, i]
  tc_latlon[1, m] = lon[j, i]
  tc_pmin[m] = slp[j, i]
  tc_vmax[m] = tc.maximum_wind(wind_speed, j, i)

np.save(workdir+casename+'/'+tstr+'/tc_center_ens', tc_center)
np.save(workdir+casename+'/'+tstr+'/tc_latlon_ens', tc_latlon)
np.save(workdir+casename+'/'+tstr+'/tc_vmax_ens', tc_vmax)
np.save(workdir+casename+'/'+tstr+'/tc_pmin_ens', tc_pmin)
