#!/usr/bin/env python3
import numpy as np
import sys
import util
import wrf_functions as wrf
import tropical_cyclone as tc

workdir = '/scratch/02135/yingyue/Patricia/'
casename = sys.argv[1] #'multiscale/fc'
tstr = sys.argv[2] #'201510211200'
nens = int(sys.argv[3]) #60
domain_id = int(sys.argv[4])
print(tstr)
if (domain_id==2):
  dx = 3 #km
if (domain_id==3):
  dx = 1
z_coord = np.arange(0, 20.1, 0.5)*1000 #m
dr = 5 #km
r_coord = np.arange(0, 201, dr)
nz = z_coord.size
nr = r_coord.size

Vt = np.zeros((nz, nr, 4, nens))
Vr = np.zeros((nz, nr, 4, nens))
for m in range(nens):
  filename = workdir+casename+'/'+tstr+'/wrfinput_d{:02d}_{:03d}'.format(domain_id, m+1)
  ua = wrf.getvar(filename, 'ua')[0, :, :, :]
  va = wrf.getvar(filename, 'va')[0, :, :, :]
  z = wrf.getvar(filename, 'z')[0, :, :, :]
  u = wrf.vert_interp(ua, z, z_coord)
  v = wrf.vert_interp(va, z, z_coord)
  center = np.load(workdir+casename+'/'+tstr+'/tc_center_ens.npy')
  tc_j = center[0, m]
  tc_i = center[1, m]
  nv, ny, nx = ua.shape
  Vt[:, :, :, m] = 0.0
  Vr[:, :, :, m] = 0.0
  count = np.zeros((nr, 4))
  count[:] = 1e-10
  for i in range(nx):
    for j in range(ny):
      rad = np.sqrt((j-tc_j)**2 + (i-tc_i)**2)
      rad_ind = int(rad*dx/dr)
      if (rad_ind < nr-1 and rad_ind > 0):
        tangent_wind = -u[:, j, i]*(j-tc_j)/rad + v[:, j, i]*(i-tc_i)/rad
        radial_wind = u[:, j, i]*(i-tc_i)/rad + v[:, j, i]*(j-tc_j)/rad
        if (j>=tc_j and i>=tc_i):
          quad_id = 0
        if (j>=tc_j and i<tc_i):
          quad_id = 1
        if (j<tc_j and i<tc_i):
          quad_id = 2
        if (j<tc_j and i>=tc_i):
          quad_id = 3
        Vt[:, rad_ind, quad_id, m] += tangent_wind
        Vr[:, rad_ind, quad_id, m] += radial_wind
        count[rad_ind, quad_id] += 1
  for k in range(nz):
    Vt[k, :, :, m] = Vt[k, :, :, m]/count
    Vr[k, :, :, m] = Vr[k, :, :, m]/count

np.save(workdir+casename+'/'+tstr+'/Vt_ens', Vt)
np.save(workdir+casename+'/'+tstr+'/Vr_ens', Vr)
