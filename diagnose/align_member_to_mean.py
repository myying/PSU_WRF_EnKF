#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
n = 360
nens = 60

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

land = wrf.getvar(workdir+'fort.80011', 'HGT')
varname = 'U10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
xb_mean = wrf.getvar(workdir+'fort.{}'.format(80011+nens), varname)
xa_mean = wrf.getvar(workdir+'fort.{}'.format(90011+nens), varname)
nt, ny, nx = xb.shape
ub1, vb1 = util.optical_flow_HS(xb[0, 0:n, 0:n], xb_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
ua1, va1 = util.optical_flow_HS(xa[0, 0:n, 0:n], xa_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
varname = 'V10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
xb_mean = wrf.getvar(workdir+'fort.{}'.format(80011+nens), varname)
xa_mean = wrf.getvar(workdir+'fort.{}'.format(90011+nens), varname)
nt, ny, nx = xb.shape
ub2, vb2 = util.optical_flow_HS(xb[0, 0:n, 0:n], xb_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
ua2, va2 = util.optical_flow_HS(xa[0, 0:n, 0:n], xa_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
ub = 0.5*(ub1+ub2)
vb = 0.5*(vb1+vb2)
ua = 0.5*(ua1+ua2)
va = 0.5*(va1+va2)

##save a copy of u,v
np.save('ub_{:03d}'.format(m), ub)
np.save('vb_{:03d}'.format(m), vb)
np.save('ua_{:03d}'.format(m), ua)
np.save('va_{:03d}'.format(m), va)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  ##align member to mean
  if(xb.ndim==4):
    xb[0, :, 0:n, 0:n] = util.warp(xb[0, :, 0:n, 0:n], -ub, -vb)
    xa[0, :, 0:n, 0:n] = util.warp(xa[0, :, 0:n, 0:n], -ua, -va)
  if(xb.ndim==3):
    xb[0, 0:n, 0:n] = util.warp(xb[0, 0:n, 0:n], -ub, -vb)
    xa[0, 0:n, 0:n] = util.warp(xa[0, 0:n, 0:n], -ua, -va)

  wrf.writevar(workdir+'fort.{}'.format(m+30010), varname, xb)
  wrf.writevar(workdir+'fort.{}'.format(m+40010), varname, xa)

