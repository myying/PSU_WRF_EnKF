#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
current_scale = int(sys.argv[2])
n = 300

workdir = '/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

# varname = 'P'
# xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
# xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
# u, v = util.optical_flow_HS(np.mean(xb[0, :, 0:n, 0:n], axis=0), np.mean(xa[0, :, 0:n, 0:n], axis=0), 5)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
  xs = wrf.getvar(workdir+'fort.{}'.format(m+60010), varname)
  xin = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  xsw = xs.copy()
  xout = xin.copy()

  if(current_scale<3):
    if(xb.ndim==4):
      u, v = util.optical_flow_HS(np.mean(xb[0, :, 0:n, 0:n], axis=0), np.mean(xa[0, :, 0:n, 0:n], axis=0), 5)
      xsw[0, :, 0:n, 0:n] = util.warp(xs[0, :, 0:n, 0:n], -u, -v)
    if(xb.ndim==3):
      u, v = util.optical_flow_HS(xb[0, 0:n, 0:n], xa[0, 0:n, 0:n], 5)
      xsw[0, 0:n, 0:n] = util.warp(xs[0, 0:n, 0:n], -u, -v)
    xout = xin + xa - xb + xsw - xs
  else:
    xout = xin + xa - xb

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xout)

