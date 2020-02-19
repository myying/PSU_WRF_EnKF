#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
current_scale = int(sys.argv[2])
n = 360

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

# varname = 'P'
# xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
# xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
# u, v = util.optical_flow_HS(np.mean(xb[0, :, 0:n, 0:n], axis=0), np.mean(xa[0, :, 0:n, 0:n], axis=0), 5)

if(current_scale<3):
  land = wrf.getvar(workdir+'fort.80011', 'HGT')
  varname = 'U10'
  xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
  nt, ny, nx = xb.shape
  u1, v1 = util.optical_flow_HS(xb[0, 0:n, 0:n], xa[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
  varname = 'V10'
  xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
  nt, ny, nx = xb.shape
  u2, v2 = util.optical_flow_HS(xb[0, 0:n, 0:n], xa[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
  u = 0.5*(u1+u2)
  v = 0.5*(v1+v2)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+50010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+70010), varname)
  # xs = wrf.getvar(workdir+'fort.{}'.format(m+60010), varname)
  xin = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  # xsw = xs.copy()
  xout = xin.copy()

  if(current_scale<3):
    if(xb.ndim==4):
      # xsw[0, :, 0:n, 0:n] = util.warp(xs[0, :, 0:n, 0:n], -u, -v)
      xout[0, :, 0:n, 0:n] = util.warp(xin[0, :, 0:n, 0:n], -u, -v)
    if(xb.ndim==3):
      # xsw[0, 0:n, 0:n] = util.warp(xs[0, 0:n, 0:n], -u, -v)
      xout[0, 0:n, 0:n] = util.warp(xin[0, 0:n, 0:n], -u, -v)
    # xout = xin + xa - xb + xsw - xs
  # else:
    # xout = xin + xa - xb

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xout)

