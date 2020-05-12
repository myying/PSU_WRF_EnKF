#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
current_scale = int(sys.argv[2])
n = 360
nens = 60
relax_coef = 0.7

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

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
  xb_mean = wrf.getvar(workdir+'fort.{}'.format(50011+nens), varname)
  xa_mean = wrf.getvar(workdir+'fort.{}'.format(70011+nens), varname)
  xin = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  xbw = xb.copy()
  xout = xin.copy()

  if(current_scale<3):
    if(xb.ndim==4):
      xout[0, :, 0:n, 0:n] = util.warp(xin[0, :, 0:n, 0:n], -u, -v)
      xbw[0, :, 0:n, 0:n] = util.warp(xb[0, :, 0:n, 0:n], -u, -v)
    if(xb.ndim==3):
      xout[0, 0:n, 0:n] = util.warp(xin[0, 0:n, 0:n], -u, -v)
      xbw[0, 0:n, 0:n] = util.warp(xb[0, 0:n, 0:n], -u, -v)
    xout = xout + xa - xbw
  else:
    xout = xin + xa - xb

  ##relaxation
  xa_pert = xa - xa_mean
  xb_pert = xb - xb_mean
  if(current_scale<3):
    if(xb.ndim==4):
      xb_pert[0, :, 0:n, 0:n] = util.warp(xb_pert[0, :, 0:n, 0:n], -u, -v)
    if(xb.ndim==3):
      xb_pert[0, 0:n, 0:n] = util.warp(xb_pert[0, 0:n, 0:n], -u, -v)
  xout = xout + relax_coef*(xb_pert - xa_pert)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xout)

