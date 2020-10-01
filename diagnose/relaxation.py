#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
nens = int(sys.argv[2]) ##60
relax_coef = float(sys.argv[3]) ##0.25
n = 360

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

land = wrf.getvar(workdir+'fort.80011', 'HGT')
varname = 'U10'
xb = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)
nt, ny, nx = xb.shape
u1, v1 = util.optical_flow_HS(xb[0, 0:n, 0:n], xa[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
varname = 'V10'
xb = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)
nt, ny, nx = xb.shape
u2, v2 = util.optical_flow_HS(xb[0, 0:n, 0:n], xa[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
u = 0.5*(u1+u2)
v = 0.5*(v1+v2)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
  xb_mean = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
  xa_mean = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)

  xbw = xb.copy()
  xbmw = xb_mean.copy()
  if(xb.ndim==4):
    xbw[0, :, 0:n, 0:n] = util.warp(xb[0, :, 0:n, 0:n], -u, -v)
    xbmw[0, :, 0:n, 0:n] = util.warp(xb_mean[0, :, 0:n, 0:n], -u, -v)
  if(xb.ndim==3):
    xbw[0, 0:n, 0:n] = util.warp(xb[0, 0:n, 0:n], -u, -v)
    xbmw[0, 0:n, 0:n] = util.warp(xb_mean[0, 0:n, 0:n], -u, -v)

  xa = xa_mean + (1-relax_coef)*(xa-xa_mean) + relax_coef*(xbw-xbmw)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xa)

