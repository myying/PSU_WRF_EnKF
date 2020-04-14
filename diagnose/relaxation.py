#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
n = 360
nens = 60
relax_coef = 0.7

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

land = wrf.getvar(workdir+'fort.80011', 'HGT')
varname = 'U10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
nt, ny, nx = xb.shape
u1, v1 = util.optical_flow_HS(xa[0, 0:n, 0:n], xb[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
varname = 'V10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
nt, ny, nx = xb.shape
u2, v2 = util.optical_flow_HS(xa[0, 0:n, 0:n], xb[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
u = 0.5*(u1+u2)
v = 0.5*(v1+v2)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
  # xb_mean = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
  # xa_mean = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  xa_relax = xa.copy()

  # xa_relax = xa_mean + (1-relax_coef)*(xa - xa_mean) + relax_coef*(xb - xb_mean)
  if(xb.ndim==4):
    xa_relax[0, :, 0:n, 0:n] = util.warp(xa[0, :, 0:n, 0:n], -relax_coef*u, -relax_coef*v)
  if(xb.ndim==3):
    xa_relax[0, 0:n, 0:n] = util.warp(xa[0, 0:n, 0:n], -relax_coef*u, -relax_coef*v)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xa_relax)

