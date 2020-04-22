#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
n = 360
nens = 60
relax_coef = 0.3

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

land = wrf.getvar(workdir+'fort.80011', 'HGT')
##displacement from prior mean to post mean
varname = 'U10'
xb_mean = wrf.getvar(workdir+'fort.{}'.format(80011+nens), varname)
xa_mean = wrf.getvar(workdir+'fort.{}'.format(90011+nens), varname)
u1, v1 = util.optical_flow_HS(xb_mean[0, 0:n, 0:n], xa_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
varname = 'V10'
xb_mean = wrf.getvar(workdir+'fort.{}'.format(80011+nens), varname)
xa_mean = wrf.getvar(workdir+'fort.{}'.format(90011+nens), varname)
u2, v2 = util.optical_flow_HS(xb_mean[0, 0:n, 0:n], xa_mean[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
u_mean = 0.5*(u1+u2)
v_mean = 0.5*(v1+v2)

##read prior member, shift mean position, then find displacement from post member to that
varname = 'U10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xb[0, 0:n, 0:n] = util.warp(xb[0, 0:n, 0:n], -u_mean, -v_mean)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
u1, v1 = util.optical_flow_HS(xa[0, 0:n, 0:n], xb[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
varname = 'V10'
xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
xb[0, 0:n, 0:n] = util.warp(xb[0, 0:n, 0:n], -u_mean, -v_mean)
xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
u2, v2 = util.optical_flow_HS(xa[0, 0:n, 0:n], xb[0, 0:n, 0:n], 5, land[0, 0:n, 0:n])
u = 0.5*(u1+u2)
v = 0.5*(v1+v2)

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
  # xb_mean = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
  # xa_mean = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)

  # xa = xa_mean + (1-relax_coef)*(xa-xa_mean) + relax_coef*(xb-xb_mean)

  if(xb.ndim==4):
    xa[0, :, 0:n, 0:n] = util.warp(xa[0, :, 0:n, 0:n], -relax_coef*u, -relax_coef*v)
  if(xb.ndim==3):
    xa[0, 0:n, 0:n] = util.warp(xa[0, 0:n, 0:n], -relax_coef*u, -relax_coef*v)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xa)

