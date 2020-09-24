#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

m = int(sys.argv[1])
n = 360
nens = 60
relax_coef_disp = 0.7
relax_coef_pert = 0.7

workdir = './' #'/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'

ub = np.load('ub_{:03d}.npy'.format(m))
vb = np.load('vb_{:03d}.npy'.format(m))
ua = np.load('ua_{:03d}.npy'.format(m))
va = np.load('va_{:03d}.npy'.format(m))

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
  xb_mean = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
  xa_mean = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)

  xb_align = wrf.getvar(workdir+'fort.{}'.format(m+30010), varname)
  xa_align = wrf.getvar(workdir+'fort.{}'.format(m+40010), varname)
  xb_align_mean = wrf.getvar(workdir+'fort.{}'.format(nens+30011), varname)
  xa_align_mean = wrf.getvar(workdir+'fort.{}'.format(nens+40011), varname)

  if(xb.ndim==4):
    nt, nz, ny, nx = xb.shape
  if(xb.ndim==3):
    nt, ny, nx = xb.shape

  ##amplitude part
  xa = xa_mean + (1-relax_coef_pert)*(xa_align-xa_align_mean) + relax_coef_pert*(xb_align-xb_align_mean)

  ##find displacement part and relax to prior
  if(xb.ndim==4):
    xa[0, :, 0:n, 0:n] = util.warp(xa[0, :, 0:n, 0:n], (1-relax_coef_disp)*ua, (1-relax_coef_disp)*va)
    xa[0, :, 0:n, 0:n] = util.warp(xa[0, :, 0:n, 0:n], relax_coef_disp*ub, relax_coef_disp*vb)
  if(xb.ndim==3):
    xa[0, 0:n, 0:n] = util.warp(xa[0, 0:n, 0:n], (1-relax_coef_disp)*ua, (1-relax_coef_disp)*va)
    xa[0, 0:n, 0:n] = util.warp(xa[0, 0:n, 0:n], relax_coef_disp*ub, relax_coef_disp*vb)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xa)

