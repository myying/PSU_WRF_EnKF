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

for varname in ('U', 'V', 'W', 'P', 'PH', 'T', 'MU', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'U10', 'V10', 'T2', 'Q2', 'PSFC'):
  xb = wrf.getvar(workdir+'fort.{}'.format(m+80010), varname)
  xa = wrf.getvar(workdir+'fort.{}'.format(m+90010), varname)
  xb_mean = wrf.getvar(workdir+'fort.{}'.format(nens+80011), varname)
  xa_mean = wrf.getvar(workdir+'fort.{}'.format(nens+90011), varname)

  xa = xa_mean + (1-relax_coef)*(xa-xa_mean) + relax_coef*(xb-xb_mean)

  wrf.writevar(workdir+'fort.{}'.format(m+90010), varname, xa)

