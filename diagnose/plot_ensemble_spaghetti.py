#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt
import wrf_functions as wrf

workdir = '/glade/scratch/mying/Patricia_multiscale/run/201510230600/enkf/d02/'
nens = 20

lat = wrf.getvar(workdir+'fort.80011', 'XLAT')[0, :, :]
lon = wrf.getvar(workdir+'fort.80011', 'XLONG')[0, :, :]
ny, nx = lon.shape
var_ens = np.zeros((nens, ny, nx))

for m in range(nens):
  var = wrf.getvar(workdir+'fort.{}'.format(90011+m), 'PSFC')
  var_ens[m, :, :] = var[0, :, :]

plt.switch_backend('Agg')
plt.figure(figsize=(10, 8))
ax = plt.axes()

for m in range(nens):
  ax.contour(lon, lat, var_ens[m, :, :], (96000,), colors='k')

plt.savefig("1.pdf")
