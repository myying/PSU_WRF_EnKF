#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

workdir = '/glade/scratch/mying/Patricia/'
casename = sys.argv[1] #'multiscale/fc'
tstr = sys.argv[2] #'201510211200'
nens = int(sys.argv[3]) #60
domain_id = int(sys.argv[4])

if (domain_id==2):
  dx = 3 #km
if (domain_id==3):
  dx = 1
z_coord = np.arange(0, 20.1, 0.5)*1000 #m
dr = 5 #km
r_coord = np.arange(0, 201, dr)

plt.switch_backend('Agg')

plt.figure(figsize=(5,3))

Vt = np.load(workdir+casename+'/'+tstr+'/Vt_ens.npy')
Vr = np.load(workdir+casename+'/'+tstr+'/Vr_ens.npy')
nz, nr, nq, nens = Vt.shape

Vt_mean = np.mean(np.mean(Vt[:, :, :, 0:nens], axis=2), axis=2)
Vr_mean = np.mean(np.mean(Vr[:, :, :, 0:nens], axis=2), axis=2)
ax = plt.subplot(111)
c = ax.contourf(Vt_mean, np.arange(10, 70, 5), cmap='rainbow')
plt.colorbar(c)
ax.set_xticks(np.arange(0, nr, 5))
ax.set_xticklabels(r_coord[::5])
ax.set_yticks(np.arange(0, nz, 10))
ax.set_yticklabels(z_coord[::10]/1000)
ax.set_xlim(0, 30)
ax.set_ylim(0, 36)

plt.savefig(workdir+casename+'/'+tstr+'/VtVr_mean.png', dpi=100)
