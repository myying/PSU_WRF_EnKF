#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

workdir = sys.argv[1] #'/glade/scratch/mying/Patricia/run/201510230600/wrf_ens/'
filename = sys.argv[2] #'wrfout_d02_2015-10-23_07:00:00'
varname = 'wind'
clevel = np.arange(10, 100, 10)

var = wrf.getvar(workdir+'/'+filename, varname)[0, :, :]

