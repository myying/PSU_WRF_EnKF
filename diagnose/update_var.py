#!/usr/bin/env python3
import numpy as np
import util
import wrf_functions as wrf
import sys

source_file = sys.argv[1]
target_file = sys.argv[2]
varname = sys.argv[3]

varin = wrf.getvar(source_file, varname)
wrf.writevar(target_file, varname, varin)
