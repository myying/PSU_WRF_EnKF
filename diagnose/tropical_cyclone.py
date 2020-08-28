###Utility function for tropical cyclone diagnostics
import numpy as np
import util

def find_center(pp):
  nj, ni = pp.shape
  min_pp = 1e10
  min_i = 0
  min_j = 0
  buff = 20
  for i in range(buff, ni-buff):
    for j in range(buff, nj-buff):
      avg_pp = np.mean(pp[j-buff:j+buff, i-buff:i+buff])
      if(avg_pp < min_pp):
        min_pp = avg_pp
        min_i = i
        min_j = j
  return min_j, min_i

def maximum_wind(wind_speed, min_j, min_i):
  nj, ni = wind_speed.shape
  smth = 1
  wind_speed = util.smooth2d(wind_speed, smth)
  buff = 30
  vmax = 0
  for i in range(max(0, min_i-buff), min(ni, min_i+buff)):
    for j in range(max(0, min_j-buff), min(nj, min_j+buff)):
      if(wind_speed[j, i] > vmax):
        vmax = wind_speed[j, i]
  return vmax
