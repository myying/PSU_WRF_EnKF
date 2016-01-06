#!/usr/bin/env python

#hroi = 30
hroi_large = 400
hroi_lhydro = 60 
thin = 1
thin_d = 1 #thin*thin_d = frequency that data used
yyyy = 2011
mm = 11
i_dd = 01
i_hh = 00
f_dd = 30
f_hh = 21
tint = 180
sat = 'mviriNOM_m07'

import numpy as np
import matplotlib.pyplot as plt
import datetime
import random

itime = datetime.datetime(yyyy,mm,i_dd,i_hh,0,0)
ftime = datetime.datetime(yyyy,mm,f_dd,f_hh,0,0)

for dom in [1]:
  time = itime
  while time <= ftime:
    filename = '/work/02135/yingyue/DYNAMO/EnKF_OSSE/201110120000/truth/BT/BT_d0'+str(dom)+'_'+time.strftime('%Y-%m-%d_%H:%M')+'.bin'
    if dom == 1:
      xmax = 333 
      ymax = 222
  
    print filename
    data = np.fromfile(filename,dtype='>f4')
    sim=data[:].reshape(4,ymax,xmax)
    lons = sim[0,:,:]
    lats = sim[1,:,:]
    Tb   = sim[2:,:,:]
  
    file_output = '/work/02135/yingyue/DYNAMO/EnKF_OSSE/Meteosat7_SO/BT_d0'+str(dom)+'_'+time.strftime('%Y%m%d%H%M')+'_so'
    f = open(file_output,'w')

    chlist = [3]

    for ch in chlist:
      ich = ch-2
  
      hroi_hy = hroi_lhydro
      hroi_d = hroi_large

      for i in range(5,xmax-4,thin):
        for j in range(5,ymax-4,thin):
           error = 3.0
           Tbout = Tb[ich,j-1,i-1]+random.gauss(0,error)
           if Tbout>160 and Tbout<330:  
              dataset = time.strftime('%Y%m%d%H%M')+"{0:>15}".format(sat)+"{0:>12}".format(ch)+"{0:12.3f}".format(float(lats[j-1,i-1]))
              dataset = dataset+"{0:12.3f}".format(float(lons[j-1,i-1]))+"{0:12.3f}".format(float(Tbout))+"{0:>12}".format(hroi_hy)+"{0:>12}".format(hroi_d)+"{0:12.3f}".format(error)+'\n'
              f.write(dataset)

    time = time + datetime.timedelta(minutes = tint) 


