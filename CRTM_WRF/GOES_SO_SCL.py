#!/usr/bin/env python

#hroi = 30
hroi_large = 400
hroi_lhydro = 60  #now setting this in enkf.f
thin = 3
thin_d = 1 #thin*thin_d = frequency that data used
yyyy = 2011
mm = 10
i_dd = 13
i_hh = 00
f_dd = 18
f_hh = 00
tint = 10
sat = 'abi_gr' #imgr_g13, abi_gr
#ch = '14'

import numpy as np
import matplotlib.pyplot as plt
import datetime
import random

itime = datetime.datetime(yyyy,mm,i_dd,i_hh,0,0)
ftime = datetime.datetime(yyyy,mm,f_dd,f_hh,0,0)

for dom in range(1,2):
  time = itime
  while time <= ftime:
    filename = '/work/02135/yingyue/DYNAMO/EnKF_OSSE/BT_bin/BT_d0'+str(dom)+'_'+time.strftime('%Y-%m-%d_%H:%M')+'.bin'
    # Load binary file
    if dom == 1:
      xmax = 333 #224 / 199
      ymax = 222 #149 / 149
    if dom == 2:
      xmax = 201 #399 / 201
      ymax = 150 #225 / 150
    if dom == 3:
      xmax = 255 #357 / 255
      ymax = 255 #357 / 255
  
    print filename
    data = np.fromfile(filename,dtype='>f4')
    sim=data[:].reshape(12,ymax,xmax)
    lons = sim[0,:,:]
    lats = sim[1,:,:]
    Tb   = sim[2:,:,:]
  
    file_output = '/work/02135/yingyue/DYNAMO/EnKF_OSSE/BT_SO/BT_d0'+str(dom)+'_'+time.strftime('%Y%m%d%H%M')+'_so'
    f = open(file_output,'w')
    if sat == 'abi_gr':
      chlist = [0,8,9,10] #[8,9,10,14] 
    if sat == 'imgr_g13':
      chlist = [3,4]

    for ich in chlist:
      ch = ich
      if sat == 'abi_gr':
        ch_tb = ch-7
        if ich == 0:
           ch = 8
           ch_tb = 1
      elif sat == 'imgr_g13':
        ch_tb = ch-1
  
      if ich == 0:
         clevs = int(thin*thin_d)
         hroi_hy = hroi_lhydro
         hroi_d = hroi_large
      else:
         clevs = thin
         hroi_hy = hroi_lhydro#hroi
         hroi_d = hroi_large#hroi #hroi 

      for i in range(2,xmax,clevs):
        for j in range(2,ymax,clevs):
          if not (ich == 8 and (i-2)%(int(thin*thin_d)) == 0 and (j-2)%(int(thin*thin_d)) == 0):
            ii = i
            jj = j
            #if ch == 8:
            #  ii = i
            #  jj = j
            #elif ch == 9:
            #  ii = i
            #  jj = j+2
            #elif ch == 10:
            #  ii = i+2
            #  jj = j+2

            if ch == 14:
              error = 5.0
            else:
              error = 3.0
            Tbout = Tb[ch_tb,jj,ii]+random.gauss(0,error)
            if sat == 'abi_gr':
              ##setting for cloud region
              if Tb[7,jj,ii]<285 and ch == 8: # 285: selecting cloud region 
                if Tbout>160 and Tbout<330:  
                  dataset = time.strftime('%Y%m%d%H%M')+"{0:>12}".format(sat)+"{0:>12}".format(ch)+"{0:12.3f}".format(float(lats[jj,ii]))
                  dataset = dataset+"{0:12.3f}".format(float(lons[jj,ii]))+"{0:12.3f}".format(float(Tbout))+"{0:>12}".format(hroi_hy)+"{0:>12}".format(hroi_d)+"{0:12.3f}".format(error)+'\n'
                  f.write(dataset)
              ##setting for in-between region
              elif (285 <= Tb[7,jj,ii] and Tb[7,jj,ii]<292) and (ch == 10 or ch == 8):
                if Tbout>160 and Tbout<330:
                  dataset = time.strftime('%Y%m%d%H%M')+"{0:>12}".format(sat)+"{0:>12}".format(ch)+"{0:12.3f}".format(float(lats[jj,ii]))
                  dataset = dataset+"{0:12.3f}".format(float(lons[jj,ii]))+"{0:12.3f}".format(float(Tbout))+"{0:>12}".format(int(hroi_hy))+"{0:>12}".format(hroi_d)+"{0:12.3f}".format(error)+'\n'
                  f.write(dataset)
              ##setting for clear region
              elif Tb[7,jj,ii] >= 292 and (ch == 10 or ch == 9 or ch == 8):
                if Tbout>160 and Tbout<330:
                  dataset = time.strftime('%Y%m%d%H%M')+"{0:>12}".format(sat)+"{0:>12}".format(ch)+"{0:12.3f}".format(float(lats[jj,ii]))
                  dataset = dataset+"{0:12.3f}".format(float(lons[jj,ii]))+"{0:12.3f}".format(float(Tbout))+"{0:>12}".format(int(hroi_hy))+"{0:>12}".format(hroi_d)+"{0:12.3f}".format(error)+'\n'
                  f.write(dataset)

    time = time + datetime.timedelta(minutes = tint) 


