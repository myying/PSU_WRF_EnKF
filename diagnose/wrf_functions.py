###WRF utility functions
import numpy as np

def ncread(filename, varname):
  import netCDF4
  f = netCDF4.Dataset(filename)
  dat = f.variables[varname]
  return np.array(dat)

def ncwrite(filename, varname, varout):
  from netCDF4 import Dataset
  f = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
  nm, ni, nj = varout.shape
  mm = f.createDimension('m', nm)
  ii = f.createDimension('i', ni)
  jj = f.createDimension('j', nj)
  dat = f.createVariable(varname, np.float32, ('m', 'i', 'j'))
  for m in range(nm):
    dat[m, :, :] = varout[m, :, :]
  f.close()


def writevar(infile, varname, var):
  import netCDF4
  f = netCDF4.Dataset(infile, 'r+')
  f.variables[varname][:] = var
  f.close()

def getvar(infile, varname):
  Rd = 287.0
  Rv = 461.6
  Cp = 1004.5
  g = 9.81
  svp1 = 0.6112
  svp2 = 17.67
  svp3 = 29.65
  ep2 = Rd/Rv
  ep3 = 0.622
  t0 = 273.15

  if (varname == 'ua'):
    dat = ncread(infile, 'U')
    nt, nz, ny, nx = dat.shape
    var = 0.5*(dat[:, :, :, 0:nx-1] + dat[:, :, :, 1:nx])

  elif (varname == 'va'):
    dat = ncread(infile, 'V')
    nt, nz, ny, nx = dat.shape
    var = 0.5*(dat[:, :, 0:ny-1, :] + dat[:, :, 1:ny, :])

  elif (varname == 'wa'):
    dat = ncread(infile, 'W')
    nt, nz, ny, nx = dat.shape
    var = 0.5*(dat[:, 0:nz-1, :, :] + dat[:, 1:nz, :, :])

  elif (varname == 'pres'):
    var = ncread(infile, 'P') + ncread(infile, 'PB')

  elif (varname == 'slp'):
    var = ncread(infile, 'P')[:, 0, :, :] + ncread(infile, 'PB')[:, 0, :, :]

  elif (varname == 'z'):
    var = (ncread(infile, 'P') + ncread(infile, 'PB')) / g

  elif (varname == 'th'):
    var = ncread(infile, 'T') + 300

  elif (varname == 'tk'):
    th = ncread(infile, 'T') + 300
    p = ncread(infile, 'P') + ncread(infile, 'PB')
    var = th * (p/100000.0) ** (Rd/Cp)

  elif (varname == 'wind'):
    dat = ncread(infile, 'U')
    nt, nz, ny, nx = dat.shape
    u = 0.5*(dat[:, :, :, 0:nx-1] + dat[:, :, :, 1:nx])
    dat = ncread(infile, 'V')
    nt, nz, ny, nx = dat.shape
    v = 0.5*(dat[:, :, 0:ny-1, :] + dat[:, :, 1:ny, :])
    # u = ncread(infile, 'U10')
    # v = ncread(infile, 'V10')
    var = np.sqrt(u**2 + v**2)

  else:
    var = ncread(infile, varname)

  return var

