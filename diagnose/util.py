###Utility function for handling qgmodel fields
import numpy as np

def advance_time(timestart, incr_minute):
  import datetime
  ccyy = int(timestart[0:4])
  mm = int(timestart[4:6])
  dd = int(timestart[6:8])
  hh = int(timestart[8:10])
  ii = int(timestart[10:12])
  t1 = datetime.datetime(ccyy, mm, dd, hh, ii, 0)
  t2 = t1 + datetime.timedelta(minutes=incr_minute)
  ccyy = t2.year
  mm = t2.month
  dd = t2.day
  hh = t2.hour
  ii = t2.minute
  timeout = "{:04d}{:02d}{:02d}{:02d}{:02d}".format(ccyy, mm, dd, hh, ii)
  return timeout

def wrf_time_string(tstr):
  ccyy = tstr[0:4]
  mm = tstr[4:6]
  dd = tstr[6:8]
  hh = tstr[8:10]
  ii = tstr[10:12]
  wrf_tstr = ccyy+'-'+mm+'-'+dd+'_'+hh+':'+ii+':00'
  return wrf_tstr

###spatial operation
def deriv_x(f):
  fx = 0.5*(np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0))
  return fx

def deriv_y(f):
  fy = 0.5*(np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1))
  return fy

def deriv_t(f):
  ft = f.copy()
  nx, ny, nt = f.shape
  ft[:, :, 1:nt] = 0.5*(np.roll(f[:, :, 1:nt], -1, axis=2) - np.roll(f[:, :, 1:nt], 1, axis=2))
  ft[:, :, 0] = f[:, :, 1] - f[:, :, 0]
  ft[:, :, nt-1] = f[:, :, nt-1] - f[:, :, nt-2]
  return ft

def deriv_xy(f):
  fxy = 0.25*(np.roll(np.roll(f, -1, axis=0), -1, axis=1) + np.roll(np.roll(f, 1, axis=0), 1, axis=1) - np.roll(np.roll(f, -1, axis=0), 1, axis=1) - np.roll(np.roll(f, 1, axis=0), -1, axis=1))
  return fxy

def laplacian(f):
  del2f = (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0) + np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1))/6 + (np.roll(np.roll(f, -1, axis=1), -1, axis=0) + np.roll(np.roll(f, -1, axis=1), 1, axis=0) + np.roll(np.roll(f, 1, axis=1), -1, axis=0) + np.roll(np.roll(f, 1, axis=1), 1, axis=0))/12 - f
  return del2f

def deriv_xx(f):
  fxx = (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0))/2 - f
  return fxx

def deriv_yy(f):
  fyy = (np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1))/2 - f
  return fyy

def laplacian_spacetime(f):
  del2f = (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0) + np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1) + np.roll(f, -1, axis=2) + np.roll(f, 1, axis=2))/14 + (
      np.roll(np.roll(f, -1, axis=1), -1, axis=0) + np.roll(np.roll(f, -1, axis=1), 1, axis=0) + np.roll(np.roll(f, 1, axis=1), -1, axis=0) + np.roll(np.roll(f, 1, axis=1), 1, axis=0) +
      np.roll(np.roll(f, -1, axis=1), -1, axis=2) + np.roll(np.roll(f, -1, axis=1), 1, axis=2) + np.roll(np.roll(f, 1, axis=1), -1, axis=2) + np.roll(np.roll(f, 1, axis=1), 1, axis=2) +
      np.roll(np.roll(f, -1, axis=2), -1, axis=0) + np.roll(np.roll(f, -1, axis=2), 1, axis=0) + np.roll(np.roll(f, 1, axis=2), -1, axis=0) + np.roll(np.roll(f, 1, axis=2), 1, axis=0))/28 + (
      np.roll(np.roll(np.roll(f, 1, axis=0), 1, axis=1), 1, axis=2) +
      np.roll(np.roll(np.roll(f, 1, axis=0), 1, axis=1), -1, axis=2) +
      np.roll(np.roll(np.roll(f, 1, axis=0), -1, axis=1), 1, axis=2) +
      np.roll(np.roll(np.roll(f, -1, axis=0), 1, axis=1), 1, axis=2) +
      np.roll(np.roll(np.roll(f, -1, axis=0), -1, axis=1), 1, axis=2) +
      np.roll(np.roll(np.roll(f, 1, axis=0), -1, axis=1), -1, axis=2) +
      np.roll(np.roll(np.roll(f, -1, axis=0), 1, axis=1), -1, axis=2) +
      np.roll(np.roll(np.roll(f, -1, axis=0), -1, axis=1), -1, axis=2))/56 - f
  return del2f

###compute optical flow in observed space using hierarchical HS algorithm
def optical_flow_HS(Im1, Im2, nlevel, mask):
  ni, nj = Im1.shape
  u = np.zeros((ni, nj))
  v = np.zeros((ni, nj))
  for lev in range(nlevel, -1, -1):
    Im1warp = warp(Im1, -u, -v)
    Im1c, off = coarsen(Im1warp, lev)
    Im2c, off = coarsen(Im2, lev)
    maskc, off = coarsen(mask, lev)

    niter = 100
    w1 = 100
    w2 = 0
    Ix = 0.5*(deriv_x(Im1c) + deriv_x(Im2c))
    Iy = 0.5*(deriv_y(Im1c) + deriv_y(Im2c))
    It = Im2c - Im1c
    du = np.zeros(Ix.shape)
    dv = np.zeros(Ix.shape)
    for k in range(niter):
      du[0,:] = 0; du[-1,:] = 0; du[:,0] = 0; du[:,-1] = 0
      dv[0,:] = 0; dv[-1,:] = 0; dv[:,0] = 0; dv[:,-1] = 0
      du[np.where(maskc>0)] = 0
      dv[np.where(maskc>0)] = 0
      ubar2 = laplacian(du) + du
      vbar2 = laplacian(dv) + dv
      ubar1 = deriv_xx(du) + du
      vbar1 = deriv_yy(dv) + dv
      uxy = deriv_xy(du)
      vxy = deriv_xy(dv)
      du = (w1*ubar2 + w2*(ubar1+vxy))/(w1+w2) - Ix*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
      dv = (w1*vbar2 + w2*(vbar1+uxy))/(w1+w2) - Iy*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)

    u += sharpen(du*2**lev, lev, off)
    v += sharpen(dv*2**lev, lev, off)
    u[0,:] = 0; u[-1,:] = 0; u[:,0] = 0; u[:,-1] = 0
    v[0,:] = 0; v[-1,:] = 0; v[:,0] = 0; v[:,-1] = 0
  return u, v

def warp(Im, u, v):
  warp_Im = Im.copy()
  if(Im.ndim==2):
    ni, nj = Im.shape
  if(Im.ndim==3):
    nz, ni, nj = Im.shape
  for i in range(ni):
    for j in range(nj):
      if(Im.ndim==2):
        warp_Im[i, j] = interp2d(Im, (i+u[i, j], j+v[i, j]))
      if(Im.ndim==3):
        warp_Im[:, i, j] = interp2d(Im, (i+u[i, j], j+v[i, j]))
  return warp_Im

def warp_spacetime(f, u, v, q):
  fw = f.copy()
  nx, ny, nt = f.shape
  for t in range(nt):
    for x in range(nx):
      for y in range(ny):
        fw[x, y, t] = interp3d(f, (x+u[x, y, t], y+v[x, y, t], t+q[x, y, t]))
  return fw

def coarsen(Im, level):
  off = np.zeros(level)
  for k in range(level):
    ni, nj = Im.shape
    if(np.mod(ni,2)==0):
      off[k] = 0
      Im1 = 0.25*(Im[0:ni:2, :][:, 0:nj:2] + Im[1:ni:2, :][:, 0:nj:2] + Im[0:ni:2, 1:nj:2] + Im[1:ni:2, 1:nj:2])
    else:
      off[k] = 1
      Im1 = 0.25*(Im[0:ni-1:2, :][:, 0:nj-1:2] + Im[1:ni:2, :][:, 0:nj-1:2] + Im[0:ni-1:2, 1:nj:2] + Im[1:ni:2, 1:nj:2])
    Im = Im1
  return Im, off

def sharpen(Im, level, off):
  for k in range(level):
    ni, nj = Im.shape
    if(off[level-1-k]==0):
      Im1 = np.zeros((ni*2, nj))
      Im1[0:ni*2:2, :] = Im
      Im1[1:ni*2:2, :] = 0.5*(np.roll(Im, -1, axis=0) + Im)
      Im2 = np.zeros((ni*2, nj*2))
      Im2[:, 0:nj*2:2] = Im1
      Im2[:, 1:nj*2:2] = 0.5*(np.roll(Im1, -1, axis=1) + Im1)
    else:
      Im1 = np.zeros((ni*2+1, nj))
      Im1[0:ni*2:2, :] = Im
      Im1[1:ni*2+1:2, :] = 0.5*(np.roll(Im, -1, axis=0) + Im)
      Im2 = np.zeros((ni*2+1, nj*2+1))
      Im2[:, 0:nj*2:2] = Im1
      Im2[:, 1:nj*2+1:2] = 0.5*(np.roll(Im1, -1, axis=1) + Im1)
    Im = Im2
  return Im

def interp2d(x, loc):
  if(x.ndim==2):
    ni, nj = x.shape
  if(x.ndim==3):
    nz, ni, nj = x.shape
  io = loc[0]
  jo = loc[1]
  io1 = int(np.floor(io)) % ni
  jo1 = int(np.floor(jo)) % nj
  io2 = int(np.floor(io+1)) % ni
  jo2 = int(np.floor(jo+1)) % nj
  di = io - np.floor(io)
  dj = jo - np.floor(jo)
  if(x.ndim==2):
    xo = (1-di)*(1-dj)*x[io1, jo1] + di*(1-dj)*x[io2, jo1] + (1-di)*dj*x[io1, jo2] + di*dj*x[io2, jo2]
  if(x.ndim==3):
    xo = (1-di)*(1-dj)*x[:, io1, jo1] + di*(1-dj)*x[:, io2, jo1] + (1-di)*dj*x[:, io1, jo2] + di*dj*x[:, io2, jo2]
  return xo

def interp3d(x, loc):
  ni, nj, nk = x.shape
  io = loc[0]
  jo = loc[1]
  ko = loc[2]
  io1 = int(np.floor(io)) % ni
  jo1 = int(np.floor(jo)) % nj
  ko1 = int(np.floor(ko)) % nk
  io2 = int(np.floor(io+1)) % ni
  jo2 = int(np.floor(jo+1)) % nj
  ko2 = int(np.floor(ko+1)) % nk   ####caution here!!!
  di = io - np.floor(io)
  dj = jo - np.floor(jo)
  dk = ko - np.floor(ko)
  xo1 = (1-di)*(1-dj)*x[io1, jo1, ko1] + di*(1-dj)*x[io2, jo1, ko1] + (1-di)*dj*x[io1, jo2, ko1] + di*dj*x[io2, jo2, ko1]
  xo2 = (1-di)*(1-dj)*x[io1, jo1, ko2] + di*(1-dj)*x[io2, jo1, ko2] + (1-di)*dj*x[io1, jo2, ko2] + di*dj*x[io2, jo2, ko2]
  xo = (1-dk)*xo1 + dk*xo2
  return xo

def regrid(x1, ni, nj):
  ni1, nj1 = x1.shape
  ii1, jj1 = np.mgrid[0:ni1:1.0*ni1/ni, 0:nj1:1.0*nj1/nj]
  x = np.zeros((ni, nj))
  for i in range(ni):
    for j in range(nj):
      x[i, j] = interp2d(x1, np.array([ii1[i, j], jj1[i, j]]))
  return x


def runningsmooth(x):
  nx = x.size
  x2 = x.copy()
  x2[1] = (x[0]+x[1]+x[2])/3.0
  x2[2] = (x[1]+x[2]+x[3]+x[4]+x[5])/5.0
  x2[nx-2] = (x[nx-3]+x[nx-2]+x[nx-1])/3.0
  x2[nx-3] = (x[nx-5]+x[nx-4]+x[nx-3]+x[nx-2]+x[nx-1])/5.0
  for i in np.arange(3, nx-3):
    x2[i] = (x[i+3]+x[i+2]+x[i+1]+x[i]+x[i-1]+x[i-2]+x[i-3])/7.0
  return x2

def smooth2d(x, smth):
  if smth > 0:
    x_smooth = np.zeros(x.shape)
    cw = 0.0
    for i in np.arange(-smth, smth, 1):
      for j in np.arange(-smth, smth, 1):
        w = np.exp(-(i**2+j**2)/(smth/2.0)**2)
        cw += w
        x_smooth += w * np.roll(np.roll(x, j, axis=0), i, axis=1)
    x_smooth = x_smooth/cw
  else:
    x_smooth = x
  return x_smooth

def smooth_spec(wn, pwr, smth):
  n = wn.size
  pwr_smth = pwr.copy()
  for m in range(1, n):
    pwr_smth[m] = np.mean(pwr[int(max(0, np.floor(m/smth))):int(min(np.ceil(m*smth), n-1))+1])
  return pwr_smth

####random fields
def generate_fft_index(n):
  nup = int(np.ceil((n+1)/2))
  if n%2 == 0:
    wn = np.concatenate((np.arange(0, nup), np.arange(2-nup, 0)))
  else:
    wn = np.concatenate((np.arange(0, nup), np.arange(1-nup, 0)))
  return wn

def gaussian_random_field(Pk, n):
  wn = generate_fft_index(n)
  kx, ky = np.meshgrid(wn, wn)
  k2d = np.sqrt(kx**2 + ky**2)
  k2d[np.where(k2d==0.0)] = 1e-10
  noise = np.fft.fft2(np.random.normal(0, 1, (n, n)))
  amplitude = Pk(k2d)
  amplitude[np.where(k2d==1e-10)] = 0.0
  noise1 = np.real(np.fft.ifft2(noise * amplitude))
  return (noise1 - np.mean(noise1))/np.std(noise1)

####diagnostics
def rmse(x, xt):
  return np.sqrt(np.mean((x-xt)**2))

def sprd(xens):
  return np.sqrt(np.mean(np.std(xens, axis=0)**2))

def skewness(xens):
  from scipy import stats
  return np.sqrt(np.mean(stats.skew(xens, axis=0)**2))

def pattern_correlation(x, xt):
  nx, ny = x.shape
  x_mean = np.mean(x)
  xt_mean = np.mean(xt)
  xp = x - x_mean
  xtp = xt - xt_mean
  cov = np.sum(xp * xtp)
  x_norm = np.sum(xp ** 2)
  xt_norm = np.sum(xtp ** 2)
  pcorr = cov/np.sqrt(x_norm * xt_norm)
  return pcorr

def sample_correlation(x1, x2):
  nens = x1.size
  x1_mean = np.mean(x1)
  x2_mean = np.mean(x2)
  x1p = x1 - x1_mean
  x2p = x2 - x2_mean
  cov = np.sum(x1p * x2p)
  x1_norm = np.sum(x1p ** 2)
  x2_norm = np.sum(x2p ** 2)
  corr = cov/np.sqrt(x1_norm * x2_norm)
  return corr

def ttest(a, b):
  from scipy import stats
  N = a.size
  var_a = a.var(ddof=1)
  var_b = b.var(ddof=1)
  s = np.sqrt((var_a + var_b)/2)
  diff = a.mean() - b.mean()
  t = diff/(s * np.sqrt(2/N))
  df = 2*N - 2
  if diff < 0:
    p_value = stats.t.cdf(t, df=df)
  else:
    p_value = 1 - stats.t.cdf(t, df=df)
  return diff, p_value
