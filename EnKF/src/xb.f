!=======================================================================================
subroutine xb_to_surface(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,xland,lu_index,znu,znw,p_top,times,xb)
use constants
use namelist_define
use obs_define
use netcdf
use wrf_tools
use map_utils
implicit none
character(len=10), intent(in)           :: inputfile
character(len=19), intent(in)           :: times
character(len=10) :: obstype
type(proj_info), intent(in)             :: proj                   ! 1st guest map info
real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
integer, intent(in)                     :: ix, jx, kx, nv, iob
real, dimension(ix, jx ), intent(in)    :: xlong, xland, lu_index
real,  intent(in)                       :: p_top
real, dimension(kx), intent(in)         :: znu
real, dimension(kx+1), intent(in)       :: znw
real, intent(out)                       :: xb
real, dimension(ix, jx, kx+1)           :: ph, phb
real, dimension(ix, jx, kx  )           :: pt, qv, qc, qr, pb
real, dimension(ix, jx      )           :: mu, mub, t2m, th2m, q2m, u10m, v10m
real, dimension(ix+1, jx, kx)           :: u
real, dimension(ix, jx+1, kx)           :: v
real, dimension(2, 2, kx+1)             :: p
real, dimension(kx)                     :: pres, ptt, qvt, ht
real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v, dir, spd
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk, obs_ii,obs_jj
integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_u, i_v
integer                                 :: i_t2, i_th2, i_q2, i_u10, i_v10
real                                    :: dx, dxm, dy, dym, dz, dzm
real                                    :: psfc, tsfc, psfcm, u10, v10, t2, q2, th2
real, dimension(ix, jx      )           :: rough

xb=-888888.

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

call roughness_from_landuse ( 'USGS', times, ix, jx, lu_index, rough ) 

!- calculate q, pressure profile on (obs_ii, obs_jj)
i_qv = 0 
i_qc = 0 
i_qr = 0 
i_mu = 0
i_mub = 0
i_ph = 0
i_phb = 0
i_pt = 0
i_u = 0
i_v = 0
i_t2 = 0
i_th2 = 0
i_q2 = 0
i_u10 = 0
i_v10 = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
   if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
   if ( enkfvar(m) == 'MU        ' ) i_mu=m
   if ( enkfvar(m) == 'MUB       ' ) i_mub=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_pt=m
   if ( enkfvar(m) == 'U         ' ) i_u=m
   if ( enkfvar(m) == 'V         ' ) i_v=m
   if ( enkfvar(m) == 'T2        ' ) i_t2=m
   if ( enkfvar(m) == 'TH2       ' ) i_th2=m
   if ( enkfvar(m) == 'Q2        ' ) i_q2=m
   if ( enkfvar(m) == 'U10       ' ) i_u10=m
   if ( enkfvar(m) == 'V10       ' ) i_v10=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_qc>0) qc (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qc )
if(i_qr>0) qr (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qr )
if(i_mu>0)  mu (i1:i1+1, j1:j1+1)   = xa( 1:2, 1:2, 1, i_mu )
if(i_mub>0) mub (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_mub )
if(i_t2>0)  t2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_t2 )
if(i_th2>0) th2m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_th2 )
if(i_q2>0)  q2m (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_q2 )
if(i_u10>0) u10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_u10 )
if(i_v10>0) v10m (i1:i1+1, j1:j1+1) = xa( 1:2, 1:2, 1, i_v10 )
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu  )
if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub )
if ( i_t2== 0  ) call get_variable2d(inputfile, 'T2        ', ix, jx, 1,    t2m )
if ( i_th2== 0 ) call get_variable2d(inputfile, 'TH2       ', ix, jx, 1,    th2m)
if ( i_q2== 0  ) call get_variable2d(inputfile, 'Q2        ', ix, jx, 1,    q2m )
if ( i_u10== 0 ) call get_variable2d(inputfile, 'U10       ', ix, jx, 1,    u10m)
if ( i_v10== 0 ) call get_variable2d(inputfile, 'V10       ', ix, jx, 1,    v10m)

qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)

!...... get total qv at obs' position from horizontal interpolation
qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
!...... get mu,mub at obs' position from horizontal interpolation
mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
!...... calculate pressure profile from qv
call cal_press_from_q( kx, znu, znw, qvt, mu1, mub1, p_top, pres )

!- calculate t (not theta) and height profile
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_phb )
if(i_pt>0) pt (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_pt )
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )

!...... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
p(1:2,1:2,1:3) = (ph(i1:i1+1, j1:j1+1, 1:3) + phb(i1:i1+1, j1:j1+1, 1:3))/g
p(1,1,1:3) = dym*(dx*p(2,1,1:3) + dxm*p(1,1,1:3)) + dy*(dx*p(2,2,1:3) + dxm*p(1,2,1:3))
ht(1:2) = 0.5*(p(1,1,1:2)+p(1,1,2:3))
!...... get ptt(theta)  at obs' position from horizontal interpolation
ptt(1) = dym*(dx*pt(i1+1,j1,  1) + dxm*pt(i1,j1,  1)) + dy*(dx*pt(i1+1,j1+1,1) + dxm*pt(i1,j1+1,1))
ptt(1) = theta_to_temp(ptt(1)+to, pres(1))
 
!-------------------------------------------------
!- calculate surface p and ground t 
   psfc = pres(1) + (1. - znu(1))*mub1
   tsfc = ptt(1) + 0.0065*(ht(1)-p(1,1,1))
!
!-------------------------------------------------
!...... get model u and v
if(i_u>0) u(i1:i1+2,j1:j1+1,1)=xa(1:3,1:2,1,i_u)
if(i_v>0) v(i1:i1+1,j1:j1+2,1)=xa(1:2,1:3,1,i_v)
if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )

!...... horizontal interp for U
if ( obs_ii-i1 == 0.500 ) then
   grid_u = u(i1+1,j1,1)*(j1+1-obs_jj) + u(i1+1,j1+1,1)*(obs_jj-j1)
else if ( obs_ii-i1 > 0.500 ) then
   grid_u = (j1+1-obs_jj)*( u(i1+1,j1  ,1)*(i1+1.5-obs_ii)+u(i1+2,j1  ,1)*(obs_ii-i1-0.5) ) + &
             (obs_jj-j1)  *( u(i1+1,j1+1,1)*(i1+1.5-obs_ii)+u(i1+2,j1+1,1)*(obs_ii-i1-0.5) )
else if ( obs_ii-i1 < 0.500 ) then
   grid_u = (j1+1-obs_jj)*( u(i1,j1  ,1)*(i1+0.5-obs_ii)+u(i1+1,j1  ,1)*(obs_ii-i1+0.5) ) + &
             (obs_jj-j1)  *( u(i1,j1+1,1)*(i1+0.5-obs_ii)+u(i1+1,j1+1,1)*(obs_ii-i1+0.5) )
endif

!...... horizontal interp for V
if ( obs_jj-j1 == 0.500 ) then
   grid_v = v(i1,j1+1,1)*(i1+1-obs_ii) + v(i1+1,j1+1,1)*(obs_ii-i1)
else if ( obs_jj-j1 > 0.500 ) then
   grid_v = (i1+1-obs_ii)*( v(i1  ,j1+1,1)*(j1+1.5-obs_jj) + v(i1  ,j1+2,1)*(obs_jj-j1-0.5) ) + &
             (obs_ii-i1)  *( v(i1+1,j1+1,1)*(j1+1.5-obs_jj) + v(i1+1,j1+2,1)*(obs_jj-j1-0.5) )
else if ( obs_jj-j1 < 0.500 ) then
   grid_v = (i1+1-obs_ii)*( v(i1  ,j1,1)*(j1+0.5-obs_jj) + v(i1  ,j1+1,1)*(obs_jj-j1+0.5) ) + &
             (obs_ii-i1)  *( v(i1+1,j1,1)*(j1+0.5-obs_jj) + v(i1+1,j1+1,1)*(obs_jj-j1+0.5) )
endif

!-------------------------------------------------
!- calculate 10m wind, 2m t and q
!call sfc_wtq( psfc, tsfc, pres(1), ptt(1), qvt(1), grid_u, grid_v,          &
!              pres(2), ptt(2), qvt(2), ht(1), rough(i1,j1), xland(i1,j1),   &
!              u10, v10, t2, q2 )  
t2  = dym*(dx*t2m (i1+1,j1) + dxm*t2m (i1,j1)) + dy*(dx*t2m (i1+1,j1+1) + dxm*t2m (i1, j1+1))
th2 = dym*(dx*th2m(i1+1,j1) + dxm*th2m(i1,j1)) + dy*(dx*th2m(i1+1,j1+1) + dxm*th2m(i1, j1+1))
q2  = dym*(dx*q2m (i1+1,j1) + dxm*q2m (i1,j1)) + dy*(dx*q2m (i1+1,j1+1) + dxm*q2m (i1, j1+1))
u10 = dym*(dx*v10m(i1+1,j1) + dxm*u10m(i1,j1)) + dy*(dx*u10m(i1+1,j1+1) + dxm*u10m(i1, j1+1))
v10 = dym*(dx*u10m(i1+1,j1) + dxm*v10m(i1,j1)) + dy*(dx*v10m(i1+1,j1+1) + dxm*v10m(i1, j1+1))


!-------------------------------------------------
!- Correct surface pressure
call da_sfc_pre ( psfcm, psfc, t2, q2, p(1,1,1), obs%sta(iob,1), obs%sta(iob,3), obs%sta(iob,4)/1000.)

!if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1  ,1), u(i1+1,j1  ,1), u(i1+2,j1  ,1)  
!if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1+1,1), u(i1+1,j1+1,1), u(i1+2,j1+1,1)  
!if( print_detail > 100 )write(*,'(3x,a,5f)')'grid_u, grid_v, u10, v10 =',grid_u, grid_v, u10, v10
!if( print_detail > 100 )write(*,'(3x,a,4f)')'station elevation, p, t, q =', obs%sta(iob,1:4)
!if( print_detail > 100 )write(*,'(3x,a,5f)')'t2, q2, model_terrain =', t2, q2, p(1,1,1)
!if( print_detail > 100 )write(*,'(3x,a,5f)')'psfc and corrected psfc =', psfc, psfcm
!-------------------------------------------------
!- get xb
!...... surface pressure
if ( obstype(10:10) == 'P' ) then
     xb = psfcm
     xb = psfc
else if ( obstype(10:10) == 'U' ) then
     xb = u10
else if ( obstype(10:10) == 'V' ) then
     xb = v10
else if ( obstype(10:10) == 'T' ) then
     xb = t2
else if ( obstype(10:10) == 'Q' ) then
     xb = q2*1000.
else if ( obstype(10:10) == 'D' ) then
     xb = mixrat_to_tdew(q2, psfcm)
else if ( obstype(10:10) == 'R' ) then
     xb = rel_humidity(q2, t2, psfcm)
else if ( obstype(10:10) == 'S' ) then   !wind speed
     long = ( xlong(i1  ,j1)*dym + xlong(i1  ,j1+1)*dy ) * dxm +          & 
            ( xlong(i1+1,j1)*dym + xlong(i1+1,j1+1)*dy ) * dx
     call gridwind_to_truewind(long, proj, grid_u, grid_v, true_u, true_v ) 
     call uv_to_dirspd(true_u,true_v,dir,spd)
     xb = spd 
     !xb = sqrt(u10**2.+v10**2.)
endif

!if( print_detail > 100 )write(*,'(3x,a,3f)')'xb_to_surface '//obstype//' obs xb, obs_ii, obs_jj :', xb, obs_ii, obs_jj

end subroutine xb_to_surface

!=======================================================================================
subroutine calc_sounding_position_k (inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,obs_kk)
use constants
use namelist_define
use obs_define
use netcdf
use wrf_tools
use map_utils
use mpi_module

implicit none
character(len=10), intent(in)           :: inputfile
type(proj_info), intent(in)             :: proj                   ! 1st guest map info
real, dimension(3,3,kx+1,nv), intent(in)  :: xa                     ! 1st guest
integer, intent(in)                     :: ix, jx, kx, nv, iob 
real, dimension(ix, jx ), intent(in)    :: xlong
real, dimension(kx+1), intent(in)       :: znw
real, dimension(kx), intent(in)         :: znu
real,  intent(in)                       :: p_top

real :: obs_ii, obs_jj
real, intent(out)   :: obs_kk
character(len=10)   :: obstype

real, dimension(ix, jx, kx+1)           :: ph, phb
real, dimension(ix, jx, kx  )           :: pt, pb, qv, qc, qr
real, dimension(ix, jx      )           :: mu, mub
real, dimension(ix+1, jx, kx)           :: u
real, dimension(ix, jx+1, kx)           :: v
real, dimension(kx+1)                   :: work, worku, workv, ph1
real, dimension(kx)                     :: pres, ptt, qvt
real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v
real                                    :: dx, dxm, dy, dym, dz, dzm
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk
integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj ) 
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

!- calculate pressure profile on (obs_ii, obs_jj)
 i_ph = 0
 i_phb = 0
 i_mu = 0
 i_mub = 0
 i_pt = 0
 i_qv = 0
 do m = 1, nv
    if ( enkfvar(m) == 'PH        ' ) i_ph=m
    if ( enkfvar(m) == 'PHB       ' ) i_phb=m
    if ( enkfvar(m) == 'T         ' ) i_pt=m
    if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
    if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
    if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
    if ( enkfvar(m) == 'MU        ' ) i_mu=m
    if ( enkfvar(m) == 'MUB       ' ) i_mub=m
 enddo
 if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
 if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
 if(i_pt>0) pt(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pt)
 if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
 if(i_qc>0) qc(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qc)
 if(i_qr>0) qr(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qr)
 if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
 if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)
 if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
 if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
 if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx, 1, pt )
 if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
 if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx, 1, qc )
 if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx, 1, qr )
 if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1, mu )
 if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1, mub)

!...... get mut mu, mub at obs' position from horizontal interpolation
 mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
 mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))

 qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
 qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
 ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))

 ph(i1:i1+1,j1:j1+1,1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
 ph1(1:kx+1) = dym*(dx*ph(i1+1,j1,1:kx+1) + dxm*ph(i1,j1,1:kx+1)) + dy*(dx*ph(i1+1,j1+1,1:kx+1) + dxm*ph(i1,j1+1,1:kx+1))

 if ( obstype(1:1) == 'H' ) then
    call to_zk(obs%position(iob,4), ph1(1:kx)/g, obs_kk,kx) 
 else if ( obstype(1:1) == 'P' ) then
    call eta_to_pres(znw(1:kx+1), mu1+mub1, qvt(1:kx), ph1(1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))
    call to_zk(obs%position(iob,4), pres(1:kx), obs_kk,kx)
 endif

 !if(obs_kk>kx .or. obs_kk<1) print *, 'obs z location error: k=',obs_kk
end subroutine calc_sounding_position_k 

subroutine xb_to_sounding (inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,xb)
use constants
use namelist_define
use obs_define
use netcdf
use wrf_tools
use map_utils
use mpi_module

implicit none
character(len=10), intent(in)           :: inputfile
type(proj_info), intent(in)             :: proj                   ! 1st guest map info
real, dimension(3,3,kx+1,nv), intent(in)  :: xa                     ! 1st guest
integer, intent(in)                     :: ix, jx, kx, nv, iob 
real, dimension(ix, jx ), intent(in)    :: xlong
real, dimension(kx+1), intent(in)       :: znw
real, dimension(kx), intent(in)         :: znu
real,  intent(in)                       :: p_top

real :: obs_ii, obs_jj, obs_kk
real, intent(out)   :: xb
character(len=10)   :: obstype

real, dimension(ix, jx, kx+1)           :: ph, phb
real, dimension(ix, jx, kx  )           :: pt, pb, qv, qc, qr
real, dimension(ix, jx      )           :: mu, mub
real, dimension(ix+1, jx, kx)           :: u
real, dimension(ix, jx+1, kx)           :: v

real, dimension(kx+1)                   :: work, worku, workv, ph1
real, dimension(kx)                     :: pres, ptt, qvt
real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v
real                                    :: dx, dxm, dy, dym, dz, dzm
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk
integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_u, i_v

xb=-888888.

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
obs_kk=obs%position(iob,3)
i1 = int( obs_ii )
j1 = int( obs_jj ) 
k1  = int( obs_kk )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj
dz  = obs_kk-real(k1)
dzm = real(k1+1)-obs_kk

 i_ph = 0
 i_phb = 0
 i_mu = 0
 i_mub = 0
 i_pt = 0
 i_qv = 0
 do m = 1, nv
    if ( enkfvar(m) == 'PH        ' ) i_ph=m
    if ( enkfvar(m) == 'PHB       ' ) i_phb=m
    if ( enkfvar(m) == 'T         ' ) i_pt=m
    if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
    if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
    if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
    if ( enkfvar(m) == 'MU        ' ) i_mu=m
    if ( enkfvar(m) == 'MUB       ' ) i_mub=m
 enddo
 if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
 if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
 if(i_pt>0) pt(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pt)
 if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
 if(i_qc>0) qc(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qc)
 if(i_qr>0) qr(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qr)
 if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
 if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)
 if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
 if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
 if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx, 1, pt )
 if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
 if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx, 1, qc )
 if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx, 1, qr )
 if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1, mu )
 if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1, mub)

!...... get mut mu, mub at obs' position from horizontal interpolation
 mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
 mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))

 qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
 qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
 ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))
 ph(i1:i1+1,j1:j1+1,1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
 ph1(1:kx+1) = dym*(dx*ph(i1+1,j1,1:kx+1) + dxm*ph(i1,j1,1:kx+1)) + dy*(dx*ph(i1+1,j1+1,1:kx+1) + dxm*ph(i1,j1+1,1:kx+1))
 call eta_to_pres(znw(1:kx+1), mu1+mub1, qvt(1:kx), ph1(1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))

!.. U, V
if ( obstype(10:10) == 'U' .or. obstype(10:10) == 'V' .or. obstype(10:10) == 'S' ) then
   i_u = 0
   i_v = 0
   do m = 1, nv
      if ( enkfvar(m) == 'U         ' ) i_u=m
      if ( enkfvar(m) == 'V         ' ) i_v=m
   enddo
   u(i1:i1+2,j1:j1+1,k1:k1+1)=xa(1:3,1:2,k1:k1+1,i_u)
   v(i1:i1+1,j1:j1+2,k1:k1+1)=xa(1:2,1:3,k1:k1+1,i_v)

   worku = -88888.
   workv = -88888.
   if ( obs_ii-i1 == 0.500 ) then
      worku(k1:k1+1) = dym*u(i1+1,j1,k1:k1+1) + dy*u(i1+1,j1+1,k1:k1+1)
   else if ( obs_ii-i1 > 0.500 ) then
      worku(k1:k1+1) = dym*( u(i1+1,j1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1,k1:k1+1)*(dx-0.5) ) + &
                       dy *( u(i1+1,j1+1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1+1,k1:k1+1)*(dx-0.5) )
   else if ( obs_ii-i1 < 0.500 ) then
      worku(k1:k1+1) = dym*( u(i1,j1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1,k1:k1+1)*(dx+0.5) ) + &
                       dy *( u(i1,j1+1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1+1,k1:k1+1)*(dx+0.5) )
   endif
   if ( obs_jj-j1 == 0.500 ) then
      workv(k1:k1+1) = v(i1,j1+1,k1:k1+1)*dxm + v(i1+1,j1+1,k1:k1+1)*dx
   else if ( obs_jj-j1 > 0.500 ) then
      workv(k1:k1+1) = dxm*( v(i1,j1+1,k1:k1+1)*(dym+0.5) + v(i1,j1+2,k1:k1+1)*(dy-0.5) ) + &
                       dx *( v(i1+1,j1+1,k1:k1+1)*(dym+0.5) + v(i1+1,j1+2,k1:k1+1)*(dy-0.5) )
   else if ( obs_jj-j1 < 0.500 ) then
      workv(k1:k1+1) = dxm*( v(i1,j1,k1:k1+1)*(dym-0.5) + v(i1,j1+1,k1:k1+1)*(dy+0.5) ) + &
                       dx *( v(i1+1,j1,k1:k1+1)*(dym-0.5) + v(i1+1,j1+1,k1:k1+1)*(dy+0.5) )
   endif

!       if( print_detail > 100 ) write(*,*)'i1,j1,k1, v=',i1,j1,k1,v(i1,j1,k1),v(i1,j1+1,k1),v(i1,j1+2,k1),v(i1+1,j1,k1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'U profile =', worku(k1:k1+1)
!       if( print_detail > 100 ) write(*,'(a,100f10.2)')'V profile =', workv(k1:k1+1)

   grid_u = dzm*worku(k1)+dz*worku(k1+1)
   grid_v = dzm*workv(k1)+dz*workv(k1+1)
!       if( print_detail > 100 ) write(*,'(a,f10.2)')'grid U =', grid_u, 'grid V =', grid_v

!------011108 changed obs. wind to gridwind
   long = ( xlong(i1  ,j1)*dym + xlong(i1  ,j1+1)*dy ) * dxm +          & 
          ( xlong(i1+1,j1)*dym + xlong(i1+1,j1+1)*dy ) * dx
   call gridwind_to_truewind( long, proj, grid_u, grid_v, true_u, true_v )
 
   if ( obstype(10:10) == 'U' ) then
       xb = grid_u
   else if ( obstype(10:10) == 'V' ) then
       xb = grid_v
   else if ( obstype(10:10) == 'S' ) then   !wind speed
       xb = sqrt(grid_u**2.+grid_v**2.)
   endif
!       if( print_detail > 100 ) then
!           write(*,'(12f10.3)')obs_kk, grid_u, grid_v, dz, dzm, long, worku(k1:k1+1), workv(k1:k1+1), true_u, true_v
!       endif

!.. T
else if ( obstype(10:10) == 'T' ) then
   do k = k1, k1+1
      work(k) = theta_to_temp(ptt(k)+to, pres(k))
   enddo
   if ( obs_kk .le. 1. ) then
      xb = work(1)
   else
      xb = dzm*work(k1)+dz*work(k1+1)
   endif

!.. TD(if used TD, not RH)
else if ( obstype(10:10) == 'D' ) then
   do k = k1, k1+1
      work(k) = mixrat_to_tdew(qvt(k), pres(k))
   enddo
   if ( obs_kk .le. 1. ) then
      xb = work(1)
   else
      xb = dzm*work(k1) + dz*work(k1+1)
   endif

!.. RH
else if ( obstype(10:10) == 'R' )then
   do k = k1, k1+1
      work(k) = theta_to_temp(ptt(k)+to, pres(k))
   enddo
   do k = k1, k1+1
      work(k) = rel_humidity(qvt(k), work(k), pres(k))
   enddo
   if ( obs_kk .le. 1. ) then
      xb = work(1)
   else
      xb = dzm*work(k1) + dz*work(k1+1)
   endif
   if ( xb > 100. )xb = 100.

!.. Q
else if ( obstype(10:10) == 'Q' )then
   if ( obs_kk .le. 1. ) then
      xb = qvt(k1)*1000.
   else
      xb = (dzm*qvt(k1)+dz*qvt(k1+1))*1000.
   endif
   if ( xb .le. 0.0 ) xb = -99999.

!.. HG
else if ( obstype(10:10) == 'H' )then
   call destag_zstag(znu, znw, kx, ph1(1:kx+1), work(1:kx))
   if ( obs_kk .le. 1. ) then
      xb = work(1)/g
   else
      xb = (dzm*work(k1) + dz*work(k1+1))/g
   endif 

!.. GPSRO refractivity
else if ( obstype(10:10) == 'N' ) then
   do k = k1, k1+1
      work(k) = theta_to_temp(ptt(k)+to, pres(k))
   end do
   do k = k1, k1+1
      work(k) = gpsref(pres(k), work(k), qvt(k))
   end do
   if ( obs_kk .le. 1. ) then
      xb = work(k1)
   else
      xb = dzm*work(k1)+dz*work(k1+1)
   endif
   !if ( xb .lt. 0.0 ) xb = -99999.

endif

   !if( print_detail > 1 )write(*,'(a,3f8.2, f10.1)')'xb_to_sounding '//obstype//' obs position :', obs_ii, obs_jj, obs%position(iob,3:4)
    
end subroutine xb_to_sounding


!=======================================================================================
subroutine calc_rv_position_k(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znw,obs_kk)
  use constants 
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils 
  use netcdf
  use wrf_tools
  use radar
  implicit none
  character (len=10), intent(in) :: inputfile
  type(proj_info), intent(in) :: proj
  integer, intent(in)         :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in)  :: xa
  real, dimension(ix, jx ), intent(in) :: xlong
  real, dimension(kx+1), intent(in)    :: znw
  real, intent(out)                    :: obs_kk 
  real  :: obs_ii, obs_jj, obs_hh, radar_ii, radar_jj, radar_elv

  real, dimension(ix+1, jx, kx)        :: u
  real, dimension(ix, jx+1, kx)        :: v
  real, dimension(ix, jx, kx+1)        :: w, ph, phb
  real, dimension(ix, jx, kx  )        :: t, qv, p, pb
  real, dimension(ix, jx      )        :: mu, mub, terrain
    
  integer :: m, ii, jj, kk, i, j, k, i1, j1
  integer :: i_u, i_v, i_w, i_t, i_qv, i_mu, i_ph, i_phb, i_mub, i_pb
  obs_ii=obs%position(iob,1)
  obs_jj=obs%position(iob,2)
  obs_hh=obs%position(iob,4)

  radar_ii=obs%sta(iob,1)
  radar_jj=obs%sta(iob,2)
  radar_elv=obs%sta(iob,4)

  i1 = int( obs_ii )
  j1 = int( obs_jj )

  i_u  = 0 ; i_v   = 0 ; i_w  = 0 ; i_t = 0 ; i_qv = 0 ; i_mu = 0; i_mub = 0;
  i_ph = 0 ; i_phb = 0 ; i_pb = 0 ;
     do m = 1, nv
        if ( enkfvar(m) == 'T         ' ) i_t=m
        if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
        if ( enkfvar(m) == 'PH        ' ) i_ph=m
        if ( enkfvar(m) == 'PHB       ' ) i_phb=m
        if ( enkfvar(m) == 'PB        ' ) i_pb=m
        if ( enkfvar(m) == 'MU        ' ) i_mu=m
        if ( enkfvar(m) == 'MUB       ' ) i_mub=m
     end do
     if(i_t>0) t(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_t)
     if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
     if(i_pb>0) pb(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pb)
     if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
     if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
     if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
     if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)
     if( i_t <1  ) call get_variable3d(inputfile, 'T         ', ix, jx, kx, 1, t )
     if( i_qv <1 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
     if( i_pb <1 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx, 1, pb )
     if( i_ph <1 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
     if( i_phb <1) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb )
     if( i_mu <1 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1, mu )
     if( i_mub <1) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1, mub )
     do j = j1, j1+1
     do i = i1, i1+1
        call cal_ph( kx, znw, t(i,j,1:kx), qv(i,j,1:kx), pb(i,j,1:kx), mu(i,j), mub(i,j), phb(i,j,1:kx+1), ph(i,j,1:kx+1) )
     enddo
     enddo

!     if (print_detail > 1000) write(*,'(a,70f8.0)') 'phb+ph = ',(phb(i1,j1,1:kx+1)+ph(i1,j1,1:kx+1))/g
!     if (print_detail > 1000) write(*,'(a,f8.0)') 'obs_hh =', obs_hh
     call calc_radar_data_point_kk ( ix, jx, kx, ph, phb, obs_ii, obs_jj, obs_hh, obs_kk )
end subroutine calc_rv_position_k

subroutine xb_to_rv(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znw,xb)
  use constants 
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils 
  use netcdf
  use wrf_tools
  use radar
  implicit none
  character (len=10), intent(in) :: inputfile
  type(proj_info), intent(in) :: proj
  integer, intent(in)         :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in)  :: xa
  real, dimension(ix, jx ), intent(in) :: xlong
  real, dimension(kx+1), intent(in)    :: znw
  real, intent(out)                    :: xb
  real  :: obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv

  real, dimension(ix+1, jx, kx)        :: u
  real, dimension(ix, jx+1, kx)        :: v
  real, dimension(ix, jx, kx+1)        :: w, ph, phb
  real, dimension(ix, jx, kx  )        :: t, qv, p, pb
  real, dimension(ix, jx      )        :: mu, mub, terrain
    
  integer :: m, ii, jj, kk, i, j, k, i1, j1
  integer :: i_u, i_v, i_w, i_t, i_qv, i_mu, i_ph, i_phb, i_mub, i_pb

  obs_ii=obs%position(iob,1)
  obs_jj=obs%position(iob,2)
  obs_kk=obs%position(iob,3)
  obs_hh=obs%position(iob,4)

  radar_ii=obs%sta(iob,1)
  radar_jj=obs%sta(iob,2)
  radar_elv=obs%sta(iob,4)

  i1 = int( obs_ii )
  j1 = int( obs_jj )
    
  i_u  = 0 ; i_v   = 0 ; i_w  = 0 ; i_t = 0 ; i_qv = 0 ; i_mu = 0; i_mub = 0;
  i_ph = 0 ; i_phb = 0 ; i_pb = 0 ;

  do m=1,nv
     if ( enkfvar(m) == 'U         ' ) i_u=m
     if ( enkfvar(m) == 'V         ' ) i_v=m
     if ( enkfvar(m) == 'W         ' ) i_w=m
  enddo
  if(i_u>0) u(i1:i1+2,j1:j1+1,1:kx)=xa(1:3,1:2,1:kx,i_u)
  if(i_v>0) v(i1:i1+1,j1:j1+2,1:kx)=xa(1:2,1:3,1:kx,i_v)
  if(i_w>0) w(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_w)
  if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
  if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )
  if ( i_w == 0 ) call get_variable3d(inputfile, 'W         ', ix, jx, kx+1, 1, w )
  

  call calculate_rv ( ix, jx, kx, u, v, w, xlong, proj, obs_ii, obs_jj, obs_kk, obs_hh, radar_ii, radar_jj, radar_elv, xb )
!  write(*,*)'i,j,k =', obs_ii, obs_jj, obs_kk
!  write(*,*)'u,v,w =', u(int(obs_ii),int(obs_jj),int(obs_kk)), v(int(obs_ii),int(obs_jj),int(obs_kk)), w(int(obs_ii),int(obs_jj),int(obs_kk))

end subroutine xb_to_rv

!=======================================================================================
subroutine xb_to_pw(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb)
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools
  implicit none
  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb
  real, dimension(2, 2)                   :: pw,temp,vtemp,zdiff
  real, dimension(ix+1, jx, kx  )         :: u
  real, dimension(ix, jx+1, kx  )         :: v
  real, dimension(ix, jx, kx  )           :: t,p,pb,z,qv
  real, dimension(ix, jx, kx+1)           :: ph,phb
  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  integer                                 :: i_t, i_u, i_v, i_ph, i_phb, i_qv, i_p, i_pb
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

i_qv = 0 
i_ph = 0
i_phb = 0
i_t = 0
i_p = 0
i_pb = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'P         ' ) i_p=m
   if ( enkfvar(m) == 'PB        ' ) i_pb=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_t=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_phb )
if(i_p>0) p (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_p )
if(i_pb>0) pb (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_pb )
if(i_t>0) t (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_t )
     
if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)
do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:',i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
pw = 0.
do k=1,kx-1
  temp=(t(i1:i1+1,j1:j1+1,k)+300.)*(p(i1:i1+1,j1:j1+1,k)/Pr)**(rd/cp)
  vtemp=(1+0.61*qv(i1:i1+1,j1:j1+1,k))*temp
  zdiff=z(i1:i1+1,j1:j1+1,k+1)-z(i1:i1+1,j1:j1+1,k)
  pw=pw+p(i1:i1+1,j1:j1+1,k)/(rd*vtemp)*zdiff*qv(i1:i1+1,j1:j1+1,k)
enddo
xb = dym*(dx*pw(2,1) + dxm*pw(1,1)) + dy*(dx*pw(2,2) + dxm*pw(1,2))
xb = xb/10

end subroutine xb_to_pw

!=======================================================================================
subroutine xb_to_slp(inputfile,xa,ix,jx,kx,nv,iob,znu,znw,xb )

!---------------------
! slp_center subroutine calculates sea level pressure and find the hurricane center
!---------------------
  use constants
  use namelist_define
  use obs_define
  use netcdf
  use wrf_tools

  implicit none

  character(len=10), intent(in)           :: inputfile
  character(len=10) :: obstype
  integer, intent(in)                     :: ix, jx, kx, nv, iob
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, intent(out)                       :: xb

  real, dimension(2, 2)                   :: slp
  real, dimension(ix+1, jx, kx  )         :: u
  real, dimension(ix, jx+1, kx  )         :: v
  real, dimension(ix, jx, kx  )           :: t
  real, dimension(ix, jx, kx+1)           :: ph
  real, dimension(ix, jx, kx+1)           :: phb
  real, dimension(ix, jx, kx  )           :: p
  real, dimension(ix, jx, kx  )           :: z
  real, dimension(ix, jx, kx  )           :: pb
  real, dimension(ix, jx, kx  )           :: qv, qc, qr
  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
  integer                                 :: i_t, i_u, i_v, i_ph, i_phb, i_qv, i_qc, i_qr, i_p, i_pb  ! variables flag
  real                                    :: obs_ii, obs_jj, dx,dxm,dy,dym

obstype=obs%type(iob)
obs_ii=obs%position(iob,1)
obs_jj=obs%position(iob,2)
i1 = int( obs_ii )
j1 = int( obs_jj )
dx  = obs_ii-real(i1)
dxm = real(i1+1)-obs_ii
dy  = obs_jj-real(j1)
dym = real(j1+1)-obs_jj

i_qv = 0 
i_qc = 0 
i_qr = 0 
i_ph = 0
i_phb = 0
i_t = 0
i_p = 0
i_pb = 0
do m = 1, nv
   if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
   if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
   if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
   if ( enkfvar(m) == 'P         ' ) i_p=m
   if ( enkfvar(m) == 'PB        ' ) i_pb=m
   if ( enkfvar(m) == 'PH        ' ) i_ph=m
   if ( enkfvar(m) == 'PHB       ' ) i_phb=m
   if ( enkfvar(m) == 'T         ' ) i_t=m
enddo
if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
if(i_qc>0) qc (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qc )
if(i_qr>0) qr (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qr )
if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_ph )
if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:kx+1) = xa( 1:2,1:2,1:kx+1,i_phb )
if(i_p>0) p (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_p )
if(i_pb>0) pb (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_pb )
if(i_t>0) t (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2,1:2,1:kx,i_t )
     
if ( i_t  == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
if ( i_p  == 0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
ph(i1:i1+1, j1:j1+1, 1:kx+1) = (ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1))/g
p(i1:i1+1, j1:j1+1, 1:kx)  = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)

do j=j1,j1+1
do i=i1,i1+1
  z(i, j, 1:kx) = (ph(i, j, 1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i, j, 2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
enddo
enddo
do j=j1,j1+1
do i=i1,i1+1
do k=1,kx
  if( p(i,j,k)/100. >= 1200. ) write(*,'(a,3i4,2f8.0)')'P error at:', i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
enddo
enddo
enddo
   call compute_seaprs(2, 2, kx, z(i1:i1+1,j1:j1+1,:), t(i1:i1+1,j1:j1+1,:), p(i1:i1+1,j1:j1+1,:), qv(i1:i1+1,j1:j1+1,:), slp, 0)

   xb = dym*(dx*slp(2,1) + dxm*slp(1,1)) + dy*(dx*slp(2,2) + dxm*slp(1,2))
   xb = xb*100.

end subroutine xb_to_slp



!=======================================================================================
subroutine xb_to_radiance(inputfile,xa,ix,jx,kx,nv,iob,xlong,xlat,znw,hgt,tsk,landmask,xb,cloud_flag)
  USE constants
  USE netcdf
  USE mpi_module
  USE CRTM_Module
  use namelist_define
  use obs_define
  use wrf_tools

  implicit none
  integer, intent(in)                      :: ix, jx, kx, nv, iob, cloud_flag
  character(len=10), intent(in)            :: inputfile
  real, dimension(3,3,kx+1,nv), intent(in) :: xa                     ! 1st guest model states
  real, intent(out)                        :: xb
  real, dimension(ix, jx ), intent(in)     :: xlong, xlat, landmask, hgt, tsk
  real, dimension(kx+1),intent(in)         :: znw
  integer                                  :: i1,j1
  integer :: i_p,i_pb,i_ph,i_phb,i_pt,i_qv,i_qr,i_qc,i_qi,i_qs,i_qg,i_psfc,i_mu,i_mub
  real                                     :: obs_ii,obs_jj, dx,dxm,dy,dym,mu1,mub1

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
  REAL, PARAMETER :: Re=6378000.0
  !====================
  !setup for GOES-ABI
   REAL, PARAMETER :: sat_h=35780000.0
   REAL, PARAMETER :: sat_lon=57.0/180.0*3.14159
   INTEGER, parameter :: n_ch=2
  !====================
  !INTEGER, intent(in) :: ix = ix  !total number of the x-grid
  !INTEGER, parameter, intent(in) :: jx = jx  !total number of the y-grid
  !INTEGER, parameter, intent(in) :: kx = kx  !level range
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1 
  !INTEGER, PARAMETER :: N_LAYERS    = kx
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
  !INTEGER, PARAMETER :: N_CLOUDS    = kx*5
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: FILE_NAME
  CHARACTER(256) :: obstype
  CHARACTER(12)  :: sat_id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels 
  INTEGER :: l, m, irec, yend, ystart, nyi
  integer :: ncid,ncrcode
  character(LEN=16) :: var_name
  character(LEN=3)  :: file_ens
  integer :: x,y,tt,v,z,n,reci,ens,n_ec,num_radgrid
  INTEGER :: ncl,icl,k1,k2
  real, dimension(ix,jx) :: psfc,mu,mub
  real, dimension(ix,jx,kx) :: p,pb,pt,qv,qc,qr,qi,qs,qg
  real, dimension(ix,jx,kx+1) :: ph,phb
  real, dimension(kx) :: delz, pres,pt1,qv1,qc1,qr1,qi1,qs1,qg1,qvt,ht,tk
  real, dimension(kx+1) :: ph1
  real :: lat,lon,hgt1,tsk1,landmask1,psfc1

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)

  ! ============================================================================
  CALL CRTM_Version( Version )

  ! Get sensor id from user
  Sensor_Id = trim(adjustl(obs%sat(iob)))
  
  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. This initializes the CRTM for the sensors
  !     predefined in the example SENSOR_ID parameter.
  !     NOTE: The coefficient data file path is hard-
  !           wired for this example.
  ! --------------------------------------------------
  Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                            ChannelInfo  , &  ! Output
                            IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                            IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                            File_Path='coefficients/')
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/obs%ch(iob)/) )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 3a. Allocate the ARRAYS
  ! -----------------------
  ! Note that only those structure arrays with a channel
  ! dimension are allocated here because we've parameterized
  ! the number of profiles in the N_PROFILES parameter.
  !
  ! Users can make the 
  ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
  ALLOCATE( RTSolution( ChannelInfo(1)%n_Channels, N_PROFILES ), STAT=Allocate_Status )

  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The input FORWARD structure
  CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, kx*5, N_AEROSOLS)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

  !Observation location x,y, find z location
  obs_ii=obs%position(iob,1)
  obs_jj=obs%position(iob,2)
  i1=nint(obs_ii)
  j1=nint(obs_jj)
  dx  = obs_ii-real(i1)
  dxm = real(i1+1)-obs_ii
  dy  = obs_jj-real(j1)
  dym = real(j1+1)-obs_jj

  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! Fill the Atm structure array.
  ! NOTE: This is an example program for illustrative purposes only.
  !       Typically, one would not assign the data as shown below,
  !       but rather read it from file
  !
  ! 4a1. Loading Atmosphere and Surface input
  ! --------------------------------
  i_p = 0
  i_pb = 0
  i_ph = 0
  i_phb = 0
  i_pt = 0
  i_qv = 0
  i_qc = 0
  i_qr = 0
  i_qi = 0
  i_qs = 0
  i_qg = 0
  i_psfc = 0
  i_mu = 0
  i_mub = 0
  do m = 1, nv
    if ( enkfvar(m) == 'P         ' ) i_p=m
    if ( enkfvar(m) == 'PB        ' ) i_pb=m
    if ( enkfvar(m) == 'PH        ' ) i_ph=m
    if ( enkfvar(m) == 'PHB       ' ) i_phb=m
    if ( enkfvar(m) == 'T         ' ) i_pt=m
    if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
    if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
    if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
    if ( enkfvar(m) == 'QICE      ' ) i_qi=m
    if ( enkfvar(m) == 'QSNOW     ' ) i_qs=m
    if ( enkfvar(m) == 'QGRAUP    ' ) i_qg=m
    if ( enkfvar(m) == 'PSFC      ' ) i_psfc=m
    if ( enkfvar(m) == 'MU        ' ) i_mu=m
    if ( enkfvar(m) == 'MUB       ' ) i_mub=m
  enddo
  if(i_p>0) p(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_p)
  if(i_pb>0) pb(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pb)
  if(i_ph>0) ph(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_ph)
  if(i_phb>0) phb(i1:i1+1,j1:j1+1,1:kx+1)=xa(1:2,1:2,1:kx+1,i_phb)
  if(i_pt>0) pt(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_pt)
  if(i_qv>0) qv(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qv)
  if(i_qc>0) qc(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qc)
  if(i_qr>0) qr(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qr)
  if(i_qi>0) qi(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qi)
  if(i_qs>0) qs(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qs)
  if(i_qg>0) qg(i1:i1+1,j1:j1+1,1:kx)=xa(1:2,1:2,1:kx,i_qg)
  if(i_psfc>0) psfc(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_psfc)
  if(i_mu>0) mu(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mu)
  if(i_mub>0) mub(i1:i1+1,j1:j1+1)=xa(1:2,1:2,1,i_mub)

  if ( i_p ==  0 ) call get_variable3d(inputfile, 'P         ', ix, jx, kx, 1, p )
  if ( i_pb == 0 ) call get_variable3d(inputfile, 'PB        ', ix, jx, kx, 1, pb)
  if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
  if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
  if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )
  if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
  if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
  if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
  if ( i_qi == 0 ) call get_variable3d(inputfile, 'QICE      ', ix, jx, kx,   1, qi )
  if ( i_qs == 0 ) call get_variable3d(inputfile, 'QSNOW     ', ix, jx, kx,   1, qs )
  if ( i_qg == 0 ) call get_variable3d(inputfile, 'QGRAUP    ', ix, jx, kx,   1, qg )
  if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu )
  if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub)

  mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
  mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
  qv1 = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
  qr1 = dym*(dx*qr(i1+1,j1,1:kx) + dxm*qr(i1,j1,1:kx)) + dy*(dx*qr(i1+1,j1+1,1:kx) + dxm*qr(i1,j1+1,1:kx))
  qc1 = dym*(dx*qc(i1+1,j1,1:kx) + dxm*qc(i1,j1,1:kx)) + dy*(dx*qc(i1+1,j1+1,1:kx) + dxm*qc(i1,j1+1,1:kx))
  qi1 = dym*(dx*qi(i1+1,j1,1:kx) + dxm*qi(i1,j1,1:kx)) + dy*(dx*qi(i1+1,j1+1,1:kx) + dxm*qi(i1,j1+1,1:kx))
  qs1 = dym*(dx*qs(i1+1,j1,1:kx) + dxm*qs(i1,j1,1:kx)) + dy*(dx*qs(i1+1,j1+1,1:kx) + dxm*qs(i1,j1+1,1:kx))
  qg1 = dym*(dx*qg(i1+1,j1,1:kx) + dxm*qg(i1,j1,1:kx)) + dy*(dx*qg(i1+1,j1+1,1:kx) + dxm*qg(i1,j1+1,1:kx))
  qvt = qv1+qr1+qc1
  pt1(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))
  ph(i1:i1+1, j1:j1+1, 1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
  ph1(1:kx+1) = dym*(dx*ph(i1+1,j1,1:kx+1) + dxm*ph(i1,j1,1:kx+1)) + dy*(dx*ph(i1+1,j1+1,1:kx+1) + dxm*ph(i1,j1+1,1:kx+1))
  p(i1:i1+1, j1:j1+1, 1:kx) = p(i1:i1+1, j1:j1+1, 1:kx) + pb(i1:i1+1, j1:j1+1, 1:kx)
  pres(1:kx) = dym*(dx*p(i1+1,j1,1:kx) + dxm*p(i1,j1,1:kx)) + dy*(dx*p(i1+1,j1+1,1:kx) + dxm*p(i1,j1+1,1:kx))
  tk = (pt1 + to) * ( (pres / P1000MB) ** (R_D/Cpd) )
  lat = dym*(dx*xlat(i1+1,j1  ) + dxm*xlat(i1,j1  )) + dy*(dx*xlat(i1+1,j1+1) + dxm*xlat(i1,j1+1))
  lon = dym*(dx*xlong(i1+1,j1  ) + dxm*xlong(i1,j1  )) + dy*(dx*xlong(i1+1,j1+1) + dxm*xlong(i1,j1+1))
  lat = lat/180.0*3.14159
  lon = lon/180.0*3.14159
  hgt1 = dym*(dx*hgt(i1+1,j1  ) + dxm*hgt(i1,j1  )) + dy*(dx*hgt(i1+1,j1+1) + dxm*hgt(i1,j1+1))
  tsk1 = dym*(dx*tsk(i1+1,j1  ) + dxm*tsk(i1,j1  )) + dy*(dx*tsk(i1+1,j1+1) + dxm*tsk(i1,j1+1))
  psfc1 = dym*(dx*psfc(i1+1,j1  ) + dxm*psfc(i1,j1  )) + dy*(dx*psfc(i1+1,j1+1) + dxm*psfc(i1,j1+1))
  landmask1 = dym*(dx*landmask(i1+1,j1  ) + dxm*landmask(i1,j1  )) + dy*(dx*landmask(i1+1,j1+1) + dxm*landmask(i1,j1+1))

  where(qv1.lt.0.0) qv1=1.0e-8
  where(qc1.lt.0.0) qc1=0.0
  where(qi1.lt.0.0) qi1=0.0
  where(qr1.lt.0.0) qr1=0.0
  where(qs1.lt.0.0) qs1=0.0
  where(qg1.lt.0.0) qg1=0.0

  if(cloud_flag==0) then  !!set hydrometeors to zero if calculating for clear-sky Tb
    qc1=0.0
    qi1=0.0
    qr1=0.0
    qs1=0.0
    qg1=0.0
  end if


  ! 4a2. Converting WRF data for CRTM structure
  ! --------------------------------
  !--- converting the data to CRTM structure

  !*******************************************************************************
  ! satellite information
  !*******************************************************************************

  sat_dis=sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon-sat_lon)*cos(lat))
  SCAN_ANGLE=180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon-sat_lon)*cos(lat))**2))
  ZENITH_ANGLE=SCAN_ANGLE+180.0/3.14159*acos(cos(lon-sat_lon)*cos(lat))

  !*******************************************************************************
  ! load WRF data into CRTM structures
  !*******************************************************************************
  !--- calcurating delz
   do z=1,kx
    if(z.eq.1) then
     delz(z) = ph1(z+1) / 9.806 - hgt1
    else
     delz(z) = (ph1(z+1)-ph1(z))/2/9.806
    endif
   enddo
  !---Atmospheric Profile
   Atm(1)%Climatology         = TROPICAL
   Atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
   Atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)
   Atm(1)%Level_Pressure(0) = (pres(kx)*3.0/2.0 - pres(kx-1)/2.0)/100.0  ! convert from Pa to hPA
   do z=kx,1,-1
     if(z.eq.1) then
       Atm(1)%Level_Pressure(kx-z+1) = psfc1/100.0  ! convert from Pa tohPA
     else
       Atm(1)%Level_Pressure(kx-z+1) = ((pres(z-1)+pres(z))/2.0)/100.0  ! convert from Pa to hPA
     endif
     Atm(1)%Pressure(kx-z+1)       = pres(z) / 100.0
     Atm(1)%Temperature(kx-z+1)    = tk(z)
     Atm(1)%Absorber(kx-z+1,1)     = qv1(z)*1000.0
   enddo
   Atm(1)%Absorber(:,2) = 5.0E-2
   !---Cloud Profile
   do z=1,kx*5
     Atm(1)%Cloud(z)%Type = 0
     Atm(1)%Cloud(z)%Effective_Radius = 0.0
     Atm(1)%Cloud(z)%Water_Content = 0.0
   enddo
   ncl = 0
   icl = 0
   !--calculating # of clouds (cloud and rain)
   do z=kx,1,-1
     if(qc1(z).gt.0.0) then
       ncl = ncl + 1
     endif
     if(qr1(z).gt.0.0) then
       ncl = ncl + 1
     endif
     if(qi1(z).gt.0.0) then
       ncl = ncl + 1
     endif
     if(qs1(z).gt.0.0) then
       ncl = ncl + 1
     endif
     if(qg1(z).gt.0.0) then
       ncl = ncl + 1
     endif
   enddo
   !--Data for cloud
   Atm(1)%n_Clouds = ncl
   if ( Atm(1)%n_Clouds > 0 ) then
     do z=kx,1,-1
       if(qc1(z).gt.0.0) then
         icl = icl + 1
         k1 = kx-z+1
         k2 = kx-z+1
         Atm(1)%Cloud(icl)%Type = WATER_CLOUD
         Atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
         Atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
             qc1(z)*pres(z)/287.2/(tk(z)+0.61*(qv1(z)/(1+qv1(z))))*delz(z)
       endif
     enddo
     do z=kx,1,-1
       if(qr1(z).gt.0.0) then
         icl = icl + 1
         k1 = kx-z+1
         k2 = kx-z+1
         Atm(1)%Cloud(icl)%Type = RAIN_CLOUD
         Atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
         Atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
             qr1(z)*pres(z)/287.2/(tk(z)+0.61*(qv1(z)/(1+qv1(z))))*delz(z)
       endif
     enddo
     do z=kx,1,-1
       if(qi1(z).gt.0.0) then
         icl = icl + 1
         k1 = kx-z+1
         k2 = kx-z+1
         Atm(1)%Cloud(icl)%Type = ICE_CLOUD
         Atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
         Atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
             qi1(z)*pres(z)/287.2/(tk(z)+0.61*(qv1(z)/(1+qv1(z))))*delz(z)
       endif
     enddo
     do z=kx,1,-1
       if(qs1(z).gt.0.0) then
         icl = icl + 1
         k1 = kx-z+1
         k2 = kx-z+1
         Atm(1)%Cloud(icl)%Type = SNOW_CLOUD
         Atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
         Atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
             qs1(z)*pres(z)/287.2/(tk(z)+0.61*(qv1(z)/(1+qv1(z))))*delz(z)
       endif
     enddo
     do z=kx,1,-1
       if(qg1(z).gt.0.0) then
         icl = icl + 1
         k1 = kx-z+1
         k2 = kx-z+1
         Atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
         Atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
         Atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
             qg1(z)*pres(z)/287.2/(tk(z)+0.61*(qv1(z)/(1+qv1(z))))*delz(z)
!!!!Virtual temperature calcualtion can be shared!!!!
       endif
     enddo
   endif

  !---Surface data
   if(landmask1.eq.1.0) then
     Sfc(1)%Water_Coverage = 0.0_fp
     Sfc(1)%Land_Coverage = 1.0_fp
     Sfc(1)%Land_Temperature = tsk1
     Sfc(1)%Soil_Temperature = tsk1
   else
     Sfc(1)%Water_Coverage = 1.0_fp
     Sfc(1)%Land_Coverage = 0.0_fp
     Sfc(1)%Water_Type = 1  ! Sea water
     Sfc(1)%Water_Temperature = tsk1
   endif

  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )


  ! 4c. Use the SOI radiative transfer algorithm
  ! --------------------------------------------
  Options%RT_Algorithm_ID = RT_SOI

  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution , &
                               Options = Options )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! ============================================================================
  ! 6. **** output ****
  !
  ! User should read the user guide or the source code of the routine
  ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
  ! select the needed variables for outputs.  These variables are contained
  ! in the structure RTSolution.
  !
  !DO m = 1, N_PROFILES
    !WRITE( *,'(//7x,"Profile ",i0," output for ",a )') n, TRIM(Sensor_Id)
    !DO l = 1, n_Channels
      !WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
      !CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
    !END DO
  !END DO

  xb=real(RTSolution(1,1)%Brightness_Temperature)
  if(xb<100 .or. xb>400) xb=-888888.

  ! ============================================================================
  ! 7. **** DESTROY THE CRTM ****
  !
  !WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

end subroutine xb_to_radiance

