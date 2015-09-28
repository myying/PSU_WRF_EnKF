!=======================================================================================
subroutine xb_to_surface(inputfile,proj,xa,ix,jx,kx,nv,iob,xland,lu_index,znu,znw,p_top,times,xb)
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
real, dimension(ix, jx ), intent(in)    :: xland, lu_index
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
real                                    :: mu1, mub1, long, grid_u, grid_v
integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk, obs_ii,obs_jj
integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_var, i_u, i_v
integer                                 :: i_t2, i_th2, i_q2, i_u10, i_v10
real                                    :: dx, dxm, dy, dym, dz, dzm
real                                    :: psfc, tsfc, psfcm, u10, v10, t2, q2, th2
real, dimension(ix, jx      )           :: rough
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
     xb = sqrt(u10**2.+v10**2.)
endif

!if( print_detail > 100 )write(*,'(3x,a,3f)')'xb_to_surface '//obstype//' obs xb, obs_ii, obs_jj :', xb, obs_ii, obs_jj

end subroutine xb_to_surface

!=======================================================================================
   subroutine xb_to_sounding (inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,xb,itruth,isimulated )
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
   integer, intent(in)                     :: itruth, isimulated

   real :: obs_ii, obs_jj, obs_kk, obs_pp !
   real, intent(out)   :: xb
   character(len=10)   :: obstype

   real, dimension(ix, jx, kx+1)           :: ph, phb
   real, dimension(ix, jx, kx  )           :: pt, pb, qv, qc, qr
   real, dimension(ix, jx      )           :: mu, mub
   real, dimension(ix+1, jx, kx)           :: u
   real, dimension(ix, jx+1, kx)           :: v
   real, dimension(kx+1)                   :: znw0
   real, dimension(kx)                     :: znu0

   real, dimension(2, 2, kx+1)             :: p
   real, dimension(kx+1)                   :: work, worku, workv

   real, dimension(kx)                     :: pres, ptt, qvt
   real                                    :: mu1, mub1, long, grid_u, grid_v, true_u, true_v
   real                                    :: dx, dxm, dy, dym, dz, dzm
   integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk
   integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_var, i_u, i_v

   obstype=obs%type(iob)
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   obs_pp=obs%position(iob,4)
   i1 = int( obs_ii )
   j1 = int( obs_jj ) 
   dx  = obs_ii-real(i1)
   dxm = real(i1+1)-obs_ii
   dy  = obs_jj-real(j1)
   dym = real(j1+1)-obs_jj

!.. calculate pressure from geopotential height
    if ( itruth == 0 )then
         call get_variable1d(inputfile, 'ZNW       ', kx+1, 1, znw0)
         call get_variable1d(inputfile, 'ZNU       ', kx  , 1, znu0)
    else
         znw0 = znw 
         znu0 = znu
    endif 

!.. Calculate obs%position(iob,3)
    if ( isimulated == 0 .or. obstype(10:10) == 'T' .or. obstype(10:10) == 'D' .or.  &
                              obstype(10:10) == 'R' .or. obstype(10:10) == 'Q' ) then
   
!..    get data from background
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
       if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )
       if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
       if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
       if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
       if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu )
       if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub)

       qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)

!...... get mut mu, mub at obs' position from horizontal interpolation
       mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
       mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))

!........ get qvt (qvapor) at obs' position from horizontal interpolation
!........ get ptt(theta)  at obs' position from horizontal interpolation
       qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
       ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))

       p(1:2,1:2,1:kx+1) = ph(i1:i1+1, j1:j1+1, 1:kx+1) + phb(i1:i1+1, j1:j1+1, 1:kx+1)
       p(1,1,1:kx+1) = dym*(dx*p(2,1,1:kx+1) + dxm*p(1,1,1:kx+1)) + dy*(dx*p(2,2,1:kx+1) + dxm*p(1,2,1:kx+1))

       call eta_to_pres(znw0(1:kx+1), mu1+mub1, qvt(1:kx), p(1,1,1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))
!   if( print_detail==3 )write(*,'(3x,a,4f)')'obs location =', obs%position(iob,1:4)
!   if( print_detail==3 )write(*,'(3x,a,f)')'          mu =', mu(1,1)
!   if( print_detail==3 )write(*,'(3x,a,10f)')'          p  =', p(1,1,1:10)
!       call cal_press_from_q( kx, znu0, znw0, qvt, mu1, mub1, p_top, pres )
!       if( print_detail > 100 ) write(*,'(a,100f10.0)')'Pressure = ', pres
!       if( print_detail > 100 ) write(*,'(a,100f10.4)')'qvt = ', qvt

       if ( isimulated == 0 ) then
          if ( obstype(1:1) == 'H' ) then
!........... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
             call to_zk(obs%position(iob,4), p(1,1,1:kx)/g, obs%position(iob,3), kx) 
          else if ( obstype(1:1) == 'P' ) then
             call to_zk(obs%position(iob,4), pres(1:kx), obs%position(iob,3), kx)
          endif
!          if ( obs%position(iob,3) .ge. real(kx-1) ) obs%position(iob,3) = real(kx-1)
          if ( obs%position(iob,3) .lt. 1. ) obs%position(iob,3) = 1.
          obs_kk = obs%position(iob,3)
       endif
    endif

!.. Get xb from background
   if ( obs_kk .gt. real(kx-1) .or. obs_kk .lt. 1. ) then
      xb = -999999.
      return
   endif

   k1  = int( obs_kk )
   dz  = obs_kk-real(k1)
   dzm = real(k1+1)-obs_kk

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

       if ( obs_kk .le. 1. ) then
          grid_u = worku(k1)
          grid_v = workv(k1)
       else
          grid_u = dzm*worku(k1)+dz*worku(k1+1)
          grid_v = dzm*workv(k1)+dz*workv(k1+1)
       endif
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
          xb = work(k1)
       else
          xb = dzm*work(k1)+dz*work(k1+1)
       endif

!.. TD(if used TD, not RH)
    else if ( obstype(10:10) == 'D' ) then
       do k = k1, k1+1
          work(k) = mixrat_to_tdew(qvt(k), pres(k))
       enddo
       if ( obs_kk .le. 1. ) then
          xb = work(k1)
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
          xb = work(k1)
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
       call destag_zstag(znu0, znw0, kx, p(1,1,1:kx+1), work(1:kx))
       if ( obs_kk .le. 1. ) then
          xb = work(k1)/g
       else
          xb = (dzm*work(k1) + dz*work(k1+1))/g
       endif 

    endif

!   if( print_detail > 100 )write(*,'(a,3f8.2, f10.1)')'xb_to_sounding '//obstype//' obs position :', obs_ii, obs_jj, obs%position(iob,3:4)
    
   end subroutine xb_to_sounding
!
!=======================================================================================
   subroutine xb_to_idealsound(inputfile,xa,ix,jx,kx,nv,iob,xb)
   use constants
   use namelist_define
   use obs_define
   use netcdf
   use wrf_tools
   use map_utils

   implicit none

   character(len=10), intent(in)           :: inputfile
   real, dimension(3,3,kx+1,nv), intent(in) :: xa
   integer, intent(in) :: ix, jx, kx, nv, iob
   real :: obs_ii, obs_jj, obs_kk
   character(len=10) :: obstype
   real, intent(out) :: xb
   integer :: i1, j1, k1, i, j, k, m, ii, jj, kk

   obstype=obs%type(iob)
   obs_ii=obs%position(iob,1)
   obs_jj=obs%position(iob,2)
   obs_kk=obs%position(iob,3)
   i1 = int( obs_ii )
   j1 = int( obs_jj )
   k1 = int( obs_kk )

   if (obstype .eq. 'idealU    ' ) then
      do m = 1, nv
         if ( enkfvar(m) == 'U         ' ) then
            xb = xa(1,1,k1,m)
         endif
      enddo
   else if (obstype .eq. 'idealV    ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'V         ' ) then
            xb = xa(1,1,k1,m)
         endif
      enddo
   else if (obstype .eq. 'idealPT   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'T         ' ) then
            xb = xa(1,1,k1,m)
         endif
      enddo
   else if (obstype .eq. 'idealQV   ' ) then
      do m = 1, nv
         if (enkfvar(m) == 'QVAPOR    ' ) then
            xb = xa(1,1,k1,m)*1000.
         endif
      enddo
   else if (obstype .eq. 'idealPH   ' ) then
      xb = 0.
      do m = 1, nv
         if ( enkfvar(m) == 'PH        ' .or. enkfvar(m) == 'PHB       ' ) then
            xb = xb + xa(1,1,k1,m)
         endif
      enddo
   endif
   end subroutine xb_to_idealsound

!=======================================================================================
  subroutine xb_to_rv(inputfile,proj,xa,ix,jx,kx,nv,iob,xlong,znw,xb,kkflag)
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
  integer, intent(in)                  :: kkflag
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
  
  if ( kkflag == 0 ) then
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
     call calc_radar_data_point_kk ( ix, jx, kx, ph, phb, obs_ii, obs_jj, obs_hh, obs%position(iob,3) )
     obs_kk=obs%position(iob,3)
!     if (print_detail > 1000) write(*,'(a, 4f8.2)')'obs_ijhk =', obs_ii, obs_jj, obs_hh, obs_kk
  endif

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
