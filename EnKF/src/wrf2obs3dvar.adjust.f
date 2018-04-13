program wrf2obs3dvar
!parse ATOVS obs and increase the observation error 1.5 times

use netcdf
use mapinfo_define
use map_utils
use wrf_tools
use constants

implicit none
character (len=125)        :: str
character (len=80)         :: obsfile,obsfile1,truthfile,times,times1
character (len=180)        :: fmt_fmt, info_fmt, srfc_fmt, each_fmt
integer                    :: synop,metar,ship,buoy,bogus,temp,amdar,airep,tamdar,pilot, &
                              satem, satob, gpspw, gpszd, gpsrf,gpsep,ssmt1,ssmt2,     &
                              tovs, qscat, profl, airsr, other, total
integer :: total1,total2,temp1,temp2,qscat1,qscat2,satob1,satob2
character (len=40)         :: name
character (len=40)         :: id
integer, dimension(7)      :: qcint
integer                    :: i,j,k,n,m,iost,obs_level,ounit
character(len=12)          :: platform
character(len=19)          :: date
real                       :: latitude, longitude, elevation
real,dimension(3)          :: slp,pw
real, dimension(7,3)       :: obs_data

type(proj_info) :: proj
integer :: ix,jx,kx, i1,j1,k1
real :: obs_ii,obs_jj,obs_kk,truth_tt,truth_td,truth_qq
real,dimension(:,:,:),allocatable :: u,v,pt,ph,phb,p,pb,zg,qv,qc,qr,qall
real,dimension(:,:),allocatable :: landmask,mu,mub
real,dimension(:),allocatable :: znw0,znu0, pres,qvt,ptt,zgt,qv1
real :: dx,dxm,dy,dym,dz,dzm, mu1,mub1, trueu,truev,gridu,gridv,dum
real,dimension(2) :: worku,workv,work
character(len=10) :: rdate,rtime,rzone
integer, dimension(8)  :: rvalues
real :: gaussdev
integer :: gridobs_ks, gridobs_ke, gridobs_int_k, k_levels
integer :: gridobs_js, gridobs_je, gridobs_is, gridobs_ie, gridobs_int_x
integer :: cyg_gridobs_js, cyg_gridobs_je, cyg_gridobs_is, cyg_gridobs_ie, cyg_gridobs_int_x
real :: obs_q,obs_pres,obs_temp,mean_qv

character (len=14) :: tstr
character (len=80) :: dt

obsfile = 'atovs.ascii'
obsfile1 = 'atovs1.ascii'
truthfile='wrfout.truth'

total=418
k_levels=44

!get truth wrfout fields
call get_ij(truthfile,ix,jx,kx)
call set_domain_proj(truthfile,proj)

allocate(u(ix+1,jx,kx))
allocate(v(ix,jx+1,kx))
allocate(pt(ix,jx,kx))
allocate(p(ix,jx,kx))
allocate(pb(ix,jx,kx))
allocate(ph(ix,jx,kx+1))
allocate(phb(ix,jx,kx+1))
allocate(zg(ix,jx,kx+1))
allocate(qall(ix,jx,kx))
allocate(qv(ix,jx,kx))
allocate(qc(ix,jx,kx))
allocate(qr(ix,jx,kx))
allocate(landmask(ix,jx))
allocate(mu(ix,jx))
allocate(mub(ix,jx))
allocate(znw0(kx+1))
allocate(znu0(kx))
allocate(pres(kx))
allocate(qvt(kx))
allocate(qv1(kx))
allocate(ptt(kx))
allocate(zgt(kx+1))

call get_variable3d(truthfile,'T         ', ix, jx, kx,   1, pt)
call get_variable3d(truthfile,'PH        ', ix, jx, kx+1,   1, ph)
call get_variable3d(truthfile,'PHB       ', ix, jx, kx+1,   1, phb)
zg=ph+phb
call get_variable3d(truthfile,'QVAPOR    ', ix, jx, kx,   1, qv)
call get_variable3d(truthfile,'QCLOUD    ', ix, jx, kx,   1, qc)
call get_variable3d(truthfile,'QRAIN     ', ix, jx, kx,   1, qr)
qall=qv+qc+qr
call get_variable3d(truthfile,'U         ', ix+1, jx, kx,   1, u)
call get_variable3d(truthfile,'V         ', ix, jx+1, kx,   1, v)
call get_variable3d(truthfile,'P         ', ix, jx, kx,   1, p)
call get_variable3d(truthfile,'PB        ', ix, jx, kx,   1, pb)
call get_variable2d(truthfile,'LANDMASK  ', ix, jx, 1,    landmask)
call get_variable2d(truthfile,'MU        ', ix, jx, 1,    mu )
call get_variable2d(truthfile,'MUB       ', ix, jx, 1,    mub)
call get_variable1d(truthfile,'ZNW       ', kx+1, 1, znw0)
call get_variable1d(truthfile,'ZNU       ', kx  , 1, znu0)

call get_times(truthfile,'T',times)

!write wrf output observations to the obs_3dvar format
open(10,file=obsfile,status='old',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ',obsfile 
open(11,file=obsfile1,status='replace',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ',obsfile1

info_fmt='(A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A40)'
srfc_fmt='(F12.3,I4,F7.2,F12.3,I4,F7.3)'
each_fmt='(3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2))'


!write observations
do n=1,total
   read(10, fmt=info_fmt)platform,date,name,obs_level,latitude,longitude,elevation,id
   read(10, fmt=srfc_fmt)slp(1),qcint(1),slp(3),pw(1),qcint(2),pw(3)

     call latlon_to_ij(proj,latitude,longitude,obs_ii,obs_jj)
     i1 = int( obs_ii )
     j1 = int( obs_jj )
     dx  = obs_ii-real(i1)
     dxm = real(i1+1)-obs_ii
     dy  = obs_jj-real(j1)
     dym = real(j1+1)-obs_jj
     mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
     mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
     qvt(1:kx) = dym*(dx*qall(i1+1,j1,1:kx) + dxm*qall(i1,j1,1:kx)) + dy*(dx*qall(i1+1,j1+1,1:kx) + dxm*qall(i1,j1+1,1:kx))
     qv1(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
     ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))
     zgt(1:kx+1) = dym*(dx*zg(i1+1,j1,1:kx+1) + dxm*zg(i1,j1,1:kx+1)) + dy*(dx*zg(i1+1,j1+1,1:kx+1) + dxm*zg(i1,j1+1,1:kx+1))
     call eta_to_pres(znw0(1:kx+1), mu1+mub1, qvt(1:kx), zgt(1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))
     write(11,fmt=info_fmt)platform,date,name,obs_level,latitude,longitude,elevation,id
     write(11,fmt=srfc_fmt)slp(1),qcint(1),slp(3),pw(1),qcint(2),pw(3)

   do k = 1, obs_level
     read(10, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
        call to_zk(obs_data(1,1), pres(1:kx), obs_kk, kx)
        if ( obs_kk .lt. 1. ) obs_kk = 1.
        k1  = int( obs_kk )
        dz  = obs_kk-real(k1)
        dzm = real(k1+1)-obs_kk

          if (obs_data(5,1).ne.-888888.) then  !T
             do m = k1,k1+1
                work(m-k1+1) = theta_to_temp(ptt(m)+to, pres(m))
             enddo
             if ( obs_kk .le. 1. ) then
               truth_tt = work(1)
             else
               truth_tt = dzm*work(1)+dz*work(2)
             endif
             obs_data(5,1) = truth_tt + 1.5*(obs_data(5,1)-truth_tt)
             obs_data(5,3) = obs_data(5,3)*1.5
          endif
          if (obs_data(6,1).ne.-888888.) then  !TD
             do m = k1,k1+1
                work(m-k1+1) = mixrat_to_tdew(qvt(m), pres(m))
             enddo
             if ( obs_kk .le. 1. ) then
               truth_td = work(1)
             else
               truth_td = dzm*work(1)+dz*work(2)
             endif
             obs_data(6,1) = truth_td + 1.5*(obs_data(6,1)-truth_td)
             obs_data(6,3) = obs_data(6,3)*1.5
          endif
          if (obs_data(7,1).ne.-888888.) then  !Q
             do m = k1,k1+1
                work(m-k1+1) = 1000*qv1(m)
             enddo
             if ( obs_kk .le. 1. ) then
               truth_qq = work(1)
             else
               truth_qq = dzm*work(1)+dz*work(2)
             endif
             obs_data(7,1) = truth_qq + 1.5*(obs_data(7,1)-truth_qq)
             if(obs_data(7,1).lt.0.0) obs_data(7,1)=0.0
             obs_data(7,3) = obs_data(7,3)*1.5
          endif

        write(11, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
   enddo
end do
close(10)
close(11)
end program wrf2obs3dvar
!==============================================================================

