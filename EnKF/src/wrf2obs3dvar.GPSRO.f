program wrf2obs3dvar

use netcdf
use mapinfo_define
use map_utils
use wrf_tools
use constants

implicit none
character (len=125)        :: str
character (len=80)         :: obsfile,truthfile,times
character (len=180)        :: fmt_fmt, info_fmt, srfc_fmt, each_fmt
integer :: total
character (len=40)         :: name
character (len=40)         :: id
integer, dimension(7)      :: qcint
integer                    :: i,j,k,n,m,iost,obs_level,ounit=10
character(len=12)          :: platform
character(len=19)          :: date
real                       :: latitude, longitude, elevation

type(proj_info) :: proj
integer :: ix,jx,kx, i1,j1,k1, k_levels
real :: obs_ii,obs_jj,obs_kk,obs_pp
real,dimension(:,:,:),allocatable :: u,v,pt,ph,phb,p,pb,zg,qv,qc,qr,qall
real,dimension(:,:),allocatable :: landmask,mu,mub
real,dimension(:),allocatable :: znw0,znu0, pres,ht,qvt,qallt,ptt,zgt,zlevels
real :: dx,dxm,dy,dym,dz,dzm, mu1,mub1, trueu,truev,gridu,gridv,dum
real,dimension(2) :: worku,workv,work
character(len=10) :: rdate,rtime,rzone
integer, dimension(8)  :: rvalues
real :: gaussdev,ran3
real :: z_error,t_error,q_error
real :: obs_q,obs_pres,obs_temp,mean_qv
real,dimension(7,3) :: obs_data

obsfile = 'gpsro.ascii'
truthfile='wrfout.truth'

total=20*25 !number of obs in domain
k_levels=85

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
allocate(mu(ix,jx),mub(ix,jx))
allocate(znw0(kx+1),znu0(kx))
allocate(pres(kx),ht(kx),qvt(kx),qallt(kx),ptt(kx))
allocate(zgt(kx+1))
allocate(zlevels(k_levels))

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

!vertical levels avail: GPSRO more vertical levels than ATOVS
! GPSRO ~100 m resolution at Lower troposphere and ~1 km in stratosphere
zlevels(1:10)=(/(i,i=50,500,50)/)
zlevels(11:55)=(/(i,i=600,5000,100)/)
zlevels(56:71)=(/(i,i=5400,11600,400)/)
zlevels(72:85)=(/(i,i=12000,25000,1000)/)

!open obsfile for writing
open(ounit,file=obsfile,status='replace',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ', obsfile

info_fmt='(A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A40)'
srfc_fmt='(F12.3,I4,F7.2,F12.3,I4,F7.3)'
each_fmt='(3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2))'

do n=1,total

  !randomly assign i j location in domain
  call date_and_time(rdate, rtime, rzone, rvalues)
  obs_ii=ran3(sum(rvalues))*(ix-1-5)+1+5
  call date_and_time(rdate, rtime, rzone, rvalues)
  obs_jj=ran3(sum(rvalues))*(jx-1-5)+1+5
  
  i1 = int( obs_ii )
  j1 = int( obs_jj )
  dx  = obs_ii-real(i1)
  dxm = real(i1+1)-obs_ii
  dy  = obs_jj-real(j1)
  dym = real(j1+1)-obs_jj
  mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
  mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
  qallt(1:kx) = dym*(dx*qall(i1+1,j1,1:kx) + dxm*qall(i1,j1,1:kx)) + dy*(dx*qall(i1+1,j1+1,1:kx) + dxm*qall(i1,j1+1,1:kx))
  qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
  ptt(1:kx) = dym*(dx*pt(i1+1,j1,1:kx) + dxm*pt(i1,j1,1:kx)) + dy*(dx*pt(i1+1,j1+1,1:kx) + dxm*pt(i1,j1+1,1:kx))
  zgt(1:kx+1) = dym*(dx*zg(i1+1,j1,1:kx+1) + dxm*zg(i1,j1,1:kx+1)) + dy*(dx*zg(i1+1,j1+1,1:kx+1) + dxm*zg(i1,j1+1,1:kx+1))
  ht(1:kx) = 0.5*(zgt(1:kx)+zgt(2:kx+1))/g
  call eta_to_pres(znw0(1:kx+1), mu1+mub1, qallt(1:kx), zgt(1:kx+1), ptt(1:kx)+to, kx, pres(1:kx))

  call ij_to_latlon(proj,obs_ii,obs_jj,latitude,longitude)
  write(ounit,fmt=info_fmt) 'FM-132 GPSRO',times,'Synthetic GPSRO sounding', &
       k_levels,latitude,longitude,-888888.,'IDEAL                            '
  write(ounit,fmt=srfc_fmt)-888888.000,-88,200.00,-888888.000,-88,0.200

  do k = 1, k_levels
     obs_data(:,1)=-888888.
     qcint(:)=-88
     obs_data(:,3)=0.0

     !pres, height
     obs_pres=interp_hght(pres(1:kx),ht,zlevels(k),kx)
     obs_data(1,1)=obs_pres
     qcint(1)=0
     obs_data(1,3)=100.0
     obs_data(4,1)=zlevels(k)
     qcint(4)=0
     obs_data(4,3)=0.0

     call to_zk(obs_pres, pres(1:kx), obs_kk, kx)
     if (obs_pres.eq.-888888. .or. obs_kk.lt.1. .or. obs_kk.gt.real(kx)) cycle

     k1  = int( obs_kk )
     dz  = obs_kk-real(k1)
     dzm = real(k1+1)-obs_kk

     !T
     do m = k1,k1+1
        work(m-k1+1) = theta_to_temp(ptt(m)+to, pres(m))
     enddo
     qcint(5)=0
     obs_data(5,1) = dzm*work(1)+dz*work(2)
     obs_data(5,3)=t_error(obs_pres)
     call date_and_time(rdate, rtime, rzone, rvalues)
     obs_data(5,1) = obs_data(5,1) + obs_data(5,3)*gaussdev(sum(rvalues))

     !Qv
     if(obs_pres.ge.20000) then
       do m = k1,k1+1
         work(m-k1+1) = 1000*qvt(m)
       enddo
       qcint(7)=0
       obs_data(7,1)=dzm*work(1)+dz*work(2)
       obs_data(7,3)=q_error(obs_pres)
       call date_and_time(rdate, rtime, rzone, rvalues)
       obs_data(7,1) = obs_data(7,1) + obs_data(7,3)*gaussdev(sum(rvalues))
     end if

     write(ounit, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
  enddo

end do

close(ounit)
end program wrf2obs3dvar



!==============================================================================
function t_error(p)
real, intent(in) :: p
real, dimension(5) :: terr_ref,p_ref
real t_error
p_ref=(/100,400,700,850,1000/)*100
terr_ref=(/1.5,1.2,1.5,2.0,2.2/)
if(p .le. p_ref(1)) t_error=terr_ref(1)
if(p .ge. p_ref(5)) t_error=terr_ref(5)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     t_error = terr_ref(i) + (terr_ref(i+1)-terr_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function t_error

function q_error(p)
real, intent(in) :: p
real, dimension(6) :: q_err_ref,p_ref 
real q_error
p_ref=(/200,300,400,500,700,850/)*100
q_err_ref=(/0.03,0.1,0.2,0.5,1.0,1.5/)
if(p .le. p_ref(1)) q_error=q_err_ref(1)
if(p .ge. p_ref(6)) q_error=q_err_ref(6)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     q_error = q_err_ref(i) + (q_err_ref(i+1)-q_err_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function q_error

 function gaussdev(idum)
! Returns a normally distributed deviate with 0 mean and unit variance,
! using ran3(idum) as the source of uniform deviates; see Num'l Recipes,
! ch 7.2.
      integer idum
      real gaussdev
      integer iset
      real fac,gset,rsq,v1,v2,ran1
      save iset, gset
      data iset/0/
      if (iset.eq.0) then
 10      v1 = 2.*ran3(idum)-1.
         v2 = 2.*ran3(idum)-1.
         rsq = v1**2 + v2**2
         if (rsq.ge.1. .or. rsq.eq.0.) goto 10
         fac = sqrt( -2.*log(rsq)/rsq )
         gset = v1*fac
         gaussdev = v2*fac
         iset= 1
      else
         gaussdev = gset
         iset=0
      end if
      return
 end function gaussdev

 function ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
 end function ran3

!==============================================================================

