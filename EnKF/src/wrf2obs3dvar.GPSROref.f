program wrf2obs3dvar
!add refractivity GPSRO obs to obs3dvar format data

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
real :: obs_ii,obs_jj,obs_kk,obs_pp,obs_tt,obs_qq,obs_nn
real,dimension(:,:,:),allocatable :: u,v,pt,ph,phb,p,pb,zg,qv,qc,qr,qall
real,dimension(:,:),allocatable :: landmask,mu,mub
real,dimension(:),allocatable :: znw0,znu0, pres,qvt,ptt,zgt
real :: dx,dxm,dy,dym,dz,dzm, mu1,mub1, trueu,truev,gridu,gridv,dum
real,dimension(2) :: worku,workv,work
character(len=10) :: rdate,rtime,rzone
integer, dimension(8)  :: rvalues
real :: gaussdev
integer :: gridobs_ks, gridobs_ke, gridobs_int_k, k_levels
integer :: gridobs_js, gridobs_je, gridobs_is, gridobs_ie, gridobs_int_x
integer :: cyg_gridobs_js, cyg_gridobs_je, cyg_gridobs_is, cyg_gridobs_ie, cyg_gridobs_int_x
real :: rh_error,z_error,spd_error,dir_error, t_error,td_error,npc_error
real :: obs_q,obs_pres,obs_temp,mean_qv,gpsref

character (len=14) :: tstr
character (len=80) :: dt

obsfile = 'gpsro.ascii'
obsfile1 = 'gpsro1.ascii'
truthfile='wrfout.truth'

total=20*25
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
allocate(mu(ix,jx))
allocate(mub(ix,jx))
allocate(znw0(kx+1))
allocate(znu0(kx))
allocate(pres(kx))
allocate(qvt(kx))
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
             do m = k1,k1+1
                work(m-k1+1) = theta_to_temp(ptt(m)+to, pres(m))
                if(qvt(m).eq.-888888) qvt(m)=0.0
                work(m-k1+1) = gpsref(pres(m),work(m-k1+1),qvt(m))
             enddo
             if ( obs_kk .le. 1. ) then
               obs_data(6,1) = work(1)
             else
               obs_data(6,1) = dzm*work(1)+dz*work(2)
             endif

        qcint(6)=0
        obs_data(6,3)=obs_data(6,1)*npc_error(obs_data(1,1))/100
        call date_and_time(rdate, rtime, rzone, rvalues)
        obs_data(6,1) = obs_data(6,1) + obs_data(6,3)*gaussdev(sum(rvalues))

        write(11, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
   enddo
end do
close(10)
close(11)
end program wrf2obs3dvar

!==============================================================================
function npc_error(p)  !!!percentage error for refractivity
real, intent(in) :: p
real, dimension(6) :: npc_err_ref,p_ref 
real npc_error
p_ref=(/100,300,500,700,850,1000/)*100
npc_err_ref=(/0.2,0.2,0.5,0.8,1,1/)
if(p .le. p_ref(1)) npc_error=npc_err_ref(1)
if(p .ge. p_ref(6)) npc_error=npc_err_ref(6)
do i=1,5
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     npc_error = npc_err_ref(i) + (npc_err_ref(i+1)-npc_err_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function npc_error

function gpsref(p,t,q)
  real, intent(in) :: p,t,q
  real pres,qv,ew,gpsref
  real,parameter :: rd=287.05, rv=461.51, c1=77.6d-6, c2=3.73d-1, rdorv=rd/rv
  pres=0.01*p
  ew=q*pres/(rdorv+(1.0-rdorv)*q)
  gpsref=c1*pres/t+c2*ew/(t**2)
  gpsref=gpsref*1.0d6
end function

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

