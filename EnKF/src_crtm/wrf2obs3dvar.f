program wrf2obs3dvar

use netcdf
use mapinfo_define
use map_utils
use wrf_tools
use constants

implicit none
character (len=125)        :: str
character (len=80)         :: obs_3dvar_file,obsfile1,obsfile2,truthfile,times,times1
character (len=180)        :: fmt_fmt, info_fmt, srfc_fmt, each_fmt
integer                    :: synop,metar,ship,buoy,bogus,temp,amdar,airep,tamdar,pilot, &
                              satem, satob, gpspw, gpszd, gpsrf,gpsep,ssmt1,ssmt2,     &
                              tovs, qscat, profl, airsr, other, total
integer :: total1,total2,temp1,temp2,satob1,satob2
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
real :: obs_ii,obs_jj,obs_kk,obs_pp
real,dimension(:,:,:),allocatable :: u,v,pt,ph,phb,p,pb,zg,qv,qc,qr,qall
real,dimension(:,:),allocatable :: mu,mub
real,dimension(:),allocatable :: znw0,znu0, pres,qvt,ptt,zgt
real :: dx,dxm,dy,dym,dz,dzm, mu1,mub1, trueu,truev,gridu,gridv
real,dimension(2) :: worku,workv,work
character(len=10) :: rdate,rtime,rzone
integer, dimension(8)  :: rvalues
real :: gaussdev
integer :: gridobs_ks, gridobs_ke, gridobs_int_k, k_levels
integer :: gridobs_js, gridobs_je, gridobs_is, gridobs_ie, gridobs_int_x
real :: rh_error,z_error,spd_error,dir_error, t_error,td_error,qpc_error
real :: obs_q,obs_pres,obs_temp,mean_qv

character (len=14) :: tstr
character (len=80) :: dt
integer*8 :: t1,t2
real, parameter :: lat_dif =  0.461
real, parameter :: lon_dif =  0.659

obs_3dvar_file = 'obs_3dvar.ascii'
obsfile1 = 'obs1.ascii'
obsfile2 = 'obs2.ascii'
truthfile='wrfout.truth'
gridobs_is=1
gridobs_ie=333
gridobs_js=1
gridobs_je=222
gridobs_ks=1
gridobs_ke=44
gridobs_int_x=10
gridobs_int_k=1

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
call get_variable2d(truthfile,'MU        ', ix, jx, 1,    mu )
call get_variable2d(truthfile,'MUB       ', ix, jx, 1,    mub)
call get_variable1d(truthfile,'ZNW       ', kx+1, 1, znw0)
call get_variable1d(truthfile,'ZNU       ', kx  , 1, znu0)

call get_times(truthfile,'T',times)

!read in obs_3dvar 1st pass: count observations
open(10,file=obs_3dvar_file,status='old',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ',obs_3dvar_file

read(10, fmt='(7x,i7)')total
read(10, fmt='(6(7x,i7,2x))') synop, metar, ship, buoy, bogus, temp
read(10, fmt='(6(7x,i7,2x))') amdar, airep, tamdar,pilot,satem, satob
read(10, fmt='(6(7x,i7,2x))') gpspw, gpszd, gpsrf, gpsep, ssmt1, ssmt2
read(10, fmt='(5(7x,i7,2x))') tovs,qscat, profl, airsr, other
synop=0; metar=0; ship=0; buoy=0; bogus=0; temp=0
amdar=0; airep=0; tamdar=0; pilot=0; satem=0; satob=0
gpspw=0; gpszd=0; gpsrf=0; gpsep=0; ssmt1=0; ssmt2=0
tovs=0; qscat=0; profl=0; airsr=0; other=0
total1=0; total2=0
temp1=0; temp2=0
satob1=0; satob2=0

do i = 1, 12
   read(10,*)
enddo
read(10,'(a180)')fmt_fmt
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
info_fmt(1:180)=' '
info_fmt(1:j-i)=fmt_fmt(i+1:j)
read(10,'(a180)')fmt_fmt
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
srfc_fmt(1:180)=' '
srfc_fmt(1:j-i)=fmt_fmt(i+1:j)
read(10,'(a180)')fmt_fmt
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
each_fmt(1:180)=' '
each_fmt(1:j-i)=fmt_fmt(i+1:j)
read(10,*)

do n = 1, total
   read(10, fmt=info_fmt) platform,date,name,obs_level,latitude,longitude,elevation,id
!!!
latitude = latitude + lat_dif
longitude = longitude + lon_dif
!!!
   read(10, fmt=srfc_fmt) slp(1),qcint(1),slp(3),pw(1),qcint(2),pw(3)
   tstr=date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)//date(18:19)
   read(tstr,'(i14)') t1
   tstr=times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16)//times(18:19)
   read(tstr,'(i14)') t2
   if ( platform(4:6).eq.'88 ' ) then
!     if(t1.lt.t2) then
     call latlon_to_ij(proj,latitude,longitude,obs_ii,obs_jj)
     if ( obs_ii > 2 .and. obs_ii < proj%nx-1 .and. obs_jj > 2 .and. obs_jj<proj%ny-1) then
       total1=total1+1
       satob1=satob1+1
     endif
!     else
!       total2=total2+1
!       satob2=satob2+1
!     endif
   endif
   if ( platform(4:6).eq.'35 ' ) then
!     if(t1.lt.t2) then
!       total1=total1+1
!       temp1=temp1+1
!     else
!       total2=total2+1
!       temp2=temp2+1
!     endif
   endif
   do k = 1, obs_level
     read(10,*)
   enddo
enddo
close(10)


!write wrf output observations to the obs_3dvar format
open(11,file=obsfile1,status='replace',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ',obsfile1
open(12,file=obsfile2,status='replace',form='formatted',iostat=iost)
if(iost .ne. 0) write(*,*) 'error opening ',obsfile2

write(11,fmt='(a,i7,a,f8.0,a)') 'TOTAL =', total1, ', MISS. =',-888888.,','
write(11, fmt='(6(a,i7,a))') 'SYNOP =',synop,', ','METAR =',metar,', ','SHIP  =',ship,', ',&
                             'BUOY  =',buoy, ', ','BOGUS =',bogus,', ','TEMP  =',temp1,', '
write(11, fmt='(6(a,i7,a))') 'AMDAR =',amdar,', ','AIREP =',airep,', ','TAMDAR=',tamdar,', ',&
                             'PILOT =',pilot,', ','SATEM =',satem,', ','SATOB =',satob1,', '
write(11, fmt='(6(a,i7,a))') 'GPSPW =',gpspw,', ','GPSZD =',gpszd,', ','GPSRF =',gpsrf,', ',&
                             'GPSEP =',gpsep,', ','SSMT1 =',ssmt1,', ','SSMT2 =',ssmt2,', '
write(11, fmt='(5(a,i7,a))') 'TOVS  =',tovs, ', ','QSCAT =',qscat,', ','PROFL =',profl,', ',&
                             'AIRSR =',airsr,', ','OTHER =',other,', '

write(12,fmt='(a,i7,a,f8.0,a)') 'TOTAL =', total2, ', MISS. =',-888888.,','
write(12, fmt='(6(a,i7,a))') 'SYNOP =',synop,', ','METAR =',metar,', ','SHIP  =',ship,', ',&
                             'BUOY  =',buoy, ', ','BOGUS =',bogus,', ','TEMP  =',temp2,', '
write(12, fmt='(6(a,i7,a))') 'AMDAR =',amdar,', ','AIREP =',airep,', ','TAMDAR=',tamdar,', ',&
                             'PILOT =',pilot,', ','SATEM =',satem,', ','SATOB =',satob2,', '
write(12, fmt='(6(a,i7,a))') 'GPSPW =',gpspw,', ','GPSZD =',gpszd,', ','GPSRF =',gpsrf,', ',&
                             'GPSEP =',gpsep,', ','SSMT1 =',ssmt1,', ','SSMT2 =',ssmt2,', '
write(12, fmt='(5(a,i7,a))') 'TOVS  =',tovs, ', ','QSCAT =',qscat,', ','PROFL =',profl,', ',&
                             'AIRSR =',airsr,', ','OTHER =',other,', '

open(10,file=obs_3dvar_file,status='old',form='formatted',iostat=iost)
!same header as obs_3dvar
do i=1,5
  read(10,*) 
end do
do i=1,16
  read(10,'(a)') str
  write(11,'(a)') str
  write(12,'(a)') str
end do

!write observations
do n=1,total
   read(10, fmt=info_fmt)platform,date,name,obs_level,latitude,longitude,elevation,id
   read(10, fmt=srfc_fmt)slp(1),qcint(1),slp(3),pw(1),qcint(2),pw(3)
   tstr=date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)//date(18:19)
   read(tstr,'(i14)') t1
   tstr=times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16)//times(18:19)
   read(tstr,'(i14)') t2
!!!
latitude = latitude + lat_dif
longitude = longitude + lon_dif
!!!

   if ( platform(4:6)=='88 ') then ! .or. platform(4:6)=='35 ') then
     call latlon_to_ij(proj,latitude,longitude,obs_ii,obs_jj)
     if ( obs_ii > 2 .and. obs_ii < proj%nx-1 .and. obs_jj > 2 .and. obs_jj <proj%ny-1) then
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
!       if (platform(4:6)=='35 ' .or. platform(4:6)=='88 ' ) then  !screening
!         if(slp(1).ne.-888888.) then !slp
!         endif
!         if(pw(1).ne.-888888.) then !pw
!         endif
!       endif
!       if(t1.lt.t2) then
         ounit=11
!       else
!         ounit=12
!       endif
       write(ounit,fmt=info_fmt)platform,date,name,obs_level,latitude,longitude,elevation,id
       write(ounit, fmt=srfc_fmt)slp(1),qcint(1),slp(3),pw(1),qcint(2),pw(3)
     endif
   endif

   do k = 1, obs_level
     read(10, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
     if ( platform(4:6)=='88 ') then ! .or. platform(4:6)=='35 ' ) then    !screening
      if ( obs_ii > 2 .and. obs_ii < proj%nx-1 .and. obs_jj > 2 .and. obs_jj <proj%ny-1) then
        call to_zk(obs_data(1,1), pres(1:kx), obs_kk, kx)
        if ( obs_kk .lt. 1. ) obs_kk = 1.
        k1  = int( obs_kk )
        dz  = obs_kk-real(k1)
        dzm = real(k1+1)-obs_kk
        if (obs_data(1,1).ne.-888888. .and. obs_kk.ge.1. .and. obs_kk.le.real(kx)) then
          if (obs_data(2,1).ne.-888888.) then  !wind
             if(obs_ii-i1==0.5) then
               worku(1:2)=dym*u(i1+1,j1,k1:k1+1) + dy*u(i1+1,j1+1,k1:k1+1)
             elseif(obs_ii-i1>0.5) then
               worku(1:2) = dym*( u(i1+1,j1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1,k1:k1+1)*(dx-0.5) ) + &
                            dy *( u(i1+1,j1+1,k1:k1+1)*(dxm+0.5)+u(i1+2,j1+1,k1:k1+1)*(dx-0.5) )
             elseif(obs_ii-i1<0.5) then
               worku(1:2) = dym*( u(i1,j1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1,k1:k1+1)*(dx+0.5) ) + &
                            dy *( u(i1,j1+1,k1:k1+1)*(dxm-0.5)+u(i1+1,j1+1,k1:k1+1)*(dx+0.5) )
             endif
             if(obs_jj-j1==0.5) then
               workv(1:2) = v(i1,j1+1,k1:k1+1)*dxm + v(i1+1,j1+1,k1:k1+1)*dx
             elseif(obs_jj-j1>0.5) then
               workv(1:2) = dxm*( v(i1,j1+1,k1:k1+1)*(dym+0.5) +v(i1,j1+2,k1:k1+1)*(dy-0.5) ) + &
                            dx *( v(i1+1,j1+1,k1:k1+1)*(dym+0.5) +v(i1+1,j1+2,k1:k1+1)*(dy-0.5) )
             elseif(obs_jj-j1<0.5) then
               workv(1:2) = dxm*( v(i1,j1,k1:k1+1)*(dym-0.5) +v(i1,j1+1,k1:k1+1)*(dy+0.5) ) + &
                            dx *( v(i1+1,j1,k1:k1+1)*(dym-0.5) +v(i1+1,j1+1,k1:k1+1)*(dy+0.5) )
             endif
             if(obs_kk.le.1.) then
               gridu=worku(1); gridv=workv(1)
             else
               gridu=dzm*worku(1)+dz*worku(2)
               gridv=dzm*workv(1)+dz*workv(2)
             endif
             call gridwind_to_truewind(longitude,proj,gridu,gridv,trueu,truev)
             call date_and_time(rdate, rtime, rzone, rvalues)
             trueu = trueu + obs_data(2,3)*gaussdev(sum(rvalues))
             truev = truev + obs_data(2,3)*gaussdev(sum(rvalues))
             call uv_to_dirspd(trueu,truev,obs_data(3,1),obs_data(2,1))
          endif
          if (obs_data(5,1).ne.-888888.) then  !T
             do m = k1,k1+1
                work(m-k1+1) = theta_to_temp(ptt(m)+to, pres(m))
             enddo
             if ( obs_kk .le. 1. ) then
               obs_data(5,1) = work(1)
             else
               obs_data(5,1) = dzm*work(1)+dz*work(2)
             endif
             call date_and_time(rdate, rtime, rzone, rvalues)
             obs_data(5,1) = obs_data(5,1) + obs_data(5,3)*gaussdev(sum(rvalues))
          endif
          if (obs_data(6,1).ne.-888888.) then  !TD
             do m = k1,k1+1
                work(m-k1+1) = mixrat_to_tdew(qvt(m), pres(m))
             enddo
             if ( obs_kk .le. 1. ) then
               obs_data(6,1) = work(1)
             else
               obs_data(6,1) = dzm*work(1)+dz*work(2)
             endif
             call date_and_time(rdate, rtime, rzone, rvalues)
             obs_data(6,1) = obs_data(6,1) + obs_data(6,3)*gaussdev(sum(rvalues))
          endif
          if (obs_data(7,1).ne.-888888.) then  !RH
             do m = k1,k1+1
                work(m-k1+1) = theta_to_temp(ptt(m)+to, pres(m))
                work(m-k1+1) = rel_humidity(qvt(m),work(m-k1+1), pres(m))
             enddo
             if ( obs_kk .le. 1. ) then
               obs_data(7,1) = work(1)
             else
               obs_data(7,1) = dzm*work(1)+dz*work(2)
             endif
             call date_and_time(rdate, rtime, rzone, rvalues)
             obs_data(7,1) = obs_data(7,1) + obs_data(7,3)*gaussdev(sum(rvalues))
             if(obs_data(7,1)<0.) obs_data(7,1)=0.
             if(obs_data(7,1)>100.) obs_data(7,1)=100.
          endif
        else
          obs_data(:,1)=-888888.
        endif
        write(ounit, fmt=each_fmt)((obs_data(i,1),qcint(i),obs_data(i,3)),i=1,7)
      endif
     endif
   enddo
end do
close(10)


close(11)
close(12)
end program wrf2obs3dvar

!==============================================================================
function z_error(p)
real, intent(in) :: p
real, dimension(10) :: zerr_ref,p_ref 
real z_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
zerr_ref=(/40,33,28,25,19,15,12,9,8,7/)
if(p .le. p_ref(1)) z_error=zerr_ref(1)
if(p .ge. p_ref(10)) z_error=zerr_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     z_error = zerr_ref(i) + (zerr_ref(i+1)-zerr_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function z_error

function spd_error(p)
real, intent(in) :: p
real, dimension(10) :: err_ref,p_ref 
real spd_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
err_ref=(/2.7,3.0,3.3,3.3,3.3,2.8,2.3,1.4,1.1,1.1/)
if(p .le. p_ref(1)) spd_error=err_ref(1)
if(p .ge. p_ref(10)) spd_error=err_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     spd_error = err_ref(i) + (err_ref(i+1)-err_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function spd_error

function dir_error(p)
real, intent(in) :: p
real, dimension(10) :: err_ref,p_ref 
real dir_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
err_ref=(/15,15,15,16,16,16,17,18,18,20/)
if(p .le. p_ref(1)) spd_error=err_ref(1)
if(p .ge. p_ref(10)) spd_error=err_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     dir_error = err_ref(i) + (err_ref(i+1)-err_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function dir_error

function t_error(p)
real, intent(in) :: p
real, dimension(10) :: terr_ref,p_ref
real t_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
terr_ref=(/1.8,1.8,1.8,1.7,1.7,1.3,1.7,1.7,2.5,2.6/)
if(p .le. p_ref(1)) t_error=terr_ref(1)
if(p .ge. p_ref(10)) t_error=terr_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     t_error = terr_ref(i) + (terr_ref(i+1)-terr_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function t_error

function td_error(p)
real, intent(in) :: p
real, dimension(10) :: terr_ref,p_ref
real td_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
terr_ref=(/3.0,3.0,3.0,3.0,4.0,4.0,4.0,5.0,4.0,3.5/)
if(p .le. p_ref(1)) td_error=terr_ref(1)
if(p .ge. p_ref(10)) td_error=terr_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     td_error = terr_ref(i) + (terr_ref(i+1)-terr_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function td_error

function rh_error(p)
real, intent(in) :: p
real, dimension(10) :: rherr_ref,p_ref 
real rh_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
rherr_ref=(/10,10,10,10,10,10,10,10,10,15/)
if(p .le. p_ref(1)) rh_error=rherr_ref(1)
if(p .ge. p_ref(10)) rh_error=rherr_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     rh_error = rherr_ref(i) + (rherr_ref(i+1)-rherr_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function rh_error

function qpc_error(p)
real, intent(in) :: p
real, dimension(10) :: qpc_err_ref,p_ref 
real qpc_error
p_ref=(/100,150,200,250,300,400,500,700,850,1000/)*100
qpc_err_ref=(/30,35,40,55,45,32,25,30,27,27/)
if(p .le. p_ref(1)) qpc_error=qpc_err_ref(1)
if(p .ge. p_ref(10)) qpc_error=qpc_err_ref(10)
do i=1,9
  if( p .gt. p_ref(i) .and. p .lt. p_ref(i+1) ) &
     qpc_error = qpc_err_ref(i) + (qpc_err_ref(i+1)-qpc_err_ref(i))*(p-p_ref(i))/(p_ref(i+1)-p_ref(i))
enddo
return
end function qpc_error

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

