subroutine enkf (wrf_file,ix,jx,kx,ni,nj,nk,nv,nm,g_comm,s_comm,gid,sid,iid,jid,xom,xo,xm,x,ind,proj,times)
use constants
use namelist_define
use obs_define
use mpi_module
use netcdf
use radar
implicit none
! variable description:
! numbers_en = number of ensemble members (n)
! nicpu, njcpu = number of cpus used in x and y direction. nicpu*njcpu = number of cpus used for a domain
! sid = subdomain id for a cpu. iid, jid= cpu id in decomposed subdomain (coordinate in x and y direction)
! gid = cpugroup id. A cpugroup is cpus with the same gid, but different sid's.
!                    the ensemble members are distributed to each cpugroups in the same way.
! nmcpu = number of cpugroups
! g_comm = MPI_COMM splitted with same sid (all cpus in the same subdomain)
! s_comm = MPI_COMM splitted with same gid (all cpus in the same cpugroup, i.e. covering all subdomains)
! ind : randomly permutated indices
integer, intent(in) :: sid,gid,iid,jid
integer, intent(in) :: ix,jx,kx,ni,nj,nk,nv,nm, g_comm,s_comm
integer, dimension(obs%num), intent(in) :: ind
character (len=10), intent(in) :: wrf_file
type (proj_info)    :: proj
character (len=10)  :: obstype
character (len=80), intent(in) :: times
! fac = 1/(n-1) used for averaging in variance calculation
! corr_coef = localization factor
! ngx, ngz = roi in horizontall and vertical
! var,cov = error variance,covariance of something
! y_hxm = y-hxm (the innovation vector, mean), hxa is the perturbation (of ensemble members).
real      :: fac,d,alpha,beta,var,cov,y_hxm,corr_coef,d_ogn
real      :: la,ka,var_a,var_b,m_d2,m_var_a,m_var_b
real,dimension(numbers_en) :: hxa
integer   :: ngx, ngz, kzdamp
integer   :: i, j, k, m, n, iob, iiob, nob, ie, iunit,ounit,ii, jj, kk, is, it, ig, iv, i1,j1, itot
integer   :: ist,ied,jst,jed,kst,ked, istart,iend,jstart,jend, uist,uied,ujst,ujed
integer   :: sist,sied,sjst,sjed, sistart,siend,sjstart,sjend,sis,sie,sjs,sje
real      :: gaussdev, error, xb
integer, dimension(8)  :: values
character (len=10) :: filename, date, time, zone, varname
character (len=20) :: format1
double precision :: timer,t0,t1,t2,t3,t4,t5
integer :: grid_id, fid, rcode, num_update_var, assimilated_obs_num, update_flag
real                    :: p_top
real, dimension(ix, jx) :: xlong, xlat, xland, lu_index, hgt, tsk
real, dimension(kx+1)   :: znw
real, dimension(kx)     :: znu
real, dimension(3)      :: center_xb
! state vectors (x): xm=ensemble mean; x=perturbation from ensemble mean
! _f=first guess; _t=truth
! ya=Hxa,Hxm: the state vector translated to observation space. ya(number of obs, number of members + mean)
! km        : kalman gain in updating x
real, dimension(ni,nj,nk,nv,nm), intent(inout) :: x,xo
real, dimension(ni,nj,nk,nv), intent(inout)    :: xm,xom
real, dimension(ni,nj,nk,nv,nm)                :: xf
real, dimension(ni,nj,nk,nv)                   :: x_t, std_x,std_xf, xmf
real, dimension(3,3,nk,nv,nm,nicpu*njcpu) :: xobsend, xob
real, dimension(3,3,nk,nv) :: tmp, tmpsend
real, dimension(obs%num,numbers_en+1) :: ya, yasend, yf,yf1
real, dimension(obs%num) :: obs_kk, obs_kk_send
integer, dimension(obs%num) :: kick_flag
real, allocatable, dimension(:,:,:)   :: km, kmsend, km1, km1send
real, allocatable, dimension(:,:,:,:) :: x1
! for satellite radiance
integer :: iob_radmin,iob_radmax
real, dimension(obs%num) :: yasend_tb, ym_radiance
real, dimension(ni,nj,nk)     :: xq_n,xq_p
real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend


read(wrf_file(6:10),'(i5)')iunit
if ( my_proc_id==0 ) then
   write(*, *)'   '
   write(*, *)'---------------------------------------------------'
   write(*, *)'.... Begin EnKF Assimilation ....'
   write(*, *)'   '
endif

! Get domain info
call get_variable2d( wrf_file, 'XLONG     ', ix, jx, 1, xlong )
call get_variable2d( wrf_file, 'XLAT      ', ix, jx, 1, xlat )
call get_variable2d( wrf_file, 'LANDMASK  ', ix, jx, 1, xland )
call get_variable2d( wrf_file, 'HGT       ', ix, jx, 1, hgt )
call get_variable2d( wrf_file, 'TSK       ', ix, jx, 1, tsk )
call get_variable2d( wrf_file, 'LU_INDEX  ', ix, jx, 1, lu_index )
call get_variable1d( wrf_file, 'ZNW       ', kx+1, 1, znw )
call get_variable1d( wrf_file, 'ZNU       ', kx, 1, znu )
call get_variable0d( wrf_file, 'P_TOP     ', 1, p_top )
call open_file(wrf_file, nf_nowrite, fid)
rcode = nf_get_att_int(fid, nf_global, 'GRID_ID', grid_id)
num_update_var = 0
do m = 1, 20
   if ( len_trim(updatevar(m))>=1 ) num_update_var=num_update_var+1
enddo

!subdomain start and end indices in full domain
istart=iid*ni+1
iend=(iid+1)*ni
jstart=jid*nj+1
jend=(jid+1)*nj
if(iid==(nicpu-1)) iend=ix+1
if(jid==(njcpu-1)) jend=jx+1
if(iend<istart .or. jend<jstart) then
  if(m==1.and.gid==0) write(*,'(a,i6,a)') '*******domain is too small to be decomposed! set nicpu,njcpu smaller.********'
  stop
endif

! I. calculate ya=Hxa,Hxm
if(my_proc_id==0) write(*,*) 'Calculating Hx...'
timer = MPI_Wtime()

!finding observation z location from xom
z_cycle: do iob=1,obs%num
  obstype=obs%type(iob)
  i1=int(obs%position(iob,1))
  j1=int(obs%position(iob,2))
  tmpsend=0.
  tmp=0.
  do j=j1,j1+2
  do i=i1,i1+2
     if(i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) then
        tmpsend(i-i1+1,j-j1+1,:,:)=xom(i-istart+1,j-jstart+1,:,:)
     endif
  enddo
  enddo
  call MPI_Allreduce(tmpsend,tmp,3*3*nk*nv,MPI_REAL,MPI_SUM,s_comm,ierr)
  write(filename,'(a5,i5.5)') wrf_file(1:5), iunit+numbers_en+1-1
  if ( obstype(1:5) == 'Radar' ) then
      call calc_rv_position_k(filename,proj,tmp,ix,jx,kx,nv,iob,xlong,znw,obs%position(iob,3))
  else if ( obstype(1:1) == 'P' .or. obstype(1:1) == 'H' ) then
      call calc_sounding_position_k(filename,proj,tmp,ix,jx,kx,nv,iob,xlong,znu,znw,p_top,obs%position(iob,3))
  else if ( obstype(1:8) == 'Radiance' ) then
      obs%position(iob,3) = 15.
  else
      obs%position(iob,3) = 1.
  endif
end do z_cycle

!print out obs if debugging
if(print_detail>4 .and. my_proc_id==0) then
   do iob=1,obs%num
     write(*,'(a,i6,3a,f10.2,a,3f10.2)') 'No.', iob,' ',obs%type(iob),'= ', obs%dat(iob), ' at (x,y,z)=', obs%position(iob,1:3)
   enddo
endif

! use cpus in the s_comm simultaneously to parallelize the calculation of Hx
! y=Hx stored in ya(obs_num, numbers_en+1)
ya=0.
yasend=0.
nob=nicpu*njcpu !number of obs in a batch, processed by cpus with different sid simotaneously (parallelized)
obs_cycle: do ig=1,int(obs%num/nob)+1
! prepare x profiles near obs position (xob) for xb calculation
   xobsend=0.
   xob=0.
   do ii=1,nob
     iob=(ig-1)*nob+ii
     if(iob<=obs%num) then
       i1=int(obs%position(iob,1))
       j1=int(obs%position(iob,2))
       do j=j1,j1+2
       do i=i1,i1+2
         if(i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) &
            xobsend(i-i1+1,j-j1+1,:,:,:,ii)=xo(i-istart+1,j-jstart+1,:,:,:)
       enddo
       enddo
     endif
   enddo
   call MPI_Allreduce(xobsend,xob,3*3*nk*nv*nm*nob,MPI_REAL,MPI_SUM,s_comm,ierr)
   iob=(ig-1)*nob+sid+1
   if(iob>obs%num) cycle obs_cycle
   obstype=obs%type(iob)
   do n = 1, nm
     ie=(n-1)*nmcpu+gid+1
     if (ie<=numbers_en+1) then
       write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+ie-1
       if ( obstype(1:5) == 'Radar' ) then
         call xb_to_rv(filename,proj,xob(:,:,:,:,n,sid+1),ix,jx,kx,nv,iob,xlong,znw,yasend(iob,ie),obs%position(iob,3),1)
       else if ( obstype == 'longitude ' .or. obstype == 'latitude  ' .or. obstype == 'min_slp   ' ) then
         call hurricane_center_assimilation(filename,ix,jx,kx,int(obs%position(iob,1)),int(obs%position(iob,2)), center_xb,znu,znw,xlong,xlat,proj)
         if ( obstype == 'longitude ' ) yasend(iob,ie) = center_xb(1)
         if ( obstype == 'latitude  ' ) yasend(iob,ie) = center_xb(2)
         if ( obstype == 'min_slp   ' ) yasend(iob,ie) = center_xb(3)*100.
       else if ( obstype(1:1) == 'P' .or. obstype(1:1) == 'H'  ) then
         call xb_to_sounding (filename,proj,xob(:,:,:,:,n,sid+1),ix,jx,kx,nv,iob,xlong,znu,znw,p_top,yasend(iob,ie))
       else if ( obstype(1:1) == 'S' ) then
         call xb_to_surface(filename,proj,xob(:,:,:,:,n,sid+1),ix,jx,kx,nv,iob,xlong,xland,lu_index,znu,znw,p_top,times,yasend(iob,ie))
       else if ( obstype(1:3) == 'slp' ) then
         call xb_to_slp(filename,xob(:,:,:,:,n,sid+1),ix,jx,kx,nv,iob,znu,znw,yasend(iob,ie))
       else if ( obstype(1:2) == 'pw' ) then
         call xb_to_pw(filename,xob(:,:,:,:,n,sid+1),ix,jx,kx,nv,iob,znu,znw,yasend(iob,ie))
       end if
     endif
   enddo
end do obs_cycle

!--ensemble loop for satellite radiance
 !---every grid is calculated in subroutine xb_to_radiance

if(raw%radiance%num.ne.0) then
  yasend_tb=0.
  do ie = 1, numbers_en
    yasend_tb = 0.0
    write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+ie-1
    !if ( my_proc_id == 0 ) write(*,*) "calculating radiance prior for member",ie
    call xb_to_radiance(filename,proj,ix,jx,kx,xlong,xlat,xland,iob_radmin,iob_radmax,yasend_tb,1)
    yasend(iob_radmin:iob_radmax,ie) = yasend_tb(iob_radmin:iob_radmax)
    yasend(iob_radmin:iob_radmax,numbers_en+1)=yasend(iob_radmin:iob_radmax,numbers_en+1)+yasend_tb(iob_radmin:iob_radmax)/real(numbers_en)  !calculate ya mean here
  enddo
endif
call MPI_Allreduce(yasend,ya,obs%num*(numbers_en+1),MPI_REAL,MPI_SUM,comm,ierr)

!print out obs/prior for debugging
if(print_detail>4 .and. my_proc_id==0) then
   do iob=1,obs%num
     write(*,'(a,i6,3a,f10.2,a,3f10.2)') 'No.', iob,' ',obs%type(iob),'=', obs%dat(iob), ' at (x,y,z)=', obs%position(iob,1:3)
   enddo
   write(format1,'(a,i4,a)') '(i5,a,',numbers_en+1,'f10.2)'
   do iob=1,obs%num
     if(print_detail>4 .and. my_proc_id==0) write(*,format1) iob, ' '//obs%type(iob)//' ya=',ya(iob,:)!numbers_en+1)
   enddo
endif

if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' Calculation of y=Hx tooks ', MPI_Wtime()-timer, ' seconds.'

!calculate ensemble perturbations
!apply prior inflation to ensemble spread
do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   if ( ie<=numbers_en+1 ) x(:,:,:,:,n)=inflate*(x(:,:,:,:,n)-xm)
enddo
xf=x   !save a copy of prior perturbation (for relaxation)
xmf=xm
yf=ya  !save a copy of observation prior (for error statistics)

!prior multiplicative inflation
! x=x*inflate

! II. assimilate obs
if(my_proc_id==0) write(*,*) 'Assimilating obs...'

timer=MPI_Wtime()
t1=0.; t2=0.; t3=0.; t4=0.;

assimilated_obs_num=0
obs_assimilate_cycle : do it = 1,obs%num
! 0. get basic information of the obs being assimilated
   iob=ind(it)
   kick_flag(iob)=0
   obstype = obs%type(iob)
   error = obs%err(iob)
   y_hxm = obs%dat(iob) - ya(iob,numbers_en+1)

   if ( my_proc_id==0 ) write(*,'(a,i6,a,f10.2,a,f10.2,a,f8.2,a,f8.2,a,i4,a,i4)') &
      'No.',iob,' '//obstype//' =',obs%dat(iob), ' ya=', ya(iob,numbers_en+1), ' y-ya=', y_hxm, &
      ' err=',error,' hroi=',obs%roi(iob,1),' vroi=',obs%roi(iob,2)

   if( abs(y_hxm)>(error*5.) .and. &
      .not.(obstype=='min_slp   ' .or. obstype=='longitude ' .or. obstype=='latitude  ' &
      .or. obstype=='slp       ' .or. obstype=='Radiance  ' &
      .or. obstype(1:5)=='Radar') ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for large error'
      kick_flag(iob)=1
      cycle obs_assimilate_cycle
   endif

   if(obs%position(iob,3)>kx .or. obs%position(iob,3)<1) then
      kick_flag(iob)=1
      cycle obs_assimilate_cycle
   endif

   if( any(yf(iob,:)==-888888.0) ) then
      if ( my_proc_id==0 ) write(*,*)' ...kicked off for invalid value'
      kick_flag(iob)=1
      cycle obs_assimilate_cycle
   endif

   assimilated_obs_num=assimilated_obs_num+1
   ngx = max(obs%roi(iob,1),max(nicpu,njcpu)/2+1)
   ngz = obs%roi(iob,2)

! fac*var=hBh'=<hxa hxa'>=<(ya-yam)(ya-yam)'>
! alpha is for Square-Root filter (Whitaker and Hamill 2002, MWR130,1913-1924)
! d=hBh'+r, for each obs, d reduces to a scalar, r=error^2.
   var = 0.
   do ie = 1, numbers_en
      hxa(ie) = ya(iob,ie)-ya(iob,numbers_en+1)
      var     = var + hxa(ie)*hxa(ie)
   enddo
   fac  = 1./real(numbers_en-1)
   d    = fac * var + error * error
   alpha = 1.0/(1.0+sqrt(error*error/d))

! cycle through variables to process, in x, 2D variables are stored in 3D form (with values only on
! k=1), when sending them among cpus, only the lowest layer (:,:,1) are sent and received.
update_x_var : do m=1,nv
t0=MPI_Wtime()
   varname=enkfvar(m)
   update_flag = 0
   do iv = 1, num_update_var
     if ( varname .eq. updatevar(iv) ) update_flag = 1
   enddo
!
!!---Observation Error Inflation & Successive Covariance Localization (Zhang et al. 2009)
!!      for Radiance assimilation by Minamide 2015.3.14
   d    = fac * var + error * error
   alpha = 1.0/(1.0+sqrt(error*error/d))
   if (obstype=='Radiance  ') then
     d = max(fac * var + error * error, y_hxm * y_hxm)
     alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
     !obs%err(iob) = sqrt(d-fac * var)
     if ( my_proc_id == 0 .and. sqrt(d-fac * var) > error .and. varname=='T         ')&
          write(*,*) 'observation-error inflated to ',sqrt(d-fac * var)
   endif
!!---OEI & SCL end

   if (obstype=='Radiance  ') then
     if (varname=='P         ' .or. varname=='PH        ' .or. varname=='MU        ' .or. varname=='PSFC      ') then
       update_flag=0
     else
       update_flag=1
     end if
   end if
   !if (obstype=='Radiance  ') then
   !  if (varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. &
   !    varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
   !    update_flag=1
   !  else
   !    update_flag=0
   !  end if
   !end if

   if ( update_flag==0 ) cycle update_x_var

! start and end indices of the update zone of the obs
     ist = max( update_is, int(obs%position(iob,1))-ngx )
     ied = min( update_ie, int(obs%position(iob,1))+ngx )
     jst = max( update_js, int(obs%position(iob,2))-ngx )
     jed = min( update_je, int(obs%position(iob,2))+ngx )
     kst = max( update_ks, int(obs%position(iob,3))-ngz )
     ked = min( update_ke, int(obs%position(iob,3))+ngz )
     call wrf_var_dimension ( wrf_file, varname, ix, jx, kx, ii, jj, kk )
     if ( kk == 1 ) then
       kst  = 1
       ked = 1
     endif
! indices of the update zone slab handled by cpu with sid.
     uist=iid*int((ied-ist+1)/nicpu)+ist
     uied=(iid+1)*int((ied-ist+1)/nicpu)+ist-1
     ujst=jid*int((jed-jst+1)/njcpu)+jst
     ujed=(jid+1)*int((jed-jst+1)/njcpu)+jst-1
     if(iid==nicpu-1) uied=ied
     if(jid==njcpu-1) ujed=jed
     if(uied<uist .or. ujed<ujst) then
       if(m==1.and.gid==0) then
          write(*,'(a,i6,a)') '*******update zone of obs #',iob,' is too small to be decomposed.********'
          write(*,*) 'update zone', uist,uied,ujst,ujed,kst,ked
          write(*,*) 'obs location', obs%position(iob,:)
          stop
       endif
     endif
     allocate(x1(uied-uist+1, ujed-ujst+1, ked-kst+1, nm))
     x1=0.

! 1. send and recv x from all sid to the slab a cpu need for later use.
!    if data needed from/to sid, we determine the indices of data sent/received by the
!    sid of dest/source cpus, and perform the send/recv
!
     do is=0,nicpu*njcpu-1
       if(sid>is) then
         ii=mod(is,nicpu)
         jj=int(is/nicpu)
         sistart=ii*ni+1
         siend=(ii+1)*ni
         sjstart=jj*nj+1
         sjend=(jj+1)*nj
         if(ii==nicpu-1) siend=ix+1
         if(jj==njcpu-1) sjend=jx+1
         if(uist<=siend .and. uied>=sistart .and. ujst<=sjend .and. ujed>=sjstart) then
            call MPI_RECV(x1(max(uist,sistart)-uist+1:min(uied,siend)-uist+1, &
                         max(ujst,sjstart)-ujst+1:min(ujed,sjend)-ujst+1, 1:ked-kst+1,:), &
                        (min(uied,siend)-max(uist,sistart)+1)*(min(ujed,sjend)-max(ujst,sjstart)+1)*(ked-kst+1)*nm, &
                         MPI_REAL, is, 0, s_comm, status, ierr )
         endif
       endif
     enddo
     do i=1,nicpu*njcpu
       is=mod(sid+i,nicpu*njcpu)
       ii=mod(is,nicpu)
       jj=int(is/nicpu)
       sist=ii*int((ied-ist+1)/nicpu)+ist
       sied=(ii+1)*int((ied-ist+1)/nicpu)+ist-1
       sjst=jj*int((jed-jst+1)/njcpu)+jst
       sjed=(jj+1)*int((jed-jst+1)/njcpu)+jst-1
       if(ii==nicpu-1) sied=ied
       if(jj==njcpu-1) sjed=jed
       if(sist<=iend .and. sied>=istart .and. sjst<=jend .and. sjed>=jstart) then
         if(sid==is) then
           x1(max(sist,istart)-sist+1:min(sied,iend)-sist+1, &
              max(sjst,jstart)-sjst+1:min(sjed,jend)-sjst+1, 1:ked-kst+1,:) = &
           x( max(sist,istart)-istart+1:min(sied,iend)-istart+1, &
              max(sjst,jstart)-jstart+1:min(sjed,jend)-jstart+1, kst:ked,m,:)
         else
           call MPI_SEND(x(max(sist,istart)-istart+1:min(sied,iend)-istart+1, &
                         max(sjst,jstart)-jstart+1:min(sjed,jend)-jstart+1,kst:ked,m,:), &
                         (min(sied,iend)-max(sist,istart)+1)*(min(sjed,jend)-max(sjst,jstart)+1)*(ked-kst+1)*nm, &
                         MPI_REAL, is, 0, s_comm, ierr)
         endif
       endif
     enddo
     do is=0,nicpu*njcpu-1
       if(sid<is) then
         ii=mod(is,nicpu)
         jj=int(is/nicpu)
         sistart=ii*ni+1
         siend=(ii+1)*ni
         sjstart=jj*nj+1
         sjend=(jj+1)*nj
         if(ii==nicpu-1) siend=ix+1
         if(jj==njcpu-1) sjend=jx+1
         if(uist<=siend .and. uied>=sistart .and. ujst<=sjend .and. ujed>=sjstart) then
            call MPI_RECV(x1(max(uist,sistart)-uist+1:min(uied,siend)-uist+1, &
                             max(ujst,sjstart)-ujst+1:min(ujed,sjend)-ujst+1, 1:ked-kst+1,:), &
                          (min(uied,siend)-max(uist,sistart)+1)*(min(ujed,sjend)-max(ujst,sjstart)+1)*(ked-kst+1)*nm, &
                          MPI_REAL, is, 0, s_comm, status, ierr )
         endif
       endif
     enddo
t1=t1+MPI_Wtime()-t0
t0=MPI_Wtime()

! 2. calculate Bh' part of kalman gain, Bh'=<xa xa'>h'=SUM(xa*hxa)*fac
!    the summation among ensemble members is done by allreduce within the g_comm group (members are
!    distributed among cpus with different gid.
     allocate ( km     (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     allocate ( kmsend (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     km=0.
     kmsend=0.
     do n=1,nm
       ie=(n-1)*nmcpu+gid+1
       if(ie<=numbers_en) then
         kmsend=kmsend+x1(:,:,:,n)*hxa(ie)
       endif
     enddo
     call MPI_Allreduce(kmsend,km,(uied-uist+1)*(ujed-ujst+1)*(ked-kst+1),MPI_REAL,MPI_SUM,g_comm,ierr)
     km = km * fac / d
t2=t2+MPI_Wtime()-t0
t0=MPI_Wtime()

! 3. calculate corr_coef=localization factor, 1 at obs position, reduce to 0 at roi
!    2-D variables are treated as at surface k=1
     do k = kst,ked
     do j = ujst,ujed
     do i = uist,uied
       call corr(real(i-obs%position(iob,1)),real(j-obs%position(iob,2)),real(k-obs%position(iob,3)),obs%roi(iob,1),obs%roi(iob,2),corr_coef)
       !if ( obstype == 'longitude ' .or. obstype == 'latitude  ' ) corr_coef = 1.0
       km(i-uist+1,j-ujst+1,k-kst+1) = km(i-uist+1,j-ujst+1,k-kst+1) * corr_coef
     enddo
     enddo
     enddo
t3=t3+MPI_Wtime()-t0
t0=MPI_Wtime()

! 4. send kalman gain to sid slabs of domain, to be finally updated to x
!    the organization of sending/receiving is similar as step 1.
     if( ist<=iend .and. ied>=istart .and. jst<=jend .and. jed>=jstart ) then
       allocate(km1(min(ied,iend)-max(ist,istart)+1, min(jed,jend)-max(jst,jstart)+1, ked-kst+1))
       km1=0.
     endif
     if( ist<=iend .and. ied>=istart .and. jst<=jend .and. jed>=jstart ) then
       sis=max(ist,istart)
       sie=min(ied,iend)
       sjs=max(jst,jstart)
       sje=min(jed,jend)
       do is=0,nicpu*njcpu-1
         if(sid>is) then
           ii=mod(is,nicpu)
           jj=int(is/nicpu)
           sist=ii*int((ied-ist+1)/nicpu)+ist
           sied=(ii+1)*int((ied-ist+1)/nicpu)+ist-1
           sjst=jj*int((jed-jst+1)/njcpu)+jst
           sjed=(jj+1)*int((jed-jst+1)/njcpu)+jst-1
           if(ii==nicpu-1) sied=ied
           if(jj==njcpu-1) sjed=jed
           if(sist<=sie .and. sied>=sis .and. sjst<=sje .and. sjed>=sjs) then
             call MPI_RECV(km1(max(sist,sis)-sis+1:min(sied,sie)-sis+1, &
                             max(sjst,sjs)-sjs+1:min(sjed,sje)-sjs+1, 1:ked-kst+1), &
                          (min(sied,sie)-max(sist,sis)+1)*(min(sjed,sje)-max(sjst,sjs)+1)*(ked-kst+1), &
                          MPI_REAl, is, 0, s_comm, status, ierr)
           endif
         endif
       enddo
     endif
     do i=1,nicpu*njcpu
       is=mod(sid+i,nicpu*njcpu)
       ii=mod(is,nicpu)
       jj=int(is/nicpu)
       sistart=ii*ni+1
       siend=(ii+1)*ni
       sjstart=jj*nj+1
       sjend=(jj+1)*nj
       if(ii==nicpu-1) siend=ix+1
       if(jj==njcpu-1) sjend=jx+1
       if( ist<=siend .and. ied>=sistart .and. jst<=sjend .and. jed>=sjstart ) then
         sis=max(ist,sistart)
         sie=min(ied,siend)
         sjs=max(jst,sjstart)
         sje=min(jed,sjend)
         if(uist<=sie .and. uied>=sis .and. ujst<=sje .and. ujed>=sjs) then
           if(sid==is) then
             km1(max(uist,sis)-sis+1:min(uied,sie)-sis+1, &
                 max(ujst,sjs)-sjs+1:min(ujed,sje)-sjs+1, 1:ked-kst+1) = &
             km (max(uist,sis)-uist+1:min(uied,sie)-uist+1, &
                 max(ujst,sjs)-ujst+1:min(ujed,sje)-ujst+1, 1:ked-kst+1)
           else
             call MPI_SEND(km(max(uist,sis)-uist+1:min(uied,sie)-uist+1, &
                            max(ujst,sjs)-ujst+1:min(ujed,sje)-ujst+1, 1:ked-kst+1), &
                           (min(uied,sie)-max(uist,sis)+1)*(min(ujed,sje)-max(ujst,sjs)+1)*(ked-kst+1), &
                           MPI_REAl, is, 0, s_comm, ierr)
           endif
         endif
       endif
     enddo
     if( ist<=iend .and. ied>=istart .and. jst<=jend .and. jed>=jstart ) then
       sis=max(ist,istart)
       sie=min(ied,iend)
       sjs=max(jst,jstart)
       sje=min(jed,jend)
       do is=0,nicpu*njcpu-1
         if(sid<is) then
           ii=mod(is,nicpu)
           jj=int(is/nicpu)
           sist=ii*int((ied-ist+1)/nicpu)+ist
           sied=(ii+1)*int((ied-ist+1)/nicpu)+ist-1
           sjst=jj*int((jed-jst+1)/njcpu)+jst
           sjed=(jj+1)*int((jed-jst+1)/njcpu)+jst-1
           if(ii==nicpu-1) sied=ied
           if(jj==njcpu-1) sjed=jed
           if(sist<=sie .and. sied>=sis .and. sjst<=sje .and. sjed>=sjs) then
             call MPI_RECV(km1(max(sist,sis)-sis+1:min(sied,sie)-sis+1, &
                             max(sjst,sjs)-sjs+1:min(sjed,sje)-sjs+1, 1:ked-kst+1), &
                            (min(sied,sie)-max(sist,sis)+1)*(min(sjed,sje)-max(sjst,sjs)+1)*(ked-kst+1), &
                            MPI_REAl, is, 0, s_comm, status, ierr)
           endif
         endif
       enddo
     endif
t4=t4+MPI_Wtime()-t0
t0=MPI_Wtime()

! 5. Update x and xm
!    for xa, update with alpha*corr_coef*Bh'*(0-hxa)/d, (Square-Root filter)
!    for xm, update with corr_coef*Bh'*(y-hxm)/d
     if( ist<=iend .and. ied>=istart .and. jst<=jend .and. jed>=jstart ) then
       do n=1,nm
         ie=(n-1)*nmcpu+gid+1
         if(ie<=numbers_en) then
             x (max(ist,istart)-istart+1:min(ied,iend)-istart+1, &
                max(jst,jstart)-jstart+1:min(jed,jend)-jstart+1, kst:ked,m,n) = &
             x (max(ist,istart)-istart+1:min(ied,iend)-istart+1, &
                max(jst,jstart)-jstart+1:min(jed,jend)-jstart+1, kst:ked,m,n) - km1 * alpha * hxa(ie)
         endif
       enddo
       xm (max(ist,istart)-istart+1:min(ied,iend)-istart+1, &
           max(jst,jstart)-jstart+1:min(jed,jend)-jstart+1, kst:ked,m) = &
       xm (max(ist,istart)-istart+1:min(ied,iend)-istart+1, &
           max(jst,jstart)-jstart+1:min(jed,jend)-jstart+1, kst:ked,m) + km1 * y_hxm

       ! remove negative Q values
       !if(varname.eq.'QVAPOR    ') then
       !  do n=1,nm
       !    ie=(n-1)*nmcpu+gid+1
       !    if(ie<=numbers_en) then
       !      where( (xm(ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m)+&
       !              x (ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m,n))<0. ) &
       !            x (ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m,n) = &
       !         0.-xm(ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m)
       !    endif
       !  enddo
       !  where( xm(ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m)<0.) &
       !         xm(ist-istart+1:ied-istart+1,jst-jstart+1:jed-jstart+1,kst:ked,m)=0.
       !endif
       deallocate(km1)
     endif
     deallocate(km)
     deallocate(kmsend)
     deallocate(x1)
t5=t5+MPI_Wtime()-t0
t0=MPI_Wtime()

   enddo update_x_var

! 6. update relevant ya
!    for members ya(:,1:n)=ya+ym (perturbation+mean), and ym=ya(:,n+1)
!    h_i B h_j'= <ya_i ya_j'> = cov(ya_i,ya_j)
!    d   = variance(ya) + obserr variance
!    ym=ym+corr_coef*hBh'(y-ym)/d
!    ya=ya+alpha*corr_coef*hBh'(0-ya)/d
!    --- basically these are the update equations of x left-multiplied by H.
   ngx = max(obs%roi(iob,1),max(nicpu,njcpu)/2+1)
   ist = max( update_is, int(obs%position(iob,1))-ngx )
   ied = min( update_ie, int(obs%position(iob,1))+ngx )
   jst = max( update_js, int(obs%position(iob,2))-ngx )
   jed = min( update_je, int(obs%position(iob,2))+ngx )
   kst = max( update_ks, int(obs%position(iob,3))-ngz )
   ked = min( update_ke, int(obs%position(iob,3))+ngz )
   update_y_cycle : do iiob=1,obs%num
! skip those ya outside update zone
      if ( obs%position(iiob,1)<ist .or. obs%position(iiob,1)>ied .or. &
           obs%position(iiob,2)<jst .or. obs%position(iiob,2)>jed .or. &
           obs%position(iiob,3)<kst .or. obs%position(iiob,3)>ked ) cycle update_y_cycle

      call corr(real(obs%position(iiob,1)-obs%position(iob,1)), real(obs%position(iiob,2)-obs%position(iob,2)), &
                real(obs%position(iiob,3)-obs%position(iob,3)), obs%roi(iob,1), obs%roi(iob,2), corr_coef)
      obstype = obs%type(iiob)
      var=0.
      cov=0.
      do ie=1,numbers_en
         hxa(ie) = ya(iob,ie)-ya(iob,numbers_en+1)
         var     = var + hxa(ie)*hxa(ie)
         cov     = cov + hxa(ie)*(ya(iiob,ie)-ya(iiob,numbers_en+1))
      enddo
      error=obs%err(iob)
      fac  = 1./real(numbers_en-1)
      d    = fac * var + error * error
      alpha= 1.0/(1.0+sqrt(error*error/d))
!!! relaxation by quality contoling
!   if ((obstype=='Radiance  ') .and. (abs(y_hxm)>max(error*3.,sqrt(fac *var))))then
!! error = oma_omb method
!    d_ogn = d
!     d    = fac * var + max(abs(y_hxm-fac*var*y_hxm/d)*abs(y_hxm),error**2)
!     alpha =1.0/(1.0+sqrt(max(abs(y_hxm-fac*var*y_hxm/d_ogn)*abs(y_hxm),error**2)/d))
!!error = y_hxm
!     d    = fac * var + abs(y_hxm) * abs(y_hxm)
!     alpha = 1.0/(1.0+sqrt(abs(y_hxm)*abs(y_hxm)/d))
!!

!!!---- OEI
   if (obstype=='Radiance  ') then
      d = max(fac * var + error * error, y_hxm * y_hxm)
      alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
   endif
!!! relaxation end


      do ie=1,numbers_en+1
         if(ie<=numbers_en) &
            ya(iiob,ie)=ya(iiob,ie)-corr_coef*alpha*fac*cov*(ya(iob,ie)-ya(iob,numbers_en+1))/d !perturbation
         ya(iiob,ie)=ya(iiob,ie)+corr_coef*fac*cov*(obs%dat(iob)-ya(iob,numbers_en+1))/d        !mean
         !if(obstype(10:10)=='Q' .and. ya(iiob,ie)<0.) ya(iiob,ie)=0.  ! remove negative values of Q
      enddo
   end do update_y_cycle

end do obs_assimilate_cycle

yf1=ya !save another copy of ya

if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 1 tooks ', t1, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 2 tooks ', t2, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 3 tooks ', t3, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 4 tooks ', t4, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' 5 tooks ', t5, ' seconds.'
if ( my_proc_id == 0 ) write(*,'(a,f7.2,a)')' Assimilation tooks ', MPI_Wtime()-timer, ' seconds.'

!---------------------------------------------------------------------------------
! IV. Diagnostics for filter performance.
! 1. assimilated obs information output to fort.10000.
m_d2=0.0
m_var_b=0.0
m_var_a=0.0
n=0
do iob=1,obs%num
  call ij_to_latlon(proj, obs%position(iob,1), obs%position(iob,2), center_xb(1), center_xb(2))
  var_a=0.0
  var_b=0.0
  do ie=1,numbers_en
    var_a=var_a+(ya(iob,ie)-ya(iob,numbers_en+1))**2
    var_b=var_b+(yf(iob,ie)-yf(iob,numbers_en+1))**2
  enddo
  var_a=var_a/real(numbers_en-1)
  var_b=var_b/real(numbers_en-1)
  if(kick_flag(iob)==0) then
    ounit=10000
  else
    ounit=10001
  endif
  !write diagnostics to fort.1000?
  if ( my_proc_id == 0 ) then
    write(ounit,'(a,6f11.2,2i4,6f11.2)') obs%type(iob), &   !type of observation
      center_xb(1), center_xb(2), &     !latitude, longitude
      obs%position(iob,4), &            !height
      obs%position(iob,1:3), &          !grid location: i,j,k
      obs%roi(iob,1:2), &               !horizontal ROI (# of grid points), vertical ROI (# of layers)
      obs%dat(iob), &                   !observation (y^o)
      obs%err(iob), &                   !observation error variance (R)
      yf(iob,numbers_en+1), &           !observation prior mean (H\bar{x}^f)
      ya(iob,numbers_en+1), &           !observation posterior mean (H\bar{x}^a)
      sqrt(var_b), &                    !observation prior spread (sqrt{H P^f H^T})
      sqrt(var_a)                       !observation posterior spread (sqrt{H P^f H^T}) (before relaxation)
  endif
  !innovation statistics (used in adaptive relaxation)
  !if(my_proc_id==0) write(*,*) iob,'=',kick_flag(iob)
  if(kick_flag(iob)==0) then
    n=n+1
    m_d2=m_d2+((obs%dat(iob)-yf(iob,numbers_en+1))/obs%err(iob))**2
    m_var_b=m_var_b+var_b/(obs%err(iob)**2)
    m_var_a=m_var_a+var_a/(obs%err(iob)**2)
  endif
end do
m_d2=m_d2/real(n)
m_var_b=m_var_b/real(n)
m_var_a=m_var_a/real(n)

!---------------------------------------------------------------------------------
! V. Relaxation
if ( my_proc_id==0 ) write(*,*)'Performing covariance relaxation...'
!------inflation factor from innovation statistics
if(my_proc_id==0) write(*,*) 'lambda=',sqrt((m_d2-1)/m_var_b), ' kappa=',(sqrt(m_var_b)-sqrt(m_var_a))/sqrt(m_var_a)
!!!adaptively determine mixing coef
if(relax_adaptive) then
  la=max(sqrt((m_d2-1.0)/m_var_b),1.0)
  ka=(sqrt(m_var_b)-sqrt(m_var_a))/sqrt(m_var_a)
  mixing=(la-1.0)/ka
  if(my_proc_id==0) write(*,*) 'd2=',m_d2
  if(my_proc_id==0) write(*,*) 'lambda=',la,' kappa=',ka,' mixing=',mixing
end if

!---option1. relax to prior perturbation (Zhang et al. 2004)
if(relax_opt==0) then
  do n = 1, nm
    ie = (n-1)*nmcpu+gid+1
    if( ie<=numbers_en ) &
      x(:,:,:,:,n)=(1.0-mixing)*x(:,:,:,:,n) + mixing*xf(:,:,:,:,n)
  enddo
end if

!---option2. relax to prior spread (Whitaker and Hamill 2012), adaptive version (ACR; Ying and Zhang 2015)
if(relax_opt==1) then
  !------ensemble spread in state space
  std_x =0.0
  std_xf=0.0
  call MPI_Allreduce(sum( x**2,5), std_x,ni*nj*nk*nv,MPI_REAL,MPI_SUM,g_comm,ierr)
  call MPI_Allreduce(sum(xf**2,5),std_xf,ni*nj*nk*nv,MPI_REAL,MPI_SUM,g_comm,ierr)
  std_x =sqrt( std_x/real(numbers_en-1))
  std_xf=sqrt(std_xf/real(numbers_en-1))
  do n=1,nm
    ie=(n-1)*nmcpu+gid+1
    if(ie==numbers_en+1) &
      x(:,:,:,:,n)=x(:,:,:,:,n)*(mixing*(std_xf-std_x)/std_x+1)
  enddo
end if

 !!!! removing analysis increment near model top - Yue Ying 2016
kzdamp=10 !!!number of damping layers from model top: check zdamp
if ( my_proc_id==0 ) write(*,*) 'removing analysis increment near model top'
do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   do k=0,kzdamp
      beta=(cos((real(k)/real(kzdamp))*4*atan(1.0)/2))**2  !beta=1 at top
      if( ie<=numbers_en ) then
        x(:,:,nk-k,:,n) = (1-beta)*x(:,:,nk-k,:,n) + beta*xf(:,:,nk-k,:,n)
      end if
      xm(:,:,nk-k,:) = (1-beta)*xm(:,:,nk-k,:) + beta*xmf(:,:,nk-k,:)
   end do
enddo


!Add the mean back to the analysis field
if ( my_proc_id==0 ) write(*,*)'Add the mean back to the analysis field...'
do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   if( ie<=numbers_en+1 )  &
      x(:,:,:,:,n) = x(:,:,:,:,n) + xm
enddo

!!! Removing negative Q-value by Minamide 2015.5.26
!if ( my_proc_id==0 ) write(*,*) 'updating negative values'
!if(raw%radiance%num.ne.0) then
!  do m=1,nv
!    varname=enkfvar(m)
!    xq_p = 0.
!    xq_n = 0.
!    xq_psend = 0.
!    xq_nsend = 0.
!    if (varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. &
!        varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
!      do n=1,nm
!        ie=(n-1)*nmcpu+gid+1
!        if(ie==numbers_en+1) write(*,*)'original xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
!        if(ie<=numbers_en) then
!          where(x(:,:,:,m,n) >= 0.) xq_psend(:,:,:,n) = x(:,:,:,m,n)
!          where(x(:,:,:,m,n) < 0.) xq_nsend(:,:,:,n) = x(:,:,:,m,n)
!        endif
!      enddo
!      call MPI_Allreduce(sum(xq_psend,4),xq_p,ni*nj*nk,MPI_REAL,MPI_SUM,comm,ierr)
!      call MPI_Allreduce(sum(xq_nsend,4),xq_n,ni*nj*nk,MPI_REAL,MPI_SUM,comm,ierr)
!      if ( my_proc_id==0 ) write(*,*) 'xq_p',minval(xq_p),'~',maxval(xq_p)
!      if ( my_proc_id==0 ) write(*,*) 'xq_n',minval(xq_n),'~',maxval(xq_n)
!      do n=1,nm
!        ie=(n-1)*nmcpu+gid+1
!        if(ie<=numbers_en) then
!          where(x(:,:,:,m,n) < 0.) x(:,:,:,m,n) = 0.
!          where(xq_p >= abs(xq_n).and.xq_p > 0.) x(:,:,:,m,n) = x(:,:,:,m,n)*(xq_p+xq_n)/xq_p
!          where(xq_p < abs(xq_n).or. xq_p == 0.) x(:,:,:,m,n) = 0.
!        endif
!        if(ie<=numbers_en+1) where(xm(:,:,:,m) < 0.) x(:,:,:,m,n) = 0.
!        if(ie==numbers_en+1) write(*,*)'non-negative xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
!      enddo
!    endif
!  enddo
!endif

end subroutine enkf
