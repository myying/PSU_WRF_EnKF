subroutine update_inflation(wrf_file,ix,jx,kx,ni,nj,nk,nv,nm,g_comm,s_comm,gid,sid,iid,jid,xom,xo,xm,x,inf_mean,inf_sd,ind,proj,times)
use constants
use namelist_define
use obs_define
use mpi_module
use netcdf
use radar
implicit none
integer, intent(in) :: sid,gid,iid,jid
integer, intent(in) :: ix,jx,kx,ni,nj,nk,nv,nm, g_comm,s_comm
integer, dimension(obs%num), intent(in) :: ind
character (len=10), intent(in) :: wrf_file
type (proj_info)    :: proj
character (len=10)  :: obstype
character (len=80), intent(in) :: times
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
real, dimension(ni,nj,nk,nv,nm), intent(inout) :: x,xo
real, dimension(ni,nj,nk,nv), intent(inout)    :: xm,xom,inf_mean,inf_sd
real, dimension(ni,nj,nk,nv,nm)                :: xf
real, dimension(ni,nj,nk,nv)                   :: x_t, std_x,std_xf, xmf, gamma_corr
real, dimension(3,3,nk,nv,nm,nicpu*njcpu) :: xobsend, xob
real, dimension(3,3,nk,nv) :: tmp, tmpsend
real, dimension(obs%num,numbers_en+1) :: ya, yasend, yf
real, dimension(obs%num) :: obs_kk, obs_kk_send
integer, dimension(obs%num) :: kick_flag
real, allocatable, dimension(:,:,:)   :: m1,m1send,m2,m2send,km,km1,km1send
real, allocatable, dimension(:,:,:,:) :: x1
integer :: iob_radmin,iob_radmax
real, dimension(obs%num) :: yasend_tb, ym_radiance
real, dimension(ni,nj,nk)     :: xq_n,xq_p
real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend
real :: rate,dist_2,sigma_p_2,sigma_o_2,inf_mean_old


read(wrf_file(6:10),'(i5)')iunit

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
do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   if ( ie<=numbers_en+1 ) x(:,:,:,:,n)=x(:,:,:,:,n)-xm
enddo
xf=x   !save a copy of prior perturbation (for relaxation)
yf=ya  !save a copy of observation prior (for error statistics)

! II. assimilate obs
if(my_proc_id==0) write(*,*) 'Assimilating obs...'

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
   varname=enkfvar(m)
   update_flag = 0
   do iv = 1, num_update_var
     if ( varname .eq. updatevar(iv) ) update_flag = 1
   enddo

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

! 2. calculate correlation
!    the summation among ensemble members is done by allreduce within the g_comm group (members are
!    distributed among cpus with different gid.
     allocate ( m1     (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     allocate ( m1send (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     allocate ( m2     (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     allocate ( m2send (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     allocate ( km     (uied-uist+1, ujed-ujst+1, ked-kst+1) )
     m1=0.
     m2=0.
     m1send=0.
     m2send=0.
     km=0.
     do n=1,nm
       ie=(n-1)*nmcpu+gid+1
       if(ie<=numbers_en) then
         m1send=m1send+x1(:,:,:,n)*hxa(ie)
         m2send=m2send+x1(:,:,:,n)**2
       endif
     enddo
     call MPI_Allreduce(m1send,m1,(uied-uist+1)*(ujed-ujst+1)*(ked-kst+1),MPI_REAL,MPI_SUM,g_comm,ierr)
     call MPI_Allreduce(m2send,m2,(uied-uist+1)*(ujed-ujst+1)*(ked-kst+1),MPI_REAL,MPI_SUM,g_comm,ierr)
     if(var<=0.) then
       km=0.
     else
       km=abs(m1/sqrt(m2*var))
     endif

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

!! 5. Update inf_mean
     gamma_corr=0.
     if( ist<=iend .and. ied>=istart .and. jst<=jend .and. jed>=jstart ) then
       gamma_corr(max(ist,istart)-istart+1:min(ied,iend)-istart+1, &
                  max(jst,jstart)-jstart+1:min(jed,jend)-jstart+1, kst:ked, m) = km1
       do i=max(ist,istart)-istart+1, min(ied,iend)-istart+1
       do j=max(jst,jstart)-jstart+1, min(jed,jend)-jstart+1
       do k=kst, ked
         if(gamma_corr(i,j,k,m)>0.) then
           dist_2=y_hxm**2
           sigma_p_2=var*fac
           sigma_o_2=1 !!error*error
           inf_mean_old=inf_mean(i,j,k,m)
           call change_GA_IG(inf_mean(i,j,k,m),inf_sd(i,j,k,m),rate)
           call linear_bayes(dist_2,sigma_p_2,sigma_o_2,inf_mean_old, &
                             gamma_corr(i,j,k,m),numbers_en,rate,inf_mean(i,j,k,m))
         endif
       enddo
       enddo
       enddo
       deallocate(km1)
     endif
     deallocate(km)
     deallocate(m1,m1send,m2,m2send)
     deallocate(x1)

   enddo update_x_var
end do obs_assimilate_cycle

!!!inflate
if ( my_proc_id==0 ) write(*,*)'Add the mean back to the analysis field...'
do n = 1, nm
   ie = (n-1)*nmcpu+gid+1
   if( ie<=numbers_en+1 )  &
      x(:,:,:,:,n) = inf_mean*x(:,:,:,:,n) + xm
enddo

end subroutine update_inflation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linear_bayes(dist_2, sigma_p_2, sigma_o_2, &
           lambda_mean, gamma_corr, ens_size, beta, new_cov_inflate)

real, intent(in)    :: dist_2, sigma_p_2, sigma_o_2, lambda_mean
real, intent(in)    :: gamma_corr, beta
real, intent(inout) :: new_cov_inflate
real :: PI=3.141592654

integer  :: ens_size
real :: theta_bar_2, like_bar, like_prime, theta_bar
real :: a, b, c, plus_root, minus_root, deriv_theta
real :: fac1, fac2, like_ratio

! Scaling factors
fac1 = (1.0 + gamma_corr * (sqrt(lambda_mean) - 1.0))**2
fac2 = -1.0 / ens_size

! Compute value of theta at current lambda_mean
if ( fac1 < abs(fac2) ) fac2 = 0.0
theta_bar_2 = (fac1+fac2) * sigma_p_2 + sigma_o_2
theta_bar   = sqrt(theta_bar_2)

! Compute constant coefficient for likelihood at lambda_bar
like_bar = exp(- 0.5 * dist_2 / theta_bar_2) / (sqrt(2.0 * PI) * theta_bar)

! If like_bar goes to 0, can't do anything, so just keep current values
! Density at current inflation value must be positive
if(like_bar <= 0.0) then
   new_cov_inflate = lambda_mean
   return
endif

! Next compute derivative of likelihood at this point
deriv_theta = 0.5 * sigma_p_2 * gamma_corr * ( 1.0 - gamma_corr + &
              gamma_corr * sqrt(lambda_mean) ) / ( theta_bar * sqrt(lambda_mean) )
like_prime  = like_bar * deriv_theta * (dist_2 / theta_bar_2 - 1.0) / theta_bar

! If like_prime goes to 0, can't do anything, so just keep current values
! We're dividing by the derivative in the quadratic equation, so this
! term better non-zero!
if(like_prime == 0. .OR. abs(like_bar) <= 0. .OR. abs(like_prime) <= 0. ) then
   new_cov_inflate = lambda_mean
   return
endif
like_ratio = like_bar / like_prime

a = 1.0 - lambda_mean / beta
b = like_ratio - 2.0 * lambda_mean
c = lambda_mean**2 - like_ratio * lambda_mean

! Use nice scaled quadratic solver to avoid precision issues
call solve_quadratic(a, b, c, plus_root, minus_root)

! Do a check to pick closest root
if(abs(minus_root - lambda_mean) < abs(plus_root - lambda_mean)) then
   new_cov_inflate = minus_root
else
   new_cov_inflate = plus_root
endif

! Do a final check on the sign of the updated factor
! Sometimes the factor can be very small (almost zero)
! From the selection process above it can be negative
! if the positive root is far away from it.
! As such, keep the current factor value
if(new_cov_inflate <= 0.0 .OR. new_cov_inflate /= new_cov_inflate) new_cov_inflate = lambda_mean

end subroutine linear_bayes


subroutine solve_quadratic(a, b, c, r1, r2)
real, intent(in)  :: a, b, c
real, intent(out) :: r1, r2
real :: scaling, as, bs, cs, disc

! Scale the coefficients to get better round-off tolerance
scaling = max(abs(a), abs(b), abs(c))
as = a / scaling
bs = b / scaling
cs = c / scaling

! Get discriminant of scaled equation
disc = sqrt(bs**2 - 4.0 * as * cs)

if(bs > 0.0) then
   r1 = (-bs - disc) / (2 * as)
else
   r1 = (-bs + disc) / (2 * as)
endif

! Compute the second root given the larger one
r2 = (cs / as) / r1

end subroutine solve_quadratic


!!----------------------------------------------
!!> Routine to change the Gaussian prior into an inverse gamma (IG).
!!> The Gaussian prior is represented by a mode (:= mean) and a variance; var
subroutine change_GA_IG(mode, var, beta)

real, intent(in)  :: mode, var
real, intent(out) :: beta

integer :: i
real :: var_p(3), mode_p(9)   ! var and mode to the Nth power
real :: AA, BB, CC, DD, EE

! Computation savers - powers are computationally expensive
var_p(1) = var
do i=2, 3
   var_p(i) = var_p(i-1)*var
enddo

mode_p(1) = mode
do i = 2, 9
  mode_p(i) = mode_p(i-1)*mode
enddo

! Calculate the rate parameter for IG distribution.
! It's a function of both the prior mean and variannce,
! obtained as a "real" solution to a cubic polynomial.
AA = mode_p(4) * sqrt((var_p(2) + 47.0*var*mode_p(2) + 3.0*mode_p(4)) / var_p(3))
BB = 75.0*var_p(2)*mode_p(5)
CC = 21.0*var*mode_p(7)
DD = var_p(3)*mode_p(3)
EE = (CC + BB + DD + mode_p(9) + 6.0*sqrt(3.0)*AA*var_p(3)) / var_p(3)

beta = (7.0*var*mode + mode_p(3))/(3.0*var) + &
       EE**(1.0/3.0)/3.0 + mode_p(2)*(var_p(2) + 14.0*var*mode_p(2) + &
       mode_p(4)) / (3.0*var_p(2)*EE**(1.0/3.0))

end subroutine change_GA_IG

