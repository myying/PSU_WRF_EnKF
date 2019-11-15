program scale_decompose
use constants
use namelist_define
use mpi_module
use mapinfo_define
use netcdf
use multiscale_utils
!!!read intermediate updated state 900??, decompose into
!!!current scale copmonent 500?? and smaller scale component 600??
integer            :: i_unit=90010, s_unit=50010, s1_unit=60010
integer            :: unit
integer            :: ie,ix,jx,kx,len,i,j,k,n,ni,nj,nk
integer            :: ii,jj,kk,nn,m,nv,nm,fid
real               :: member_per_cpu
character (len=10) :: wrf_file
real, allocatable, dimension(:,:) :: dat2d
real, allocatable, dimension(:,:,:) :: dat3d
real, allocatable, dimension(:,:,:,:,:) :: x,xout,xout1
type(proj_info)    :: proj
character (len=80) :: times
!----------------------------------------------------------------
! Initialize parallel stuff
call parallel_start()
! Get WRF model information from the 1th member
! ( reference to MM5_to_GrADS: module_wrf_to_grads_util.f)
write(wrf_file,'(a5,i5.5)')'fort.',i_unit+1
call get_wrf_info(wrf_file, ix, jx, kx, times, proj) ! saved wrf info to proj
! Initilize and Read in namelist information
call read_namelist(ix, jx, kx)

!!only using nmcpu here, no domain decomposition (because fft2 requires whole domain)
member_per_cpu=real(numbers_en+1)/nmcpu
nm=ceiling(member_per_cpu)

nv = 0
do m=1,30
  if(len_trim(enkfvar(m))>=1) nv=nv+1
enddo
ni=ix+1
nj=jx+1
nk=kx+1
allocate(x(ni,nj,nk,nv,nm),xout(ni,nj,nk,nv,nm),xout1(ni,nj,nk,nv,nm))
x=0.

!!!krange should have ns-1 entries, each is the wavenumber separating scales.
if ( my_proc_id == 0 ) then
  write(*,*) '---------------------------------------------------'
  print*, num_scales
  print*, krange(1:num_scales-1)
endif

!!! Read in initial Ensemble in x
if ( my_proc_id == 0 ) write(*,*)'Read in initial Ensemble in x'
do n=1,nm
  ie=(n-1)*my_proc_id+1
  if(ie<=numbers_en) then
    write(wrf_file,'(a5,i5.5)') 'fort.',i_unit+ie
    do m=1,nv
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      allocate(dat3d(ii,jj,kk))
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,dat3d)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,dat3d)
      endif
      x(1:ii,1:jj,1:kk,m,n)=dat3d(1:ii,1:jj,1:kk)
      deallocate(dat3d)
    enddo
  endif
enddo

!!! Get scale component
xout=x
xout1=x
do n=1,nm
  ie=(n-1)*my_proc_id+1
  if(ie<=numbers_en) then
    do m=1,nv
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      nn=min(ii,jj)
      do k=1,kk
        call scale_bandpass(x(1:nn,1:nn,k,m,n),krange(1:num_scales-1), &
                            current_scale,xout(1:nn,1:nn,k,m,n),xout1(1:nn,1:nn,k,m,n))
      enddo
    enddo
  endif
enddo

!!! output netcdf file
do n=1,nm
  ie=(n-1)*my_proc_id+1
  if(ie<=numbers_en) then
    write(wrf_file,'(a5,i5.5)') 'fort.',s_unit+ie
    call open_file(wrf_file,nf_write,fid)
    do m=1,nv
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      allocate(dat3d(ii,jj,kk))
      dat3d=0.
      dat3d(1:ii,1:jj,1:kk)=xout(1:ii,1:jj,1:kk,m,n)
      if(kk>1) then
        call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,dat3d)
      else if (kk==1) then
        call write_variable2d(fid,enkfvar(m),ii,jj,1,dat3d)
      endif
      deallocate(dat3d)
    enddo
    call close_file(fid)
  endif
enddo
do n=1,nm
  ie=(n-1)*my_proc_id+1
  if(ie<=numbers_en) then
    write(wrf_file,'(a5,i5.5)') 'fort.',s1_unit+ie
    call open_file(wrf_file,nf_write,fid)
    do m=1,nv
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      allocate(dat3d(ii,jj,kk))
      dat3d=0.
      dat3d(1:ii,1:jj,1:kk)=xout1(1:ii,1:jj,1:kk,m,n)
      if(kk>1) then
        call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,dat3d)
      else if (kk==1) then
        call write_variable2d(fid,enkfvar(m),ii,jj,1,dat3d)
      endif
      deallocate(dat3d)
    enddo
    call close_file(fid)
  endif
enddo

if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of scale_decompose.exe'

deallocate(x,xout,xout1)
call parallel_finish()

end program scale_decompose
