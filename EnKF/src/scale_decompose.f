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
integer            :: ie,ix,jx,kx,len,i,j,k,n
integer            :: ii,jj,kk,nn,m,nv,fid
character (len=10) :: wrf_file
real, allocatable, dimension(:,:,:) :: x,xout,xout1
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

!!only using nmcpu=nens here, no domain decomposition (because fft2 requires whole domain)
!!ignore namelist parallel options
if(nprocs .ne. numbers_en) then
  if(my_proc_id==0) print*,'ERROR: nprocs/=numbers_en'
endif

nv = 0
do m=1,30
  if(len_trim(enkfvar(m))>=1) nv=nv+1
enddo

!!!krange should have ns-1 entries, each is the wavenumber separating scales.
if ( my_proc_id == 0 ) then
  write(*,*) '---------------------------------------------------'
  print*, 'num_scales=',num_scales
  print*, 'current_scale=',current_scale
  print*, 'krange=',krange(1:num_scales-1)
endif

ie=my_proc_id+1
if(ie<=numbers_en) then
  do m=1,nv
    if(my_proc_id==0) print*, 'processing ', enkfvar(m)

    write(wrf_file,'(a5,i5.5)') 'fort.',i_unit+ie
    call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
    allocate(x(ii,jj,kk),xout(ii,jj,kk),xout1(ii,jj,kk))
    x=0.

    !!! Read in initial Ensemble in x
    if(kk>1) then
      call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,x)
    else if (kk==1) then
      call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,x)
    endif

    !!! Get scale component
    xout=x
    xout1=x
    nn=min(ii,jj)
    do k=1,kk
      call scale_bandpass(x(1:nn,1:nn,k),krange(1:num_scales-1), &
                          current_scale,xout(1:nn,1:nn,k),xout1(1:nn,1:nn,k))
    enddo

    !!!output netcdf file
    write(wrf_file,'(a5,i5.5)') 'fort.',s_unit+ie
    call open_file(wrf_file,nf_write,fid)
    if(kk>1) then
      call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,xout)
    else if (kk==1) then
      call write_variable2d(fid,enkfvar(m),ii,jj,1,xout)
    endif
    call close_file(fid)
    write(wrf_file,'(a5,i5.5)') 'fort.',s1_unit+ie
    call open_file(wrf_file,nf_write,fid)
    if(kk>1) then
      call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,xout1)
    else if (kk==1) then
      call write_variable2d(fid,enkfvar(m),ii,jj,1,xout1)
    endif
    call close_file(fid)
    deallocate(x,xout,xout1)
  enddo
endif

if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of scale_decompose.exe'

call parallel_finish()

end program scale_decompose
