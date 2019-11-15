program alignment
use constants
use namelist_define
use mpi_module
use mapinfo_define
use netcdf
use multiscale_utils
!!!find displacement vectors q from scale-s updates (500??->700??)
!!!warp smaller scale component (600??) with q
!!!assemble increments = warped 600?? - 600?? + 700?? - 500??
!!!  and add to final output 900??
integer            :: b_unit=50010,a_unit=70010,s_unit=60010
integer            :: i_unit=80010,o_unit=90010
integer            :: ie,ix,jx,kx,len,i,j,k,n
integer            :: ii,jj,kk,nn,m,nv,nm,fid
real               :: member_per_cpu
character (len=10) :: wrf_file
real, allocatable, dimension(:,:) :: u,v
real, allocatable, dimension(:,:,:) :: xb,xa,xs,xsw,xout
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

!!only using nmcpu here, no domain decomposition (because requires whole domain)
member_per_cpu=real(numbers_en+1)/nmcpu
nm=ceiling(member_per_cpu)
!if ( my_proc_id == 0 ) then
!  write(*,*) '---------------------------------------------------'
!  print*, num_scales
!  print*, krange
!endif

!-- allocate wrf variables
nv = 0
do m=1,30
  if(len_trim(enkfvar(m))>=1) nv=nv+1
enddo

!!!!!BETTER PARALLEL:
!!!! loop over members and variables????
do n=1,nm
  ie=(n-1)*my_proc_id+1
  if(ie<=numbers_en) then
    do m=1,nv
      write(wrf_file,'(a5,i5.5)') 'fort.',b_unit+ie
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      allocate(u(ii,jj),v(ii,jj))
      allocate(xb(ii,jj,kk),xa(ii,jj,kk),xs(ii,jj,kk),xsw(ii,jj,kk),xout(ii,jj,kk))
      xb=0.; xa=0.; xs=0.; xsw=0.; xout=0.

      !!! Read in prior ensemble at current scale
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xb)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xb)
      endif

      !!! Read in posterior ensemble at current scale
      write(wrf_file,'(a5,i5.5)') 'fort.',a_unit+ie
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xa)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xa)
      endif

      !!!read smaller scale component xs
      write(wrf_file,'(a5,i5.5)') 'fort.',s_unit+ie
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xs)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xs)
      endif

      do k=1,kk
        !!! Compute displacement vectors using Horn Schunck
        call optical_flow_HS(xb(:,:,k),xa(:,:,k),100.,u,v)
        !!! Warp smaller-scale field xs -> xsw
        call warp_state(xs(:,:,k),xsw(:,:,k),u,v)
      enddo

      !!!!read current intermediate state
      write(wrf_file,'(a5,i5.5)') 'fort.',o_unit+ie
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xout)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xout)
      endif

      !!!!add increment and write back
      xout=xout+xsw-xs!+xa-xb

      call open_file(wrf_file,nf_write,fid)
      if(kk>1) then
        call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,xout)
      else if (kk==1) then
        call write_variable2d(fid,enkfvar(m),ii,jj,1,xout)
      endif
      call close_file(fid)

      deallocate(u,v,xa,xb,xs,xsw,xout)
    enddo
  endif
enddo
if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of alignment.exe'
call parallel_finish()

end program alignment

