program alignment
use constants
use namelist_define
use mpi_module
use mapinfo_define
use map_utils
use netcdf
use multiscale_utils
!!!find displacement vectors q from scale-s updates (500??->700??)
!!!warp smaller scale component (600??) with q
!!!assemble increments = warped 600?? - 600?? + 700?? - 500??
!!!  and add to final output 900??
integer            :: b_unit=50010,a_unit=70010,s_unit=60010,o_unit=90010
integer            :: num_dom=2, dm
integer            :: ie,ix,jx,kx,len,i,j,k,n
integer            :: ii,jj,kk,nn,m,nv,fid
integer            :: ii1,jj1,kk1,ix1,jx1,kx1
real               :: i_start,j_start,grid_ratio,u1,v1,i0,j0,i1,j1,i1_new,j1_new
character (len=10) :: wrf_file
character (len=17) :: wrf_file1
real, allocatable, dimension(:,:) :: u,v
real, allocatable, dimension(:,:,:) :: xb,xa,xs,xsw,xin,xout,delx,x1,x1w
type(proj_info)    :: proj,proj1
character (len=80) :: times
!----------------------------------------------------------------
! Initialize parallel stuff
call parallel_start()
! Get WRF model information from the 1th member
! ( reference to MM5_to_GrADS: module_wrf_to_grads_util.f)
write(wrf_file,'(a5,i5.5)')'fort.',o_unit+1
call get_wrf_info(wrf_file, ix, jx, kx, times, proj)
! Initilize and Read in namelist information
call read_namelist(ix, jx, kx)

!!only using nmcpu=nens here, no domain decomposition (because requires whole domain)
!!ignore namelist parallel options
if(nprocs .ne. numbers_en) then
  if(my_proc_id==0) print*,'ERROR: nprocs/=numbers_en'
endif

!-- allocate wrf variables
nv = 0
do m=1,30
  if(len_trim(enkfvar(m))>=1) nv=nv+1
enddo
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

    write(wrf_file,'(a5,i5.5)') 'fort.',b_unit+ie
    call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
    allocate(u(ii,jj),v(ii,jj))
    allocate(xb(ii,jj,kk),xa(ii,jj,kk),xs(ii,jj,kk),xsw(ii,jj,kk))
    allocate(xin(ii,jj,kk),xout(ii,jj,kk),delx(ii,jj,kk))
    xb=0.; xa=0.; xs=0.; xsw=0.; xin=0.; xout=0.; delx=0.

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

    !!!!read current intermediate state
    write(wrf_file,'(a5,i5.5)') 'fort.',o_unit+ie
    if(kk>1) then
      call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xin)
    else if (kk==1) then
      call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xin)
    endif


    if(current_scale<num_scales) then
      !!!read smaller scale component xs
      write(wrf_file,'(a5,i5.5)') 'fort.',s_unit+ie
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,xs)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,xs)
      endif

      !!! Compute displacement vectors using Horn Schunck
      call optical_flow_HS(sum(xb,3)/real(kk),sum(xa,3)/real(kk),100.,u,v)

      !!! Warp smaller-scale field xs -> xsw
      xsw=xs
      buffer=4
      do k=1,kk
        do i=1+buffer,ii-buffer
        do j=1+buffer,jj-buffer
          xsw(i,j,k) = interp2d(xs(:,:,k),real(i)-u(i,j),real(j)-v(i,j))
        enddo
        enddo
      enddo
      if(num_dom>1) then
        do dm=2,num_dom
          write(wrf_file1,'(a5,i1,a6,i5.5)') '../d0',dm,'/fort.',o_unit+ie
          call get_wrf_info(wrf_file1, ix1, jx1, kx1, times, proj1)
          call wrf_var_dimension(wrf_file1,enkfvar(m),ix1,jx1,kx1,ii1,jj1,kk1)
          allocate(x1(ii1,jj1,kk),x1w(ii1,jj1,kk))
          if(kk>1) then
            call get_variable3d(wrf_file1,enkfvar(m),ii1,jj1,kk,1,x1)
          else if (kk==1) then
            call get_variable2d(wrf_file1,enkfvar(m),ii1,jj1,1,x1)
          endif
          x1w=x1
          !!find nested domain location
          call latlon_to_ij(proj,proj1%lat1,proj1%lon1,i_start,j_start)
          grid_ratio=real(proj%dx)/real(proj1%dx)
          do k=1,kk
            do i1=1,ii1
            do j1=1,jj1
              !!!find displacement vector at i,j
              i0=i_start+i1/grid_ratio
              j0=j_start+j1/grid_ratio
              u1=interp2d(u,i0,j0)*grid_ratio
              v1=interp2d(v,i0,j0)*grid_ratio
              i1_new=i1-u1
              j1_new=j1-v1
              !!!find warped state
              if(i1_new<1 .or. i1_new>ii1 .or. j1_new<1 .or. j1_new>jj1) then
                x1w(i1,j1,k)=interp2d(xin(:,:,k),i0-u(int(i0),int(j0)),j0-v(int(i0),int(j0)))
              else
                x1w(i1,j1,k)=interp2d(x1(:,:,k),i1_new,j1_new)
              endif
            end do
            end do
          end do
          !!write warped file back
          call open_file(wrf_file1,nf_write,fid)
          if(kk>1) then
            call write_variable3d(fid,enkfvar(m),ii1,jj1,kk,1,x1w)
          else if (kk==1) then
            call write_variable2d(fid,enkfvar(m),ii1,jj1,1,x1w)
          endif
          call close_file(fid)
          deallocate(x1,x1w)
        end do
      end if
    end if

    !!!!add increment and write back
    delx=xa-xb
    if(run_alignment .and. current_scale<num_scales) delx=delx+xsw-xs
    xout=xin+delx

    write(wrf_file,'(a5,i5.5)') 'fort.',o_unit+ie
    call open_file(wrf_file,nf_write,fid)
    if(kk>1) then
      call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,xout)
    else if (kk==1) then
      call write_variable2d(fid,enkfvar(m),ii,jj,1,xout)
    endif
    call close_file(fid)

    deallocate(u,v,xa,xb,xs,xsw,xin,xout,delx)
  enddo
endif
if ( my_proc_id == 0 ) write(*,'(a)')' Successful completion of alignment.exe'
call parallel_finish()

end program alignment

