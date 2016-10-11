Program dualresolution_enkf

  use constants
  use namelist_define
  use mpi_module
  use mapinfo_define
  use netcdf
  use obs_define
  use map_utils
  use wrf_tools
  use radar
  
!----------------------------------------------------------------
  implicit none

  integer            :: xfl_unit=90001, xa_unit=90010, xal_unit
  integer            :: xfh_unit=90002

  integer            :: ixl, jxl, kxl, ixh, jxh, kxh 
  integer            :: iil, jjl, kkl, iih, jjh, kkh, n, fid, rcode
  integer            :: i_parent_start, j_parent_start, i_parent_end, j_parent_end, parent_grid_ratio
  real               :: dx, dy

  character (len=10)                    :: xfl_file, xfh_file, xal_file
  real, allocatable, dimension(:,:,:)   :: xfl, xal, xfh, xah

  integer            :: i, j, k, il, jl, il1, jl1, m, num_var 


!-----------------------------------------------------------------------------
!  get domain dimension
  write(xfl_file,'(a5,i5.5)')'fort.',xfl_unit
!  call get_var_ij ( xfl_file, 'T         ', ixl, jxl, kxl )
   call get_ij ( xfl_file, ixl, jxl, kxl )

  write(xfh_file,'(a5,i5.5)')'fort.',xfh_unit
!  call get_var_ij ( xfh_file, 'T         ', ixh, jxh, kxh )
   call get_ij ( xfh_file, ixh, jxh, kxh )
  call open_file ( xfh_file, nf_nowrite, fid )
  rcode = nf_get_att_int(fid, nf_global, 'I_PARENT_START', i_parent_start)
  rcode = nf_get_att_int(fid, nf_global, 'J_PARENT_START', j_parent_start)
  rcode = nf_get_att_int(fid, nf_global, 'PARENT_GRID_RATIO', parent_grid_ratio)
  call close_file ( fid )
  i_parent_end = (ixh-1)/parent_grid_ratio+i_parent_start
  j_parent_end = (jxh-1)/parent_grid_ratio+j_parent_start
  

!----------------------------------------------------------------
! Initilize and Read in namelist information
  call read_namelist(ixl, jxl, kxl)
  xal_unit = xa_unit + numbers_en
  write(xal_file,'(a5,i5.5)')'fort.',xal_unit

!----------------------------------------------------------------
! variables loop
  num_var = 0
  do m = 1, 20
     if ( len_trim( enkfvar(m) ) >= 1 ) num_var = num_var +1
  enddo
  do n = 1, num_var

     call wrf_var_dimension ( xal_file, enkfvar(n), ixl, jxl, kxl, iil, jjl, kkl )
     call wrf_var_dimension ( xal_file, enkfvar(n), ixh, jxh, kxh, iih, jjh, kkh )
   
     allocate ( xfl ( iil, jjl, kkl ) )
     allocate ( xal ( iil, jjl, kkl ) )
     allocate ( xfh ( iih, jjh, kkh ) )
     allocate ( xah ( iih, jjh, kkh ) )
     
     if ( kkl .ne. kkh ) stop ' vertical leves not match!!!'
     if ( kkl > 1 ) then
        call get_variable3d( xfl_file, enkfvar(n), iil, jjl, kkl, 1, xfl )
        call get_variable3d( xal_file, enkfvar(n), iil, jjl, kkl, 1, xal )
        call get_variable3d( xfh_file, enkfvar(n), iih, jjh, kkh, 1, xfh )
     else if ( kkl == 1 ) then
        call get_variable2d( xfl_file, enkfvar(n), iil, jjl, 1, xfl )
        call get_variable2d( xal_file, enkfvar(n), iil, jjl, 1, xal )
        call get_variable2d( xfh_file, enkfvar(n), iih, jjh, 1, xfh )
     endif
     xal = xal - xfl

     do k = 1, kkl
        do j = 1, jjh
        do i = 1, iih
             il = int((i-1)/parent_grid_ratio) + i_parent_start
             jl = int((j-1)/parent_grid_ratio) + j_parent_start
             il1 = il + 1
             jl1 = jl + 1
             dx = real(i-1)/real(parent_grid_ratio) - real(il-i_parent_start) 
             dy = real(j-1)/real(parent_grid_ratio) - real(jl-j_parent_start) 
             xah(i,j,k) = (1.-dy)*((1.0-dx)*xal(il,jl  ,k) + dx*xal(il+1,jl  ,k))    &
                            + dy *((1.0-dx)*xal(il,jl+1,k) + dx*xal(il+1,jl+1,k))
             xah(i,j,k) = xah(i,j,k) + xfh(i,j,k)
        enddo
        enddo
     enddo

     call open_file ( xfh_file, nf_write, fid )
     if ( kkl > 1 ) then
        call write_variable3d (fid, enkfvar(n), iih, jjh, kkh, 1, xah )
     else if ( kkl == 1 ) then
        call write_variable2d (fid, enkfvar(n), iih, jjh,  1, xah )
     endif
     call close_file ( fid )
 
     deallocate ( xfl )
     deallocate ( xal )
     deallocate ( xfh )
     deallocate ( xah )

  enddo   !do n = 1, num_var
     
end     
!========================================================================================


