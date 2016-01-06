!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   module netcdf - module that contains several subroutines used to 
!                   read and write netcdf data
!
!     created Oct. 2004 Ryan Torn, U. Washington
!     add some subroutines 07/2005 Zhiyong Meng TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      MODULE NETCDF

        include 'netcdf.inc'

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   close_file - subroutine that closes a netcdf file
!
!      fid - integer file id
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine close_file(fid)

      integer, intent(in) :: fid

      integer :: rcode

      rcode = nf_close(fid)
      if ( rcode .ne. 0 ) then
         call netcdf_error('File close', 'close', rcode, 0, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_bdy_lev - subroutine that reads one vertical level of WRF 
!                 boundary data
!
!      fid - integer file id
!      var - character name of variable
!      ixy - number of grid points in x or y
!       iz - number of grid points in z
!      wid - number of grid points in width
!       it - time to get
!      dat - data read from file
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bdy_lev(fid, var, ixy, iz, wid, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ixy, iz, wid, it 
      real, intent(out)             :: dat(ixy,wid)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write data
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ixy
      istart(2) = iz
      iend(2) = 1
      istart(3) = 1
      iend(3) = wid
      istart(4) = it
      iend(4) = 1

      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode .ne. 0) then
        call netcdf_error('get_bdy_lev', var, rcode, ixy, iz, wid, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_dimlen - function that reading the length of a dimension
!
!       fid - integer file id
!   dimname - name of dimension
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function get_dimlen(fid, dimname)

      character (len=*), intent(in) :: dimname
      integer, intent(in)           :: fid

      integer :: did, rcode, get_dimlen

      rcode = nf_inq_dimid(fid, dimname, did)
      rcode = nf_inq_dimlen(fid, did, get_dimlen)

      return
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable0d - subroutine that reads netcdf data for a real    
!
!     fid - integer file id
!     var - name of the variable to read
!      it - time to read the data from
!     dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable0d(file, var, it, dat)

      character (len=*), intent(in) :: file, var
      integer, intent(in)           ::  it
      real, intent(out)             :: dat
      integer                       :: fid
      integer :: istart, iend, varid, rcode
      
      call open_file (file, nf_nowrite, fid)
      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart = 1
      iend = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable0d', var, rcode, 1, 0, 0, it)
      endif

      call close_file(fid)

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d - subroutine that reads netcdf data for a 1    
!                    dimensional quantity
!
!     fid - integer file id
!     var - name of the variable to read
!      iz - z dimension size
!      it - time to read the data from
!     dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d(file, var, iz, it, dat)

      character (len=*), intent(in) :: file, var
      integer, intent(in)           ::  iz, it
      real, intent(out)             :: dat(iz)
      integer                       :: fid
      integer :: istart(2), iend(2), varid, rcode
      
      call open_file (file, nf_nowrite, fid)
      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = iz
      istart(2) = it
      iend(2) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable1d', var, rcode, iz, 0, 0, it)
      endif

      call close_file(fid)

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d_local - subroutine that reads netcdf data for a 1    
!                          dimensional quantity over an interval
!
!       fid - integer file id
!       var - name of the variable to read
!       izs - z dimension start
!       ize - z dimension end
!        it - time to read the data from
!       dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d_local(fid, var, izs, ize, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, izs, ize, it
      real, intent(out)             :: dat(ize-izs+1)

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = izs
      iend(1) = ize-izs+1
      istart(2) = it
      iend(2) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'izs', istart(1)
        call netcdf_error('get_variable1d_local', var, rcode, 
     &                         iend(1), 0, 0, it) 
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d_int - subroutine that reads netcdf integer data 
!                        for a 1 dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       iz - z dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d_int(file, var, iz, it, dat)

      character (len=*), intent(in) :: file, var
      integer, intent(in)           :: iz, it
      integer, intent(out)          :: dat(iz)

      integer :: fid, istart(2), iend(2), varid, rcode

      call open_file (file, nf_nowrite, fid)
      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = iz
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_get_vara_int(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable1d_int', var, rcode, iz, 0, 0,it)
      endif
      call close_file(fid)

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable2d - subroutine that reads netcdf data for a 2    
!                    dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - x dimension size
!       i2 - y dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable2d(file, var, i1, i2, it, dat)

      character (len=*), intent(in) :: file, var
      integer, intent(in)           :: i1, i2, it
      real, intent(out)             :: dat(i1,i2)
      integer                       :: fid

      integer :: istart(3), iend(3), varid, rcode

      call open_file (file, nf_nowrite, fid)
      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = it
      iend(3) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable2d', var, rcode, i1, i2, 0, it)
      endif

      call close_file(fid)

      return
      end subroutine get_variable2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable2d_local - subroutine that reads netcdf data for a 2    
!                          dimensional quantity over a patch
!
!     fid - integer file id
!     var - name of the variable to read
!     i1s - x dimension start
!     i1e - x dimension end
!     i2s - x dimension start
!     i2e - x dimension end
!      it - time to read the data from
!     dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable2d_local(fid, var, 
     &                           i1s, i1e, i2s, i2e, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1s, i1e, i2s, i2e, it
      real, intent(out)             :: dat(i1e-i1s+1,i2e-i2s+1)

      integer :: istart(3), iend(3), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = i1s
      iend(1) = i1e-i1s+1
      istart(2) = i2s
      iend(2) = i2e-i2s+1
      istart(3) = it
      iend(3) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'i1s i2s', istart(1:2)
        call netcdf_error('get_variable2d_local', var, rcode, 
     &                             iend(1), iend(2), 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable3d - subroutine that reads netcdf data for a 3 
!                    dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - x dimension size
!       i2 - y dimension size
!       i3 - z dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable3d(file, var, i1, i2, i3, it, dat)

      character (len=*), intent(in) :: file, var
      integer, intent(in)           :: i1, i2, i3, it
      real, intent(out)             :: dat(i1,i2,i3)
      integer                       :: fid, length
      character :: var_name(NF_MAX_NAME)
      integer :: istart(4), iend(4), varid, rcode

      !length = len_trim(var)
      call open_file(file, nf_nowrite, fid)

      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = 1
      iend(3) = i3
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable3d', var, rcode, i1, i2, i3, it)
      endif

      call close_file(fid)
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable3d_local - subroutine that reads netcdf data for a
!                          patch of a 3 dimensional quantity
!
!       fid - integer file id
!       var - name of the variable to read
!       i1s - x dimension start
!       i1e - x dimension end
!       i2s - y dimension start
!       i2e - y dimension end
!       i3s - x dimension start
!       i3e - x dimension end
!        it - time to read the data from
!       dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable3d_local(fid, var, 
     &                          i1s, i1e, i2s, i2e, i3s, i3e, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid,i1s,i1e,i2s,i2e,i3s,i3e,it
      real, intent(out)           :: dat(i1e-i1s+1,i2e-i2s+1,i3e-i3s+1)

      integer :: istart(4), iend(4), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = i1s
      iend(1) = i1e-i1s+1
      istart(2) = i2s
      iend(2) = i2e-i2s+1
      istart(3) = i3s
      iend(3) = i3e-i3s+1
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'i1s i2s i3s', istart(1:3)
        call netcdf_error('get_variable3d_local', var, rcode, 
     &                       iend(1), iend(2), iend(3), it)
      endif
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_lev - subroutine that reads netcdf data for one vertical 
!                      level of a 3-D quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       ix - x dimension size
!       iy - y dimension size
!       iz - z level
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_lev(fid, var, ix, iy, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, iy, iz, it
      real, intent(out)             :: dat(ix,iy)

      integer :: istart(4), iend(4), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = ix
      istart(2) = 1
      iend(2) = iy
      istart(3) = iz
      iend(3) = 1
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_lev', var, rcode, ix, iy, iz,it)
      endif
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_vec - subroutine that reads netcdf data for a vector    
!                      of data
!
!       fid - integer file id
!       var - name of the variable to read
!        iz - length of the vector of data
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_vec(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz
      real, intent(out)             :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart = 1
      iend = iz

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_veci - subroutine that reads netcdf data for a vector    
!                       of integer data
!
!       fid - integer file id
!       var - name of the variable to read
!        iz - length of the vector of data
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_veci(fid, var, iz, dat)

      character (len=*), intent(in)  :: var
      integer, intent(in)            :: fid, iz
      integer, intent(out)           :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart = 1
      iend = iz

      ! read the data
      rcode = nf_get_vara_int(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_vec_local - subroutine that reads netcdf data for a     
!                            vector of data over a limited area
!
!       fid - integer file id
!       var - name of the variable to read
!    istart - index of beginning of data to read
!        ie - index of end of data to read
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_vec_local(fid, var, istart, ie, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, istart, ie
      real, intent(out)             :: dat(ie-istart+1)

      integer :: iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      iend = ie-istart+1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec_local',var,rcode,istart,
     &                            ie,0,0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_text - subroutine that reads netcdf text data for 1 dimensional 
!              quantity 
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - length of string to read
!       it - time to read the data from
!  outtext - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_text(fid, var, i1, it, outtxt)

      character (len=*), intent(in)   :: var
      integer, intent(in)             :: fid, i1, it
      character (len=i1), intent(out) :: outtxt

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_get_vara_text(fid, varid, istart, iend, outtxt)

      if (rcode.ne.0) then 
        call netcdf_error('get_text', var, rcode, i1, 0, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   netcdf_error - subroutine that prints an error and stops a program 
!                  if a netcdf routine throws an error
!
!     place - program where the netcdf error happens
!       var - variable of error
!        rc - integer code of error
!        ix - dimension 1 of error
!        iy - dimension 2 of error
!        iz - dimension 3 of error
!        it - time dimension of error
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine netcdf_error(place, var, rc, ix, iy, iz, it)

      character (len=*), intent(in) :: place, var
      integer, intent(in)           :: rc, ix, iy, iz, it

      write(6,*) 'NETCDF ERROR !!!!!!!!'
      write(6,*)            
      write(6,*) 'Error occurs at ',place
      write(6,*) 'With variable ','"'//var//'"'
      write(6,*) 'Netcdf message ',nf_strerror(rc)
      write(6,*) 'ix iy iz it',ix,iy,iz,it

      open(25, file='enkf_error',status='unknown')
      write(25,*) nf_strerror(rc); close(25)

      stop
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   open_file - subroutine that opens a netcdf file 
!
!   filename - name of file to open
!    permiss - permissions of file to open
!        fid - integer file id
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine open_file(filename, permiss, fid)

      character (len=*), intent(in) :: filename
      integer, intent(in)           :: permiss
      integer, intent(out)          :: fid

      integer :: rcode

      rcode = nf_open(filename, permiss, fid)
      if ( rcode .ne. 0 ) then
         call netcdf_error(filename, 'open', rcode, 0, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_bdy_lev - subroutine that writes one vertical level of 
!                   boundary data
!
!     fid - integer file id
!     var - name of variable
!     ixy - number of grid points in x or y
!      iz - vertical level to write out
!     wid - width of boundary data
!      it - time to write
!     dat - data to write out
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_bdy_lev(fid, var, ixy, iz, wid, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ixy, iz, wid, it
      real, intent(in)              :: dat(ixy,wid)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write data
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ixy
      istart(2) = iz
      iend(2) = 1
      istart(3) = 1
      iend(3) = wid
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)

      if (rcode .ne. 0) then
        call netcdf_error('write_bdy_lev', var, rcode, ixy, iz, wid, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_text - subroutine that writes netcdf text data for 1 
!                dimensional quantity 
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - length of string to write
!       it - time to write the data from
!  outtext - array containing data to be written
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_text(fid, var, i1, it, outtxt)

      character (len=*), intent(in)   :: var
      integer, intent(in)             :: fid, i1, it
      character (len=i1), intent(in ) :: outtxt

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_put_vara_text(fid, varid, istart, iend, outtxt)

      if (rcode.ne.0) then 
        call netcdf_error('write_text', var, rcode, i1, 0, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable1d - subroutine that writes netcdf dat for a 1
!                          dimensional quantity
!
!       fid - integer file id
!       var - name of the variable to write
!        ix - x dimension size
!        it - time to write out data
!       dat - array containing dat to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable1d(fid, var, ix, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, it
      real, intent(in)              :: dat(ix)

      integer :: istart(2), iend(2), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = it
      iend(2) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable1d', var, rcode, ix, 0, 0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable1d_int - subroutine that writes netcdf dat for a 1
!                          dimensional integer quantity
!
!     fid - integer file id
!     var - name of the variable to write
!      ix - x dimension size
!      it - time to write out data
!     dat - array containing dat to write
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable1d_int(fid, var, ix, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, it, dat(ix)
      
      integer :: istart(2), iend(2), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = it
      iend(2) = 1

      rcode = nf_put_vara_int(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable1d_int', var, rcode, ix, 0, 
     &                         0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable2d - subroutine that writes netcdf dat for a 2
!                      dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - x dimension size
!       i2 - y dimension size
!       it - time to write out data
!      dat - array containing dat to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable2d(fid, var, i1, i2, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, it
      real, intent(in)              :: dat(i1,i2)

      integer :: istart(3), iend(3), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = it
      iend(3) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable2d', var, rcode, i1, i2, 0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable3d - subroutine that writes netcdf data for a 3 
!                      dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - x dimension size
!       i2 - y dimension size
!       i3 - z dimension size
!       it - time to write data
!      dat - array containing data to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable3d(fid, var, i1, i2, i3, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, i3, it
      real, intent(in)              :: dat(i1,i2,i3)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = 1
      iend(3) = i3
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable3d', var, rcode, i1, i2, i3,it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_var_lev - subroutine that writes netcdf data for one vertical 
!                   level of a 3D field
!
!     fid - integer file id
!     var - name of the variable to write
!      ix - x dimension size
!      iy - y dimension size
!      iz - vertical level to write
!      it - time to write out data
!     dat - array containing data to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_var_lev(fid, var, ix, iy, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, iy, iz, it
      real, intent(in)              :: dat(ix,iy)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = 1
      iend(2) = iy
      istart(3) = iz
      iend(3) = 1
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_var_lev', var, rcode, ix, iy, iz,it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_vec -  subroutine that writes netcdf data for a 
!                         vector of data
!
!      fid - integer file id
!      var - name of the variable to write
!       iz - x dimension size
!      dat - array containing dat to write
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_vec(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz
      real, intent(in)              :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart = 1
      iend   = iz

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_vec -  subroutine that writes netcdf data for a 
!                         vector of integer data
!
!      fid - integer file id
!      var - name of the variable to write
!       iz - x dimension size
!      dat - array containing dat to write
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_veci(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz, dat(iz)

      integer :: istart, iend, varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart = 1
      iend   = iz

      rcode = nf_put_vara_int(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    get ix jx from a netcdf file
!    created 07/2005 Zhiyong Meng TAMU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_ij (file,ix,iy,iz)
      character (len=*), intent(in) :: file
      integer :: fid, rcode, mapp, ix, iy, iz

      call open_file(file, nf_nowrite, fid)
          rcode = nf_get_att_int(fid, nf_global,
     &                            'WEST-EAST_GRID_DIMENSION', ix)
          rcode = nf_get_att_int(fid, nf_global,
     &                            'SOUTH-NORTH_GRID_DIMENSION', iy)
          rcode = nf_get_att_int(fid, nf_global,
     &                            'BOTTOM-TOP_GRID_DIMENSION', iz)
      call close_file(fid)
      ix = ix -1
      iy = iy - 1
      iz = iz - 1
      end subroutine get_ij
!!!!!!
!--------get soil variable levels
!!!!!!
      subroutine get_soilkk(file,k)
      character (len=*), intent(in) :: file
      integer :: fid, rcode, dimid,k 
      call open_file(file, nf_nowrite, fid)
      rcode=nf_inq_dimid(fid,'soil_layers_stag',dimid)
      rcode=nf_inq_dimlen(fid,dimid,k)
      call close_file(fid)
      end subroutine get_soilkk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    get any variable in 3D format 
!    created 07/2005 Zhiyong Meng TAMU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_var (file, plot_var, imx,imy,imz,data)
      implicit none
      character (LEN =  *)        :: file
      integer                      :: ix,iy,iz,i,j,k,imx,imy,imz
      integer                      :: fid,rcode,id_var
      real                         :: data(imx,imy,imz)
      character (len=*)           :: plot_var
      integer istart(4), iend(4)
     
      istart=1
      iend  =1
      iend(1)=imx
      iend(2)=imy
      iend(3)=imz
      call open_file (file, nf_nowrite, fid)
      rcode = nf_inq_varid ( fid, plot_var, id_var )
      call ncvgt( fid,id_var,istart,iend,data,rcode)
      end subroutine get_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    put variable
!    created 07/2005 Zhiyong Meng TAMU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine put_var (file, plot_var, imx,imy,imz,data)
      implicit none
      character (LEN =  *)        :: file
      integer                      :: ix,iy,iz,i,j,k,imx,imy,imz
      integer                      :: fid,rcode,id_var
      real                         :: data(imx,imy,imz)
      character (len=*)           :: plot_var
      integer istart(4), iend(4)

      istart=1
      iend  =1
      iend(1)=imx
      iend(2)=imy
      iend(3)=imz
      call open_file (file, nf_write, fid)
      rcode = nf_inq_varid ( fid, plot_var, id_var )
      call ncvpt( fid,id_var,istart,iend,data,rcode)
      end subroutine put_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    get the times of the data
!    created 07/2005 Zhiyong Meng  TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_times (file, plot_var, time)
      implicit none
      character (len=*) :: plot_var
      integer            :: i,j,k
      integer            :: fid,rcode,id_var
      integer            :: ndims, natts, ivtype,id_time, n_times,unlimDimID
      integer            :: dimids(10),dims(4)
      integer            :: istart_t(2), iend_t(2), itimes
      character (LEN =  *)        :: file
      character (len=80) :: times(100),print_time,varnam, time

      call open_file (file, nf_nowrite, fid)
      id_time = ncvid( fid, 'Times', rcode )
      rcode = nf_inq_var(fid,id_time,varnam,ivtype,ndims,dimids,natts)
      do i=1,ndims
        rcode = nf_inq_dimlen( fid, dimids(i), dims(i) )
      enddo

      n_times = dims(2)
      do i=1,dims(2)
        istart_t(1) = 1
        iend_t(1) = dims(1)
        istart_t(2) = i
        iend_t(2) = 1
        rcode = NF_GET_VARA_TEXT(fid,id_time,istart_t,iend_t,times(i))
      enddo

      itimes = n_times
      print_time = times(itimes)
      time=times(1)
    
!      print*, time

      call close_file(fid)
    
      end subroutine get_times 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    get the hour of the data
!    created 07/2005 Zhiyong Meng  TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_hour (file, plot_var, hour)
      implicit none
      character (len=*) :: plot_var
      integer            :: i,j,k
      integer            :: fid,rcode,id_var
      integer            :: ndims, natts, ivtype,id_time, n_times,unlimDimID
      integer            :: dimids(10),dims(4)
      integer            :: istart_t(2), iend_t(2), itimes
      character (len=2)  :: hour
      character (LEN =  *)        :: file
      character (len=80) :: times(100),print_time,varnam

      call open_file (file, nf_nowrite, fid)
      id_time = ncvid( fid, 'Times', rcode )
      rcode = nf_inq_var(fid,id_time,varnam,ivtype,ndims,dimids,natts)
      do i=1,ndims
        rcode = nf_inq_dimlen( fid, dimids(i), dims(i) )
      enddo

      n_times = dims(2)
      do i=1,dims(2)
        istart_t(1) = 1
        iend_t(1) = dims(1)
        istart_t(2) = i
        iend_t(2) = 1
        rcode = NF_GET_VARA_TEXT(fid,id_time,istart_t,iend_t,times(i))
      enddo

      itimes = n_times
      print_time = times(itimes)
      hour = print_time(12:13)
    
      print*, times(1)

      call close_file(fid)
    
      end subroutine get_hour 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    get the i,j,k of a variable
!    created July 2005  Zhiyong Meng  TAMU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_var_ij (file, plot_var, ix, iy, iz)
      implicit none
      character (LEN =  *)        :: file
      integer                      :: ix,iy,iz,i
      integer                      :: fid,rcode,id_var,dimlen
      integer :: dimids(NF_MAX_VAR_DIMS),dims(4)
      integer ndims, natts, ivtype
      character (len=*) :: plot_var
      character :: var_name(NF_MAX_NAME)
      integer istart(4), iend(4)

      call open_file(file, nf_nowrite, fid)
      write(*,*)' file =', file, ' fid =', fid

      rcode = nf_inq_varid(fid, plot_var, id_var)
      write(*,*)' plot_var =', plot_var, ' id_var=', id_var

      rcode = nf_inq_var(fid,id_var,var_name,ivtype,ndims,dimids,natts)
      write(*,*)' ndims =', ndims

      do i=1,ndims
         rcode = nf_inq_dimlen( fid, dimids(i), dims(i) )
      enddo

      istart  = 1
      iend    = 1
      do i = 1,ndims-1
         iend(i) = dims(i)
      enddo

      ix=iend(1)
      iy=iend(2)
      iz=iend(3)

      call close_file(fid)
 
      end subroutine get_var_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE NETCDF
