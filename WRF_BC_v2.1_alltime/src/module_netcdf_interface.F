MODULE module_netcdf_interface

CONTAINS

!--------------------------------------------------------------------

  subroutine get_times_cdf( file, times, n_times, max_times, debug )
        
  implicit none

  include 'netcdf.inc'

  integer, intent(in)  ::  max_times  
  integer, intent(out) ::  n_times
  character (len=80), intent(in) :: file
  character (len=80), intent(out) :: times(max_times)
  logical, intent(in ) :: debug

  integer cdfid, rcode, id_time
  character (len=80) :: varnam, time1
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode == 0) then
    if(debug) write(6,*) ' open netcdf file ', trim(file)
  else
    write(6,*) ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_time = ncvid( cdfid, 'Times', rcode )

  rcode = nf_inq_var( cdfid, id_time, varnam, ivtype, ndims, dimids, natts )
  if(debug) then
    write(6,*) ' number of dims for Time ',ndims
  endif
  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  get the times
  
  n_times = idims(2)
  do i=1,idims(2)
    istart(1) = 1
    iend(1) = idims(1)
    istart(2) = i
    iend(2) = 1

    rcode = NF_GET_VARA_TEXT  ( cdfid, id_time,  &
                                istart, iend,    &
                                times(i)          )
    time1 = times(i)

    if(debug) write(6,*) trim(file), time1(1:19)
    times(i) = time1(1:19)
  enddo

  write(6,*) ' exiting get_times_cdf '

  call ncclos(cdfid,rcode)

  end subroutine get_times_cdf

!-------------------------------------------------------------------------------

  subroutine get_dims_cdf( file, var, dims, ndims, debug )
        
  implicit none

  include 'netcdf.inc'

  character (len=80), intent(in) :: file
  character (len=*), intent(in) :: var
  logical, intent(in ) :: debug
  integer, intent(out), dimension(4) :: dims
  integer, intent(out) :: ndims

  integer cdfid, rcode, id_time
  character (len=80) :: varnam, time1
  integer :: natts, istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )
  
  if( rcode == 0) then
    if(debug) write(6,*) ' open netcdf file ', trim(file)
  else
    write(6,*) ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_time = ncvid( cdfid, var, rcode )

  rcode = nf_inq_var( cdfid, id_time, varnam, ivtype, ndims, dimids, natts )
  if(debug) then
    write(6,*) ' number of dims for ',var,' ',ndims
  endif
  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), dims(i) )
    if(debug) write(6,*) ' dimension ',i,dims(i)
  enddo

  call ncclos(cdfid,rcode)

  end subroutine get_dims_cdf

!-------------------------------------------------------------------------------

  subroutine get_gl_att_int_cdf( file, att_name, value, debug )
        
  implicit none

  include 'netcdf.inc'

  character (len=80), intent(in) :: file
  character (len=*), intent(in) :: att_name
  logical, intent(in ) :: debug
  integer, intent(out) :: value

  integer :: cdfid, rcode

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode == 0) then
    if(debug) write(6,*) ' open netcdf file ', trim(file)
  else
    write(6,*) ' error openiing netcdf file ', trim(file)
    stop
  end if

  rcode = NF_GET_ATT_INT(cdfid, nf_global, att_name, value )

  call ncclos(cdfid,rcode)

  if(debug) write(6,*) ' global attribute ',att_name,' is ',value

  end subroutine get_gl_att_int_cdf

!-------------------------------------------------------------------------------

  subroutine get_gl_att_real_cdf( file, att_name, value, debug )
        
    implicit none

    include 'netcdf.inc'

    character (len=80), intent(in) :: file
    character (len=*), intent(in) :: att_name
    logical, intent(in ) :: debug
    real,    intent(out) :: value
    real(kind=8)         :: tmp

    integer :: cdfid, rcode, ivtype

    cdfid = ncopn(file, NCNOWRIT, rcode )

    if( rcode == 0) then
      if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
      write(6,*) ' error openiing netcdf file ', trim(file)
      stop
    end if

    rcode = NF_INQ_ATTTYPE( cdfid, nf_global, att_name, ivtype )

    write(unit=*, fmt='(a, i6)') &
         'ivtype:', ivtype, &
         'NF_REAL=', NF_REAL, &
         'NF_DOUBLE=', NF_DOUBLE, &
         'kind(value)=', kind(value)

    if((ivtype == NF_REAL) .and. (kind(value) == 4)) then
       rcode = NF_GET_ATT_REAL(cdfid, nf_global, att_name, value )
    else if((ivtype == NF_DOUBLE) .and. (kind(value) == 4)) then
       rcode = NF_GET_ATT_REAL(cdfid, nf_global, att_name, tmp )
       value = tmp
    else if((ivtype == NF_DOUBLE) .and. (kind(value) == 8)) then
       rcode = NF_GET_ATT_REAL(cdfid, nf_global, att_name, value )
    else
       write(unit=*, fmt='(a, i6)') &
            'Unrecognizable ivtype:', ivtype
       stop
    endif

    call ncclos(cdfid,rcode)

    if(debug) write(6,*) ' global attribute ',att_name,' is ',value

  end subroutine get_gl_att_real_cdf

!--------------------------------------------------------------------

  subroutine get_var_3d_real_cdf( file, var, data, &
                                  i1, i2, i3, time, debug )
        
  implicit none

  include 'netcdf.inc'

  integer, intent(in)  ::  i1, i2, i3, time
  character (len=80), intent(in) :: file
  logical, intent(in ) :: debug
  character (len=*), intent(in) :: var
  real, dimension(i1,i2,i3), intent(out) :: data
  real(kind=8), dimension(i1,i2,i3) :: tmp

  character (len=80) :: varnam, time1

  integer :: cdfid, rcode, id_data
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode /= 0) then
    write(6,*) ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_data = ncvid( cdfid, var, rcode )

  rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )

  if(debug) then
    write(unit=*, fmt='(3a,i6)') ' get_var_3d_real_cdf: dims for ',var,' ',ndims
    write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
    write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
    write(unit=*, fmt='(a,i6)') ' kind(data)=', kind(data)
  endif

  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) write(unit=*, fmt='(a,2i6)') ' dimension ',i,idims(i)
  enddo

!  check the dimensions

   if( (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (i3 /= idims(3)) .or.  &
       (time > idims(4))     )  then

     write(6,*) ' error in 3d_var_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) i3, idims(3)
     write(6,*) time, idims(4)
     write(6,*) ' error stop '
     stop

   end if

!  get the data
  
    istart(1) = 1
    iend(1) = i1
    istart(2) = 1
    iend(2) = i2
    istart(3) = 1
    iend(3) = i3
    istart(4) = time
    iend(4) = 1

    if((ivtype == NF_REAL) .and. (kind(data) == 4)) then
       call ncvgt( cdfid,id_data,istart,iend,data,rcode)
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 4)) then
       call ncvgt( cdfid,id_data,istart,iend,tmp,rcode)
       data = tmp
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 8)) then
       call ncvgt( cdfid,id_data,istart,iend,data,rcode)
    else
       write(unit=*, fmt='(a, i6)') &
            'Unrecognizable ivtype:', ivtype
       stop
    endif

    if(debug) then
       write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1,1)
    endif

    call ncclos(cdfid,rcode)

  end subroutine get_var_3d_real_cdf

!--------------------------------------------------------------------

  subroutine get_var_2d_real_cdf( file, var, data, &
                                  i1, i2, time, debug )
        
  implicit none

  include 'netcdf.inc'

  integer, intent(in)  ::  i1, i2, time
  character (len=80), intent(in) :: file
  logical, intent(in ) :: debug
  character (len=*), intent(in) :: var
  real, dimension(i1,i2), intent(out) :: data
  real(kind=8), dimension(i1,i2) :: tmp

  integer cdfid, rcode, id_data
  character (len=80) :: varnam, time1
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode /= 0) then
    write(unit=*, fmt='(2a)') ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_data = ncvid( cdfid, var, rcode )

  rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )

  if(debug) then
    write(unit=*, fmt='(3a,i6)') ' get_var_2d_real_cdf: dims for ',var,' ',ndims
  endif

  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) then
      write(unit=*, fmt='(a,2i6)') ' dimension ',i,idims(i)
      write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
      write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
    endif
  enddo

!  check the dimensions

   if( (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (time > idims(3))     )  then

     write(6,*) ' error in 2d_var_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) time, idims(4)
     write(6,*) ' error stop '
     stop

   end if

!  get the data
  
    istart(1) = 1
    iend(1) = i1
    istart(2) = 1
    iend(2) = i2
    istart(3) = time
    iend(3) = 1

    if((ivtype == NF_REAL) .and. (kind(data) == 4)) then
       call ncvgt( cdfid,id_data,istart,iend,data,rcode)
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 8)) then
       call ncvgt( cdfid,id_data,istart,iend,data,rcode)
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 4)) then
       call ncvgt( cdfid,id_data,istart,iend,tmp,rcode)
       data = tmp
    else
       write(unit=*, fmt='(a, i6)') &
            'Unrecognizable ivtype:', ivtype
       stop
    endif

    if(debug) then
       write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1)
    endif

    call ncclos(cdfid,rcode)

  end subroutine get_var_2d_real_cdf

!--------------------------------------------------------------------

  subroutine put_var_3d_real_cdf( file, var, data, &
                                  i1, i2, i3, time, debug )
        
  implicit none

  include 'netcdf.inc'

  integer, intent(in)  ::  i1, i2, i3, time
  character (len=80), intent(in) :: file
  logical, intent(in ) :: debug
  character (len=*), intent(in) :: var
  real, dimension(i1,i2,i3), intent(in) :: data
  real(kind=8), dimension(i1,i2,i3) :: tmp

  integer cdfid, rcode, id_data
  character (len=80) :: varnam, time1
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCWRITE, rcode )

  if( rcode /= 0) then
    write(unit=*, fmt='(2a)') ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_data = ncvid( cdfid, var, rcode )

  rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )

  if(debug) then
    write(unit=*, fmt='(3a,i6)') ' put_var_3d_real_cdf: dims for ',var,' ',ndims
  endif

  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  check the dimensions

   if( (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (i3 /= idims(3)) .or.  &
       (time > idims(4))     )  then

     write(6,*) ' error in 3d_var_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) i3, idims(3)
     write(6,*) time, idims(4)
     write(6,*) ' error stop '
     stop

   end if

!  get the data
  
    istart(1) = 1
    iend(1) = i1
    istart(2) = 1
    iend(2) = i2
    istart(3) = 1
    iend(3) = i3
    istart(4) = time
    iend(4) = 1

    if((ivtype == NF_REAL) .and. (kind(data) == 4)) then
       call ncvpt( cdfid,id_data,istart,iend,data,rcode)
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 8)) then
       tmp = data
       call ncvpt( cdfid,id_data,istart,iend,tmp,rcode)
    else if((ivtype == NF_DOUBLE) .and. (kind(data) == 4)) then
       tmp = data
       call ncvpt( cdfid,id_data,istart,iend,tmp,rcode)
    else
       write(unit=*, fmt='(a, i6)') &
            'Unrecognizable ivtype:', ivtype
       stop
    endif

    call ncclos(cdfid,rcode)

  end subroutine put_var_3d_real_cdf

!--------------------------------------------------------------------

  subroutine put_var_2d_real_cdf( file, var, data, &
                                  i1, i2, time, debug )
        
    implicit none

    include 'netcdf.inc'

    integer, intent(in)  ::  i1, i2, time
    character (len=80), intent(in) :: file
    logical, intent(in ) :: debug
    character (len=*), intent(in) :: var
    real, dimension(i1,i2), intent(in) :: data
    real(kind=8), dimension(i1,i2) :: tmp

    integer :: cdfid, rcode, id_data
    character (len=80) :: varnam, time1
    integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
    integer :: i, ivtype

    cdfid = ncopn(file, NCWRITE, rcode )

    if( rcode == 0) then
      if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
      write(6,*) ' error openiing netcdf file ', trim(file)
      stop
    end if

    id_data = ncvid( cdfid, var, rcode )

    rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )
    if(debug) then
      write(6,*) ' number of dims for ',var,' ',ndims
    endif
    do i=1,ndims
      rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
      if(debug) write(6,*) ' dimension ',i,idims(i)
    enddo

!---check the dimensions

    if((i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (time > idims(3))     )  then

       write(6,*) ' error in 3d_var_real read, dimension problem '
       write(6,*) i1, idims(1)
       write(6,*) i2, idims(2)
       write(6,*) time, idims(3)
       write(6,*) ' error stop '
       stop
     end if

!----get the data
  
     istart(1) = 1
     iend(1) = i1
     istart(2) = 1
     iend(2) = i2
     istart(3) = time
     iend(3) = 1

     if((ivtype == NF_REAL) .and. (kind(data) == 4)) then
        call ncvpt( cdfid,id_data,istart,iend,data,rcode)
     else if((ivtype == NF_DOUBLE) .and. (kind(data) == 8)) then
        tmp = data
        call ncvpt( cdfid,id_data,istart,iend,tmp,rcode)
     else if((ivtype == NF_DOUBLE) .and. (kind(data) == 4)) then
        tmp = data
        call ncvpt( cdfid,id_data,istart,iend,tmp,rcode)
     else
        write(unit=*, fmt='(a, i6)') &
            'Unrecognizable ivtype:', ivtype
        stop
     endif

     call ncclos(cdfid,rcode)

  end subroutine put_var_2d_real_cdf

!--------------------------------------------------------------------

  subroutine get_var_2d_int_cdf(file, var, data, i1, i2, time, debug)
        
  implicit none

  include 'netcdf.inc'

  integer, intent(in)  ::  i1, i2, time
  character (len=80), intent(in) :: file
  logical, intent(in ) :: debug
  character (len=*), intent(in) :: var
  integer, dimension(i1,i2), intent(out) :: data

  integer cdfid, rcode, id_data
  character (len=80) :: varnam, time1
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode /= 0) then
    write(unit=*, fmt='(2a)') ' error openiing netcdf file ', trim(file)
    stop
  end if

  id_data = ncvid( cdfid, var, rcode )

  rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )

  if(debug) then
    write(unit=*, fmt='(3a,i6)') ' get_var_2d_real_cdf: dims for ',var,' ',ndims
  endif

  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) then
      write(unit=*, fmt='(a,2i6)') ' dimension ',i,idims(i)
      write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
      write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
    endif
  enddo

!  check the dimensions

   if( (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (time > idims(3))     )  then

     write(6,*) ' error in 2d_var_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) time, idims(4)
     write(6,*) ' error stop '
     stop

   end if

!  get the data
  
    istart(1) = 1
    iend(1) = i1
    istart(2) = 1
    iend(2) = i2
    istart(3) = time
    iend(3) = 1

    call ncvgt( cdfid,id_data,istart,iend,data,rcode)

    if(debug) then
       write(unit=*, fmt='(a, i8)') ' Sample data=', data(1,1)
    endif

    call ncclos(cdfid,rcode)

  end subroutine get_var_2d_int_cdf

!-------------------------------------------------------------------------------

  subroutine get_att_cdf( file, var, debug )

     implicit none

     include 'netcdf.inc'

     character (len=80), intent(in) :: file
     character (len=*), intent(in) :: var
     logical, intent(in ) :: debug

     integer :: cdfid, status, varid, n, natts
     character (len=256) :: att_name

     status = NF_OPEN(file, NF_NOWRITE, cdfid )

     status = NF_INQ_VARID( cdfid, var, varid )

     if( status == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
     else
       write(6,*) ' error openiing netcdf file ', trim(file)
       stop
     end if

     status = NF_INQ_VARNATTS(cdfid, varid, natts )

     do n=1, natts
        status = NF_INQ_ATTNAME(cdfid, varid, n, att_name )
   
        write(unit=*, fmt='(a,i2,2a)') &
          'att_name(',n,')=', trim(att_name)
     enddo

     status = NF_CLOSE(cdfid)

  end subroutine get_att_cdf

!-------------------------------------------------------------------------------

  subroutine put_att_cdf( file, var, att_name, text, debug )

     implicit none

     include 'netcdf.inc'

     character (len=80), intent(in) :: file
     character (len=*), intent(in) :: var, att_name, text
     logical, intent(in ) :: debug

     integer :: cdfid, status, varid, n, natts
     character (len=256) :: loc_att_name

     status = NF_OPEN(file, NF_WRITE, cdfid )

     if( status == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
     else
       write(6,*) ' error openiing netcdf file ', trim(file)
       stop
     end if

     status = NF_INQ_VARID( cdfid, var, varid )

     status = NF_INQ_VARNATTS(cdfid, varid, natts )

     do n=1, natts
        status = NF_INQ_ATTNAME(cdfid, varid, n, loc_att_name )

        write(unit=*, fmt='(a,i2,2a)') &
          'loc_att_name(',n,')=', trim(loc_att_name)

        if(trim(loc_att_name) == trim(att_name)) then
           write(unit=*, fmt='(2a)') &
             'att_name=', trim(att_name)

           status = NF_PUT_ATT_TEXT(cdfid, varid, trim(att_name), len(text), trim(text))

           if(status == 0) then
              if(debug) then
                 write(unit=*, fmt='(4a)') &
                      'write ', trim(att_name), 'to netcdf file ', trim(file)
              endif
           else
              write(unit=*, fmt='(a, i8)') &
                     'Status= ', status

              write(unit=*, fmt='(4a)') &
                     'Failed to write ', trim(att_name), 'to netcdf file ', trim(file)

!             if(status /= NF_NOERR) call handle_err(status)
              stop
           endif

           exit
        endif
     enddo

     status = NF_CLOSE(cdfid)

  end subroutine put_att_cdf

  subroutine get_bdyfrq(this_datestr, next_datestr, bdyfrq, debug )
        
  implicit none

  include 'netcdf.inc'

  character(len=80), intent(in ) :: this_datestr, next_datestr
  real,              intent(out) :: bdyfrq
  logical,           intent(in ) :: debug

  real    :: this_bdy_time, next_bdy_time

  CALL atotime(this_datestr, this_bdy_time)
  CALL atotime(next_datestr, next_bdy_time)

  bdyfrq = next_bdy_time - this_bdy_time

  if(debug) then
     write(unit=*, fmt='(a, f20.1)') &
          'next_bdy_time=', next_bdy_time, &
          'this_bdy_time=', this_bdy_time, &
          'bdyfrq       =', bdyfrq
  endif

  end subroutine get_bdyfrq

  subroutine get_bdytimestr_cdf(file, time_flag, bdytimestr, max_times, debug)
        
  implicit none

  include 'netcdf.inc'

  integer,           intent(in ) :: max_times
  logical,           intent(in ) :: debug
  character(len=80), intent(in ) :: file
  character(len=*),  intent(in ) :: time_flag
  character(len=80), intent(out) :: bdytimestr(max_times)

  character(len=80) :: varnam, time1

  integer :: cdfid, rcode, id_time
  integer :: ndims, natts, idims(10), istart(10),iend(10), dimids(10)
  integer :: i, ivtype

  cdfid = ncopn(file, NCNOWRIT, rcode )

  if( rcode == 0) then
    if(debug) write(6,*) ' open netcdf file ', trim(file)
  else
    write(6,*) ' error openiing netcdf file ', trim(file)
    stop
  end if

  varnam = 'md___' // trim(time_flag) // 'e_x_t_d_o_m_a_i_n_m_e_t_a_data_'

  id_time = ncvid( cdfid, varnam, rcode )

  rcode = nf_inq_var( cdfid, id_time, varnam, ivtype, ndims, dimids, natts )

  if(debug) then
    write(6,*) ' number of dims for bdytime ',ndims
  endif

  do i=1,ndims
    rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
    if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  get the bdytime
  
  do i=1,idims(2)
    istart(1) = 1
    iend(1) = idims(1)
    istart(2) = i
    iend(2) = 1

    rcode = NF_GET_VARA_TEXT  ( cdfid, id_time,  &
                                istart, iend,    &
                                bdytimestr(i) )
    time1 = bdytimestr(i)

    if(debug) write(6,*) trim(file), time1(1:19)
    bdytimestr(i) = time1(1:19)
  enddo

  call ncclos(cdfid,rcode)

  end subroutine get_bdytimestr_cdf

  SUBROUTINE atotime(date_char, st)

   implicit none

   character(len=80), intent(in ) :: date_char
   real,              intent(out) :: st
  
   integer                :: ccyy,mo,dd,hh,mi,ss,i
   integer, dimension(12) :: mmday

   character(len=1)       :: dash, underscore, colon
  
   mmday=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  
   read(date_char(1:19),'(i4,a,4(i2,a),i2)') &
        ccyy, dash, &
          mo, dash, &
          dd, underscore, &
          hh, colon, &
          mi, colon, &
          ss
  
   if(mod(ccyy,4) == 0) then
      mmday(2) = 29
      if(mod(ccyy,400) == 0) then
         mmday(2) = 29
      else if (mod(ccyy,100) == 0) then
         mmday(2) = 28
      endif
   endif
  
   dd=dd+365*(ccyy-2000)

   do i=1,mo-1
      dd=dd+mmday(i)
   enddo
  
   st = real(ss) &
      + 60.0*(real(mi) &
      + 60.0*(real(hh) &
      + 24.0* real(dd)))

END SUBROUTINE atotime

END MODULE module_netcdf_interface

