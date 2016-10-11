!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program replace

! 1, in shell script (out of this program), copy gfs_initial file to each output members (fort.900??);
! 2, link gfs_initial as fort.70010; and link ensemble members perturbed by 3dvar as fort.800??;
! 3, for each variable, 
!    -- calculate the average of perturbation xm, and perturbation xa=xa-xm;
!    -- replace xm by gfs initial xm=gfs;
!    -- calculate new variable value: xa=xm+xa; 
!    -- output to fort.900??.
  
!----------------------------------------------------------------
  use netcdf
  USE mpi_module
  implicit none

  integer            :: i_unit=80010, o_unit=90010, gfs_unit=70010
  integer            :: unit
  integer            :: numbers_en 
  integer            :: ix, jx, kx, i, j, k, len 
  integer            :: ii, jj, kk, m, n, fid, iargc, nmi,mstart, mend, ie
  character (len=10) :: wrf_file, filec, wrf_ifile
  character (len=5)  :: dum
  integer, parameter :: n_cdf = 1000
  real, dimension(n_cdf-1,2) :: erfinv

!----------------------------------------------------------------
! variables will be replaced
  character (len=10), allocatable, dimension(:) :: varname
  integer                                       :: varnum
  character (len=10)                            :: var

  real, allocatable, dimension(:,:,:,:) :: xa      !!! ensemble members
  real, allocatable, dimension(:,:,:)   :: xm      !!! ensemble mean
  real, allocatable, dimension(:,:,:)   :: gfs     !!! initial field
  real, allocatable, dimension(:,:,:)   :: dat3d, work,xq_n,xq_p
  real, allocatable, dimension(:,:)     :: cen_loc
  real, allocatable, dimension(:,:,:,:) :: xq_nsend,xq_psend


  call parallel_start()
!----------------------------------------------------------------
! variables will be replaced
  varnum = 23
  allocate( varname (varnum) )
!----------------------------------------------------------------

! Get ensemble number
  if (iargc() .lt. 1) then
     write(6,*) ' usage: replace_mean ensember_number '
     stop
  else
     call getarg(1,dum)
     read( dum, '(i5)' ) numbers_en
     write(*,*)'numbers_en: ', numbers_en
  endif 

  varname = (/'U         ', 'V         ', 'W         ', 'PH        ', 'P         ', &
              'T         ', 'MU        ', 'MUB       ', 'PHB       ', 'PB        ', &
              'Q2        ', 'T2        ', 'TH2       ', 'PSFC      ', &
              'U10       ', 'V10       ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', &
              'QICE      ', 'QSNOW     ', 'QGRAUP    ', 'TSK       '/)


!----------------------------------------------------------------
! Get WRF model dimension from the 1th member
  call get_var_ij ( 'fort.70010', 'T         ', ix, jx, kx )
  write(*,*)'wrf dimension: ', ix, jx, kx

!----------------------------------------------------------------
  if(mod(numbers_en,nprocs).eq.0) then
    nmi=numbers_en/nprocs
  else
    nmi=numbers_en/nprocs+1
  endif
  mstart=my_proc_id*nmi+1
  mend=(my_proc_id+1)*nmi

 do_wrf_var  : do m = 1, varnum
 var = varname(m)
 if(my_proc_id==0) write(*,*)var
!.... get dimensions
   call wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
   allocate( xa    ( ii, jj, kk, nmi ) )
   allocate( xm    ( ii, jj, kk ) )
   allocate( gfs   ( ii, jj, kk ) )
   allocate( dat3d ( ii, jj, kk ) )
   allocate( work  ( ii, jj, kk ) )
   allocate( xq_p  ( ii, jj, kk ) )
   allocate( xq_n  ( ii, jj, kk ) )
   allocate( xq_psend    ( ii, jj, kk, nmi ) )
   allocate( xq_nsend    ( ii, jj, kk, nmi ) )

!.... get gfs initial 
   if ( kk > 1 ) then
     call get_variable3d( 'fort.70010', var, ii, jj, kk, 1, gfs )
   else if ( kk == 1 ) then
     call get_variable2d( 'fort.70010', var, ii, jj, 1, gfs )
   endif

  work = 0.
  xa = 0.
  do_ensemble_member : do n=1,nmi   !mstart,mend
  ie = (n-1)*nprocs+my_proc_id+1
  !ie = nmi*my_proc_id+n 
  if (ie <= numbers_en) then
!.... get ensemble and calculate average
      write(wrf_file,'(a5,i5.5)')'fort.',i_unit+ie
!....... get data and sum
      if ( kk > 1 ) then
        call get_variable3d( wrf_file, var, ii, jj, kk, 1, dat3d )
      else if ( kk == 1 ) then
        call get_variable2d( wrf_file, var, ii, jj, 1, dat3d )
      endif
      xa(:,:,:,n)=dat3d(:,:,:)
  endif
  end do do_ensemble_member

  CALL MPI_Allreduce(sum(xa,4),work,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
  xm = work/float(numbers_en)
  if ( my_proc_id == 0 ) write(*,*)'nmi=',nmi,'xq-mean ',minval(xm),'~',maxval(xm)
  work = 0.

!.... replace xm by gfs
!.... output to ensemble
  do n=1,nmi !mstart,mend
  ie = (n-1)*nprocs+my_proc_id+1
  if (ie <= numbers_en) then
    work=xa(:,:,:,n)-xm+gfs
  endif
  end do

!!! Removing negative Q-value by Minamide 2015.5.26 
  if (var=='QCLOUD    ' .or. var=='QRAIN     ' .or. var=='QICE      ' .or. &
      var=='QGRAUP    ' .or. var=='QSNOW     ') then 
   xq_p = 0.
   xq_n = 0.
   xq_psend = 0.
   xq_nsend = 0.
   do n=1,nmi !mstart,mend
    ie = (n-1)*nprocs+my_proc_id+1
    if (ie <= numbers_en) then
     where(work >= 0.) xq_psend(:,:,:,n) = work
     where(work < 0.) xq_nsend(:,:,:,n) = work
    endif
   enddo
   call MPI_Allreduce(sum(xq_psend,4),xq_p,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
   call MPI_Allreduce(sum(xq_nsend,4),xq_n,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
   if( my_proc_id == 0 ) write(*,*)'original xq value',minval(work),'~',maxval(work)
   do n=1,nmi
    ie = (n-1)*nprocs+my_proc_id+1
    if(ie<=numbers_en) then
      where(work < 0.) work = 0.
      where(xq_p >= abs(xq_n).and.xq_p > 0.) work = work*(xq_p+xq_n)/xq_p
      where(xq_p < abs(xq_n).or. xq_p == 0.)  work = 0.
    endif
   enddo
   if( my_proc_id == 0 ) write(*,*)'non-negative xq value',minval(work),'~',maxval(work)
  endif

  do n=1,nmi !mstart,mend
  ie = (n-1)*nprocs+my_proc_id+1
  if (ie <= numbers_en) then
    write(wrf_file,'(a5,i5.5)')'fort.',o_unit+ie
    if(my_proc_id==0) write(*,*)'output to ', wrf_file
    call open_file( wrf_file, nf_write, fid )

    if ( kk > 1 ) then
       call write_variable3d(fid, var, ii, jj, kk, 1, work )
    else if ( kk == 1 ) then
       call write_variable2d(fid, var, ii, jj,  1, work )
    endif
    call close_file( fid )
  endif
  end do 
 
  deallocate( xa    )
  deallocate( xm    )
  deallocate( gfs   )
  deallocate( dat3d )
  deallocate( work  )
  deallocate( xq_p  )
  deallocate( xq_n  )
  deallocate( xq_psend  )
  deallocate( xq_nsend  )
  call mpi_barrier(comm,ierr)

 end do do_wrf_var

  call parallel_finish()

  if(my_proc_id==0) write(*,*)'!!! Successful completion of replace_mean.exe !!!'

end program replace

!==============================================================================
subroutine wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )

   character(len=10), intent(in)   :: var
   integer, intent(in)             :: ix, jx, kx
   integer, intent(out)            :: ii, jj, kk

   ii = ix
   jj = jx
   kk = kx
   if      ( var == 'U         ' ) then
      ii = ix + 1
   else if ( var == 'V         ' ) then
      jj = jx + 1
   else if ( var == 'W         ' .or. var == 'PH        ' .or. var == 'PHB       ' ) then
      kk = kx + 1
   else if ( var == 'MU        ' .or. var == 'MUB       ' .or. var == 'Q2        '  &
        .or. var == 'T2        ' .or. var == 'TH2       ' .or. var == 'PSFC      '  &
        .or. var == 'SST       ' .or. var == 'TSK       ' .or. var == 'XICE      '  &
        .or. var == 'SFROFF    ' .or. var == 'UDROFF    ' .or. var == 'IVGTYP    '  &
        .or. var == 'ISLTYP    ' .or. var == 'VEGFRA    ' .or. var == 'GRDFLX    '  &
        .or. var == 'SNOW      ' .or. var == 'SNOWH     ' .or. var == 'CANWAT    '  &
        .or. var == 'SST       ' .or. var == 'MAPFAC_M  ' .or. var == 'F         '  &
        .or. var == 'E         ' .or. var == 'SINALPHA  ' .or. var == 'COSALPHA  '  &
        .or. var == 'HGT       ' .or. var == 'TSK       ' .or. var == 'RAINC     '  &
        .or. var == 'RAINNC    ' .or. var == 'SWDOWN    ' .or. var == 'GLW       '  &
        .or. var == 'XLAT      ' .or. var == 'XLONG     ' .or. var == 'TMN       '  &
        .or. var == 'XLAND     ' .or. var == 'PBLH      ' .or. var == 'HFX       '  &
        .or. var == 'QFX       ' .or. var == 'LH        ' .or. var == 'SNOWC     '  &
        .or. var == 'SR        ' .or. var == 'POTEVP    ' .or. var == 'U10       '  &
        .or. var == 'V10       ' .or. var == 'MU0       ' ) then
      kk = 1
   else if ( var == 'MAPFAC_U  ' ) then
      kk = 1
      ii = ix + 1
   else if ( var == 'MAPFAC_V  ' ) then
      kk = 1
      jj = jx + 1
   else if ( var == 'FNM       ' .or. var == 'FNP       '  &
        .or. var == 'RDNW      ' .or. var == 'RDN       '  &
        .or. var == 'DNW       ' .or. var == 'DN        '  &
        .or. var == 'ZNU       '                          ) then
      ii = 1
      jj = 1
   else if ( var == 'ZNW       '                          ) then
      ii = 1
      jj = 1
      kk = kx + 1
   endif

end subroutine wrf_var_dimension

