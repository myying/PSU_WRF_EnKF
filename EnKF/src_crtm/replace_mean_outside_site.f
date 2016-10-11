!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program replace_mean

! 1, in shell script (out of this program), link enkf members (fort.900??);
! 2, link wrfinput_d0? as fort.70010;
! 3, link best track file: YYYYMMDDHH.bal092011.dat
! 4, for each variable,  xa = xf + w*xm + (1-w)*gfsm
  
!----------------------------------------------------------------
  use netcdf
  use mapinfo_define
  use map_utils
  use wrf_tools

  implicit none
  type(proj_info)    :: proj

  integer            :: i_unit=80010, o_unit=90010, gfs_unit=70010
  integer            :: unit
  integer            :: numbers_en 
  integer            :: ix, jx, kx, i, j, k, len 
  integer            :: ii, jj, kk, m, n, fid, iargc
  character (len=10) :: wrf_file, filec
  character (len=10)  :: dum
  character (len=80) :: times
!----------------------------------------------------------------
! variables will be replaced
  character (len=10), allocatable, dimension(:) :: varname
  integer                                       :: varnum
  character (len=10)                            :: var

  real, allocatable, dimension(:,:,:)   :: xa      !!! ensemble members
  real, allocatable, dimension(:,:,:)   :: xm      !!! ensemble mean
  real, allocatable, dimension(:,:,:)   :: gfsm    !!! initial field
  real, allocatable, dimension(:,:,:)   :: work

  real                 :: si,sj,slat,slon

  real               :: Rmax, Rmin, dx_rate, Rmax_i, Rmin_i
  integer            :: d_ii, d_jj
  integer            :: rep_is, rep_ie, rep_js, rep_je
  real               :: r, alpha

!----------------------------------------------------------------
! variables will be replaced
  varnum = 17
  allocate( varname (varnum) )
!----------------------------------------------------------------
! Get ensemble number
  if (iargc() .lt. 3) then
     write(6,*) ' usage: replace_mean_outside_site.exe site_lat site_lon numbers_en'
     stop
  else
     call getarg(1, dum)
     read(dum,'(f7.2)')slat
     call getarg(2, dum)
     read(dum,'(f7.2)')slon
     call getarg(3, dum)
     read(dum,'(i4)')numbers_en
  endif 

  write(*,*)'numbers_en =', numbers_en
  write(*,*)'site lat/lon =', slat, slon

  varname = (/'U         ', 'V         ', 'W         ', 'PH        ', &
              'T         ', 'MU        ', 'MUB       ',  &
              'Q2        ', 'T2        ', 'TH2       ', 'PSFC      ', &
              'U10       ', 'V10       ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', &
              'TSK       '/)

!----------------------------------------------------------------
! Get WRF model dimension from the 1th member
  call set_domain_proj('fort.70010', proj )
!  call get_var_ij ( 'fort.70010', 'T         ', ix, jx, kx )
  call get_ij ( 'fort.70010', ix, jx, kx )
  call get_times ( 'fort.70010', 'T         ', times )

  write(*,*)'wrf dimension: ', ix, jx, kx

!----------------------------------------------------------------
  ! 3, calculate storm center's position according to wrf domain grid
  call latlon_to_ij( proj, slat,slon, si, sj) 
  Rmax=300.*1000. 
  Rmin=200.*1000.  
  Rmax_i = Rmax/proj%dx
  Rmin_i = Rmin/proj%dx
  rep_is = int(si - Rmax_i)
  rep_ie = int(si + Rmax_i)
  rep_js = int(sj - Rmax_i)
  rep_je = int(sj + Rmax_i)
  if ( rep_is .le. 1 ) rep_is=1
  if ( rep_ie .ge. ix) rep_ie=ix
  if ( rep_js .le. 1 ) rep_js=1
  if ( rep_je .ge. jx) rep_je=jx

  write(*,*)'Rmax_i, Rmin_i =',Rmax_i, Rmin_i
  write(*,*)'rep_is, rep_ie, rep_js, rep_je =', rep_is, rep_ie, rep_js, rep_je

!----------------------------------------------------------------
  do_wrf_var  : do m = 1, varnum

!.... get dimensions
      var = varname(m)
      write(*,*)var
      call wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
      allocate( xa    ( ii, jj, kk ) )
      allocate( xm    ( ii, jj, kk ) )
      allocate( gfsm  ( ii, jj, kk ) )
      allocate( work  ( ii, jj, kk ) )

!.... get gfs initial 
      write(wrf_file,'(a5,i5.5)')'fort.',gfs_unit
!....... get data and sum
      if ( kk > 1 ) then
         call get_variable3d( wrf_file, var, ii, jj, kk, 1, gfsm )
      else if ( kk == 1 ) then
         call get_variable2d( wrf_file, var, ii, jj, 1, gfsm )
      endif

!.... get ensemble and calculate average
      xm = 0.0
      do_ensemble_member : do n = 1, numbers_en
         write(wrf_file,'(a5,i5.5)')'fort.',i_unit+n
!....... get data and sum
         if ( kk > 1 ) then
            call get_variable3d( wrf_file, var, ii, jj, kk, 1, xa )
         else if ( kk == 1 ) then
            call get_variable2d( wrf_file, var, ii, jj, 1, xa )
         endif
         xm = xm + xa
      end do do_ensemble_member
      xm = xm/float(numbers_en)

!.... mixing and output to ensemble
      do n = 1, numbers_en

         write(wrf_file,'(a5,i5.5)')'fort.',i_unit+n
         if ( kk > 1 ) then
            call get_variable3d( wrf_file, var, ii, jj, kk, 1, xa )
         else if ( kk == 1 ) then
            call get_variable2d( wrf_file, var, ii, jj, 1, xa )
         endif 

         write(wrf_file,'(a5,i5.5)')'fort.',o_unit+n
         write(*,*)'output to ', wrf_file
         call open_file( wrf_file, nf_write, fid )

         work(:,:,:)=(xa(:,:,:)-xm(:,:,:))+gfsm(:,:,:)
         do k = 1, kk
            do j = rep_js, rep_je
            do i = rep_is, rep_ie
               r = sqrt((i-si)**2 + (j-sj)**2)
               if ( r < Rmin_i ) then
                  alpha = 1.0
               else if ( r > Rmax_i ) then
                  alpha = 0.0
               else if ( r > Rmin_i .and. r <= Rmax_i ) then
                  !alpha = ((Rmax_i-Rmin_i)**2 - (r-Rmin_i)**2)/((Rmax_i-Rmin_i)**2 + (r-Rmin_i)**2)
                  alpha = (Rmax_i-r)/(Rmax_i-Rmin_i)
               endif

               work(i,j,k) = (xa(i,j,k)-xm(i,j,k)) + alpha*xm(i,j,k) + (1.0-alpha)*gfsm(i,j,k)
            end do
            end do
         end do

         if ( kk > 1 ) then
            call write_variable3d(fid, var, ii, jj, kk, 1, work )
         else if ( kk == 1 ) then
            call write_variable2d(fid, var, ii, jj,  1, work )
         endif
         call close_file( fid )
      enddo 

      write(wrf_file,'(a5,i5.5)')'fort.',o_unit+numbers_en+1
      write(*,*)'output to ', wrf_file
      call open_file( wrf_file, nf_write, fid )
      work=gfsm
      do k = 1, kk
         do j = rep_js, rep_je
         do i = rep_is, rep_ie
            r = sqrt((i-si)**2 + (j-sj)**2)
            if ( r < Rmin_i ) then
               alpha = 1.0
            else if ( r > Rmax_i ) then
               alpha = 0.0
            else if ( r > Rmin_i .and. r <= Rmax_i ) then
               alpha = (Rmax_i-r)/(Rmax_i-Rmin_i)
            endif
            work(i,j,k) = alpha*xm(i,j,k) +(1.0-alpha)*gfsm(i,j,k)
         end do
         end do
      end do
      if ( kk > 1 ) then
         call write_variable3d(fid, var, ii, jj, kk, 1, work)
      else if ( kk == 1 ) then
         call write_variable2d(fid, var, ii, jj,  1, work)
      endif
      call close_file( fid )

      deallocate( xa    )
      deallocate( xm    )
      deallocate( gfsm  )
      deallocate( work  )

  end do do_wrf_var
  write(*,*)'!!! Successful completion of replace_mean.exe !!!'

end program replace_mean

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

