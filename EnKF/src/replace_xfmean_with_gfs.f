!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program replace_xfmean_with_gfs

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
  character (len=4)  :: numbers_enc

  character (len=80) :: times
  character (len=26) :: finish_file_flg
!----------------------------------------------------------------
! variables will be replaced
  character (len=10), allocatable, dimension(:) :: varname
  integer                                       :: varnum
  character (len=10)                            :: var

  real, allocatable, dimension(:,:,:)   :: xa      !!! ensemble members
  real, allocatable, dimension(:,:,:)   :: xm      !!! ensemble mean
  real, allocatable, dimension(:,:,:)   :: gfsm    !!! initial field
  real, allocatable, dimension(:,:,:)   :: work

  character(len=24)    :: btkfilename
  character(len= 8)    :: cstorm
  character(len=10)    :: cdate
  integer              :: year, month, day, hour
  real                 :: tc_ix, tc_jy
  real, dimension(4)   :: tc_center    !lat, lon, pmin, vmax

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
  if (iargc() .lt. 2) then
     write(6,*) ' usage: replace_xfmean_with_gfs.exe stormid ensember_number'
     stop
  else
     call getarg(1, cstorm)
     call getarg(2, numbers_enc)
     read(numbers_enc,'(i)')numbers_en
  endif 

  write(*,*)'numbers_en =', numbers_en

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
! Get btk file (Best Track File)
  cdate(1:10)=times(1:4)//times(6:7)//times(9:10)//times(12:13)
  read(cdate,'(i4,3i2)')year, month, day, hour
  btkfilename=cdate(1:10)//'.b'//cstorm//'.dat'
  call rd_interp_besttrack(cdate, btkfilename, tc_center)

  ! 3, calculate storm center's position according to wrf domain grid
  call latlon_to_ij( proj, tc_center(1), tc_center(2), tc_ix, tc_jy)
  write(*,*)'TC position=', tc_ix, tc_jy
  Rmax=800.*1000.   !200km
  Rmin=200.*1000.   !100km
  Rmax_i = Rmax/proj%dx
  Rmin_i = Rmin/proj%dx
  rep_is = int(tc_ix - Rmax_i)
  rep_ie = int(tc_ix + Rmax_i)
  rep_js = int(tc_jy - Rmax_i)
  rep_je = int(tc_jy + Rmax_i)
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

         ! r>800km xg=xf'+gfs
         work(:,:,:)=(xa(:,:,:)-xm(:,:,:))+gfsm(:,:,:)
         ! r<200km xg=xa
         !work(int(tc_ix-Rmin_i):int(tc_ix+Rmin_i), int(tc_jy-Rmin_i):int(tc_jy+Rmin_i), :) = xa(int(tc_ix-Rmin_i):int(tc_ix+Rmin_i), int(tc_jy-Rmin_i):int(tc_jy+Rmin_i), :)
         
         do k = 1, kk
            do j = rep_js, rep_je
            do i = rep_is, rep_ie
               r = sqrt((i-tc_ix)**2 + (j-tc_jy)**2)
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

      deallocate( xa    )
      deallocate( xm    )
      deallocate( gfsm  )
      deallocate( work  )

  end do do_wrf_var

  finish_file_flg = times(1:4)//times(6:7)//times(9:10)//times(12:13)//times(15:16)//'.finish_mixing'
  do k = 1, 12
     if ( finish_file_flg(k:k) .eq. " " )finish_file_flg(k:k) ='0'
  enddo
  open(10,file=finish_file_flg)
  write(10,*)times
  !write(10,*)'xa = xm + 0.8*(xa-xm) + 0.2*(gfs-gfsm)'
  !write(10,*)'work(:,:,:)=xa(:,:,:,n)+(gfs(:,:,:,n)-gfsm(:,:,:))*0.15' 
  write(10,*)'xg = (xf-xfmean) + w*xfmean +(1-w)*xgfs'
  close(10)

end program replace_xfmean_with_gfs

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

!==============================================================================
!=======================================================================================
   subroutine rd_interp_besttrack ( times, filename, center )

   implicit none
   character (len=80), intent(in)        :: times
   character (len=*), intent(in)         :: filename
   real, dimension(4), intent(out)       :: center

   integer            :: length, iost, hours_from_date
   integer, dimension(3) :: year, month, day, hour, dhour     !1=current, 2=pre, 3=post
   real, dimension(3) :: lat, lon, wsp, slp
   integer            :: iyear, imonth, iday, ihour, ilat, ilon, iwsp, islp, diff_hour
   character (len=1)  :: clat, clon

   dhour(1:3)=-999999
   read(times, '(i4, 3i2)')year(1), month(1), day(1), hour(1)
   dhour(1) = hours_from_date ( year(1), month(1), day(1), hour(1))

   length = len_trim(filename)
   open(39, file=filename(1:length), status='old', form = 'formatted', iostat = iost )
   if( iost .ne. 0 ) then
     center = -9999999.0
     write(*,*)'CANNOT find best-track file!!!!!!'
     return
   end if

   do_get_track_loop : do
      read(39, '(8x, i4,3i2, 16x, i4, a1, 1x, i5, a1, 1x, i4, 1x, i5)', iostat = iost)  &
              iyear, imonth, iday, ihour, ilat, clat, ilon, clon, iwsp, islp
      if( iost .ne. 0 ) exit
      !write(*,'(i4,3i2,i4, a1, 1x, i5, a1, 1x, i4, 1x, i5)')iyear, imonth, iday, ihour, ilat, clat, ilon, clon, iwsp, islp
      if (clat == 'S') ilat = -ilat
      if (clon == 'W') ilon = -ilon
      diff_hour = hours_from_date ( iyear, imonth, iday, ihour)
      if (diff_hour == dhour(1) ) then
         lat(2) = ilat/10.
         lon(2) = ilon/10.
         wsp(2) = iwsp*0.514444
         slp(2) = islp*1.0
         dhour(2) = 99999
         lat(3) = ilat/10.
         lon(3) = ilon/10.
         wsp(3) = iwsp*0.514444
         slp(3) = islp*1.0
         dhour(3) = 99999
         exit
      else if (diff_hour < dhour(1) ) then
         year(2) = iyear
         month(2) = imonth
         day(2) = iday
         hour(2) = ihour
         lat(2) = ilat/10.
         lon(2) = ilon/10.
         wsp(2) = iwsp*0.514444
         slp(2) = islp*1.0
         dhour(2) = dhour(1) - diff_hour
      else
         year(3) = iyear
         month(3) = imonth
         day(3) = iday
         hour(3) = ihour
         lat(3) = ilat/10.
         lon(3) = ilon/10.
         wsp(3) = iwsp*0.514444
         slp(3) = islp*1.0
         dhour(3) = diff_hour - dhour(1)
         exit
      end if
   end do do_get_track_loop
   close(39)
   if (dhour(2) < -9999 .or. dhour(3) < -9999 ) then
      center(1:3) = -999999.0
   else
      center(1) = lat(2)+(lat(3)-lat(2))*abs(dhour(2))/(abs(dhour(2))+abs(dhour(3)))
      center(2) = lon(2)+(lon(3)-lon(2))*abs(dhour(2))/(abs(dhour(2))+abs(dhour(3)))
      center(3) = slp(2)+(slp(3)-slp(2))*abs(dhour(2))/(abs(dhour(2))+abs(dhour(3)))
      center(4) = wsp(2)+(wsp(3)-wsp(2))*abs(dhour(2))/(abs(dhour(2))+abs(dhour(3)))
   end if
   write(*,*)'hurricane center from obs is:',center
   return

   end subroutine rd_interp_besttrack
!==============================================================================
   function hours_from_date(iyear, imonth, iday, ihour)
   integer :: iyear, imonth, iday, ihour
   integer :: hours_from_date
   integer :: y, m, d, h

   m = mod(imonth+9, 12)
   y = iyear - int(m/10)
   hours_from_date = 365*y + int(y/4) + int(y/100) + int(y/400) + int((m*306+5)/10) + iday - 1
   hours_from_date = hours_from_date*24 + ihour
   return
   end function hours_from_date

!=======================================================================================
