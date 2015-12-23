!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  program so_to_3dvar

  
!----------------------------------------------------------------
! Definition of WRF model structure
! Definition of EnKF-specific variables
! Definition of Radar Observation structure
! definition of namelist
!  use constants
  use namelist_define
  use mapinfo_define
  use mpi_module
  use netcdf
  use obs_define
  use map_utils
  use wrf_tools
  use radar
!----------------------------------------------------------------
  implicit none
  
  integer            :: e_unit=70010     ! pure ensemble forecast
  integer            :: i_unit=80010     ! enkf forecast
  integer            :: o_unit=90010     ! enkf update
  integer            :: unit
  integer            :: ie, ix, jx, kx, ntot, itot, i, j, k, len, nr
  integer            :: ii, jj, kk, m, num_var
  character (len=10) :: wrf_file, filec
  real, allocatable, dimension(:,:) :: xa0      !!! ensemble members 
  real, allocatable, dimension(:,:) :: xm0      !!! ensemble mean

  integer, allocatable, dimension(:) :: var_start  ! every variables start point in xa0

!----------------------------------------------------------------
! definition of WRF-specific variables
!----------------------------------------------------------------
  type(proj_info)              :: proj
  real                         :: ptop
  character (len=80)           :: times, char_total_Radar, char_ned
!----------------------------------------------------------------
! Define observation variables
  integer                          :: iobs, num_rv      
  real, allocatable, dimension(:)  :: ob, obs_ii, obs_jj, obs_kk, obs_hh, obs_az, obs_ra, obs_el, obs_lat, obs_lon  

! just for radar
  real, allocatable, dimension(:)  :: radar_ii, radar_jj, radar_elv

  integer   :: iob, iunit, lev, nu, i1, j1
  real      :: di, di1, dj, dj1
  real, allocatable, dimension(:, :)    :: xlat, xlon, work2d
  integer   :: isign, length_sign, length_outfile
  character (len=120) :: outfile, sign
  integer, external :: iargc

!----------------------------------------------------------------
! Initialize parallel stuff
!----------------------------------------------------------------
  call parallel_start()

!----------------------------------------------------------------
! Get argc 
  if (iargc() .lt. 2) then
     write(6,*) ' usage: verify_v isign verification_point_vr_output_file'
     write(6,*) ' isign: 1, assimilation data(start from record 1 of SO)'
     write(6,*) '        2, verification data(start from record 2 of SO)'
     write(6,*) '        elae,  all so data'
     write(6,*) '  '
     write(6,*) ' example: verify.exe 2 verify_ahgx_vhgxsr1_nx10z12_err3_0212_d02'
     stop
  else
     call getarg(1,sign)
     read(sign(1:1),'(i1)')isign
     call getarg(2, outfile)
     length_outfile = len_trim(outfile)
  endif
      
!----------------------------------------------------------------
  write(wrf_file,'(a5,i5.5)')'fort.',e_unit+1
  call get_wrf_info(wrf_file, ix, jx, kx, times, proj) ! saved wrf info to proj 

!----------------------------------------------------------------
! Initilize and Read in namelist information
!----------------------------------------------------------------
  call read_namelist(ix, jx, kx)

!----------------------------------------------------------------
! Calculate the Radar position reference to WRF Model grid coordinate
!  if ( use_radar_rv ) then

!... get radar information from radar_data.info and saved to obs type data array
     call get_radar_info ( "radar_data.info" ) 

!... compute the position of every radar data points in wrf domain
     do nr = 1, raw%radar_stn_num
        call radar_position ( proj, nr ) 
     enddo

!  endif

! Get Radar data and process them( delete the data out of model range, data thinning, qc)
! (if used simulated data, first get the truth and then interp to radar data grid)
  
!   call real_obser ( ix, jx, kx, proj, times, isign )
   call get_wsr88d_radar ( ix, jx, kx, proj, times )

! calculate radar rv obser points
  num_rv = 0
  do nr  = 1, raw%radar_stn_num
     num_rv = num_rv + raw%radar(nr)%radar_stn%numObs
  enddo

! get obser point info
  allocate( ob        ( num_rv ) )
  allocate( obs_ii    ( num_rv ) )
  allocate( obs_jj    ( num_rv ) )
  allocate( obs_lat   ( num_rv ) )
  allocate( obs_lon   ( num_rv ) )
  allocate( obs_kk    ( num_rv ) )
  allocate( obs_hh    ( num_rv ) )
  allocate( obs_az    ( num_rv ) )
  allocate( obs_ra    ( num_rv ) )
  allocate( obs_el    ( num_rv ) )
  allocate( radar_ii  ( num_rv ) )
  allocate( radar_jj  ( num_rv ) )
  allocate( radar_elv ( num_rv ) )
  allocate( xlat      (ix, jx  ) )
  allocate( xlon      (ix, jx  ) )

  call get_variable2d( wrf_file, 'XLAT      ', ix, jx, 1, xlat  )
  call get_variable2d( wrf_file, 'XLONG     ', ix, jx, 1, xlon  )

  open(99, file=outfile(1:length_outfile)//'_so.3dvar',form='formatted')

  char_total_Radar( 1:40)='                                            '
  char_total_Radar(41:80)='                                            '
  char_total_Radar(1:14)=outfile(1:14)
  write(char_total_Radar(15:17),'(i3)')raw%radar_stn_num
  write(99,'(a80)')char_total_Radar

  write(99,'(a)')'  '

  iobs = 0
  do nr  = 1, raw%radar_stn_num

     write(99,'(a)')'  '

           char_ned( 1:40) = '                                        '
           char_ned(41:80) = '                                        '
           char_ned( 1: 4) = raw%radar(nr)%radar_stn%idc
           char_ned( 5: 7) = '   '
           char_ned( 8:19) = raw%radar(nr)%radar_stn%name(1:12)
     write(char_ned(20:27),'(f8.3)')raw%radar(nr)%radar_stn%lon
     write(char_ned(30:37),'(f8.3)')raw%radar(nr)%radar_stn%lat
     write(char_ned(40:47),'(f8.3)')raw%radar(nr)%radar_stn%elv
     write(char_ned(50:68),'( a19)')times(1:19)
     write(char_ned(69:74),'(  i6)')raw%radar(nr)%radar_stn%numObs
     write(char_ned(75:80),'(  i6)')1
     write(99,'(a80)')char_ned(1:80)

     write(99,'(a)')'  '
     write(99,'(a)')'  '

     reports: do lev = 1, raw%radar(nr)%radar_stn%numObs

              iobs = iobs + 1
              ob  ( iobs ) = raw%radar(nr)%rv(lev)
              obs_ii(iobs) = raw%radar(nr)%radar_point(lev)%ii
              obs_jj(iobs) = raw%radar(nr)%radar_point(lev)%jj
              obs_hh(iobs) = raw%radar(nr)%radar_point(lev)%hgt*1000.
              obs_az(iobs) = raw%radar(nr)%radar_point(lev)%azim
              obs_el(iobs) = raw%radar(nr)%radar_point(lev)%elev
              obs_ra(iobs) = raw%radar(nr)%radar_point(lev)%rdis

              radar_ii(iobs)  = raw%radar(nr)%radar_stn%i_radar
              radar_jj(iobs)  = raw%radar(nr)%radar_stn%j_radar
              radar_elv(iobs) = raw%radar(nr)%radar_stn%elv

              i = int(obs_ii(iobs))
              i1 = i+1
              j = int(obs_jj(iobs))
              j1 = j+1
              di = obs_ii(iobs) - real(i)
              di1 = real(i1) - obs_ii(iobs) 
              dj = obs_jj(iobs) - real(j)
              dj1 = real(j1) - obs_jj(iobs)

              obs_lat(iobs) = dj1*( di1*xlat(i,j )+di*xlat(i1,j) ) +  &
                               dj*( di1*xlat(i,j1)+di*xlat(i1,j1) )
              obs_lon(iobs) = di1*( dj1*xlon(i,j )+dj*xlon(i,j1) ) +  &
                               di*( dj1*xlon(i1,j)+dj*xlon(i1,j1) )

              write(99,'(A3,I3,I6,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)' ) &
                         raw%radar(nr)%radar_stn%idc(2:4), 128, lev,   &
                         times(1:19),                                  &
                         obs_lat(iobs), obs_lon(iobs), obs_hh(iobs), 1 

              write(99,'( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )' )  &
                         obs_hh(iobs),                          &
                         ob  ( iobs ), 0, raw%radar(nr)%radar_stn%err_rv, &
                         -99999.000, 9, -99999.000
                         
     end do reports
                         
  enddo
  close(99)

     

  open(99, file=outfile(1:length_outfile)//'.dat',form='formatted')
  write(99,*)'azimuth elevation radius   i       j       h      rv    lat      lon ' 
  do j=1,iobs
     write(99,'(9f8.3)')obs_az(j)*180./3.1415927, obs_el(j)*180./3.1415927, &
           obs_ra(j), obs_ii(j), obs_jj(j), obs_hh(j), &
           ob(j), obs_lat(j), obs_lon(j)
  enddo
  close(99)

    call parallel_finish()
end program so_to_3dvar
