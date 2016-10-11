subroutine get_all_obs( wrf_file, ix, jx, kx, times, proj )
!  PURPOSE: Get all Observations and stored in ob(obs%num), obs%type(obs%num), obs%err(obs%num),
!     obs%position(obs%num,4) and obs%sta(obs%num,4)
!  Origin:        02/19/2007 (Yonghui Weng)
use constants
use namelist_define
use mpi_module
use mapinfo_define
use obs_define
use map_utils
use netcdf
use wrf_tools
use radar
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
character (len=80), intent(in)       :: times

integer                              :: iobs, grid_id, fid,rcode

integer                              :: n, nr, lev, nvar, i, j, k, lvl, unit
real,dimension(3)                    :: center
real                                 :: center_io, center_jo

character (len=10)                   :: truthfile
character (len=80)                   :: filename
real, dimension(ix+1, jx, kx  )      :: u
real, dimension(ix, jx+1, kx  )      :: v
real, dimension(ix, jx, kx  )        :: t, qv
real, dimension(ix, jx, kx+1)        :: ph, phb
real, dimension(ix, jx, kx)          :: p, pb
character (len=6)                    :: pfile         !! out_##

!*******************************************************************************
! Get Model pressure and lat and long
!*******************************************************************************
  call get_variable3d( wrf_file, 'P         ', ix, jx, kx, 1, p )
  call get_variable3d( wrf_file, 'PB        ', ix, jx, kx, 1, pb)
  p = p + pb 
  call get_variable3d( wrf_file, 'PH        ', ix, jx, kx+1, 1, ph )
  call get_variable3d( wrf_file, 'PHB       ', ix, jx, kx+1, 1, phb)
  ph = (ph + phb)/9.8 
  call open_file(wrf_file, nf_nowrite, fid)
  rcode = nf_get_att_int(fid, nf_global, 'GRID_ID', grid_id)

!*******************************************************************************
! Get all raw observations into raw (type)
!*******************************************************************************
  if ( .not. use_ideal_obs ) then
!   call get_gtsobs_bufr ( times, ix, jx, kx, p, ph, proj )
   call get_gtsobs_3dvar ( times, ix, jx, kx, p, ph, proj )
!      filename = 'hurricane_best_track                                                                '
!      if ( expername .eq. 'hurricane ' ) &
!      call gtsobs_time_shift ( times, filename )
  endif
! Hurricane center position and minimum sea level pressure
if ( use_hurricane_PI ) then
   call rd_interp_besttrack ( times, 'hurricane_best_track', center )
   call latlon_to_ij( proj, center(1), center(2), center_io, center_jo )
   if(my_proc_id==0) write(*,'(a,5f8.2)')'Observation hurricane center =', center_io, center_jo, center(1:3)
endif
! Radar data 
if ( use_radar_rv ) then

!... get radar information from radar_data.info and saved to obs type data array
   call get_radar_info ( "radar_data.info" )

!... compute the position of every radar data points in wrf domain
   do nr = 1, raw%radar_stn_num      
      call radar_position ( proj, nr )
   enddo

!... get WSR-88D radar SO data
   if ( .not. use_ideal_obs  )call get_wsr88d_radar ( ix, jx, kx, proj, times )

!... ideal WSR-88D radar SO data 
   if ( use_ideal_obs ) call wsr88d_ideal_data( wrf_file, ix, jx, kx, proj )
  
endif 

! Airborne radar data
if ( use_airborne_rv ) then
  if ( .not. use_ideal_obs  ) call get_airborne ( ix, jx, kx, proj, times )
endif


!---edited by Minamide 2014.11.18
! Radiance data (satellite observation)
if ( use_radiance ) then
  if ( .not. use_ideal_obs  ) call get_radiance ( ix, jx, kx, proj, times )
endif


!  if ( use_radar_rv ) call output_simulated_rv ( "asimulated_rv" )
!   if ( use_simulated ) call simulated_obser ( wrf_file, ix, jx, kx, proj )
! Add all obs type observations into obs%dat(num_rv)


!*******************************************************************************
!. allocate variables
!*******************************************************************************
iobs = 0               ! iobs for all observation, not only Rv

!. calculate obs_gts records
if ( .not. use_ideal_obs ) then
   do n = 1, raw%gts%num
      if (raw%gts%slp(n,1).gt.80000. .and. raw%gts%slp(n,1).lt.105000.)iobs = iobs + 1
      if (raw%gts%pres(n,1,1).gt.40000. .and. raw%gts%pres(n,1,1).lt.105000.)iobs = iobs + 1
      do j = 1, raw%gts%levels(n)
         if (raw%gts%spd(n,j,1).ge.0. .and. raw%gts%spd(n,j,1).lt.200. .and.  &
             raw%gts%wd(n,j,1).ge.0. .and. raw%gts%wd(n,j,1).le.360. ) iobs = iobs + 2
!......u,v
!            if ( abs(raw%gts%spd(n,j,1)) .lt. 200. ) iobs = iobs + 1
!            if ( abs(raw%gts%wd (n,j,1)) .le. 360. ) iobs = iobs + 1
         if (raw%gts%t(n,j,1).gt.0. .and. raw%gts%t(n,j,1).lt.350. ) iobs = iobs + 1
         if (raw%gts%td(n,j,1).gt.0. .and. raw%gts%td(n,j,1).lt.350. .or.    &
             raw%gts%rh(n,j,1).gt.0. .and. raw%gts%rh(n,j,1).le.100. ) iobs = iobs + 1
      enddo
   enddo
endif

!. calculate WSR-88D radar RV records
if ( use_radar_rv ) then
   do nr  = 1, raw%radar_stn_num
      iobs = iobs + raw%radar(nr)%radar_stn%numObs
   enddo
endif

!. calculate airborne radar records
if ( use_airborne_rv) iobs = iobs + raw%airborne%num

!---edited by Minamide 2014.11.18
!. calculate satellite radiance records
if ( use_radiance) iobs = iobs + raw%radiance%num

!. calculate tracker records
if ( use_hurricane_PI ) iobs = iobs + 3            ! for hurricane center lat and lon and slp

if(my_proc_id==0) write(*,*)iobs,' observation data will be loaded'

allocate( obs%dat      ( iobs    ) )
allocate( obs%type     ( iobs    ) )
allocate( obs%err      ( iobs    ) )
allocate( obs%position ( iobs, 4 ) )
allocate( obs%sta      ( iobs, 4 ) )
allocate( obs%roi      ( iobs, 3 ) )
allocate( obs%sat      ( iobs    ) )
allocate( obs%ch       ( iobs    ) )

obs%dat        = -888888.
obs%type       = '          '
obs%err        = -888888.
obs%position   = -888888.
obs%sta        = -888888.
obs%roi        = -888888
obs%sat        = '            '
obs%ch         = -888888

!----------------------------------------------------------------------------
obs%num = 0               ! obs%num for all observation, not only Rv
!*******************************************************************************
!. save obs data to obs varry
!....... Surface U, V, T, RH, SLP
!*******************************************************************************
   if ( use_surface ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_surface_data( wrf_file, ix, jx, kx, proj, 'surface ', &
              datathin_surface, hroi_surface, vroi_surface, grid_id )
      endif
   endif


!....... Metar U, V, T, RH, SLP
   if ( use_metar   ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_surface_data( wrf_file, ix, jx, kx, proj, 'metar   ', &
              datathin_metar  , hroi_metar  , vroi_metar  , grid_id )
      endif
   endif

!....... Buoy and Ship U, V, T, RH, SLP
   if ( use_sfcshp  ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_surface_data( wrf_file, ix, jx, kx, proj, 'sfcshp  ', &
              datathin_sfcshp , hroi_sfcshp , vroi_sfcshp , grid_id )
      endif
   endif

!....... SSMI U, V, T, RH, SLP
   if ( use_spssmi  ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_surface_data( wrf_file, ix, jx, kx, proj, 'spssmi  ', &
              datathin_spssmi , hroi_spssmi , vroi_spssmi , grid_id )
      endif
   endif

!....... Sounding U, V, T, RH, SLP
   if ( use_sounding ) then
      if ( use_ideal_obs ) then 
         call ideal_sounding_data( wrf_file, ix, jx, kx )
      else
         call sort_sounding_data( wrf_file, ix, jx, kx, proj, 'sounding', &
              datathin_sounding, hroi_sounding, vroi_sounding, grid_id)
      endif
   endif

!....... Get Profiler U, V
   if ( use_profiler ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_upperair_data( wrf_file, ix, jx, kx, proj, 'profiler', &
              datathin_profiler, hroi_profiler, vroi_profiler, grid_id)
      endif
   endif

!....... Get Aircft U, V, T
   if ( use_aircft ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_upperair_data( wrf_file, ix, jx, kx, proj, 'aircft  ', &
              datathin_aircft  , hroi_aircft  , vroi_aircft  , grid_id)
      endif
   endif

!....... Get Atovs T and Td
   if ( use_atovs ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_upperair_data( wrf_file, ix, jx, kx, proj, 'atovs   ', &
              datathin_atovs   , hroi_atovs   , vroi_atovs   , grid_id)
      endif
   endif

!....... Get SatwndU, V, T
   if ( use_satwnd ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_upperair_data( wrf_file, ix, jx, kx, proj, 'satwnd  ', &
              datathin_satwnd  , hroi_satwnd  , vroi_satwnd  , grid_id)
      endif
   endif

!.......Get GPSPW
   if ( use_gpspw ) then
      if ( use_ideal_obs ) then
         write(*,*)' Please choose use_sounding and use_ideal '
         stop
      else
         call sort_surface_data( wrf_file, ix, jx, kx, proj, 'gpspw   ', &
              datathin_gpspw, hroi_gpspw, vroi_gpspw, grid_id )
      endif
   endif

!....... Get Rv
   if ( use_radar_rv ) then
      call sort_radarRV_data( wrf_file, ix, jx, kx, proj, 'RadarRV   ', &
              datathin_radar  , hroi_radar  , vroi_radar  , grid_id )
   endif

!....... Get Airbone Rv
   if ( use_airborne_rv ) then
      call sort_radarRV_data( wrf_file, ix, jx, kx, proj, 'AirborneRV', &
              datathin_airborne  , hroi_airborne  , vroi_airborne  , grid_id )
   endif

!---edited by Minamide 2014.11.18
!....... Get Satellite radiance
   if ( use_radiance ) then
      call sort_radiance_data( wrf_file, ix, jx, kx, proj, 'radiance', &
              datathin_radiance, hroi_radiance, vroi_radiance, grid_id )
   endif

!... hurricane center lat, lon and minimum sea-level-pressure
  if ( use_hurricane_PI ) then 
     if( center(2) .ge. -180. .and. center(2) .le. 360. .and. &
         center(1) .ge.  -90. .and. center(1) .le.  90. .and. &
         center(3) .ge.  800. .and. center(3) .le. 1200. )then

         obs%num                 = obs%num + 1
         obs%dat(obs%num)      = center(3) - pr/100.
         obs%type(obs%num) = 'min_slp   '
         obs%err(obs%num) = 5.
         obs%position(obs%num,1)  = center_io
         obs%position(obs%num,2)  = center_jo
         obs%position(obs%num,3)  = 1
         obs%roi     (obs%num,1)  = hroi_hurricane_PI 
         obs%roi     (obs%num,2)  = vroi_hurricane_PI

      else
         if(my_proc_id==0) then
           write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*)'No hurricane center data will be assimilateed. '
           write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*)
         endif
      endif
   endif

   if(my_proc_id==0) write(*,*)obs%num,' observation data will be assimilated'

end subroutine get_all_obs
!=======================================================================================

subroutine rd_hurricane_track ( times, filename, center )
implicit none
character (len=80), intent(in)        :: times
character (len=*), intent(in)        :: filename
real, dimension(3), intent(out)       :: center
character          :: year*4, month*2, day*2, hour*2, minute*2
real, dimension(3) :: dat
integer            :: length, iost

length = len_trim(filename)
open(39, file=filename(1:length), status='old', form = 'formatted', iostat = iost )
if( iost .ne. 0 ) stop 'Cannot find the hurricane track file '

do_get_track_loop : do
   read(39, '(a2,1x,a2,3f7.1)', iostat = iost)day, hour, dat(1:3)
   if( iost .ne. 0 ) exit
   if( day == times(9:10) .and. hour== times(12:13) ) then
       center = dat
       exit
   endif
end do do_get_track_loop
close(39)
write(*,*)'hurricane center from obs is:',center

end subroutine rd_hurricane_track
!=======================================================================================

subroutine rd_interp_besttrack ( times, filename, center )

   use mpi_module
   implicit none
   character (len=80), intent(in)        :: times
   character (len=*), intent(in)         :: filename
   real, dimension(3), intent(out)       :: center

   integer            :: length, iost, hours_from_date
   integer, dimension(3) :: year, month, day, hour, dhour     !1=current, 2=pre, 3=post
   real, dimension(3) :: lat, lon, wsp, slp
   integer            :: iyear, imonth, iday, ihour, ilat, ilon, iwsp, islp, diff_hour
   real               :: awsp
   character (len=1)  :: clat, clon
   character (len=3)  :: nhc

   dhour(1:3)=-999999
   read(times, '(i4, 1x, i2, 1x, i2, 1x, i2)')year(1), month(1), day(1), hour(1)
   dhour(1) = hours_from_date ( year(1), month(1), day(1), hour(1))

   length = len_trim(filename)
   open(39, file=filename(1:length), status='old', form = 'formatted', iostat = iost )
   if( iost .ne. 0 ) then
     center = -9999999.0
     if ( my_proc_id == 0 ) write(*,*)'CANNOT find best-track file!!!!!!'
     return
   end if

   read(39, '(a3)', iostat = iost)nhc
   backspace(39)
   if (nhc .eq. "NHC" ) then   !tcvital
      read(39, '(32x,i4, a1, i5, a1, 4x, i4, i5)', iostat = iost)  &
          ilat, clat, ilon, clon, iwsp, islp
      if ( iost .ne. 0 ) return
      if (clat == 'S') ilat = -ilat
      if (clon == 'W') ilon = -ilon
      center(1) = ilat/10.
      center(2) = ilon/10.
      center(3) = islp*1.0

   else if (nhc .eq. "AL," ) then   ! best track
        do_get_track_loop : do
           !-----Best Track --b-deck
           read(39, '(8x, i4,3i2, 16x, i4, a1, 1x, i5, a1, 1x, i4, 1x, i5)', iostat = iost)  &
                   iyear, imonth, iday, ihour, ilat, clat, ilon, clon, iwsp, islp
           awsp = iwsp*0.514444
           !-----tcvitals
           !read(39, '(19x, i4,2i2,1x,i2,2x,i4, a1, i5, a1,9x,i4, 10x, i3)', iostat = iost)  &
           !        iyear, imonth, iday, ihour, ilat, clat, ilon, clon, islp, iwsp 
           !awsp = iwsp*1.0
           if( iost .ne. 0 ) exit
           if (clat == 'S') ilat = -ilat
           if (clon == 'W') ilon = -ilon
           diff_hour = hours_from_date ( iyear, imonth, iday, ihour)
           if (diff_hour == dhour(1) ) then
              lat(2) = ilat/10.
              lon(2) = ilon/10.
              wsp(2) = awsp
              slp(2) = islp*1.0
              dhour(2) = 99999
              lat(3) = ilat/10.
              lon(3) = ilon/10.
              wsp(3) = awsp
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
              wsp(2) = awsp
              slp(2) = islp*1.0
              dhour(2) = dhour(1) - diff_hour
           else
              year(3) = iyear
              month(3) = imonth
              day(3) = iday
              hour(3) = ihour
              lat(3) = ilat/10.
              lon(3) = ilon/10.
              wsp(3) = awsp
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
        end if

   end if
   if ( my_proc_id == 0 )write(*,*)'hurricane center from obs is:',center
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
!TODO
!subroutine get_gtsobs_bufr (times,ix,jx,kx,p,ph,proj)
!!get obs in bufr format
!use constants
!use namelist_define
!use obs_define
!use mapinfo_define
!use map_utils
!use mpi_module
!implicit none
!
!end subroutine get_gtsobs_bufr
!=======================================================================================

subroutine get_gtsobs_3dvar ( times, ix, jx, kx, p, ph, proj )
!get obs in little_r format
use constants
use namelist_define
use obs_define
use mapinfo_define
use map_utils
use mpi_module
implicit none
type(proj_info)                          :: proj
character (len=80), intent(in)           :: times
character (len=80)                       :: obs_3dvar_file
integer, intent(in)                      :: ix, jx, kx
real, dimension(ix, jx, kx), intent(in)  :: p
real, dimension(ix, jx, kx+1), intent(in):: ph
character (len=180)                      :: fmt_fmt, info_fmt, srfc_fmt, each_fmt
integer                                  :: synop,metar,ship,buoy,bogus,temp,amdar,airep,tamdar,pilot, &
                                            satem, satob, gpspw, gpszd, gpsrf,gpsep,ssmt1,ssmt2,     &
                                            tovs, qscat, profl, airsr, other, total 
character (len=40)                       :: name
character (len= 5)                       :: id
real, dimension(3)                       :: pw
real, dimension(7, 3, 3000)              :: obs_data     !7=pres, spd, wd, height, t, td, rh
                                                         !3=data, qc, error
                                                         !3000=levels
integer                                  :: obs_level
integer, dimension(7)                    :: qcint
integer                                  :: i, j, k, n, io, jo, kk, ksund, iost
real                                     :: aio,ajo,p_diff,z_diff,po_diff,zo_diff
integer, parameter                       :: numlevs=300

if ( my_proc_id == 0 ) then
write(6, *)'   '
write(6, *)'---------------------------------------------------'
write(6, *)'.... Getting GTS (3dvar format) Data ....'
endif
!  get data
obs_3dvar_file = "                                                                                "
obs_3dvar_file = 'obs_3dvar_'// times(1:4)//times(6:7)//times(9:10)//times(12:13)//           &
                 times(15:16)//times(18:19)
i = len_trim( obs_3dvar_file )
if ( my_proc_id == 0 ) write(*, *)'.... ', obs_3dvar_file(1:i),'  ....'
open (10, file = obs_3dvar_file(1:i), status = 'old', form = 'formatted', iostat = iost )
if( iost .ne. 0 ) then
   if ( my_proc_id == 0 ) then
     write(*,*)obs_3dvar_file(1:i),' does not exist, please check it.'
     write(*,*)'OBS(in 3dvar format) data will be not assimilated !'
   endif
      use_surface  = .false.
      use_metar    = .false.
      use_profiler = .false.
      use_aircft   = .false.
      use_spssmi   = .false.
      use_satwnd   = .false.
      use_sfcshp   = .false.
      use_sounding = .false.
      return
endif

!... get the data number
read(10, fmt='(7x,i7)')total                                                 
read(10, fmt='(6(7x,i7,2x))') synop, metar, ship, buoy, bogus, temp
read(10, fmt='(6(7x,i7,2x))') amdar, airep, tamdar,pilot,satem, satob
read(10, fmt='(6(7x,i7,2x))') gpspw, gpszd, gpsrf, gpsep, ssmt1, ssmt2
read(10, fmt='(5(7x,i7,2x))') tovs,qscat, profl, airsr, other
raw%gts%num        = total
!....  platforms
!
!   Name    WMO Codes     WMO Code names
!   synop    12,14       'SYNOP','SYNOP MOBIL'                      ! Land synoptic reports
!   ship     13          'SHIP'                                     ! Ship and moored buoy synoptic reports
!   metar    15,16       'METAR','SPECI'
!   buoy     18          'BUOY'                                     ! Drifting buoy reports
!   pilot    32,33,34    'PILOT','PILOT SHIP','PILOT MOBIL'         ! Upper air wind profiles: UG UH UI UP UQ  UY
!   sound    35,36,37,38 'TEMP','TEMP SHIP, 'TEMP DROP','TEMP MOBIL'! Upper air temperature, humidity and wind profiles: UE UF UK UM US
!   amdar    42          'AMDAR'                                    ! Aircraft Meteorological Data And Reporting Relay: UD
!   satem    86          'SATEM'
!   satob    88          'SATOB'                                    ! satwind: JMA satob-coded winds
!   airep    96,97       'AIREP'                                    ! Aircraft reports: UA UB
!   gpspw    111         'GPSPW'
!   gpsztd   114         'GPSZD'
!   gpsref   116         'GPSRF'
!   gpseph   118         'GPSEP'
!   ssmt1    121         'SSMT1'
!   ssmt2    122         'SSMT2'
!   ssmi     125,126     'SSMI'                                     ! Micro-wave radiances from Special Sensor Microwave Imager on DMSP satellites
!   tovs     131         'TOVS'                                     ! TIROS Operational Vertical Sounder
!   qscat    281         'Quikscat'                                 ! Seawinds data obtained from QuikSCAT
!   profl    132         'Profilers'                                ! 
!   bogus    135         'BOGUS'
!   airs     133         'AIRSRET'
!   other Any other code 'UNKNOWN'

!.... reduced to BUFRLIB classes:
!   sounding  32,33,34,35,36,37,38,                  UPPER-AIR (RAOB, PIBAL, RECCO, DROPS) REPORTS 
!   surface   12,14,                                 SURFACE LAND (SYNOPTIC) REPORTS
!   metar     15,16,                                 SURFACE LAND (METAR) REPORTS
!   profiler  132                                    WIND PROFILER REPORTS
!   aircft    42,96,97                               AIREP/PIREP, AMDAR (ASDAR/ACARS), E-ADAS (AMDAR BUFR) AIRCRAFT REPORTS
!   sfcshp    13,18,                                 SURFACE MARINE (SHIP, BUOY, C-MAN PLATFORM) REPORTS
!   satwnd    88,                                    SATELLITE-DERIVED WIND REPORTS
!   spssmi    125,126                                DMSP SSM/I RETRIEVAL PRODUCTS (REPROCESSED WIND SPEED, TPW)

!... skip blank line
do i = 1, 12
   read(10,*)
enddo

!... get format information
read(10,'(a180)')fmt_fmt 
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
info_fmt(1:180)=' '
info_fmt(1:j-i)=fmt_fmt(i+1:j)

read(10,'(a180)')fmt_fmt 
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
srfc_fmt(1:180)=' '
srfc_fmt(1:j-i)=fmt_fmt(i+1:j)

read(10,'(a180)')fmt_fmt 
i = index(fmt_fmt,'=')
j = len_trim(fmt_fmt)
each_fmt(1:180)=' '
each_fmt(1:j-i)=fmt_fmt(i+1:j)

read(10,*)

!... allocate
! JPOTERJOY: use numlevels to allocate arrays
allocate( raw%gts%platform ( total ) ) 
allocate( raw%gts%date     ( total ) ) 
allocate( raw%gts%latitude ( total ) ) 
allocate( raw%gts%longitude ( total ) ) 
allocate( raw%gts%elevation ( total ) ) 
allocate( raw%gts%slp      ( total, 3 ) ) 
allocate( raw%gts%pw      ( total, 3 ) ) 
allocate( raw%gts%levels   ( total ) ) 
allocate( raw%gts%pres     ( total, numlevs, 3 ) ) 
allocate( raw%gts%spd      ( total, numlevs, 3 ) ) 
allocate( raw%gts%wd       ( total, numlevs, 3 ) ) 
allocate( raw%gts%height   ( total, numlevs, 3 ) ) 
allocate( raw%gts%t        ( total, numlevs, 3 ) ) 
allocate( raw%gts%td       ( total, numlevs, 3 ) ) 
allocate( raw%gts%rh       ( total, numlevs, 3 ) ) 

!... initialize
raw%gts%platform = '            '
raw%gts%date     = '                   '
raw%gts%latitude = -888888.
raw%gts%longitude= -888888.
raw%gts%elevation= -888888.
raw%gts%slp      = -888888.
raw%gts%pw      = -888888.
raw%gts%levels   = 0
raw%gts%pres     = -888888.
raw%gts%spd      = -888888.
raw%gts%wd       = -888888.
raw%gts%height   = -888888.
raw%gts%t        = -888888.
raw%gts%td       = -888888.
raw%gts%rh       = -888888.

!... Get data and sort to raw%gts
do n = 1, total
   read(10, fmt=info_fmt)raw%gts%platform(n), raw%gts%date(n), name, &
                         obs_level, raw%gts%latitude(n),     &
                         raw%gts%longitude(n), raw%gts%elevation(n), id
   read(10, fmt=srfc_fmt)raw%gts%slp(n,1),qcint(1),raw%gts%slp(n,3), &
                         raw%gts%pw(n,1), qcint(2),raw%gts%pw(n,3)
   raw%gts%slp(n,2) = real(qcint(1))
   raw%gts%pw(n,2) = real(qcint(2))
   
   do k = 1, obs_level
      read(10, fmt=each_fmt)((obs_data(i,1,k),qcint(i),obs_data(i,3,k)),i=1,7) 
      obs_data(1:7,2,k) = real(qcint(1:7))
      if (raw%gts%platform(n) .eq. "FM-12 MTCS s") then
         !obs_data(2,3,k)=obs_data(2,3,k)*10.0
         !obs_data(3,3,k)=obs_data(3,3,k)*10.0
         obs_data(2,3,k)=7.0
         obs_data(3,3,k)=7.0
      endif
   enddo

!!.Selecting the data closed to model level.
   call latlon_to_ij( proj, raw%gts%latitude(n), raw%gts%longitude(n), aio, ajo )
   io = int(aio)
   jo = int(ajo)
   if( io.gt.0 .and. io.lt.ix .and. jo.gt.0 .and. jo.lt.jx) then
     ksund = 0
     do k = 1, kx-1
        p_diff = p(io,jo,k) - p(io,jo,k+1)
        z_diff = ph(io,jo,k+1) - ph(io,jo,k)
        kk = -99
        do i = 1, obs_level
           if ( obs_data(1,1,i) .gt. 1000. .and. obs_data(1,1,i) .lt. 105000. ) then   ! decided by pressure
              po_diff=obs_data(1,1,i)-p(io,jo,k)
              if (po_diff.gt.0 .and. po_diff.lt.p_diff ) then
                 p_diff = po_diff
                 kk = i
              endif
           end if
           if ( obs_data(4,1,i) .ge. 0. .and. obs_data(4,1,i) .lt. 90000. )then   ! decided by height
              zo_diff=obs_data(4,1,i)-ph(io,jo,k)
              if (zo_diff.gt.0 .and. zo_diff.lt.z_diff ) then
                 z_diff = zo_diff
                 kk = i
              endif
           endif
        end do
        if ( kk > 0 .and. kk < obs_level+1 ) then
           ksund = ksund + 1
           do j = 1, 3
              raw%gts%pres(n, ksund, j) = obs_data(1, j, kk)
              raw%gts%spd(n, ksund, j) = obs_data(2, j, kk)
              raw%gts%wd(n, ksund, j) = obs_data(3, j, kk)
              raw%gts%height(n, ksund, j) = obs_data(4, j, kk)
              raw%gts%t(n, ksund, j) = obs_data(5, j, kk)
              raw%gts%td(n, ksund, j) = obs_data(6, j, kk)
              raw%gts%rh(n, ksund, j) = obs_data(7, j, kk)
           enddo
        endif
     enddo    !do k = 1, kx-1 
     raw%gts%levels(n) = ksund
   end if
!!.Without selecting obs_level close to model_level.
!  raw%gts%levels(n) = obs_level
!  do i = 1, 3
!     raw%gts%pres(n, 1:obs_level, i) = obs_data(1, i, 1:obs_level)
!     raw%gts%spd(n, 1:obs_level, i) = obs_data(2, i, 1:obs_level)
!     raw%gts%wd(n, 1:obs_level, i) = obs_data(3, i, 1:obs_level)
!     raw%gts%height(n, 1:obs_level, i) = obs_data(4, i, 1:obs_level)
!     raw%gts%t(n, 1:obs_level, i) = obs_data(5, i, 1:obs_level)
!     raw%gts%td(n, 1:obs_level, i) = obs_data(6, i, 1:obs_level)
!     raw%gts%rh(n, 1:obs_level, i) = obs_data(7, i, 1:obs_level)
!  enddo
enddo   !do n = 1, total
close(10)

end subroutine get_gtsobs_3dvar
!=======================================================================================

  subroutine get_wsr88d_radar ( ix, jx, kx, proj, times )
  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use mpi_module
  use map_utils
  use netcdf
  use wrf_tools
  use radar
  implicit none
  type(proj_info), intent(in)         :: proj
  integer, intent(in)                 :: ix, jx, kx
  real, dimension(ix, jx      )       :: xlong
  real, dimension(ix, jx, kx+1)       :: ph
  real, dimension(ix, jx, kx+1)       :: phb
  character (len=80)                  :: radar_file
  character (len=80)                  :: times
  integer                             :: i, n, iost, num
  real                                :: az, el, ra, rw
  real                                :: s, h, rx, ry, ir1, jr1

 if ( my_proc_id == 0 ) then
  write(6, *)'   '
  write(6, *)'---------------------------------------------------'
  write(6, *)'.... Getting Real Observation Data  ....'
 endif
! get data
  if ( use_radar_rv ) then

  radar_file = "                                                                            "
  do_radar_stn        : do n = 1, raw%radar_stn_num

     i = len_trim(raw%radar(n)%radar_stn%idc)
     radar_file(2:4) = raw%radar(n)%radar_stn%idc(1:i)
     radar_file = raw%radar(n)%radar_stn%idc(1:i)//'_'//times(1:4)//times(6:7)//         &
                  times(9:10)//times(12:13)//times(15:16)//times(18:19)//'_so'
     i = len_trim( radar_file )
     if(my_proc_id==0) write(6, *)'.... ', radar_file(1:i),'  ....'
     open (10, file = radar_file(1:i), status = 'old', form = 'formatted', iostat = iost )
     if( iost .ne. 0 ) then
        if(my_proc_id==0) write(*,*)radar_file(1:i),' does not exsit, please check it.'
!           stop 'get_wsr88d_radar'
        cycle do_radar_stn
     endif

!...... get the data number of every radar
     num = 0
     do_get_raw_data_loop_read : do
        read(10, '(4f12.3)', iostat = iost ) az, el, ra, rw
        if( iost .ne. 0 ) exit
        num = num + 1
     end do do_get_raw_data_loop_read

!...... allocate
     allocate( raw%radar(n)%radar_point( num ) )
     allocate( raw%radar(n)%rf( num ) )
     allocate( raw%radar(n)%rv( num ) )

!...... get data
     rewind(10)
     num = 0
     do_get_raw_data_loop : do

!......... Enkf with odd data, and verify with even data
        read(10, '(4f12.3)', iostat = iost ) az, el, ra, rw
        if( iost .ne. 0 ) exit
        az = az * rad_per_deg
        el = el * rad_per_deg
  
!......... compute the surface distance and ASL height
        call calc_radar_data_point_s_h (s, h, ra, el, raw%radar(n)%radar_stn%elv)
        if ( h >= 0.05 .and. h <= 15. ) then

!............ compute data pint's position in model domain -- i, j
           rx = s * sin( az )
           ry = s * cos( az )
           ir1 = raw%radar(n)%radar_stn%i_radar + ( rx / proj%dx *1000. )
           jr1 = raw%radar(n)%radar_stn%j_radar + ( ry / proj%dx *1000. )

!............ evaluate
           if ( ir1 > 2 .and. ir1 < proj%nx-1 .and. jr1 > 2 .and. jr1 < proj%ny-1 ) then
              num = num + 1
              raw%radar(n)%radar_point(num)%ii   = ir1
              raw%radar(n)%radar_point(num)%jj   = jr1
              raw%radar(n)%radar_point(num)%hgt  = h
              raw%radar(n)%radar_point(num)%rdis = ra
              raw%radar(n)%radar_point(num)%azim = az
              raw%radar(n)%radar_point(num)%elev = el
              raw%radar(n)%rv( num )             = rw
           endif
        endif

     end do do_get_raw_data_loop 
     raw%radar(n)%radar_stn%numObs = num
     close (10)

  end do do_radar_stn

  endif   !if ( use_radar_rv ) then
end subroutine get_wsr88d_radar

!=======================================================================================
  subroutine get_airborne ( ix, jx, kx, proj, times )
  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils
  use netcdf
  use wrf_tools
  use radar
  implicit none
  type(proj_info), intent(in)         :: proj
  integer, intent(in)                 :: ix, jx, kx
  real, dimension(ix, jx      )       :: xlong
  real, dimension(ix, jx, kx+1)       :: ph
  real, dimension(ix, jx, kx+1)       :: phb
  character (len=80)                  :: radar_file
  character (len=80)                  :: times
  character (len=12)                  :: so_time
  integer                             :: i, n, iost, num
  real                                :: lat, lon, height, az, el, ra, rw
  real                                :: s, h, rx, ry, ir1, jr1, is, js

  if ( my_proc_id == 0 ) then
    write(6, *)'   '
    write(6, *)'---------------------------------------------------'
    write(6, *)'.... Getting Airbone Data  ....'
  endif
! get data
  if ( use_airborne_rv ) then

  radar_file = "                                                                            "
  radar_file = 'airborne_'//times(1:4)//times(6:7)//         &
               times(9:10)//times(12:13)//times(15:16)//'_so'
  i = len_trim( radar_file )
  if ( my_proc_id == 0 ) write(6, *)'.... ', radar_file(1:i),'  ....'
  open (10, file = radar_file(1:i), status = 'old', form = 'formatted', iostat = iost )
  if( iost .ne. 0 ) then
      if ( my_proc_id == 0 ) write(*,*)'============================================'
      if ( my_proc_id == 0 ) write(*,*)radar_file(1:i),' does not exist, please check it.'
      if ( my_proc_id == 0 ) write(*,*)'============================================'
!         stop 'get_airborne'
  endif

!...... get the data number
     num = 0
     do_get_raw_data_loop_read : do
        read(10, '(a12,7f10.3)', iostat = iost ) so_time, lat, lon, height, az, el, ra, rw
        if( iost .ne. 0 ) exit
        num = num + 1
     end do do_get_raw_data_loop_read

!...... allocate
     allocate( raw%airborne%lat( num ) )
     allocate( raw%airborne%lon( num ) )
     allocate( raw%airborne%height( num ) )
     allocate( raw%airborne%radar_ii( num ) )
     allocate( raw%airborne%radar_jj( num ) )

     allocate( raw%airborne%azim( num ) )
     allocate( raw%airborne%elev( num ) )
     allocate( raw%airborne%range( num ) )
     allocate( raw%airborne%rv( num ) )
     allocate( raw%airborne%ii( num ) )
     allocate( raw%airborne%jj( num ) )
     allocate( raw%airborne%hh( num ) )

!...... get data
     rewind(10)
     num = 0
     do_get_raw_data_loop : do

!......... Enkf with odd data, and verify with even data
        read(10, '(a12,7f10.3)', iostat = iost ) so_time, lat, lon, height, az, el, ra, rw
        if( iost .ne. 0 ) exit
        az = az * rad_per_deg
        el = el * rad_per_deg

!......... compute the surface distance and ASL height
        call calc_radar_data_point_s_h (s, h, ra, el, height)
        if ( h >= 0.05 .and. h <= 15. ) then

!............ calculate radar center's position according to wrf domain grid
           call latlon_to_ij( proj, lat, lon, is, js )

!............ compute data pint's position in model domain -- i, j
           rx = s * cos( az )
           ry = s * sin( az )
           ir1 = is + ( rx / proj%dx *1000. )
           jr1 = js + ( ry / proj%dx *1000. )

!............ evaluate
           if ( ir1 > 2 .and. ir1 < proj%nx-1 .and. jr1 > 2 .and. jr1 < proj%ny-1 ) then
              num = num + 1
              raw%airborne%lat(num) = lat
              raw%airborne%lon(num) = lon
              raw%airborne%height(num) = height
              raw%airborne%radar_ii(num) = is
              raw%airborne%radar_jj(num) = js

              raw%airborne%azim(num) = az
              raw%airborne%elev(num) = el
              raw%airborne%range(num) = ra
              raw%airborne%rv(num) = rw
              raw%airborne%ii(num) = ir1
              raw%airborne%jj(num) = jr1
              raw%airborne%hh(num) = h
           endif
        endif

     end do do_get_raw_data_loop
     raw%airborne%num = num
     close (10)

  endif   !if ( use_airborne_rv ) then
end subroutine get_airborne
!=======================================================================================


!---edited by Minamide 2014.11.18 (for satellite data)
!=======================================================================================
  subroutine get_radiance ( ix, jx, kx, proj, times )
  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils
  use netcdf
  use wrf_tools
  use radar
  implicit none
  type(proj_info), intent(in)         :: proj
  integer, intent(in)                 :: ix, jx, kx
  real, dimension(ix, jx      )       :: xlong
  real, dimension(ix, jx, kx+1)       :: ph
  real, dimension(ix, jx, kx+1)       :: phb
  character (len=80)                  :: radiance_file
  character (len=80)                  :: times
  character (len=12)                  :: so_time
  character (len=12)                  :: sat_id
  integer                             :: i, n, iost, num
  integer                             :: ch_info,hroi_rad,hroi_drad
  real                                :: lat, lon, tb, err
  real                                :: s, h, rx, ry, ir1, jr1, is, js

  if ( my_proc_id == 0 ) then
    write(6, *)'   '
    write(6, *)'---------------------------------------------------'
    write(6, *)'.... Getting Satellite Radiance Data  ....'
  endif
! get data
  if ( use_radiance ) then

  radiance_file = ""
  radiance_file = 'radiance_'//times(1:4)//times(6:7)//         &
               times(9:10)//times(12:13)//times(15:16)//'_so'
  i = len_trim( radiance_file )
  if ( my_proc_id == 0 ) write(6, *)'.... ', radiance_file(1:i),'  ....'
  open (10, file = radiance_file(1:i), status = 'old', form = 'formatted', iostat =iost )
  if( iost .ne. 0 ) then
      if ( my_proc_id == 0 )write(*,*)'============================================'
      if ( my_proc_id == 0 )write(*,*)radiance_file(1:i),' does not exist, please check it.'
      if ( my_proc_id == 0 )write(*,*)'============================================'
!         stop 'get_radiance'
  endif

!...... get the data number
     num = 0
     do_get_raw_data_loop_read : do
        read(10, '(2a12,i12,3f12.3)', iostat = iost ) so_time, sat_id, ch_info, lat, lon, tb
        if( iost .ne. 0 ) exit
        num = num + 1
     end do do_get_raw_data_loop_read

!...... allocate
     allocate( raw%radiance%lat( num ) )
     allocate( raw%radiance%lon( num ) )
     allocate( raw%radiance%platform( num ) )
     allocate( raw%radiance%ch( num ) )
     allocate( raw%radiance%ii( num ) )
     allocate( raw%radiance%jj( num ) )
     allocate( raw%radiance%tb( num ) )
     allocate( raw%radiance%hroi( num ) )
     allocate( raw%radiance%hroi_d( num ) )
     allocate( raw%radiance%err( num ) )

!...... get data
     rewind(10)
     num = 0
     do_get_raw_data_loop : do

!......... Enkf with odd data, and verify with even data
        read(10, '(2a12,i12,3f12.3,2i12,f12.3)', iostat = iost ) so_time, sat_id, ch_info, lat, lon, tb, hroi_rad,hroi_drad,err
        if( iost .ne. 0 ) exit
!......... calculate radar center's position according to wrf domain grid
        call latlon_to_ij( proj, lat, lon, is, js )
!......... evaluate
           if ( is > 2 .and. is < proj%nx-1 .and. js > 2 .and. js <proj%ny-1) then
              num = num + 1
              raw%radiance%lat(num) = lat
              raw%radiance%lon(num) = lon
              raw%radiance%platform(num) = sat_id
              raw%radiance%ch(num) = ch_info
              raw%radiance%ii(num) = is
              raw%radiance%jj(num) = js
              raw%radiance%tb(num) = tb
              raw%radiance%hroi(num) = (hroi_rad*1000)/proj%dx
              raw%radiance%hroi_d(num) = (hroi_drad*1000)/proj%dx
              raw%radiance%err(num) = err
           else
           endif

     end do do_get_raw_data_loop
     raw%radiance%num = num
     close (10)

  endif   !if ( use_radiance ) then
end subroutine get_radiance


!=======================================================================================
  subroutine output_simulated_rv ( filename )
  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils
  use netcdf
  use wrf_tools
  implicit none
  character (len=*), intent(in)      :: filename
  character (len=8) :: stid
  integer       :: i, nr, num, k, length, nlev, nflag
  real          :: tim

! output in GrADS format, x: distance, y: azimuth, z: elevation angle levels,
!                         t: radar number;    rv: rv(x,y,z)
  length = len_trim( filename )

! each radar a file
  do nr = 1, raw%radar_stn_num

!... output grd file
  i = len_trim(raw%radar(nr)%radar_stn%idc)
  open(10, file = filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.grd',  &
           form='unformatted', access='direct', recl=9*4)
  do num = 1, raw%radar(nr)%radar_stn%numObs
     tim = 0.0
     nlev = 1
     nflag = 1
     write(stid,'(i8.8)')num
     write(10)stid,raw%radar(nr)%radar_point(num)%jj,raw%radar(nr)%radar_point(num)%ii,tim,nlev,nflag, &
              raw%radar(nr)%rv(num),raw%radar(nr)%radar_point(num)%hgt
  enddo
  nlev = 0
  write(10)stid,raw%radar(nr)%radar_point(num-1)%jj,raw%radar(nr)%radar_point(num-1)%ii,tim,nlev,nflag
  close(10)

!... output ctl file
  open(10, file = filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.ctl',  &
           form='formatted' )
  write(10,'(a)')'dset ^'//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.grd'
  write(10,'(a)')'dtype  station'
  write(10,'(a)')'stnmap '//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.map'
  write(10,'(a)')'undef -888888.'
  write(10,'(a)')'title '//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)
  write(10,'(a)')'tdef         1  linear 12z02aug2005  1hr '
!     write(10,'(a,i5)')'vars   ', raw%radar(nr)%radar_stn%levels*2
!     do k = 1, raw%radar(nr)%radar_stn%levels
  write(10,'(a,i5)')'vars   ', 2
  do k = 1, 1
     write(10,'(a3,i2.2,5x,a,f8.3,a)')'vr_',k, '0 99 radial velocity ',raw%radar(nr)%radar_stn%elv_ang(k),' elevation'
     write(10,'(a4,i2.2,4x,a,f8.3,a)')'hgt_',k,'0 99 height above SL ',raw%radar(nr)%radar_stn%elv_ang(k),' elevation'
  enddo
  write(10,'(a)')'endvars'
  close(10)

  enddo 

  return

end subroutine output_simulated_rv

!=======================================================================================
  subroutine output_verify_result ( filename, iobs, ob, fo3dm, obs_ii, obs_jj )

  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use map_utils
  use netcdf
  use wrf_tools

  implicit none

  character (len=*), intent(in)      :: filename
  integer, intent(in)                :: iobs
  real, dimension(iobs), intent(in)  :: ob, obs_ii, obs_jj
  real, dimension(iobs,3), intent(in):: fo3dm

  character (len=8) :: stid

  integer       :: i, nr, num, k, length, nlev, nflag
  real          :: tim

!------------------------------------------------------------------------------
! output in GrADS format, x: distance, y: azimuth, z: elevation angle levels,
!                         t: radar number;    rv: rv(x,y,z)

  length = len_trim( filename )

! each radar a file
  do nr = 1, raw%radar_stn_num

!... output grd file
  i = len_trim(raw%radar(nr)%radar_stn%idc)
  open(10, file = filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.grd',  &
           form='unformatted', access='direct', recl=10*4)
  do num = 1, iobs
     tim = 0.0
     nlev = 1
     nflag = 1
     write(stid,'(i8.8)')num
     write(10)stid,obs_jj(num),obs_ii(num),tim,nlev,nflag, &
              fo3dm(num,1:3)
  enddo
  nlev = 0
  write(10)stid,obs_jj(num),obs_ii(num),tim,nlev,nflag
  close(10)

!... output ctl file
  open(10, file = filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.ctl',  &
           form='formatted' )
  write(10,'(a)')'dset ^'//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.grd'
  write(10,'(a)')'dtype  station'
  write(10,'(a)')'stnmap '//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)//'.map'
  write(10,'(a)')'undef -888888.'
  write(10,'(a)')'title '//filename(1:length)//'_'//raw%radar(nr)%radar_stn%idc(1:i)
  write(10,'(a)')'tdef         1  linear 12z02aug2005  1hr '
  write(10,'(a,i5)')'vars   ', 3
  write(10,'(a)')'vr_en 0 99 pure ensemble forecast'
  write(10,'(a)')'vr_fc 0 99 enkf ensemble forecast'
  write(10,'(a)')'vr_up 0 99 enkf update'
  write(10,'(a)')'endvars'
  close(10)

  enddo

  return

end subroutine output_verify_result

!=======================================================================================
subroutine ideal_sounding_data( wrf_file, ix, jx, kx )

use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
!----------------------------------------------------------------------------
implicit none

character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx

integer                              :: i, j, k, nvar, unit
character (len=10)                   :: truthfile
real, dimension(ix+1, jx, kx  )      :: u
real, dimension(ix, jx+1, kx  )      :: v
real, dimension(ix, jx, kx  )        :: t, qv
real, dimension(ix, jx, kx+1)        :: ph, phb

!----------------------------------------------------------------------------
read(wrf_file(6:10),'(i5)')unit
write(truthfile,'(a5,i5.5)')wrf_file(1:5),unit-1

call get_variable3d( truthfile, 'U         ', ix+1, jx, kx, 1, u )
call get_variable3d( truthfile, 'V         ', ix, jx+1, kx, 1, v )
call get_variable3d( truthfile, 'T         ', ix, jx, kx, 1, t )
call get_variable3d( truthfile, 'QVAPOR    ', ix, jx, kx, 1, qv )
call get_variable3d( truthfile, 'PH        ', ix, jx, kx+1, 1, ph )
call get_variable3d( truthfile, 'PHB       ', ix, jx, kx+1, 1, phb )
qv = qv*1000.
ph = ph + phb

do k = gridobs_ks, gridobs_ke, gridobs_int_k
do j = gridobs_js, gridobs_je, gridobs_int_x
do i = gridobs_is, gridobs_ie, gridobs_int_x
do nvar = 1, 4   !U, V, Theta, Qv, P
   obs%num                 = obs%num + 1
   obs%position(obs%num,1) = i
   obs%position(obs%num,2) = j
   obs%position(obs%num,3) = k
   if ( nvar == 1 ) then
      obs%dat( obs%num ) = u(i,j,k)
      obs%type(obs%num) = 'idealU    '
      obs%err(obs%num) = 2.0
   elseif ( nvar == 2 ) then
      obs%dat( obs%num ) = v(i,j,k)
      obs%type(obs%num) = 'idealV    '
      obs%err(obs%num) = 2.0
   elseif ( nvar == 3 ) then
      obs%dat( obs%num ) = t(i,j,k)
      obs%type(obs%num) = 'idealPT   '
      obs%err(obs%num) = 1.0
   elseif ( nvar == 4 ) then
      obs%dat( obs%num ) = qv(i,j,k)
      obs%type(obs%num) = 'idealQV   '
      obs%err(obs%num) = 1.0
   elseif ( nvar == 5 ) then
      obs%dat( obs%num ) = ph(i,j,k)
      obs%type(obs%num) = 'idealPH   '
      obs%err(obs%num) = 150.0
   endif
enddo
enddo
enddo
enddo

  return

end subroutine ideal_sounding_data

!=======================================================================================
subroutine sort_sounding_data( wrf_file, ix, jx, kx, proj, instrument, datathin, hroi, vroi, grid_id)
use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
use map_utils
use wrf_tools
!----------------------------------------------------------------------------
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
integer, intent(in)                  :: datathin, hroi, vroi, grid_id
character (len=8), intent(in)       :: instrument
integer                              :: i, j, k, n, ista,sta, level,numlevs
real                                 :: x, y, gridu, gridv, trueu, truev
real, allocatable,dimension(:,:)     :: slp
integer, allocatable, dimension(:)   :: levels
real, allocatable, dimension(:,:,:,:):: data
real, allocatable, dimension(:)      :: obs_lat, obs_lon, obs_elv
character(len=12)                    :: fm
character(len=6)                     :: pfile
integer                              :: start_data, inter_data, iroi, ngxn


!. get data
sta = 0
numlevs = 0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)
   if( (instrument.eq.'sounding' .and. (fm(4:6).eq.'32 ' .or. fm(4:6).eq.'33 ' .or.  &
                                        fm(4:6).eq.'34 ' .or. fm(4:6).eq.'35 ' .or.  &
                                        fm(4:6).eq.'36 ' .or. fm(4:6).eq.'37 ' .or.  &
                                        fm(4:6).eq.'38 ') )      .or.                &
       (instrument.eq.'profiler' .and.  fm(4:6).eq.'132' )       .or.                &
       (instrument.eq.'aircft  ' .and. (fm(4:6).eq.'42 ' .or. fm(4:6).eq.'96 ' .or.  &
                                        fm(4:6).eq.'97 ') )      .or.                &
       (instrument.eq.'atovs   ' .and.  fm(4:6).eq.'131' )       .or.                &
       (instrument.eq.'satwnd  ' .and. fm(4:6).eq.'88 ' )      )then
         sta = sta + 1
         if(raw%gts%levels(n).gt.numlevs) numlevs=raw%gts%levels(n)
   endif
enddo
allocate( data(sta,numlevs,3,7) )
allocate(slp(sta,3))
allocate(levels(sta))
allocate(obs_lat(sta))
allocate(obs_lon(sta))
allocate(obs_elv(sta))

ista = 0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)

   if( (instrument.eq.'sounding' .and. (fm(4:6).eq.'32 ' .or. fm(4:6).eq.'33 ' .or.  &
                                        fm(4:6).eq.'34 ' .or. fm(4:6).eq.'35 ' .or.  &
                                        fm(4:6).eq.'36 ' .or. fm(4:6).eq.'37 ' .or.  &
                                        fm(4:6).eq.'38 ') )      .or.                &
       (instrument.eq.'profiler' .and.  fm(4:6).eq.'132' )       .or.                &
       (instrument.eq.'aircft  ' .and. (fm(4:6).eq.'42 ' .or. fm(4:6).eq.'96 ' .or.  &
                                        fm(4:6).eq.'97 ') )      .or.                &
       (instrument.eq.'atovs   ' .and.  fm(4:6).eq.'131' )       .or.                &
       (instrument.eq.'satwnd  ' .and. fm(4:6).eq.'88 ' )      )then

         ista = ista + 1
         obs_lat(ista) = raw%gts%latitude(n)
         obs_lon(ista) = raw%gts%longitude(n)
         obs_elv(ista) = raw%gts%elevation(n)
         slp(ista,1:3) = raw%gts%slp(n,1:3)
         levels(ista) = raw%gts%levels(n)
         do k =1, levels(ista)
            do i = 1, 3
               data(ista,k,i,1) = raw%gts%pres(n,k,i)
               data(ista,k,i,2) = raw%gts%spd(n,k,i)
               data(ista,k,i,3) = raw%gts%wd(n,k,i)
               data(ista,k,i,4) = raw%gts%height(n,k,i)
               data(ista,k,i,5) = raw%gts%t(n,k,i)
               data(ista,k,i,6) = raw%gts%td(n,k,i)
            enddo
!................ convert rh to q and calculate q_error from rh_error
!                 data(ista,k,i,7) = raw%gts%rh(n,k,i)
               call rel_humidity_to_q(raw%gts%pres(n,k,1), raw%gts%t(n,k,1), raw%gts%rh(n,k,1), data(ista,k,1,7), &
                                         raw%gts%t(n,k,3), raw%gts%rh(n,k,3), data(ista,k,3,7) ) 
               data(ista,k,2,7) = raw%gts%rh(n,k,2)
         enddo
   endif
enddo

start_data = 1
inter_data = 1
if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)

if(print_detail > 20) then
write(*,*)'start_data, ista, inter_data =',start_data, ista, inter_data
end if
iroi = 0

!----------------------------------------------------------------------------
!. data thinning   
! JPOTERJOY: perform thinning on vertical levels, not soundings
!   do_reports : do n = start_data, ista, inter_data
do_reports : do n = start_data, ista

  iroi = iroi + 1
  call cal_hroi ( instrument, grid_id, iroi, ngxn )
  call latlon_to_ij( proj, obs_lat(n), obs_lon(n), x, y ) 
  if(print_detail > 20) write(*,'(i5,2x,a,4f10.3)')iroi,instrument, obs_lat(n), obs_lon(n), x, y
  if ( x.lt.1. .or. x.gt.real(ix) .or. y.lt.1. .or. y.gt.real(jx) ) cycle do_reports

!... slp
  if ( slp(n,1) .ge. 80000. .and. slp(n,1).lt.120000. ) then
       obs%num                 = obs%num + 1
       obs%dat     (obs%num  ) = slp(n,1)
       obs%type    (obs%num  ) = 'slp       '
       obs%err     (obs%num  ) = slp(n,3)
       obs%position(obs%num,1) = x
       obs%position(obs%num,2) = y
       obs%position(obs%num,3) = 1.
       obs%roi     (obs%num,1) = hroi * ngxn
       obs%roi     (obs%num,2) = vroi
   endif

  if ( instrument.eq.'sounding' ) then
!... SurfacePressure
   do k = 1, 1
      if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .lt. 105000. .and.  &
          abs(data(n,k,2,1)).lt. 5. .and. abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. ) then
          obs%num                 = obs%num + 1
          obs%dat     (obs%num  ) = data(n,k,1,1)
          obs%type (obs%num  ) = 'S'//instrument//'P'
          obs%err     (obs%num  ) = data(n,k,3,1)
          obs%position(obs%num,1) = x
          obs%position(obs%num,2) = y
          obs%position(obs%num,4) = data(n,k,1,1)
          obs%sta     (obs%num,1) = obs_elv(n)
          obs%sta     (obs%num,2) = data(n,1,1,1)
          obs%sta     (obs%num,3) = data(n,1,1,5)
          obs%sta     (obs%num,4) = data(n,1,1,7)
          obs%roi     (obs%num,1) = hroi * ngxn
          obs%roi     (obs%num,2) = vroi
      endif
   enddo
   endif

!..... T
! JPOTERJOY: perform thinning on vertical levels, not soundings
!      do k = 1, levels(n)
   do k = 1, levels(n), inter_data
      if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .le. 100000. ) then     !
          if ( data(n,k,1,5).ge.100. .and. data(n,k,1,5) .le. 350. .and. abs(data(n,k,2,5)).lt. 90. ) then 
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,5)
               obs%type    (obs%num  ) = 'P'//instrument//'T'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'T' 
               obs%err     (obs%num  ) = data(n,k,3,5)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,5).ge.100. .and. data(n,k,1,5) .le. 350. .and. abs(data(n,k,2,5)).lt. 90. ) then
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,5)
               obs%type    (obs%num  ) = 'H'//instrument//'T'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'T'
               obs%err     (obs%num  ) = data(n,k,3,5)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      endif
   enddo

!................ u, v
! JPOTERJOY: perform thinning on vertical levels, not soundings
!      do k = 1, levels(n)
   do k = 1, levels(n), inter_data
     if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .lt. 105000. ) then
          if ( data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200. .and.   &
               data(n,k,1,3).ge.0. .and. data(n,k,1,3).lt.360. .and.   & 
               abs(data(n,k,2,2)).lt. 90.  ) then 
               call dir_spd2xy(data(n,k,1,3), data(n,k,1,2), trueu, truev)
               call truewind_to_gridwind(obs_lon(n), proj, trueu, truev, gridu, gridv)
!................ u
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridu
               obs%type    (obs%num  ) = 'P'//instrument//'U'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'U' 
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
!................ v
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridv
               obs%type    (obs%num  ) = 'P'//instrument//'V'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'V' 
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

          if ( (data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200.) .and.   &
               (data(n,k,1,3).lt.0. .or.  data(n,k,1,3).gt.360.) .and.   &
                abs(data(n,k,2,2)).lt. 90.  ) then
!................ wind speed
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,2)
               obs%type    (obs%num  ) = 'P'//instrument//'S'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'S'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
               

     else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200. .and.   &
               data(n,k,1,3).ge.0. .and. data(n,k,1,3).lt.360. .and.   &
               abs(data(n,k,2,2)).lt. 90.  ) then
               call dir_spd2xy(data(n,k,1,3), data(n,k,1,2), trueu, truev)
               call truewind_to_gridwind(obs_lon(n), proj, trueu, truev, gridu, gridv)
!................ u
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridu
               obs%type    (obs%num  ) = 'H'//instrument//'U'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'U'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
!................ v
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridv
               obs%type    (obs%num  ) = 'H'//instrument//'V'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'V'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

          if ( (data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200.) .and.   &
               (data(n,k,1,3).lt.0. .or.  data(n,k,1,3).gt.360.) .and.   &
                abs(data(n,k,2,2)).lt. 90.  ) then
!................ wind speed
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,2)
               obs%type    (obs%num  ) = 'H'//instrument//'S'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'S'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

     endif
   enddo

!................ Q(g/kg)
! JPOTERJOY: perform thinning on vertical levels, not soundings
!      do k = 1, levels(n)
   do k = 1, levels(n), inter_data
      if ( data(n,k,1,1) .gt. 25000. .and. data(n,k,1,1) .lt. 100000. ) then     !
          if ( data(n,k,1,7).ge.0.01 .and. data(n,k,1,7) .le. 100. .and. abs(data(n,k,2,7)).lt. 5.) then 
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,7)
               obs%type    (obs%num  ) = 'P'//instrument//'Q'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'Q' 
               obs%err     (obs%num  ) = data(n,k,3,7)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,7).ge.0.01 .and. data(n,k,1,7) .le. 100. .and. abs(data(n,k,2,7)).lt. 5.) then
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,7)
               obs%type    (obs%num  ) = 'H'//instrument//'Q'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'Q'
               obs%err     (obs%num  ) = data(n,k,3,7)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      endif
   enddo

  end do do_reports

  return

end subroutine sort_sounding_data
!=======================================================================================
subroutine sort_upperair_data( wrf_file, ix, jx, kx, proj, instrument, datathin, hroi, vroi, grid_id)
use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
use map_utils
use wrf_tools
!----------------------------------------------------------------------------
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
integer, intent(in)                  :: datathin, hroi, vroi, grid_id
character (len=8), intent(in)       :: instrument
integer                              :: i, j, k, n, ista,sta, level,numlevs
real                                 :: x, y, gridu, gridv, trueu, truev
real, allocatable,dimension(:,:)   :: slp
integer, allocatable,dimension(:)   :: levels
real, allocatable, dimension(:,:,:,:) :: data
real, allocatable, dimension(:)     :: obs_lat, obs_lon, obs_elv
character(len=12)                    :: fm
character(len=6)                     :: pfile
integer                              :: start_data, inter_data, iroi, ngxn

!. get data
sta = 0
numlevs=0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)
   if( (instrument.eq.'sounding' .and. (fm(4:6).eq.'32 ' .or. fm(4:6).eq.'33 ' .or.  &
                                        fm(4:6).eq.'34 ' .or. fm(4:6).eq.'35 ' .or.  &
                                        fm(4:6).eq.'36 ' .or. fm(4:6).eq.'37 ' .or.  &
                                        fm(4:6).eq.'38 ') )      .or.                &
       (instrument.eq.'profiler' .and.  fm(4:6).eq.'132' )       .or.                &
       (instrument.eq.'aircft  ' .and. (fm(4:6).eq.'42 ' .or. fm(4:6).eq.'96 ' .or.  &
                                        fm(4:6).eq.'97 ') )      .or.                &
       (instrument.eq.'atovs   ' .and.  fm(4:6).eq.'131' )       .or.                &
       (instrument.eq.'satwnd  ' .and. fm(4:6).eq.'88 ' )      )then
         sta = sta + 1
         if(raw%gts%levels(n).gt.numlevs) numlevs=raw%gts%levels(n)
   endif
enddo
allocate( data(sta,numlevs,3,7) )
allocate(slp(sta,3))
allocate(levels(sta))
allocate(obs_lat(sta))
allocate(obs_lon(sta))
allocate(obs_elv(sta))

ista = 0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)

   if( (instrument.eq.'sounding' .and. (fm(4:6).eq.'32 ' .or. fm(4:6).eq.'33 ' .or.  &
                                        fm(4:6).eq.'34 ' .or. fm(4:6).eq.'35 ' .or.  &
                                        fm(4:6).eq.'36 ' .or. fm(4:6).eq.'37 ' .or.  &
                                        fm(4:6).eq.'38 ') )      .or.                &
       (instrument.eq.'profiler' .and.  fm(4:6).eq.'132' )       .or.                &
       (instrument.eq.'aircft  ' .and. (fm(4:6).eq.'42 ' .or. fm(4:6).eq.'96 ' .or.  &
                                        fm(4:6).eq.'97 ') )      .or.                &
       (instrument.eq.'atovs   ' .and.  fm(4:6).eq.'131' )       .or.                &
       (instrument.eq.'satwnd  ' .and. fm(4:6).eq.'88 ' )      )then

         ista = ista + 1
         obs_lat(ista) = raw%gts%latitude(n)
         obs_lon(ista) = raw%gts%longitude(n)
         obs_elv(ista) = raw%gts%elevation(n)
         slp(ista,1:3) = raw%gts%slp(n,1:3)
         levels(ista) = raw%gts%levels(n)
         do k =1, levels(ista)
            do i = 1, 3
               data(ista,k,i,1) = raw%gts%pres(n,k,i)
               data(ista,k,i,2) = raw%gts%spd(n,k,i)
               data(ista,k,i,3) = raw%gts%wd(n,k,i)
               data(ista,k,i,4) = raw%gts%height(n,k,i)
               data(ista,k,i,5) = raw%gts%t(n,k,i)
               data(ista,k,i,6) = raw%gts%td(n,k,i)
            enddo
!................ convert rh to q and calculate q_error from rh_error
!                 data(ista,k,i,7) = raw%gts%rh(n,k,i)
               call rel_humidity_to_q(raw%gts%pres(n,k,1), raw%gts%t(n,k,1), raw%gts%rh(n,k,1), data(ista,k,1,7), &
                                         raw%gts%t(n,k,3), raw%gts%rh(n,k,3), data(ista,k,3,7) ) 
               data(ista,k,2,7) = raw%gts%rh(n,k,2)
         enddo
   endif
enddo

!----------------------------------------------------------------------------
!. data thinning   
start_data = 1
inter_data = 1
if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)

if(print_detail > 20) then
write(*,*)'start_data, ista, inter_data =',start_data, ista, inter_data
end if
iroi = 0
do_reports : do n = start_data, ista, inter_data

  iroi = iroi + 1
  call cal_hroi ( instrument, grid_id, iroi, ngxn )
  call latlon_to_ij( proj, obs_lat(n), obs_lon(n), x, y ) 
  if(print_detail > 20) write(*,'(i5,2x,a,4f10.3)')iroi,instrument, obs_lat(n), obs_lon(n), x, y
  if ( x.lt.1. .or. x.gt.real(ix) .or. y.lt.1. .or. y.gt.real(jx) ) cycle do_reports

!... slp
  if ( slp(n,1) .ge. 80000. .and. slp(n,1).lt.120000. ) then
       obs%num                 = obs%num + 1
       obs%dat     (obs%num  ) = slp(n,1)
       obs%type    (obs%num  ) = 'slp       '
       obs%err     (obs%num  ) = slp(n,3)
       obs%position(obs%num,1) = x
       obs%position(obs%num,2) = y
       obs%position(obs%num,3) = 1.
       obs%roi     (obs%num,1) = hroi * ngxn
       obs%roi     (obs%num,2) = vroi
   endif

  if ( instrument.eq.'sounding' ) then
!... SurfacePressure
   do k = 1, 1
      if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .lt. 105000. .and.  &
          abs(data(n,k,2,1)).lt. 5. .and. abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. ) then
          obs%num                 = obs%num + 1
          obs%dat     (obs%num  ) = data(n,k,1,1)
          obs%type (obs%num  ) = 'S'//instrument//'P'
          obs%err     (obs%num  ) = data(n,k,3,1)
          obs%position(obs%num,1) = x
          obs%position(obs%num,2) = y
          obs%position(obs%num,4) = data(n,k,1,1)
          obs%sta     (obs%num,1) = obs_elv(n)
          obs%sta     (obs%num,2) = data(n,1,1,1)
          obs%sta     (obs%num,3) = data(n,1,1,5)
          obs%sta     (obs%num,4) = data(n,1,1,7)
          obs%roi     (obs%num,1) = hroi * ngxn
          obs%roi     (obs%num,2) = vroi
      endif
   enddo
   endif

!..... T
   do k = 1, levels(n)
      if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .le. 100000. ) then     !
          if ( data(n,k,1,5).ge.100. .and. data(n,k,1,5) .le. 350. .and. abs(data(n,k,2,5)).lt. 90. ) then 
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,5)
               obs%type    (obs%num  ) = 'P'//instrument//'T'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'T' 
               obs%err     (obs%num  ) = data(n,k,3,5)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,5).ge.100. .and. data(n,k,1,5) .le. 350. .and. abs(data(n,k,2,5)).lt. 90. ) then
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,5)
               obs%type    (obs%num  ) = 'H'//instrument//'T'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'T'
               obs%err     (obs%num  ) = data(n,k,3,5)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      endif
   enddo

!................ u, v
   do k = 1, levels(n)
     if ( data(n,k,1,1) .gt. 100. .and. data(n,k,1,1) .lt. 105000. ) then
          if ( data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200. .and.   &
               data(n,k,1,3).ge.0. .and. data(n,k,1,3).lt.360. .and.   & 
               abs(data(n,k,2,2)).lt. 90.  ) then 
               call dir_spd2xy(data(n,k,1,3), data(n,k,1,2), trueu, truev)
               call truewind_to_gridwind(obs_lon(n), proj, trueu, truev, gridu, gridv)
!................ u
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridu
               obs%type    (obs%num  ) = 'P'//instrument//'U'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'U' 
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
!................ v
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridv
               obs%type    (obs%num  ) = 'P'//instrument//'V'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'V' 
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

          if ( (data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200.) .and.   &
               (data(n,k,1,3).lt.0. .or.  data(n,k,1,3).gt.360.) .and.   &
                abs(data(n,k,2,2)).lt. 90.  ) then
!................ wind speed
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,2)
               obs%type    (obs%num  ) = 'P'//instrument//'S'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'S'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
               

     else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200. .and.   &
               data(n,k,1,3).ge.0. .and. data(n,k,1,3).lt.360. .and.   &
               abs(data(n,k,2,2)).lt. 90.  ) then
               call dir_spd2xy(data(n,k,1,3), data(n,k,1,2), trueu, truev)
               call truewind_to_gridwind(obs_lon(n), proj, trueu, truev, gridu, gridv)
!................ u
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridu
               obs%type    (obs%num  ) = 'H'//instrument//'U'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'U'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
!................ v
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = gridv
               obs%type    (obs%num  ) = 'H'//instrument//'V'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'V'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

          if ( (data(n,k,1,2).ge.0. .and. data(n,k,1,2).le.200.) .and.   &
               (data(n,k,1,3).lt.0. .or.  data(n,k,1,3).gt.360.) .and.   &
                abs(data(n,k,2,2)).lt. 90.  ) then
!................ wind speed
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,2)
               obs%type    (obs%num  ) = 'H'//instrument//'S'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'S'
               obs%err     (obs%num  ) = data(n,k,3,2)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif

     endif
   enddo

!................ Q(g/kg)
   do k = 1, levels(n)
      if ( data(n,k,1,1) .gt. 25000. .and. data(n,k,1,1) .lt. 100000. ) then     !
          if ( data(n,k,1,7).ge.0.01 .and. data(n,k,1,7) .le. 100. .and. abs(data(n,k,2,7)).lt. 5.) then 
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,7)
               obs%type    (obs%num  ) = 'P'//instrument//'Q'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'Q' 
               obs%err     (obs%num  ) = data(n,k,3,7)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,1)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      else if ( data(n,k,1,1) .lt. 0. .or. data(n,k,1,1) .gt. 200000. ) then
          if ( data(n,k,1,7).ge.0.01 .and. data(n,k,1,7) .le. 100. .and. abs(data(n,k,2,7)).lt. 5.) then
               obs%num                 = obs%num + 1
               obs%dat     (obs%num  ) = data(n,k,1,7)
               obs%type    (obs%num  ) = 'H'//instrument//'Q'
               if( abs(data(n,k,1,4)-obs_elv(n)) .le. 5. .and. obs_elv(n).gt.50. )obs%type(obs%num  ) = 'S'//instrument//'Q'
               obs%err     (obs%num  ) = data(n,k,3,7)
               obs%position(obs%num,1) = x
               obs%position(obs%num,2) = y
               obs%position(obs%num,4) = data(n,k,1,4)
               obs%sta     (obs%num,1) = obs_elv(n)
               obs%sta     (obs%num,2) = data(n,1,1,1)
               obs%sta     (obs%num,3) = data(n,1,1,5)
               obs%sta     (obs%num,4) = data(n,1,1,7)
               obs%roi     (obs%num,1) = hroi * ngxn
               obs%roi     (obs%num,2) = vroi
          endif
      endif
   enddo

  end do do_reports

  return

end subroutine sort_upperair_data
!=======================================================================================
subroutine sort_surface_data( wrf_file, ix, jx, kx, proj, instrument, datathin, hroi, vroi, grid_id )
use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
use map_utils
use wrf_tools
!----------------------------------------------------------------------------
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
integer, intent(in)                  :: datathin, hroi, vroi, grid_id
character (len=8), intent(in)       :: instrument

integer                              :: i, j, k, n, ista, sta

real                                 :: x, y, trueu, truev, gridu, gridv
real, allocatable,dimension(:,:)     :: slp,pw
real, allocatable,dimension(:,:,:)   :: data
real, allocatable,dimension(:)       :: obs_lat, obs_lon, obs_elv
character(len=12)                    :: fm

integer                              :: start_data, inter_data, iroi, ngxn

!----------------------------------------------------------------------------
!. get data
sta = 0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)
   if( (instrument.eq.'surface ' .and. (fm(4:6).eq.'12 ' .or. fm(4:6).eq.'14 ') ) .or. &
       (instrument.eq.'metar   ' .and. (fm(4:6).eq.'15 ' .or. fm(4:6).eq.'16 ') ) .or. &
       (instrument.eq.'sfcshp  ' .and. (fm(4:6).eq.'13 ' .or. fm(4:6).eq.'18 ') ) .or. &
       (instrument.eq.'spssmi  ' .and. (fm(4:6).eq.'125' .or. fm(4:6).eq.'126') ) .or. &
       (instrument.eq.'gpspw   ' .and. (fm(4:6).eq.'111') ) )then
         sta = sta + 1
   endif
enddo
allocate(slp(sta,3))
allocate(pw(sta,3))
allocate(data(sta,3,7))
allocate(obs_lat(sta))
allocate(obs_lon(sta))
allocate(obs_elv(sta))

ista = 0
do n = 1, raw%gts%num
   fm = raw%gts%platform(n)

   if( (instrument.eq.'surface ' .and. (fm(4:6).eq.'12 ' .or. fm(4:6).eq.'14 ') ) .or. &
       (instrument.eq.'metar   ' .and. (fm(4:6).eq.'15 ' .or. fm(4:6).eq.'16 ') ) .or. &
       (instrument.eq.'sfcshp  ' .and. (fm(4:6).eq.'13 ' .or. fm(4:6).eq.'18 ') ) .or. &
       (instrument.eq.'spssmi  ' .and. (fm(4:6).eq.'125' .or. fm(4:6).eq.'126') ) .or. &
       (instrument.eq.'gpspw   ' .and. (fm(4:6).eq.'111') ) )then

         ista = ista + 1
         obs_lat(ista) = raw%gts%latitude(n)
         obs_lon(ista) = raw%gts%longitude(n)
         obs_elv(ista) = raw%gts%elevation(n)
         do i = 1, 3
            slp(ista,i) = raw%gts%slp(n,i)
            pw(ista,i) = raw%gts%pw(n,i)
            data(ista,i,1) = raw%gts%pres(n,1,i)
            data(ista,i,2) = raw%gts%spd(n,1,i)
            data(ista,i,3) = raw%gts%wd(n,1,i)
            data(ista,i,4) = raw%gts%height(n,1,i)
            data(ista,i,5) = raw%gts%t(n,1,i)
            data(ista,i,6) = raw%gts%td(n,1,i)
         enddo
!.......... convert rh to q and calculate q_error from rh_error
!           data(ista,i,7) = raw%gts%rh(n,1,i)
         call rel_humidity_to_q(raw%gts%pres(n,1,1), raw%gts%t(n,1,1), raw%gts%rh(n,1,1), data(ista,1,7), &
                                 raw%gts%t(n,1,3), raw%gts%rh(n,1,3), data(ista,3,7) ) 
         data(ista,2,7) = raw%gts%rh(n,1,2)
   endif
enddo

!----------------------------------------------------------------------------
!. data thinning   
start_data = 1
inter_data = 1
if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)

iroi = 0
do_reports : do n = start_data, ista, inter_data

  iroi = iroi + 1
  call cal_hroi ( instrument, grid_id, iroi, ngxn )
  call latlon_to_ij( proj, obs_lat(n), obs_lon(n), x, y ) 
  if ( x.lt.1. .or. x.gt.real(ix) .or. y.lt.1. .or. y.gt.real(jx) ) cycle do_reports

!... slp
  if ( slp(n,1) .ge. 80000. .and. slp(n,1).lt.120000. .and. abs(slp(n,2)).lt. 5. ) then
       obs%num                 = obs%num + 1
       obs%dat     (obs%num  ) = slp(n,1)
       obs%type    (obs%num  ) = 'slp       '
       obs%err     (obs%num  ) = slp(n,3)
       obs%position(obs%num,1) = x
       obs%position(obs%num,2) = y
       obs%position(obs%num,3) = 1.
       obs%sta     (obs%num,1) = obs_elv(n)
       obs%sta     (obs%num,2) = data(n,1,1)
       obs%sta     (obs%num,3) = data(n,1,5)
       obs%sta     (obs%num,4) = data(n,1,7)
       obs%roi     (obs%num,1) = hroi * ngxn
       obs%roi     (obs%num,2) = vroi
   endif

!... pw
   if ( pw(n,1) .ge. 0.0 .and. pw(n,1) .lt. 100.0 .and. abs(pw(n,2)).lt.5.) then
       obs%num                 = obs%num + 1
       obs%dat     (obs%num  ) = pw(n,1)
       obs%type    (obs%num  ) = 'pw        '
       obs%err     (obs%num  ) = pw(n,3)
       obs%position(obs%num,1) = x
       obs%position(obs%num,2) = y
       obs%position(obs%num,3) = 1.
       obs%sta     (obs%num,1) = obs_elv(n)
       obs%sta     (obs%num,2) = data(n,1,1)
       obs%sta     (obs%num,3) = data(n,1,5)
       obs%sta     (obs%num,4) = data(n,1,7)
       obs%roi     (obs%num,1) = hroi * ngxn
       obs%roi     (obs%num,2) = vroi
   endif

!......... T
   if ( data(n,1,5).ge.100. .and. data(n,1,5) .le. 350. .and. abs(data(n,2,5)).lt. 5. ) then 
        obs%num                 = obs%num + 1
        obs%dat     (obs%num  ) = data(n,1,5)
        obs%type    (obs%num  ) = 'S'//instrument//'T'
        obs%err     (obs%num  ) = data(n,3,5)
        obs%position(obs%num,1) = x
        obs%position(obs%num,2) = y
        obs%position(obs%num,3) = 1.
        obs%sta     (obs%num,1) = obs_elv(n)
        obs%sta     (obs%num,2) = data(n,1,1)
        obs%sta     (obs%num,3) = data(n,1,5)
        obs%sta     (obs%num,4) = data(n,1,7)
        obs%roi     (obs%num,1) = hroi * ngxn
        obs%roi     (obs%num,2) = vroi
   endif

   if ( data(n,1,2).ge.0. .and. data(n,1,2).le.200. .and.   &
        data(n,1,3).ge.0. .and. data(n,1,3).lt.360. .and.   &
        abs(data(n,2,2)).lt. 5. ) then
        call dir_spd2xy(data(n,1,3), data(n,1,2), trueu, truev) 
        call truewind_to_gridwind(obs_lon(n), proj, trueu, truev, gridu, gridv) 
!......... u
        obs%num                 = obs%num + 1
        obs%dat     (obs%num  ) = gridu
        obs%type    (obs%num  ) = 'S'//instrument//'U'
        obs%err     (obs%num  ) = data(n,3,2)
        obs%position(obs%num,1) = x
        obs%position(obs%num,2) = y
        obs%position(obs%num,3) = 1
        obs%sta     (obs%num,1) = obs_elv(n)
        obs%sta     (obs%num,2) = data(n,1,1)
        obs%sta     (obs%num,3) = data(n,1,5)
        obs%sta     (obs%num,4) = data(n,1,7)
        obs%roi     (obs%num,1) = hroi * ngxn
        obs%roi     (obs%num,2) = vroi
!......... v
        obs%num                 = obs%num + 1
        obs%dat     (obs%num  ) = gridv
        obs%type    (obs%num  ) = 'S'//instrument//'V'
        obs%err     (obs%num  ) = data(n,3,2)
        obs%position(obs%num,1) = x
        obs%position(obs%num,2) = y
        obs%position(obs%num,3) = 1.
        obs%sta     (obs%num,1) = obs_elv(n)
        obs%sta     (obs%num,2) = data(n,1,1)
        obs%sta     (obs%num,3) = data(n,1,5)
        obs%sta     (obs%num,4) = data(n,1,7)
        obs%roi     (obs%num,1) = hroi * ngxn
        obs%roi     (obs%num,2) = vroi
   endif

!......... Q(g/kg)
   if ( data(n,1,1) .gt. 25000. .and. data(n,1,1) .lt. 100000. ) then     !
     if ( data(n,1,7).ge.0.01 .and. data(n,1,7) .le. 100. .and. abs(data(n,2,7)).lt. 5.) then 
        obs%num                 = obs%num + 1
        obs%dat     (obs%num  ) = data(n,1,7)
        obs%type    (obs%num  ) = 'S'//instrument//'Q'
        obs%err     (obs%num  ) = data(n,3,7)
        obs%position(obs%num,1) = x
        obs%position(obs%num,2) = y
        obs%position(obs%num,3) = 1.
        obs%sta     (obs%num,1) = obs_elv(n)
        obs%sta     (obs%num,2) = data(n,1,1)
        obs%sta     (obs%num,3) = data(n,1,5)
        obs%sta     (obs%num,4) = data(n,1,7)
        obs%roi     (obs%num,1) = hroi * ngxn
        obs%roi     (obs%num,2) = vroi
     endif
   endif

end do do_reports

  return

end subroutine sort_surface_data


!=======================================================================================
subroutine sort_radarRV_data( wrf_file, ix, jx, kx, proj, instrument, datathin, hroi, vroi, grid_id )

use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
use map_utils
!----------------------------------------------------------------------------
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
integer, intent(in)                  :: datathin, hroi, vroi, grid_id
character (len=10), intent(in)       :: instrument

integer                              :: i, j, k, nr, n, ista, sta

real                                 :: x, y, u, v
real, allocatable, dimension(:,:)    :: data

integer                              :: start_data, inter_data, iroi, ngxn

!----------------------------------------------------------------------------
!. get data
sta = 0
if ( instrument == 'RadarRV   ' )then
   do nr = 1, raw%radar_stn_num
   do n  = 1, raw%radar(nr)%radar_stn%numObs
      sta = sta + 1
   enddo
   enddo
else if ( instrument == 'AirborneRV' )then
   do nr = 1, raw%airborne%num
      sta = sta + 1
   enddo
endif
allocate(data(sta,8))

ista = 0
if ( instrument == 'RadarRV   ' )then
   do nr = 1, raw%radar_stn_num
   do n  = 1, raw%radar(nr)%radar_stn%numObs
      ista = ista + 1
       data(ista,1) = raw%radar(nr)%rv(n)
       data(ista,2) = raw%radar(nr)%radar_stn%err_rv
       data(ista,3) = raw%radar(nr)%radar_point(n)%ii
       data(ista,4) = raw%radar(nr)%radar_point(n)%jj
       data(ista,5) = raw%radar(nr)%radar_point(n)%hgt
       data(ista,6) = raw%radar(nr)%radar_stn%i_radar
       data(ista,7) = raw%radar(nr)%radar_stn%j_radar
       data(ista,8) = raw%radar(nr)%radar_stn%elv
   enddo
   enddo
else if ( instrument == 'AirborneRV' )then
   do nr = 1, raw%airborne%num
      ista = ista + 1
       data(ista,1) = raw%airborne%rv(nr)
       data(ista,2) = 3.0
       data(ista,3) = raw%airborne%ii(nr)
       data(ista,4) = raw%airborne%jj(nr)
       data(ista,5) = raw%airborne%hh(nr)
       data(ista,6) = raw%airborne%radar_ii(nr)
       data(ista,7) = raw%airborne%radar_jj(nr)
       data(ista,8) = raw%airborne%height(nr)
   enddo
endif

!----------------------------------------------------------------------------
!. data thinning
start_data = 1
inter_data = 1
if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)

iroi = 0
do_reports : do n = start_data, ista, inter_data

  iroi = iroi + 1
  call cal_hroi ( 'Radar   ', grid_id, iroi, ngxn )
  x = data(n,3)
  y = data(n,4)
  if(x.le.1. .and. x.ge.real(ix) .and. y.le.1. .and. y.ge.real(jx) ) cycle do_reports
  obs%num                 = obs%num + 1
  obs%dat( obs%num )      = data(n,1) 
  obs%type(obs%num)       = 'RadarRV   '
  obs%err(obs%num)        = data(n,2)
  obs%position(obs%num,1) = data(n,3)
  obs%position(obs%num,2) = data(n,4)
  obs%position(obs%num,4) = data(n,5)
  obs%sta(obs%num,1)      = data(n,6)
  obs%sta(obs%num,2)      = data(n,7)
  obs%sta(obs%num,4)      = data(n,8)
  obs%roi(obs%num,1)      = hroi*ngxn
  obs%roi(obs%num,2)      = vroi

end do do_reports

!   deallocate ( raw%radar )
!   deallocate ( raw%airborne%rv )
!   deallocate ( raw%airborne%ii )
!   deallocate ( raw%airborne%jj )
!   deallocate ( raw%airborne%hh )
!   deallocate ( raw%airborne%radar_ii )
!   deallocate ( raw%airborne%radar_jj )
!   deallocate ( raw%airborne%height )

  return

end subroutine sort_radarRV_data   

!========================================================================================


!---edited by Minamide 2014.11.18 (for satellite data)
!=======================================================================================
subroutine sort_radiance_data( wrf_file, ix, jx, kx, proj, instrument, datathin,hroi, vroi, grid_id )

use constants
use namelist_define
use mpi_module
use obs_define
use netcdf
use map_utils
!----------------------------------------------------------------------------
implicit none

type(proj_info)                      :: proj
character (len=10), intent(in)       :: wrf_file
integer, intent(in)                  :: ix, jx, kx
integer, intent(in)                  :: datathin, hroi, vroi, grid_id
character (len=10), intent(in)       :: instrument

integer                              :: i, j, k, nr, n, ista, sta

real                                 :: x, y, u, v
real, allocatable, dimension(:,:)    :: data
character (len=12), allocatable, dimension(:) :: data_sat
integer, allocatable, dimension(:)   :: data_ch, data_hroi, data_hroi_d


integer                              :: start_data, inter_data, iroi, ngxn

!----------------------------------------------------------------------------
!. get data

sta = 0
do nr = 1, raw%radiance%num
   sta = sta + 1
enddo
allocate(data(sta,4))
allocate(data_sat(sta))
allocate(data_ch(sta))
allocate(data_hroi(sta))
allocate(data_hroi_d(sta))

ista = 0
do nr = 1, raw%radiance%num
    ista = ista + 1
    data(ista,1) = raw%radiance%tb(nr)
    data(ista,2) = raw%radiance%err(nr)
    data(ista,3) = raw%radiance%ii(nr)
    data(ista,4) = raw%radiance%jj(nr)
    data_sat(ista) = raw%radiance%platform(nr)
    data_ch(ista)  = raw%radiance%ch(nr)
    data_hroi(ista)  = raw%radiance%hroi(nr)
    data_hroi_d(ista)  = raw%radiance%hroi_d(nr)
enddo
!----------------------------------------------------------------------------
!. data thinning
start_data = 1
inter_data = 1
if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)

iroi = 0
!inter_data does not work  2014.11.20 Minamide
do_reports : do n = start_data, ista,inter_data
  iroi = iroi + 1
  call cal_hroi ( instrument, grid_id, iroi, ngxn )
  x = data(n,3)
  y = data(n,4)

  if(x.le.1. .and. x.ge.real(ix) .and. y.le.1. .and. y.ge.real(jx) ) cycle do_reports
  obs%num                 = obs%num + 1
  obs%dat( obs%num )      = data(n,1)
  obs%type(obs%num)       = 'Radiance  '
  obs%err(obs%num)        = data(n,2)
  obs%position(obs%num,1) = data(n,3)
  obs%position(obs%num,2) = data(n,4)
  obs%sat(obs%num)        = data_sat(n)
  obs%ch(obs%num)         = data_ch(n)
  obs%roi(obs%num,1)      = data_hroi(n)*ngxn
  obs%roi(obs%num,3)      = data_hroi_d(n)*ngxn
!  obs%roi(obs%num,1)      = 1
  obs%roi(obs%num,2)      = vroi
end do do_reports

  return

end subroutine sort_radiance_data
!========================================================================================

subroutine bcast_obs(root)

!----------------------------------------------------------------------------
!  PURPOSE: broadcast obs structure from root to all other processors
!        10/03/2012 Jon Poterjoy
!----------------------------------------------------------------------------
use constants
use mpi_module
use obs_define
!----------------------------------------------------------------------------
implicit none

integer, intent(in)               :: root
integer                           :: m, n
real, allocatable, dimension(:)   :: buff     !!! array buffer
character(len=10)                 :: buff2    !!! string buffer

  ! Send data from root
call MPI_BCAST(obs%num,          1, MPI_INTEGER, root, comm, ierr)
allocate( buff  (obs%num) )

!. allocate and initialize arrays
if ( my_proc_id .ne. root ) then
   allocate( obs%dat      ( obs%num    ) )
   allocate( obs%type     ( obs%num    ) )
   allocate( obs%err      ( obs%num    ) )
   allocate( obs%position ( obs%num, 4 ) )
   allocate( obs%sta      ( obs%num, 4 ) )
   allocate( obs%roi      ( obs%num, 2 ) )
   allocate( obs%sat      ( obs%num    ) )
   allocate( obs%ch       ( obs%num    ) )

   obs%dat        = -888888.
   obs%type       = '          '
   obs%err        = -888888.
   obs%position   = -888888.
   obs%sta        = -888888.
   obs%roi        = -888888
   obs%sat        = '            '
   obs%ch         = -888888
endif

if ( my_proc_id == root ) buff = obs%dat
call MPI_BCAST(buff,      obs%num,    MPI_REAL, root, comm, ierr)
obs%dat = buff

do m = 1, obs%num
   if ( my_proc_id == root ) buff2 = obs%type(m)
   call MPI_BCAST(buff2,      10,    MPI_CHARACTER, root, comm, ierr)
   obs%type(m) = buff2
enddo

if ( my_proc_id == root ) buff = obs%err
call MPI_BCAST(buff,      obs%num,    MPI_REAL, root, comm, ierr)
obs%err = buff
 
do m = 1, 4
   if ( my_proc_id == root ) buff = obs%position(:,m)
   call MPI_BCAST(buff,      obs%num,    MPI_REAL, root, comm, ierr)
   do n = 1, obs%num
     obs%position(n,m) = buff(n)
   enddo
enddo
 
do m = 1, 4
   if ( my_proc_id == root ) buff = obs%sta(:,m)
   call MPI_BCAST(buff,      obs%num,    MPI_REAL, root, comm, ierr)
   do n = 1, obs%num
     obs%sta(n,m) = buff(n)
   enddo
enddo
 
do m = 1, 2
   if ( my_proc_id == root ) buff = obs%roi(:,m)
   call MPI_BCAST(buff,      obs%num,    MPI_REAL, root, comm, ierr)
   do n = 1, obs%num
     obs%roi(n,m) = buff(n)
   enddo
enddo

call MPI_BARRIER(comm, ierr)

  return

end subroutine bcast_obs
!=======================================================================================
