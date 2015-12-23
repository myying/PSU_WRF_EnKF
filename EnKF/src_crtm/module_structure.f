module constants

!-----------------------------------------------------------------------------
!  PURPOSE: Common reference point for constants.
!
!  METHOD:  Straightforward definitions.
!------------------------------------------------------------------------------

   IMPLICIT NONE

!------------------------------------------------------------------------------
!  [1.0] Physical parameter constants (all NIST standard values):
!------------------------------------------------------------------------------
   real, parameter :: g=9.81,                       &! gravitational constant
                      to=300.,                      &! base pot. temp.
                      cp=1004.5,                    &! specific heat
                      rd=287.,                      &! dry air gas constant
                      rv=461.6,                     &! moist air gas constant
                      lhv=2.50*10**6,               &! latent heat of vap.
                      R_earth=6374.,                &! radius of the earth
                      Pr=100000.,                   &! reference pressure
                      Po=101325.,                   &! altimeter reference P
                      Talt=288.15,                  &! attimeter reference T
                      pbrh=611.,                    &! vapor pressure at 0 C (Pa)
                      T_frez=273.15,                &! water freezing point
                      RvRd=Rv/Rd,                   &! Rd/Rv
                      kappa=2./7.,                  &! kappa for pot. temp
                      pid=3.1415/180.,              &! radians to degrees
                      om_earth=(2.*3.1415)/(24.*3600.),&! earth rot
                      lap_std_atmos=0.0065,         &! standard atmosphere
                      th_denom=0.03727594,          &! Pr**(-kappa)
                      missing=-888888.,             &! missing value
                      es_alpha = 611.2,             &! saturation vapor pressure
                      es_beta = 17.67,              &!
                      es_gamma = 243.5
  real, parameter ::  pi = 3.1415927,               &!
                      deg_per_rad = 180./pi,        &!
                      rad_per_deg = pi / 180.,      &!
                      earth_radius_m = 6371200.
  integer,parameter :: PROJ_LATLON = 0,             &!
                       PROJ_MERC   = 1,             &!
                       PROJ_LC     = 3,             &!
                       PROJ_PS     = 5

end module constants

!==============================================================================

module namelist_define

   implicit none

!------------------------------------------------------------------------------
!  [1.0] Namelist parameters:
!------------------------------------------------------------------------------
!-- enkf_parameter
   integer       :: numbers_en          ! ensemble size
   character(len=10) :: expername       ! experiment name, such as hurricane, clear-air
   character(len=10), dimension(30) :: enkfvar   ! include the variables which will be updated and be used to calculate the XB
   character(len=10), dimension(20) :: updatevar ! updated variables
   integer       :: update_is, update_ie, update_js, update_je, update_ks, update_ke  ! domain to be updated
   real          :: inflate             !
   real          :: mixing              !
   integer       :: print_detail        ! used to debug the code
!-- parallel
   logical       :: manual_parallel, random_order
   integer       :: nmcpu, nicpu, njcpu

!-- use_osse
   logical       :: use_ideal_obs       ! .true. : use ideal observation, no real obs. position. Just for use_sounding or use_groundbase_radar
   integer       :: gridobs_is, gridobs_ie, gridobs_js, gridobs_je, gridobs_ks, gridobs_ke  ! domain to pick up data as gridobservation, for use_ideal_obs
   integer       :: gridobs_int_x, gridobs_int_k     ! horizontal and vertical interval grid for picking up gridobs, for use_ideal_obs
   logical       :: use_simulated       ! .true. : use simulated radar data

!-- use_hurricane_position_intensity
   logical       :: use_hurricane_PI    ! .true. : assimilated hurricane position and intensity 
   integer       :: hroi_hurricane_PI   ! horizontal radius of influence for hurricane PI
   integer       :: vroi_hurricane_PI   ! vertical radius of influence for hurricane PI

!-- use_surface_obs
   logical       :: use_surface         ! .true. : assimilated SURFACE LAND (SYNOPTIC, METAR) REPORTS
   integer       :: datathin_surface    ! 0=all data, 2=1/2 data, 10=1/10 data
                                        ! 2: get the 1st, 3rd, 5th ... data
                                        !-2: get the 2nd, 4th, 6th ... data
   integer       :: hroi_surface        ! horizontal radius of influence for surface
   integer       :: vroi_surface        ! vertical radius of influence for surface

!-- use_sounding_obs
   logical       :: use_sounding        ! .true. : assimilated UPPER-AIR (RAOB, PIBAL, RECCO, DROPS) REPORTS
   integer       :: datathin_sounding   ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_sounding       ! horizontal radius of influence for sounding 
   integer       :: vroi_sounding       ! vertical radius of influence for sounding 

!-- use_profiler_obs
   logical       :: use_profiler        ! .true. : assimilated WIND PROFILER REPORTS
   integer       :: datathin_profiler   ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_profiler       ! horizontal radius of influence for profiler 
   integer       :: vroi_profiler       ! vertical radius of influence for profiler 

!-- use_aircft_obs
   logical       :: use_aircft          ! .true. : assimilated AIREP/PIREP, AMDAR (ASDAR/ACARS), E-ADAS (AMDAR BUFR) AIRCRAFT REPORTS
   integer       :: datathin_aircft     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_aircft         ! horizontal radius of influence for aircft   
   integer       :: vroi_aircft         ! vertical radius of influence for aircft   

!-- use_metar_obs
   logical       :: use_metar           ! .true. : assimilated auto-meteo-station report
   integer       :: datathin_metar      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_metar          ! horizontal radius of influence for sfcshp   
   integer       :: vroi_metar          ! vertical radius of influence for sfcshp   

!-- use_sfcshp_obs
   logical       :: use_sfcshp          ! .true. : assimilated SURFACE MARINE (SHIP, BUOY, C-MAN PLATFORM) REPORTS
   integer       :: datathin_sfcshp     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_sfcshp         ! horizontal radius of influence for sfcshp   
   integer       :: vroi_sfcshp         ! vertical radius of influence for sfcshp   

!-- use_spssmi_obs
   logical       :: use_spssmi          ! .true. : assimilated DMSP SSM/I RETRIEVAL PRODUCTS (REPROCESSED WIND SPEED, TPW)
   integer       :: datathin_spssmi     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_spssmi         ! horizontal radius of influence for spssmi   
   integer       :: vroi_spssmi         ! vertical radius of influence for spssmi   

!-- use_atovs_obs
   logical       :: use_atovs           ! .true. : assimilated WIND PROFILER REPORTS
   integer       :: datathin_atovs      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_atovs          ! horizontal radius of influence for atovs    
   integer       :: vroi_atovs          ! vertical radius of influence for atovs    

!-- use_satwnd_obs
   logical       :: use_satwnd          ! .true. : assimilated SATELLITE-DERIVED WIND REPORTS
   integer       :: datathin_satwnd     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_satwnd         ! horizontal radius of influence for satwnd   
   integer       :: vroi_satwnd         ! vertical radius of influence for satwnd   

!-- use_gpspw_obs
   logical       :: use_gpspw          ! .true. : assimilated SATELLITE-DERIVED WIND REPORTS
   integer       :: datathin_gpspw     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_gpspw         ! horizontal radius of influence for satwnd   
   integer       :: vroi_gpspw         ! vertical radius of influence for satwnd   

!-- use_groundbase_radar_obs
   integer       :: radar_number        ! radar number
   logical       :: use_radar_rf        ! .true. : assimilated radar_rf 
   logical       :: use_radar_rv        ! .true. : assimilated radar_rv 
   integer       :: datathin_radar      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_radar          ! horizontal radius of influence for radar    
   integer       :: vroi_radar          ! vertical radius of influence for radar    

!-- use_airborne_radar
   logical       :: use_airborne_rf     ! .true. : assimilated airborne_rf 
   logical       :: use_airborne_rv     ! .true. : assimilated airborne_rv 
   integer       :: datathin_airborne   ! 0=all data, 2=1/2 data, 10=1/10 airborne
   integer       :: hroi_airborne       ! horizontal radius of influence for airborne    
   integer       :: vroi_airborne       ! vertical radius of influence for airborne    

!-- use_radiance
   logical       :: use_radiance        ! .true. : assimilated radiance 
   integer       :: datathin_radiance   ! 0=all data, 2=1/2 data, 10=1/10 radiance
   integer       :: hroi_radiance       ! horizontal radius of influence for radiance    
   integer       :: vroi_radiance       ! vertical radius of influence for radiance  



!-- Namelist contents :

   namelist /enkf_parameter / numbers_en, expername, enkfvar, updatevar,                          &
                                 update_is, update_ie, update_js, update_je, update_ks, update_ke,   &
                                 inflate, mixing, random_order, print_detail
   namelist /parallel       / manual_parallel, nmcpu, nicpu, njcpu
   namelist /osse           / use_ideal_obs, gridobs_is, gridobs_ie, gridobs_js, gridobs_je,      &
                                 gridobs_ks, gridobs_ke, gridobs_int_x, gridobs_int_k, use_simulated
   namelist /hurricane_PI   / use_hurricane_PI, hroi_hurricane_PI, vroi_hurricane_PI
   namelist /surface_obs    / use_surface, datathin_surface, hroi_surface, vroi_surface
   namelist /sounding_obs   / use_sounding, datathin_sounding, hroi_sounding, vroi_sounding
   namelist /profiler_obs   / use_profiler, datathin_profiler, hroi_profiler, vroi_profiler
   namelist /aircft_obs     / use_aircft, datathin_aircft, hroi_aircft, vroi_aircft
   namelist /metar_obs      / use_metar , datathin_metar , hroi_metar , vroi_metar 
   namelist /sfcshp_obs     / use_sfcshp, datathin_sfcshp, hroi_sfcshp, vroi_sfcshp
   namelist /spssmi_obs     / use_spssmi, datathin_spssmi, hroi_spssmi, vroi_spssmi
   namelist /atovs_obs      / use_atovs, datathin_atovs, hroi_atovs, vroi_atovs
   namelist /satwnd_obs     / use_satwnd, datathin_satwnd, hroi_satwnd, vroi_satwnd
   namelist /gpspw_obs      / use_gpspw, datathin_gpspw, hroi_gpspw, vroi_gpspw
   namelist /radar_obs      / radar_number, use_radar_rf, use_radar_rv, datathin_radar, hroi_radar, vroi_radar
   namelist /airborne_radar / use_airborne_rf, use_airborne_rv, datathin_airborne, hroi_airborne, vroi_airborne
   namelist /radiance / use_radiance, datathin_radiance, hroi_radiance, vroi_radiance


end module namelist_define

!==============================================================================

module mapinfo_define

 IMPLICIT NONE 

 ! Define data structures to define various projections

      TYPE proj_info

      INTEGER    :: code     ! Integer code for projection type
      REAL       :: lat1    ! SW latitude (1,1) in degrees (-90->90N)
      REAL       :: lon1    ! SW longitude (1,1) in degrees (-180->180E)
      REAL       :: dx       ! Grid spacing in meters at truelats, used
                             ! only for ps, lc, and merc projections
      REAL       :: dlat     ! Lat increment for lat/lon grids
      REAL       :: dlon     ! Lon increment for lat/lon grids
      REAL       :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      REAL       :: truelat1 ! First true latitude (all projections)
      REAL       :: truelat2 ! Second true lat (LC only)
      REAL       :: hemi     ! 1 for NH, -1 for SH
      REAL       :: cone     ! Cone factor for LC projections
      REAL       :: polei    ! Computed i-location of pole point
      REAL       :: polej    ! Computed j-location of pole point
      REAL       :: rsw      ! Computed radius to SW corner
      REAL       :: rebydx   ! Earth radius divided by dx
      LOGICAL    :: init     ! Flag to indicate if this struct is
                                 ! ready for use
      INTEGER    :: nx
      INTEGER    :: ny
      END TYPE proj_info
end module mapinfo_define

!==============================================================================
module obs_define

!                    _
!                    | radar_stn   % numObs, levels, lat, lon, elv, i_radar, j_radar 
! obs % radar(200) % | radar_point % ii, jj, hgt 
!                    | rf           
!                    | rv
!                    - 

   implicit none

   TYPE Radar_stn_type
        INTEGER                :: idn           ! station id number
        CHARACTER (LEN = 4)    :: idc           ! station id character
        CHARACTER (LEN = 50)   :: name          ! Station name
        CHARACTER (LEN = 19)   :: date_char     ! CCYY-MM-DD_HH:MM:SS date
        INTEGER                :: numObs        ! number of Obs
        INTEGER                :: levels        ! number of levels
        REAL                   :: lat           ! Latitude in degree
        REAL                   :: lon           ! Longitude in degree
        REAL                   :: elv           ! Elevation in m
        real, dimension(20)    :: elv_ang       ! elevation angles
        real                   :: i_radar       ! radar base position according to wrf model domain
        real                   :: j_radar       ! i_radar: west-east; j_radar: south-north
        real                   :: min_dis       ! distance between 1st data and radar
        real                   :: max_dis       ! valid data distance
        real                   :: d_dis         ! resolution in radial
        real                   :: d_azim        ! resolution in azimuth
        real                   :: err_rf        ! reflective observation error 
        real                   :: err_rv        ! radial velocity observation error
   END TYPE Radar_stn_type

   type radar_data_point_type
        real                    :: ii          ! each valid radar data point according to wrf domain
        real                    :: jj          ! ii: west-east; jj: south-north
        real                    :: kk          ! kk: half eta level
        real                    :: hgt         ! each valid radar data height above sea level
        real                    :: rdis        ! radial distance   ! just output for plot
        real                    :: azim        ! azimuth           ! just output for plot
        real                    :: elev        ! elevation angles  ! just output for plot
   end type radar_data_point_type

   type radar_data_type
        type ( Radar_stn_type )        :: radar_stn
        type ( radar_data_point_type ), allocatable,dimension(:) :: radar_point 
        real, allocatable,dimension(:)             :: rf    ! radar reflective 
        real, allocatable,dimension(:)             :: rv    ! radar radial velocity
   end type radar_data_type
 
   type airborne_data_type
        integer                                    :: num
        real, allocatable,dimension(:)             :: lat, lon, height, radar_ii, radar_jj   !for radar location
        real, allocatable,dimension(:)             :: azim, elev, range, rf, rv, ii, jj, hh
   end type airborne_data_type 

   type gts_data_type
        integer                                    :: num
        character(len=12),allocatable,dimension(:) :: platform
        character(len=19),allocatable,dimension(:) :: date
        real, allocatable,dimension(:)             :: latitude, longitude, elevation
        real, allocatable,dimension(:,:)           :: slp                              !slp(:,3), 3: data, qc, error
        real, allocatable,dimension(:,:)           :: pw                               !pw(:,3), 3: data, qc, error
        integer, allocatable,dimension(:)          :: levels
        real, allocatable,dimension(:,:,:)         :: pres, spd, wd, height, t, td, rh !(:,:,3), 3: data, qc, error
   end type gts_data_type

   type radiance_data_type
        integer                                    :: num
        character(len=12),allocatable,dimension(:) :: platform
        real, allocatable,dimension(:)             :: lat, lon,ii, jj, tb, err
        integer, allocatable,dimension(:)          :: ch, hroi, hroi_d
   end type radiance_data_type


   type raw_type
        integer                                  :: radar_stn_num
        type ( Radar_data_type ), allocatable,dimension( : )  :: radar
        type ( airborne_data_type  )             :: airborne
        type ( gts_data_type      )              :: gts
        type ( radiance_data_type      )         :: radiance
   end type raw_type

   type obs_type
       integer                                   :: num          !! observation number
       real, allocatable, dimension(:)           :: dat              !! observation data
       character(len=10), allocatable, dimension(:) :: type  !! observation type
       real, allocatable, dimension(:)           :: err         !! observation error
       real, allocatable, dimension(:,:)         :: position    !! (ob_num,4)
                                                                    !! observation data position in wrf domain (ii,jj,kk,hh)
                                                                    !! for radar, it is km; for sounding, it's mb;
       real, allocatable, dimension(:,:)         :: sta         !! (ob_num,4)
                                                                !! station attribution. For radar, it includes
                                                                !! radar station's (ii,jj,kk,hh) in wrf domain;
                                                                !! for surface obs, it include (elevation, station pressure, t, q)
       integer, allocatable, dimension(:,:)      :: roi         !! (ob_num,3) : 1=horizontal, 2=vertical, 3=h.(non-Q in radiance)
       character(len=12), allocatable, dimension(:) :: sat      !! Name of the Satellite
       integer, allocatable, dimension(:)        :: ch          !! channel of the satellite


   end type obs_type

   type ( raw_type )   :: raw
   type ( obs_type )   :: obs

!---------------------------------------------------------------------------------------------------------------------

end module obs_define
!==============================================================================
