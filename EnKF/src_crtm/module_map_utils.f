!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

       MODULE map_utils

! Module that defines constants, data structures, and
! subroutines used to convert grid indices to lat/lon
! and vice versa.   
!
! SUPPORTED PROJECTIONS
! ---------------------
! Cylindrical Lat/Lon (code = PROJ_LATLON)
! Mercator (code = PROJ_MERC)
! Lambert Conformal (code = PROJ_LC)
! Polar Stereographic (code = PROJ_PS)
!
! REMARKS
! -------
! The routines contained within were adapted from routines
! obtained from the NCEP w3 library.  The original NCEP routines were less
! flexible (e.g., polar-stereo routines only supported truelat of 60N/60S)
! than what we needed, so modifications based on equations in Hoke, Hayes, and
! Renninger (AFGWC/TN/79-003) were added to improve the flexibility.  
! Additionally, coding was improved to F90 standards and the routines were
! combined into this module.  
!
! ASSUMPTIONS
! -----------
!  Grid Definition:
!    For mercator, lambert conformal, and polar-stereographic projections,
!    the routines within assume the following:
!
!       1.  Grid is dimensioned (i,j) where i is the East-West direction, 
!           positive toward the east, and j is the north-south direction, 
!           positive toward the north.  
!       2.  Origin is at (1,1) and is located at the southwest corner,
!           regardless of hemispere.
!       3.  Grid spacing (dx) is always positive.
!       4.  Values of true latitudes must be positive for NH domains
!           and negative for SH domains.
!
!     For the latlon projection, the grid origin may be at any of the
!     corners, and the deltalat and deltalon values can be signed to 
!     account for this using the following convention:
!       Origin Location        Deltalat Sign      Deltalon Sign
!       ---------------        -------------      -------------
!        SW Corner                  +                   +
!        NE Corner                  -                   -
!        NW Corner                  -                   +
!        SE Corner                  +                   -
!       
!  Data Definitions:
!       1. Any arguments that are a latitude value are expressed in 
!          degrees north with a valid range of -90 -> 90
!       2. Any arguments that are a longitude value are expressed in
!          degrees east with a valid range of -180 -> 180.
!       3. Distances are in meters and are always positive.
!       4. The standard longitude (stdlon) is defined as the longitude
!          line which is parallel to the y-axis (j-direction), along
!          which latitude increases (NOT the absolute value of latitude, but
!          the actual latitude, such that latitude increases continuously
!          from the south pole to the north pole) as j increases.  
!       5. One true latitude value is required for polar-stereographic and
!          mercator projections, and defines at which latitude the 
!          grid spacing is true.  For lambert conformal, two true latitude
!          values must be specified, but may be set equal to each other to
!          specify a tangent projection instead of a secant projection.  
!       
! USAGE
! -----
! To use the routines in this module, the calling routines must have the 
! following statement at the beginning of its declaration block:
!   USE map_utils
! 
! The use of the module not only provides access to the necessary routines,
! but also defines a structure of TYPE (proj_info) that can be used
! to declare a variable of the same type to hold your map projection
! information.  It also defines some integer parameters that contain
! the projection codes so one only has to use those variable names rather
! than remembering the acutal code when using them.  The basic steps are
! as follows:
!  
!   1.  Ensure the "USE map_utils" is in your declarations.
!   2.  Declare the projection information structure as type(proj_info):
!         TYPE(proj_info) :: proj
!   3.  Populate your structure by calling the map_set routine:
!         CALL map_set(code,lat1,lon1,dx,stdlon,truelat1,truelat2,nx,ny,proj)
!       where:
!         code (input) = one of PROJ_LATLON, PROJ_MERC, PROJ_LC, or PROJ_PS
!         lat1 (input) = Latitude of grid origin point (i,j)=(1,1) 
!                         (see assumptions!)
!         lon1 (input) = Longitude of grid origin 
!         dx (input) = grid spacing in meters (ignored for LATLON projections)
!         stdlon (input) = Standard longitude for PROJ_PS and PROJ_LC, 
!               deltalon (see assumptions) for PROJ_LATLON, 
!               ignored for PROJ_MERC
!         truelat1 (input) = 1st true latitude for PROJ_PS, PROJ_LC, and
!                PROJ_MERC, deltalat (see assumptions) for PROJ_LATLON
!         truelat2 (input) = 2nd true latitude for PROJ_LC, 
!                ignored for all others.
!         nx = number of points in east-west direction
!         ny = number of points in north-south direction
!         proj (output) = The structure of type (proj_info) that will be fully 
!                populated after this call
!
!   4.  Now that the proj structure is populated, you may call any 
!       of the following routines:
!       
!       latlon_to_ij(proj, lat, lon, i, j)
!       ij_to_latlon(proj, i, j, lat, lon)
!       truewind_to_gridwind(lon, proj, ugrid, vgrid, utrue, vtrue)
!       gridwind_to_truewind(lon, proj, utrue, vtrue, ugrid, vgrid)
!       compare_projections(proj1, proj2, same_proj)
!
!       It is incumbent upon the calling routine to determine whether or
!       not the values returned are within your domain bounds.  All values
!       of i, j, lat, and lon are REAL values.
!
!
! REFERENCES
! ----------
!  Hoke, Hayes, and Renninger, "Map Preojections and Grid Systems for
!       Meteorological Applications." AFGWC/TN-79/003(Rev), Air Weather
!       Service, 1985.
!
!  NCAR MM5v3 Modeling System, REGRIDDER program, module_first_guess_map.F
!  NCEP routines w3fb06, w3fb07, w3fb08, w3fb09, w3fb11, w3fb12
!
! HISTORY
! -------
! 27 Mar 2001 - Original Version
!               Brent L. Shaw, NOAA/FSL (CSU/CIRA)
! 02 Apr 2001 - Added routines to rotate winds from true to grid
!               and vice versa.
!               Brent L. Shaw, NOAA/FSL (CSU/CIRA)
! 09 Apr 2001 - Added compare_projections routine to compare two
!               sets of projection parameters.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use constants
      use mapinfo_define

      IMPLICIT NONE
   
      ! Define data structures to define various projections

      CONTAINS

      SUBROUTINE map_init(proj)
      ! Initializes the map projection structure to missing values

      IMPLICIT NONE
      TYPE(proj_info), INTENT(INOUT)  :: proj

      proj%lat1 =    -999.9
      proj%lon1 =    -999.9
      proj%dx    =    -999.9
      proj%stdlon =   -999.9
      proj%truelat1 = -999.9
      proj%truelat2 = -999.9
      proj%hemi     = 0.0
      proj%cone     = -999.9
      proj%polei    = -999.9
      proj%polej    = -999.9
      proj%rsw      = -999.9
      proj%init     = .FALSE.
      proj%nx       = -99
      proj%ny       = -99 
      END SUBROUTINE map_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE map_set(proj_code,lat1,lon1,dx,stdlon,
     &      truelat1,truelat2,idim,jdim,proj)
      ! Given a partially filled proj_info structure, this routine computes
      ! polei, polej, rsw, and cone (if LC projection) to complete the 
      ! structure.  This allows us to eliminate redundant calculations when
      ! calling the coordinate conversion routines multiple times for the
      ! same map.
      ! This will generally be the first routine called when a user wants
      ! to be able to use the coordinate conversion routines, and it
      ! will call the appropriate subroutines based on the 
      ! proj%code which indicates which projection type  this is.

      IMPLICIT NONE
    
      ! Declare arguments
      INTEGER, INTENT(IN)               :: proj_code
      REAL, INTENT(IN)                  :: lat1
      REAL, INTENT(IN)                  :: lon1
      REAL, INTENT(IN)                  :: dx
      REAL, INTENT(IN)                  :: stdlon
      REAL, INTENT(IN)                  :: truelat1
      REAL, INTENT(IN)                  :: truelat2
      INTEGER, INTENT(IN)               :: idim
      INTEGER, INTENT(IN)               :: jdim
      TYPE(proj_info), INTENT(OUT)      :: proj

      ! Local variables


      ! Executable code

      ! First, check for validity of mandatory variables in proj
      IF ( ABS(lat1) .GT. 90.001 ) THEN
        PRINT '(A)', 'Latitude of origin corner required as follows:'
        PRINT '(A)', '    -90N <= lat1 < = 90.N'
        STOP 'MAP_INIT'
      ENDIF
      IF ( ABS(lon1) .GT. 180.) THEN
        PRINT '(A)', 'Longitude of origin required as follows:'
        PRINT '(A)', '   -180E <= lon1 <= 180W'
        STOP 'MAP_INIT'
      ENDIF
      IF ((dx .LE. 0.).AND.(proj_code .NE. PROJ_LATLON)) THEN
        PRINT '(A)', 'Require grid spacing (dx) in meters be positive!'
        STOP 'MAP_INIT'
      ENDIF
      IF ((ABS(stdlon) .GT. 180.).AND.(proj_code .NE. PROJ_MERC)) THEN
        PRINT '(A)', 'Need orientation longitude (stdlon) as: '
        PRINT '(A)', '   -180E <= lon1 <= 180W' 
        STOP 'MAP_INIT'
      ENDIF
      IF (ABS(truelat1).GT.90.) THEN
        PRINT '(A)', 'Set true latitude 1 for all projections!'
        STOP 'MAP_INIT'
      ENDIF
   
      CALL map_init(proj) 
      proj%code  = proj_code
      proj%lat1 = lat1
      proj%lon1 = lon1
      proj%dx    = dx
      proj%stdlon = stdlon
      proj%truelat1 = truelat1
      proj%truelat2 = truelat2
      proj%nx = idim
      proj%ny = jdim
      IF (proj%code .NE. PROJ_LATLON) THEN
        proj%dx = dx
        IF (truelat1 .LT. 0.) THEN
          proj%hemi = -1.0 
        ELSE
          proj%hemi = 1.0
        ENDIF
        proj%rebydx = earth_radius_m / dx
      ENDIF
      pick_proj: SELECT CASE(proj%code)

      CASE(PROJ_PS)
        !PRINT '(A)', 'Setting up POLAR STEREOGRAPHIC map...'
        CALL set_ps(proj)

      CASE(PROJ_LC)
        !PRINT '(A)', 'Setting up LAMBERT CONFORMAL map...'
        IF (ABS(proj%truelat2) .GT. 90.) THEN
         PRINT '(A)', 'Second true latitude not set,assuming atangent'
         PRINT '(A,F10.3)', 'projection at truelat1: ', proj%truelat1
          proj%truelat2=proj%truelat1
        ELSE 
          ! Ensure truelat1 < truelat2
          proj%truelat1 = MIN(truelat1,truelat2)
          proj%truelat2 = MAX(truelat1,truelat2)
        ENDIF
        CALL set_lc(proj)
   
      CASE (PROJ_MERC)
        !PRINT '(A)', 'Setting up MERCATOR map...'
        CALL set_merc(proj)
   
      CASE (PROJ_LATLON)
        !PRINT '(A)', 'Setting up CYLINDRICAL EQUIDISTANT LATLON map...'
        ! Convert lon1 to 0->360 notation
        IF (proj%lon1 .LT. 0.) proj%lon1 = proj%lon1 + 360.
        proj%dlat = truelat1
        proj%dlon = stdlon 
   
      CASE DEFAULT
        PRINT '(A,I2)', 'Unknown projection code: ', proj%code
        STOP 'MAP_INIT'
    
      END SELECT pick_proj
      proj%init = .TRUE.
      RETURN
      END SUBROUTINE map_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE latlon_to_ij(proj, lat, lon, i, j)
      ! Converts input lat/lon values to the cartesian (i,j) value
      ! for the given projection. 

      IMPLICIT NONE
      TYPE(proj_info), INTENT(IN)          :: proj
      REAL, INTENT(IN)                     :: lat
      REAL, INTENT(IN)                     :: lon
      REAL, INTENT(OUT)                    :: i
      REAL, INTENT(OUT)                    :: j

      IF (.NOT.proj%init) THEN
        PRINT '(A)', 'You have not called map_set for this proj!'
        STOP 'LATLON_TO_IJ'
      ENDIF

      SELECT CASE(proj%code)
 
      CASE(PROJ_LATLON)
        CALL llij_latlon(lat,lon,proj,i,j)

      CASE(PROJ_MERC)
        CALL llij_merc(lat,lon,proj,i,j)

      CASE(PROJ_PS)
        CALL llij_ps(lat,lon,proj,i,j)
      
      CASE(PROJ_LC)
        CALL llij_lc(lat,lon,proj,i,j)

      CASE DEFAULT
        PRINT '(A,I2)', 'Unrecognized map projection code: '
     &                  , proj%code
        STOP 'LATLON_TO_IJ'
 
      END SELECT
      RETURN
      END SUBROUTINE latlon_to_ij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ij_to_latlon(proj, i, j, lat, lon)
      ! Computes geographical latitude and longitude for a given (i,j) point
      ! in a grid with a projection of proj

      IMPLICIT NONE
      TYPE(proj_info),INTENT(IN)          :: proj
      REAL, INTENT(IN)                    :: i
      REAL, INTENT(IN)                    :: j
      REAL, INTENT(OUT)                   :: lat
      REAL, INTENT(OUT)                   :: lon

      IF (.NOT.proj%init) THEN
       PRINT '(A)', 'You have not called map_set for this projection!'
       STOP 'IJ_TO_LATLON'
      ENDIF
      SELECT CASE (proj%code)

      CASE (PROJ_LATLON)
        CALL ijll_latlon(i, j, proj, lat, lon)

      CASE (PROJ_MERC)
        CALL ijll_merc(i, j, proj, lat, lon)

      CASE (PROJ_PS)
        CALL ijll_ps(i, j, proj, lat, lon)

      CASE (PROJ_LC)
        CALL ijll_lc(i, j, proj, lat, lon)

      CASE DEFAULT
        PRINT '(A,I2)', 'Unrecognized map projection code: ', 
     &               proj%code
        STOP 'IJ_TO_LATLON'

      END SELECT
      RETURN
      END SUBROUTINE ij_to_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE set_ps(proj)
      ! Initializes a polar-stereographic map projection from the partially
      ! filled proj structure. This routine computes the radius to the
      ! filled proj structure. This routine computes the radius to the
      ! southwest corner and computes the i/j location of the pole for use
      ! in llij_ps and ijll_ps.
      IMPLICIT NONE
 
      ! Declare args
      TYPE(proj_info), INTENT(INOUT)    :: proj

      ! Local vars
      REAL                              :: ala1
      REAL                              :: alo1
      REAL                              :: reflon
      REAL                              :: scale_top

      ! Executable code
      reflon = proj%stdlon + 90.
  
      ! Cone factor
      proj%cone = 1.0

      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Compute radius to lower-left (SW) corner
      ala1 = proj%lat1 * rad_per_deg
      proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.+proj%hemi*
     &                  SIN(ala1))

      ! Find the pole point
      alo1 = (proj%lon1 - reflon) * rad_per_deg
      proj%polei = 1. - proj%rsw * COS(alo1)
      proj%polej = 1. - proj%hemi * proj%rsw * SIN(alo1)
      !PRINT '(A,2F10.1)', 'Computed (I,J) of pole point: ',proj%polei,proj%polej
      RETURN
      END SUBROUTINE set_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE llij_ps(lat,lon,proj,i,j)
      ! Given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      IMPLICIT NONE

      ! Delcare input arguments
      REAL, INTENT(IN)               :: lat
      REAL, INTENT(IN)               :: lon
      TYPE(proj_info),INTENT(IN)     :: proj

      ! Declare output arguments     
      REAL, INTENT(OUT)              :: i !(x-index)
      REAL, INTENT(OUT)              :: j !(y-index)

      ! Declare local variables
    
      REAL                           :: reflon
      REAL                           :: scale_top
      REAL                           :: ala
      REAL                           :: alo
      REAL                           :: rm

      ! BEGIN CODE

  
      reflon = proj%stdlon + 90.
   
      ! Compute numerator term of map scale factor

      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Find radius to desired point
      ala = lat * rad_per_deg
      rm = proj%rebydx * COS(ala) * scale_top/(1. + proj%hemi *
     &                          SIN(ala))
      alo = (lon - reflon) * rad_per_deg
      i = proj%polei + rm * COS(alo)
      j = proj%polej + proj%hemi * rm * SIN(alo)
 
      RETURN
      END SUBROUTINE llij_ps
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ijll_ps(i, j, proj, lat, lon)

      ! This is the inverse subroutine of llij_ps.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  

      IMPLICIT NONE

      ! Declare input arguments
      REAL, INTENT(IN)                    :: i    ! Column
      REAL, INTENT(IN)                    :: j    ! Row
      TYPE (proj_info), INTENT(IN)        :: proj
    
      ! Declare output arguments
      REAL, INTENT(OUT)                   :: lat     ! -90 -> 90 North
      REAL, INTENT(OUT)                   :: lon     ! -180 -> 180 East

      ! Local variables
      REAL                                :: reflon
      REAL                                :: scale_top
      REAL                                :: xx,yy
      REAL                                :: gi2, r2
      REAL                                :: arccos

      ! Begin Code

      ! Compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj%stdlon + 90.
   
      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)

      ! Compute radius to point of interest
      xx = i - proj%polei
      yy = (j - proj%polej) * proj%hemi
      r2 = xx**2 + yy**2

      ! Now the magic code
      IF (r2 .EQ. 0.) THEN 
        lat = proj%hemi * 90.
        lon = reflon
      ELSE
        gi2 = (proj%rebydx * scale_top)**2.
        lat = deg_per_rad * proj%hemi * ASIN((gi2-r2)/(gi2+r2))
        arccos = ACOS(xx/SQRT(r2))
        IF (yy .GT. 0) THEN
          lon = reflon + deg_per_rad * arccos
        ELSE
          lon = reflon - deg_per_rad * arccos
        ENDIF
      ENDIF
  
      ! Convert to a -180 -> 180 East convention
      IF (lon .GT. 180.) lon = lon - 360.
      IF (lon .LT. -180.) lon = lon + 360.
      RETURN
  
      END SUBROUTINE ijll_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE set_lc(proj)
      ! Initialize the remaining items in the proj structure for a
      ! lambert conformal grid.

      IMPLICIT NONE
      
      TYPE(proj_info), INTENT(INOUT)     :: proj

      REAL                               :: arg
      REAL                               :: deltalon1
      REAL                               :: tl1r
      REAL                               :: ctl1r

      ! Compute cone factor
      CALL lc_cone(proj%truelat1, proj%truelat2, proj%cone)
      ! PRINT '(A,F8.6)', 'Computed cone factor: ', proj%cone
      ! Compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj%lon1 - proj%stdlon
      IF (deltalon1 .GT. +180.) deltalon1 = deltalon1 - 360.
      IF (deltalon1 .LT. -180.) deltalon1 = deltalon1 + 360.

      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)

      ! Compute the radius to our known lower-left (SW) corner
      proj%rsw = proj%rebydx * ctl1r/proj%cone * 
     &        (TAN((90.*proj%hemi-proj%lat1)*rad_per_deg/2.) / 
     &         TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))
     &              **proj%cone

      ! Find pole point
      arg = proj%cone*(deltalon1*rad_per_deg)
      proj%polei = 1. - proj%hemi * proj%rsw * SIN(arg)
      proj%polej = 1. + proj%rsw * COS(arg)  
      !PRINT '(A,2F10.3)', 'Computed pole i/j = ', proj%polei, proj%polej
      RETURN
      END SUBROUTINE set_lc                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE lc_cone(truelat1, truelat2, cone)

  ! Subroutine to compute the cone factor of a Lambert Conformal projection

      IMPLICIT NONE
      
      ! Input Args
      REAL, INTENT(IN)             :: truelat1  ! (-90 -> 90 degrees N)
      REAL, INTENT(IN)             :: truelat2  !   "   "  "   "     "

      ! Output Args
      REAL, INTENT(OUT)            :: cone

      ! Locals

      ! BEGIN CODE

      ! First, see if this is a secant or tangent projection.  For tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the 
      ! Earth surface at this latitude.  For secant projections, the cone
      ! intersects the Earth surface at each of the distinctly different
      ! latitudes
      IF (ABS(truelat1-truelat2) .GT. 0.1) THEN

        ! Compute cone factor following:
        cone=(ALOG(COS(truelat1*rad_per_deg))-ALOG(COS(truelat2*
     &        rad_per_deg))) / 
     &    (ALOG(TAN((90.-ABS(truelat1))*rad_per_deg*0.5 ))- 
     &     ALOG(TAN((90.-ABS(truelat2))*rad_per_deg*0.5 )) )
      ELSE
         cone = SIN(ABS(truelat1)*rad_per_deg )  
      ENDIF
      RETURN
      END SUBROUTINE lc_cone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ijll_lc( i, j, proj, lat, lon)

  ! Subroutine to convert from the (i,j) cartesian coordinate to the 
  ! geographical latitude and longitude for a Lambert Conformal projection.

  ! History:
  ! 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
  ! 
      IMPLICIT NONE

      ! Input Args
      REAL, INTENT(IN)              :: i    ! Cartesian X coordinate
      REAL, INTENT(IN)              :: j    ! Cartesian Y coordinate
      TYPE(proj_info),INTENT(IN)    :: proj ! Projection info structure

      ! Output Args                 
      REAL, INTENT(OUT)             :: lat  ! Latitude (-90->90 deg N)
      REAL, INTENT(OUT)             :: lon  ! Longitude (-180->180 E)

      ! Locals 
      REAL                          :: inew
      REAL                          :: jnew
      REAL                          :: r
      REAL                          :: chi,chi1,chi2
      REAL                          :: r2
      REAL                          :: xx
      REAL                          :: yy

      ! BEGIN CODE

      chi1 = (90. - proj%hemi*proj%truelat1)*rad_per_deg
      chi2 = (90. - proj%hemi*proj%truelat2)*rad_per_deg

      ! See if we are in the southern hemispere and flip the indices
      ! if we are. 
      IF (proj%hemi .EQ. -1.) THEN 
        inew = -i + 2.
        jnew = -j + 2.
      ELSE
        inew = i
        jnew = j
      ENDIF

      ! Compute radius**2 to i/j location
      xx = inew - proj%polei
      yy = proj%polej - jnew
      r2 = (xx*xx + yy*yy)
      r = SQRT(r2)/proj%rebydx
   
      ! Convert to lat/lon
      IF (r2 .EQ. 0.) THEN
        lat = proj%hemi * 90.
        lon = proj%stdlon
      ELSE
         
        ! Longitude
        lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/
     &           proj%cone
        lon = AMOD(lon+360., 360.)

        ! Latitude.  Latitude determined by solving an equation adapted 
        ! from:
        !  Maling, D.H., 1973: Coordinate Systems and Map Projections
        ! Equations #20 in Appendix I.  
          
        IF (chi1 .EQ. chi2) THEN
          chi = 2.0*ATAN( ( r/TAN(chi1) )**(1./proj%cone) * 
     &          TAN(chi1*0.5) )
        ELSE
          chi = 2.0*ATAN( (r*proj%cone/SIN(chi1))**(1./proj%cone) * 
     &             TAN(chi1*0.5)) 
        ENDIF
        lat = (90.0-chi*deg_per_rad)*proj%hemi

      ENDIF

      IF (lon .GT. +180.) lon = lon - 360.
      IF (lon .LT. -180.) lon = lon + 360.
      RETURN
      END SUBROUTINE ijll_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE llij_lc( lat, lon, proj, i, j)

  ! Subroutine to compute the geographical latitude and longitude values
  ! to the cartesian x/y on a Lambert Conformal projection.
      
      IMPLICIT NONE

      ! Input Args
      REAL, INTENT(IN)              :: lat   ! Latitude (-90->90 deg N)
      REAL, INTENT(IN)              :: lon   ! Longitude (-180->180 E)
      TYPE(proj_info),INTENT(IN)      :: proj! Projection info structure

      ! Output Args                 
      REAL, INTENT(OUT)             :: i     ! Cartesian X coordinate
      REAL, INTENT(OUT)             :: j     ! Cartesian Y coordinate

      ! Locals 
      REAL                          :: arg
      REAL                          :: deltalon
      REAL                          :: tl1r
      REAL                          :: rm
      REAL                          :: ctl1r
      

      ! BEGIN CODE
      
      ! Compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj%stdlon
      IF (deltalon .GT. +180.) deltalon = deltalon - 360.
      IF (deltalon .LT. -180.) deltalon = deltalon + 360.
      
      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)     
   
      ! Radius to desired point
      rm = proj%rebydx * ctl1r/proj%cone * 
     &      (TAN((90.*proj%hemi-lat)*rad_per_deg/2.) / 
     &       TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))**
     &         proj%cone

      arg = proj%cone*(deltalon*rad_per_deg)
      i = proj%polei + proj%hemi * rm * SIN(arg)
      j = proj%polej - rm * COS(arg)

      ! Finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the SW corner
      ! (what we assume) which is different than the original NCEP
      ! algorithms which used the NE corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      IF (proj%hemi .EQ. -1.) THEN
        i = 2. - i  
        j = 2. - j
      ENDIF
      RETURN
      END SUBROUTINE llij_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE set_merc(proj)
  
      ! Sets up the remaining basic elements for the mercator projection

      IMPLICIT NONE
      TYPE(proj_info), INTENT(INOUT)       :: proj
      REAL                                 :: clain


      !  Preliminary variables

      clain = COS(rad_per_deg*proj%truelat1)
      proj%dlon = proj%dx / (earth_radius_m * clain)

      ! Compute distance from equator to origin, and store in the 
      ! proj%rsw tag.

      proj%rsw = 0.
      IF (proj%lat1 .NE. 0.) THEN
        proj%rsw = (ALOG(TAN(0.5*((proj%lat1+90.)*rad_per_deg))))
     &          /proj%dlon
      ENDIF
      RETURN
      END SUBROUTINE set_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE llij_merc(lat, lon, proj, i, j)

      ! Compute i/j coordinate from lat lon for mercator projection
  
      IMPLICIT NONE
      REAL, INTENT(IN)              :: lat
      REAL, INTENT(IN)              :: lon
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL,INTENT(OUT)              :: i
      REAL,INTENT(OUT)              :: j
      REAL                          :: deltalon

      deltalon = lon - proj%lon1
      IF (deltalon .LT. -180.) deltalon = deltalon + 360.
      IF (deltalon .GT. 180.) deltalon = deltalon - 360.
      i = 1. + (deltalon/(proj%dlon*deg_per_rad))
      j = 1. + (ALOG(TAN(0.5*((lat + 90.) * rad_per_deg)))) /
     &        proj%dlon - proj%rsw

      RETURN
      END SUBROUTINE llij_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ijll_merc(i, j, proj, lat, lon)

      ! Compute the lat/lon from i/j for mercator projection

      IMPLICIT NONE
      REAL,INTENT(IN)               :: i
      REAL,INTENT(IN)               :: j    
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL, INTENT(OUT)             :: lat
      REAL, INTENT(OUT)             :: lon 


      lat = 2.0*ATAN(EXP(proj%dlon*(proj%rsw + j-1.)))*
     &               deg_per_rad - 90.
      lon = (i-1.)*proj%dlon*deg_per_rad + proj%lon1
      IF (lon.GT.180.) lon = lon - 360.
      IF (lon.LT.-180.) lon = lon + 360.
      RETURN
      END SUBROUTINE ijll_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE llij_latlon(lat, lon, proj, i, j)
 
      ! Compute the i/j location of a lat/lon on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: lat
      REAL, INTENT(IN)             :: lon
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: i
      REAL, INTENT(OUT)            :: j

      REAL                         :: deltalat
      REAL                         :: deltalon
      REAL                         :: lon360
      REAL                         :: latinc
      REAL                         :: loninc


      ! Compute deltalat and deltalon as the difference between the input 
      ! lat/lon and the origin lat/lon

      deltalat = lat - proj%lat1

      ! To account for issues around the dateline, convert the incoming
      ! longitudes to be 0->360.
      IF (lon .LT. 0) THEN 
        lon360 = lon + 360. 
      ELSE 
        lon360 = lon
      ENDIF    
      deltalon = lon360 - proj%lon1      
      IF (deltalon .LT. 0) deltalon = deltalon + 360.
 
      ! Compute i/j
      i = deltalon/proj%dlon + 1.
      j = deltalat/proj%dlat + 1.
      RETURN
      END SUBROUTINE llij_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ijll_latlon(i, j, proj, lat, lon)
 
      ! Compute the lat/lon location of an i/j on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: i
      REAL, INTENT(IN)             :: j
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: lat
      REAL, INTENT(OUT)            :: lon

      REAL                         :: deltalat
      REAL                         :: deltalon
      REAL                         :: lon360


      ! Compute deltalat and deltalon 

      deltalat = (j-1.)*proj%dlat
      deltalon = (i-1.)*proj%dlon
      lat = proj%lat1 + deltalat
      lon = proj%lon1 + deltalon

      IF ((ABS(lat) .GT. 90.).OR.(ABS(deltalon) .GT.360.)) THEN
        ! Off the earth for this grid
        lat = -999.
        lon = -999.
      ELSE
        lon = lon + 360.
        lon = AMOD(lon,360.)
        IF (lon .GT. 180.) lon = lon -360.
      ENDIF

      RETURN
      END SUBROUTINE ijll_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE gridwind_to_truewind(lon,proj,ugrid,vgrid,utrue,
     &        vtrue)
  
      ! Subroutine to convert a wind from grid north to true north.

      IMPLICIT NONE
      
      ! Arguments
      REAL, INTENT(IN)          :: lon   ! Longitude of point in degrees
      TYPE(proj_info),INTENT(IN):: proj  ! Projection info structure 
      REAL, INTENT(IN)          :: ugrid ! U-component, grid-relative
      REAL, INTENT(IN)          :: vgrid ! V-component, grid-relative
      REAL, INTENT(OUT)         :: utrue ! U-component, earth-relative
      REAL, INTENT(OUT)         :: vtrue ! V-component, earth-relative

      ! Locals
      REAL                      :: alpha
      REAL                      :: diff

      IF ((proj%code .EQ. PROJ_PS).OR.(proj%code .EQ. PROJ_LC))THEN
        diff = lon - proj%stdlon
        IF (diff .GT. 180.) diff = diff - 360.
        IF (diff .LT.-180.) diff = diff + 360.

        alpha = diff * proj%cone * rad_per_deg * 
     &         SIGN(1.,proj%truelat1)

        utrue = vgrid * SIN(alpha) + ugrid * COS(alpha)
        vtrue = vgrid * COS(alpha) - ugrid * SIN(alpha)
      ELSEIF ((proj%code .EQ. PROJ_MERC).OR.(proj%code .EQ. 
     &          PROJ_LATLON))THEN
        utrue = ugrid
        vtrue = vgrid
      ELSE  
        PRINT '(A)', 'Unrecognized projection.'
        STOP 'GRIDWIND_TO_TRUEWIND'
      ENDIF

      RETURN
      END SUBROUTINE gridwind_to_truewind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE truewind_to_gridwind(lon, proj, utrue, vtrue, 
     &     ugrid, vgrid)
 
      ! Subroutine to compute grid-relative u/v wind components from the earth-
      ! relative values for a given projection.

      IMPLICIT NONE
      
      ! Arguments
      REAL, INTENT(IN)                 :: lon
      TYPE(proj_info),INTENT(IN)       :: proj 
      REAL, INTENT(IN)                 :: utrue
      REAL, INTENT(IN)                 :: vtrue
      REAL, INTENT(OUT)                :: ugrid
      REAL, INTENT(OUT)                :: vgrid

      ! Locals
      REAL                             :: alpha
      REAL                             :: diff

      IF ((proj%code .EQ. PROJ_PS).OR.(proj%code .EQ. PROJ_LC))THEN
        
        diff = proj%stdlon - lon
        IF (diff .GT. 180.) diff = diff - 360.
        IF (diff .LT.-180.) diff = diff + 360.
        alpha = diff * proj%cone * rad_per_deg * SIGN(1.,proj
     &          %truelat1)
        ugrid = vtrue * SIN(alpha) + utrue * COS(alpha)
        vgrid = vtrue * COS(alpha) - utrue * SIN(alpha)

      ELSEIF((proj%code .EQ. PROJ_MERC).OR.(proj%code .EQ. 
     &             PROJ_LATLON)) THEN
        ugrid = utrue
        vgrid = vtrue
      ELSE
        PRINT '(A)', 'Unrecognized map projection.'
        STOP 'TRUEWIND_TO_GRIDWIND'
      ENDIF
      RETURN
      END SUBROUTINE truewind_to_gridwind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE compare_projections(proj1, proj2, same_proj)
  
      ! Subroutine to compare two proj_info structures to determine if the
      ! maps are the same.  

      IMPLICIT NONE

      ! Arguments
      TYPE(proj_info), INTENT(IN)   :: proj1
      TYPE(proj_info), INTENT(IN)   :: proj2
      LOGICAL, INTENT(OUT)          :: same_proj

      ! Locals
   
      same_proj = .false.

      ! Make sure both structures are initialized

      IF ((.NOT. proj1%init).OR.(.NOT.proj2%init)) THEN
        PRINT '(A)', 'COMPARE_PROJECTIONS: Map_Set not called yet!'
        same_proj = .false.
      ELSE  
        same_proj = .true.
      ENDIF
          
      ! Check projection type
      IF (same_proj) THEN
        IF (proj1%code .NE. proj2%code) THEN
          PRINT '(A)', 'COMPARE_PROJEC: Different projection type.'
          same_proj = .false.
        ELSE
          same_proj = .true.
        ENDIF
      ENDIF

      ! Check corner lat/lon
      IF (same_proj) THEN
        IF ( (ABS(proj1%lat1-proj2%lat1) .GT. 0.001) .OR. 
     &         (ABS(proj1%lon1-proj2%lon1) .GT. 0.001) ) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different lat1/lon1'
          same_proj = .false.
        ENDIF
      ENDIF

      ! Compare dx
      IF ((same_proj).AND.(proj1%code .NE. PROJ_LATLON)) THEN
        IF ( ABS(proj1%dx-proj2%dx).GT.1. ) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different dx'
          same_proj = .false.
        ENDIF
      ENDIF
      ! Compare dimensions
      IF ((same_proj).AND.( (proj1%nx .NE. proj2%nx).OR. 
     &                       (proj1%ny .NE. proj2%ny) ) ) THEN
        PRINT '(A)', 'COMPARE_PROJECTIONS: Different dimensions'
        same_proj = .false.
      ENDIF        

      IF ((same_proj).AND.(proj1%code .EQ. PROJ_LATLON)) THEN
        IF ( (proj1%dlat .NE. proj2%dlat).OR. 
     &        (proj1%dlon .NE. proj2%dlon)) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different dlat/dlon'
          same_proj = .false.
        ENDIF
      ENDIF

      ! Compare stdlon for Polar and LC projections
      IF ( (same_proj).AND. (  (proj1%code .EQ. PROJ_PS).OR. 
     &                          (proj1%code .EQ. PROJ_LC) ) ) THEN
        IF (proj1%stdlon .NE. proj2%stdlon) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different stdlon.'
          same_proj = .false.
        ENDIF
      ENDIF
      ! Compare true latitude1 if polar stereo, mercator, or lambert
      IF ( (same_proj).AND. ( (proj1%code .EQ. PROJ_PS) .OR. 
     &                         (proj1%code .EQ. PROJ_LC) .OR. 
     &                         (proj1%code .EQ. PROJ_MERC) ) ) THEN
        IF (proj1%truelat1 .NE. proj2%truelat1) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different truelat1'
          same_proj = .false.
        ENDIF
      ENDIF

      ! Compare true latitude2 if LC
      IF ( (same_proj).AND.(proj1%code .EQ. PROJ_LC)) THEN
        IF (proj1%truelat2 .NE. proj2%truelat2) THEN
          PRINT '(A)', 'COMPARE_PROJECTIONS: Different truelat2'
          same_proj = .false.
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE compare_projections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE compute_msf_lc(lat,truelat1,truelat2,msf)
  
      ! Computes the map scale factor for a Lambert Conformal 
      ! grid at a given
      ! latitude.

      IMPLICIT NONE
      REAL, INTENT(IN)            :: lat  ! latitude where msf is
      REAL, INTENT(IN)            :: truelat1
      REAL, INTENT(IN)            :: truelat2
      REAL, INTENT(OUT)           :: msf   

      REAL                        :: cone
      REAL                        :: psi1, psix, pole

      CALL lc_cone(truelat1,truelat2, cone)
      IF (truelat1 .GE. 0.) THEN
        psi1 = (90. - truelat1) * rad_per_deg
        pole =90.
      ELSE
        psi1 = (90. + truelat1) * rad_per_deg
        pole = -90.
      ENDIF
      psix = (pole - lat)*rad_per_deg
      msf = (SIN(psi1)/SIN(psix)) 
     &     * ((TAN(psix*0.5)/TAN(psi1*0.5))**cone)
      RETURN
      END SUBROUTINE compute_msf_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE compute_msf_ps(lat,truelat1,msf)

      ! Computes the map scale factor for a Polar Stereographic grid at
      ! latitude.

      IMPLICIT NONE
      REAL, INTENT(IN)            :: lat  ! latitude where msf is requested
      REAL, INTENT(IN)            :: truelat1
      REAL, INTENT(OUT)           :: msf

      REAL                        :: psi1, psix, pole

      IF (truelat1 .GE. 0.) THEN
        psi1 = (90. - truelat1) * rad_per_deg
        pole =90.
      ELSE
        psi1 = (90. + truelat1) * rad_per_deg
        pole = -90.
      ENDIF
      psix = (pole - lat)*rad_per_deg
      msf = ((1.+COS(psi1))/(1.0 + COS(psix)))
      RETURN
      END SUBROUTINE compute_msf_ps                                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      END MODULE map_utils
