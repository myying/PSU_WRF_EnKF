!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   module wrf_tools - module containing several subroutines and 
!                      functions to calculate commen meteorological 
!                      fields.
!
!     created Oct. 2004 Ryan Torn, U. Washington
!     add some new subroutines Zhiyong Meng 07/2005 TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE WRF_TOOLS

      use constants
      use mapinfo_define
      use map_utils
      use netcdf

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   altimeter - function that calculates the altimeter value at a 
!               grid point
!
!      pres - surface pressure value
!      hght - surface height value
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function altimeter(pres, hght)

        real, parameter  :: pow = (Rd*lap_std_atmos)/g,               &
                           invpow = g/(Rd*lap_std_atmos),               &
                           Patoinhg = 0.0002953
        real, intent(in) :: pres, hght
        real             :: altimeter

        altimeter = (pres * (1 + ((Po/pres)**pow) *                     &
                     (lap_std_atmos*hght/Talt)) ** invpow)               

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   colunm_int - function that calculates the column integral of a 
!                quantity such as water vapor.
!
!    infld - input field to take the column integral of
!     dens - dry air density at each vertical level
!     hght - staggered geopotential height values
!       iz - number of grid points in z direction
!
!     created June 2005 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function column_int(infld, dens, hght, iz)

        integer, intent(in) :: iz
        real, intent(in)    :: infld(iz), dens(iz), hght(iz+1)

        integer             :: k
        real                :: column_int

        column_int = 0.
        do k = 1, iz
          column_int = column_int + infld(k) * dens(k) *           &
                            ( hght(k+1) - hght(k) )
        enddo

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   couple - subroutine that couples quantities needed in wrf boundary
!            files since variables are coupled in the model
!
!   varble - data input to be coupled, then output
!     name - name of variable to be coupled (one letter)
!     mass - mass field
!      msf - map scale factor (needed for u and v fields)
!       ix - x dimension size of variable
!       iy - y dimension size of variable
!      ixm - x dimension size of mass field
!      iym - y dimension size of mass field
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine couple(varble, vname, mass, msf, ix, iy, ixm, iym)

      character (len=*), intent(in) :: vname
      integer, intent(in)           :: ix, iy, ixm, iym
      real, intent(in)              :: mass(ixm,iym), msf(ix,iy)
      real, intent(inout)           :: varble(ix,iy)

      integer :: i, j
      real :: mauv(ix,iy)

      select case (vname)
        case ('U') ! couple the u field
    
          ! First, put mass on u points
          do i = 2, ixm; do j = 1, iym
            mauv(i,j) = 0.5 * ( mass(i,j) + mass(i-1,j) )
          enddo ; enddo

          do j = 1, iym
            mauv(1,j) = mass(1,j)
            mauv(ixm+1,j) = mass(ixm,j)
          enddo

          ! couple to u
          do i = 1, ix; do j = 1, iy
            varble(i,j) = (varble(i,j) * mauv(i,j)) / msf(i,j)
          enddo ; enddo

        case ('V') ! couple the v field

          do i = 1, ixm ; do j = 2, iym
            mauv(i,j) = 0.5 * ( mass(i,j) + mass(i,j-1) )
          enddo ; enddo

          do i = 1, ixm
            mauv(i,1) = mass(i,1)
            mauv(i,iym+1) = mass(i,iym)
          enddo

          do i = 1, ix; do j = 1, iy
            varble(i,j) = (varble(i,j) * mauv(i,j)) / msf(i,j)
          enddo ; enddo

        case default  ! couple all other fields

          do i = 1, ix; do j = 1, iy
            varble(i,j) = varble(i,j) * mass(i,j)
          enddo ; enddo

      end select

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   destag_zstag - subroutine that destaggers staggered quantities in 
!                  the vertical direction
!
!     zus - unstaggered eta values
!      zs - staggered eta values
!      iz - number of unstaggered vertical grid points
!     vin - input staggered field
!    vout - output unstaggered field
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine destag_zstag(zus, zs, iz, vin, vout)

        real, parameter     :: shght = (Rd*256.)/g
        integer, intent(in) :: iz
        real, intent(in)    :: zus(iz), zs(iz+1), vin(iz+1)
        real, intent(out)   :: vout(iz)

        integer :: k
        real :: znfac, vinl(iz+1)

        do k = 1, (iz+1)
          vinl(k) = exp( - vin(k) / shght )
        enddo

        do k = 1, iz
          znfac = (zus(k) - zs(k)) / (zs(k+1)-zs(k))
          vout(k) = znfac * vinl(k+1) + (1.-znfac) * vinl(k)
          vout(k) = - shght * log(vout(k))
        enddo

        return
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   dir_spd2xy - subroutine that transforms a vector speed and 
!                direction to components in the x and y direction
!
!     dir - vector direction ( 0 is true north )
!     spd - vector magnitude
!    xdir - component in the x direction
!    ydir - component in the y direction
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine dir_spd2xy(dir, spd, xdir, ydir)

        real, intent(in)  :: dir, spd
        real, intent(out) :: xdir, ydir

        xdir = spd * cos(pid*(90.+dir))
        ydir = -spd * sin(pid*(90.+dir))

        return
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   dry_density - function that gives the dry air density at a mass 
!                 grid point.
!
!     eta - vertical eta levels bounding a grid box
!    hght - geopotential bounding a grid box
!    mass - dry-air mass at a grid point
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function dry_density(eta, hght, mass)

        real, intent(in) :: eta(2), hght(2), mass

        real :: dry_density

        dry_density = -mass * ((eta(2) - eta(1)) / (hght(2) - hght(1)))

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   earth_dist - function that calculates the earth distance between 
!                two points.
!
!    lon1 - longitude of point 1
!    lat1 - latitude of point 1
!    lon2 - longitude of point 2
!    lat2 - latitude of point 2
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function earth_dist(lon1, lat1, lon2, lat2)

        real, intent(in) :: lon1, lat1, lon2, lat2
        real :: earth_dist, pro

        pro =    sin(lat1*pid) * sin(lat2*pid) +           &
                cos(lat1*pid) * cos(lat2*pid) *            &
                cos((lon2-lon1)*pid)

        if ( abs(pro) .gt. 1.001 ) then
          print*,'In earth_dist, pro > 1, returning',pro
          earth_dist = missing
        elseif ( abs(pro) .gt. 1 ) then
          earth_dist = 0.
        else
          earth_dist = R_earth * acos(pro)
        endif

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   eta_to_pres - subroutine that calculates the pressure of a 
!                 vertical column of grid points
!
!     eta - vertical eta levels
!    mass - dry-air mass for column
!   vapor - column of water vapor mixing ratio
!    hght - column of geopotential values
!    thta - column of potential temperatures
!      iz - number of vertical grid points
!    pres - output vertical pressure column
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine eta_to_pres(eta, mass, vapor, hght, thta, iz, pres)

        integer, intent(in) :: iz
        real, intent(in)    :: eta(iz+1), mass, vapor(iz), hght(iz+1),    &
                              thta(iz)
        real, intent(out)   :: pres(iz)

        integer :: i, is, ie
        real    :: rho_dry, kapdiv

        kapdiv = 1./(1.-kappa)
        do i = 1, iz
          is=i  ;  ie = is+1
          rho_dry = dry_density(eta(is:ie),hght(is:ie),mass)
          pres(i) = (rho_dry * Rd * th_denom * thta(i) *                  &
                        (1. + vapor(i) * RvRd )) ** kapdiv
        enddo

        return
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   far2kelvin - function that converts temperatures from degrees F to 
!                Kelvin
!
!     tin - input temperature in deg. F
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function far2kelvin(tin)

        real, intent(in)  :: tin

        real :: far2kelvin

        far2kelvin = (5./9.)*(tin-32) + 273.15

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   find_eta - function that finds the model level closest to the 
!              specified eta value.
!
!        etain - model eta levels
!     eta_want - desired eta level
!           iz - number of eta levels
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function find_eta(etain, eta_want, iz)

        integer, intent(in) :: iz
        real, intent(in)    :: etain(iz), eta_want

        integer :: find_eta, k
        real :: dist, edist

        find_eta = 0
        edist = 1.
        do k = 1, iz
          dist = abs(etain(k) - eta_want)
          if ( dist .lt. edist ) then
            find_eta = k
            edist = dist
          endif
        enddo

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   horiz_interp - function that horizontally interpolates a field
!
!    cdat - input 2x2 data to interpolate
!      di - weighting coefficient in x direction
!      dj - weighting coefficient in y direction
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function horiz_interp(cdat, di, dj)

        real, intent(in) :: cdat(2,2), di, dj
        real :: horiz_interp

        if ( (cdat(1,1).ne.missing) .AND. (cdat(2,1).ne.missing) .AND.       &
            (cdat(1,2).ne.missing) .AND. (cdat(2,2).ne.missing) ) then

          horiz_interp = (1.-di)*((1.-dj)*cdat(1,1) + dj*cdat(1,2)) +       &
                           (di)*((1.-dj)*cdat(2,1) + dj*cdat(2,2))

        else
         
          horiz_interp = missing

        endif

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_hght - function that interpolates a variable in the
!                 vertical based on height.
!
!    varble - input variable to interpolate
!      hght - vector of geopotential heights
!     level - vertical level desired
!        iz - number of vertical grid points
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function interp_hght(varble, hght, level, iz)

      integer, intent(in) :: iz
      real, intent(in)    :: varble(iz), hght(iz), level

      integer :: k, klev
      real :: m, interp_hght

      if ((hght(1) .le. level) .AND. (hght(iz) .ge. level)) then

        ! search for appropriate level
        do k = 2, iz
          if ( hght(k) > level ) then
            klev = k - 1
            exit
          endif
        enddo

        ! linearly interpolate
        m = (varble(klev+1) - varble(klev)) / (hght(klev+1)-hght(klev))
        interp_hght = m * (level - hght(klev)) + varble(klev)

      else  ! level outside domain

        interp_hght = missing

      endif

      return
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_pres - function that interpolates a field in pressure
!
!    varble - vertical column of variable to interpolate
!      pres - vertical column of pressure
!  pres_lev - desired pressure level to interpolate to
!        iz - number of vertical grid points
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function interp_pres(varble, pres, pres_lev, iz)

        integer, intent(in) :: iz
        real, intent(in)    :: varble(iz), pres(iz), pres_lev

        integer :: k, klev
        real    :: interp_pres, m

        if ((pres(1) .ge. pres_lev) .AND. (pres(iz) .le. pres_lev)) then

          do k = 2, iz
            if (pres(k) < pres_lev) then
              klev = k-1
              exit
            endif
          enddo

          m = ( varble(klev+1) - varble(klev) ) /          &
              ( log(pres(klev+1)) - log(pres(klev)) )
!               ( pres(klev+1) - pres(klev) )

          interp_pres = m * (log(pres_lev) - log(pres(klev))) +          &
!          interp_pres = m * (pres_lev - pres(klev)) +
                          varble(klev)

        else
          
          interp_pres = missing

        endif

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_thta_pres - function that vertically interpolates the 
!                      potential temperature to a desired pressure 
!                      level.  Special care is taken to account for 
!                      tropopause
!
!      varble - vertical column of theta to interpolate
!        pres - vertical column of pressure
!    pres_lev - desired pressure level to interpolate to
!          iz - number of vertical grid points
!
!     created May 2005 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function interp_thta_pres(varble, pres, pres_lev, iz)

        integer, intent(in) :: iz
        real, intent(in)    :: varble(iz), pres(iz), pres_lev

        integer :: k, ka, kb
        real    :: interp_thta_pres, m, mb, ma, bb, ba, pint, lap, bint

        if ((pres(1) .ge. pres_lev) .AND. (pres(iz) .le. pres_lev)) then

          if ( pres(2) .lt. pres_lev ) then ! level just above surface

            m = (varble(2) - varble(1)) / (pres(2) - pres(1))
            interp_thta_pres = m * (pres_lev - pres(1)) + varble(1)

          elseif (pres(iz-1).gt.pres_lev) then ! level just below top

            m = (varble(iz) - varble(iz-1)) / (pres(iz) - pres(iz-1))
            interp_thta_pres = m * (pres_lev-pres(iz-1)) + varble(iz-1)

          else ! anywhere in between

            ! find the levels where desired pressure level exists      
            do k = 3, (iz-1)
              if (pres(k) < pres_lev) then
                kb = k-1
                ka = k
                exit
              endif
            enddo

            ! calculate slope and intercept above and below point
            mb = (varble(kb) - varble(kb-1)) / (pres(kb) - pres(kb-1))
            bb = varble(kb-1) - mb * pres(kb-1)
            ma = (varble(ka+1) - varble(ka)) / (pres(ka+1) - pres(ka))
            ba = varble(ka) - ma * pres(ka)
            pint = (ba-bb) / (mb-ma)

            ! interpolate linearly in pressure
            if ( pres_lev .ge. pint ) then
              lap = mb ; bint = bb
            else
              lap = ma ; bint = ba 
            endif
            interp_thta_pres = lap * pres_lev + bint
 
          endif

        else
          
          interp_thta_pres = missing

        endif

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ll_in_domain - function that decides whether a lat-lon pair is 
!                  within a WRF domain
!
!     lat - latitude of point
!     lon - longitude of point
!    proj - projection structure of the domain
!      ix - number of grid points in x
!      iy - number of grid points in y
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function ll_in_domain(lat, lon, proj, ix, iy)

        integer, intent(in)         :: ix, iy
        real, intent(in)            :: lat, lon
        type(proj_info), intent(in) :: proj

        integer :: ii, jj
        logical :: ll_in_domain
        real :: ir, jr

        call latlon_to_ij(proj, lat, lon, ir, jr)
        ii = floor(ir)  ;  jj = floor(jr)

        if ( (ii .ge. 1) .AND. (ii .le. (ix-1)) .AND.        &
             (jj .ge. 1) .AND. (jj .le. (iy-1)) ) then

          ll_in_domain = .true.

        else
         
          ll_in_domain = .false.

        endif

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   mixrat_to_tdew - function that converts from a mixing ratio to a 
!                    dew point temperature.
!
!    qvap - input water vapor mixing ratio
!    pres - input pressure in Pa
!
!     created April 2005 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function mixrat_to_tdew(qvap, pres)

        real, intent(in) :: qvap, pres

        real :: evap, mixrat_to_tdew

        evap = qvap * RvRd * pres
        mixrat_to_tdew = 1./((1./T_frez) - (Rv/lhv)*log(evap/pbrh))

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   mp_use - subroutine that reads the microphysics method and sets 
!            logicals for each microphysics type. Valid for WRF version 
!            2.0.2
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine mp_use(imp, qcu, qru, qiu, qsu, qgu)

        integer, intent(in)  :: imp
        logical, intent(out) :: qcu, qru, qiu, qsu, qgu

        select case(imp)
          case (0)  ! No Microphysics
            qcu = .FALSE.
            qru = .FALSE.
            qsu = .FALSE.
            qiu = .FALSE.
            qgu = .FALSE.
          case (1)  ! Kessler warm rain
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .FALSE.
            qiu = .FALSE.
            qgu = .FALSE.
          case (2)  ! Lin et. al. scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .TRUE.
            qiu = .TRUE.
            qgu = .TRUE.
          case (3)  ! WSM 3 class simple ice scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .FALSE.
            qiu = .FALSE.
            qgu = .FALSE.
          case (4)  ! NCEP 5 class scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .TRUE.
            qiu = .TRUE.
            qgu = .FALSE.
          case (5)  ! Eta Ferrier 2-class scheme
            qcu = .TRUE.
            qru = .FALSE.
            qsu = .FALSE.
            qiu = .FALSE.
            qgu = .FALSE.
          case (6)  ! WSM 6 class scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .TRUE.
            qiu = .TRUE.
            qgu = .TRUE.
          case (98)  ! NCEP 3 class simple ice scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .FALSE.
            qiu = .FALSE.
            qgu = .FALSE.
          case (99)  ! NCEP 5 class scheme
            qcu = .TRUE.
            qru = .TRUE.
            qsu = .TRUE.
            qiu = .TRUE.
            qgu = .FALSE.
        end select

        return
        end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   pres_alt_to_pres - function that converts from pressure altitude 
!                      to pressure.  (used in decoding ACARS)
!
!   hght - input pressure altitude
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function pres_alt_to_pres(hght)

        real, parameter  :: C1 = g / (lap_std_atmos * Rd)
        real, intent(in) :: hght

        real :: pres_alt_to_pres

        pres_alt_to_pres = Po * exp( C1 * log( 1 -         &
                             (lap_std_atmos * hght) / Talt) )

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   radar_velocity - function that calculates the velocity parallel 
!                    to the radar beam
!
!    uin - input u velocity at a point
!    vin - input v velocity at a point
!   plat - latitude of velocity point
!   plon - longitude of velocity point
!   rlat - latitude of radar
!   rlon - longitude of radar
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function radar_velocity(uin, vin, plat, plon, rlat, rlon)

        real, intent(in) :: uin, vin, plat, plon, rlat, rlon

        real :: qrain, pos_vec(2), radar_velocity

        pos_vec(1) = plon - rlon ; pos_vec(2) = plat - rlat ;
        pos_vec(:) = pos_vec(:) / sqrt(sum(pos_vec(:)*pos_vec(:)))

        radar_velocity = uin * pos_vec(1) + vin * pos_vec(2)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   reflectivity_qi - function that calculates radar reflctivity for 
!                     a microphysics package that contains ice.  
!                     Coefficients taken from WRF RIP converter
!
!    rho_dry - dry-air density at a grid point
!      qrain - rain water mixing ratio at grid point
!       qice - ice mixing ratio at a grid point
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function reflectivity_qi(rho_dry, qrain, qice)

        real, intent(in) :: rho_dry, qrain, qice

        real :: zfac, reflectivity_qi

        zfac = 3.630803e-9 * 1.e18 * ((rho_dry * qrain)**1.75)
        zfac = zfac + 2.18500e-10 * 1.e18 * ((rho_dry * qice)**1.75)
        
        reflectivity_qi = max(10. * log(zfac) / log(10.), 0.)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   reflectivity_qg - function that calculates radar reflectivity for 
!                     a microphysics package that contains graupel.
!                     Coefficients taken from WRF RIP converter
!
!    rho_dry - dry-air density at a grid point
!      qrain - rain water mixing ratio at a grid point
!       qice - ice mixing ratio at a grid point
!     qgraup - graupel mixing ratio at a grid point
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function reflectivity_qg(rho_dry, qrain, qice, qgraup)

        real, intent(in) :: rho_dry, qrain, qice, qgraup

        real :: zfac, reflectivity_qg    

        zfac = 3.630803e-9 * 1.e18 * ((rho_dry * qrain)**1.75)
        zfac = zfac + 2.18500e-10 * 1.e18 * ((rho_dry * qice)**1.75)
        zfac = zfac + 1.033267e-9 * 1.e18 * ((rho_dry * qgraup)**1.75)

        reflectivity_qg = max(10. * log(zfac) / log(10.), 0.)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   reflectivity_qr - function that calculates radar reflectivity for 
!                     a microphysics package that contains rain.
!                     Coefficients taken from WRF RIP converter
!
!   rho_dry - dry-air density at a grid point
!      temp - temperature at a grid point
!     qrain - rain water mixing ratio at a grid point
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function reflectivity_qr(rho_dry, temp, qrain)

        real, intent(in) :: rho_dry, temp, qrain


        real :: zfac, reflectivity_qr

        if ( temp > T_frez ) then
          zfac = 3.630803e-9 * 1.e18 * ((rho_dry * qrain)**1.75)
        else 
          zfac = 2.18500e-10 * 1.e18 * ((rho_dry * qrain)**1.75)
        endif

        reflectivity_qr = max(10. * log(zfac) / log(10.), 0.)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   rel_humidity - function that calculates the relative humidity at 
!                  a grid point
!
!    qvap - water vapor mixing ratio
!    temp - temperature 
!    pres - pressure
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function rel_humidity(qvap, temp, pres)

        real, intent(in) :: qvap, temp, pres

        real :: es, ws, rel_humidity

        es = pbrh * exp((lhv/Rv)*(1/T_frez - 1/temp))
        ws = (es / (pres-es)) / RvRd
        rel_humidity = ( qvap / ws ) * 100.

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   rel_humidity_td - function that gives relative humidity given the 
!                     temperature, dew point and pressure
!
!    a_temp - air temperature
!   td_temp - dew point temperature
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function rel_humidity_td(t_air, t_dew)

        real, intent(in) :: t_air, t_dew

        real :: e, es, rel_humidity_td

        e  = pbrh * exp((lhv/Rv)*(1/T_frez - 1/t_dew))
        es = pbrh * exp((lhv/Rv)*(1/T_frez - 1/t_air))

        rel_humidity_td = (e / es)*100.

        return
        end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   rel_humidity_to_q - subroutine that convert rh to q and get q_error from rh_error and t_error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rel_humidity_to_q ( pres, temp, rh, q, t_error, rh_error, q_error )

   implicit none

   real, intent(in)          :: pres        ! Observed pressure.
   real, intent(in)          :: temp        ! Observed temperature.
   real, intent(in)          :: rh       ! Observed Relative humidity.
   real, intent(in)          :: t_error  ! Temperature observation error.
   real, intent(in)          :: rh_error ! RH observation error.
   real, intent(out)         :: q        ! Observed specific humidity.
   real, intent(out)         :: q_error  ! q observation error.

   real                      :: es, es0    ! Saturation vapor pressure.
   real                      :: qs, qs0    ! Saturation specific humidity.
   real                      :: pmb, rh1

   if ( abs(pres - missing) > 1. .and. abs(temp - missing) > 1. .and. &
        abs(rh - missing) > 1. ) then

        rh1 = rh
        if( rh > 100. ) rh1 = 100.

!.......Pa to hPa
        pmb = pres/100.             !pa --> hPa

!...... vapor pressure in mb
        es = 6.112*17.67*243.5*t_error/((temp-273.15+243.5)*(temp-273.15+243.5))*     &
             exp(17.67*(temp-273.15)/(temp-273.15+243.5))
        es0 = 6.112*exp(17.67*(temp-273.15)/(temp-273.15+243.5))

!...... saturation mixing ratio in kg/kg
        qs = 0.622*(pmb*es)/((pmb-es0)*(pmb-es0))
        qs0 = 0.622*es0/(pmb-es0)

!...... mixing ratio in g/kg
        q_error = 0.01*(rh1*qs + rh_error*qs0)*1000.
        q       = 0.01*rh1*qs0*1000.

!...... mixing ratio must be positive
        if ( q < 0. ) q = 0.

    else
        q = missing
        q_error = missing
    endif

    end subroutine rel_humidity_to_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  set_domain_proj - subroutine that gets a structure to use for 
!                    getting the i, j point later
!
!     file - file to use to get domain parameters
!    projt - structure that contains projection information
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine set_domain_proj(file, projt)

        character (len=*), intent(in) :: file
        type(proj_info), intent(out)  :: projt

        integer :: fid, rcode, mapp, ix, iy
        real :: out2d(1,1), la1, lo1, cen_lon, truelat1, truelat2,    &
                 grid_space
        
        call open_file(file, nf_nowrite, fid)
          rcode = nf_get_att_int(fid, nf_global,    &
                                  'WEST-EAST_GRID_DIMENSION', ix)
          rcode = nf_get_att_int(fid, nf_global,    &
                                  'SOUTH-NORTH_GRID_DIMENSION', iy)
          ix = ix - 1  ;  iy = iy - 1
          rcode = nf_get_att_real(fid, nf_global,'DX', grid_space)
          rcode = nf_get_att_real(fid, nf_global,'STAND_LON', cen_lon)
          rcode = nf_get_att_real(fid, nf_global,'TRUELAT1', truelat1)
          rcode = nf_get_att_real(fid, nf_global,'TRUELAT2', truelat2)
          rcode = nf_get_att_int(fid, nf_global, 'MAP_PROJ', mapp)
          call get_variable2d_local(fid, 'XLAT', 1, 1, 1, 1, 1, out2d)
          la1 = out2d(1,1)
          call get_variable2d_local(fid, 'XLONG', 1, 1, 1, 1, 1, out2d)
          lo1 = out2d(1,1)
        call close_file(fid)

        select case (mapp)
          case(1)
            call map_set(PROJ_LC, la1, lo1, grid_space, cen_lon,    &
                             truelat1, truelat2, ix, iy, projt)
          case(2)
            call map_set(PROJ_PS, la1, lo1, grid_space, cen_lon,   &
                             truelat1, truelat2, ix, iy, projt)
          case(3)
            call map_set(PROJ_MERC, la1, lo1, grid_space, cen_lon,   &
                             truelat1, truelat2, ix, iy, projt)
        end select

        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   slp_standard_atmos - function that calculates the sea level 
!                        pressure by assuming a US standard atmos for 
!                        temperature profile
!
!     temp - surface temperature
!     pres - surface pressure
!     qvap - surface water vapor mixing ratio
!     topo - surface elevation
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function slp_standard_atmos(temp, pres, qvap, topo)

        real, intent(in) :: temp, pres, qvap, topo

        real :: Tv, scale, slp_standard_atmos 

        Tv = temp * ( 1. + 0.608 * qvap )
        slp_standard_atmos = pres * exp((g * topo) / (Rd * Tv))

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  slp_ben - function that calculates the sea level pressure using the 
!            Benjamin tecnique
!
!     tmp - column of temperatures on mass levels
!       p - pressure on lowest mass level
!       q - column of water vapor mixing ratio
!       h - column of geopotential heights on mass levels
!      iz - number of grid points in z
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function slp_ben(tmp, p, q, h, iz)

        integer, intent(in) :: iz
        real, intent(in)    :: tmp(iz), p, q(iz), h(iz)

        integer :: i, j
        real :: tsurf, scale, t2km, q2km, topp2, slp_ben

        topp2 = h(1) + 2000.
        t2km = interp_hght(tmp, h, topp2, iz)
        q2km = interp_hght(q, h, topp2, iz)
        t2km = t2km*(1 + 0.608 * q2km)
        tsurf = t2km + lap_std_atmos * 2000. 
        scale = (g * h(1)) / (Rd*(2*tsurf + lap_std_atmos*h(1))/2)
        slp_ben = p * exp(scale)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   static_stability - function that calculates the static 
!                      stability given the potential temperature at
!                      two vertical levels and density.
!
!   
!     created May 2005 Ryan Torn U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function static_stability(thta1, thta2, pdiff, dens)

        real, intent(in) :: thta1, thta2, pdiff, dens
        
        real :: static_stability

        static_stability = (log(thta2) - log(thta1)) / pdiff / dens

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   temp_to_theta - function that converts from temperature to theta
!
!    temp - input temperature 
!    pres - input pressure 
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function temp_to_theta( temp, pres )

        real, intent(in) :: temp, pres

        real :: temp_to_theta

        temp_to_theta = temp / ( th_denom * (pres**kappa) )

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   theta_to_temp - function that converts theta to temperature
!
!    theta - input potential temperature
!     pres - input pressure
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function theta_to_temp( theta, pres )

        real, intent(in) :: theta, pres

        real :: theta_to_temp

        theta_to_temp = theta * th_denom * (pres**kappa)

        return
        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   uv_to_dirspd - subroutine that converts from a u and v component of 
!                  the wind to the direction and speed.
!
!    uin - input u component of the wind
!    vin - input v component of the wind
!    dir - wind direction
!    spd - wind speed
!
!     created April 2005 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine uv_to_dirspd(uin, vin, dir, spd)

        real, intent(in)  :: uin, vin
        real, intent(out) :: dir, spd

        spd = sqrt(uin*uin + vin*vin)
        dir = atan2(uin,vin) / pid + 180.

        return
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! toGrid -- This subroutine transfers obs. x to grid i and calculate its
!           distance to grid i and i+1
!           adaped from 3DVAR/WRF Zhiyong Meng 07/2005 TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine toGrid (x, ib, ie, i, dx, dxm)
        implicit none
        real,                     intent(in)  :: x
        integer,                  intent(in)  :: ib, ie
        real,                     intent(out) :: dx, dxm
        integer,                  intent(out) :: i

        i = int (x)

        if (i < ib) then
            i = ib
        else if (i >= ie) then
            i = ie - 1
        endif

        dx = x - real(i)
        dxm= 1.0 - dx

        end subroutine toGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  to_zk --- This subroutine transfers obs. p to grid k and calculate its
!           distance to grid k and k+1
!         adapted from 3DVAR/WRF Zhiyong Meng 07/2005 TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
       subroutine to_zk(obs_v, mdl_v, zk, iz )

       implicit none

       integer                                     :: k, iz
       real,                           intent(in)  :: obs_v
       real, dimension(1:iz),          intent(in)  :: mdl_v
       real,                           intent(out) :: zk

       if (obs_v > mdl_v(1) ) then
          zk = 1.
          return
       else if (obs_v < mdl_v(iz)) then
          zk = iz+iz
          return
       else 
         do k = 1,iz-1
            if(obs_v <= mdl_v(k) .and. obs_v >= mdl_v(k+1)) then
             zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
             exit
            endif
         enddo
       endif

       end subroutine to_zk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   zstag - subroutine that staggers destaggered quantities in
!                  the vertical direction
!
!     zus - unstaggered eta values
!      zs - staggered eta values
!      iz - number of unstaggered vertical grid points
!     vin - input staggered field
!    vout - output unstaggered field
!
!     adaped from Torn code Zhiyong Meng 07/2005 TAMU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine zstag(zus, zs, iz, vin, vout)

        real, parameter     :: shght = (Rd*256.)/g
        integer, intent(in) :: iz
        real, intent(in)    :: zus(iz), zs(iz+1), vin(iz)
        real, intent(out)   :: vout(iz+1)

        integer :: k
        real :: znfac, vinl(iz)

         vout(1) = 0.
        do k = 1, iz
          vinl(k) = exp( - vin(k) / shght )
        enddo

        do k = 2, iz
          znfac = (zs(k) - zus(k-1)) / (zus(k)-zus(k-1))
          vout(k) = znfac * vinl(k) + (1.-znfac) * vinl(k-1)
          vout(k) = - shght * log(vout(k))
!          print*,vout(k)
        enddo

!          znfac = (zs(iz+1) - zus(iz-1)) / (zus(iz)-zus(iz-1))
!          print*,znfac
!          vout(iz+1) = znfac * vinl(iz) + (1.-znfac) * vinl(iz-1)
           vout(iz+1) = 2 * vout(iz)- vout(iz-1)
!          vout(iz+1) = - shght * log(vout(iz+1))
!          print*,vout(iz+1)

        return
        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ENDMODULE WRF_TOOLS
