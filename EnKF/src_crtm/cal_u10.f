!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the  10m wind, 2m temperature and moisture based on the
! similarity theory/
!
!  The unit for pressure   : psfc, ps, ps2     is Pa.
!  The unit for temperature: tg, ts, ts2, t2   is K.
!  The unit for moisture   : qs, qs2, q2       is kg/kg.
!  The unit for wind       : us, vs, u10, v10  is m/s.
!  The unit for height     : hs, roughness     is m.
!  xland and regime are dimensionless.
!
! Reference:
! ---------
!  Detail MM5/physics/pbl_sfc/mrfpbl/MRFPBL.F
!
!  Input Variables:
! 
!   psfc, tg               : surface pressure and ground temperature
!   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma leve1
!   ps2, ts2, qs2          : model variables at the second lowest half 
!                            sigma level
!
!
!  Constants:
!
!   hs                     : height at the lowest half sigma level
!   roughness              : roughness
!   xland                  : land-water-mask (=2 water, =1 land)
!
!  Output Variables:
!
!   regime                 : PBL regime
!   u10, v10               : 10-m high observed wind components
!   t2 , q2                : 2-m high observed temperature and mixing ratio
!
!-----------------------------------------------------------------------------!
!  
!                      psim  : mechanical psi at lowlest sigma leve1
!                      psim2 : mechanical psi at 2m 
!                      psimz : mechanical psi at 10m 
!
!-----------------------------------------------------------------------------!
!  adapted from 3dVAR/WRF July 2005  Zhiyong Meng TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SFC_WTQ ( psfc, tg, ps, ts, qs, us, vs, ps2, ts2, qs2, &
                           hs, roughness, xland, &
                           u10, v10, t2, q2, regime)
      USE WRF_TOOLS

      IMPLICIT NONE

      REAL, INTENT (in)                :: ps , ts , qs , us, vs
      REAL, INTENT (in)                :: ps2, ts2, qs2, psfc, tg
      REAL, INTENT (in)                :: hs, roughness, xland

      REAL,    INTENT (out)            :: regime
      REAL   , INTENT (out)            :: u10, v10, t2, q2

! Maximum number of iterations in computing psim, psih

      INTEGER, PARAMETER :: k_iteration = 10
!      INTEGER, PARAMETER :: k_iteration = 1

! h10 is the height of 10m where the wind observed
! h2  is the height of 2m where the temperature and 
!        moisture observed.

      REAL, PARAMETER :: h10 = 10., h2 = 2.
!
! Default roughness over the land

      REAL, PARAMETER :: zint0 = 0.01
!
! Von Karman constant

      REAL, PARAMETER :: k_kar = 0.4
!
! Working variables

      REAL :: Vc2, Va2, V2
      REAL :: rib, rr, rz, r2, rcp, xx, yy, cc
      REAL :: psiw, psiz, mol, ust, hol, holz, hol2
      REAL :: psim, psimz, psim2, psih, psihz, psih2
      REAL :: psit, psitz, psit2, psiq, psiqz, psiq2
      REAL :: gzsoz0, gz10oz0, gz2oz0
      REAL :: eg, qg, tvg, tvs, tvs2, tv
      REAL :: ths, thg, thvs, thvg, thvs2
      REAL :: dz, arg, e1, zq0, z0
      INTEGER :: k, nn, nz, n2

      REAL, PARAMETER :: ka = 2.4E-5

!-----------------------------------------------------------------------------!
      rcp = Rd/Cp

! 1 Compute the roughness length based upon season and land use 
! =====================================

! 1.1 Define the rouhness length
!     -----------------

      z0 = roughness

      if( z0 < 0.0001) z0 = 0.0001

! 1.2 Define the rouhgness length for moisture
!     -----------------

      if ( xland .ge. 1.5 ) then
        zq0 = z0
       else 
        zq0 =  zint0
      endif 

! 1.3 Define the some constant variable for psi
!     -----------------

      gzsoz0 = log(hs/z0)

      gz10oz0 = log(h10/z0) 

      gz2oz0 = log(h2/z0)


! 2. Calculate the virtual temperature
! =====================================

! 2.1 Compute Virtual temperature on the lowest half sigma level
!     ---------------------------------------------------------

      tvs  = ts * (1. + 0.608 * qs)

! 2.2 Compute Virtual temperature on the second lowest half sigma level
!     -----------------------------------------------------------------

      tvs2  = ts2 * (1. + 0.608 * qs2)

! 2.3 Convert ground virtual temperature assuming it's saturated
!     ----------------------------------------------------------

      call DA_TP_To_Qs( tg, psfc, eg, qg )
      tvg  = tg * (1. + 0.608 * qg)

! 3.  Compute the potential temperature
! ======================================

! 3.1 Potential temperature on the lowest half sigma level
!     ----------------------------------------------------

      ths  = ts * (1000. / (ps/100.)) ** rcp

! 3.2 Potential temperature at the ground
!     -----------------------------------

      thg  = tg * (1000. / (psfc/100.)) ** rcp


! 4. Virtual potential temperature
! ==================================

! 4.1 Virtual potential temperature on the lowest half sigma level
!     ------------------------------------------------------------

      thvs = tvs * (1000. / (ps/100.)) ** rcp

! 4.2 Virtual potential temperature on the second lowest half sigma level
!     ------------------------------------------------------------------

      thvs2 = tvs2 * (1000. / (ps2/100.)) ** rcp

! 4.3 Virtual potential temperature at ground
!     ---------------------------------------

      thvg = tvg * (1000. / (psfc/100.)) ** rcp


! 5.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
! =================================================

! 5.1 Velocity
!     --------
!
!     Wind speed:

      Va2 =   us*us + vs*vs
!  
!     Convective velocity:

      if ( thvg >= thvs ) then
         Vc2 = 4. * (thvg - thvs)
        else
         Vc2 = 0.
      endif
!
      V2  = Va2 + Vc2

! 5.2 Bulk richardson number
!     ----------------------

      rib = (g * hs / ths) * (thvs - thvg) / V2

! 6.  CALCULATE PSI BASED UPON REGIME
! =======================================

! 6.1 Stable conditions (REGIME 1)
!     ---------------------------

      IF       (rib .GE. 0.2) THEN
         regime = 1.1
         psim = -10.*gzsoz0
         psimz = -10.*gz10oz0
         psim2 = -10.*gz2oz0
         psim = max(psim,-10.)
         psimz = max(psimz,-10.)
         psim2 = max(psim2,-10.)
         psih = psim
         psihz = psimz
         psih2 = psim2

! 6.2 Mechanically driven turbulence (REGIME 2)
!     ------------------------------------------

      ELSE IF ((rib .LT. 0.2) .AND. (rib .GT. 0.)) THEN

         regime = 2.1
         psim = ( -5. * rib ) * gzsoz0 / (1.1 - 5.*rib)
         psimz = ( -5. * rib ) * gz10oz0 / (1.1 - 5.*rib)
         psim2 = ( -5. * rib ) * gz2oz0 / (1.1 - 5.*rib)

         psim = max(psim,-10.)
         psimz = max(psimz,-10.)
         psim2 = max(psim2,-10.)
         psih = psim
         psihz = psimz
         psih2 = psim2

! 6.3 Unstable Forced convection (REGIME 3)
!     -------------------------------------

      ELSE IF ((rib .EQ. 0.) .or. (rib.LT.0.0 .and. thvs2.GT.thvs)) THEN
         regime = 3.1
         psim = 0.
         psimz = 0.
         psim2 = 0. 
         psih = psim
         psihz = psimz
         psih2 = psim2


! 6.4 Free convection (REGIME 4)
!     --------------------------

      ELSE
        regime = 4.1

!      Calculate psi m and pshi h using iteration method
       
        psim = 0.
        psih = 0.
        cc = 2. * atan(1.0)

!        do k = 1 , k_iteration

! 6.4.1  Calculate   ust, m/L (mol), h/L (hol)
!        --------------------------

!       Friction speed

          ust = k_kar * sqrt(v2) /( gzsoz0 - psim)

!       Heat flux factor

          mol = k_kar * (ths - thg )/( gzsoz0 - psih)

!       Ratio of PBL height to Monin-Obukhov length

          if ( ust .LT. 0.01 ) then
             hol = rib * gzsoz0
           else
             hol = k_kar * g * hs * mol / ( ths * ust * ust )
          endif

! 6.4.2  Calculate n, nz, R, Rz
!        --------------------------

          hol = min(hol,0.)
          hol = max(hol,-10.)

          holz = (h10 / hs) * hol
          holz = min(holz,0.)
          holz = max(holz,-10.)

          hol2 = (h2 / hs) * hol
          hol2 = min(hol2,0.)
          hol2 = max(hol2,-10.)

! 6.4.3 Calculate Psim & psih
!        --------------------------

!       Using the look-up table:
!          nn = int( -100. * hol )
!          rr = ( -100. * hol ) - nn
!          psim = psimtb(nn) + rr * ( psimtb(nn+1) - psimtb(nn))
!          psih = psihtb(nn) + rr * ( psihtb(nn+1) - psihtb(nn))
!       Using the continuous function:
          xx = (1. - 16. * hol) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psim = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psih = 2. * yy

!       Using the look-up table:
!          nz = int( -100. * holz )
!          rz = ( -100. * holz ) - nz
!          psimz = psimtb(nz) + rz * ( psimtb(nz+1) - psimtb(nz))
!          psihz = psihtb(nz) + rz * ( psihtb(nz+1) - psihtb(nz))
!       Using the continuous function:
          xx = (1. - 16. * holz) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psimz = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psihz = 2. * yy

!       Using the look-up table:
!          n2 = int( -100. * hol2 )
!          r2 = ( -100. * hol2 ) - n2
!          psim2 = psimtb(n2) + r2 * ( psimtb(n2+1) - psimtb(n2))
!          psih2 = psihtb(n2) + r2 * ( psihtb(n2+1) - psihtb(n2))
!       Using the continuous function:
          xx = (1. - 16. * hol2) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psim2 = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psih2 = 2. * yy

!      enddo 

! 6.4.4 Define the limit value for psim & psih
!        --------------------------

       psim = min(psim,0.9*gzsoz0)
       psimz = min(psimz,0.9*gz10oz0)
       psim2 = min(psim2,0.9*gz2oz0)
       psih = min(psih,0.9*gzsoz0)
       psihz = min(psihz,0.9*gz10oz0)
       psih2 = min(psih2,0.9*gz2oz0)

      ENDIF  ! Regime

! 7.  CALCULATE PSI FOR WIND, TEMPERATURE AND MOISTURE
! =======================================

      psiw = gzsoz0 - psim
      psiz = gz10oz0 - psimz
      psit = gzsoz0 - psih 
      psit2 = gz2oz0 - psih2

!     Friction speed
      ust = k_kar * sqrt(v2) /( gzsoz0 - psim)

      psiq  = log(k_kar*ust*hs/ka + hs / zq0 ) - psih
      psiq2 = log(k_kar*ust*h2/ka + h2 / zq0 ) - psih2

! 8.  CALCULATE 10M WIND, 2M TEMPERATURE AND MOISTURE
! =======================================

      u10 = us * psiz / psiw
      v10 = vs * psiz / psiw
      t2 = ( thg + ( ths - thg )*psit2/psit)*((psfc/100.)/1000.)**rcp
      q2 = qg + (qs - qg)*psiq2/psiq  

    END SUBROUTINE SFC_WTQ



! ---------------
      SUBROUTINE cal_10mu ( psfc, tg, ps, ts, qs, us, vs, ps2, ts2, qs2, &
                           hs, roughness, xland, &
                           u10, v10, t2, q2, regime)
      USE WRF_TOOLS

      IMPLICIT NONE

      REAL, INTENT (in)                :: ps , ts , qs , us, vs
      REAL, INTENT (in)                :: ps2, ts2, qs2, psfc, tg
      REAL, INTENT (in)                :: hs, roughness, xland

      REAL,    INTENT (out)            :: regime
      REAL   , INTENT (out)            :: u10, v10, t2, q2

! Maximum number of iterations in computing psim, psih

      INTEGER, PARAMETER :: k_iteration = 10
!      INTEGER, PARAMETER :: k_iteration = 1

! h10 is the height of 10m where the wind observed
! h2  is the height of 2m where the temperature and 
!        moisture observed.

      REAL, PARAMETER :: h10 = 10., h2 = 2.
!
! Default roughness over the land

      REAL, PARAMETER :: zint0 = 0.01
!
! Von Karman constant

      REAL, PARAMETER :: k_kar = 0.4
!
! Working variables

      REAL :: Vc2, Va2, V2
      REAL :: rib, rr, rz, r2, rcp, xx, yy, cc
      REAL :: psiw, psiz, mol, ust, hol, holz, hol2
      REAL :: psim, psimz, psim2, psih, psihz, psih2
      REAL :: psit, psitz, psit2, psiq, psiqz, psiq2
      REAL :: gzsoz0, gz10oz0, gz2oz0
      REAL :: eg, qg, tvg, tvs, tvs2, tv
      REAL :: ths, thg, thvs, thvg, thvs2
      REAL :: dz, arg, e1, zq0, z0
      INTEGER :: k, nn, nz, n2

      REAL, PARAMETER :: ka = 2.4E-5

!-----------------------------------------------------------------------------!
      rcp = Rd/Cp

! 1 Compute the roughness length based upon season and land use 
! =====================================

! 1.1 Define the rouhness length
!     -----------------

      z0 = roughness

      if( z0 < 0.0001) z0 = 0.0001

! 1.2 Define the rouhgness length for moisture
!     -----------------

      if ( xland .ge. 1.5 ) then
        zq0 = z0
       else 
        zq0 =  zint0
      endif 

! 1.3 Define the some constant variable for psi
!     -----------------

      gzsoz0 = log(hs/z0)

      gz10oz0 = log(h10/z0) 

      gz2oz0 = log(h2/z0)


! 2. Calculate the virtual temperature
! =====================================

! 2.1 Compute Virtual temperature on the lowest half sigma level
!     ---------------------------------------------------------

      tvs  = ts * (1. + 0.608 * qs)

! 2.2 Compute Virtual temperature on the second lowest half sigma level
!     -----------------------------------------------------------------

      tvs2  = ts2 * (1. + 0.608 * qs2)

! 2.3 Convert ground virtual temperature assuming it's saturated
!     ----------------------------------------------------------

      call DA_TP_To_Qs( tg, psfc, eg, qg )
      tvg  = tg * (1. + 0.608 * qg)

! 3.  Compute the potential temperature
! ======================================

! 3.1 Potential temperature on the lowest half sigma level
!     ----------------------------------------------------

      ths  = ts * (1000. / (ps/100.)) ** rcp

! 3.2 Potential temperature at the ground
!     -----------------------------------

      thg  = tg * (1000. / (psfc/100.)) ** rcp


! 4. Virtual potential temperature
! ==================================

! 4.1 Virtual potential temperature on the lowest half sigma level
!     ------------------------------------------------------------

      thvs = tvs * (1000. / (ps/100.)) ** rcp

! 4.2 Virtual potential temperature on the second lowest half sigma level
!     ------------------------------------------------------------------

      thvs2 = tvs2 * (1000. / (ps2/100.)) ** rcp

! 4.3 Virtual potential temperature at ground
!     ---------------------------------------

      thvg = tvg * (1000. / (psfc/100.)) ** rcp


! 5.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
! =================================================

! 5.1 Velocity
!     --------
!
!     Wind speed:

      Va2 =   us*us + vs*vs
!  
!     Convective velocity:

      if ( thvg >= thvs ) then
         Vc2 = 4. * (thvg - thvs)
        else
         Vc2 = 0.
      endif
!
      V2  = Va2 + Vc2

! 5.2 Bulk richardson number
!     ----------------------

      rib = (g * hs / ths) * (thvs - thvg) / V2

! 6.  CALCULATE PSI BASED UPON REGIME
! =======================================

! 6.1 Stable conditions (REGIME 1)
!     ---------------------------

      IF       (rib .GE. 0.2) THEN
         regime = 1.1
         psim = -10.*gzsoz0
         psimz = -10.*gz10oz0
         psim2 = -10.*gz2oz0
         psim = max(psim,-10.)
         psimz = max(psimz,-10.)
         psim2 = max(psim2,-10.)
         psih = psim
         psihz = psimz
         psih2 = psim2

! 6.2 Mechanically driven turbulence (REGIME 2)
!     ------------------------------------------

      ELSE IF ((rib .LT. 0.2) .AND. (rib .GT. 0.)) THEN

         regime = 2.1
         psim = ( -5. * rib ) * gzsoz0 / (1.1 - 5.*rib)
         psimz = ( -5. * rib ) * gz10oz0 / (1.1 - 5.*rib)
         psim2 = ( -5. * rib ) * gz2oz0 / (1.1 - 5.*rib)

         psim = max(psim,-10.)
         psimz = max(psimz,-10.)
         psim2 = max(psim2,-10.)
         psih = psim
         psihz = psimz
         psih2 = psim2

! 6.3 Unstable Forced convection (REGIME 3)
!     -------------------------------------

      ELSE IF ((rib .EQ. 0.) .or. (rib.LT.0.0 .and. thvs2.GT.thvs)) THEN
         regime = 3.1
         psim = 0.
         psimz = 0.
         psim2 = 0. 
         psih = psim
         psihz = psimz
         psih2 = psim2


! 6.4 Free convection (REGIME 4)
!     --------------------------

      ELSE
        regime = 4.1

!      Calculate psi m and pshi h using iteration method
       
        psim = 0.
        psih = 0.
        cc = 2. * atan(1.0)

!        do k = 1 , k_iteration

! 6.4.1  Calculate   ust, m/L (mol), h/L (hol)
!        --------------------------

!       Friction speed

          ust = k_kar * sqrt(v2) /( gzsoz0 - psim)

!       Heat flux factor

          mol = k_kar * (ths - thg )/( gzsoz0 - psih)

!       Ratio of PBL height to Monin-Obukhov length

          if ( ust .LT. 0.01 ) then
             hol = rib * gzsoz0
           else
             hol = k_kar * g * hs * mol / ( ths * ust * ust )
          endif

! 6.4.2  Calculate n, nz, R, Rz
!        --------------------------

          hol = min(hol,0.)
          hol = max(hol,-10.)

          holz = (h10 / hs) * hol
          holz = min(holz,0.)
          holz = max(holz,-10.)

          hol2 = (h2 / hs) * hol
          hol2 = min(hol2,0.)
          hol2 = max(hol2,-10.)

! 6.4.3 Calculate Psim & psih
!        --------------------------

!       Using the look-up table:
!          nn = int( -100. * hol )
!          rr = ( -100. * hol ) - nn
!          psim = psimtb(nn) + rr * ( psimtb(nn+1) - psimtb(nn))
!          psih = psihtb(nn) + rr * ( psihtb(nn+1) - psihtb(nn))
!       Using the continuous function:
          xx = (1. - 16. * hol) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psim = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psih = 2. * yy

!       Using the look-up table:
!          nz = int( -100. * holz )
!          rz = ( -100. * holz ) - nz
!          psimz = psimtb(nz) + rz * ( psimtb(nz+1) - psimtb(nz))
!          psihz = psihtb(nz) + rz * ( psihtb(nz+1) - psihtb(nz))
!       Using the continuous function:
          xx = (1. - 16. * holz) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psimz = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psihz = 2. * yy

!       Using the look-up table:
!          n2 = int( -100. * hol2 )
!          r2 = ( -100. * hol2 ) - n2
!          psim2 = psimtb(n2) + r2 * ( psimtb(n2+1) - psimtb(n2))
!          psih2 = psihtb(n2) + r2 * ( psihtb(n2+1) - psihtb(n2))
!       Using the continuous function:
          xx = (1. - 16. * hol2) ** 0.25
          yy = log((1.+xx*xx)/2.)
          psim2 = 2. * log((1.+xx)/2.) + yy - 2. * atan(xx) + cc
          psih2 = 2. * yy

!      enddo 

! 6.4.4 Define the limit value for psim & psih
       psim = min(psim,0.9*gzsoz0)
       psimz = min(psimz,0.9*gz10oz0)

      ENDIF  ! Regime

! 7.  CALCULATE PSI FOR WIND, TEMPERATURE AND MOISTURE
      psiw = gzsoz0 - psim
      psiz = gz10oz0 - psimz

! 8.  CALCULATE 10M WIND, 2M TEMPERATURE AND MOISTURE
! =======================================

      u10 = us * psiz / psiw
     endD SUBROUTINE cal_10mu
