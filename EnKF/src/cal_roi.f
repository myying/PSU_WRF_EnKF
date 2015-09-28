!========================================================================================
   subroutine cal_hroi ( instrument, grid_id, iroi, ngxn ) 

   implicit none

   character (len=8), intent(in)          :: instrument
   integer, intent(in)                    :: grid_id, iroi
   integer, intent(inout)                 :: ngxn

!   if ( instrument == 'Radar   ' .or. instrument == 'aircft  ' .or.  &
!        instrument == 'metar   ' .or. instrument == 'satwnd  ' ) then
!!Successive covariance localization (SCL) technique:
!!http://hfip.psu.edu/fuz4/2011/WengZhang2011MWR.pdf
   if ( instrument == 'Radar   ' ) then
        if ( grid_id == 1 ) then
             ngxn = 1
        else if ( grid_id == 2 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
        else if ( grid_id == 3 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
             if( mod(iroi,9) == 1 ) ngxn = 9
        else if ( grid_id == 4 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
             if( mod(iroi,9) == 1 ) ngxn = 9
             if( mod(iroi,27) == 1) ngxn = 27
        endif
   else
!        ngxn = 3**(grid_id-1)
! need to adjust HROI in namelist.enkf to match the (hroi_n_grid)*(grid spacing)=hroi_km
! namelist.enkf sets hroi_n_grid
       ngxn = 1
   endif

   end subroutine cal_hroi
!========================================================================================
   subroutine corr(dx,dy,dz,ngx,ngz,corr_coef)
! This is alternative routine to calculate corr_coef, if schur_matrix
! requires too much memory to be calculated before hand.

  implicit none

  real, intent(in)    :: dx,dy,dz        !dx: x distance(grids) from model grid to obs
  integer, intent(in) :: ngx, ngz        !ngx: horrizontal cutting off distances
  integer :: i,j,k
  real, intent(out)   :: corr_coef

  integer :: horradi
  real :: k1
  real :: comp_cov_factor
  real :: horrad, distance, distanceh

! calculate horizontal radius at height k
     horrad = (real(ngx)/real(ngz))*sqrt(real(ngz**2. - dz**2.))
     horradi = int(horrad)   ! grid points within the radius

! equivalence of k in terms of dx  ! added FZ 2004/09/09
     k1 = dz * real(ngx) / real(ngz)

        distanceh = sqrt(real(dx**2. + dy**2.))   ! hor. distance from z-axis

        if ((dx==0.) .and. (dy==0.) .and. (dz==0.)) then
           corr_coef = 1.

        else if (distanceh<=horrad) then

           distance  = sqrt( real(dx**2. + dy**2. + k1**2.))   ! 3-d distance from obs

           corr_coef  = comp_cov_factor(dble(distance),dble(ngx/2.))

        else
           corr_coef  = 0.

        end if

   end subroutine corr
!========================================================================================
   subroutine corr_matrix(nx,nz,cmatr)

! written by Altug Aksoy, 05/20/2003

! this subroutine computes the coefficient matrix to be used for
! compact correlation function calculations. the matrix is computed
! once and stored so that it can be refered to when calculating the gain.

  implicit none

  integer, intent(in) :: nx,nz
  integer :: i,j,k
  integer :: centerx, centerz, horradi, totradi

  real :: k1
  real , intent(out), dimension(2*nx+1,2*nx+1,2*nz+1) :: cmatr
  real :: comp_cov_factor
  real :: horrad, distance, distanceh, term1, term2, totrad, schur

  cmatr = 0.

  centerx = nx+1   ! location of origin (obs) within the matrix
  centerz = nz+1

  do k = -nz,nz

! calculate horizontal radius at height k
     horrad = (real(nx)/real(nz))*sqrt(real(nz**2. - k**2.))
     horradi = int(horrad)   ! grid points within the radius

! equivalence of k in terms of dx  ! added FZ 2004/09/09
     k1 = real(k) * real(nx) / real(nz)

     do j = -horradi,horradi
     do i = -horradi,horradi

        distanceh = sqrt(real(i**2. + j**2.))   ! hor. distance from z-axis

        if ((i==0) .and. (j==0) .and. (k==0)) then
           cmatr(centerx,centerx,centerz) = 1.

        else if (distanceh<=horrad) then

           distance  = sqrt( real(i**2. + j**2. + k1**2.))   ! 3-d distance from obs
!           distance  = sqrt( real(i**2. + j**2. ))   ! 2-d distance from obs

           schur  = comp_cov_factor(dble(distance),dble(nx/2.))
           cmatr(centerx+i, centerx+j,centerz+k) = schur

        end if

      enddo
      enddo

   enddo

   end subroutine corr_matrix
!========================================================================================
   subroutine corr_matrix_h(nx,cmatr)
! this subroutine computes the coefficient matrix to be used for
! compact correlation function calculations. the matrix is computed
! once and stored so that it can be refered to when calculating the gain.

  implicit none

  integer, intent(in) :: nx
  integer :: i,j,k
  integer :: centerx, horradi, totradi

  real , intent(out), dimension(2*nx+1,2*nx+1) :: cmatr
  real :: comp_cov_factor
  real :: horrad, distance, distanceh, term1, term2, totrad, schur

  cmatr = 0.

  centerx = nx+1   ! location of origin (obs) within the matrix
     do j = -nx, nx
     do i = -nx, nx
       distanceh = sqrt(real(i**2. + j**2.))   ! hor. distance from z-axis

        if ((i==0) .and. (j==0)) then
           cmatr(centerx,centerx) = 1.

        else if (distanceh<=real(nx)) then
           distance  = sqrt( real(i**2. + j**2. ))   ! 2-d distance from obs

           schur  = comp_cov_factor(dble(distance),dble(nx/2.))
           cmatr(centerx+i, centerx+j) = schur

        end if

     enddo
     enddo

  end subroutine corr_matrix_h
!==============================================================================
   function comp_cov_factor(z_in, c)

   implicit none

   real comp_cov_factor
   double precision z_in, c
   double precision z, r

!  Computes a covariance cutoff function from Gaspari and Cohn
!  (their eqn. 4.10) QJRMS, 125, 723-757.

!  z_in is the distance while c is the cutoff distance.
!  For distances greater than 2c, the cov_factor returned goes to 0.

   z = dabs(z_in)
   r = z / c

   if(z >= 2*c) then
      comp_cov_factor = 0.0
   else if(z >= c .and. z < 2*c) then
      comp_cov_factor =                                               &
          ( ( ( ( r/12.  -0.5 )*r  +0.625 )*r +5./3. )*r  -5. )*r     &
                                                 + 4. - 2./(3.*r)
   else
      comp_cov_factor =                                               &
          ( ( ( -0.25*r +0.5 )*r +0.625 )*r  -5./3. )*r**2 + 1.
   endif

   end function comp_cov_factor
!==============================================================================
