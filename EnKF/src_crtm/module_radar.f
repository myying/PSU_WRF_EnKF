module radar
!------------------------------------------------------------------------------
!--  This module 
!--
!------------------------------------------------------------------------------

  use constants
  use namelist_define
  use mapinfo_define
  use obs_define
  use mpi_module
  use map_utils
  use netcdf
  use wrf_tools

  contains

!=======================================================================================
  subroutine get_radar_info ( filename )

  implicit none

  character(len=*), intent(in) :: filename
  character(len=400)           :: dat

  integer                      :: i, j, k, len, nr

  integer                      :: elva, hgt
  
  if ( my_proc_id == 0 ) then
     write(6, *)'   '
     write(6, *)'---------------------------------------------------'
     write(6, *)'.... Radar Station Information from: ', filename
  endif
  len = len_trim( filename )
  open(10, file=filename(1:len), status='old' ) 

  do nr = 1, raw%radar_stn_num

     read(10,'(a)')dat
  
     i = index( dat, "station idn=" )
     read( dat(i+13:i+18), '(i6)')  raw%radar(nr)%radar_stn%idn 
    
     i = index( dat, "id=" )
     raw%radar(nr)%radar_stn%idc = dat(i+4:i+7)

     i = index( dat, "name=" )
     j = index( dat(i:i+49), """ " )
     raw%radar(nr)%radar_stn%name = dat(i+6:i+j-2)

     i = index( dat, "lat=" )
     j = index( dat(i:i+15), """ " )
     read ( dat(i+5:i+j-2), '(f)') raw%radar(nr)%radar_stn%lat

     i = index( dat, "lon=" )
     j = index( dat(i:i+15), """ " )
     read ( dat(i+5:i+j-2), '(f)') raw%radar(nr)%radar_stn%lon

     i = index( dat, "elev=" )
     j = index( dat(i:i+15), """ " )
     read ( dat(i+6:i+j-2), '(i)') elva

    read(10,'(a)')dat

     i = index( dat, "hgt=" ) 
     j = index( dat(i:i+15), """ " )
     read ( dat(i+5:i+j-2), '(i)') hgt

     raw%radar(nr)%radar_stn%elv = real(elva) + real(hgt) 

     i = index( dat, "scan=" )
     j = index( dat(i:i+15), """ " )
     read ( dat(i+6:i+j-2), '(i)') raw%radar(nr)%radar_stn%levels
     
     i = index( dat, "elevation=" )
     j = index( dat(i:i+100), """ " )
     read ( dat(i+11:i+j-2), *)   &
          raw%radar(nr)%radar_stn%elv_ang(1:raw%radar(nr)%radar_stn%levels)

     i = index( dat, "mindis=" ) 
     j = index( dat(i:i+20), """ " )
     read ( dat(i+8:i+j-2), '(f)') raw%radar(nr)%radar_stn%min_dis

     i = index( dat, "maxdis=" ) 
     j = index( dat(i:i+20), """ " )
     read ( dat(i+8:i+j-2), '(f)') raw%radar(nr)%radar_stn%max_dis

     i = index( dat, "ddis=" )
     j = index( dat(i:i+20), """ " )
     read ( dat(i+6:i+j-2), '(f)') raw%radar(nr)%radar_stn%d_dis

     i = index( dat, "dazm=" )
     j = index( dat(i:i+20), """ " )
     read ( dat(i+6:i+j-2), '(f)') raw%radar(nr)%radar_stn%d_azim

     i = index( dat, "rf_err=" )
     j = index( dat(i:i+20), """ " )
     read ( dat(i+8:i+j-2), '(f)') raw%radar(nr)%radar_stn%err_rf

     i = index( dat, "rv_err=" )
     j = index( dat(i:i+20), """ " )
     read ( dat(i+8:i+j-2), '(f)') raw%radar(nr)%radar_stn%err_rv

     i = len_trim(raw%radar(nr)%radar_stn%idc)
     if ( my_proc_id == 0 ) then
        write(6,'(a,i6.6,a)')" NEXDAR:",raw%radar(nr)%radar_stn%idn,":"//raw%radar(nr)%radar_stn%idc(1:i)//":"//raw%radar(nr)%radar_stn%name
        write(6,'(a,f7.2,a,f7.2)')" lat = " , raw%radar(nr)%radar_stn%lat, " lon = ", raw%radar(nr)%radar_stn%lon
        write(6,'(a,f6.1,a)')" The height of the radar antenna above sea level : ", real(raw%radar(nr)%radar_stn%elv), "(m)"
        write(6,'(a,20f5.1)')" Elevation : ", raw%radar(nr)%radar_stn%elv_ang(1:raw%radar(nr)%radar_stn%levels)
     endif
      
  end do

     close(10)


  end subroutine get_radar_info

!=======================================================================================
  subroutine radar_position ( proj, nr )

  implicit none
!
  type(proj_info), intent(in)         :: proj
  integer, intent(in)                 :: nr     ! the nrth radar

  real                                :: xla, xlo, is, js

!------------------------------------------------------------------------------
!-- compute radar position according to wrf domain grid
  xla = raw%radar(nr)%radar_stn%lat
  xlo = raw%radar(nr)%radar_stn%lon

  call latlon_to_ij( proj, xla, xlo, is, js ) 

  raw%radar(nr)%radar_stn%i_radar = is
  raw%radar(nr)%radar_stn%j_radar = js

  if ( my_proc_id == 0 ) write(6,'(a,i1,a,f7.2,a,f7.2,a)')" The ",nr,"th Radar Position according to Model : (", &
           is,",",js,")"

!------------------------------------------------------------------------------

  end subroutine radar_position  
   
!=======================================================================================
  subroutine calc_radar_data_point_s_h (s, h, r, elev, h0)

  implicit none

  real, intent(out) :: s           ! surface distance from radar to data point
  real, intent(out) :: h           ! vertical height above sea level for data point
  real, intent(in)  :: r           ! range from radar to data point
  real, intent(in)  :: elev        ! elevation angle
  real, intent(in)  :: h0          ! radar base height above sea level (unit: m)
  real              :: ke

  ke = 4./3.

!------------------------------------------------------------------------------
! 
  h = sqrt( r**2. + ( ke*R_earth )**2. + 2.*r*ke*R_earth*sin(elev) ) &
           - ke*R_earth + h0 / 1000. 

  s = ke*R_earth*asin( (r*cos(elev))/(ke*R_earth+h) ) 

!------------------------------------------------------------------------------

  end subroutine calc_radar_data_point_s_h

!=======================================================================================
  subroutine wsr88d_ideal_data_point( proj )

  implicit none

  type(proj_info), intent(in)         :: proj

!.work arrays and variables
  integer     :: n, k, kl, nr, na, i, j, obs_num, allocation_status
  real        :: r_dis, alpha2, alpha0, rand, theta, s, h 
  real, dimension(:), allocatable  :: obs_range, obs_azimuth, obs_elev, obs_irp, obs_jrp, obs_hrp
!  real, dimension(5000)  :: obs_range, obs_azimuth, obs_elev, obs_irp, obs_jrp, obs_hrp
  real        :: rx, ry, ir1, jr1, random
  integer     :: rand00

!------------------------------------------------------------------------------
  do_radar_stn           :  do n = 1, raw%radar_stn_num

!------------------------------------------------------------------------------
!... allocatable varies
     kl = raw%radar(n)%radar_stn%levels
     nr = int( ( raw%radar(n)%radar_stn%max_dis -    &
                 raw%radar(n)%radar_stn%min_dis ) /  &
                 raw%radar(n)%radar_stn%d_dis + 1   )     
     na = int( (359. - 0.) / raw%radar(n)%radar_stn%d_azim + 1 )
     allocate (obs_range (kl*nr*na), stat = allocation_status)
     allocate (obs_azimuth (kl*nr*na), stat = allocation_status)
     allocate (obs_elev (kl*nr*na), stat = allocation_status)
     allocate (obs_irp (kl*nr*na), stat = allocation_status)
     allocate (obs_jrp (kl*nr*na), stat = allocation_status)
     allocate (obs_hrp (kl*nr*na), stat = allocation_status)
     
!... pick out data points 
     obs_num = 0
     do_elevation_angles_loop :  do k = 1, raw%radar(n)%radar_stn%levels

        do_radial_range_loop  : do j = 1, nr

!......... define radial position
           r_dis = raw%radar(n)%radar_stn%min_dis +       &
                   raw%radar(n)%radar_stn%d_dis *         &
                   ( float(j) - 0.5 - float(k) / float(kl) * (-1.)**k ) 

!......... define azimuth angle
           alpha2 = raw%radar(n)%radar_stn%d_azim *       &
                    raw%radar(n)%radar_stn%max_dis / r_dis
           alpha2 = min(alpha2, 89.)
           na = int( 360. / alpha2 )
           alpha0 = alpha2*(1.-random(rand00))
           alpha0 = alpha0*(1.-float(k)/float(2*kl+1))*(-1.)**(k+1)
           alpha0 = alpha0*(1.+random(rand00))

           do_azimuth_loop     : do i = 1, na

              theta = alpha0 + alpha2 * float(i-1)
              theta = theta * rad_per_deg

!............ compute the surface distance and ASL height 
              call calc_radar_data_point_s_h (s, h, r_dis, &
                   raw%radar(n)%radar_stn%elv_ang(k)*rad_per_deg, raw%radar(n)%radar_stn%elv)

!............ compute data pint's position in model domain -- i, j
              rx = s * sin( theta )
              ry = s * cos( theta )
              ir1 = raw%radar(n)%radar_stn%i_radar + ( rx / proj%dx *1000. )
              jr1 = raw%radar(n)%radar_stn%j_radar + ( ry / proj%dx *1000. )

!............ screen output
              if( k == raw%radar(n)%radar_stn%levels .and. j == nr/5 .and. i == na ) then
                 write(6,'(a,i1,a,f7.1,a,f7.1,a,f4.1,a,f8.2,a,f8.2,a)')                     &
                     " The ",n,"th Radar's data point(", r_dis, "km,", theta, "deg,",       &
                     raw%radar(n)%radar_stn%elv_ang(k),                                     &
                     "deg)'s position according to Model : (", ir1,",",jr1,")"
              endif

!............ save to temp arries
              if ( (ir1 .ge. 2. .and. ir1 .le. (proj%nx -1) ) .and.                &
                   (jr1 .ge. 2. .and. jr1 .le. (proj%ny -1) ) .and.                &
                   (h .ge. 0.001 .and. h .le. 5.)       ) then
                 obs_num = obs_num + 1
                 obs_range(obs_num) = r_dis
                 obs_azimuth(obs_num) = theta
                 obs_elev(obs_num) = raw%radar(n)%radar_stn%elv_ang(k)
                 obs_irp(obs_num) = ir1
                 obs_jrp(obs_num) = jr1
                 obs_hrp(obs_num) = h
              endif

           end do do_azimuth_loop
        end do do_radial_range_loop
     end do do_elevation_angles_loop

!------------------------------------------------------------------------------
!... save position info. to raw%radar(n) type 
     raw%radar(n)%radar_stn%numObs = obs_num
     !allocate( raw%radar(n)%radar_point( obs_num ) )
     !allocate( raw%radar(n)%rf( obs_num ) )
     !allocate( raw%radar(n)%rv( obs_num ) )

     do i = 1, obs_num
        raw%radar(n)%radar_point(i)%ii   = obs_irp(i)
        raw%radar(n)%radar_point(i)%jj   = obs_jrp(i)
        raw%radar(n)%radar_point(i)%hgt  = obs_hrp(i)
        raw%radar(n)%radar_point(i)%rdis = obs_range(i)
        raw%radar(n)%radar_point(i)%azim = obs_azimuth(i)
        raw%radar(n)%radar_point(i)%elev = obs_elev(i)
     enddo

     deallocate ( obs_range )
     deallocate ( obs_azimuth )
     deallocate ( obs_elev )
     deallocate ( obs_irp )
     deallocate ( obs_jrp )
     deallocate ( obs_hrp )
!------------------------------------------------------------------------------
  end do do_radar_stn

  end subroutine wsr88d_ideal_data_point

!=======================================================================================

  subroutine wsr88d_ideal_data ( filename, ix, jx, kx, proj )

  implicit none

  type(proj_info), intent(in)         :: proj
  character (len=10), intent(in)      :: filename
  integer, intent(in)                 :: ix, jx, kx

  real, dimension(ix+1, jx, kx)       :: u
  real, dimension(ix, jx+1, kx)       :: v
  real, dimension(ix, jx, kx+1)       :: w
  real, dimension(ix, jx      )       :: xlong
  real, dimension(ix, jx, kx+1)       :: ph
  real, dimension(ix, jx, kx+1)       :: phb
  real                                :: rv

  character (len=10)                  :: wrf_file
  real                                :: ii, jj, kk, hh, i_radar, j_radar, elv
  integer                             :: unit, n, l

!------------------------------------------------------------------------------
  if ( my_proc_id == 0 ) then
      write(6, *)'   '
      write(6, *)'---------------------------------------------------'
      write(6, *)'.... Use Model Output to Creat Simulated Observation  ....'
  endif

!------------------------------------------------------------------------------
! compute radar data point's i, j, h, azimuth, range, elevation

  call wsr88d_ideal_data_point( proj )

!------------------------------------------------------------------------------
  read( filename, '(5x, i5)' ) unit
  write( wrf_file, '(a5, i5.5)' ) filename(1:5), unit-1 

!--------------------------------------------------------------------------------------
! get every radar and every data point's ii, jj, kk in 1st guest domain
  if ( use_radar_rv .or. use_radar_rf ) then
     call get_variable3d ( wrf_file, 'PH        ', ix, jx, kx+1, 1, ph  ) 
     call get_variable3d ( wrf_file, 'PHB       ', ix, jx, kx+1, 1, phb ) 

     do_radar_stn        : do n = 1, raw%radar_stn_num 
     do_radar_data_point : do l = 1, raw%radar(n)%radar_stn%numObs
        ii = raw%radar(n)%radar_point(l)%ii
        jj = raw%radar(n)%radar_point(l)%jj
        hh = raw%radar(n)%radar_point(l)%hgt 
        call calc_radar_data_point_kk ( ix, jx, kx, ph, phb, ii, jj, hh, kk )
        raw%radar(n)%radar_point(l)%kk = kk

        if ( n == raw%radar_stn_num .and. l == raw%radar(n)%radar_stn%numObs ) then
           if ( my_proc_id == 0 )write(6, '(a,i2,a,i7,a,3f7.2)')' The ', n, 'th radar ',   &
                 l, 'th data point''s ii jj kk : ', ii, jj, kk 
        endif

     end do do_radar_data_point
     end do do_radar_stn

  endif
 
! so every obser data point position ( ii, jj, kk ) is known

!--------------------------------------------------------------------------------------
  if ( use_radar_rv ) then

!    get model output u, v, w
     call get_variable3d ( wrf_file, 'U         ', ix+1, jx  , kx  , 1, u )     
     call get_variable3d ( wrf_file, 'V         ', ix, jx+1  , kx  , 1, v )     
     call get_variable3d ( wrf_file, 'W         ', ix, jx  , kx+1  , 1, w )     
     call get_variable2d ( wrf_file, 'XLONG     ', ix, jx  , 1,         xlong )     

     do n = 1, raw%radar_stn_num
     do l = 1, raw%radar(n)%radar_stn%numObs
        ii = raw%radar(n)%radar_point(l)%ii
        jj = raw%radar(n)%radar_point(l)%jj
        kk = raw%radar(n)%radar_point(l)%kk
        hh = raw%radar(n)%radar_point(l)%hgt
        i_radar = raw%radar(n)%radar_stn%i_radar
        j_radar = raw%radar(n)%radar_stn%j_radar
        elv     = raw%radar(n)%radar_stn%elv
 
        call calculate_rv ( ix, jx, kx, u, v, w, xlong, proj, ii, jj, kk,      &
                            hh, i_radar, j_radar, elv, rv)
        raw%radar(n)%rv(l) = rv

        if ( n == raw%radar_stn_num .and. l == raw%radar(n)%radar_stn%numObs ) then
           if ( my_proc_id == 0 )write(6, '(a,i2,a,i7,a,f7.2)')' The ', n, 'th radar ',   &
                 l, 'th data point''s radial velocity : ', rv 
        endif

     enddo
     enddo
  endif
  
!  if ( use_groundbase_radar .and. use_radar_rf  ) then
!     do n = 1, raw%radar_stn_num
!        call simulated_rf ( wrf_file, ix, jx, kx, proj, n )
!     enddo
!  endif

   end subroutine wsr88d_ideal_data

!=======================================================================================
  subroutine calc_radar_data_point_kk ( ix, jx, kx, ph, phb, ii, jj, hh, kk )

  implicit none

  integer, intent(in)                         :: ix, jx, kx
  real, intent(in)                            :: ii, jj, hh
  real, intent(out)                           :: kk

  real, dimension(ix, jx, kx+1), intent(in)   :: ph
  real, dimension(ix, jx, kx+1), intent(in)   :: phb
  
  real, dimension(2, 2, kx)                   :: z
  real, dimension(2, 2, kx+1)                 :: zw
  real, dimension(kx)                         :: zm
  real, dimension(kx+1)                       :: p

  integer                                     :: i, j, k, i1, j1, k0, k1, k2
  real                                        :: ire, jre, kre
  real                                        :: zm_min


!------------------------------------------------------------------------------
! get the leftdown point height around obser data point
  i1 = int ( ii )
  j1 = int ( jj )
  ire = ii - i1
  jre = jj - j1

  p (1:kx+1) = ( ph(i1,j1,1:kx+1) + phb(i1,j1,1:kx+1) ) / g
  zm(1:kx  ) = 0.5*( p(1:kx) + p(2:kx+1) )
  zm         = zm/1000.
  
!------------------------------------------------------------------------------
! find out the most close level around obser data point

!  zm(1:kx) = zm(1:kx) - hh    
  k0 = 1
  k1 = 1
  k2 = 2
!  zm_min = hh - zm(1)
!  do k = 2, kx
!     if ( zm(k) < zm_min ) then
!        zm_min = abs( zm(k) )
!        k0 = k
!     endif
!  enddo
!  if ( zm(k0) <= 0. )then
!     k1 = k0
!     k2 = k0 + 1
!  else if ( zm(k0) > 0. ) then
!     k2 = k0                         ! k1 is under k2
!     k1 = k0 - 1
!  endif 

  if ( zm(1) >= hh ) then
     kk = 1
     return
  else if ( zm(kx) <= hh ) then
     kk = kx
     return
  else
     k0 = 1
     do k =2, kx-1
        if ( zm(k) < hh .and. zm(k+1) >= hh ) then
           k0 = k-1
        endif
     enddo
     k1 = k0
     k2 = k1 + 1
  endif
  

!  write(*,*)'k1, k2 =', k1,k2

!------------------------------------------------------------------------------
! get 8 points' z value  
  kk_value   : do

     do_k_cycle : do k = k1, k2
     do_j_cycle : do j = j1, j1+1
     do_i_cycle : do i = i1, i1+1
        p(k:k+1) = ( ph(i,j,k:k+1) + phb(i,j,k:k+1) ) / g
        z(i-i1+1,j-j1+1,k) = 0.5*( p(k) + p(k+1) ) / 1000.
     end do do_i_cycle  
     end do do_j_cycle   

        zm(k) = ( z(1,1,k)*(1-ire) + z(2,1,k)*ire )*(1-jre) +         &
                ( z(1,2,k)*(1-ire) + z(2,2,k)*ire )*jre     
     end do do_k_cycle

     if ( hh >= zm(k1) .and. hh <= zm(k2) ) then
        kk = k1 + (hh - zm(k1))/(zm(k2)-zm(k1))
        exit kk_value
     else if ( hh < zm(k1) ) then
        k1 = k1 - 1
        k2 = k1 + 1
     else if ( hh > zm(k2) ) then
        k2 = k2 + 1
        k1 = k2 - 1
     endif

     if ( k1 < 1 ) then
        kk = 1
        exit kk_value
     endif
    
     if ( k2 > kx ) then
        kk = kx
        exit kk_value
     endif 

  end do kk_value
!------------------------------------------------------------------------------

  end subroutine calc_radar_data_point_kk

!=======================================================================================
  subroutine calculate_rv ( ix, jx, kx, ugrid, vgrid, wgrid, xlong, proj, ii, jj, kk,  &
                            hh, i_radar, j_radar, elv, rv)
  
  implicit none

  type(proj_info), intent(in)         :: proj                   ! 1st guest map info
  integer, intent(in)                 :: ix, jx, kx             ! 1st guest dimension
  real,  intent(in)                   :: ii, jj, kk, hh         ! radar data point position and height
  real,  intent(in)                   :: i_radar, j_radar, elv  ! radar station position and elv
  real                                :: r                      ! radar data point distance

  real, dimension(ix+1, jx, kx), intent(in) :: ugrid
  real, dimension(ix, jx+1, kx), intent(in) :: vgrid
  real, dimension(ix, jx, kx+1), intent(in) :: wgrid
  real, dimension(ix, jx      ), intent(in) :: xlong

  real, intent(out)                   :: rv

  integer                             :: n, l, i1, j1, k1, i, j, k  
  real                                :: ire, jre, kre   ! residue for ii, jj, kk
  real, dimension(3,  3)              :: wind1    ! for 3 levels, 3 for u, v, w
  real, dimension(3)                  :: wind     ! for obser point u, v, w
  real, dimension(2)                  :: windt    ! 2: ut, vt
  real                                :: long     ! 
  real                                :: x, y, z


!------------------------------------------------------------------------------
! calculate the rv

     i1  = int( ii )
     j1  = int( jj )
     k1  = int( kk )
     ire = ii - i1
     jre = jj - j1
     kre = kk - k1

!------------------------------------------------------------------------------
!  horizontal interpolate u,v,w 
     do k = k1, k1+1

!       for 2 levels adjacent obser point u
        if ( ire == 0.500 ) then
           wind1(k-k1+1,1) = ugrid(i1+1,j1,k)*(1-jre)+ugrid(i1+1,j1+1,k)*jre
        else if ( ire > 0.500 ) then
           wind1(k-k1+1,1) = ( ugrid(i1+1,j1,k)*(1.5-ire) + ugrid(i1+2,j1,k)*(ire-0.5) )*(1-jre)+  &
                             ( ugrid(i1+1,j1+1,k)*(1.5-ire) + ugrid(i1+2,j1+1,k)*(ire-0.5) )*jre
        else if ( ire < 0.500 ) then
           wind1(k-k1+1,1) = ( ugrid(i1,j1,k)*(0.5-ire) + ugrid(i1+1,j1,k)*(0.5+ire) )*(1-jre)+  &
                             ( ugrid(i1,j1+1,k)*(0.5-ire) + ugrid(i1+1,j1+1,k)*(0.5+ire) )*jre
        endif

!       for 2 levels adjacent obser point v
        if ( jre == 0.500 ) then
           wind1(k-k1+1,2) = vgrid(i1,j1+1,k)*(1-ire)+vgrid(i1+1,j1+1,k)*ire
        else if ( jre > 0.500 ) then
           wind1(k-k1+1,2) = ( vgrid(i1,j1+1,k)*(1.5-jre) + vgrid(i1,j1+2,k)*(jre-0.5) )*(1-ire)+  &
                             ( vgrid(i1+1,j1+1,k)*(1.5-jre) + vgrid(i1+1,j1+2,k)*(jre-0.5) )*ire
        else if ( jre < 0.500 ) then
           wind1(k-k1+1,2) = ( vgrid(i1,j1,k)*(0.5-jre) + vgrid(i1,j1+1,k)*(0.5+jre) )*(1-ire)+  &
                             ( vgrid(i1+1,j1,k)*(0.5-jre) + vgrid(i1+1,j1+1,k)*(0.5+jre) )*ire
        endif

    enddo

!       for 3 levels adjacent obser point w
    do k = k1, k1+2
       wind1(k-k1+1,3) = ( wgrid(i1,j1,k)*(1-ire) + wgrid(i1+1,j1,k)*ire )*(1-jre) + &
                         ( wgrid(i1,j1+1,k)*(1-ire) + wgrid(i1+1,j1+1,k)*ire )*jre
    enddo

! vertical interpolate u,v,w
    do i = 1, 2     ! for u,v
       wind(i) = wind1(1,i)*(1-kre) + wind1(2,i)*kre
    enddo  
    if ( kre == 0.500 ) then
       wind(3) = wind1(k1,3)
    else if ( kre < 0.500 ) then
       wind(3) = wind1(1,3)*(0.5-kre) + wind1(2,3)*(0.5+kre)
    else if ( kre > 0.500 ) then
       wind(3) = wind1(2,3)*(1.5-kre) + wind1(3,3)*(kre-0.5)
    endif
               
!------------------------------------------------------------------------------
! calculate obser point longtitude
!    long = ( xlong(i1,j1)*(1-jre) + xlong(i1,j1+1)*jre ) * (1-ire) +          &
!           ( xlong(i1+1,j1)*(1-jre) + xlong(i1+1,j1+1)*jre ) * ire
  
! convert wind from grid north to true north 
!    call gridwind_to_truewind( long, proj, wind(1), wind(2),    &
!         windt(1), windt(2) )
!------------------------------------------------------------------------------
! calculate rv: vr=(x/r)u+(y/r)v+(z/r)w
     x = ( ii - i_radar ) * proj%dx / 1000.
     y = ( jj - j_radar ) * proj%dx / 1000.
     z =  hh - elv / 1000.
     r = sqrt( x**2. + y**2. + z**2. )
!    if ( use_simulated ) then
       rv = (x/r)*wind(1) + (y/r)*wind(2) + (z/r)*wind(3) 
!    else
!       rv = (x/r)*windt(1) + (y/r)*windt(2) + (z/r)*wind(3) 
!    endif
!    write(*,'(8f10.4)')windt(1),windt(2),wind(3),x,y,z,r,rv
 
!------------------------------------------------------------------------------

  end subroutine calculate_rv 

!=======================================================================================

!------------------------------------------------------------------------------
end module radar
