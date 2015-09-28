Program produce_simulated_dropsonde

  use constants
  use mapinfo_define
  use netcdf
  use map_utils
  use wrf_tools

!----------------------------------------------------------------
  implicit none

  type(proj_info)                       :: proj
  character (len=80)                    :: times, filename
  character (len=10)                    :: truth_file
  integer                               :: ix, jx, kx

  integer                               :: pattern, num_obs
  real                                  :: flight_dis, flight_spd, drop_interval, fight_hgt
  real                                  :: flight_dis1, flight_spd1, drop_interval1, fight_hgt1    !the 1st stage flight state, for p3,p4,p5,p6
  real, allocatable, dimension(:)       :: obs_x, obs_y, obs_i, obs_j
  real, allocatable, dimension(:)       :: obs_lat, obs_lon
  integer, allocatable, dimension(:)    :: obs_lev
  character(len=19),allocatable,dimension(:) :: obs_time
  real, allocatable, dimension(:,:)     :: obs_p, obs_z, obs_t, obs_td, obs_spd, obs_dir, obs_ph, obs_rh, obs_theta, obs_u, obs_v, obs_qv
 
  real, allocatable, dimension(:,:)     :: xlong, xlat
  real, allocatable, dimension(:)       :: znu, znw
  real, allocatable, dimension(:,:,:)   :: u, v
  real, allocatable, dimension(:,:,:)   :: ph, phb, p, pb, qv, qc, qr, theta

  real,dimension(4)                     :: center    !lat, lon, slp, wsp
  real                                  :: hio, hjo, dx, dy, time_interval

  integer                               :: fid, rcode, i, j, k, n, ii, jj, kk
  real                                  :: dis, di, dj, missing_r, dis0, ddis, alpha, utrue, vtrue
  CHARACTER (LEN = 120)                 :: fmt_info, fmt_srfc, fmt_each

  missing_r = -888888.
!----------------------------------------------------------------
  pattern = 3    ! 1, 2, 3, 4, 5, 6, 71, 72
!  pattern = 71
  fight_hgt = 6000.        !m
! Observations pattern 1: dropsondes every 5 minutes from an aircraft starting 27/12 
!                         traveling at a speed of 110m/s straight from East to West 
!                         for 200km centered on the hurricane eye, then goes 45 degrees 
!                         to the southeast for 100*sqrt(2)km, then goes straight from south 
!                         to the north for another 200km.
! Observations pattern 2: dropsondes every 8 minutes from an aircraft starting 27/15 
!                         traveling at a speed of 110m/s straight from East to West for 
!                         200km centered on the hurricane eye, then turns right 60 degrees 
!                         to the southeast for 100km, which turn right 60degree for another 
!                         200km and then turn right for another 60 degree for 100km and turns 
!                         right for 60degree for 200km. 
  if ( pattern == 1 ) then
     flight_dis = (400.+100.*sqrt(2.))*1000   !km
     flight_spd = 110.        !m/s
     drop_interval = 5.*60.   !5 min = 300s
  else if ( pattern == 2 ) then
     flight_dis = 800.*1000   !800km
     flight_spd = 110.        !m/s
     drop_interval = 8.*60.   !8 min = 480s
  else if ( pattern == 3 ) then
     flight_dis = (400.+100.*sqrt(2.))*1000   !km
     flight_spd = 110.        !m/s
     drop_interval = 5.*60.   !5 min = 300s
  else if ( pattern == 4 ) then
     flight_dis = 800.*1000   !800km
     flight_spd = 110.        !m/s
     drop_interval = 8.*60.   !8 min = 480s
  else if ( pattern == 5 ) then
     flight_dis = 2.*(400.+100.*sqrt(2.))*1000   !km
     flight_spd = 110.        !m/s
     drop_interval = 5.*60.   !5 min = 300s
  else if ( pattern == 6 ) then
     flight_dis = 1600.*1000   !800km
     flight_spd = 110.        !m/s
     drop_interval = 8.*60.   !8 min = 480s
  else if ( pattern == 71 .or. pattern == 72 ) then
     flight_dis = (400.+100.*sqrt(2.))*1000   !km
     flight_spd = 110.        !m/s
     drop_interval = 5.*60.   !5 min = 300s
  endif

  num_obs = int(flight_dis/(flight_spd*drop_interval))
  allocate( obs_x( num_obs ))
  allocate( obs_y( num_obs ))
  allocate( obs_i( num_obs ))
  allocate( obs_j( num_obs ))

!.... calculate the horizontal position in km centered by hurricane eye
  do n = 1, num_obs 
     dis = real(n-1)*flight_spd*drop_interval/1000.

     if ( pattern == 1 ) then
        if ( dis <= 200. ) then
           obs_x(n) = 100. - dis
           obs_y(n) = 0.
        else if ( dis > 200. .and. dis <= (200.+100.*2.0**0.5) ) then
           obs_x(n) =  (dis-200.)*cos(45./180.*3.1415926)-100.
           obs_y(n) = -(dis-200.)*sin(45./180.*3.1415926) 
        else if ( dis > (200.+100.*2.0**0.5) .and. dis <= (400.+100.*2.0**0.5) ) then
           obs_x(n) = 0.
           obs_y(n) = -100. + (dis-200.-100.*2.0**0.5)
        endif
     else if ( pattern == 2 ) then
        if ( dis <= 200. ) then
           obs_x(n) = 100. - dis
           obs_y(n) = 0.
        else if ( dis > 200. .and. dis <= 300. ) then
           obs_x(n) =  (dis-200.)*cos(60./180.*3.1415926)-100.
           obs_y(n) = -(dis-200.)*sin(60./180.*3.1415926)
        else if ( dis > 300. .and. dis <= 500. ) then
           obs_x(n) = -100.*cos(60./180.*3.1415926)+(dis-300.)*sin(30./180.*3.1415926) 
           obs_y(n) = -100.*sin(60./180.*3.1415926)+(dis-300.)*cos(30./180.*3.1415926) 
        else if ( dis > 500. .and. dis <= 600. ) then
           obs_x(n) = -(dis-500.)+100.*sin(30./180.*3.1415926) 
           obs_y(n) = 100.*cos(30./180.*3.1415926) 
        else if ( dis > 600. .and. dis <= 800. ) then
           obs_x(n) = (dis-700.)*sin(30./180.*3.1415926)
           obs_y(n) = (700.-dis)*sin(60./180.*3.1415926) 
        endif
     else if ( pattern == 3 ) then   ! rotate 45degree and continue for p1
        if ( dis <= 200. ) then
           ddis = dis - 100.
           obs_x(n) =  ddis*cos(45./180.*3.1415926)
           obs_y(n) = -ddis*sin(45./180.*3.1415926)
        else if ( dis > 200. .and. dis <= (200.+100.*sqrt(2.)) ) then
           ddis = dis - (200.+100.*sqrt(2.)/2.)
           obs_x(n) = 100.*sqrt(2.)/2.
           obs_y(n) = ddis*sin(45./180.*3.1415926) 
        else if ( dis > (200.+100.*sqrt(2.)) .and. dis <= (400.+100.*2.0**0.5) ) then
           ddis = 200.+100.*sqrt(2.)+100. - dis
           obs_x(n) = ddis*cos(45./180.*3.1415926)
           obs_y(n) = ddis*sin(45./180.*3.1415926)
        endif 
     else if ( pattern == 4 ) then
        if ( dis <= 200. ) then
           ddis = dis - 100.
           obs_x(n) = ddis*sin(15./180.*3.1415926)
           obs_y(n) = ddis*cos(15./180.*3.1415926)
        else if ( dis > 200. .and. dis <= 300. ) then
           ddis = dis - 200.
           obs_x(n) = 100.*sin(15./180.*3.1415926) - ddis*cos(15./180.*3.1415926)
           obs_y(n) = 100.*cos(15./180.*3.1415926) - ddis*sin(15./180.*3.1415926)
        else if ( dis > 300. .and. dis <= 500. ) then
           ddis = dis - 400.
           obs_x(n) = ddis*cos(45./180.*3.1415926)
           obs_y(n) =-ddis*sin(45./180.*3.1415926)
        else if ( dis > 500. .and. dis <= 600. ) then
           ddis = dis - 500.
           obs_x(n) = 100.*cos(45./180.*3.1415926) + ddis*cos(75./180.*3.1415926)
           obs_y(n) =-100.*sin(45./180.*3.1415926) + ddis*sin(75./180.*3.1415926)
        else if ( dis > 600. .and. dis <= 800. ) then
           ddis = 700. - dis
           obs_x(n) = ddis*cos(15./180.*3.1415926)
           obs_y(n) = ddis*sin(15./180.*3.1415926)
        endif
     else if ( pattern == 5 ) then
        if ( dis <= 400. ) then
           obs_x(n) = 200. - dis
           obs_y(n) = 0.
        else if ( dis > 400. .and. dis <= (400.+200.*2.0**0.5) ) then
           obs_x(n) =  (dis-400.)*cos(45./180.*3.1415926)-200.
           obs_y(n) = -(dis-400.)*sin(45./180.*3.1415926) 
        else if ( dis > (400.+200.*2.0**0.5) .and. dis <= (800.+200.*2.0**0.5) ) then
           obs_x(n) = 0.
           obs_y(n) = -200. + (dis-400.-200.*2.0**0.5)
        endif
     else if ( pattern == 6 ) then
        if ( dis <= 400. ) then
           obs_x(n) = 200. - dis
           obs_y(n) = 0.
        else if ( dis > 400. .and. dis <= 600. ) then
           obs_x(n) =  (dis-400.)*cos(60./180.*3.1415926)-200.
           obs_y(n) = -(dis-400.)*sin(60./180.*3.1415926)
        else if ( dis > 600. .and. dis <= 1000. ) then
           obs_x(n) = -200.*cos(60./180.*3.1415926)+(dis-600.)*sin(30./180.*3.1415926)
           obs_y(n) = -200.*sin(60./180.*3.1415926)+(dis-600.)*cos(30./180.*3.1415926)
        else if ( dis > 1000. .and. dis <= 1200. ) then
           obs_x(n) = -(dis-1000.)+200.*sin(30./180.*3.1415926)
           obs_y(n) = 200.*cos(30./180.*3.1415926)
        else if ( dis > 1200. .and. dis <= 1600. ) then
           obs_x(n) = (dis-1400.)*sin(30./180.*3.1415926)
           obs_y(n) = (1400.-dis)*sin(60./180.*3.1415926)
        endif
     else if ( pattern == 71 .or. pattern == 72) then
        if (pattern == 72 ) dis = dis + 400.+100.*2.0**0.5
        if ( dis <= 400. ) then
           obs_x(n) = 200. - dis
           obs_y(n) = 0.
        else if ( dis > 400. .and. dis <= (400.+200.*2.0**0.5) ) then
           obs_x(n) =  (dis-400.)*cos(45./180.*3.1415926)-200.
           obs_y(n) = -(dis-400.)*sin(45./180.*3.1415926)
        else if ( dis > (400.+200.*2.0**0.5) .and. dis <= (800.+200.*2.0**0.5) ) then
           obs_x(n) = 0.
           obs_y(n) = -200. + (dis-400.-200.*2.0**0.5)
        endif
        dis0 = sqrt(obs_x(n)**2 + obs_y(n)**2)
        if ( dis0 <= 100. )fight_hgt = 3000.
        if ( dis0 >  100. )fight_hgt = 6000.
     endif

  enddo 

!----------------------------------------------------------------
  truth_file = 'fort.80010'
  call get_wrf_info( truth_file, ix, jx, kx, times, proj)

  allocate( xlong ( ix, jx ) )
  allocate( xlat  ( ix, jx ) )
  allocate( znu   ( kx ) )
  allocate( znw   ( kx+1 ) )
  allocate( u     ( ix+1, jx, kx  ) )
  allocate( v     ( ix, jx+1, kx  ) )
  allocate( theta ( ix  , jx, kx  ) )
  allocate( ph    ( ix  , jx, kx+1) )
  allocate( phb   ( ix  , jx, kx+1) )
  allocate( p     ( ix  , jx, kx  ) )
  allocate( pb    ( ix  , jx, kx  ) )
  allocate( qv    ( ix  , jx, kx  ) )
  allocate( qc    ( ix  , jx, kx  ) )
  allocate( qr    ( ix  , jx, kx  ) )

  call get_variable2d( truth_file, 'XLONG     ', ix, jx, 1, xlong )
  call get_variable2d( truth_file, 'XLAT      ', ix, jx, 1, xlat  )
  call get_variable1d( truth_file, 'ZNU       ', kx, 1, znu )
  call get_variable1d( truth_file, 'ZNW       ', kx+1, 1, znw )
  call open_file ( truth_file, nf_nowrite, fid )
  rcode = nf_get_att_real(fid, nf_global, 'DX', dx )
  rcode = nf_get_att_real(fid, nf_global, 'DY', dy )
  call close_file ( fid )

!----------------------------------------------------------------
! calculate hurricane eye position
  call hurricane_center_wind ( truth_file, proj, ix, jx, kx, xlong, xlat, znu, znw, &
                                    ix/2, jx/2, center )
  hio = center(1)
  hjo = center(2)

!----------------------------------------------------------------
! calculate obs position in wrf-domain grid
  write(*,*)hio, hjo, dx, dy
  do n = 1, num_obs
     obs_i(n) = obs_x(n)*1000./dx + hio
     obs_j(n) = obs_y(n)*1000./dy + hjo
  enddo

  allocate( obs_lat( num_obs ))
  allocate( obs_lon( num_obs ))
  allocate( obs_time( num_obs ))
  allocate( obs_lev( num_obs ))
  allocate( obs_p( num_obs,kx ))
  allocate( obs_z( num_obs,kx ))
  allocate( obs_t( num_obs,kx ))
  allocate( obs_td( num_obs,kx ))
  allocate( obs_rh( num_obs,kx ))
  allocate( obs_spd( num_obs,kx ))
  allocate( obs_dir( num_obs,kx ))

  allocate( obs_theta( num_obs,kx ))
  allocate( obs_u( num_obs,kx ))
  allocate( obs_v( num_obs,kx ))
  allocate( obs_qv( num_obs,kx ))

  allocate( obs_ph( num_obs,kx+1 ))
!----------------------------------------------------------------
! get background
  call get_variable3d( truth_file, 'U         ', ix+1, jx, kx, 1, u )
  call get_variable3d( truth_file, 'V         ', ix, jx+1, kx, 1, v )
  call get_variable3d( truth_file, 'T         ', ix, jx, kx, 1, theta )
  call get_variable3d( truth_file, 'PH        ', ix, jx, kx+1, 1, ph )
  call get_variable3d( truth_file, 'PHB       ', ix, jx, kx+1, 1, phb )
  call get_variable3d( truth_file, 'P         ', ix, jx, kx, 1, p )
  call get_variable3d( truth_file, 'PB        ', ix, jx, kx, 1, pb )
  call get_variable3d( truth_file, 'QVAPOR    ', ix, jx, kx, 1, qv )
  call get_variable3d( truth_file, 'QCLOUD    ', ix, jx, kx, 1, qc )
  call get_variable3d( truth_file, 'QRAIN     ', ix, jx, kx, 1, qr )
  ph = (ph + phb)/9.81
  p  = p + pb
  qv = (qv + qc + qr)*1000.

!----------------------------------------------------------------
! calculate obs
  do n = 1, num_obs

     call ij_to_latlon(proj, obs_i(n), obs_j(n), obs_lat(n), obs_lon(n))

     time_interval = real(n-1)*drop_interval
     call cal_time(times(1:19), int(time_interval), obs_time(n))

     ii = int(obs_i(n))
     jj = int(obs_j(n))
     di = obs_i(n) - real(ii)
     dj = obs_j(n) - real(jj)

     do k = 1, kx+1
        obs_ph(n,k)=horiz_interp(ph(ii:ii+1,jj:jj+1,k), di, dj)
     enddo
     do k = 1, kx
        obs_z(n,k)=(obs_ph(n,k)+obs_ph(n,k+1))/2.
     enddo 

     obs_lev(n) = 0
     do k = 1, kx
        if ( obs_z(n,k) <= fight_hgt ) then
           obs_lev(n) = obs_lev(n) + 1
           obs_p(n,k)=horiz_interp( p(ii:ii+1,jj:jj+1,k), di, dj)        

           obs_u(n,k)=horiz_interp( u(ii:ii+1,jj:jj+1,k), di, dj)        
           obs_v(n,k)=horiz_interp( v(ii:ii+1,jj:jj+1,k), di, dj)        
           call gridwind_to_truewind(obs_lon(n), proj, obs_u(n,k), obs_v(n,k), utrue, vtrue)           
           call uv_to_dirspd(utrue, vtrue, obs_dir(n,k), obs_spd(n,k))

           obs_theta(n,k)=horiz_interp( theta(ii:ii+1,jj:jj+1,k), di, dj)        
           obs_t(n,k)=theta_to_temp(obs_theta(n,k)+to, obs_p(n,k)) 

           obs_qv(n,k)=horiz_interp( qv(ii:ii+1,jj:jj+1,k), di, dj)        
!           obs_td(n,k) = mixrat_to_tdew(obs_qv(n,k), obs_p(n,k)) + T_frez 
           obs_td(n,k) = missing_r
!           obs_rh(n,k) = rel_humidity(obs_qv(n,k), obs_t(n,k), obs_p(n,k))
           obs_rh(n,k) = missing_r
           if ( obs_rh(n,k) .gt.100.) obs_rh(n,k)=100.
        endif
     enddo
     write(*,*)obs_x(n), obs_y(n), obs_i(n), obs_j(n)
  enddo 

!----------------------------------------------------------------
! output to 3d-VAR obs format  

! 1. OPEN FILE FOR VALID OBSERVATIONS OUTPUT
! ==========================================
! 1.1 Name of output file
!     -------------------
  filename = 'ADPUPA_p1h6d2_'//times(1:4)//times(6:7)//times(9:10)//times(12:13)
  write(filename(9:9),'(i1.1)')pattern
  write(filename(11:11),'(i1.1)')int(fight_hgt/1000.)
  if(pattern == 71 .or. pattern == 72) then
     write(filename(9:9),'(i1.1)')pattern/10
     write(filename(11:11),'(i1.1)')pattern-pattern/10*10
  endif
! 1.2 OPEN FILE
  open(99, file=filename, FORM   = 'FORMATTED', ACCESS = 'SEQUENTIAL', STATUS = 'REPLACE')

! 2.  WRITE OBSERVATIONS
! ======================

! 2.1  Loop over stations
!      ------------------
  stations: DO n = 1, num_obs
     
! 2.4 Write station info
!     ------------------
  write(99,'(a8,1x,2f7.2,1x,f7.3,1x,2f8.1,1x,f7.1,i8)')  &
            '99999   ',obs_lon(n),obs_lat(n),1.,-88888.,120.,11.0,obs_lev(n)
  write(99,'(a72)')'       P       Q     QOE       T     TOE       Z       U       V     WOE'
  do k = 1, obs_lev(n)
     write(99,'(9f8.1)')obs_p(n,k)/100.,obs_qv(n,k)*1000.,2.0,obs_t(n,k)-273.15,1.0,obs_z(n,k),obs_u(n,k),obs_v(n,k),2.0
  enddo

  enddo stations
  close(99)

  end
 



!==============================================================================
  subroutine cal_time(time_in, second1, time_out)

  implicit none
  character(len=19), intent(in)   :: time_in
  integer, intent(in)             :: second1
  character(len=19), intent(out)  :: time_out

  integer                         :: year, month, day, hour, minute, second
  integer                         :: i, j, k

  read(time_in,'(i4,5(1x,i2))')year, month, day, hour, minute, second

  second = second + second1
  if ( second .lt. 0 ) then
     minute = minute + second/60 - 1
     second = -(second/60 - 1)*60 + second
  else
     minute = minute + second/60
     second = second - second/60*60
  endif
  if ( minute .lt. 0 ) then
     hour = hour + minute/60 - 1
     minute = -(minute/60 - 1)*60 + minute
  else
     hour = hour + minute/60
     minute = minute - minute/60*60
  endif
  if ( hour .lt. 0 ) then
     day = day + hour/24 - 1
     hour = -(hour/24 -1)*24 + hour
  else
     day = day + hour/24
     hour = hour - hour/24*24
  endif

  do_day_loop : do
     if ( day .le. 0 ) then
        month = month - 1
        if ( month.eq.2 .and. mod(year,4).eq.0 ) day = 29 + day   !day=0-->day=29; day=-1-->day=28
        if ( month.eq.2 .and. mod(year,4).ne.0 ) day = 28 + day
        if ( month.eq.1 .or. month.eq.3 .or. month.eq.5 .or. month.eq.7 .or. &
             month.eq.8 .or. month.eq.10 .or. month.eq.12 )day = 31 + day
        if ( month.eq.4 .or. month.eq.6 .or. month.eq.9 .or. month.eq.11 )day = 30 + day
     else
        if ( month.eq.2 )then
           if (mod(year,4).eq.0 .and. day.gt.29 )then
              month = month + 1
              day = day - 29
           else if (mod(year,4).ne..0 .and. day.gt.8 )then
              month = month + 1
              day = day - 28
           else
              exit do_day_loop
           endif
        else if ( month.eq.1 .or. month.eq.3 .or. month.eq.5 .or. month.eq.7 .or. &
                  month.eq.8 .or. month.eq.10 .or. month.eq.12 ) then
           if ( day .gt. 31 )then
              month = month + 1
              day = day - 31
           else
              exit do_day_loop
           endif
        else
           if ( day .gt. 30 )then
              month = month + 1
              day = day - 30
           else
              exit do_day_loop
           endif
        endif
     endif
  enddo do_day_loop

  do_month_loop : do
     if ( month.le.0 ) then
        year = year - 1
        month = month + 12
     else if ( month.gt.12 ) then
        year = year + 1
        month = month - 12
     else
        exit do_month_loop
     endif
  enddo do_month_loop

  time_out = time_in
  write(time_out,'(i4,5(1x,i2))')year, month, day, hour, minute, second

  end subroutine cal_time
!==============================================================================
