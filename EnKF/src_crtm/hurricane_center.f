!=======================================================================================
   subroutine hurricane_center_shift(filename,ix,jx,kx,xlong,xlat,znu,znw,proj,times)
!  hurricane center shift (real data relative to real hurricane center, 
!  simulated data relative to truth's hurricane center)

   use constants
   use namelist_define
   use obs_define
   use mapinfo_define
   use map_utils
   use netcdf
   implicit none

   character (len=10), intent(in)          :: filename
   type (proj_info)                        :: proj 
   character (len=80), intent(in)          :: times
   integer, intent(in) :: ix, jx, kx
   real, dimension(kx), intent(in)         :: znu
   real, dimension(kx+1), intent(in)       :: znw
   real, dimension(ix, jx), intent(in)     :: xlong, xlat

   real,dimension(3)      :: center, center_xb
   real                   :: center_io, center_jo, xb
   integer                :: iob

!.... obs relative to hurricane center, so move obs' position around truth hurricane center
      if( print_detail > 3)write(*,*)'1st obs position', obs%position(1,1),obs%position(1,2)


!......... get real hurricane center position and change to model position
           call rd_hurricane_track ( times, 'hurricane_best_track', center )
           call latlon_to_ij ( proj, center(1), center(2), center_io, center_jo )
           if( print_detail > 100)write(*,*)' hurricane center in model: ', center_io, center_jo

!......... get truth's hurricane center position
           if( print_detail > 100)write(*,*)'call hurricane_center_assimilation'
           call hurricane_center_assimilation (filename,ix,jx,kx,int(center_io),int(center_jo),center_xb, &
                                               znu,znw,xlong,xlat,proj)

!......... shift simulated obs
           if( print_detail > 1)write(*,*)'hurricane obs position ', center(1), center(2), center_io, center_jo
           center(1) = center_xb(1) - center_io
           center(2) = center_xb(2) - center_jo
           if( print_detail > 1)write(*,*)'truth hurricane position ', center_xb(1:2)
           do iob = 1, obs%num
              if( print_detail > 3)write(*,*)'real  obs position ', obs%position(iob,1),obs%position(iob,2)
              if (obs%type(iob) == 'longtitude' )then
                 obs%dat(iob)      = center_xb(1)
                 obs%position(iob,1)  = center_xb(1)
                 obs%position(iob,2)  = center_xb(2)
              else if (obs%type(iob) == 'latitude  ' )then
                 obs%dat(iob)      = center_xb(2)
                 obs%position(iob,1)  = center_xb(1)
                 obs%position(iob,2)  = center_xb(2)
              else if (obs%type(iob) == 'min_slp   ' ) then
                 obs%dat(iob)      = center_xb(3)
                 obs%position(iob,1)  = center_xb(1)
                 obs%position(iob,2)  = center_xb(2)
              else
                 if ( use_simulated ) then
                    obs%position(iob,1) = obs%position(iob,1) + center(1)
                    obs%position(iob,2) = obs%position(iob,2) + center(2)
                    obs%sta(iob,1)      = obs%sta(iob,1)      + center(1)
                    obs%sta(iob,2)      = obs%sta(iob,2)      + center(2)
                 endif
              endif
              if( print_detail > 3)write(*,*)'truth obs position ', obs%position(iob,1),obs%position(iob,2)
           enddo

   return

   end subroutine hurricane_center_shift
!==============================================================================
subroutine hurricane_center_assimilation(inputfile,ix,jx,kx,io,jo,center,znu,znw,xlong,xlat,proj)
  use constants
  use namelist_define
  use netcdf
  use wrf_tools
  implicit none
  character(len=10), intent(in) :: inputfile
  type(proj_info) :: proj
  integer, intent(in) :: ix,jx,kx,io,jo
  real, dimension(kx), intent(in)         :: znu
  real, dimension(kx+1), intent(in)       :: znw
  real, dimension(ix, jx), intent(in)     :: xlong, xlat 
  real, intent(out), dimension(3) :: center 
  real, dimension(ix, jx)                 :: slp, slpsmth
  real, dimension(ix+1, jx, kx  )         :: u
  real, dimension(ix, jx+1, kx  )         :: v
  real, dimension(ix, jx, kx  )           :: t
  real, dimension(ix, jx, kx+1)           :: ph
  real, dimension(ix, jx, kx+1)           :: phb
  real, dimension(ix, jx, kx  )           :: p
  real, dimension(ix, jx, kx  )           :: z
  real, dimension(ix, jx, kx  )           :: pb
  real, dimension(ix, jx, kx  )           :: qv, qc, qr
  integer :: m, i, j, k, i1, j1, i2, j2, search_radiu
  real    :: dvdx, dudy, vor, vor1, wsp, wsp1

  call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
  call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
  call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
  call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
  call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
  call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
  call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
  call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

  ph = (ph + phb)/g
  p  = p + pb

!.... interp PH to unstaggered eta level
   do j = 1, jx
   do i = 1, ix
      z(i,j,1:kx) = (ph(i,j,1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i,j,2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
      do k=1,kx
      if( p(i,j,k)/100. >= 1200. ) then
          write(*,'(a,3i4,2f8.0)')'P error at:', i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
      endif
      enddo
   enddo
   enddo

!.... calculate slp
   call compute_seaprs(ix, jx, kx, z, t, p, qv, slp, 0)

   slpsmth = slp
   call smooth(slpsmth, ix, jx, 305)

!.... find the minumum slp
   center(3) = 200000.    !pa
   search_radiu = int(700./proj%dx*1000.)
   i1 = max(io-search_radiu, 5)
   i2 = min(io+search_radiu, ix-5)
   j1 = max(jo-search_radiu, 5)
   j2 = min(jo+search_radiu, jx-5)
   do j = j1, j2
   do i = i1, i2
         if( slpsmth(i,j) < center(3) )then
             center(3) = slpsmth(i,j)
             center(1) = i*1.
             center(2) = j*1.
         endif
   enddo
   enddo
   center(3) = slp(int(center(1)), int(center(2)))

!.... when minslp > 995mb, then find the max curvature of wind;
!.... that is to calculate the vor for (u,v)/(u**2+v**2)**0.5
   if (center(3) >= 995. ) then
       vor1 = 0.
       search_radiu = int(600./proj%dx*1000.)
       i1 = max(io-search_radiu, 5)
       i2 = min(io+search_radiu, ix-5)
       j1 = max(jo-search_radiu, 5)
       j2 = min(jo+search_radiu, jx-5)
       do j = j1, j2
       do i = i1, i2
          dvdx = 0.5*( v(i+1,j+1,8)/sqrt(u(i+1,j+1,8)*u(i+1,j+1,8)+v(i+1,j+1,8)*v(i+1,j+1,8))  &
                      -v(i  ,j+1,8)/sqrt(u(i  ,j+1,8)*u(i  ,j+1,8)+v(i  ,j+1,8)*v(i  ,j+1,8))  &
                      +v(i+1,j  ,8)/sqrt(u(i+1,j  ,8)*u(i+1,j  ,8)+v(i+1,j  ,8)*v(i+1,j  ,8))  &
                      -v(i  ,j  ,8)/sqrt(u(i  ,j  ,8)*u(i  ,j  ,8)+v(i  ,j  ,8)*v(i  ,j  ,8))  &
                     )/proj%dx
          dudy = 0.5*( u(i+1,j+1,8)/sqrt(u(i+1,j+1,8)*u(i+1,j+1,8)+v(i+1,j+1,8)*v(i+1,j+1,8))  &
                      -u(i+1,j  ,8)/sqrt(u(i+1,j  ,8)*u(i+1,j  ,8)+v(i+1,j  ,8)*v(i+1,j  ,8))  &
                      +u(i  ,j+1,8)/sqrt(u(i  ,j+1,8)*u(i  ,j+1,8)+v(i  ,j+1,8)*v(i  ,j+1,8))  &
                      -u(i  ,j  ,8)/sqrt(u(i  ,j  ,8)*u(i  ,j  ,8)+v(i  ,j  ,8)*v(i  ,j  ,8))  &
                     )/proj%dx
          vor  = (dvdx-dudy)*(1030.-slp(i,j))**4
          if( vor .gt. vor1 )then
              center(1) = i*1.
              center(2) = j*1.
              vor1 = vor
          endif
      enddo
      enddo
      center(3) = slp(int(center(1)), int(center(2)))
   endif

!.... find the minumum slp
!   search_radiu = int(150./proj%dx*1000.)
!   i1 = max(int(center(1))-search_radiu, 5)
!   i2 = min(int(center(1))+search_radiu, ix-5)
!   j1 = max(int(center(2))-search_radiu, 5)
!   j2 = min(int(center(2))+search_radiu, jx-5)
!   do j = j1, j2
!   do i = i1, i2
!      if( slpsmth(i,j) < center(3) )then
!          center(3) = slpsmth(i,j)
!          center(1) = i*1.
!          center(2) = j*1.
!      endif
!   enddo
!   enddo
!   center(3) = slp(int(center(1)), int(center(2)))
   center(3) = slp(io, jo)

  end subroutine hurricane_center_assimilation

!==============================================================================
! subroutine hurricane_center_wind ( wrf_file, proj, ix, jx, kx, xlong, xlat, znu, znw, &
!                                    io, jo, center )
!
!---------------------
! slp_center subroutine calculates sea level pressure and find the hurricane center
!---------------------
!
!  use constants
!  use namelist_define
!  use netcdf
!  use wrf_tools
!
!  implicit none
!   
!  type(proj_info)                         :: proj
!  integer, intent(in)                     :: ix, jx, kx       ! 1st guest dimension
!  real, dimension(kx), intent(in)         :: znu
!  real, dimension(kx+1), intent(in)       :: znw
!  integer, intent(in)                     :: io, jo    ! center position of observation
!  real, dimension(ix, jx), intent(in)     :: xlong, xlat
!  character(len=10), intent(in)           :: wrf_file
!
!  real, intent(out), dimension(4)         :: center   !1-lon, 2-lat, 3-slp, 4-wsp
!
!  real, dimension(ix, jx)                 :: slp, slpsmth
!  real, dimension(ix+1, jx, kx  )         :: u
!  real, dimension(ix, jx+1, kx  )         :: v
!  real, dimension(ix, jx, kx  )           :: t
!  real, dimension(ix, jx, kx+1)           :: ph
!  real, dimension(ix, jx, kx+1)           :: phb
!  real, dimension(ix, jx, kx  )           :: p
!  real, dimension(ix, jx, kx  )           :: pb
!  real, dimension(ix, jx, kx  )           :: qv, qc, qr
!  real, dimension(ix, jx )                :: u10, v10, u101, v101
!  real, dimension(ix, jx )                :: height300, height400, height500, height
!  integer                                 :: itot, m, i, j, k, i1, j1, i2, j2, search_radiu
!  real                                    :: dx, heightmin, dvdx, dudy, vor, vor1, wsp, wsp1
!  integer                                 :: iloc, jloc
!  integer                                 :: isearch_r, isearchs, isearche, jsearchs, jsearche
!
!  integer                                 :: rcode, fid, varid
!    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! 1st, get ph, phb, t, qv, p, pb from xa(ntot)
!
!   call get_variable3d( wrf_file, 'U         ', ix+1, jx, kx, 1, u )
!   call get_variable3d( wrf_file, 'V         ', ix, jx+1, kx, 1, v )
!   call get_variable3d( wrf_file, 'PH        ', ix, jx, kx+1, 1, ph )
!   call get_variable3d( wrf_file, 'PHB       ', ix, jx, kx+1, 1, phb )
!   call get_variable3d( wrf_file, 'P         ', ix, jx, kx, 1, p )
!   call get_variable3d( wrf_file, 'PB        ', ix, jx, kx, 1, pb )
!   ph = (ph + phb)/g    
!   p  = p + pb
!....... interp PH to unstaggered eta level
!   do j = 1, jx
!   do i = 1, ix
!   do k = 1, kx
!      ph(i,j,k) = (ph(i,j,k)*(znw(k+1)-znu(k))+ph(i,j,k+1)*(znu(k)-znw(k)))/(znw(k+1)-znw(k))
!      if( p(i,j,k)/100. >= 1200. ) then
!          write(*,'(a,3i4,2f8.0)')'P error at:', i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
!      endif
!   enddo
!   enddo
!   enddo
!
!.......compute 500, 400 and 300mb height
!   do j = 1, jx
!   do i = 1, ix
!      height300(i,j)=interp_thta_pres(ph(i,j,1:kx), p(i,j,1:kx)/100., 300., kx) 
!      height400(i,j)=interp_thta_pres(ph(i,j,1:kx), p(i,j,1:kx)/100., 400., kx) 
!      height500(i,j)=interp_thta_pres(ph(i,j,1:kx), p(i,j,1:kx)/100., 500., kx) 
!      height(i,j)=height300(i,j)+height400(i,j)+height500(i,j)
!      height(i,j)=height500(i,j)
!   enddo
!   enddo
!
!....... calculate slp
!   call open_file( wrf_file, nf_nowrite, fid )
!   rcode = nf_inq_varid(     fid, 'PMSL      ', varid )
!   if ( rcode .ne. 0 ) then
!       call get_variable3d( wrf_file, 'T         ', ix, jx, kx, 1, t )
!       call get_variable3d( wrf_file, 'QVAPOR    ', ix, jx, kx, 1, qv )
!       call get_variable3d( wrf_file, 'QCLOUD    ', ix, jx, kx, 1, qc )
!       call get_variable3d( wrf_file, 'QRAIN     ', ix, jx, kx, 1, qr )
!       qv = qv + qc + qr
!       call compute_seaprs(ix, jx, kx, ph, t, p, qv, slp, 0)
!   else
!       call get_variable2d( wrf_file, 'PMSL      ', ix, jx, 1, slp)
!   endif
!
!   slpsmth = slp 
!
!.... smooth
!   call smooth(slpsmth, ix, jx, 305)
!   call smooth(height, ix, jx, 305)
!
!...  find the minumum height and the position (iloc, jloc)
!   rcode = nf_get_att_real(fid, nf_global, 'DX', dx )
!   dx=dx/1000.
!   isearch_r=int(750./dx)
!   if ( io .lt. 1 .or. jo .lt. 1 .or. io .ge. ix .or. jo .ge. jx ) then
!      isearchs = 5
!      isearche = ix-5
!      jsearchs = 5
!      jsearche = jx-5
!   else
!      isearchs = max(io - isearch_r, 5)
!      isearche = min(io + isearch_r, ix-5)
!      jsearchs = max(jo - isearch_r, 5)
!      jsearche = min(jo + isearch_r, jx-5)
!   endif
!
!   iloc=-99
!   jloc=-99
!   heightmin=99999999.
!   do j = jsearchs, jsearche
!   do i = isearchs, isearche
!      if ( height(i,j) < heightmin ) then
!           heightmin = height(i,j)
!           iloc = i
!           jloc = j
!      endif
!   enddo
!   enddo
!   write(*,*)io, jo, iloc, jloc, isearch_r, isearchs, isearche, jsearchs, jsearche
!   iloc=io
!   jloc=jo
!
!.... setup the search radis 
!   isearch_r=int(150./dx)
!   if ( iloc .lt. 1 .or. jloc .lt. 1 .or. iloc .ge. ix .or. jloc .ge. jx ) then
!      isearchs = 5
!      isearche = ix-5
!      jsearchs = 5
!      jsearche = jx-5
!   else
!      isearchs = max(iloc - isearch_r, 5)
!      isearche = min(iloc + isearch_r, ix-5)
!      jsearchs = max(jloc - isearch_r, 5)
!      jsearche = min(jloc + isearch_r, jx-5)
!   endif
!
!.... find the minumum slp
!   center(3) = 200000.    !pa
!   do j = jsearchs, jsearche
!   do i = isearchs, isearche
!      if( abs(xlong(i,j)-xlong(io,jo)).le.3. .and. abs(xlat(i,j)-xlat(io,jo)).le.3. ) then
!         if( slpsmth(i,j) < center(3) )then
!             center(3) = slpsmth(i,j)
!             center(1) = i*1.
!             center(2) = j*1.
!         endif
!      endif
!   enddo
!   enddo
!   center(3) = slp(int(center(1)), int(center(2)))
!
!.... when minslp > 995mb, then find the max curvature of wind;
!.... that is to calculate the vor for ((u,v)/(u**2+v**2)**0.5)*(1030.-slp(i,j))**4
!   if (center(3) >= 995. ) then
!       vor1 = 0. 
!       do j = jsearchs, jsearche
!       do i = isearchs, isearche
!          dvdx = 0.5*( v(i+1,j+1,8)/sqrt(u(i+1,j+1,8)*u(i+1,j+1,8)+v(i+1,j+1,8)*v(i+1,j+1,8))  &
!                      -v(i  ,j+1,8)/sqrt(u(i  ,j+1,8)*u(i  ,j+1,8)+v(i  ,j+1,8)*v(i  ,j+1,8))  &
!                      +v(i+1,j  ,8)/sqrt(u(i+1,j  ,8)*u(i+1,j  ,8)+v(i+1,j  ,8)*v(i+1,j  ,8))  &
!                      -v(i  ,j  ,8)/sqrt(u(i  ,j  ,8)*u(i  ,j  ,8)+v(i  ,j  ,8)*v(i  ,j  ,8))  &
!                     )/proj%dx
!          dudy = 0.5*( u(i+1,j+1,8)/sqrt(u(i+1,j+1,8)*u(i+1,j+1,8)+v(i+1,j+1,8)*v(i+1,j+1,8))  &
!                      -u(i+1,j  ,8)/sqrt(u(i+1,j  ,8)*u(i+1,j  ,8)+v(i+1,j  ,8)*v(i+1,j  ,8))  &
!                      +u(i  ,j+1,8)/sqrt(u(i  ,j+1,8)*u(i  ,j+1,8)+v(i  ,j+1,8)*v(i  ,j+1,8))  &
!                      -u(i  ,j  ,8)/sqrt(u(i  ,j  ,8)*u(i  ,j  ,8)+v(i  ,j  ,8)*v(i  ,j  ,8))  &
!                     )/proj%dx
!          vor  = (dvdx-dudy)*(1030.-slp(i,j))**4
!          if( vor .gt. vor1 )then
!              center(1) = i*1.
!              center(2) = j*1.
!              vor1 = vor
!          endif
!      enddo
!      enddo
!      center(3) = slp(int(center(1)), int(center(2)))
!   endif
!
!!.... find the minumum slp
!   search_radiu = int(150./proj%dx*1000.)
!   i1 = max(int(center(1))-search_radiu, 5)
!   i2 = min(int(center(1))+search_radiu, ix-5)
!   j1 = max(int(center(2))-search_radiu, 5)
!   j2 = min(int(center(2))+search_radiu, jx-5)
!   do j = j1, j2
!   do i = i1, i2
!      if( slpsmth(i,j) < center(3) )then
!          center(3) = slpsmth(i,j)
!          center(1) = i*1.
!          center(2) = j*1.
!      endif
!   enddo
!   enddo
!   center(3) = slp(int(center(1)), int(center(2)))
!
!.... smooth
!   call open_file( wrf_file, nf_nowrite, fid )
!   rcode = nf_inq_varid(     fid, 'U10       ', varid )
!   if ( rcode .eq. 0 ) then
!      call get_variable2d( wrf_file, 'U10       ', ix, jx, 1, u10 )
!      call get_variable2d( wrf_file, 'V10       ', ix, jx, 1, v10 )
!   else
!      u10(1:ix,1:jx)=u(1:ix,1:jx,1)
!      v10(1:ix,1:jx)=v(1:ix,1:jx,1)
!   endif
!   call smooth(u10, ix, jx, 301)
!   call smooth(v10, ix, jx, 301)
!
!   wsp1 = 1.**2.
!   search_radiu = int(500./proj%dx*1000.)
!   i1 = max(int(center(1))-search_radiu, 5)
!   i2 = min(int(center(1))+search_radiu, ix-5)
!   j1 = max(int(center(2))-search_radiu, 5)
!   j2 = min(int(center(2))+search_radiu, jx-5)
!   do j = j1, j2 
!   do i = i1, i2
!   do j = jsearchs, jsearche
!   do i = isearchs, isearche
!
!      wsp  = u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j)
!      if( wsp .gt.  wsp1 )then  
!          wsp1 = wsp
!      endif
!   enddo
!   enddo
!   center(4) = sqrt(wsp1)
!
!  end subroutine hurricane_center_wind
!
!==============================================================================
