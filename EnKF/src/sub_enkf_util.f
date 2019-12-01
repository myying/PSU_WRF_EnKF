!=====================================================================================
subroutine  read_namelist(ix, jx, kx)

!----------------------------------------------------------------------------
!  PURPOSE:       Read in EnKF input parameters from namelist.
!  METHOD:        Open namelist input file and read in values.
!
!  Origin:        03/09/2006 (Yonghui Weng)
!----------------------------------------------------------------------------

   use constants
   use namelist_define
   use mpi_module
   use mapinfo_define
   use obs_define
   implicit none

!  Local scalars:

   integer, intent(in) :: ix, jx, kx
   character(len=80)   :: namelist_file       ! Input namelist filename.
   integer, parameter  :: namelist_unit=7   ! Input namelist unit.
   integer             :: iost               ! Error code.
   integer             :: i

!------------------------------------------------------------------------------
!  [1.0] namelist initialize
!------------------------------------------------------------------------------
!-- enkf_parameter
   numbers_en   = 999
   expername    = '          '
   enkfvar      = '          '
   updatevar    = '          '
   update_is    = 999
   update_ie    = 999
   update_js    = 999 
   update_je    = 999 
   update_ks    = 999 
   update_ke    = 999 
   inflate      = 0.0
   relax_opt    = 0
   relax_adaptive = .false.
   mixing       = 0.0
   print_detail  = 2
!--parallel
   manual_parallel = .false.
   nmcpu       = 1
   nicpu       = 1
   njcpu       = 1
   random_order = .true.

!-- multiscale
   num_scales = 1
   krange(1) = 1
   current_scale = 1
   run_alignment = .false.

!-- use_hurricane_position_intensity
   use_hurricane_PI = .false.      ! hurricane position(center lat and lon) and inteneity(min SLP)
   hroi_hurricane_PI = 999
   vroi_hurricane_PI = 999

!-- use_surface_obs
   use_surface      = .false.      ! SURFACE LAND (SYNOPTIC, METAR) REPORTS
   datathin_surface = 999
   hroi_surface     = 999
   vroi_surface     = 999

!-- use_sounding_obs
   use_sounding      = .false. 
   datathin_sounding = 999
   datathin_sounding_vert = 999
   hroi_sounding     = 999
   vroi_sounding     = 999

!-- use_profiler_obs
   use_profiler      = .false. 
   datathin_profiler = 999
   datathin_profiler_vert = 999
   hroi_profiler     = 999
   vroi_profiler     = 999

!-- use_aircft_obs
   use_aircft      = .false. 
   datathin_aircft = 999
   hroi_aircft     = 999
   vroi_aircft     = 999

!-- use_metar_obs
   use_metar       = .false. 
   datathin_metar  = 999
   hroi_metar      = 999
   vroi_metar      = 999

!-- use_sfcshp_obs
   use_sfcshp      = .false. 
   datathin_sfcshp = 999
   hroi_sfcshp     = 999
   vroi_sfcshp     = 999

!-- use_spssmi_obs
   use_spssmi      = .false. 
   datathin_spssmi = 999
   hroi_spssmi     = 999
   vroi_spssmi     = 999

!-- use_atovs_obs
   use_atovs         = .false. 
   datathin_atovs    = 999
   hroi_atovs        = 999
   vroi_atovs        = 999

!-- use_satwnd_obs
   use_satwnd      = .false. 
   datathin_satwnd = 999
   hroi_satwnd     = 999
   vroi_satwnd     = 999

!-- use_gpspw_obs
   use_gpspw      = .false. 
   datathin_gpspw = 999
   hroi_gpspw     = 999
   vroi_gpspw     = 999

!-- use_groundbase_radar
   radar_number   = 999
   use_radar_rf   = .false. 
   use_radar_rv   = .false. 
   datathin_radar = 999
   hroi_radar     = 999
   vroi_radar     = 999

!-- use_airborne_radar
   use_airborne_rf   = .false. 
   use_airborne_rv   = .false. 
   datathin_airborne = 999
   hroi_airborne     = 999
   vroi_airborne     = 999
  
!-- use_radiance
   use_radiance   = .false.
   datathin_radiance = 999
   hroi_radiance     = 999
   vroi_radiance     = 999

!-- use_seawind
   use_seawind   = .false.
   datathin_seawind = 999
   hroi_seawind     = 999
   vroi_seawind     = 999


 
!------------------------------------------------------------------------
!  [2.0] read namelist 
!------------------------------------------------------------------------

   namelist_file = 'namelist.enkf'
   iost = 0
   open ( file = namelist_file, unit = namelist_unit,                   &
          status = 'old', access = 'sequential', iostat = iost ) 
   if( iost .ne. 0 ) then
       write(*,*)'namelist.enkf does not exist, please check it.'
       stop 'read_namelist'
   endif

!-- enkf_parameter
   iost = 0
   read ( unit = namelist_unit, nml = enkf_parameter, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'enkf_parameter, please check it.'
       stop 'read_namelist enkf_parameter'
   endif
   read ( unit = namelist_unit, nml = parallel, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'parallel, please check it.'
       stop 'read_namelist enkf_parameter'
   endif
   read ( unit = namelist_unit, nml = multiscale, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'multiscale, please check it.'
       stop 'read_namelist enkf_parameter'
   endif
   read ( unit = namelist_unit, nml = hurricane_PI, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'hurricane_PI, please check it.'
       stop 'read_namelist hurricane_PI'
   endif
   read ( unit = namelist_unit, nml = surface_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'surface_obs, please check it.'
       stop 'read_namelist surface_obs'
   endif
   read ( unit = namelist_unit, nml = sounding_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'sounding_obs, please check it.'
       stop 'read_namelist sounding_obs'
   endif
   read ( unit = namelist_unit, nml = profiler_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'profiler_obs, please check it.'
       stop 'read_namelist profiler_obs'
   endif
   read ( unit = namelist_unit, nml = aircft_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'aircft_obs, please check it.'
       stop 'read_namelist aircft_obs'
   endif
   read ( unit = namelist_unit, nml = metar_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'metar_obs, please check it.'
       stop 'read_namelist metar_obs'
   endif
   read ( unit = namelist_unit, nml = sfcshp_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'sfcshp_obs, please check it.'
       stop 'read_namelist sfcshp_obs'
   endif
   read ( unit = namelist_unit, nml = spssmi_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'spssmi_obs, please check it.'
       stop 'read_namelist spssmi_obs'
   endif
   read ( unit = namelist_unit, nml = atovs_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'atovs_obs, please check it.'
       stop 'read_namelist atovs_obs'
   endif
   read ( unit = namelist_unit, nml = satwnd_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'satwnd_obs, please check it.'
       stop 'read_namelist satwnd_obs'
   endif
   read ( unit = namelist_unit, nml = seawind_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'seawind_obs, please check it.'
       stop 'read_namelist seawind_obs'
   endif
   read ( unit = namelist_unit, nml = gpspw_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'gpspw_obs, please check it.'
       stop 'read_namelist gpspw_obs'
   endif
   read ( unit = namelist_unit, nml = radar_obs, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'radar_obs, please check it.'
       stop 'read_namelist radar_obs'
   endif
   read ( unit = namelist_unit, nml = airborne_radar, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'airborne_radar, please check it.'
       stop 'read_namelist airborne_radar'
   endif
   read ( unit = namelist_unit, nml = radiance, iostat = iost )
   if( iost .ne. 0 ) then
       write(*,*)'radiance, please check it.'
       stop 'read_namelist radiance'
   endif

   if( update_is <= 1 .or. update_is >=ix ) update_is = 1
   if( update_ie <= 1 .or. update_ie >=ix ) update_ie = ix
   if( update_js <= 1 .or. update_js >=jx ) update_js = 1
   if( update_je <= 1 .or. update_je >=jx ) update_je = jx
   if( update_ks <= 1 .or. update_ks >=kx ) update_ks = 1
   if( update_ke <= 1 .or. update_ke >=kx ) update_ke = kx

   if ( my_proc_id == 0 ) then
      write(6, *)'   '
      write(6, *)'---------------------------------------------------'
      write(6, *)'.... namelist info ....'
      write(6,'(a,3x,i5    )') ' numbers_en = ',numbers_en 
      write(6,'(a,3x,a10   )') ' expername  = ',expername
      write(6,'(a,3x,20a   )') ' enkfvar = ',enkfvar
      write(6,'(a,3x,20a   )') ' updatevar = ',updatevar
      write(6,'(a18,i3,5(a1,i3))') ' update domain  : ',update_is,':',update_ie,'; ',update_js,':',update_je,'; ',update_ks,':',update_ke
      write(6,'(a,3x,f5.2   )') ' inflate    = ',inflate  
      write(6,'(a,3x,f5.2   )') ' mixing     = ',mixing  
      write(6,'(a,3x,i5    )') ' print_detail    = ',print_detail

      if ( use_hurricane_PI ) then
         write(6,'(a      )') ' ===== assimilate hurricane position and intensity ===== '
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_hurricane_PI, vroi_hurricane_PI
      endif

      if ( use_surface ) then
         write(6,'(a      )') ' ===== assimilate SURFACE LAND (SYNOPTIC, METAR) REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_surface
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_surface, vroi_surface
      endif

      if ( use_sounding ) then
         write(6,'(a      )') ' ===== assimilate UPPER-AIR (RAOB, PIBAL, RECCO, DROPS) REPORTS ===== '
         write(6,'(a,2i4   )') '       data thinning: ',datathin_sounding, datathin_sounding_vert
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_sounding, vroi_sounding
      endif

      if ( use_profiler ) then
         write(6,'(a      )') ' ===== assimilate PROFILER REPORTS ===== '
         write(6,'(a,2i4   )') '       data thinning: ',datathin_profiler, datathin_profiler_vert
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_profiler, vroi_profiler
      endif

      if ( use_aircft ) then
         write(6,'(a      )') ' ===== assimilate AIREP/PIREP, AMDAR (ASDAR/ACARS), E-ADAS (AMDAR BUFR) AIRCRAFT REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_aircft
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_aircft, vroi_aircft
      endif

      if ( use_metar ) then
         write(6,'(a      )') ' ===== assimilate AUTO-METEO-STATION REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_metar
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_metar , vroi_metar 
      endif

      if ( use_sfcshp ) then
         write(6,'(a      )') ' ===== assimilate SURFACE MARINE (SHIP, BUOY, C-MAN PLATFORM) REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_sfcshp
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_sfcshp, vroi_sfcshp
      endif

      if ( use_spssmi ) then
         write(6,'(a      )') ' ===== assimilate DMSP SSM/I RETRIEVAL PRODUCTS (REPROCESSED WIND SPEED, TPW) ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_spssmi
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_spssmi, vroi_spssmi
      endif

      if ( use_atovs ) then
         write(6,'(a      )') ' ===== assimilate ATOVS retrievals ===== '
         write(6,'(a,2i4   )') '       data thinning: ',datathin_atovs, datathin_atovs_vert
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_atovs, vroi_atovs
      endif

      if ( use_satwnd ) then
         write(6,'(a      )') ' ===== assimilate SATELLITE-DERIVED WIND REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_satwnd
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_satwnd, vroi_satwnd
      endif

      if ( use_seawind ) then
         write(6,'(a      )') ' ===== assimilate SEA SURFACE WIND ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_seawind
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_seawind, vroi_seawind
      endif

      if ( use_gpspw ) then
         write(6,'(a      )') ' ===== assimilate GPS PRECIP WATER REPORTS ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_gpspw
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_gpspw, vroi_gpspw
      endif

      if ( use_radiance ) then
         write(6,'(a      )') ' ===== assimilate RADIANCE ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_radiance
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_radiance, vroi_radiance
      endif

      if ( use_radar_rv .or. use_radar_rf ) then
         if ( use_radar_rf ) write(*,'(a)') ' ===== assimilate NEXRAD reflectivity ===== '
         if ( use_radar_rv ) write(*,'(a)') ' ===== assimilate NEXRAD radial velocity ===== '
         write(6,'(a,i2,a )') '       ',radar_number,' radar data will be assimilated '
         write(6,'(a,i4   )') '       data thinning: ',datathin_radar
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_radar, vroi_radar
      endif

      if ( use_airborne_rv .or. use_airborne_rf ) then
         if ( use_airborne_rf ) write(*,'(a)') ' ===== assimilate airborne radar reflectivity ===== '
         if ( use_airborne_rv ) write(*,'(a)') ' ===== assimilate airborne radar radial velocity ===== '
         write(6,'(a,i4   )') '       data thinning: ',datathin_airborne
         write(6,'(a,2i4  )') '       ROI for horizonal and vertical is:', hroi_airborne, vroi_airborne
      endif

   endif

   raw%radar_stn_num             = radar_number
   allocate( raw%radar( radar_number ) )

   close ( unit = namelist_unit )
   return
end subroutine  read_namelist 

!========================================================================================
subroutine get_wrf_info(wrf_file, ix, jx, kx, times, proj)

!-----------------------------------------------------------------------------
!  get wrf output file's information such as ix,jx,kx, eta levels, proj
!-----------------------------------------------------------------------------

   use constants
   use namelist_define
   use mpi_module
   use mapinfo_define
   use netcdf
   use wrf_tools
   implicit none

   type(proj_info)              :: proj
   character (len=10), intent(in)  :: wrf_file
   character (len=80), intent(out) :: times
   integer, intent(out)            :: ix, jx, kx

   if ( my_proc_id == 0 ) then
      write(6, *)'   '
      write(6, *)'---------------------------------------------------'
      write(6, *)'.... Model Domain Information ....'
   endif
!-----------------------------------------------------------------------------
!  get domain dimension
!-----------------------------------------------------------------------------
   !call get_var_ij ( wrf_file, 'T', ix, jx, kx )
   call get_ij ( wrf_file, ix, jx, kx )
   if ( my_proc_id == 0 ) write(6,'(3(a7,i3))')"   ix :",ix,"   jx :",jx,"   kx:",kx

!-----------------------------------------------------------------------------
!  get times
   call get_times ( wrf_file, 'T', times )
   if ( my_proc_id == 0 ) write(6,'(a,a19)')" 1st guest field time : ",times(1:19)

!-----------------------------------------------------------------------------
!  get map info
   call set_domain_proj( wrf_file, proj )
   if ( my_proc_id == 0 ) then
      write(6,'(a,i1)')" Domain projection type : ", proj%code
      write(6,'(a,f7.3,a,f8.3,a)')" left_bottom(1,1) : (", proj%lat1, ",",proj%lon1,")"
      write(6,'(a,f7.0)')" Grid spacing in meters : ", proj%dx
      write(6,'(a,f8.3)')" Standard Longitude : ", proj%stdlon
      write(6,'(a,2f8.3)')" True Latitude      : ", proj%truelat1, proj%truelat2
   endif

end subroutine get_wrf_info

!========================================================================================
subroutine rd_truth_wrf(wrf_file,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,x)
   use namelist_define
   use wrf_tools
   use mpi_module
   implicit none

   character (len=10), intent(in) :: wrf_file
   integer, intent(in)  :: ix, jx, kx,ni,nj,nk,nv,nm,gid,sid,iid,jid 
   real, dimension(ni,nj,nk,nv), intent(out) :: x
   integer  :: ii, jj, kk, m, istart,iend,jstart,jend
   real, allocatable, dimension(:,:,:) :: dat3d

   do m=1,nv
      call wrf_var_dimension ( wrf_file, enkfvar(m), ix, jx, kx, ii, jj, kk )
      istart=iid*ni+1
      iend=(iid+1)*ni
      jstart=jid*nj+1
      jend=(jid+1)*nj
      if(iid==(nicpu-1)) iend=ii
      if(jid==(njcpu-1)) jend=jj
      allocate( dat3d ( ii, jj, kk ) )
      if ( kk > 1 ) then
         call get_variable3d( wrf_file, enkfvar(m), ii, jj, kk, 1, dat3d )
      else if ( kk == 1 ) then
         call get_variable2d( wrf_file, enkfvar(m), ii, jj, 1, dat3d )
      endif
      x(1:iend-istart+1,1:jend-jstart+1,1:kk,m) = dat3d(istart:iend,jstart:jend,1:kk)
      deallocate( dat3d )
   end do

end subroutine rd_truth_wrf

!========================================================================================
subroutine read_ensemble(i_unit,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid,x)
use namelist_define
use wrf_tools
use mpi_module
implicit none
integer, intent(in)  :: i_unit,ix,jx,kx,ni,nj,nk,nv,nm,gid,sid,iid,jid
real, dimension(ni,nj,nk,nv,nm), intent(inout) :: x
integer :: ii,jj,kk,ie,n,m,istart,iend,jstart,jend
character(len=10) :: wrf_file
real, allocatable, dimension(:,:,:) :: dat3d

if ( my_proc_id == 0 ) then
  write(*, *)'   '
  write(*, *)'---------------------------------------------------'
  write(*, *)'.... Getting initial ensemble filed ....'
endif
do n=1,nm
  ie=(n-1)*nmcpu+gid+1
  if(ie<=numbers_en) then
    write(wrf_file,'(a5,i5.5)') 'fort.',i_unit+ie
    do m=1,nv
      call wrf_var_dimension(wrf_file,enkfvar(m),ix,jx,kx,ii,jj,kk)
      istart=iid*ni+1
      iend=(iid+1)*ni
      jstart=jid*nj+1
      jend=(jid+1)*nj
      if(iid==(nicpu-1)) iend=ii
      if(jid==(njcpu-1)) jend=jj
      allocate(dat3d(ii,jj,kk))
      if(kk>1) then
        call get_variable3d(wrf_file,enkfvar(m),ii,jj,kk,1,dat3d)
      else if (kk==1) then
        call get_variable2d(wrf_file,enkfvar(m),ii,jj,1,dat3d)
      endif
      x(1:iend-istart+1,1:jend-jstart+1,1:kk,m,n)=dat3d(istart:iend,jstart:jend,1:kk)
      deallocate(dat3d)
    enddo
  endif
enddo
end subroutine read_ensemble

!========================================================================================
!subroutine check_obs_location(ix,jx,obs_ii,obs_jj,roi,indomain)
! check if obs is within (1) cpu's subdomain and (2) update_var domain defined by update_is,ie,js,je
! return indomain=1(True) or 0(False)
!   use namelist_define 
!   use mpi_module
!   implicit none
!   integer, intent(in)           :: ix,jx,nicpu,njcpu,micpu,mjcpu,mkcpu,roi
!   integer                       :: sid,iid,jid,gid,ncpugroups
!   integer                       :: deli,delj,istart,iend,jstart,jend
!   real, intent(in)              :: obs_ii,obs_jj
!   integer, intent(out)          :: indomain
!
!   indomain=1
!
!   if( obs_ii < real(update_is-roi) .or. obs_ii > real(update_ie+roi) .or.  &
!       obs_jj < real(update_js-roi) .or. obs_jj > real(update_je+roi) ) then
!      indomain=0
!   endif
!
!   gid=int(my_proc_id/(nicpu*njcpu*micpu*mjcpu*mkcpu))
!   ncpugroups=nprocs/(nicpu*njcpu*micpu*mjcpu*mkcpu)
!   sid=mod(my_proc_id, nicpu*njcpu*micpu*mjcpu*mkcpu)
!   jid=int(int(sid/(micpu*mjcpu*mkcpu))/nicpu)
!   iid=mod(int(sid/(micpu*mjcpu*mkcpu)),nicpu)

!   deli=int(ix/nicpu)
!   delj=int(jx/njcpu)
!   istart=iid*deli+1
!   iend=(iid+1)*deli
!   jstart=jid*delj+1
!   jend=(jid+1)*delj
!   if(iid==(nicpu-1)) iend=ix
!   if(jid==(njcpu-1)) jend=jx
!   if( obs_ii < real(istart-roi) .or. obs_ii > real(iend+roi) .or. &
!       obs_jj < real(jstart-roi) .or. obs_jj > real(jend+roi) ) then
!      indomain=0
!   endif
!
!end subroutine check_obs_location

!========================================================================================
subroutine output(o_unit,ix,jx,kx,ni,nj,nk,nv,nm,s_comm,sid,iid,jid,xout,times)
use namelist_define
use netcdf
use wrf_tools
use mpi_module
implicit none
integer, intent(in) :: o_unit,ix,jx,kx,ni,nj,nk,nv,nm,sid,iid,jid,s_comm
real,dimension(ni,nj,nk,nv), intent(in) :: xout
character (len=80), intent(in) :: times
character (len=10) :: filename
real, allocatable, dimension(:,:,:) :: dat3d,dat3drecv
integer :: ii,jj,kk,ie,n,m,fid,istart,iend,jstart,jend
write(filename,'(a5,i5.5)') 'fort.',o_unit
if(sid==0) &
  call open_file(filename,nf_write,fid)
do m=1,nv
  call wrf_var_dimension(filename,enkfvar(m),ix,jx,kx,ii,jj,kk)
  allocate(dat3d(ii,jj,kk))
  allocate(dat3drecv(ii,jj,kk))
  dat3d=0.
  dat3drecv=0.
  istart=iid*ni+1
  iend=(iid+1)*ni
  jstart=jid*nj+1
  jend=(jid+1)*nj
  if(iid==(nicpu-1)) iend=ii
  if(jid==(njcpu-1)) jend=jj
  dat3d(istart:iend,jstart:jend,1:kk)=xout(1:iend-istart+1,1:jend-jstart+1,1:kk,m)
  call MPI_Allreduce(dat3d,dat3drecv,ii*jj*kk,MPI_REAL,MPI_SUM,s_comm,ierr)
  if(sid==0) then
    if(kk>1) then
      call write_variable3d(fid,enkfvar(m),ii,jj,kk,1,dat3drecv)
    else if (kk==1) then
      call write_variable2d(fid,enkfvar(m),ii,jj,1,dat3drecv)
    endif
  endif
  deallocate(dat3d)
  deallocate(dat3drecv)
enddo
if(sid==0) &
  call close_file(fid)
end subroutine output
!==============================================================================
subroutine wrf_var_dimension ( wrf_file, var, ix, jx, kx, ii, jj, kk )
   use netcdf

   character (len=*), intent(in) :: wrf_file
   character(len=10), intent(in)   :: var
   integer, intent(in)             :: ix, jx, kx
   integer, intent(out)            :: ii, jj, kk

   ii = ix
   jj = jx
   kk = kx
   if      ( var == 'U         ' ) then
      ii = ix + 1
   else if ( var == 'V         ' ) then
      jj = jx + 1
   else if ( var == 'W         ' .or. var == 'PH        ' .or. var == 'PHB       ' ) then
      kk = kx + 1
   else if ( var == 'MU        ' .or. var == 'MUB       ' .or. var == 'Q2        '  &
        .or. var == 'T2        ' .or. var == 'TH2       ' .or. var == 'PSFC      '  &
        .or. var == 'SST       ' .or. var == 'TSK       ' .or. var == 'XICE      '  &
        .or. var == 'SFROFF    ' .or. var == 'UDROFF    ' .or. var == 'IVGTYP    '  &
        .or. var == 'ISLTYP    ' .or. var == 'VEGFRA    ' .or. var == 'GRDFLX    '  &
        .or. var == 'SNOW      ' .or. var == 'SNOWH     ' .or. var == 'CANWAT    '  &
        .or. var == 'SST       ' .or. var == 'MAPFAC_M  ' .or. var == 'F         '  &
        .or. var == 'E         ' .or. var == 'SINALPHA  ' .or. var == 'COSALPHA  '  &
        .or. var == 'HGT       ' .or. var == 'TSK       ' .or. var == 'RAINC     '  &
        .or. var == 'RAINNC    ' .or. var == 'SWDOWN    ' .or. var == 'GLW       '  &
        .or. var == 'XLAT      ' .or. var == 'XLONG     ' .or. var == 'TMN       '  &
        .or. var == 'XLAND     ' .or. var == 'PBLH      ' .or. var == 'HFX       '  &
        .or. var == 'QFX       ' .or. var == 'LH        ' .or. var == 'SNOWC     '  &
        .or. var == 'SR        ' .or. var == 'POTEVP    ' .or. var == 'U10       '  &
        .or. var == 'V10       ' ) then
      kk = 1
   else if ( var == 'MAPFAC_U  ' ) then
      kk = 1
      ii = ix + 1
   else if ( var == 'MAPFAC_V  ' ) then
      kk = 1
      jj = jx + 1
   else if ( var == 'FNM       ' .or. var == 'FNP       '  &
        .or. var == 'RDNW      ' .or. var == 'RDN       '  &
        .or. var == 'DNW       ' .or. var == 'DN        '  &
        .or. var == 'ZNU       '                          ) then
      ii = 1
      jj = 1
   else if ( var == 'ZNW       '                          ) then
      ii = 1
      jj = 1
      kk = kx + 1
   endif
      
   if( var == 'TSLB' .or. var == 'SMOIS' ) then
         call get_soilkk(wrf_file,kk)
   end if

end subroutine wrf_var_dimension

!==============================================================================
 function gaussdev(idum)
! Returns a normally distributed deviate with 0 mean and unit variance,
! using ran3(idum) as the source of uniform deviates; see Num'l Recipes,
! ch 7.2.
      integer idum
      real gaussdev

!- local variables
      integer iset
      real fac,gset,rsq,v1,v2,ran1
      save iset, gset

      data iset/0/

      if (iset.eq.0) then
 10      v1 = 2.*ran3(idum)-1.
         v2 = 2.*ran3(idum)-1.

         rsq = v1**2 + v2**2
         if (rsq.ge.1. .or. rsq.eq.0.) goto 10

         fac = sqrt( -2.*log(rsq)/rsq )
         gset = v1*fac
         gaussdev = v2*fac

         iset= 1

      else
         gaussdev = gset
         iset=0

      end if
      return
 end function gaussdev

!==============================================================================
 function ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
 end function ran3 

!==============================================================================
  subroutine get_variable_from_xa( xa, varname, ntot, ix, jx, kx,    &
                                 var_start, num_var, work3d, i_var )
!---------------------
! get_variable_from_xa subroutine get one variable from xa
!---------------------

  use namelist_define

  implicit none
  integer, intent(in)                     :: ix, jx, kx, ntot       ! 1st guest dimension
  real, dimension(ntot), intent(in)       :: xa
  integer, intent(in)                     :: num_var
  integer, dimension(num_var), intent(in) :: var_start
  character(len=10), intent(in)           :: varname

  integer, intent(out)                    :: i_var
  real, dimension(ix,jx,kx)               :: work3d

  integer                                 :: m, itot, i, j, k

  do_variable_loop : do m = 1, num_var
    
     if ( enkfvar(m) == varname ) then
        do k = 1, kx
        do j = 1, jx
        do i = 1, ix
           itot = ix*jx*(k-1) + ix*(j-1) + i - 1
           work3d(i, j, k) = xa( itot + var_start(m) )
        enddo
        enddo
        enddo
    
        i_var = 1
        exit do_variable_loop
     endif

  end do do_variable_loop

  end subroutine get_variable_from_xa

!==============================================================================
  subroutine get_variable_from_xa_point( xa, varname, ntot, ix, jx, kx,    &
                                 var_start, num_var, work1d, i_var, iob, job )
  use namelist_define

  implicit none
  integer, intent(in)                     :: ix, jx, kx, ntot, iob, job
  real, dimension(ntot), intent(in)       :: xa
  integer, intent(in)                     :: num_var
  integer, dimension(num_var), intent(in) :: var_start
  character(len=10), intent(in)           :: varname
  integer, intent(out)                    :: i_var
  real, dimension(kx)                     :: work1d

  integer                                 :: m, itot, i, j, k

  do_variable_loop : do m = 1, num_var

     if ( enkfvar(m) .eq. varname ) then
        do k = 1, kx
           itot = ix*jx*(k-1) + ix*(job-1) + iob-1
           work1d(k) = xa( itot + var_start(m) )
        enddo

        i_var = 1
        exit do_variable_loop
     endif

  end do do_variable_loop

  end subroutine get_variable_from_xa_point

!======================================================================================= 
!  subroutine seaprsnh(t, ter, ps, psfc, p, znu, ix, jx, kx, ptop, slp) 
      SUBROUTINE SEAPRSNH(T,TER,PS,psfc,PP,SIGMA,IMX,JMX,KX,PTOP, pslv ) 
!     SECTION  DIAGNOSTIC
!     PURPOSE  COMPUTES SEA LEVEL PRESSURE FROM THE RULE
!              T1/T2=(P1/P2)**(GAMMA*R/G).
!
!     INPUT       T        TEMPERATURE                CROSS    3D
!                 TER      TERRAIN                    CROSS    2D
!                 PS       P STAR = PSFC-PTOP         CROSS    2D
!                 SFP      SURFACE PRESSURE           CROSS    2D
!                 SIG      HALF SIGMA LEVELS                   1D
!                 IMX      DOT POINT DIMENSION N-S
!                 JMX      DOT POINT DIMENSION E-W
!                 KX       NUMBER OF VERTICAL LEVELS
!                 PTOP     PRESSURE AT TOP OF MODEL
!
!     OUTPUT      SLP      SEA LEVEL PRESSURE         CROSS    2D
!
      DIMENSION T(IMX,JMX,KX),PP(IMX,JMX,KX),                     &
                PS(IMX,JMX)  ,psfc(IMX,JMX) ,                     &
                TER(IMX,JMX) ,SIGMA(KX)
      DIMENSION PL(IMX,JMX),T0(IMX,JMX),TS(IMX,JMX),              &
                XKLEV(IMX,JMX)
      DIMENSION pslv(IMX,JMX)

      LOGICAL L1,L2,L3,L4

      PARAMETER (TC=273.15+17.5) ! T CRITICAL IN PSFC/PSLV
      PARAMETER (PCONST=100.)    ! HEIGHT (MB ABOVE SURFACE)
!                                  TO BEGIN EXTRAPOLATION
      PARAMETER ( R     = 287.04,                                 &
     &            CP    = 1004.0,                                 &
     &            ROVCP = R/CP,                                   &
     &            CPOVR = CP/R,                                   &
     &            GAMMA = 6.5E-3,                                 &
     &            G     = 9.81,                                   &
     &            EPS   = 0.622,                                  &
     &            XLVC1 = 3.1484E6,                               &
     &            XLVC2 = 2.37E3,                                 &
     &            SVP1  = 0.6112,                                 &
     &            SVP2  = 17.67,                                  &
     &            SVP3  = 29.65,                                  &
     &            SVPT0 = 273.15)
!
!     ... SEA LEVEL PRESSURE
!
      XTERM=GAMMA*R/G
!
!     ... COMPUTE PRESSURE AT PCONST MB ABOVE SURFACE (PL)
!
      KUPTO=KX/2
99    CONTINUE
      DO 100 J=1,JMX-1
      DO 100 I=1,IMX-1
         PL(I,J)=PSFC(I,J)-PCONST
         XKLEV(I,J)=0.
100   CONTINUE
!
!     ... FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J
!
      DO 150 J=1,JMX-1
      DO 150 I=1,IMX-1
         DO 125 K=KX-1,KUPTO,-1
            IF(((SIGMA(K  )*PS(I,J)+PTOP+PP(I,J,K)).GE.PL(I,J)) .AND.  &
            ((SIGMA(K+1)*PS(I,J)+PTOP+PP(I,J,K+1)).LT.PL(I,J)))        &
            XKLEV(I,J)=FLOAT(K)
125      CONTINUE
         IF(XKLEV(I,J).LT.1.) THEN
            PRINT *,'ERROR FINDING PRESSURE LEVEL ',PCONST,' MB ',     &
                    'ABOVE THE SURFACE'
            PRINT *,'LAST K LEVEL =',KUPTO
            IF(KUPTO.NE.1) THEN
               PRINT *,'TRYING AGAIN WITH KUPTO=1'
               KUPTO=1
               GOTO 99
            ELSE
               PRINT *,'I,J=',I,J
               PRINT *,'PL=',PL(I,J)
               PRINT *,'PSFC=',PSFC(I,J)
               PRINT *,'PP=',PP(I,J,1:KX)
               PRINT *,'SIGMA=',SIGMA(1:KX)
               PRINT *,SIGMA(K  )*PS(I,J)+PTOP+PP(I,J,K),SIGMA(K+1)*PS(I,J)+PTOP+PP(I,J,K+1)
               CALL ABORT
            END IF
         END IF
150   CONTINUE
!
!     ... GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)
!         AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE
!
      DO 200 J=1,JMX-1
      DO 200 I=1,IMX-1
!         KLO=NINT(XKLEV(I,J))+1
!         KHI=NINT(XKLEV(I,J))
         KLO=NINT(XKLEV(I,J))
         KHI=NINT(XKLEV(I,J))+1
         PLO=SIGMA(KLO)*PS(I,J)+PTOP+PP(I,J,KLO)
         PHI=SIGMA(KHI)*PS(I,J)+PTOP+PP(I,J,KHI)
         TLO=T(I,J,KLO)
         THI=T(I,J,KHI)
         TL=THI-(THI-TLO)*ALOG(PL(I,J)/PHI)/ALOG(PLO/PHI)
         TS(I,J)=TL*(PSFC(I,J)/PL(I,J))**XTERM
         TBAR=(TS(I,J)+TL)*0.5
         HL=TER(I,J)-R/G*ALOG(PL(I,J)/PSFC(I,J))*TBAR
         T0(I,J)=TL+GAMMA*HL
200   CONTINUE
!
!     ... CORRECT SEA LEVEL TEMPERATURE IF TOO HOT
!
      DO 400 J=1,JMX-1
      DO 400 I=1,IMX-1
         L1=T0(I,J).LT.TC
         L2=TS(I,J).LE.TC
         L3=.NOT.L1
         IF(L2.AND.L3)THEN
            T0(I,J)=TC
         ELSE IF((.NOT. L2).AND. L3) THEN
            T0(I,J)=TC-0.005*(TS(I,J)-TC)**2
         ENDIF
400   CONTINUE
!
!     ... COMPUTE SEA LEVEL PRESSURE
!
      DO 600 J=1,JMX-1
      DO 600 I=1,IMX-1
         pslv(I,J)=PSFC(I,J)*EXP(2.*G*TER(I,J)/(R*(TS(I,J)+T0(I,J))))
600   CONTINUE
      RETURN
      END

!=======================================================================================
   subroutine compute_seaprs(nx, ny, nz, z, t, p, q, sea_level_pressure, debug) 
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave

!      subroutine compute_seaprs ( nx , ny , nz  ,         &
!                                  z, t , p , q ,          &
!                                  sea_level_pressure,debug)
!     &                            t_sea_level, t_surf, level )
      IMPLICIT NONE
!     Estimate sea level pressure.
      INTEGER nx , ny , nz
      REAL    z(nx,ny,nz)
      REAL    t(nx,ny,nz) , p(nx,ny,nz) , q(nx,ny,nz)
!     The output is the 2d sea level pressure.
      REAL    sea_level_pressure(nx,ny)
      INTEGER level(nx,ny)
      REAL t_surf(nx,ny) , t_sea_level(nx,ny)
      LOGICAL debug

!     Some required physical constants:

      REAL R, G, GAMMA
      PARAMETER (R=287.04, G=9.81, GAMMA=0.0065)

!     Specific constants for assumptions made in this routine:

      REAL    TC, PCONST
      PARAMETER (TC=273.16+17.5, PCONST = 10000)
      LOGICAL ridiculous_mm5_test
      PARAMETER (ridiculous_mm5_test = .TRUE.)
!      PARAMETER (ridiculous_mm5_test = .false.)

!     Local variables:
      INTEGER i , j , k
      INTEGER klo , khi


      REAL plo , phi , tlo, thi , zlo , zhi
      REAL p_at_pconst , t_at_pconst , z_at_pconst
      REAL z_half_lowest

      REAL    , PARAMETER :: cp           = 7.*R/2.
      REAL    , PARAMETER :: rcp          = R/cp
      REAL    , PARAMETER :: p1000mb      = 100000.

      LOGICAL  l1 , l2 , l3, found

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

      t(:,:,:) = (t(:,:,:)+300.)*(p(:,:,:)/p1000mb)**rcp

      DO j = 1 , ny
         DO i = 1 , nx
            level(i,j) = -1

            k = 1
            found = .false.
            do while( (.not. found) .and. (k.le.nz))
               IF ( p(i,j,k) .LT. p(i,j,1)-PCONST ) THEN
                  level(i,j) = k
                  found = .true.
               END IF
               k = k+1
            END DO

            IF ( level(i,j) .EQ. -1 ) THEN
            PRINT '(A,I4,A)','Troubles finding level ',   &
                        NINT(PCONST)/100,' above ground.'
            PRINT '(A,I4,A,I4,A)',                        &
                  'Problems first occur at (',i,',',j,')'
            PRINT '(A,E13.5,A)',                           &
                  'Surface pressure = ',p(i,j,1)/100,' hPa.'
            STOP 'Error_in_finding_100_hPa_up'
         END IF
         END DO
      END DO

!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

      DO j = 1 , ny
         DO i = 1 , nx

            klo = MAX ( level(i,j) - 1 , 1      )
            khi = MIN ( klo + 1        , nz - 1 )

            IF ( klo .EQ. khi ) THEN
               PRINT '(A)','Trapping levels are weird.'
               PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
                            ': and they should not be equal.'
               STOP 'Error_trapping_levels'
            END IF

         plo = p(i,j,klo)
         phi = p(i,j,khi)
         tlo = t(i,j,klo)*(1. + 0.608 * q(i,j,klo) )
         thi = t(i,j,khi)*(1. + 0.608 * q(i,j,khi) )
!         zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
!         zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         zlo = z(i,j,klo)
         zhi = z(i,j,khi)

         p_at_pconst = p(i,j,1) - pconst
         t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
         z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

         t_surf(i,j) = t_at_pconst*(p(i,j,1)/p_at_pconst)**(gamma*R/g)
         t_sea_level(i,j) = t_at_pconst+gamma*z_at_pconst

         END DO
      END DO

!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

      IF ( ridiculous_mm5_test ) THEN
         DO j = 1 , ny
            DO i = 1 , nx
               l1 = t_sea_level(i,j) .LT. TC
               l2 = t_surf     (i,j) .LE. TC
               l3 = .NOT. l1
               IF ( l2 .AND. l3 ) THEN
                  t_sea_level(i,j) = TC
               ELSE
                  t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
               END IF
            END DO
         END DO
      END IF

!     The grand finale: ta da!

      DO j = 1 , ny
      DO i = 1 , nx
!         z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
         z_half_lowest=z(i,j,1)
         sea_level_pressure(i,j) = p(i,j,1) *              &
                               EXP((2.*g*z_half_lowest)/   &
                                   (R*(t_sea_level(i,j)+t_surf(i,j))))

         sea_level_pressure(i,j) = sea_level_pressure(i,j)*0.01

      END DO
      END DO

    if (debug) then
      print *,'sea pres input at weird location i=20,j=1,k=1'
      print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
      print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
      print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
      print *,'slp=',sea_level_pressure(20,1),     &
               sea_level_pressure(20,2),sea_level_pressure(20,3)
    endif
!      print *,'t=',t(10:15,10:15,1),t(10:15,2,1),t(10:15,3,1)
!      print *,'z=',z(10:15,1,1),z(10:15,2,1),z(10:15,3,1)
!      print *,'p=',p(10:15,1,1),p(10:15,2,1),p(10:15,3,1)
!      print *,'slp=',sea_level_pressure(10:15,10:15),     &
!         sea_level_pressure(10:15,10:15),sea_level_pressure(20,10:15)

      end subroutine compute_seaprs
!=======================================================================================
  subroutine smooth(pslab, njx, niy, numpas)
!
!   This is a smoothing routine, with several choices:
!
!   If numpas is between 1 and 99, then a 9-point weighted smoother is
!   applied numpas times.  The smoother follows equation 11-107 in
!   Haltiner and Williams. One pass completely removes 2-delta-x waves
!   on the interior.  On the outer row and column, and near missing
!   data points, smoothing is carried out in a manner that preserves
!   the domain average value of the field.
!
!   If numpas is between 101 and 199, then a smoother-desmoother is
!   applied (numpas-100) times.  One pass removes a large fraction
!   of the 2-delta-x component, but is not as harsh on longer
!   wavelengths as the 9-point smoother
!
!   If numpas is between 201 and 299, then the smoother-desmoother is
!   applied (numpas-200) times, and, after each pass, the data field
!   is forced to be non-negative.
!
!   If numpas is between 301 and 399, then a weighted
!   smoother is applied, in which the smoothed value
!   is given by a weighted average of values at
!   surrounding grid points.  The weighting function
!   is the Cressman weighting function:
!
!               w = ( D**2 - d**2 ) / ( D**2 + d**2 )
!
!   In the above, d is the distance (in grid increments)
!   of the neighboring point to the smoothing point, and
!   D is the radius of influence [in grid increments,
!   given by (numpas-300)].
!
!   If numpas is between 401 and 499, then the smoothing
!   is similar for numpas=301-399, except the weighting
!   function is the circular apperture diffraction function
!   (following a suggestion of Barnes et al. 1996):
!
!               w = bessel(3.8317*d/D)/(3.8317*d/D)
!
!   If numpas is between 501 and 599, then the smoothing
!   is similar for numpas=301-399, except the weighting
!   function is the product of the rectangular
!   apperture diffraction function in the x and y directions
!   (the function used in Barnes et al. 1996):
!
!               w = [sin(pi*x/D)/(pi*x/D)]*[sin(pi*y/D)/(pi*y/D)]
!
!   Note, the first index of pslab varies along the abcissa
!   (or x), and the second index varies along the ordinate (or y).
!
      parameter(beszero=3.8317)
!
      dimension pslab(njx,niy),work(njx,niy),xnu(2),fprint(150,150)

      rmsg=9.0e+9       ! indicates missing data or specification
      if (mod(numpas,100).eq.0) return
!
      if (numpas.le.99) then   ! 9-point smoother
!
      do ipas=1,numpas
!
      do i=1,niy
      do j=1,njx
         work(j,i)=0.
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            work(j,i)=rmsg
         else
            totgive=0.
            if (i.gt.1) then
               if (pslab(j,i-1).ne.rmsg) then
                  give=.125*pslab(j,i)
                  work(j,i-1)=work(j,i-1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i-1).ne.rmsg) then
                     give=.0625*pslab(j,i)
                     work(j-1,i-1)=work(j-1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i-1).ne.rmsg) then
                     give=.0625*pslab(j,i)
                     work(j+1,i-1)=work(j+1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (i.lt.niy) then
               if (pslab(j,i+1).ne.rmsg) then
                  give=.125*pslab(j,i)
                  work(j,i+1)=work(j,i+1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i+1).ne.rmsg) then
                     give=.0625*pslab(j,i)
                     work(j-1,i+1)=work(j-1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i+1).ne.rmsg) then
                     give=.0625*pslab(j,i)
                     work(j+1,i+1)=work(j+1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (j.gt.1) then
               if (pslab(j-1,i).ne.rmsg) then
                  give=.125*pslab(j,i)
                  work(j-1,i)=work(j-1,i)+give
                  totgive=totgive+give
               endif
            endif
            if (j.lt.njx) then
               if (pslab(j+1,i).ne.rmsg) then
                  give=.125*pslab(j,i)
                  work(j+1,i)=work(j+1,i)+give
                  totgive=totgive+give
               endif
            endif
            work(j,i)=work(j,i)+pslab(j,i)-totgive
         endif
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
!
      enddo
!
      elseif (numpas.le.299) then   ! smoother-desmoother
!
      if (numpas.ge.200) then
         nump=numpas-200
         inn=1
      else
         nump=numpas-100
         inn=0
      endif
!
      if (nump.lt.1) return
!
!   Check if any data is missing.
!
      imsg=0
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            imsg=1
            goto 15
         endif
      enddo
      enddo
 15   continue
!
      if (imsg.eq.1) then
!
!   Get average value of pslab.
!
      nval=0
      total=0.
      do 10 i=1,niy
      do 10 j=1,njx
         if (pslab(j,i).ne.rmsg) then
            total=total+pslab(j,i)
            nval=nval+1
         endif
   10 continue
      if (nval.eq.0) then
         write(*,*)'  All elements of this pslab are rmsg.'
         return
      endif
      avgval=total/nval
!
!   Set each element that is currently rmsg to avgval, and
!   keep track of them with work.
!
      do 20 i=1,niy
      do 20 j=1,njx
         if (pslab(j,i).eq.rmsg) then
            pslab(j,i)=avgval
            work(j,i)=1.
         else
            work(j,i)=0.
         endif
   20 continue
!
      endif
!
!     *** Do calculation and put into pslab array.
!
      xnu(1) = 0.50
      xnu(2) = -0.52
      je = njx - 1
      ie = niy - 1
      do 100 ipass = 1,nump*2
         kp=2-mod(ipass,2)
!
!        *** First, smooth in the njx direction.
!
         do 60 j = 2,je
            asv = pslab(j,1)
            do 50 i = 2,ie
               aplus = pslab(j,i+1)
               cell = pslab(j,i)
               pslab(j,i)= pslab(j,i) + xnu(kp)*     &
     &            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   50       continue
   60    continue
!
!        *** Now, smooth in the niy direction.
!
         do 80 i = 2,ie
            asv = pslab(1,i)
            do 70 j = 2,je
               aplus = pslab(j+1,i)
               cell = pslab(j,i)
               pslab(j,i) = pslab(j,i) + xnu(kp)*    &
     &            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   70       continue
   80    continue
!
      if (inn.eq.1) then
!
!      Make non-negative.
!
         do i=1,niy
         do j=1,njx
            pslab(j,i)=max(0.,pslab(j,i))
         enddo
         enddo
      endif
!
  100 continue
!
      if (imsg.eq.1) then
!
!      Set rmsg elements back to rmsg
!
         do 200 i=1,niy
         do 200 j=1,njx
            pslab(j,i)=work(j,i)*rmsg + (1.-work(j,i))*pslab(j,i)
  200    continue
      endif
!
      elseif (numpas.le.599) then   ! weighted smoother
!
      idist=mod(numpas,100)
      if (idist.eq.0) return
      nfp=1+2*idist
      npsq=idist*idist
      if (numpas.le.399) then  ! Cressman function
         do i=1,nfp
         do j=1,nfp
            distsq=(i-idist-1.)**2+(j-idist-1.)**2
            fprint(j,i)=max((npsq-distsq)/(npsq+distsq),0.0)
         enddo
         enddo
      elseif (numpas.le.499) then   ! Circular diffraction function
         do i=1,nfp
         do j=1,nfp
            dist=beszero/idist*sqrt((i-idist-1.)**2+(j-idist-1.)**2)
            if (i.eq.idist+1.and.j.eq.idist+1) then
               fprint(j,i)=.5
            else
               fprint(j,i)=max(0.,bes(dist)/dist)
            endif
         enddo
         enddo
      elseif (numpas.le.599) then   ! Rect. diffraction function
         do i=1,nfp
         do j=1,nfp
            if (j.eq.idist+1) then
               xfac=1.
            else
               xdist=pi/idist*(j-idist-1.)
               xfac=sin(xdist)/xdist
            endif
            if (i.eq.idist+1) then
               yfac=1.
            else
               ydist=pi/idist*(i-idist-1.)
               yfac=sin(ydist)/ydist
            endif
            fprint(j,i)=xfac*yfac
         enddo
         enddo
      endif
!
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).ne.rmsg) then
            tot=0.
            totwt=0.
            is=max(1,i-idist)
            ie=min(niy,i+idist)
            js=max(1,j-idist)
            je=min(njx,j+idist)
            do ireg=is,ie
               ifp=ireg-i+idist+1
            do jreg=js,je
               jfp=jreg-j+idist+1
               if (pslab(jreg,ireg).ne.rmsg) then
                  totwt=totwt+fprint(jfp,ifp)
                  tot=tot+fprint(jfp,ifp)*pslab(jreg,ireg)
               endif
            enddo
            enddo
            work(j,i)=tot/totwt
         else
            work(j,i)=rmsg
         endif
      enddo
      enddo
!
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
!
      endif
!
      return
      end
!=======================================================================================
      function bes(x)
      rint=0.
      do i=1,1000
         u=i*.001-.0005
         rint=rint+sqrt(1.-u*u)*cos(x*u)*.001
      enddo
      bes=2.*x*rint/(4.*atan(1.))
      return
      end
!=======================================================================================
   subroutine gtsobs_time_shift ( times, filename )

   use constants
   use namelist_define
   use obs_define

   implicit none

   character (len=80), intent(in)        :: times
   character (len=80), intent(in)        :: filename

   character          :: year*4, month*2, day*2, hour*2, minute*2
   real, dimension(3) :: dat
   integer            :: length, iost, j

   integer, dimension(2,3) :: tc_time    !2: day, hour
   real, dimension(2)      :: tc_lat, tc_lon
   real                    :: tc_lat_speed, tc_lon_speed   !unit:(degree/minute)
   real                    :: delt_time

   length = len_trim(filename)
   open(39, file=filename(1:length), status='old', form = 'formatted', iostat = iost )
   if( iost .ne. 0 ) then
      write(*,*) 'Cannot find the hurricane track file '//filename(1:length)//'='
      stop
   endif

   do_get_track_loop : do
      read(39, '(a2,1x, a2,3f7.1)', iostat = iost)day, hour, dat(1:3)
      if( iost .ne. 0 ) exit
      if( day == times(9:10) .and. hour== times(12:13) ) then
          backspace(39)
          backspace(39)
          read(39,'(i2, 1x, i2, 2f7.1)')tc_time(1:2,1),tc_lat(1),tc_lon(1)
          read(39,'(i2, 1x, i2, 2f7.1)')tc_time(1:2,2),tc_lat(2),tc_lon(2)
          delt_time = (tc_time(2,2)-tc_time(2,1))*60. +    &               !min
                      (tc_time(1,2)-tc_time(1,1))*60.*24.
          tc_lat_speed = (tc_lat(2)-tc_lat(1))/delt_time
          tc_lon_speed = (tc_lon(2)-tc_lon(1))/delt_time 
          exit
      endif
   end do do_get_track_loop
   close(39)

   read(times,'(8x,i2,1x,i2)')tc_time(1:2,1)
   do j = 1, raw%gts%num
      !read(raw%gts%date(j),'(i4,4(1x,i2))')tc_time(1:5,1)
      read(raw%gts%date(j),'(8x,i2,1x,i2)')tc_time(1:2,1)
      delt_time = (tc_time(1,1) - tc_time(1,3))*24.*60. +  &
                  (tc_time(2,1) - tc_time(2,3))*60.     
      raw%gts%latitude(j) = raw%gts%latitude(j) - delt_time*tc_lat_speed
      raw%gts%longitude(j) = raw%gts%longitude(j) - delt_time*tc_lon_speed
   enddo

   end subroutine gtsobs_time_shift
!=======================================================================================
  subroutine cal_ph( kx, znw, t, qv, pb, mu, mub, phb, ph )

! This subroutine calculates perturbation geopotential
!      along a profile line ( i, j, 1:kx )

  use constants

  integer, intent(in)                  :: kx       ! eta half levels
  real, dimension(kx+1), intent(in)    :: znw      ! eta values on full (w) levels
  real, dimension(kx  ), intent(in)    :: t        ! theta-t0 (k)
  real, dimension(kx  ), intent(in)    :: qv       ! water vapor mixing ratio(kg/kg)
  real, dimension(kx  ), intent(in)    :: pb       ! base state pressure(Pa)
  real, intent(in)                     :: mu
  real, intent(in)                     :: mub
  real, dimension(kx+1), intent(in)    :: phb
  real, dimension(kx+1), intent(out)   :: ph

  real, dimension(kx  )                :: dnw, alb, p, al, alt
  real                                 :: qvf, qvf1, qvf2
  real, dimension(kx+1)   :: pht
!------------------------------------------------------------------------------
!
  dnw(1:kx) = znw( 2:kx+1 ) - znw ( 1:kx )

!  alb(1:kx) = (rd/pr)*(t(1:kx)+to)*(pb(1:kx)/pr)**(-(cp-rd)/cp)

  k    = kx
  qvf1 = 0.5*(qv(k)+qv(k))
  qvf2 = 1./(1.+qvf1)
  qvf1 = qvf1*qvf2

  p(k) = - 0.5*( mu + qvf1*mub )*dnw(k)/qvf2
  qvf  = 1. + (rv/rd)*qv(k)
  alt(k) = (rd/pr)*(t(k)+to)*qvf*(((p(k)+pb(k))/pr)**(-(cp-rd)/cp))
  al(k)  = alt(k) - alb(k)

!  Now, integrate down the column to compute the pressure perturbation, and diagnose the two
!  inverse density fields (total and perturbation).

  do k = kx-1, 1, -1
     qvf1 = 0.5*(qv(k)+qv(k+1))
     qvf2 = 1./(1.+qvf1)
     qvf1 = qvf1*qvf2
     p(k) = p(k+1) - ( mu + qvf1*mub )*( 0.5*(dnw(k)+dnw(k+1)) )/qvf2
     qvf  = 1. + (rv/rd)*qv(k)
     alt(k) = (rd/pr)*(t(k)+to)*qvf*(((p(k)+pb(k))/pr)**(-(cp-rd)/cp))
!     al(k)  = alt(k) - alb(k)
  enddo

!  This is the hydrostatic equation used in the model after the small timesteps.
!  In the model, al (inverse density) is computed from the geopotential.

!  ph(1) = 0.
  pht(1) = phb(1)
  do k = 2, kx+1
!     ph(k) = ph(k-1) - dnw(k-1) * ( (mub+mu)*al(k-1) + mu*alb(k-1) )
     pht(k) = pht(k-1) - dnw(k-1) * alt(k-1)*(mu+mub)
     ph(k) = pht(k) - phb(k)
     if( ph(k) >= 30000000. ) then
       stop 'error in calculating ph'
     endif
  enddo

  return

  end subroutine cal_ph
!=======================================================================================
  subroutine cal_press_from_q( kx, znu, znw, qv, mu, mub, p_top, pres )

! This subroutine calculates perturbation pressure
!      along a profile line ( i, j, 1:kx )

  use constants
  implicit none
  integer, intent(in)                  :: kx       ! eta half levels
  real, dimension(kx  ), intent(in)    :: znu      ! eta values on half levels
  real, dimension(kx+1), intent(in)    :: znw      ! eta values on full (w) levels
  real, dimension(kx  ), intent(in)    :: qv       ! water vapor mixing ratio(kg/kg)
  real, intent(in)                     :: mu
  real, intent(in)                     :: mub
  real, intent(in)                     :: p_top
  real, dimension(kx  ), intent(out)   :: pres

  real, dimension(kx  )                :: dnw, alb
  real                                 :: qvf, qvf1, qvf2
  integer                              :: k
!------------------------------------------------------------------------------
!
  dnw(1:kx) = znw( 2:kx+1 ) - znw ( 1:kx )

!  alb(1:kx) = (rd/pr)*(t(1:kx)+to)*(pb(1:kx)/pr)**(-(cp-rd)/cp)

  k    = kx
  qvf1 = 0.5*(qv(k)+qv(k))
  qvf2 = 1./(1.+qvf1)
  qvf1 = qvf1*qvf2

  pres(k) = - 0.5*( mu + qvf1*mub )*dnw(k)/qvf2

!  Now, integrate down the column to compute the pressure perturbation, and diagnose the two
!  inverse density fields (total and perturbation).

  do k = kx-1, 1, -1
     qvf1 = 0.5*(qv(k)+qv(k+1))
     qvf2 = 1./(1.+qvf1)
     qvf1 = qvf1*qvf2
     pres(k) = pres(k+1) - ( mu + qvf1*mub )*( 0.5*(dnw(k)+dnw(k+1)) )/qvf2
  enddo

  do k = 1, kx
     pres(k) = pres(k) + znu(k)*mub + p_top
  enddo

  return

  end subroutine cal_press_from_q
!=======================================================================================
  subroutine roughness_from_landuse ( mminlu, times, ix, jx, lu_index, rough )

!  calculate rough

  use constants
  implicit none
  integer, intent(in)                   :: ix,jx
  character (len=4)   ,   intent(in)    :: mminlu
  character (len=19)  ,   intent(in)    :: times
  real, dimension(ix,jx),   intent(in)    :: lu_index
  real, dimension(ix,jx),   intent(out)   :: rough

  integer                               :: LS, LC, LI, LUCATS, LUSEAS, &
                                           LUMATCH, year, month, day,  &
                                           julday, Ict, isn, io_error, &
                                           m1, m2, n1, n2
  real                                  :: albd, slmo, sfem
  real(kind=4), dimension(50,2)         :: sfz0
  character (len=4)                     :: LUTYPE
  integer                               :: iost, ltbl

  read(times,'(I4,1x,I2,1X,I2)') year, month, day
  call julian_day (year,month,day,julday, 1)
  isn = 1
  IF(JULDAY < 105 .OR. JULDAY > 288) isn=2

  ltbl = 10
  open(ltbl, file = 'LANDUSE.TBL', status='old', form = 'formatted', iostat = iost )
  if( iost .ne. 0 ) stop 'Cannot find the  landuse file : LANDUSE.TBL !!!!!!'

  LUMATCH=0

  DO
      READ (ltbl,'(A4)', IOSTAT=io_error) LUTYPE
      if (io_error /= 0) exit
      READ (ltbl,*, IOSTAT=io_error) LUCATS,LUSEAS

      IF(LUTYPE == MMINLU) LUMATCH=1

      DO LS=1,LUSEAS
         READ (ltbl,*)
         DO LC=1,LUCATS
            IF(LUTYPE == MMINLU)THEN
               READ (ltbl,*) LI, ALBD, SLMO, SFEM, SFZ0(LC,LS)
               IF(LC /= LI) STOP 'MISSING LANDUSE: LC'
            ELSE
               READ (ltbl,*)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   close (ltbl)

   IF(LUMATCH == 0)THEN
      PRINT *,'LANDUSE IN INPUT FILE DOES NOT MATCH LUTABLE'
      STOP 'INCONSISTENT OR MISSING LUTABLE FILE'
   ENDIF

!  write(0,'(''dimensions of lu_index: m1,m2,n1,n2='',4I4)') m1,m2,n1,n2

   do lc = 1,ix
   do ls = 1,jx
      Li = int(lu_index(lc,ls)+0.001)
      rough(lc,ls) =  sfz0(Li,isn)/100.
   enddo
   enddo

  end subroutine roughness_from_landuse
!=======================================================================================
  subroutine julian_day(NY,NM,ND,JD,METHOD)
   IMPLICIT NONE

   integer,  intent(IN)    :: METHOD, NY
   integer,  intent(INOUT) :: NM, ND, JD

   integer                 :: LOOP
   integer, DIMENSION(12)  :: MDAY = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   IF(MOD(NY,4) == 0) then
      MDAY(2)=29
      IF(MOD(NY,100) == 0) then
         MDAY(2)=28
         IF(MOD(NY,400) == 0) then
            MDAY(2)=29
         ENDIF
      ENDIF
   ENDIF

   IF(METHOD == 1) THEN
      JD=ND
      JuDAY: DO LOOP=1,NM-1
         JD=JD+MDAY(LOOP)
      ENDDO JuDAY
   ELSE IF(METHOD == 2) THEN
      NM=1
      ND=JD
      NYEAR: DO LOOP=1,11
         IF(ND <= MDAY(LOOP)) EXIT NYEAR

         ND=ND-MDAY(LOOP)
         NM=NM+1
      ENDDO NYEAR
   END IF

  end subroutine julian_day
!=======================================================================================
  subroutine da_sfc_pre (psfcm, psm, tsm, qsm, hsm, ho, t0, qvo)
!-----------------------------------------------------------------------------!
!
! Correct pressure between two levels.
!
! Reference: make use of the hydrosatic equation:
!
!  P2 = P1 * exp [-G/R * (z2-z1) / (tv1 + tv2)/2)
!
! Where:
!  z1  = height at level 1
!  z1  = height at level 2
!  tv1 = temperature at level 1
!  tv2 = temperature at level 2
!  P1  = Pressure at level 1
!  P2  = Pressure at level 2
!-----------------------------------------------------------------------------!

  use constants
      IMPLICIT NONE

      REAL, INTENT (out)   :: psfcm   ! model pressure at ho
      REAL, INTENT (in)    :: psm, tsm, qsm

      REAL, INTENT (in)           :: hsm, ho
      REAL, INTENT (in), OPTIONAL :: t0, qvo

      REAL                 :: tvo, tvsm, tv, dz, arg0, arg

!-----------------------------------------------------------------------------!

! 1.  MODEL AND OBSERVATION VIRTUAL TEMPERATURE
! ---------------------------------------------

      tvsm = tsm  * (1. + 0.608 * qsm)
      if (present(t0) .and. present(qvo)) then
        tvo = t0  * (1. + 0.608 * qvo)
      else if (present(t0) .and. .not.present(qvo)) then
        tvo = t0
      else
        tvo = tvsm
      endif

      tv  = 0.5 * (tvsm + tvo)

! 2. HEIGHT DIFFERENCE BEWTEEN MODEL SURFACE AND OBSERVATIONS
! ------------------------------------------------------------

      dz = hsm - ho
      arg0 = dz * g / Rd

! 3.  EXTRAPOLATE PRESSURE OBS TO MODEL SURFACE
! ---------------------------------------------

! ---------------------------------------------

      arg = arg0    / tv

      psfcm = psm * exp (arg)

  end subroutine da_sfc_pre
!=======================================================================================
  subroutine sfc_wtq ( psfc, tg, ps, ts, qs, us, vs, ps2, ts2, qs2, &
                       hs, roughness, xland, u10, v10, t2, q2)

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
!
!                      psim  : mechanical psi at lowlest sigma leve1
!                      psim2 : mechanical psi at 2m
!                      psimz : mechanical psi at 10m
!
!-----------------------------------------------------------------------------!
!  adapted from 3dVAR/WRF July 2005  Zhiyong Meng TAMU

  use constants
      USE WRF_TOOLS

      IMPLICIT NONE

      REAL, INTENT (in)                :: ps , ts , qs , us, vs
      REAL, INTENT (in)                :: ps2, ts2, qs2, psfc, tg
      REAL, INTENT (in)                :: hs, roughness, xland
      REAL, INTENT (out)               :: u10, v10, t2, q2

      REAL                             :: regime            
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

      call da_tp_to_qs( tg, psfc, eg, qg )
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

  end subroutine sfc_wtq
!=======================================================================================
  subroutine da_tp_to_qs( t, p, es, qs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  PURPOSE: Convert T/p to saturation specific humidity.
!
!  METHOD: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
!          Use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) ).
!
!  HISTORY: 10/03/2000 - Creation of F90 version.           Dale Barker
!  MODIFIED: 10/01/2002                                 Wei Huang
!  adapted from 3dvar/WRF  Zhiyong Meng 07/2005 TAMU
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use constants
   IMPLICIT NONE

   real, intent(in)  :: t, p
   real, intent(out) :: es, qs

   REAL              :: t_c              ! T in degreesC.

!------------------------------------------------------------------------------
!  [1.0] Initialise:
!------------------------------------------------------------------------------

   t_c = t - T_frez

!------------------------------------------------------------------------------
!  [2.0] Calculate saturation vapour pressure:
!------------------------------------------------------------------------------

   es = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )

!------------------------------------------------------------------------------
!  [3.0] Calculate saturation specific humidity:
!------------------------------------------------------------------------------

   qs = (1/RvRd) * es / ( p - (1-1/RvRd) * es )

  end subroutine da_tp_to_qs
!=======================================================================================

SUBROUTINE quicksort(n,x,ind)
 
IMPLICIT NONE

REAL, INTENT(IN)  :: x(n)
INTEGER, INTENT(IN OUT)   :: ind(n)
INTEGER, INTENT(IN)    :: n

!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL 
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.
! (OK for N up to about a billon)

!*********************************************************************

INTEGER   :: iu(21), il(21)
INTEGER   :: m, i, j, k, l, ij, it, itt, indx
REAL      :: r
REAL      :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
  ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
  r = r + .0390625
ELSE
  r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + r*(j-i)
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90
END SUBROUTINE quicksort
