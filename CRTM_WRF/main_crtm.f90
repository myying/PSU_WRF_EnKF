!
! parallel CRTM main code for WRF output
!  
!---please set Parameters for WRf output, input file(FILE_NAME) in 3.5 & output file in  6.5.
!!


!!!!CHANGE SATELLITE: n_ch, sat_lon, Sensor_Id


PROGRAM crtm 

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE netcdf
  USE mpi_module
  USE CRTM_Module

  ! Disable all implicit typing
  !  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
  ! ----------
  ! Parameters for WRf output
  ! ----------
   INTEGER, parameter :: xmax=333      !d01:199/224 !d02:201/399 !d03:255/357!d04:429
   INTEGER, parameter :: ymax=222      !d01:149/149 !d02:150/225 !d03:255/357!d04:429
   INTEGER, parameter :: zmax=44       !level range
   INTEGER, parameter :: n_ch=2       !number of channles (10 for GOES-R ABI)
   !INTEGER, parameter :: i_ens=0       !ensemble first number, 0:mean
   !INTEGER, parameter :: t_ens=60      !ensemble last number  
   INTEGER, parameter :: i_yy=2011     !initial year
   INTEGER, parameter :: i_mm=10        !initial month
   INTEGER, parameter :: i_dd=30       !initial day
   INTEGER, parameter :: i_hh=0       !initial hour
   INTEGER, parameter :: i_mn=0       !initial minute
   INTEGER, parameter :: f_yy=2011     !final year
   INTEGER, parameter :: f_mm=10        !final month
   INTEGER, parameter :: f_dd=31       !final day
   INTEGER, parameter :: f_hh=00       !final hour
   INTEGER, parameter :: f_mn=00       !final minutes
   INTEGER, parameter :: hint=3,tint=60       !interval time [minutes] < 1 hour
   INTEGER, parameter :: dom =1        !WRF output domain
   CHARACTER(*), PARAMETER ::DATA_DIR=&
  '/scratch/02135/yingyue/EnKF_OSSE/ERA_201110120000'
  !'/scratch/02135/yingyue/EnKF_OSSE/201110120000/Meteosat7/fc'
!   CHARACTER(*), PARAMETER ::OUTPUT_DIR=DATA_DIR//'/radiance'
   CHARACTER(*), PARAMETER ::OUTPUT_DIR=&
  '/work/02135/yingyue/DYNAMO/EnKF_OSSE/201110120000/truth/BT'
!'/scratch/03154/tg824524/simulation/da_goes/Cq_sc_oei-d1h_mslp_rc05_goes_10minr30t4-300t6/radiance/dtfct/201009170600'

   REAL, PARAMETER :: P1000MB=100000.D0
   REAL, PARAMETER :: R_D=287.D0
   REAL, PARAMETER :: CP=7.D0*R_D/2.D0
   REAL, PARAMETER :: Re=6378000.0
   REAL, PARAMETER :: sat_h=35780000.0
   REAL, PARAMETER :: sat_lon=57.0/180.0*3.14159

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS EXAMPLE ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1  ! 11934=117*102
  INTEGER, PARAMETER :: N_LAYERS    = zmax
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
  INTEGER, PARAMETER :: N_CLOUDS    = zmax*5
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  ! ...but only ONE Sensor at a time
  INTEGER, PARAMETER :: N_SENSORS = 1
 
  ! Test GeometryInfo angles. The test scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis
!  REAL(fp), PARAMETER :: ZENITH_ANGLE = 20.0_fp
!  REAL(fp), PARAMETER :: SCAN_ANGLE   = 10.782   !atan((6378km*PI*20deg/180deg)/(2*35780km)) !26.37293341421_fp
  ! ============================================================================

  ! Variables
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: FILE_NAME
  CHARACTER(256) :: FILE_OUTPUT
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, m, irec

  ! Variables for WRF --- by minamide
  !character(LEN=18)  :: file_date
  !character(LEN=12)  :: dir_date
  !character(LEN=3)   :: file_ens
  INTEGER :: k1, k2
  integer :: x, y, tt, v, z, n, reci, ens, n_ec,yy,mm,dd,hh,mn
  INTEGER :: ncl,icl
  real :: xlat(xmax,ymax)  ! latitude
  real :: xlong(xmax,ymax) ! longitude
  real :: lat(xmax,ymax)   ! in radian
  real :: lon(xmax,ymax)   ! in radian
  real :: p(xmax,ymax,zmax)
  real :: pb(xmax,ymax,zmax)
  real :: pres(xmax,ymax,zmax)
  real :: ph(xmax,ymax,zmax+1)
  real :: phb(xmax,ymax,zmax+1)
  real :: delz(zmax)
  real :: t(xmax,ymax,zmax)
  real :: tk(xmax,ymax,zmax)
  real :: qvapor(xmax,ymax,zmax)
  real :: qcloud(xmax,ymax,zmax)
  real :: qrain(xmax,ymax,zmax)
  real :: qice(xmax,ymax,zmax)
  real :: qsnow(xmax,ymax,zmax)
  real :: qgraup(xmax,ymax,zmax)
  real :: psfc(xmax,ymax)
  real :: hgt(xmax,ymax)
  real :: tsk(xmax,ymax)
  real :: landmask(xmax,ymax)
  real :: Tbsend(xmax,ymax,n_ch)
  real :: Tb(xmax,ymax,n_ch)

  ! ============================================================================

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
  ! ============================================================================

  call parallel_start()

  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  if(my_proc_id==0)  write(*,*) "CRTM ver.",TRIM(Version) 

  ! Get sensor id from user
  ! -----------------------
  !WRITE( *,'(/5x,"Enter sensor id [hirs4_n18, amsua_metop-a, or mhs_n18]:")',ADVANCE='NO' )
  !READ( *,'(a)' ) Sensor_Id
  !Sensor_Id = ADJUSTL(Sensor_Id)
  Sensor_Id = ADJUSTL('mviriNOM_m07') 
  !Sensor_Id = ADJUSTL('abi_gr') 
  !WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)


  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. This initializes the CRTM for the sensors
  !     predefined in the example SENSOR_ID parameter.
  !     NOTE: The coefficient data file path is hard-
  !           wired for this example.
  ! --------------------------------------------------
  if(my_proc_id==0) WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                            ChannelInfo  , &  ! Output
                            IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                            IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                            File_Path='/work/02135/yingyue/code/CRTM/crtm_wrf/coefficients/')
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  ! ============================================================================



  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 3a. Allocate the ARRAYS
  ! -----------------------
  ! Note that only those structure arrays with a channel
  ! dimension are allocated here because we've parameterized
  ! the number of profiles in the N_PROFILES parameter.
  !
  ! Users can make the 
  ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The input FORWARD structure
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 3.5. **** ensember menbers **** by Minamide
  !
  !DATA_DIR = ADJUSTL(DATA_DIR)
  !do ens=i_ens, t_ens
  !write(dir_date,'(a6,2i2.2,a2)')'201009',i_dd,i_hh,'00'
  do yy=i_yy,f_yy
   mm_min = 1
   mm_max = 12
   if(yy.eq.f_yy) mm_max = f_mm
   if(yy.eq.i_yy) mm_min = i_mm
   do mm=mm_min,mm_max
    dd_min = 1
    if(mm.eq.1) dd_max = 31
    if(mm.eq.2) dd_max = 28
    if(mm.eq.3) dd_max = 31
    if(mm.eq.4) dd_max = 30
    if(mm.eq.5) dd_max = 31
    if(mm.eq.6) dd_max = 30
    if(mm.eq.7) dd_max = 31
    if(mm.eq.8) dd_max = 31
    if(mm.eq.9) dd_max = 30
    if(mm.eq.10) dd_max = 31
    if(mm.eq.11) dd_max = 30
    if(mm.eq.12) dd_max = 31
    if(mm.eq.f_mm) dd_max = f_dd
    if(mm.eq.i_mm) dd_min = i_dd
    do dd=dd_min,dd_max
     hh_min = 0
     hh_max = 23
     if(dd .eq. f_dd) hh_max = f_hh
     if(dd .eq. i_dd) hh_min = i_hh
     do hh=hh_min,hh_max,hint
      mn_min = 0
      mn_max = 59
      if((hh .eq. f_hh).and.(dd .eq. f_dd)) mn_max = f_mn
      if((hh .eq. i_hh).and.(dd .eq. i_dd)) mn_min = i_mn
      do mn=mn_min,mn_max,tint
      write(FILE_NAME,'(a,a,i1,a1,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a)')DATA_DIR,'/wrfout_d0',dom,'_',yy,'-',mm,'-',dd,'_',hh,':',mn,':00'
      !write(FILE_NAME,'(a,a,i4,4i2.2,a)')DATA_DIR,'/',yy,mm,dd,hh,mn,'/wrfinput_d01_mean'
      FILE_NAME = ADJUSTL(FILE_NAME)
      if(my_proc_id==0)write(*,*) FILE_NAME

  !============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! Fill the Atm structure array.
  ! NOTE: This is an example program for illustrative purposes only.
  !       Typically, one would not assign the data as shown below,
  !       but rather read it from file

  ! 4a. Atmosphere and Surface input
  ! --------------------------------
  !  CALL Load_wrfAtm_Data()
  !  CALL Load_wrfAtm_Data()

  ! 4a1. Loading Atmosphere and Surface input
  ! --------------------------------
  !   CALL Load_wrf_Data()
  call get_variable2d(FILE_NAME,'XLAT',xmax,ymax,1,xlat)
  call get_variable2d(FILE_NAME,'XLONG',xmax,ymax,1,xlong)
  call get_variable3d(FILE_NAME,'P',xmax,ymax,zmax,1,p)
  call get_variable3d(FILE_NAME,'PB',xmax,ymax,zmax,1,pb)
  call get_variable3d(FILE_NAME,'PH',xmax,ymax,zmax+1,1,ph)
  call get_variable3d(FILE_NAME,'PHB',xmax,ymax,zmax+1,1,phb)
  call get_variable3d(FILE_NAME,'T',xmax,ymax,zmax,1,t)
  call get_variable3d(FILE_NAME,'QVAPOR',xmax,ymax,zmax,1,qvapor)
  call get_variable3d(FILE_NAME,'QCLOUD',xmax,ymax,zmax,1,qcloud)
  call get_variable3d(FILE_NAME,'QRAIN',xmax,ymax,zmax,1,qrain)
  call get_variable3d(FILE_NAME,'QICE',xmax,ymax,zmax,1,qice)
  call get_variable3d(FILE_NAME,'QSNOW',xmax,ymax,zmax,1,qsnow)
  call get_variable3d(FILE_NAME,'QGRAUP',xmax,ymax,zmax,1,qgraup)
  call get_variable2d(FILE_NAME,'PSFC',xmax,ymax,1,psfc)
  call get_variable2d(FILE_NAME,'TSK',xmax,ymax,1,tsk)
  call get_variable2d(FILE_NAME,'HGT',xmax,ymax,1,hgt)
  call get_variable2d(FILE_NAME,'LANDMASK',xmax,ymax,1,landmask)
  lat = xlat/180.0*3.14159
  lon = xlong/180.0*3.14159
  pres = P + PB
  tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/CP) )
  where(qvapor.lt.0.0) qvapor=1.0e-8
  where(qcloud.lt.0.0) qcloud=0.0
  where(qice.lt.0.0) qice=0.0
  where(qrain.lt.0.0) qrain=0.0
  where(qsnow.lt.0.0) qsnow=0.0
  where(qgraup.lt.0.0) qgraup=0.0


  ! 4a2. Converting WRF data for CRTM structure
  ! --------------------------------
  !--- calculating for every grid

  if(mod(ymax,nprocs).eq.0) then
     nyi=ymax/nprocs
  else
     nyi=ymax/nprocs+1
  endif
  ystart=my_proc_id*nyi+1
  yend=min(ymax,(my_proc_id+1)*nyi)
  do y=ystart,yend
  do x=1, xmax
  CALL Convert_wrf_crtm()

!---------------------
!for Debug by Minamide
!---------------------
   !write(*,*) 'lpres',atm(1)%Level_Pressure
   !write(*,*) 'Pres',atm(1)%Pressure
   !write(*,*) 'Temp', atm(1)%Temperature
   !write(*,*) 'H2O', atm(1)%Absorber(:,1)
   !write(*,*) 'delz',delz
   !write(*,*) 'hgt',hgt(x,y)
   !write(*,*) 'ph',ph(x,y,:)
   !write(*,*) 'phb',phb(x,y,:)
   !write(*,*) 'qcloud',qcloud(x,y,:)
   !write(*,*) 'qice',qice(x,y,:)
   !write(*,*) 'qsnow',qsnow(x,y,:)
   !write(*,*) 'qrain',qrain(x,y,:)
   !write(*,*) 'qgraup',qgraup(x,y,:)
   !do z=1,ncl
   !write(*,*)'cloud',atm(1)%Cloud(z)%Type,minval(atm(1)%Cloud(z)%Water_Content),'~',maxval(atm(1)%Cloud(z)%Water_Content)
   !enddo
!--------------------

  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )


  ! 4c. Use the SOI radiative transfer algorithm
  ! --------------------------------------------
  Options%RT_Algorithm_ID = RT_SOI
  ! ============================================================================

  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution , &
                               Options = Options )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================



  ! ============================================================================
  ! 6. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  ! User should read the user guide or the source code of the routine
  ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
  ! select the needed variables for outputs.  These variables are contained
  ! in the structure RTSolution.
  !
  !DO m = 1, N_PROFILES
    !WRITE( *,'(//7x,"Profile ",i0," output for ",a )') n, TRIM(Sensor_Id)
    !DO l = 1, n_Channels
      !WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
      !CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
    !END DO
  !END DO

  !---for file output, edited 2014.9.26
  do l = 1, n_Channels
    do m = 1, N_PROFILES
      Tbsend(x,y,l) = real(RTSolution(l,m)%Brightness_Temperature)
   enddo
  enddo
   !WRITE(*,'(7x,"Profile (",i0,", ",i0,") finished Tb = ",f6.2)')x,y,Tbsend(x,y,2)
  ! ============================================================================

!--- end of x,y-loop
  end do
  end do

  CALL MPI_Allreduce(Tbsend,Tb,xmax*ymax*n_ch,MPI_REAL,MPI_SUM,comm,ierr)

  ! ============================================================================
  !6.5  **** writing the output ****
  !
  if(my_proc_id==0) then
    write(FILE_OUTPUT,'(a,a,i1,a1,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a)')OUTPUT_DIR,'/BT_d0',dom,'_',yy,'-',mm,'-',dd,'_',hh,':',mn,'.bin'
    open(10,file=FILE_OUTPUT,&
           form='unformatted',access='direct',recl=4)
    irec = 0
    do y = 1, ymax
    do x = 1, xmax
      irec= irec +1
      write( 10, rec=irec) xlong(x,y)
    enddo
    enddo
    do y = 1, ymax
    do x = 1, xmax
      irec= irec +1
      write( 10, rec=irec) xlat(x,y)
    enddo
    enddo

    do l = 1, n_ch
      do y = 1, ymax
      do x = 1, xmax
        irec= irec +1
        write( 10, rec=irec) Tb(x,y,l)
      enddo
      enddo
    enddo
    close (10)
  !  initializing the Tbsend fields for Bcast
    Tbsend = 0.0
  endif

  ! ============================================================================
  !  **** initializing all Tb and Tbsend fields ****
  !
  Tb = 0.0
  CALL MPI_BCAST(Tbsend,xmax*ymax*n_ch,MPI_REAL,0,comm,ierr)


  !---end of mn & dd & hh & mm & yy loop
  enddo 
  enddo 
  enddo 
  enddo 
  enddo
  !---end of ens-loop 
  !enddo



  ! ============================================================================
  ! 7. **** DESTROY THE CRTM ****
  !
  if(my_proc_id==0) WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

  call parallel_finish()

  ! ============================================================================
   !---for debug by Minamide
   !write(*,*) 'lpres',atm(1)%Level_Pressure
   !write(*,*) 'Pres',atm(1)%Pressure
   !write(*,*) 'Temp', atm(1)%Temperature
   !write(*,*) 'H2O', atm(1)%Absorber(:,1)
   !write(*,*) 'delz',delz
   !write(*,*) 'hgt',hgt(x,y)
   !write(*,*) 'ph',ph(x,y,:)
   !write(*,*) 'phb',phb(x,y,:)
   !write(*,*) 'qcloud',qcloud(x,y,:)
   !write(*,*) 'qice',qice(x,y,:)
   !write(*,*) 'qsnow',qsnow(x,y,:)
   !write(*,*) 'qrain',qrain(x,y,:)
   !write(*,*) 'qgraup',qgraup(x,y,:)
   !do z=1,ncl
   !write(*,*)
   !'cloud',atm(1)%Cloud(z)%Type,minval(atm(1)%Cloud(z)%Water_Content),'~',maxval(atm(1)%Cloud(z)%Water_Content)
   !enddo
  ! ============================================================================


CONTAINS

  INCLUDE "Convert_wrf_crtm.inc"

END PROGRAM crtm

